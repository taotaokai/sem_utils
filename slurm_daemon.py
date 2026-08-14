import asyncio
import json
import logging
import subprocess
import time

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    handlers=[
        logging.StreamHandler(),
    ],
)

class JobNode:
    def __init__(self, name: str, script_path: str, max_retries: int = 2):
        self.name = name
        self.script_path = script_path
        self.max_retries = max_retries
        self.retries = 0
        self.slurm_id = None
        self.status = "UNSUBMITTED"  # UNSUBMITTED, SUBMITTED, RUNNING, COMPLETED, FAILED, CANCELLED
        self.dependencies: set[str] = set()
        self.failure_reason = "UNKNOWN"


class SlurmDaemon:
    def __init__(
        self,
        workflow_config_path: str,
    ):
        self.nodes: dict[str, JobNode] = {}

        # Central registry tracking banned nodes -> expiration unix timestamp (float)
        self.global_banned_nodes: dict[str, float] = {}

        with open(workflow_config_path, "r") as f:
            config = json.load(f)
        settings = config["settings"]
        self.max_retries = settings.get("max_retries", 2)
        self.check_interval_seconds = settings.get("check_interval_seconds", 60)
        self.max_concurrent_submissions = settings.get("max_concurrent_submissions", 5)
        self.ban_duration_seconds = settings.get("ban_duration_seconds", 1800)
        self.log_file = settings.get("log_file", "slurm_daemon.log")

        self.logger = logging.getLogger(__name__)
        log_format = logging.Formatter('%(asctime)s [%(levelname)s] %(message)s')
        file_handler = logging.FileHandler(self.log_file)
        file_handler.setFormatter(log_format)
        self.logger.addHandler(file_handler)

        self.load_and_validate_workflow(workflow_config_path)

    def load_and_validate_workflow(self, config_path: str):
        with open(config_path, "r") as f:
            config = json.load(f)
        for job_cfg in config["jobs"]:
            node = JobNode(
                name=job_cfg["name"],
                script_path=job_cfg["script"],
                max_retries=job_cfg.get("max_retries", self.max_retries),
            )
            node.dependencies = set(job_cfg.get("dependencies", []))
            self.nodes[node.name] = node

        for name, node in self.nodes.items():
            for dep in node.dependencies:
                if dep not in self.nodes:
                    raise ValueError(
                        f"Configuration Error: '{name}' depends on missing job '{dep}'."
                    )

    def get_currently_banned_nodes(self) -> list[str]:
        """Cleans out expired records and returns currently blocked cluster nodes."""
        current_time = time.time()

        # Clean expired records out of the dictionary map
        self.global_banned_nodes = {
            node: expire_time
            for node, expire_time in self.global_banned_nodes.items()
            if expire_time > current_time
        }
        return list(self.global_banned_nodes.keys())

    def submit_with_dependencies(self, node: JobNode) -> str:
        cmd = ["sbatch", "--parsable"]

        valid_dep_ids = [
            self.nodes[p].slurm_id for p in node.dependencies if self.nodes[p].slurm_id
        ]
        if valid_dep_ids:
            cmd.append(f"--dependency=afterok:{':'.join(valid_dep_ids)}")

        # Dynamically evaluate the active ban list at submission time
        active_bans = self.get_currently_banned_nodes()
        if active_bans:
            exclude_str = ",".join(active_bans)
            cmd.append(f"--exclude={exclude_str}")
            self.logger.warning(
                f"Job '{node.name}' submission excluding active banned nodes: {exclude_str}"
            )

        cmd.append(node.script_path)

        try:
            result = subprocess.run(cmd, capture_output=True, text=True, check=True)
            job_id = result.stdout.strip()
            self.logger.info(f"Submitted '{node.name}' -> Slurm ID: {job_id}")
            return job_id
        except subprocess.CalledProcessError as e:
            self.logger.error(
                f"Failed to submit sbatch command for '{node.name}': {e.stderr}"
            )
            return None

    def get_submission_order(self) -> list[str]:
        in_degree = {name: len(node.dependencies) for name, node in self.nodes.items()}
        queue = [name for name, count in in_degree.items() if count == 0]
        order = []
        while queue:
            current = queue.pop(0)
            order.append(current)
            for name, node in self.nodes.items():
                if current in node.dependencies:
                    in_degree[name] -= 1
                    if in_degree[name] == 0:
                        queue.append(name)
        return order

    async def update_job_statuses(self):
        tracked_ids = {
            node.slurm_id: node for node in self.nodes.values() if node.slurm_id
        }
        if not tracked_ids:
            return

        id_str = ",".join(tracked_ids.keys())
        cmd = [
            "sacct",
            "-j",
            id_str,
            "--format=JobID,State,NodeList,ExitCode",
            "-n",
            "-P",
        ]

        try:
            proc = await asyncio.create_subprocess_exec(
                *cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE
            )
            stdout, _ = await proc.communicate()

            for line in stdout.decode().strip().split("\n"):
                if not line or "|" not in line:
                    continue
                parts = line.split("|")
                raw_id, state, node_list, exit_code = (
                    parts[0],
                    parts[1],
                    parts[2],
                    parts[3],
                )

                base_id = raw_id.split(".")[0]
                if base_id in tracked_ids:
                    node = tracked_ids[base_id]
                    self.map_state(node, state)

                    if node.status == "FAILED":
                        self.triage_failure_signature(node, state, node_list, exit_code)
        except Exception as e:
            self.logger.error(f"Error querying Slurm states: {e}")

    def map_state(self, node: JobNode, slurm_state: str):
        if "COMPLETED" in slurm_state:
            node.status = "COMPLETED"
        elif any(
            x in slurm_state
            for x in ["FAILED", "TIMEOUT", "NODE_FAIL", "OUT_OF_MEMORY"]
        ):
            node.status = "FAILED"
        elif any(x in slurm_state for x in ["RUNNING", "COMPLETING"]):
            node.status = "RUNNING"
        elif any(x in slurm_state for x in ["CANCELLED", "REVOKED"]):
            node.status = "CANCELLED"

    def triage_failure_signature(
        self, node: JobNode, state: str, node_list: str, exit_code: str
    ):
        clean_node = (
            node_list.split(",").strip()
            if node_list and node_list not in ["None", "None assigned", ""]
            else None
        )

        if any(x in state for x in ["NODE_FAIL", "TIMEOUT"]):
            node.failure_reason = "CLUSTER_ERROR"
            if clean_node:
                # Add node to central registry with expiration timestamp
                expiration = time.time() + self.ban_duration_seconds
                self.global_banned_nodes[clean_node] = expiration
                self.logger.warning(
                    f"Temporary ban issued for node '{clean_node}' until timestamp {expiration}. Reason: {state}"
                )
            return

        if exit_code and ":" in exit_code:
            try:
                app_code = int(exit_code.split(":")[0])
                if app_code != 0:
                    node.failure_reason = "APPLICATION_ERROR"
                    self.logger.error(
                        f"Application crash caught for '{node.name}' (Exit Code: {app_code}). Bypassing cluster node ban."
                    )
                    return
            except ValueError:
                pass

        node.failure_reason = "CLUSTER_ERROR"
        if clean_node:
            expiration = time.time() + self.ban_duration_seconds
            self.global_banned_nodes[clean_node] = expiration

    def update_downstream_dependencies(self, failed_node: JobNode, new_job_id: str):
        for node in self.nodes.values():
            if failed_node.name in node.dependencies and node.status in [
                "SUBMITTED",
                "RUNNING",
            ]:
                valid_dep_ids = [
                    self.nodes[p].slurm_id
                    for p in node.dependencies
                    if self.nodes[p].slurm_id
                ]
                dep_str = f"afterok:{':'.join(valid_dep_ids)}"
                cmd = [
                    "scontrol",
                    "update",
                    f"JobId={node.slurm_id}",
                    f"Dependency={dep_str}",
                ]
                subprocess.run(cmd, capture_output=True, text=True)

    def can_submit_node(self, node: JobNode) -> bool:
        for parent_name in node.dependencies:
            if self.nodes[parent_name].status == "UNSUBMITTED":
                return False
        return True

    async def run_loop(self, poll_interval: int = 15):
        submission_queue = self.get_submission_order()
        self.logger.info(
            f"Starting scheduling daemon. Default ban time: {self.ban_duration_seconds} seconds."
        )

        while True:
            await self.update_job_statuses()
            active_slots = sum(
                1 for n in self.nodes.values() if n.status in ["SUBMITTED", "RUNNING"]
            )

            # Fill out open sliding window execution channels
            for name in list(submission_queue):
                if active_slots >= self.max_concurrent_submissions:
                    break
                node = self.nodes[name]
                if node.status == "UNSUBMITTED" and self.can_submit_node(node):
                    job_id = self.submit_with_dependencies(node)
                    if job_id:
                        node.slurm_id = job_id
                        node.status = "SUBMITTED"
                        active_slots += 1
                        submission_queue.remove(name)

            # Cycle evaluation loop
            active_pipeline = False
            for name, node in self.nodes.items():
                if node.status == "PENDING":
                    parent_failed = any(
                        self.nodes[p].status in ["FAILED", "CANCELLED"]
                        for p in node.dependencies
                    )
                    if parent_failed:
                        self.logger.warning(
                            f"Job '{node.name}' cancelled due to broken upstream execution tracks."
                        )
                        node.status = "CANCELLED"
                        if name in submission_queue:
                            submission_queue.remove(name)
                        continue

                if node.status == "FAILED":
                    if node.retries < node.max_retries:
                        node.retries += 1
                        self.logger.warning(
                            f"Retrying '{node.name}' (Attempt {node.retries}/{node.max_retries}). Reason: {node.failure_reason}"
                        )

                        node.slurm_id = self.submit_with_dependencies(node)
                        node.status = "SUBMITTED"

                        if node.slurm_id:
                            self.update_downstream_dependencies(node, node.slurm_id)

                    else:
                        self.logger.error(
                            f"Job '{node.name}' permanently halted. Hard failure threshold reached."
                        )
                        continue

                if node.status in ["UNSUBMITTED", "SUBMITTED", "RUNNING"]:
                    active_pipeline = True

            if not active_pipeline and not submission_queue:
                self.logger.info(
                    "All pipeline processes completed successfully. Exiting daemon loop."
                )
                break

            await asyncio.sleep(self.check_interval_seconds)


if __name__ == "__main__":
    # Example initialization: Ban problematic infrastructure nodes for exactly 30 minutes (1800 seconds)
    daemon = SlurmDaemon("workflow.json")
    asyncio.run(daemon.run_loop())

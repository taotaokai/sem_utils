import argparse
import asyncio
import json
import logging
import os
import re
import subprocess
import time
from enum import Enum

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    handlers=[
        logging.StreamHandler(),
    ],
)


class JobStatus(Enum):
    UNSUBMITTED = 1
    SUBMITTED = 2
    RUNNING = 3
    COMPLETED = 4
    FAILED = 5
    CANCELLED = 6


class JobNode:
    def __init__(self, name: str, script_path: str, max_retries: int = 2):
        self.name = name
        self.script_path = script_path
        self.max_retries = max_retries
        self.retries = 0
        self.slurm_id = None
        self.status = (
            JobStatus.UNSUBMITTED
        )  # UNSUBMITTED, SUBMITTED, RUNNING, COMPLETED, FAILED, CANCELLED
        self.dependencies: set[str] = set()
        self.failure_reason = "UNKNOWN"

    def to_dict(self) -> dict:
        return {
            "name": self.name,
            "script_path": self.script_path,
            "max_retries": self.max_retries,
            "retries": self.retries,
            "slurm_id": self.slurm_id,
            "status": self.status.name,
            "dependencies": list(self.dependencies),
            "failure_reason": self.failure_reason,
        }

    def load_dict(self, data: dict):
        self.retries = data["retries"]
        self.slurm_id = data["slurm_id"]
        self.status = JobStatus[data["status"]]
        self.failure_reason = data["failure_reason"]


class SlurmDaemon:
    def __init__(self, workflow_config_path: str, reset_retries: bool = False):
        self.nodes: dict[str, JobNode] = {}

        # Central registry tracking banned nodes -> expiration unix timestamp (float)
        self.global_banned_nodes: dict[
            str, float
        ] = {}  # node_name -> expiration_timestamp
        self.node_failure_tracker: dict[
            str, int
        ] = {}  # node_name -> total_strike_count

        with open(workflow_config_path, "r") as f:
            config = json.load(f)
        settings = config["settings"]
        self.max_retries = settings.get("max_retries", 2)
        self.check_interval_seconds = settings.get("check_interval_seconds", 60)
        self.max_concurrent_submissions = settings.get("max_concurrent_submissions", 5)
        self.ban_duration_seconds = settings.get("ban_duration_seconds", 1800)
        self.max_node_strikes = settings.get("max_node_strikes", 3)
        self.log_file = settings.get("log_file", "slurm_daemon.log")
        self.checkpoint_path = settings.get("checkpoint_path", "checkpoint.json")

        self.logger = logging.getLogger(__name__)
        log_format = logging.Formatter("%(asctime)s [%(levelname)s] %(message)s")
        file_handler = logging.FileHandler(self.log_file)
        file_handler.setFormatter(log_format)
        self.logger.addHandler(file_handler)

        self.load_and_validate_workflow(workflow_config_path, reset_retries)

    def load_and_validate_workflow(self, config_path: str, reset_retries: bool):
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
            if not os.path.isfile(node.script_path):
                raise FileExistsError(f"'{node.script_path}' does not exist.")

        # Step 2: Overlay active state records if checkpoint file exists
        if os.path.exists(self.checkpoint_path):
            self.logger.info(
                f"Checkpoint file detected at '{self.checkpoint_path}'. Restoring pipeline state..."
            )
            try:
                with open(self.checkpoint_path, "r") as f:
                    checkpoint = json.load(f)

                # Restore global temporary infrastructure ban metadata
                self.global_banned_nodes = checkpoint.get("global_banned_nodes", {})
                self.node_failure_tracker = checkpoint.get("node_failure_tracker", {})

                # Restore execution state matrices
                for node_name, node_data in checkpoint.get("nodes", {}).items():
                    if node_name in self.nodes:
                        self.nodes[node_name].load_dict(node_data)
                        if reset_retries:
                            self.nodes[node_name].retries = 0
                        # re-submit previously cancelled jobs
                        if self.nodes[node_name].status == JobStatus.CANCELLED:
                            self.nodes[node_name].status = JobStatus.UNSUBMITTED
                self.logger.info("State pipeline recovery successfully completed.")
            except Exception as e:
                self.logger.error(
                    f"Failed to cleanly load checkpoint file: {e}. Defaulting to fresh execution parameters."
                )
        else:
            self.logger.info(
                "No active checkpoint found. Starting execution graph from scratch."
            )
            self.save_checkpoint()

    def save_checkpoint(self):
        """Atomic write operation preserving the active execution graph blueprint."""
        checkpoint_data = {
            "timestamp": time.time(),
            "global_banned_nodes": self.global_banned_nodes,
            "node_failure_tracker": self.node_failure_tracker,
            "nodes": {name: node.to_dict() for name, node in self.nodes.items()},
        }
        # Writing to temporary file first before overriding prevents corrupted data if the daemon drops mid-write
        temp_path = f"{self.checkpoint_path}.tmp"
        with open(temp_path, "w") as f:
            json.dump(checkpoint_data, f, indent=2)
        os.replace(temp_path, self.checkpoint_path)

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

        if node.dependencies:
            dep_ids = [self.nodes[p].slurm_id for p in node.dependencies]
            # Double check defensive guard: ensure no None values slipped through
            if any(sid is None for sid in dep_ids):
                self.logger.error(
                    f"Attempted to submit '{node.name}' before all parents had Slurm IDs."
                )
                return None
            cmd.append(f"--dependency=afterok:{':'.join(dep_ids)}")

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

    #    def get_submission_order(self) -> list[str]:
    #        in_degree = {name: len(node.dependencies) for name, node in self.nodes.items()}
    #        queue = [name for name, count in in_degree.items() if count == 0]
    #        order = []
    #        while queue:
    #            current = queue.pop(0)
    #            order.append(current)
    #            for name, node in self.nodes.items():
    #                if current in node.dependencies:
    #                    in_degree[name] -= 1
    #                    if in_degree[name] == 0:
    #                        queue.append(name)
    #        return order

    def get_submission_order(self) -> list[str]:
        """Calculates topological sorting execution order using Depth-First Search (DFS)."""
        visited = set()
        temp_stack = set()
        order = []

        def dfs(node_name: str):
            if node_name in temp_stack:
                raise ValueError(
                    f"Configuration Error: Circular dependency detected involving '{node_name}'!"
                )
            if node_name not in visited:
                temp_stack.add(node_name)
                # Find downstream dependent jobs (children)
                for child_name, child_node in self.nodes.items():
                    if node_name in child_node.dependencies:
                        dfs(child_name)

                temp_stack.remove(node_name)
                visited.add(node_name)
                # Insert at the beginning to reverse the post-order collection
                order.insert(0, node_name)

        # Begin DFS from independent root jobs first to build the deepest paths out systematically
        root_nodes = [
            name for name, node in self.nodes.items() if not node.dependencies
        ]
        for root in root_nodes:
            if root not in visited:
                dfs(root)

        # Catch any leftover orphaned structural networks or unvisited islands
        for name in self.nodes:
            if name not in visited:
                dfs(name)

        # Filter submission array down exclusively to items that still need processing
        return [
            name for name in order if self.nodes[name].status == JobStatus.UNSUBMITTED
        ]

    async def update_job_statuses(self):
        # delimiters to split jobID returned from sacct
        delimiters = [".", "_"]  # array jobID is in jobID_taskID
        pattern = "|".join(map(re.escape, delimiters))

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
            "-X",
        ]

        try:
            proc = await asyncio.create_subprocess_exec(
                *cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE
            )
            stdout, _ = await proc.communicate()

            checkpoint_needed = False
            for line in stdout.decode().strip().split("\n"):
                if not line or "|" not in line:
                    continue
                parts = line.split("|")
                if len(parts) < 4:
                    continue
                raw_id, state, node_list, exit_code = (
                    parts[0],
                    parts[1],
                    parts[2],
                    parts[3],
                )

                # base_id = raw_id.split(".")[0]
                # base_id = raw_id.split("_")[0]
                base_id = re.split(pattern, raw_id)[0]
                if base_id in tracked_ids:
                    node = tracked_ids[base_id]
                    old_status = node.status
                    self.map_state(node, state)

                    if (
                        node.status == JobStatus.FAILED
                        and old_status != JobStatus.FAILED
                    ):
                        self.triage_failure_signature(node, state, node_list, exit_code)

                    # Save state immediately if a task crosses into a terminal configuration
                    if old_status != node.status:
                        checkpoint_needed = True

            if checkpoint_needed:
                self.save_checkpoint()

        except Exception as e:
            self.logger.error(f"Error querying Slurm states: {e}")

    def map_state(self, node: JobNode, slurm_state: str):
        if "COMPLETED" in slurm_state:
            node.status = JobStatus.COMPLETED
        elif any(
            x in slurm_state
            for x in ["FAILED", "TIMEOUT", "NODE_FAIL", "OUT_OF_MEMORY"]
        ):
            node.status = JobStatus.FAILED
        elif any(x in slurm_state for x in ["RUNNING", "COMPLETING"]):
            node.status = JobStatus.RUNNING
        elif any(x in slurm_state for x in ["CANCELLED", "REVOKED"]):
            node.status = JobStatus.CANCELLED
        elif any(
            x in slurm_state
            for x in [
                "PENDING",
            ]
        ):
            node.status = JobStatus.SUBMITTED
        else:
            raise ValueError(
                f"unkown slurm state '{slurm_state}' for node '{node.name}'!"
            )

    # def triage_failure_signature(
    #     self, node: JobNode, state: str, node_list: str, exit_code: str
    # ):
    #     clean_node = (
    #         node_list.split(",").strip()
    #         if node_list and node_list not in ["None", "None assigned", ""]
    #         else None
    #     )

    #     if any(x in state for x in ["NODE_FAIL", "TIMEOUT"]):
    #         node.failure_reason = "CLUSTER_ERROR"
    #         if clean_node:
    #             # Add node to central registry with expiration timestamp
    #             expiration = time.time() + self.ban_duration_seconds
    #             self.global_banned_nodes[clean_node] = expiration
    #             self.logger.warning(
    #                 f"Temporary ban issued for node '{clean_node}' until timestamp {expiration}. Reason: {state}"
    #             )
    #         return

    #     if exit_code and ":" in exit_code:
    #         try:
    #             app_code = int(exit_code.split(":")[0])
    #             if app_code != 0:
    #                 node.failure_reason = "APPLICATION_ERROR"
    #                 self.logger.error(
    #                     f"Application crash caught for '{node.name}' (Exit Code: {app_code}). Bypassing cluster node ban."
    #                 )
    #                 return
    #         except ValueError:
    #             pass

    #     node.failure_reason = "CLUSTER_ERROR"
    #     if clean_node:
    #         expiration = time.time() + self.ban_duration_seconds
    #         self.global_banned_nodes[clean_node] = expiration

    def triage_failure_signature(
        self, node_node: JobNode, state: str, node_list: str, exit_code: str
    ):
        """Analyzes failures, expanding Slurm node ranges into individual node strikes."""

        # 1. Isolate verified application errors first
        if exit_code and ":" in exit_code:
            try:
                app_code = int(exit_code.split(":")[0])
                if app_code != 0:
                    node_node.failure_reason = "APPLICATION_ERROR"
                    self.logger.error(
                        f"Application crash (Exit:{app_code}) on '{node_node.name}'. Skipping node tracking."
                    )
                    return
            except ValueError:
                pass

        # 2. Treat as Cluster/Infrastructure Error
        node_node.failure_reason = "CLUSTER_ERROR"
        if not node_list or node_list in ["None", "None assigned", ""]:
            return

        # 3. Expand bracketed Slurm node lists using 'scontrol show hostnames'
        try:
            result = subprocess.run(
                ["scontrol", "show", "hostnames", node_list.strip()],
                capture_output=True,
                text=True,
                check=True,
            )
            expanded_nodes = [
                line.strip()
                for line in result.stdout.strip().split("\n")
                if line.strip()
            ]
        except subprocess.CalledProcessError as e:
            # Fallback to splitting by commas if scontrol encounters an unexpected issue
            self.logger.error(
                f"scontrol failed to parse node list '{node_list}': {e.stderr}"
            )
            expanded_nodes = [n.strip() for n in node_list.split(",") if n.strip()]

        # 4. Process strike counts for each individual node involved in the failure
        for individual_node in expanded_nodes:
            self.node_failure_tracker[individual_node] = (
                self.node_failure_tracker.get(individual_node, 0) + 1
            )
            strikes = self.node_failure_tracker[individual_node]
            self.logger.warning(
                f"Infrastructure failure recorded on node '{individual_node}'. Strike count: {strikes}/{self.max_node_strikes}"
            )

            # If threshold breached, enact time-bounded ban on the specific node
            if strikes > self.max_node_strikes:
                expiration = time.time() + self.ban_duration
                self.global_banned_nodes[individual_node] = expiration
                self.node_failure_tracker[individual_node] = (
                    0  # Reset counter for this node
                )
                self.logger.error(
                    f"Node '{individual_node}' breached failure threshold (> {self.max_node_strikes}). Banned for {self.ban_duration}s."
                )

    # def update_downstream_dependencies(self, failed_node: JobNode, new_job_id: str):
    #     for node in self.nodes.values():
    #         if failed_node.name in node.dependencies and node.status in [
    #             "SUBMITTED",
    #             "RUNNING",
    #         ]:
    #             valid_dep_ids = [
    #                 self.nodes[p].slurm_id
    #                 for p in node.dependencies
    #                 if self.nodes[p].slurm_id
    #             ]
    #             dep_str = f"afterok:{':'.join(valid_dep_ids)}"
    #             cmd = [
    #                 "scontrol",
    #                 "update",
    #                 f"JobId={node.slurm_id}",
    #                 f"Dependency={dep_str}",
    #             ]
    #             subprocess.run(cmd, capture_output=True, text=True)

    def reset_downstream_chain(self, parent_name: str, status: JobStatus):
        """resets jobs to the specified status."""
        for node in self.nodes.values():
            if parent_name in node.dependencies and node.status != status:
                self.logger.warning(
                    f"Resetting downstream job '{node.name}' to status '{status.name}'."
                )
                node.status = JobStatus.UNSUBMITTED
                node.slurm_id = None
                # Cascade recursively down to reset multi-tier children
                self.reset_downstream_chain(node.name)

    def can_submit_node(self, node: JobNode) -> bool:
        for parent_name in node.dependencies:
            parent_node = self.nodes[parent_name]
            if not parent_node.slurm_id or parent_node.status in [
                JobStatus.UNSUBMITTED,
                JobStatus.FAILED,
                JobStatus.CANCELLED,
            ]:
                return False
        return True

    async def run_loop(self):
        # submission_queue = self.get_submission_order()
        self.logger.info(
            f"Starting scheduling daemon. Default ban time: {self.ban_duration_seconds} seconds."
        )

        while True:
            await self.update_job_statuses()
            # Recalculate DFS queue order dynamically to account for rescued/reset jobs
            submission_queue = self.get_submission_order()

            active_slots = sum(
                1
                for n in self.nodes.values()
                if n.status in [JobStatus.SUBMITTED, JobStatus.RUNNING]
            )

            # 1. Slide window allocation deployments
            for name in list(submission_queue):
                if active_slots >= self.max_concurrent_submissions:
                    break
                node = self.nodes[name]
                if node.status == JobStatus.UNSUBMITTED and self.can_submit_node(node):
                    job_id = self.submit_with_dependencies(node)
                    if job_id:
                        node.slurm_id = job_id
                        node.status = JobStatus.SUBMITTED
                        active_slots += 1
                        submission_queue.remove(name)
                        self.save_checkpoint()  # Save checkpoint on successful submissions

            # 2. Process active faults or terminations
            active_pipeline = False
            for name, node in self.nodes.items():
                # if node.status == "PENDING":
                #     parent_failed = any(
                #         self.nodes[p].status in ["FAILED", "CANCELLED"]
                #         for p in node.dependencies
                #     )
                #     if parent_failed:
                #         self.logger.warning(
                #             f"Job '{node.name}' cancelled; parent pipeline broke."
                #         )
                #         node.status = "CANCELLED"
                #         if name in submission_queue:
                #             submission_queue.remove(name)
                #             self.save_checkpoint()
                #         continue

                if node.status == JobStatus.FAILED:
                    if node.retries < node.max_retries:
                        node.retries += 1
                        self.logger.warning(
                            f"Retrying '{node.name}' ({node.retries}/{node.max_retries}). Reason: {node.failure_reason}"
                        )
                        node.slurm_id = self.submit_with_dependencies(node)
                        if node.slurm_id:
                            # Rescue downstream jobs cancelled by Slurm
                            self.reset_downstream_chain(
                                node.name, JobStatus.UNSUBMITTED
                            )
                            node.status = JobStatus.SUBMITTED
                            self.save_checkpoint()
                    else:
                        self.logger.error(
                            f"Job '{node.name}' permanently halted. Retries exhausted."
                        )
                        node.status = JobStatus.CANCELLED
                        self.reset_downstream_chain(node.name, JobStatus.CANCELLED)
                        self.save_checkpoint()
                        continue

                if node.status in [
                    JobStatus.UNSUBMITTED,
                    JobStatus.SUBMITTED,
                    JobStatus.RUNNING,
                ]:
                    active_pipeline = True

            self.save_checkpoint()

            if not active_pipeline and not submission_queue:
                self.logger.info(
                    "All workflow steps finalized. Cleaning up checkpoint file."
                )
                num_failed = sum(
                    [
                        1
                        for node in self.nodes.values()
                        if node.status != JobStatus.COMPLETED
                    ]
                )
                if num_failed == 0 and os.path.exists(self.checkpoint_path):
                    os.remove(self.checkpoint_path)  # Clean up on successful run finish
                break

            await asyncio.sleep(self.check_interval_seconds)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("workflow", type=str, help="json file of workflow definition")
    parser.add_argument(
        "--reset_retries",
        action="store_true",
        help="reset #retries when reading checkpoint",
    )
    args = parser.parse_args()

    daemon = SlurmDaemon(args.workflow, reset_retries=args.reset_retries)
    asyncio.run(daemon.run_loop())

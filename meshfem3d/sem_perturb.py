#!/usr/bin/env python
# -*- coding: utf-8 -*-
import os
import argparse
import numpy as np
from mpi4py import MPI

from meshfem3d_utils import read_gll_file, write_gll_file

# MPI initialization
mpi_comm = MPI.COMM_WORLD
mpi_size = mpi_comm.Get_size()
mpi_rank = mpi_comm.Get_rank()

PERTURB_METHODS = ["absolute", "relative", "exponential"]


def parse_arguments():
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(description="Perturb model")

    parser.add_argument("nproc", type=int, help="number of model slices")
    parser.add_argument(
        "model_dir",
        help="Directory containing input model GLL files (proc*_reg1_[model_tag].bin)",
    )
    parser.add_argument(
        "dmodel_dir",
        help="Directory containing input dmodel GLL files (proc*_reg1_[dmodel_tag].bin)",
    )
    parser.add_argument(
        "out_dir",
        help="Directory for output perturbed model files (proc*_reg1_[model_tag].bin)",
    )
    parser.add_argument(
        "--model_tags",
        nargs="+",
        default=["vp", "vs"],
        help="model tags as in model_dir/proc*_reg1_[model_tag].bin",
    )
    parser.add_argument(
        "--dmodel_tags",
        nargs="+",
        default=["dvp", "dvs"],
        help="dmodel tags as in dmodel_dir/proc*_reg1_[dmodel_tag].bin",
    )
    parser.add_argument(
        "--max_values",
        nargs="+",
        default=None,
        type=float,
        help="output models are clipped above the maximum values",
    )
    parser.add_argument(
        "--min_values",
        nargs="+",
        default=None,
        type=float,
        help="output model are clipped below the minimum values",
    )
    parser.add_argument(
        "--scale",
        default=1.0,
        type=float,
        help="scale factor applied to dmodel",
    )
    parser.add_argument(
        "--method",
        choices=PERTURB_METHODS,
        default="absolute",
        help="perturbation methods: absolute (m + dm), relative (m * (1 + dm)), exponential (m * exp(dm))",
    )

    return parser.parse_args()


def apply_perturbation(
    model_gll: np.ndarray,
    dm_gll: np.ndarray,
    method: str,
    min_value=None,
    max_value=None,
) -> np.ndarray:
    """Apply perturbation based on specified type."""
    if method == "absolute":
        model_out = model_gll + dm_gll
    elif method == "relative":
        model_out = model_gll * (1.0 + dm_gll)
    elif method == "exponential":
        model_out = model_gll * np.exp(dm_gll)
    else:
        raise ValueError(f"Invalid perturbation method: {method}")
    model_out = np.clip(model_out, min_value, max_value)
    return model_out


def process(
    nproc: int,
    model_dir: str,
    dmodel_dir: str,
    out_dir: str,
    model_tags: list,
    dmodel_tags: list,
    scale: float = 1.0,
    method: str = "absolute",
    min_values: list | None = None,
    max_values: list | None = None,
) -> None:
    """
    Process and perturb GLL files.

    Args:
        nproc: Number of model slices
        model_dir: Directory with input model files
        dmodel_dir: Directory with perturbation files
        out_dir: Output directory
        model_tags: List of model parameter tags
        dmodel_tags: List of perturbation parameter tags
        dmodel_scale: Scaling factor for perturbations
        method: Method to apply perturbation
    """

    for iproc in range(mpi_rank, nproc, mpi_size):
        for i, (model_tag, dmodel_tag) in enumerate(zip(model_tags, dmodel_tags)):
            model_gll = read_gll_file(model_dir, model_tag, iproc)
            dm_gll = read_gll_file(dmodel_dir, dmodel_tag, iproc)
            dm_gll = scale * dm_gll  # scale dmodel
            if min_values is not None:
                min_value = min_values[i]
            else:
                min_value = None
            if max_values is not None:
                max_value = max_values[i]
            else:
                max_value = None
            dm_gll = apply_perturbation(
                model_gll, dm_gll, method, min_value=min_value, max_value=max_value
            )
            out_gll = apply_perturbation(model_gll, dm_gll, method)
            write_gll_file(out_dir, model_tag, iproc, out_gll)


def main():
    args = parse_arguments()
    if mpi_rank == 0:
        print(args)

    # Validate arguments
    if not os.path.exists(args.model_dir):
        raise FileNotFoundError(f"model directory not found: {args.model_dir}")
    if not os.path.exists(args.dmodel_dir):
        raise FileNotFoundError(f"dmodel directory not found: {args.dmodel_dir}")

    assert len(args.model_tags) == len(args.dmodel_tags)

    if args.max_values is not None:
        assert len(args.model_tags) == len(args.max_values)
    if args.min_values is not None:
        assert len(args.model_tags) == len(args.min_values)

    # Create output directory if it doesn't exist
    os.makedirs(args.out_dir, exist_ok=True)

    try:
        process(
            args.nproc,
            args.model_dir,
            args.dmodel_dir,
            args.out_dir,
            args.model_tags,
            args.dmodel_tags,
            scale=args.scale,
            method=args.method,
            min_values=args.min_values,
            max_values=args.max_values,
        )
    except Exception as e:
        print(f"Error: {e}")
        raise SystemExit(1)


if __name__ == "__main__":
    main()

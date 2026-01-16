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

PERTURB_TYPES = ["absolute", "relative", "exponential"]


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
        "--model_groups",
        nargs="+",
        default=["vpv,vph", "vsv,vsh", "eta,rho"],
        help="groups of comma separated parameters as in model_dir/proc*_reg1_[model_name].bin",
    )
    parser.add_argument(
        "--group_scales",
        nargs="+",
        default=[1.0,],
        type=float,
        help="scaling factors of dmodel applied to each model group",
    )
    parser.add_argument(
        "--dmodel_tag",
        default="_dmodel",
        help="dmodel tags as in dmodel_dir/proc*_reg1_[model][dmodel_tag].bin",
    )
    parser.add_argument(
        "--models",
        nargs="+",
        default=[],
        help="model names",
    )
    parser.add_argument(
        "--model_min_values",
        nargs="+",
        type=float,
        default=[],
        help="minumum values for each model",
    )
    parser.add_argument(
        "--model_max_values",
        nargs="+",
        type=float,
        default=[],
        help="maximum values for each model",
    )
    parser.add_argument(
        "--method",
        choices=PERTURB_TYPES,
        default="absolute",
        help="perturbation methods: absolute (m + dm), relative (m * (1 + dm)), exponential (m * exp(dm))",
    )
    parser.add_argument(
        "--overwrite_ok",
        action="store_true",
        help="overwrite output files if they already exist",
    )

    return parser.parse_args()


def apply_perturbation(
    model_gll: np.ndarray, dm_gll: np.ndarray, method: str
) -> np.ndarray:
    """Apply perturbation based on specified type."""
    if method == "relative":
        return model_gll * (1.0 + dm_gll)
    elif method == "absolute":
        return model_gll + dm_gll
    elif method == "exponential":
        return model_gll * np.exp(dm_gll)
    else:
        raise ValueError(f"Invalid perturbation method: {method}")


def process(
    nproc,
    model_dir,
    dmodel_dir,
    out_dir,
    model_groups,
    dmodel_scales,
    dmodel_tag="_dmodel",
    method="absolute",
    model_value_ranges={},
    overwrite_ok=False
) -> None:
    """
    Process and perturb GLL files.

    Args:
        nproc: Number of model slices
        model_dir: Directory with input model files
        dmodel_dir: Directory with perturbation files
        out_dir: Output directory
        model_groups: List of model parameter tags
        dmodel_scales: List of perturbation
        dmodel_tag: Tag for dmodel files
        method: Method to apply perturbation
        model_value_ranges: Dictionary of model parameter ranges
    """
    for iproc in range(mpi_rank, nproc, mpi_size):
        original_models = {}
        dmodels = {}
        dmodels_scaled_sum = {}
        for model_names, dmodel_scale in zip(model_groups, dmodel_scales):
            for model_name in model_names:
                if model_name not in original_models:
                    original_models[model_name] = read_gll_file(
                        model_dir, model_name, iproc
                    )
                if model_name not in dmodels:
                    dmodels[model_name] = read_gll_file(dmodel_dir, model_name + dmodel_tag, iproc)
                # cumulative sum in case a model name is present in multiple groups
                if model_name not in dmodels_scaled_sum:
                    dmodels_scaled_sum[model_name] = dmodels[model_name] * dmodel_scale
                else:
                    dmodels_scaled_sum[model_name] += dmodels[model_name] * dmodel_scale

        for model_name in dmodels_scaled_sum:
            perturbed_model = apply_perturbation(
                original_models[model_name], dmodels_scaled_sum[model_name], method
            )
            if model_name in model_value_ranges:
                min_value, max_value = model_value_ranges[model_name]
                perturbed_model = np.clip(perturbed_model, min_value, max_value)
            write_gll_file(out_dir, model_name, iproc, perturbed_model, overwrite=overwrite_ok)


def main():
    args = parse_arguments()
    if mpi_rank == 0:
        print(args)

    # Validate arguments
    if not os.path.exists(args.model_dir):
        raise FileNotFoundError(f"model directory not found: {args.model_dir}")
    if not os.path.exists(args.dmodel_dir):
        raise FileNotFoundError(f"dmodel directory not found: {args.dmodel_dir}")

    group_scales = args.group_scales
    if len(args.group_scales) == 1 and len(args.model_groups) > 1:
        group_scales = args.group_scales * len(args.model_groups)
    assert len(args.model_groups) == len(group_scales)

    assert len(args.models) == len(args.model_min_values) == len(args.model_max_values)
    model_value_ranges = {}
    for model_name, model_min_value, model_max_value in zip(args.models, args.model_min_values, args.model_max_values):
        if model_min_value > model_max_value:
            raise ValueError(f"Invalid model range: {model_name=}, {model_min_value=}, {model_max_value=}")
        model_value_ranges[model_name] = (model_min_value, model_max_value)

    # separate tags in comma separated strings
    # e.g. ["vpv,vph", "vsv,vsh"] -> [["vpv", "vph"], ["vsv", "vsh"]]
    model_groups = [s.split(",") for s in args.model_groups]

    # Create output directory if it doesn't exist
    os.makedirs(args.out_dir, exist_ok=True)

    try:
        process(
            args.nproc,
            args.model_dir,
            args.dmodel_dir,
            args.out_dir,
            model_groups,
            group_scales,
            dmodel_tag=args.dmodel_tag,
            method=args.method,
            model_value_ranges=model_value_ranges,
            overwrite_ok=args.overwrite_ok
        )
    except Exception as e:
        print(f"Error: {e}")
        raise SystemExit(1)


if __name__ == "__main__":
    main()

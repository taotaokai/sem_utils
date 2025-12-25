#!/usr/bin/env python
# -*- coding: utf-8 -*-
import os
import argparse
import numpy as np
from mpi4py import MPI

from meshfem3d_utils import (
    sem_mesh_read,
    sem_mesh_get_vol_gll,
    read_gll_file,
    write_gll_file,
)

# MPI initialization
mpi_comm = MPI.COMM_WORLD
mpi_size = mpi_comm.Get_size()
mpi_rank = mpi_comm.Get_rank()


def parse_arguments():
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(description="Perturb model")

    parser.add_argument("nproc", type=int, help="number of model slices")
    parser.add_argument(
        "mesh_dir", help="directory of mesh files, proc*_reg1_solver_data.bin"
    )
    parser.add_argument(
        "model_dir",
        help="directory of model GLL files, proc*_reg1_[model_tag].bin",
    )
    parser.add_argument(
        "model_tag",
        help="tag as in model_dir/proc*_reg1_[model_tag].bin",
    )
    parser.add_argument(
        "--nbin",
        type=int,
        default=100,
        help="number of bins for histogram",
    )
    parser.add_argument(
        "--exponential_base",
        type=float,
        default=None,
        help="i-th and i+1-th bin edge is at max_amp * exp_base**(nbin-i)",
    )
    parser.add_argument(
        "--cdf_threshold",
        type=float,
        default=None,
        help="amplitude threshold for CDF",
    )
    parser.add_argument(
        "--out_dir",
        default=".",
        help="output directory of truncated GLL files",
    )
    parser.add_argument(
        "--out_tag",
        default=None,
        help="tag for truncated GLL file as proc*_reg1_[out_tag].bin",
    )
    parser.add_argument(
        "--out_hist",
        default="histogram.txt",
        help="output histogram file name",
    )

    return parser.parse_args()


def process(
    nproc,
    mesh_dir,
    model_dir,
    model_tag,
    nbin=100,
    exponential_base=None,
    out_hist="histogram.txt",
    cdf_threshold=None,
    out_dir=".",
    out_tag=None,
):
    """Process and scale GLL files to scaled amplitude."""

    max_amp_l = 0
    for iproc in range(mpi_rank, nproc, mpi_size):
        model_gll = read_gll_file(model_dir, model_tag, iproc)
        max_amp_l = max(max_amp_l, np.max(np.abs(model_gll)))
    max_amp = mpi_comm.allreduce(max_amp_l, op=MPI.MAX)

    bin_edges = np.zeros(nbin + 1)
    if exponential_base is not None:
        bin_edges[1:] = max_amp * np.power(exponential_base, np.arange(nbin)[::-1])
    else:
        bin_edges[:] = np.linspace(0, max_amp, nbin + 1)
    # ensure all values in the first bin are included since (, ] is used
    bin_edges[0] = -1

    pdf_l = np.zeros(nbin)
    for iproc in range(mpi_rank, nproc, mpi_size):
        mesh_file = os.path.join(mesh_dir, f"proc{iproc:06d}_reg1_solver_data.bin")
        mesh_data = sem_mesh_read(mesh_file)
        vol_gll = sem_mesh_get_vol_gll(mesh_data)
        vol_gll = vol_gll.reshape(-1)
        model_gll = read_gll_file(model_dir, model_tag, iproc)
        for ibin in range(nbin):
            mask = (model_gll > bin_edges[ibin]) & (model_gll <= bin_edges[ibin + 1])
            pdf_l[ibin] += np.sum(vol_gll[mask])

    pdf = np.zeros_like(pdf_l)
    mpi_comm.Allreduce(
        [pdf_l, MPI.DOUBLE],  # Send buffer (data and MPI datatype)
        [pdf, MPI.DOUBLE],  # Receive buffer (data and MPI datatype)
        op=MPI.SUM,  # Operation to perform
    )
    pdf /= np.sum(pdf)
    cdf = np.cumsum(pdf)

    if mpi_rank == 0:
        np.savetxt(out_hist, np.column_stack([bin_edges[1:], pdf, cdf]))

    if out_tag is None:
        out_tag = model_tag

    if cdf_threshold is not None:
        if cdf_threshold > (max_cdf := max(cdf)):
            cdf_threshold = max_cdf
        ind = np.where(cdf >= cdf_threshold)[0][0]
        truncate_value = bin_edges[ind + 1]
        for iproc in range(mpi_rank, nproc, mpi_size):
            model_gll = read_gll_file(model_dir, model_tag, iproc)
            model_gll = np.clip(model_gll, -truncate_value, truncate_value)
            write_gll_file(out_dir, model_tag, iproc, model_gll)


def main():
    args = parse_arguments()
    if mpi_rank == 0:
        print(args)

    if args.exponential_base is not None:
        assert args.exponential_base > 0 and args.exponential_base < 1

    try:
        process(
            args.nproc,
            args.mesh_dir,
            args.model_dir,
            args.model_tag,
            nbin=args.nbin,
            exponential_base=args.exponential_base,
            out_hist=args.out_hist,
            cdf_threshold=args.cdf_threshold,
            out_dir=args.out_dir,
            out_tag=args.out_tag,
        )
    except Exception as e:
        print(f"Error: {e}")
        raise SystemExit(1)


if __name__ == "__main__":
    main()

#!/usr/bin/env python
# -*- coding: utf-8 -*-
import os
import argparse
import numpy as np
from scipy.io import FortranFile
from mpi4py import MPI

from meshfem3d_utils import sem_mesh_read, sem_mesh_get_vol_gll

# MPI initialization
mpi_comm = MPI.COMM_WORLD
mpi_size = mpi_comm.Get_size()
mpi_rank = mpi_comm.Get_rank()


def parse_arguments():
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(
        description="Truncate value based on cumulative distribution of GLL values"
        " sum(m * dV, where abs(m) <= m_trunc) / sum(m * dV) = trucate_percentage "
    )

    parser.add_argument("nproc", type=int, help="number of model slices")
    parser.add_argument(
        "mesh_dir", help="directory of mesh files, proc*_reg1_solver_data.bin"
    )
    parser.add_argument(
        "model_dir", help="directory of model files, proc*_reg1_[model_tag].bin"
    )
    parser.add_argument(
        "model_tag", help="tag of model file as proc*_reg1_[model_tag].bin"
    )
    parser.add_argument("out_dir", help="output dir for truncated files")
    parser.add_argument(
        "--nbins",
        type=int,
        help="number of bins for cumulative distribution function (default 100)",
        default=100,
    )
    parser.add_argument(
        "--truncate_percentage",
        type=float,
        default=0.95,
        help="threshold percentage for truncation (default 0.95)",
    )

    return parser.parse_args()


def read_mesh_file(mesh_dir, iproc):
    mesh_file = "%s/proc%06d_reg1_solver_data.bin" % (mesh_dir, iproc)
    mesh_data = sem_mesh_read(mesh_file)

    return mesh_data


def read_gll_file(model_dir, model_tag, iproc, dtype="f4"):
    """Read a Fortran unformatted file and return the data."""
    filename = os.path.join(model_dir, f"proc{iproc:06d}_reg1_{model_tag}.bin")
    if not os.path.exists(filename):
        raise FileNotFoundError(f"File not found: {filename}")

    with FortranFile(filename, "r") as f:
        data = f.read_reals(dtype=dtype)

    return data


def write_gll_file(model_dir, model_tag, iproc, data, dtype="f4"):
    """Write data to a Fortran unformatted file."""
    filename = os.path.join(model_dir, f"proc{iproc:06d}_reg1_{model_tag}.bin")

    with FortranFile(filename, "w") as f:
        f.write_record(np.array(data, dtype=dtype))


def process(nproc, mesh_dir, model_dir, model_tag, nbins=100, truncate_percentage=0.95):
    """Process and truncate GLL files based on cumulative distribution."""

    # get maximum amplitude
    max_amp = -1
    for iproc in range(mpi_rank, nproc, mpi_size):
        model_gll = read_gll_file(model_dir, model_tag, iproc)
        max_amp = max(max_amp, np.max(np.abs(model_gll)))
    max_amp = mpi_comm.allreduce(max_amp, op=MPI.MAX)

    # cumulative distribution by bins
    z_bins = np.linspace(0, max_amp, nbins)
    z_bins[0] = -1
    z_bins[-1] = max_amp + 1
    cdf = np.zeros(nbins)
    for iproc in range(mpi_rank, nproc, mpi_size):
        cdf_local = np.zeros(nbins)
        mesh_data = read_mesh_file(mesh_dir, iproc)
        model_gll = read_gll_file(model_dir, model_tag, iproc)
        model_gll = np.abs(model_gll)
        vol_gll = sem_mesh_get_vol_gll(mesh_data)
        for ibin in range(1, nbins):
            mask = (model_gll >= z_bins[ibin - 1]) & (model_gll < z_bins[ibin + 1])
            cdf_local[ibin] = cdf_local[ibin - 1] + np.sum(
                model_gll[mask] * vol_gll[mask]
            )
        cdf += cdf_local
    mpi_comm.Allreduce(MPI.IN_PLACE, cdf, op=MPI.SUM)
    cdf /= np.sum(cdf)  # normalize
    ind = (cdf >= truncate_percentage).nonzero()[0][0]
    z_cutoff = z_bins[ind]

    # write out threshold gll files
    for iproc in range(mpi_rank, nproc, mpi_size):
        model_gll = read_gll_file(model_dir, model_tag, iproc)
        mask = np.abs(model_gll) > z_cutoff
        model_gll[mask] = z_cutoff * np.sign(model_gll[mask]) # truncate amplitude 
        write_gll_file(model_dir, model_tag, iproc, model_gll)


def main():
    args = parse_arguments()

    # Create output directory if it doesn't exist
    os.makedirs(args.out_dir, exist_ok=True)

    try:
        process(
            args.nproc,
            args.mesh_dir,
            args.model_dir,
            args.model_tag,
            nbins=args.nbins,
            truncate_percentage=args.truncate_percentage,
        )
    except Exception as e:
        print(f"Error: {e}")
        raise SystemExit(1)


if __name__ == "__main__":
    main()

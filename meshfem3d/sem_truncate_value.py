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
    parser.add_argument(
        "--cdf_file",
        default="cdf.txt",
        help="filename to write out cumulative distribution function",
    )
    parser.add_argument(
        "--cdf_only", action="store_true", help="only output cumulative distribution function"
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


def process(nproc, mesh_dir, model_dir, model_tag, out_dir, nbins=100, 
            truncate_percentage=0.95, cdf_file="cdf.txt", cdf_only=False):
    """Process and truncate GLL files based on cumulative distribution."""

    # get maximum amplitude
    max_amp = -1
    for iproc in range(mpi_rank, nproc, mpi_size):
        model_gll = read_gll_file(model_dir, model_tag, iproc)
        max_amp = max(max_amp, np.max(np.abs(model_gll)))
    max_amp = mpi_comm.allreduce(max_amp, op=MPI.MAX)

    # cumulative distribution by bins
    # bin_edges = np.linspace(0, max_amp, nbins + 1)
    bin_edges = np.geomspace(1, max_amp+1, nbins + 1) - 1
    bin_edges[0] -= 1 # ensure all values in the first bin are included since (, ] is used
    print(f"DEBUG: {mpi_rank=}, {bin_edges=}")
    pdf = np.zeros(nbins)
    for iproc in range(mpi_rank, nproc, mpi_size):
        pdf_local = np.zeros(nbins)
        mesh_data = read_mesh_file(mesh_dir, iproc)
        model_gll = read_gll_file(model_dir, model_tag, iproc)
        model_gll = np.abs(model_gll) # amplitude of the model
        vol_gll = sem_mesh_get_vol_gll(mesh_data)
        vol_gll = vol_gll.flatten() # since model_gll is an 1-D array
        for ibin in range(nbins):
            # sum(m * dV), where m in (bin_edges[ibin], bin_edges[ibin+1]]
            mask = (model_gll > bin_edges[ibin]) & (model_gll <= bin_edges[ibin + 1])
            # pdf_local[ibin] = np.sum(model_gll[mask] * vol_gll[mask])
            pdf_local[ibin] = np.sum(vol_gll[mask])
        pdf += pdf_local
    mpi_comm.Allreduce(MPI.IN_PLACE, pdf, op=MPI.SUM)
    pdf /= np.sum(pdf)  # normalize
    cdf = np.cumsum(pdf)
    ind = (cdf >= truncate_percentage).nonzero()[0][0]
    truncate_value = bin_edges[ind+1]

    if mpi_rank == 0:
        np.savetxt(cdf_file, cdf)

    # write out threshold gll files
    if not cdf_only:
        for iproc in range(mpi_rank, nproc, mpi_size):
            model_gll = read_gll_file(model_dir, model_tag, iproc)
            mask = np.abs(model_gll) > truncate_value
            model_gll[mask] = truncate_value * np.sign(model_gll[mask]) # truncate amplitude 
            write_gll_file(out_dir, model_tag, iproc, model_gll)


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
            args.out_dir,
            nbins=args.nbins,
            truncate_percentage=args.truncate_percentage,
            cdf_file=args.cdf_file,
            cdf_only=args.cdf_only,
        )
    except Exception as e:
        print(f"Error: {e}")
        raise SystemExit(1)


if __name__ == "__main__":
    main()

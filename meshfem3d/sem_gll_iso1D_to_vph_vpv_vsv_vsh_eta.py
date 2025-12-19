#!/usr/bin/env python
# -*- coding: utf-8 -*-
import argparse
import os
import numpy as np
from mpi4py import MPI
import numba

from meshfem3d_utils import sem_mesh_read, read_gll_file, write_gll_file
# MPI initialization
comm_world = MPI.COMM_WORLD
mpi_size = comm_world.Get_size()
mpi_rank = comm_world.Get_rank()

def parse_arguments():
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(
        description="Convert 1D isotropy model vph, vpv, vsv, vsh models (alpha, beta, phi, xi = 0; eta = 1)"
    )

    parser.add_argument("nproc", type=int, help="number of mesh slices")
    parser.add_argument(
        "reference_dir", help="directory with 1D reference model proc*_reg1_[vp,vs].bin"
    )
    parser.add_argument(
        "out_dir", help="output dir for proc*_reg1_[vph, vpv, vsv, vsh, eta]_kernel.bin"
    )
    return parser.parse_args()


def read_model_ref(iproc, reference_dir):
    """Read velocity models from a directory.
    note: velocities in km/s
    """
    vp = read_gll_file(reference_dir, "vp", iproc)
    vs = read_gll_file(reference_dir, "vs", iproc)
    return vp, vs


def process(iproc, reference_dir, out_dir):
    """Process model data for a single processor slice."""
    # Read reference velocity models (km/s)
    vp0, vs0 = read_model_ref(iproc, reference_dir)

    # limit model value range
    alpha = np.zeros_like(vp0)
    beta = np.zeros_like(vp0)
    phi = np.zeros_like(vp0)
    xi = np.zeros_like(vp0)
    eta = np.ones_like(vp0)

    # voigt average
    vp = (1+alpha)*vp0
    vs = (1+beta)*vs0

    vpv = (1-4/5*phi)**0.5 * vp
    vph = (1+1/5*phi)**0.5 * vp
    vsv = (1-1/3*xi)**0.5 * vs
    vsh = (1+2/3*xi)**0.5 * vs

    # Write each kernel to a separate file
    write_gll_file(out_dir, "vpv", iproc, vpv)
    write_gll_file(out_dir, "vph", iproc, vph)
    write_gll_file(out_dir, "vsh", iproc, vsh)
    write_gll_file(out_dir, "vsv", iproc, vsv)
    write_gll_file(out_dir, "eta", iproc, eta)


def main():
    """Main function to orchestrate the model conversion process."""
    args = parse_arguments()

    if mpi_rank == 0:
        print(args)

    os.makedirs(args.out_dir, exist_ok=True)

    # Process each processor slice
    for iproc in range(mpi_rank, args.nproc, mpi_size):

        try:
            process(iproc, args.reference_dir, args.out_dir,)
        except Exception as e:
            print(f"Error processing iproc {iproc}: {e}")
            raise SystemExit(1)


if __name__ == "__main__":
    main()

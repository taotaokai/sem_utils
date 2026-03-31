#!/usr/bin/env python
# -*- coding: utf-8 -*-
import argparse
import os
from mpi4py import MPI

from meshfem3d_utils import (
    read_gll_file,
    write_gll_file,
    sem_VTI_model_vpv_vph_vsv_vsh_to_vp0_vs0,
)

# MPI initialization
comm_world = MPI.COMM_WORLD
mpi_size = comm_world.Get_size()
mpi_rank = comm_world.Get_rank()


def parse_arguments():
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(
        description="Convert vph, vpv, vsv, vsh to isotropic vp, vs"
    )
    parser.add_argument("nproc", type=int, help="number of mesh slices")
    parser.add_argument(
        "model_dir", help="directory with proc*_reg1_[vpv,vph,vsv,vsh,eta,rho].bin"
    )
    parser.add_argument(
        "out_dir",
        help="output dir for proc*_reg1_[vp,vs].bin",
    )
    parser.add_argument(
        "--overwrite_ok",
        action="store_true",
        help="allow overwirte existing files",
    )


    return parser.parse_args()


def read_model_tiso(iproc, model_dir):
    """Read velocity models from a directory.
    note: velocities in km/s, density in g/cm^3
    """
    vpv = read_gll_file(model_dir, "vpv", iproc)
    vph = read_gll_file(model_dir, "vph", iproc)
    vsv = read_gll_file(model_dir, "vsv", iproc)
    vsh = read_gll_file(model_dir, "vsh", iproc)
    return vpv, vph, vsv, vsh


def process(iproc, model_dir, out_dir, overwrite_ok):
    """Process model data for a single processor slice."""
    print(f"# iproc {iproc}")

    # velocity models (velocities in km/s, density in g/cm^3)
    vpv, vph, vsv, vsh = read_model_tiso(iproc, model_dir)

    vp0, vs0 = sem_VTI_model_vpv_vph_vsv_vsh_to_vp0_vs0(vpv, vph, vsv, vsh)

    # Write each kernel to a separate file
    write_gll_file(out_dir, "vp", iproc, vp0, overwrite=overwrite_ok)
    write_gll_file(out_dir, "vs", iproc, vs0, overwrite=overwrite_ok)


def main():
    """Main function to orchestrate the model conversion process."""
    args = parse_arguments()

    if mpi_rank == 0:
        print(args)

    os.makedirs(args.out_dir, exist_ok=True)

    # Process each processor slice
    for iproc in range(mpi_rank, args.nproc, mpi_size):
        try:
            process(iproc, args.model_dir, args.out_dir, args.overwrite_ok)
        except Exception as e:
            print(f"Error processing iproc {iproc}: {e}")
            raise SystemExit(1)


if __name__ == "__main__":
    main()

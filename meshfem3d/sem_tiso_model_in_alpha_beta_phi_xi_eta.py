#!/usr/bin/env python
# -*- coding: utf-8 -*-
import argparse
import os
from mpi4py import MPI

from meshfem3d_utils import read_gll_file, write_gll_file

# MPI initialization
comm_world = MPI.COMM_WORLD
mpi_size = comm_world.Get_size()
mpi_rank = comm_world.Get_rank()


def parse_arguments():
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(
        description="Convert TISO model from vpv,vph,vsv,vsh to alpha,beta,phi,xi"
    )

    parser.add_argument("nproc", type=int, help="number of mesh slices")
    parser.add_argument(
        "model_dir", help="directory with proc*_reg1_[vpv,vph,vsv,vsh,eta,rho].bin"
    )
    parser.add_argument(
        "reference_dir", help="directory with reference isotropy model proc*_reg1_[vp,vs].bin"
    )
    parser.add_argument(
        "out_dir", help="output dir for proc*_reg1_[alpha,beta,phi,xi].bin"
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
    # eta = read_gll_file(model_dir, "eta", iproc)
    # rho = read_gll_file(model_dir, "rho", iproc)
    return vpv, vph, vsv, vsh #, eta, rho


def read_model_ref(iproc, reference_dir):
    """Read velocity models from a directory.
    note: velocities in km/s
    """
    vp = read_gll_file(reference_dir, "vp", iproc)
    vs = read_gll_file(reference_dir, "vs", iproc)
    return vp, vs


def convert_model(iproc, model_dir, reference_dir, out_dir):
    """Process model files for a single processor slice."""
    print(f"# iproc {iproc}")

    # Read reference velocity models (km/s)
    vp0, vs0 = read_model_ref(iproc, reference_dir)

    # velocity models (velocities in km/s, density in g/cm^3)
    vpv, vph, vsv, vsh = read_model_tiso(iproc, model_dir)

    # convert to relative velocity perturbation
    vp = ((vpv**2 + 4 * vph**2) / 5) ** 0.5
    vs = ((2 * vsv**2 + vsh**2) / 3) ** 0.5
    alpha = vp / vp0 - 1
    beta = vs / vs0 - 1
    phi = (vph**2 - vpv**2) / vp**2
    xi = (vsh**2 - vsv**2) / vs**2

    # Write each kernel to a separate file
    write_gll_file(out_dir, "alpha", iproc, alpha)
    write_gll_file(out_dir, "beta", iproc, beta)
    write_gll_file(out_dir, "phi", iproc, phi)
    write_gll_file(out_dir, "xi", iproc, xi)


def main():
    """Main function to orchestrate the model conversion process."""
    args = parse_arguments()

    if mpi_rank == 0:
        print(args)

    # Create output directory if it doesn't exist
    os.makedirs(args.out_dir, exist_ok=True)

    # Process each processor slice
    for iproc in range(mpi_rank, args.nproc, mpi_size):

        try:
            convert_model(
                iproc,
                args.model_dir,
                args.refernce_dir,
                args.out_dir,
            )
        except Exception as e:
            print(f"Error processing iproc {iproc}: {e}")
            raise SystemExit(1)


if __name__ == "__main__":
    main()

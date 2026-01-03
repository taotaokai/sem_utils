#!/usr/bin/env python
# -*- coding: utf-8 -*-
import argparse
import os
from mpi4py import MPI
import numpy as np

from meshfem3d_utils import (
    read_gll_file,
    write_gll_file,
    sem_VTI_alpha_beta_phi_xi_to_vpv_vph_vsv_vsh,
    sem_VTI_vpv_vph_vsv_vsh_to_alpha_beta_phi_xi,
)

# MPI initialization
comm_world = MPI.COMM_WORLD
mpi_size = comm_world.Get_size()
mpi_rank = comm_world.Get_rank()


def parse_arguments():
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(
        description="Convert TISO model from vpv,vph,vsv,vsh to alpha,beta,phi,xi,vp,vs"
    )

    parser.add_argument("nproc", type=int, help="number of mesh slices")
    parser.add_argument(
        "reference_dir",
        help="directory with reference isotropy model proc*_reg1_[vp,vs].bin",
    )
    parser.add_argument(
        "model_dir", help="directory with proc*_reg1_[vpv,vph,vsv,vsh].bin"
    )
    parser.add_argument(
        "out_dir", help="output dir for proc*_reg1_[alpha,beta,phi,xi,vp,vs].bin"
    )
    parser.add_argument(
        "--reverse",
        action="store_true",
        help="flag to reverse conversion from alpha,beta,phi,xi to vpv,vph,vsv,vsh",
    )

    return parser.parse_args()


def read_vpv_vph_vsv_vsh(iproc, model_dir):
    """Read velocity models from a directory.
    note: velocities in km/s, density in g/cm^3
    """
    vpv = read_gll_file(model_dir, "vpv", iproc)
    vph = read_gll_file(model_dir, "vph", iproc)
    vsv = read_gll_file(model_dir, "vsv", iproc)
    vsh = read_gll_file(model_dir, "vsh", iproc)
    # eta = read_gll_file(model_dir, "eta", iproc)
    # rho = read_gll_file(model_dir, "rho", iproc)
    return vpv, vph, vsv, vsh  # , eta, rho


def read_alpha_beta_phi_xi(iproc, model_dir):
    """Read velocity models from a directory.
    note: velocities in km/s, density in g/cm^3
    """
    alpha = read_gll_file(model_dir, "alpha", iproc)
    beta = read_gll_file(model_dir, "beta", iproc)
    phi = read_gll_file(model_dir, "phi", iproc)
    xi = read_gll_file(model_dir, "xi", iproc)
    return alpha, beta, phi, xi


def read_model_ref(iproc, reference_dir):
    """Read velocity models from a directory.
    note: velocities in km/s
    """
    vp = read_gll_file(reference_dir, "vp", iproc)
    vs = read_gll_file(reference_dir, "vs", iproc)
    return vp, vs


def convert_model(iproc, model_dir, reference_dir, out_dir, reverse=False):
    """Process model files for a single processor slice."""
    print(f"# iproc {iproc}")

    # Read reference velocity models (km/s)
    vp0, vs0 = read_model_ref(iproc, reference_dir)

    # velocity models (velocities in km/s, density in g/cm^3)
    if reverse:
        alpha, beta, phi, xi = read_alpha_beta_phi_xi(iproc, model_dir)

        vpv, vph, vsv, vsh, vp, vs = sem_VTI_alpha_beta_phi_xi_to_vpv_vph_vsv_vsh(
            alpha, beta, phi, xi, vp0, vs0, output_iso=True
        )

        # Write each model to a separate file
        write_gll_file(out_dir, "vpv", iproc, vpv)
        write_gll_file(out_dir, "vph", iproc, vph)
        write_gll_file(out_dir, "vsv", iproc, vsv)
        write_gll_file(out_dir, "vsh", iproc, vsh)
        write_gll_file(out_dir, "vp", iproc, vp)
        write_gll_file(out_dir, "vs", iproc, vs)

    else:
        vpv, vph, vsv, vsh = read_vpv_vph_vsv_vsh(iproc, model_dir)

        alpha, beta, phi, xi, vp, vs = sem_VTI_vpv_vph_vsv_vsh_to_alpha_beta_phi_xi(
            vpv, vph, vsv, vsh, vp0, vs0, output_iso=True
        )

        # Write each model to a separate file
        write_gll_file(out_dir, "alpha", iproc, alpha)
        write_gll_file(out_dir, "beta", iproc, beta)
        write_gll_file(out_dir, "phi", iproc, phi)
        write_gll_file(out_dir, "xi", iproc, xi)
        write_gll_file(out_dir, "vp", iproc, vp)
        write_gll_file(out_dir, "vs", iproc, vs)


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
                args.reference_dir,
                args.out_dir,
                reverse=args.reverse,
            )
        except Exception as e:
            print(f"Error processing iproc {iproc}: {e}")
            raise SystemExit(1)


if __name__ == "__main__":
    main()

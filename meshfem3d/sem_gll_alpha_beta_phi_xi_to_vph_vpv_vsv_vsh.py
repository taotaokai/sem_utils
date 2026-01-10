#!/usr/bin/env python
# -*- coding: utf-8 -*-
import argparse
import os
from mpi4py import MPI

from meshfem3d_utils import (
    read_gll_file,
    write_gll_file,
    sem_VTI_model_alpha_beta_phi_xi_to_vpv_vph_vsv_vsh,
)

# MPI initialization
comm_world = MPI.COMM_WORLD
mpi_size = comm_world.Get_size()
mpi_rank = comm_world.Get_rank()


def parse_arguments():
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(
        description="Convert alpha, beta, xi, phi to vph, vpv, vsv, vsh models"
    )

    parser.add_argument("nproc", type=int, help="number of mesh slices")
    parser.add_argument(
        "reference_dir", help="directory with reference model proc*_reg1_[vp,vs].bin"
    )
    parser.add_argument(
        "model_dir", help="directory with proc*_reg1_[alpha,beta,phi,xi].bin"
    )
    parser.add_argument(
        "out_dir", help="output dir for proc*_reg1_[vpv,vph,vsv,vsh].bin"
    )
    parser.add_argument("min_alpha", type=float)
    parser.add_argument("max_alpha", type=float)
    parser.add_argument("min_beta", type=float)
    parser.add_argument("max_beta", type=float)
    parser.add_argument("min_phi", type=float)
    parser.add_argument("max_phi", type=float)
    parser.add_argument("min_xi", type=float)
    parser.add_argument("max_xi", type=float)
    parser.add_argument("min_eta", type=float)
    parser.add_argument("max_eta", type=float)

    return parser.parse_args()


def read_model(iproc, model_dir):
    # Read velocity perturbation models from a directory.
    alpha = read_gll_file(model_dir, "alpha", iproc)
    beta = read_gll_file(model_dir, "beta", iproc)
    phi = read_gll_file(model_dir, "phi", iproc)
    xi = read_gll_file(model_dir, "xi", iproc)
    eta = read_gll_file(model_dir, "eta", iproc)
    return alpha, beta, phi, xi, eta


def read_model_ref(iproc, reference_dir):
    """Read velocity models from a directory.
    note: velocities in km/s
    """
    vp = read_gll_file(reference_dir, "vp", iproc)
    vs = read_gll_file(reference_dir, "vs", iproc)
    return vp, vs


def process(
    iproc,
    reference_dir,
    model_dir,
    out_dir,
    min_alpha,
    max_alpha,
    min_beta,
    max_beta,
    min_phi,
    max_phi,
    min_xi,
    max_xi,
    min_eta,
    max_eta,
):
    """Process model data for a single processor slice."""
    # velocity perturbation models
    alpha, beta, phi, xi, eta = read_model(iproc, model_dir)

    # limit model value range
    alpha[alpha < min_alpha] = min_alpha
    alpha[alpha > max_alpha] = max_alpha
    beta[beta < min_beta] = min_beta
    beta[beta > max_beta] = max_beta
    phi[phi < min_phi] = min_phi
    phi[phi > max_phi] = max_phi
    xi[xi < min_xi] = min_xi
    xi[xi > max_xi] = max_xi
    eta[eta < min_eta] = min_eta
    eta[eta > max_eta] = max_eta

    # Read reference velocity models (km/s)
    vp0, vs0 = read_model_ref(iproc, reference_dir)

    vpv, vph, vsv, vsh = sem_VTI_model_alpha_beta_phi_xi_to_vpv_vph_vsv_vsh(
        alpha, beta, phi, xi, vp0, vs0
    )

    # Write each kernel to a separate file
    write_gll_file(out_dir, "vpv", iproc, vpv)
    write_gll_file(out_dir, "vph", iproc, vph)
    write_gll_file(out_dir, "vsv", iproc, vsv)
    write_gll_file(out_dir, "vsh", iproc, vsh)
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
            process(
                iproc,
                args.reference_dir,
                args.model_dir,
                args.out_dir,
                args.min_alpha,
                args.max_alpha,
                args.min_beta,
                args.max_beta,
                args.min_phi,
                args.max_phi,
                args.min_xi,
                args.max_xi,
                args.min_eta,
                args.max_eta,
            )
        except Exception as e:
            print(f"Error processing iproc {iproc}: {e}")
            raise SystemExit(1)


if __name__ == "__main__":
    main()

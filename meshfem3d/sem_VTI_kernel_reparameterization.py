#!/usr/bin/env python
# -*- coding: utf-8 -*-
import argparse
import os
from mpi4py import MPI

from meshfem3d_utils import (
    read_gll_file,
    write_gll_file,
    sem_VTI_kernel_cijkl_rho_to_alpha_beta_phi_xi_eta_rho,
    sem_VTI_kernel_cijkl_rho_to_beta_kappa_phi_xi_eta_rho,
)

# MPI initialization
comm_world = MPI.COMM_WORLD
mpi_size = comm_world.Get_size()
mpi_rank = comm_world.Get_rank()

PARAMETERIZATION_TYPES = [
    "alpha_beta_phi_xi",
    "beta_kappa_phi_xi",
]


def parse_arguments():
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(
        description="Convert cijkl and rho kernels to TISO (transversely isotropic) kernels"
    )

    parser.add_argument("nproc", type=int, help="number of mesh slices")
    parser.add_argument(
        "reference_dir", help="directory with reference model proc*_reg1_[vp,vs].bin"
    )
    parser.add_argument(
        "model_dir", help="directory with proc*_reg1_[alpha,beta,phi,xi,eta,rho].bin"
    )
    parser.add_argument(
        "kernel_dir", help="directory with proc*_reg1_[cijkl,rho]_kernel.bin"
    )
    parser.add_argument(
        "out_dir",
        help="output dir for proc*_reg1_[alpha,beta,phi,xi,eta,rho]_kernel.bin",
    )
    parser.add_argument(
        "--type",
        choices=PARAMETERIZATION_TYPES,
        default="alpha_beta_phi_xi",
        help="type of parameterization",
    )

    return parser.parse_args()


def process_kernel_for_alpha_beta_phi_xi(
    iproc, reference_dir, model_dir, kernel_dir, out_dir
):  # , mask=None):
    """Process kernel data for a single processor slice."""

    # Read cijkl,rho kernel
    cijkl_kernel = read_gll_file(kernel_dir, "cijkl_kernel", iproc).reshape((-1, 21))
    rho_kernel = read_gll_file(kernel_dir, "rho_kernel", iproc)

    # Read reference velocity models (km/s)
    vp0 = read_gll_file(reference_dir, "vp", iproc)
    vs0 = read_gll_file(reference_dir, "vs", iproc)

    # read models (velocities in km/s, density in g/cm^3)
    alpha = read_gll_file(model_dir, "alpha", iproc)
    beta = read_gll_file(model_dir, "beta", iproc)
    phi = read_gll_file(model_dir, "phi", iproc)
    xi = read_gll_file(model_dir, "xi", iproc)
    eta = read_gll_file(model_dir, "eta", iproc)
    rho = read_gll_file(model_dir, "rho", iproc)

    K_alpha, K_beta, K_phi, K_xi, K_eta, K_rho = (
        sem_VTI_kernel_cijkl_rho_to_alpha_beta_phi_xi_eta_rho(
            cijkl_kernel, rho_kernel, alpha, beta, phi, xi, eta, rho, vp0, vs0
        )
    )

    # Write each kernel to a separate file
    write_gll_file(out_dir, "alpha_kernel", iproc, K_alpha)
    write_gll_file(out_dir, "beta_kernel", iproc, K_beta)
    write_gll_file(out_dir, "phi_kernel", iproc, K_phi)
    write_gll_file(out_dir, "xi_kernel", iproc, K_xi)
    write_gll_file(out_dir, "eta_kernel", iproc, K_eta)
    # NOTE maybe use rhonotprime to distinguish from original rho kernel?
    write_gll_file(out_dir, "rho_kernel", iproc, K_rho)


def process_kernel_for_beta_kappa_phi_xi(
    iproc, reference_dir, model_dir, kernel_dir, out_dir
):  # , mask=None):
    """Process kernel data for a single processor slice."""

    # Read cijkl,rho kernel
    cijkl_kernel = read_gll_file(kernel_dir, "cijkl_kernel", iproc).reshape((-1, 21))
    rho_kernel = read_gll_file(kernel_dir, "rho_kernel", iproc)

    # Read reference velocity models (km/s)
    vs0 = read_gll_file(reference_dir, "vs", iproc)

    # read models (velocities in km/s, density in g/cm^3)
    beta = read_gll_file(model_dir, "beta", iproc)
    kappa = read_gll_file(model_dir, "kappa", iproc)
    phi = read_gll_file(model_dir, "phi", iproc)
    xi = read_gll_file(model_dir, "xi", iproc)
    eta = read_gll_file(model_dir, "eta", iproc)
    rho = read_gll_file(model_dir, "rho", iproc)

    K_beta, K_kappa, K_phi, K_xi, K_eta, K_rho = (
        sem_VTI_kernel_cijkl_rho_to_beta_kappa_phi_xi_eta_rho(
            cijkl_kernel, rho_kernel, beta, kappa, phi, xi, eta, rho, vs0
        )
    )

    # Write each kernel to a separate file
    write_gll_file(out_dir, "beta_kernel", iproc, K_beta)
    write_gll_file(out_dir, "kappa_kernel", iproc, K_kappa)
    write_gll_file(out_dir, "phi_kernel", iproc, K_phi)
    write_gll_file(out_dir, "xi_kernel", iproc, K_xi)
    write_gll_file(out_dir, "eta_kernel", iproc, K_eta)
    # NOTE maybe use rhonotprime to distinguish from original rho kernel?
    write_gll_file(out_dir, "rho_kernel", iproc, K_rho)


def main():
    """Main function to orchestrate the kernel conversion process."""
    args = parse_arguments()

    if mpi_rank == 0:
        print(args)

    # Create output directory if it doesn't exist
    os.makedirs(args.out_dir, exist_ok=True)

    # Process each processor slice
    for iproc in range(mpi_rank, args.nproc, mpi_size):
        try:
            if args.type == "alpha_beta_phi_xi":
                process_kernel_for_alpha_beta_phi_xi(
                    iproc,
                    args.reference_dir,
                    args.model_dir,
                    args.kernel_dir,
                    args.out_dir,
                    # mask=mask_gll,
                )
            elif args.type == "beta_kappa_phi_xi":
                process_kernel_for_beta_kappa_phi_xi(
                    iproc,
                    args.reference_dir,
                    args.model_dir,
                    args.kernel_dir,
                    args.out_dir,
                    # mask=mask_gll,
                )
        except Exception as e:
            print(f"Error processing iproc {iproc}: {e}")
            raise SystemExit(1)


if __name__ == "__main__":
    main()

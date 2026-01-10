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
        description="Convert cijkl and rho kernels to TISO (transversely isotropic) kernels"
    )

    parser.add_argument("nproc", type=int, help="number of mesh slices")
    parser.add_argument(
        "model_dir", help="directory with proc*_reg1_[alpha,beta,phi,xi,eta,rho].bin"
    )
    parser.add_argument(
        "reference_dir", help="directory with reference model proc*_reg1_[vp,vs].bin"
    )
    parser.add_argument(
        "kernel_dir", help="directory with proc*_reg1_[cijkl,rho]_kernel.bin"
    )
    parser.add_argument(
        "out_dir",
        help="output dir for proc*_reg1_[alpha,beta,phi,xi,eta,rho]_kernel.bin",
    )

    return parser.parse_args()


def read_model(iproc, model_dir):
    # Read velocity perturbation models from a directory.
    alpha = read_gll_file(model_dir, "alpha", iproc)
    beta = read_gll_file(model_dir, "beta", iproc)
    phi = read_gll_file(model_dir, "phi", iproc)
    xi = read_gll_file(model_dir, "xi", iproc)
    eta = read_gll_file(model_dir, "eta", iproc)
    rho = read_gll_file(model_dir, "rho", iproc)
    return alpha, beta, phi, xi, eta, rho


def read_model_ref(iproc, reference_dir):
    """Read velocity models from a directory.
    note: velocities in km/s
    """
    vp = read_gll_file(reference_dir, "vp", iproc)
    vs = read_gll_file(reference_dir, "vs", iproc)
    return vp, vs


def process_kernel(
    iproc, model_dir, reference_dir, kernel_dir, out_dir
):  # , mask=None):
    """Process kernel data for a single processor slice."""
    print(f"# iproc {iproc}")

    # About the kernel dimension
    # base on specfem3D_glob/src/specfem3D/compute_kernels.F90:
    # subroutine save_kernels_crust_mantle_ani()
    # ! kernel unit [ s / km^3 ]                ! for objective function
    # ! For anisotropic kernels

    # ! final unit : [s km^(-3) GPa^(-1)]       ! for cijkl_kernel
    # ! final unit : [s km^(-3) (kg/m^3)^(-1)]  ! for rho_kernel

    # Read cijkl kernel data
    cijkl_kernel = read_gll_file(kernel_dir, "cijkl_kernel", iproc).reshape((-1, 21))

    # Read rho kernel data
    rho_kernel = read_gll_file(kernel_dir, "rho_kernel", iproc)
    rho_kernel = rho_kernel * 1.0e3  # convert to [s km^(-3) (g/cm^3)^(-1)]

    # Read reference velocity models (km/s)
    vp0, vs0 = read_model_ref(iproc, reference_dir)

    # velocity models (velocities in km/s, density in g/cm^3)
    alpha, beta, phi, xi, eta, rho = read_model(iproc, model_dir)

    # convert to vph, vpv, vsh, vsv
    vpv, vph, vsv, vsh = sem_VTI_model_alpha_beta_phi_xi_to_vpv_vph_vsv_vsh(
        alpha, beta, phi, xi, vp0, vs0
    )

    # voigt-averaged isotropic velocities
    vp = (1.0 + alpha) * vp0
    vs = (1.0 + beta) * vs0

    # Extract specific components
    K_C11 = cijkl_kernel[:, 0]  # dChi/dC11
    K_C12 = cijkl_kernel[:, 1]  # dChi/dC12
    K_C13 = cijkl_kernel[:, 2]  # dChi/dC13
    K_C22 = cijkl_kernel[:, 6]  # dChi/dC22
    K_C23 = cijkl_kernel[:, 7]  # dChi/dC23
    K_C33 = cijkl_kernel[:, 11]  # dChi/dC33
    K_C44 = cijkl_kernel[:, 15]  # dChi/dC44
    K_C55 = cijkl_kernel[:, 18]  # dChi/dC55
    K_C66 = cijkl_kernel[:, 20]  # dChi/dC66

    # Convert to relative velocity using chain rule: dChi/dA = dChi/dCij * dCij/dA

    K_alpha = (
        (
            (K_C11 + K_C22 + K_C12) * (1.0 + 1.0 / 5.0 * phi)
            + K_C33 * (1.0 - 4.0 / 5.0 * phi)
            + (K_C13 + K_C23) * (1.0 + 1.0 / 5.0 * phi) * eta
        )
        * 2.0
        * rho
        * vp0**2
        * (1 + alpha)
    )

    K_beta = (
        (
            K_C12 * (-2.0 * (1.0 + 2.0 / 3.0 * xi))
            + (K_C13 + K_C23) * (-2.0 * eta * (1.0 - 1.0 / 3.0 * xi))
            + (K_C44 + K_C55) * (1.0 - 1.0 / 3.0 * xi)
            + K_C66 * (1.0 + 2.0 / 3.0 * xi)
        )
        * 2.0
        * rho
        * vs0**2
        * (1 + beta)
    )

    K_phi = (
        (
            (K_C11 + K_C22 + K_C12) * 1.0 / 5.0
            + K_C33 * (-4.0 / 5.0)
            + (K_C13 + K_C23) * 1.0 / 5.0 * eta
        )
        * vp**2
        * rho
    )

    K_xi = (
        (
            (K_C44 + K_C55) * (-1.0 / 3.0)
            + K_C66 * (2.0 / 3.0)
            + K_C12 * (-4.0 / 3.0)
            + (K_C13 + K_C23) * (2.0 / 3.0 * eta)
        )
        * vs**2
        * rho
    )

    K_eta = (K_C13 + K_C23) * (vph**2 - 2.0 * vsv**2) * rho

    K_rho = (
        rho_kernel
        + (K_C11 + K_C22) * vph**2
        + K_C33 * vpv**2
        + K_C12 * (vph**2 - 2 * vsh**2)
        + (K_C13 + K_C23) * eta * (vph**2 - vsv**2)
        + (K_C44 + K_C55) * vsv**2
        + K_C66 * vsh**2
    )

    # Write each kernel to a separate file
    write_gll_file(out_dir, "alpha_kernel", iproc, K_alpha)
    write_gll_file(out_dir, "beta_kernel", iproc, K_beta)
    write_gll_file(out_dir, "phi_kernel", iproc, K_phi)
    write_gll_file(out_dir, "xi_kernel", iproc, K_xi)
    write_gll_file(out_dir, "eta_kernel", iproc, K_eta)
    #NOTE maybe use rhonotprime to distinguish from original rho kernel?
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
            process_kernel(
                iproc,
                args.model_dir,
                args.reference_dir,
                args.kernel_dir,
                args.out_dir,
                # mask=mask_gll,
            )
        except Exception as e:
            print(f"Error processing iproc {iproc}: {e}")
            raise SystemExit(1)


if __name__ == "__main__":
    main()

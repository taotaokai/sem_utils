#!/usr/bin/env python
# -*- coding: utf-8 -*-
import argparse
import os
from mpi4py import MPI

from meshfem3d_utils import (
    sem_mesh_read,
    read_gll_file,
    write_gll_file,
    sem_VTI_model_alpha_beta_phi_xi_to_vpv_vph_vsv_vsh,
    sem_VTI_model_vpv_vph_vsv_vsh_to_alpha_beta_phi_xi,
    sem_VTI_model_beta_kappa_phi_xi_to_vpv_vph_vsv_vsh,
    sem_VTI_model_vpv_vph_vsv_vsh_to_beta_kappa_phi_xi,
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
        description="Re-parameterize VTI model from vpv,vph,vsv,vsh to alpha,beta,phi,xi or beta,kappa,phi,xi"
    )

    parser.add_argument(
        "--nproc", type=int, required=True, help="number of mesh slices"
    )
    parser.add_argument(
        "--reference_dir",
        required=True,
        help="directory with reference isotropy model proc*_reg1_[vp,vs].bin",
    )
    parser.add_argument(
        "--model_dir",
        required=True,
        help="directory with proc*_reg1_[vpv,vph,vsv,vsh].bin",
    )
    parser.add_argument(
        "--out_dir",
        required=True,
        help="output dir for proc*_reg1_[alpha,beta,phi,xi,vp,vs].bin",
    )
    parser.add_argument(
        "--reverse",
        action="store_true",
        help="flag to reverse conversion from beta,kappa,phi,xi to vpv,vph,vsv,vsh",
    )
    parser.add_argument(
        "--type",
        choices=PARAMETERIZATION_TYPES,
        default="alpha_beta_phi_xi",
        help="type of parameterization",
    )
    parser.add_argument(
        "--overwrite_ok",
        action="store_true",
        help="overwrite output files if they already exist",
    )
    parser.add_argument(
        "--mesh_dir", default=None, help="mesh directory of proc*_reg1_solver_data.bin"
    )

    return parser.parse_args()


def read_vpv_vph_vsv_vsh_eta(iproc, model_dir, shape=None):
    """Read velocity models from a directory.
    note: velocities in km/s, density in g/cm^3
    """
    vpv = read_gll_file(model_dir, "vpv", iproc, shape=shape)
    vph = read_gll_file(model_dir, "vph", iproc, shape=shape)
    vsv = read_gll_file(model_dir, "vsv", iproc, shape=shape)
    vsh = read_gll_file(model_dir, "vsh", iproc, shape=shape)
    eta = read_gll_file(model_dir, "eta", iproc, shape=shape)
    # rho = read_gll_file(model_dir, "rho", iproc)
    return vpv, vph, vsv, vsh, eta  # , rho


def read_beta_kappa_phi_xi_eta(iproc, model_dir, shape=None):
    """Read velocity models from a directory.
    note: velocities in km/s, density in g/cm^3
    """
    beta = read_gll_file(model_dir, "beta", iproc, shape=shape)
    kappa = read_gll_file(model_dir, "kappa", iproc, shape=shape)
    phi = read_gll_file(model_dir, "phi", iproc, shape=shape)
    xi = read_gll_file(model_dir, "xi", iproc, shape=shape)
    eta = read_gll_file(model_dir, "eta", iproc, shape=shape)
    return beta, kappa, phi, xi, eta


def read_alpha_beta_phi_xi_eta(iproc, model_dir, shape=None):
    """Read velocity models from a directory.
    note: velocities in km/s, density in g/cm^3
    """
    alpha = read_gll_file(model_dir, "alpha", iproc, shape=shape)
    beta = read_gll_file(model_dir, "beta", iproc, shape=shape)
    phi = read_gll_file(model_dir, "phi", iproc, shape=shape)
    xi = read_gll_file(model_dir, "xi", iproc, shape=shape)
    eta = read_gll_file(model_dir, "eta", iproc, shape=shape)
    return alpha, beta, phi, xi, eta


def convert_model_for_beta_kappa_phi_xi(
    iproc,
    model_dir,
    reference_dir,
    out_dir,
    reverse=False,
    overwrite_ok=False,
    mesh_dir=None,
):
    """Process model files for a single processor slice."""
    print(f"# iproc {iproc}")

    # read mesh
    gll_dims = None
    ispec_is_iso = None
    if mesh_dir is not None:
        mesh_file = os.path.join(mesh_dir, f"{iproc:06d}_reg1_solver_data.bin")
        mesh = sem_mesh_read(mesh_file)
        gll_dims = mesh["gll_dims"]
        ispec_is_iso = ~(mesh["ispec_is_tiso"].astype(bool))

    # Read reference velocity models (km/s)
    vs0 = read_gll_file(reference_dir, "vs", iproc, shape=gll_dims)

    # velocity models (velocities in km/s, density in g/cm^3)
    if reverse:
        beta, kappa, phi, xi, eta = read_beta_kappa_phi_xi_eta(
            iproc, model_dir, shape=gll_dims
        )

        if mesh_dir:
            phi[ispec_is_iso, ...] = 0.0
            xi[ispec_is_iso, ...] = 0.0
            eta[ispec_is_iso, ...] = 1.0

        vpv, vph, vsv, vsh, vp, vs = sem_VTI_model_beta_kappa_phi_xi_to_vpv_vph_vsv_vsh(
            beta, kappa, phi, xi, vs0, output_iso=True
        )

        # Write each model to a separate file
        write_gll_file(out_dir, "vpv", iproc, vpv, overwrite=overwrite_ok)
        write_gll_file(out_dir, "vph", iproc, vph, overwrite=overwrite_ok)
        write_gll_file(out_dir, "vsv", iproc, vsv, overwrite=overwrite_ok)
        write_gll_file(out_dir, "vsh", iproc, vsh, overwrite=overwrite_ok)
        write_gll_file(out_dir, "eta", iproc, eta, overwrite=overwrite_ok)
        write_gll_file(out_dir, "vp", iproc, vp, overwrite=overwrite_ok)
        write_gll_file(out_dir, "vs", iproc, vs, overwrite=overwrite_ok)

    else:
        vpv, vph, vsv, vsh, eta = read_vpv_vph_vsv_vsh_eta(
            iproc, model_dir, shape=gll_dims
        )

        beta, kappa, phi, xi, vp, vs = (
            sem_VTI_model_vpv_vph_vsv_vsh_to_beta_kappa_phi_xi(
                vpv, vph, vsv, vsh, vs0, output_iso=True
            )
        )

        if mesh_dir:
            phi[ispec_is_iso, ...] = 0.0
            xi[ispec_is_iso, ...] = 0.0
            eta[ispec_is_iso, ...] = 1.0

        # Write each model to a separate file
        write_gll_file(out_dir, "beta", iproc, beta, overwrite=overwrite_ok)
        write_gll_file(out_dir, "kappa", iproc, kappa, overwrite=overwrite_ok)
        write_gll_file(out_dir, "phi", iproc, phi, overwrite=overwrite_ok)
        write_gll_file(out_dir, "xi", iproc, xi, overwrite=overwrite_ok)
        write_gll_file(out_dir, "eta", iproc, eta, overwrite=overwrite_ok)
        write_gll_file(out_dir, "vp", iproc, vp, overwrite=overwrite_ok)
        write_gll_file(out_dir, "vs", iproc, vs, overwrite=overwrite_ok)


def convert_model_for_alpha_beta_phi_xi(
    iproc,
    model_dir,
    reference_dir,
    out_dir,
    reverse=False,
    overwrite_ok=False,
    mesh_dir=None,
):
    """Process model files for a single processor slice."""
    print(f"# iproc {iproc}")

    # read mesh
    gll_dims = None
    ispec_is_iso = None
    if mesh_dir is not None:
        mesh_file = os.path.join(mesh_dir, f"{iproc:06d}_reg1_solver_data.bin")
        mesh = sem_mesh_read(mesh_file)
        gll_dims = mesh["gll_dims"]
        ispec_is_iso = ~(mesh["ispec_is_tiso"].astype(bool))

    # Read reference velocity models (km/s)
    vp0 = read_gll_file(reference_dir, "vp", iproc, shape=gll_dims)
    vs0 = read_gll_file(reference_dir, "vs", iproc, shape=gll_dims)

    # velocity models (velocities in km/s, density in g/cm^3)
    if reverse:
        alpha, beta, phi, xi, eta = read_alpha_beta_phi_xi_eta(
            iproc, model_dir, shape=gll_dims
        )

        if mesh_dir:
            phi[ispec_is_iso, ...] = 0.0
            xi[ispec_is_iso, ...] = 0.0
            eta[ispec_is_iso, ...] = 1.0

        vpv, vph, vsv, vsh, vp, vs = sem_VTI_model_alpha_beta_phi_xi_to_vpv_vph_vsv_vsh(
            alpha, beta, phi, xi, vp0, vs0, output_iso=True
        )

        # Write each model to a separate file
        write_gll_file(out_dir, "vpv", iproc, vpv, overwrite=overwrite_ok)
        write_gll_file(out_dir, "vph", iproc, vph, overwrite=overwrite_ok)
        write_gll_file(out_dir, "vsv", iproc, vsv, overwrite=overwrite_ok)
        write_gll_file(out_dir, "vsh", iproc, vsh, overwrite=overwrite_ok)
        write_gll_file(out_dir, "eta", iproc, eta, overwrite=overwrite_ok)
        write_gll_file(out_dir, "vp", iproc, vp, overwrite=overwrite_ok)
        write_gll_file(out_dir, "vs", iproc, vs, overwrite=overwrite_ok)

    else:
        vpv, vph, vsv, vsh, eta = read_vpv_vph_vsv_vsh_eta(iproc, model_dir)

        alpha, beta, phi, xi, vp, vs = (
            sem_VTI_model_vpv_vph_vsv_vsh_to_alpha_beta_phi_xi(
                vpv, vph, vsv, vsh, vp0, vs0, output_iso=True
            )
        )

        if mesh_dir:
            phi[ispec_is_iso, ...] = 0.0
            xi[ispec_is_iso, ...] = 0.0
            eta[ispec_is_iso, ...] = 1.0

        # Write each model to a separate file
        write_gll_file(out_dir, "alpha", iproc, alpha, overwrite=overwrite_ok)
        write_gll_file(out_dir, "beta", iproc, beta, overwrite=overwrite_ok)
        write_gll_file(out_dir, "phi", iproc, phi, overwrite=overwrite_ok)
        write_gll_file(out_dir, "xi", iproc, xi, overwrite=overwrite_ok)
        write_gll_file(out_dir, "eta", iproc, eta, overwrite=overwrite_ok)
        write_gll_file(out_dir, "vp", iproc, vp, overwrite=overwrite_ok)
        write_gll_file(out_dir, "vs", iproc, vs, overwrite=overwrite_ok)


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
            if args.type == "alpha_beta_phi_xi":
                convert_model_for_alpha_beta_phi_xi(
                    iproc,
                    args.model_dir,
                    args.reference_dir,
                    args.out_dir,
                    reverse=args.reverse,
                    overwrite_ok=args.overwrite_ok,
                    mesh_dir=args.mesh_dir,
                )
            elif args.type == "beta_kappa_phi_xi":
                convert_model_for_beta_kappa_phi_xi(
                    iproc,
                    args.model_dir,
                    args.reference_dir,
                    args.out_dir,
                    reverse=args.reverse,
                    overwrite_ok=args.overwrite_ok,
                    mesh_dir=args.mesh_dir,
                )
        except Exception as e:
            print(f"Error processing iproc {iproc}: {e}")
            raise SystemExit(1)


if __name__ == "__main__":
    main()

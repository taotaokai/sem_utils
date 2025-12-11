#!/usr/bin/env python
# -*- coding: utf-8 -*-
import argparse
import os
import numpy as np
from scipy.io import FortranFile
from mpi4py import MPI
import numba

from meshfem3d_constants import R_EARTH_KM
from meshfem3d_utils import sem_mesh_read, read_gll_file, write_gll_file

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
        "model_dir", help="directory with proc*_reg1_[vpv,vph,vsv,vsh,eta,rho].bin"
    )
    parser.add_argument(
        "reference_dir", help="directory with reference model proc*_reg1_[vp0,vs0].bin"
    )
    parser.add_argument(
        "kernel_dir", help="directory with proc*_reg1_[cijkl,rho]_kernel.bin"
    )
    parser.add_argument(
        "out_dir", help="output dir for proc*_reg1_[vpv,vph,vsv,vsh,eta,rho]_kernel.bin"
    )
    parser.add_argument(
        "--with_mask",
        action="store_true",
        help="flag to apply Gaussian mask around given points",
    )
    parser.add_argument(
        "--mesh_dir",
        help="directory with proc*_reg1_solver_data.bin",
        default="DATABASES_MPI",
    )
    parser.add_argument(
        "--mask_list",
        default="mask_points.txt",
        help="list of x,y,z,sigma_km, where x,y,z are non-dimensionalized by R_EARTH (e.g. 6371 km)",
    )

    return parser.parse_args()


def read_mesh_file(mesh_dir, iproc):
    mesh_file = os.path.join(mesh_dir, f"proc{iproc:06d}_reg1_solver_data.bin")
    mesh_data = sem_mesh_read(mesh_file)

    return mesh_data


def read_list(point_list):
    import pandas as pd

    df = pd.read_csv(point_list, header=None, delimiter=r"\s+")
    xyz_mask = df.iloc[:, :3].to_numpy(dtype=np.float32)
    sigma_mask = df.iloc[:, 3].to_numpy(dtype=np.float32)
    sigma_mask /= R_EARTH_KM  # non-dimensionalized
    return xyz_mask, sigma_mask


@numba.jit(nopython=True, nogil=True)
def make_gaussian_mask(xyz_glob, xyz_mask, sigma_mask):
    nglob = xyz_glob.shape[0]
    nmask = xyz_mask.shape[0]
    weight_mask = np.ones(nglob, dtype=np.float32)

    sigma_squared = sigma_mask**2

    for iglob in range(nglob):
        weight = 1
        for imask in range(nmask):
            weight *= 1.0 - np.exp(
                -0.5
                * sum((xyz_glob[iglob, :] - xyz_mask[imask, :]) ** 2)
                / sigma_squared[imask]
            )
        weight_mask[iglob] = weight

    return weight_mask


def read_model_tiso(iproc, model_dir):
    """Read velocity models from a directory.
    note: velocities in km/s, density in g/cm^3
    """
    vpv = read_gll_file(model_dir, "vpv", iproc)
    vph = read_gll_file(model_dir, "vph", iproc)
    vsv = read_gll_file(model_dir, "vsv", iproc)
    vsh = read_gll_file(model_dir, "vsh", iproc)
    eta = read_gll_file(model_dir, "eta", iproc)
    rho = read_gll_file(model_dir, "rho", iproc)
    return vpv, vph, vsv, vsh, eta, rho


def read_model_tiso_perturb(iproc, model_dir):
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


def process_kernel(iproc, model_dir, reference_dir, kernel_dir, out_dir, mask=None):
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
    vpv, vph, vsv, vsh, eta, rho = read_model_tiso(iproc, model_dir)
    alpha, beta, phi, xi, eta = read_model_tiso_perturb(iproc, model_dir)

    # convert to relative velocity perturbation 
    vp = (1 + alpha) * vp0 
    vs = (1 + beta) * vs0

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

    K_eta = (
        (K_C13 + K_C23) 
        * (vph**2 - 2.0 * vsv**2)
        * rho
    )

    K_rho = (
        rho_kernel
        + (K_C11 + K_C22) * vph**2
        + K_C33 * vpv**2
        + K_C12 * (vph**2 - 2 * vsh**2)
        + (K_C13 + K_C23) * eta * (vph**2 - vsv**2)
        + (K_C44 + K_C55) * vsv**2
        + K_C66 * vsh**2
    )


    if mask is not None:
        mask = mask.flatten()
        K_alpha *= mask
        K_beta *= mask
        K_phi *= mask
        K_xi *= mask
        K_eta *= mask
        K_rho *= mask

    # Write each kernel to a separate file
    write_gll_file(out_dir, "alpha_kernel", iproc, K_alpha)
    write_gll_file(out_dir, "beta_kernel", iproc, K_beta)
    write_gll_file(out_dir, "phi_kernel", iproc, K_phi)
    write_gll_file(out_dir, "xi_kernel", iproc, K_xi)
    write_gll_file(out_dir, "eta_kernel", iproc, K_eta)
    write_gll_file(out_dir, "rho_kernel", iproc, K_rho)


def main():
    """Main function to orchestrate the kernel conversion process."""
    args = parse_arguments()

    if mpi_rank == 0:
        print(args)

    # if mpi_size != args.nproc:
    #     raise ValueError(f"{mpi_size=} != {args.nproc=}")

    # Create output directory if it doesn't exist
    os.makedirs(args.out_dir, exist_ok=True)

    mask_gll = None
    if args.with_mask:
        xyz_mask, sigma_mask = read_list(args.mask_list)

    # Process each processor slice
    for iproc in range(mpi_rank, args.nproc, mpi_size):

        if args.with_mask:
            mesh_data = read_mesh_file(args.mesh_dir, iproc)
            xyz_glob = mesh_data["xyz_glob"]
            ibool = mesh_data["ibool"]
            mask_glob = make_gaussian_mask(xyz_glob, xyz_mask, sigma_mask)
            mask_gll = mask_glob[ibool]

        try:
            process_kernel(
                iproc,
                args.model_dir,
                args.reference_dir,
                args.kernel_dir,
                args.out_dir,
                mask=mask_gll,
            )
        except Exception as e:
            print(f"Error processing iproc {iproc}: {e}")
            raise SystemExit(1)


if __name__ == "__main__":
    main()

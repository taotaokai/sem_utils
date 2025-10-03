#!/usr/bin/env python
# -*- coding: utf-8 -*-
import argparse
import os
import numpy as np
from scipy.io import FortranFile
from mpi4py import MPI
import numba

from meshfem3d_constants import R_EARTH_KM
from meshfem3d_utils import sem_mesh_read

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
        "kernel_dir", help="directory with proc*_reg1_[cijkl,rho]_kernel.bin"
    )
    parser.add_argument(
        "out_dir", help="output dir for proc*_reg1_[vpv,vph,vsv,vsh,eta,rho]_kernel.bin"
    )
    parser.add_argument("--with_mask", action="store_true", help="flag to apply Gaussian mask around given points")
    parser.add_argument("--mesh_dir", help="directory with proc*_reg1_solver_data.bin", default="DATABASES_MPI")
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
    sigma_mask /= R_EARTH_KM # non-dimensionalized
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
                * sum((xyz_glob[iglob, :] - xyz_mask[imask, :])**2)
                / sigma_squared[imask]
            )
        weight_mask[iglob] = weight

    return weight_mask

def read_fortran_file(filename, dtype="f4"):
    """Read a Fortran unformatted file and return the data."""
    with FortranFile(filename, "r") as f:
        data = f.read_reals(dtype=dtype)
    return data

def write_fortran_file(filename, data, dtype="f4"):
    """Write data to a Fortran unformatted file."""
    with FortranFile(filename, "w") as f:
        f.write_record(np.array(data, dtype=dtype))

def read_model_tiso(iproc, model_dir):
    """Read velocity models from a directory.
       note: velocities in km/s, density in g/cm^3
    """
    vpv = read_fortran_file(os.path.join(model_dir, f"proc{iproc:06d}_reg1_vpv.bin"))
    vph = read_fortran_file(os.path.join(model_dir, f"proc{iproc:06d}_reg1_vph.bin"))
    vsv = read_fortran_file(os.path.join(model_dir, f"proc{iproc:06d}_reg1_vsv.bin"))
    vsh = read_fortran_file(os.path.join(model_dir, f"proc{iproc:06d}_reg1_vsh.bin"))
    eta = read_fortran_file(os.path.join(model_dir, f"proc{iproc:06d}_reg1_eta.bin"))
    rho = read_fortran_file(os.path.join(model_dir, f"proc{iproc:06d}_reg1_rho.bin"))
    return vpv, vph, vsv, vsh, eta, rho

def process_kernel(iproc, model_dir, kernel_dir, out_dir, mask=None):
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
    cijkl_file = os.path.join(kernel_dir, f"proc{iproc:06d}_reg1_cijkl_kernel.bin")
    cijkl_kernel = read_fortran_file(cijkl_file, "f4").reshape((-1, 21))

    # Read rho kernel data
    rho_file = os.path.join(kernel_dir, f"proc{iproc:06d}_reg1_rho_kernel.bin")
    rho_kernel = read_fortran_file(rho_file, "f4")
    rho_kernel = rho_kernel * 1.0e3  # convert to [s km^(-3) (g/cm^3)^(-1)]

    # velocity models (velocities in km/s, density in g/cm^3)
    vpv, vph, vsv, vsh, eta, rho = read_model_tiso(iproc, model_dir)

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

    # Convert to Love parameters using chain rule: dChi/dA = dChi/dCij * dCij/dA
    # These represent the sensitivity kernels for the Love parameters
    # C11 = C22 = A, C33 = C
    # C44 = C55 = L, C66 = N
    # C12 = A - 2 * N, C13 = C23 = F
    K_A = K_C11 + K_C22 + K_C12
    K_C = K_C33
    K_N = K_C66 - 2 * K_C12
    K_L = K_C44 + K_C55
    K_F = K_C13 + K_C23

    # Love parameters (GPa)
    A = rho * vph**2
    C = rho * vpv**2
    N = rho * vsh**2
    L = rho * vsv**2
    F = eta * (A - 2 * L)  # F = eta * rho * (vph**2 - 2 * vsv**2)

    # Convert to kernels for relative perturbation in vph, vpv, vsh, vsv, eta and rho
    # e.g. dvph = vph * alpha_h
    K_dvph = 2 * A * (K_A + eta * K_F)
    K_dvpv = 2 * C * K_C
    K_dvsh = 2 * N * K_N
    K_dvsv = 2 * L * (K_L - 2 * eta * K_F)
    K_deta = F * K_F
    K_drho = rho * rho_kernel + A * K_A + C * K_C + N * K_N + L * K_L + F * K_F

    if mask is not None:
        mask = mask.flatten()
        K_dvph *= mask
        K_dvpv *= mask
        K_dvsh *= mask
        K_dvsv *= mask
        K_deta *= mask
        K_drho *= mask

    # Write each kernel to a separate file
    write_fortran_file(os.path.join(out_dir, f"proc{iproc:06d}_reg1_dvph_kernel.bin"), K_dvph)
    write_fortran_file(os.path.join(out_dir, f"proc{iproc:06d}_reg1_dvpv_kernel.bin"), K_dvpv)
    write_fortran_file(os.path.join(out_dir, f"proc{iproc:06d}_reg1_dvsh_kernel.bin"), K_dvsh)
    write_fortran_file(os.path.join(out_dir, f"proc{iproc:06d}_reg1_dvsv_kernel.bin"), K_dvsv)
    write_fortran_file(os.path.join(out_dir, f"proc{iproc:06d}_reg1_deta_kernel.bin"), K_deta)
    write_fortran_file(os.path.join(out_dir, f"proc{iproc:06d}_reg1_drho_kernel.bin"), K_drho)

def main():
    """Main function to orchestrate the kernel conversion process."""
    args = parse_arguments()
    print(args)

    if mpi_size != args.nproc:
        raise ValueError(f"{mpi_size=} != {args.nproc=}")

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
            process_kernel(iproc, args.model_dir, args.kernel_dir, args.out_dir, mask=mask_gll)
        except Exception as e:
            print(f"Error processing iproc {iproc}: {e}")
            raise SystemExit(1)

if __name__ == "__main__":
    main()
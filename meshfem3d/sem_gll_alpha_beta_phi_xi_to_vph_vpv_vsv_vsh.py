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
        description="Convert alpha, beta, xi, phi to vph, vpv, vsv, vsh models"
    )

    parser.add_argument("nproc", type=int, help="number of mesh slices")
    parser.add_argument(
        "reference_dir", help="directory with reference model proc*_reg1_[vp,vs].bin"
    )
    parser.add_argument(
        "model_dir", help="directory with proc*_reg1_[vpv,vph,vsv,vsh,eta,rho].bin"
    )
    parser.add_argument(
        "out_dir", help="output dir for proc*_reg1_[alpha,beta,phi,xi,eta,rho]_kernel.bin"
    )
    parser.add_argument("min_alpha")
    parser.add_argument("max_alpha")
    parser.add_argument("min_beta")
    parser.add_argument("max_beta")
    parser.add_argument("min_phi")
    parser.add_argument("max_phi")
    parser.add_argument("min_xi")
    parser.add_argument("max_xi")
    parser.add_argument("min_eta")
    parser.add_argument("max_eta")

    return parser.parse_args()


def read_fortran_file(filename, dtype="f4"):
    """Read a Fortran unformatted file and return the data."""
    with FortranFile(filename, "r") as f:
        data = f.read_reals(dtype=dtype)
    return data


def write_fortran_file(filename, data, dtype="f4"):
    """Write data to a Fortran unformatted file."""
    with FortranFile(filename, "w") as f:
        f.write_record(np.array(data, dtype=dtype))


def read_perturb_model_tiso(iproc, model_dir):
    # Read velocity perturbation models from a directory.
    alpha = read_fortran_file(os.path.join(model_dir, f"proc{iproc:06d}_reg1_alpha.bin"))
    beta = read_fortran_file(os.path.join(model_dir, f"proc{iproc:06d}_reg1_beta.bin"))
    phi = read_fortran_file(os.path.join(model_dir, f"proc{iproc:06d}_reg1_phi.bin"))
    xi = read_fortran_file(os.path.join(model_dir, f"proc{iproc:06d}_reg1_xi.bin"))
    eta = read_fortran_file(os.path.join(model_dir, f"proc{iproc:06d}_reg1_eta.bin"))
    return alpha, beta, phi, xi, eta


def read_model_ref(iproc, reference_dir):
    """Read velocity models from a directory.
    note: velocities in km/s
    """
    vp = read_fortran_file(os.path.join(reference_dir, f"proc{iproc:06d}_reg1_vp.bin"))
    vs = read_fortran_file(os.path.join(reference_dir, f"proc{iproc:06d}_reg1_vs.bin"))
    return vp, vs


def process_kernel(
        iproc, reference_dir, model_dir, out_dir,
        min_alpha, max_alpha, min_beta, max_beta, min_phi,
        max_phi, min_xi, max_xi, min_eta, max_eta, mask=None
    ):
    """Process model data for a single processor slice."""
    print(f"# iproc {iproc}")

    # velocity perturbation models
    alpha, beta, phi, xi, eta= read_perturb_model_tiso(iproc, model_dir)
    # Read reference velocity models (km/s)
    vp0, vs0 = read_model_ref(iproc, reference_dir)

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

    # voigt average
    vp = (1+alpha)*vp0
    vs = (1+beta)*vs0

    vpv = (1-4/5*phi)**0.5 * vp
    vph = (1+1/5*phi)**0.5 * vp
    vsv = (1-1/3*xi)**0.5 * vs
    vsh = (1+2/3*xi)**0.5 * vs

    # Write each kernel to a separate file
    write_fortran_file(
        os.path.join(out_dir, f"proc{iproc:06d}_reg1_vpv.bin"), vpv
    )
    write_fortran_file(
        os.path.join(out_dir, f"proc{iproc:06d}_reg1_vph.bin"), vph
    )
    write_fortran_file(
        os.path.join(out_dir, f"proc{iproc:06d}_reg1_vsh.bin"), vsh
    )
    write_fortran_file(
        os.path.join(out_dir, f"proc{iproc:06d}_reg1_vsv.bin"), vsv
    )
    write_fortran_file(
        os.path.join(out_dir, f"proc{iproc:06d}_reg1_eta.bin"), eta
    )


def main():
    """Main function to orchestrate the model conversion process."""
    args = parse_arguments()

    if mpi_rank == 0:
        print(args)

    os.makedirs(args.out_dir, exist_ok=True)

    # Process each processor slice
    for iproc in range(mpi_rank, args.nproc, mpi_size):

        try:
            process_kernel(
                iproc, args.reference_dir,args.model_dir, args.out_dit,
                args.min_alpha, args.max_alpha,
                args.min_beta, args.max_beta,
                args.min_phi, args.max_phi,
                args.min_xi, args.max_xi,
                args.min_eta, args.max_eta,
                )
        except Exception as e:
            print(f"Error processing iproc {iproc}: {e}")
            raise SystemExit(1)


if __name__ == "__main__":
    main()
#!/usr/bin/env python
# -*- coding: utf-8 -*-
import numpy as np
import argparse
import mpi4py.MPI as MPI

from meshfem3d_utils import read_gll_file, write_gll_file

# MPI initialization
mpi_comm = MPI.COMM_WORLD
mpi_size = mpi_comm.Get_size()
mpi_rank = mpi_comm.Get_rank()


def parse_arguments():
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(
        description="Evaluate GLL models by math expression"
    )
    parser.add_argument(
        "--nproc", type=int, required=True, help="number of model slices"
    )
    parser.add_argument(
        "--model_dirs",
        nargs="+",
        required=True,
        help="directories of model gll files, proc*_reg1_[model_tag].bin",
    )
    parser.add_argument(
        "--model_tags",
        nargs="+",
        required=True,
        help="model tags, e.g. rho vp",
    )
    parser.add_argument(
        "--math_expr",
        required=True,
        help="math expression, e.g v[0] * v[1]",
    )
    parser.add_argument("--out_dir", required=True, help="output directory for results")
    parser.add_argument(
        "--out_tag",
        required=True,
        help="tag for output GLL file as proc*_reg1_[out_tag].bin",
    )
    parser.add_argument(
        "--overwrite_ok",
        action="store_true",
        help="allow overwrites existing file",
    )

    return parser.parse_args()


def main():
    args = parse_arguments()
    if mpi_rank == 0:
        print(args)

    model_dirs = args.model_dirs
    model_tags = args.model_tags
    ntags = len(model_tags)

    if len(model_dirs) == 1 and ntags > 1:
        model_dirs = model_dirs * len(model_tags)
    elif len(model_tags) == 1 and len(model_dirs) > 1:
        model_tags = model_tags * len(model_dirs)
    elif len(model_dirs) != len(model_tags):
        raise ValueError("number of model directories must equal number of model tags")

    # for iproc in range(args.nproc):
    for iproc in range(mpi_rank, args.nproc, mpi_size):
        v = []
        for i, tag in enumerate(model_tags):
            tag = model_tags[i]
            gll = read_gll_file(model_dirs[i], tag, iproc)
            v.append(gll)

        out_gll = eval(args.math_expr)
        print(
            "[iproc%03d] %s: min %.3e max %.3e mean %.3e std %.3e"
            % (
                iproc,
                args.math_expr,
                np.min(out_gll),
                np.max(out_gll),
                np.mean(out_gll),
                np.std(out_gll),
            )
        )

        if args.out_dir is not None:
            out_tag = args.out_tag if args.out_tag is not None else model_tags[0]
            write_gll_file(
                args.out_dir, out_tag, iproc, out_gll, overwrite=args.overwrite_ok
            )


if __name__ == "__main__":
    main()

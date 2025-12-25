#!/usr/bin/env python
# -*- coding: utf-8 -*-
import argparse

import numpy as np
import numexpr as ne

from meshfem3d_utils import read_gll_file, write_gll_file


def parse_arguments():
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(
        description="Compare GLL files from two directories using a math expression"
    )
    parser.add_argument("nproc", type=int, help="number of model slices")
    parser.add_argument("model_dir1", help="first directory of model files, as a")
    parser.add_argument("model_dir2", help="second directory of model files, as b")
    parser.add_argument(
        "model_name", help="tag for GLL file as proc*_reg1_[model_name].bin"
    )
    parser.add_argument(
        "--math_expr", default="a - b", help="math expression, e.g a - b"
    )
    parser.add_argument("--out_dir", default=None, help="output directory for results")
    parser.add_argument(
        "--out_name",
        default=None,
        help="tag for output GLL file as proc*_reg1_[out_name].bin",
    )

    return parser.parse_args()


def main():
    args = parse_arguments()

    for iproc in range(args.nproc):

        a = read_gll_file(args.model_dir1, args.model_name, iproc)
        b = read_gll_file(args.model_dir2, args.model_name, iproc)
        c = ne.evaluate(args.math_expr)

        print("#proc", iproc)
        print("a: min %10.3e max %10.3e" % (np.min(a), np.max(a)))
        print("b: min %10.3e max %10.3e" % (np.min(b), np.max(b)))
        print(f"{args.math_expr}: ", "min %10.3e max %10.3e" % (np.min(c), np.max(c)))

        if args.out_dir is not None:
            out_name = args.out_name
            if out_name is None:
                out_name = args.model_name
            write_gll_file(args.out_dir, out_name, iproc, c)


if __name__ == "__main__":
    main()

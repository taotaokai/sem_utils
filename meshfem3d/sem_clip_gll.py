#!/usr/bin/env python
# -*- coding: utf-8 -*-
import argparse
import numpy as np

from meshfem3d_utils import read_gll_file, write_gll_file


def parse_arguments():
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(
        description="Clip GLL files by min/max thresholds"
    )

    parser.add_argument("nproc", type=int, help="number of model slices")
    parser.add_argument(
        "model_dir",
        help="Directory of model files",
    )
    parser.add_argument(
        "model_name", help="tag for GLL file as proc*_reg1_[model_name].bin"
    )
    parser.add_argument("out_dir", help="output dir for proc*_reg1_[model_name].bin")
    parser.add_argument("--min", type=float, default=None, help="minimum threshold for GLL file")
    parser.add_argument("--max", type=float, default=None, help="maximum threshold for GLL file")

    return parser.parse_args()


def main():
    args = parse_arguments()

    for iproc in range(args.nproc):
        print("# iproc", iproc)

        model_gll = read_gll_file(args.model_dir, args.model_name, iproc)

        out_gll = np.clip(model_gll, args.min, args.max)
        write_gll_file(args.out_dir, args.model_name, iproc, out_gll)


if __name__ == "__main__":
    main()

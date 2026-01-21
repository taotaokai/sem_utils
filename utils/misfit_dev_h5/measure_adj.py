#!/usr/bin/env python
# -*- coding: utf-8 -*-
import sys
import argparse
from datetime import datetime
from misfit import Misfit

parser = argparse.ArgumentParser()

parser.add_argument("misfit_h5file" )  # "misfit.h5"
parser.add_argument("misfit_parfile")  # "misfit.yaml"
parser.add_argument("solver_parfile")  # "DATA/Par_file"
parser.add_argument("cmt_file"      )  # "DATA/CMTSOLUTION"
parser.add_argument("channel_file"  )  # "data/C202006211907A/channel.txt"
parser.add_argument("data_h5file"   )  # "data/C202006211907A/data.h5"
parser.add_argument("syn_sac_dir"   )  # "mesh_hdf5/OUTPUT_FILES/sac"
parser.add_argument("out_adj_dir"   )  # "adj"
parser.add_argument("--cmt_in_ECEF", action="store_true")
parser.add_argument("--syn_is_grn", action="store_true")
parser.add_argument("--nproc", type=int, default=10)
parser.add_argument("--window_yaml", default=None)

args = parser.parse_args()

if args.nproc < 2:
    parser.error("nproc must >= 2!")

# with Misfit(args.misfit_h5file, 'w') as misfit:
if __name__ == "__main__":
    misfit = Misfit(args.misfit_h5file)
    print(f"======= [{datetime.now()}] read_config_file", flush=True)
    misfit.read_config_file(args.misfit_parfile)
    print(f"======= [{datetime.now()}] read_solver_parfile", flush=True)
    misfit.read_solver_parfile(args.solver_parfile)
    print(f"======= [{datetime.now()}] read_cmtsolution", flush=True)
    misfit.read_cmtsolution(args.cmt_file, ECEF=args.cmt_in_ECEF)
    print(f"======= [{datetime.now()}] read_channel_file", flush=True)
    misfit.read_channel_file(args.channel_file)
    print(f"======= [{datetime.now()}] read_data_h5", flush=True)
    misfit.read_data_h5(args.data_h5file)
    print(f"======= [{datetime.now()}] read_syn_sac", flush=True)
    misfit.read_syn_sac(args.syn_sac_dir, is_grn=args.syn_is_grn)
    print(f"======= [{datetime.now()}] setup_windows", flush=True)
    misfit.setup_windows(window_yaml=args.window_yaml)
    print(f"======= [{datetime.now()}] measure_adj", flush=True)
    # misfit.measure_adj()
    nproc = args.nproc - 1
    misfit.measure_adj_multiprocess(nproc=nproc)
    print(f"======= [{datetime.now()}] output_adj", flush=True)
    misfit.output_adj(args.out_adj_dir)

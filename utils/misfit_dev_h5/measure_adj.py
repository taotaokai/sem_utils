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

args = parser.parse_args()

with Misfit(args.misfit_h5file, 'w') as misfit:
    print(f"======= [{datetime.now()}] read_config_file")
    misfit.read_config_file(args.misfit_parfile)
    print(f"======= [{datetime.now()}] read_solver_parfile")
    misfit.read_solver_parfile(args.solver_parfile)
    print(f"======= [{datetime.now()}] read_cmtsolution")
    misfit.read_cmtsolution(args.cmt_file, ECEF=args.cmt_in_ECEF)
    print(f"======= [{datetime.now()}] read_channel_file")
    misfit.read_channel_file(args.channel_file)
    print(f"======= [{datetime.now()}] read_data_h5")
    misfit.read_data_h5(args.data_h5file)
    print(f"======= [{datetime.now()}] read_syn_sac")
    misfit.read_syn_sac(args.syn_sac_dir, is_grn=args.syn_is_grn)
    print(f"======= [{datetime.now()}] setup_windows")
    misfit.setup_windows()
    print(f"======= [{datetime.now()}] measure_adj")
    misfit.measure_adj()
    print(f"======= [{datetime.now()}] output_adj")
    misfit.output_adj(args.out_adj_dir)

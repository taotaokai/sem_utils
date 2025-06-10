import sys
import argparse
from datetime import datetime
from misfit import Misfit

parser = argparse.ArgumentParser()

parser.add_argument("misfit_h5file")  # "misfit.h5"
parser.add_argument("syn_sac_dir")  # "mesh_hdf5/OUTPUT_FILES/sac"
parser.add_argument("tag_diff", help="tag_diff")
parser.add_argument("--syn_is_grn", help="green\'s function", action="store_true")

args = parser.parse_args()

misfit = Misfit(args.misfit_h5file) 
misfit.read_syn_sac(args.syn_sac_dir, is_grn=args.syn_is_grn, is_diff=True, tag_diff=args.tag_diff)
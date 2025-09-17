#!/usr/bin/env python
# -*- coding: utf-8 -*-
import sys
import argparse

import numpy as np
from scipy.io import FortranFile
               
#====== user input
parser = argparse.ArgumentParser()

parser.add_argument("nproc", type=int, help="number of mesh slices")
parser.add_argument("model_dir", help="kernel dir")
parser.add_argument("out_dir", help="output dir")
parser.add_argument("--model_tag", default="cijkl_kernel", help="e.g. cijkl_kerenl")

args = parser.parse_args()
print(args)

nproc = args.nproc
model_dir = args.model_dir
model_tag = args.model_tag
out_dir = args.out_dir

#====== read in gll file
for iproc in range(nproc):

  print("# iproc", iproc)

  input_file = "%s/proc%06d_reg1_%s.bin"%(model_dir, iproc, model_tag)
  with FortranFile(input_file, 'r') as f:
    gll = f.read_reals(dtype='f4')
  gll = gll.reshape((-1, 21)) 

  C11 = gll[:, 0]
  C12 = gll[:, 1]
  C13 = gll[:, 2]
  C22 = gll[:, 6]
  C23 = gll[:, 7]
  C33 = gll[:, 11]
  C44 = gll[:, 15]
  C55 = gll[:, 18]
  C66 = gll[:, 20]

  # convert to Love parameters
  love_param = {
    'A': (3*C11 + 3*C22 + 2*C12 + 4*C66) / 8.0,
    'C': C33,
    'N': (C11 + C22 - 2*C12 + 4*C66) / 8.0,
    'L': (C44 + C55) / 2.0,
    'F': (C13 + C23) / 2.0,
  }

  for tag in love_param:
    output_file = "%s/proc%06d_reg1_%s_kernel.bin"%(out_dir, iproc, tag)
    with FortranFile(output_file, 'w') as f:
      f.write_record(np.array(love_param[tag], dtype='f4'))
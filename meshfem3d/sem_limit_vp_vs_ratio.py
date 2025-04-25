#!/usr/bin/env python
# -*- coding: utf-8 -*-
import sys

import numpy as np
from scipy.io import FortranFile
               
#====== user input
procnum_begin = int(sys.argv[1])
procnum_end = int(sys.argv[2])
vs_dir = str(sys.argv[3])
vs_tag = str(sys.argv[4])
vp_dir = str(sys.argv[5])
vp_tag = str(sys.argv[6])
min_vp_vs_ratio = float(sys.argv[7])
out_dir = str(sys.argv[8])

#====== read in gll file
for iproc in range(procnum_begin, procnum_end):

  print("# iproc", iproc)

  input_file = "%s/proc%06d_reg1_%s.bin"%(vs_dir, iproc, vs_tag)
  with FortranFile(input_file, 'r') as f:
    vs = f.read_reals(dtype='f4')

  input_file = "%s/proc%06d_reg1_%s.bin"%(vp_dir, iproc, vp_tag)
  with FortranFile(input_file, 'r') as f:
    vp = f.read_reals(dtype='f4')

  idx = vp/vs < min_vp_vs_ratio
  vp[idx] = min_vp_vs_ratio * vs[idx]

  output_file = "%s/proc%06d_reg1_%s.bin"%(out_dir, iproc, vp_tag)
  with FortranFile(output_file, 'w') as f:
    f.write_record(np.array(vp, dtype='f4'))

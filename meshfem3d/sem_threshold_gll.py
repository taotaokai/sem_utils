#!/usr/bin/env python
# -*- coding: utf-8 -*-
import sys

import numpy as np
from scipy.io import FortranFile
               
#====== user input
procnum_begin = int(sys.argv[1])
procnum_end = int(sys.argv[2])
model_dir = str(sys.argv[3])
model_ref_dir = str(sys.argv[4])
model_name = str(sys.argv[5])
out_dir = str(sys.argv[6])

#max_dgll = 1.5
max_dgll = 1.0
#min_dgll = -1

#====== read in gll file
for iproc in range(procnum_begin, procnum_end):

  print("# iproc", iproc)

  input_file = "%s/proc%06d_%s.bin"%(model_dir, iproc, model_name)
  with FortranFile(input_file, 'r') as f:
    gll = f.read_reals(dtype='f4')

  input_file = "%s/proc%06d_%s.bin"%(model_ref_dir, iproc, model_name)
  with FortranFile(input_file, 'r') as f:
    gll_ref = f.read_reals(dtype='f4')

  dgll = gll - gll_ref 

  idx = dgll>max_dgll
  gll[idx] = max_dgll + gll_ref[idx]
  print(np.sum(idx))

  #idx = dgll<min_dgll
  #gll[idx] = min_dgll + gll_ref[idx]
  #print(np.sum(idx))

  output_file = "%s/proc%06d_%s.bin"%(out_dir, iproc, model_name)
  with FortranFile(output_file, 'w') as f:
    f.write_record(np.array(gll, dtype='f4'))

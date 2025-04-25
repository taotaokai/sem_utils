#!/usr/bin/env python
# -*- coding: utf-8 -*-
import sys

import numpy as np
import numexpr as ne
from scipy.io import FortranFile
               
#====== user input
procnum_begin = int(sys.argv[1])
procnum_end = int(sys.argv[2])
model_dir1 = str(sys.argv[3])
model_dir2 = str(sys.argv[4])
model_name = str(sys.argv[5])
math_expr = str(sys.argv[6]) # e.g. "a - b"
out_dir = str(sys.argv[7])

#====== read in gll file
for iproc in range(procnum_begin, procnum_end):

  input_file = "%s/proc%06d_%s.bin"%(model_dir1, iproc, model_name)
  with FortranFile(input_file, 'r') as f:
    a = f.read_reals(dtype='f4')

  input_file = "%s/proc%06d_%s.bin"%(model_dir2, iproc, model_name)
  with FortranFile(input_file, 'r') as f:
    b = f.read_reals(dtype='f4')

  c = ne.evaluate(math_expr)
 
  print("# proc", iproc)
  print("gll1: min %10.1e max %10.1e"%(np.min(a), np.max(a)))
  print("gll2: min %10.1e max %10.1e"%(np.min(b), np.max(b)))
  print("dgll: min %10.1e max %10.1e"%(np.min(c), np.max(c)))

  output_file = "%s/proc%06d_%s.bin"%(out_dir, iproc, model_name)
  with FortranFile(output_file, 'w') as f:
    f.write_record(np.array(c, dtype='f4'))

#!/usr/bin/env python
# -*- coding: utf-8 -*-
import sys

import numpy as np
from scipy.io import FortranFile
               
#====== user input
procnum_begin = int(sys.argv[1])
procnum_end = int(sys.argv[2])
model_dir = str(sys.argv[3])
model_name = str(sys.argv[4])

#====== read in gll file
for iproc in range(procnum_begin, procnum_end):

  input_file = "%s/proc%06d_%s.bin"%(model_dir, iproc, model_name)
  
  with FortranFile(input_file, 'r') as f:
    gll = f.read_reals(dtype='f4')
  
  #print("n= %d"%(gll.size))
  print("iproc %06d min %10.1e max %10.1e"%(iproc, np.min(gll), np.max(gll)))

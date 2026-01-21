#!/usr/bin/env python
# -*- coding: utf-8 -*-
import sys

import numpy as np
from meshfem3d import read_gll_file
               
#====== user input
procnum_begin = int(sys.argv[1])
procnum_end = int(sys.argv[2])
model_dir = str(sys.argv[3])
model_name = str(sys.argv[4])

#====== read in gll file
for iproc in range(procnum_begin, procnum_end):

  gll = read_gll_file(model_dir, model_name, iproc)
  #print("n= %d"%(gll.size))
  print("iproc %06d min %10.1e max %10.1e"%(iproc, np.min(gll), np.max(gll)))
#!/usr/bin/env python
# -*- coding: utf-8 -*-
import sys

import numpy as np
# import numexpr as ne
from scipy.io import FortranFile

#====== user input
nproc = int(sys.argv[1])
model_dir = str(sys.argv[2])
model_tags = str(sys.argv[3])
math_expr = str(sys.argv[4]) # e.g. "rho * vsv**2"

#====== read in gll file
# model_dirs = model_dirs.split(",")
model_tags = model_tags.split(",")
ntags = len(model_tags)
# print(model_tags)

# assert(len(model_dirs) == len(model_tags))

for iproc in range(nproc):
    # print("# proc", iproc)
    for i in range(ntags):
        tag = model_tags[i]
        input_file = "%s/proc%06d_reg1_%s.bin"%(model_dir, iproc, tag)
        with FortranFile(input_file, 'r') as f:
            x = f.read_reals(dtype='f4')
            locals()[tag] = x
            # print("proc%03d %s min %.3f max %.3f"%(iproc, tag, np.min(x), np.max(x)))
            # a.append(x)

    # # c = ne.evaluate(math_expr)
    c = eval(math_expr)
    print("[proc%03d] %s: min %.3f max %.3f"%(iproc, math_expr, np.min(c), np.max(c)))

    # output_file = "%s/proc%06d_%s.bin"%(out_dir, iproc, out_tag)
    # with FortranFile(output_file, 'w') as f:
    #     f.write_record(np.array(c, dtype='f4'))

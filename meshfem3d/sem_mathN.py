#!/usr/bin/env python
# -*- coding: utf-8 -*-
import sys

import numpy as np
# import numexpr as ne
from scipy.io import FortranFile

#====== user input
procnum_begin = int(sys.argv[1])
procnum_end = int(sys.argv[2])
model_dir = str(sys.argv[3])
model_tags = str(sys.argv[4])
out_dir = str(sys.argv[5])
out_tag = str(sys.argv[6])
math_expr = str(sys.argv[7]) # e.g. "a[0] - a[1]"

#====== read in gll file
# model_dirs = model_dirs.split(",")
model_tags = model_tags.split(",")
ntags = len(model_tags)
print(model_tags)

# assert(len(model_dirs) == len(model_tags))

for iproc in range(procnum_begin, procnum_end):
    print("# proc", iproc)
    a = []
    for i in range(ntags):
        input_file = "%s/proc%06d_%s.bin"%(model_dir, iproc, model_tags[i])
        with FortranFile(input_file, 'r') as f:
            x = f.read_reals(dtype='f4')
            print("tag%d: min %10.1e max %10.1e"%(i, np.min(x), np.max(x)))
            a.append(x)

    # c = ne.evaluate(math_expr)
    c = eval(math_expr)
    print("out:  min %10.1e max %10.1e"%(np.min(c), np.max(c)))

    output_file = "%s/proc%06d_%s.bin"%(out_dir, iproc, out_tag)
    with FortranFile(output_file, 'w') as f:
        f.write_record(np.array(c, dtype='f4'))

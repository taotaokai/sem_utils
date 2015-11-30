#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Process misfit
"""
import sys
from misfit import Misfit
import json
import numpy as np

# read command line args
misfit_dir = str(sys.argv[1])

#------
print "\ninitialize\n"
misfit = Misfit()

#------
print "\nload data\n"
misfit.load(filename="%s/misfit.json" % (misfit_dir))

event_id = [ key for key in misfit.data['events'] ][0]

#------
print "\nrelocate\n"
window_id_list = ['F.p,P', 'F.s,S']
fix_depth = True
out_cmt_file = "%s/CMTSOLUTION.reloc" % (misfit_dir)
misfit.relocate_1d(event_id, 
        window_id_list=window_id_list,
        fix_depth=fix_depth,
        out_cmt_file=out_cmt_file)

print misfit.data['events'][event_id]['relocate']
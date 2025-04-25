#!/usr/bin/env python
# -*- coding: utf-8 -*-
""" create horizontal slice of SEM model at a given depth
"""
import sys

import numpy as np
from scipy.io import FortranFile

from meshfem3d_utils import sem_mesh_read

#====== 

mesh_file = str(sys.argv[1]) # <mesh_dir>/proc******_external_mesh.bin
out_file = str(sys.argv[2]) # proc*_idoubling.bin

mesh_data = sem_mesh_read(mesh_file)
idoubling = mesh_data['idoubling']
idoubling_ext = np.zeros(mesh_data['ibool'].shape) + idoubling

with FortranFile(out_file, 'w') as f:
  f.write_record(np.array(np.ravel(idoubling_ext, order='F'), dtype='f4'))

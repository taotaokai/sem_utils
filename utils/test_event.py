#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Process misfit
"""
import sys
from event import Event
import numpy as np

# read command line args
cmt_file = "DATA/CMTSOLUTION.init"

event = Event()

event.read_cmtsolution(cmt_file, isECEF=True)

event.make_cmtsolution("CMTSOLUTION.test", isECEF=True)
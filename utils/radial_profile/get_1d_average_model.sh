#!/bin/bash

mkdir 1d_model
ls radial_model/*.txt | awk -F"/" '{printf "./get_1d_average_model.py radial_model %s 1d_model \n",$2}' > tmp_sh
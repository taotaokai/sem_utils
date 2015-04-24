#!/bin/bash

databasedir=DATABASES_MPI
nproc=144
iregion=1

section=section
xyz=${section}.xyz

for modelname in vpv vph vsv vsh rho qmu eta
do

	echo interpolate $modelname ...
	./xinterp_mesh $databasedir $modelname $nproc $iregion < $xyz > model/${section}_${modelname}.xyz

done

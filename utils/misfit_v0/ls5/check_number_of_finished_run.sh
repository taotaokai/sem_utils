#!/bin/bash

job=${1:?[arg]need job name,e.g. syn,kernel,perturb}

ls */output_${job}/output_solver.txt | xargs grep "End of the simulation" | wc -l
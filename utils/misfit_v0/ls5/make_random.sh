#!/bin/bash

# randomly select events for linearized seismogram

event_list=${1:?[arg]need event_list}
ratio=${2:?[arg]need ratio e.g. 0.1}

awk -F"|" 'BEGIN{srand()}; NF&&$1!~/#/&&rand()<a' a=$ratio $event_list
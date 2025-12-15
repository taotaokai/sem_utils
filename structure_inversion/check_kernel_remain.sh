#!/bin/bash

echo "Working dir: $PWD"

control_file=${1:?[arg]need control_file}
source $control_file

list_file="event_within_mesh.lst"
output_file="event_within_mesh_remain.lst"

[ -f "$output_file" ] && rm "$output_file"

while read -r line; do
    [ -z "$line" ] && continue

    event_id=$(echo "$line" | awk '{print $1}')
    kernel_dir=$SEM_iter_dir/events/$event_id/output_kernel/kernel
    bin_num=$(ls $kernel_dir/*.bin | wc -l)

    if [ "$bin_num" -eq 24 ]; then
        echo "$event_id kernel and mask jobs done successfully."
    else
        echo "$event_id kernel and mask jobs failed......"
        echo "$line" >> "$output_file"
    fi

done < "$list_file"

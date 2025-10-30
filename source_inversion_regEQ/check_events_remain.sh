#!/bin/bash

cd ../../
echo "$PWD"
list_file="event_within_mesh.lst"
folder="./stage00.source/iter09/CMTSOLUTION_updated/"
output_file="event_within_mesh_remain.lst"

[ -f "$output_file" ] && rm "$output_file"

while read -r line; do

    [ -z "$line" ] && continue

    filename=$(echo "$line" | awk '{print $1}')

    if [ ! -e "$folder/$filename.cmt" ]; then
        echo "$line" >> "$output_file"
        echo "$filename does not exit in CMTSOLUTION_updated."
    fi
done < "$list_file"

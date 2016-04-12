#!/bin/bash

# use do-loop to repeate operations on sac files

sac_dir=${1:?[arg] need sac_dir}
sac_fileID=${2:?[arg] need sac_fileID}
sac_macro=${3:?[arg] need sac_macro}

tmpfile=$(mktemp -p .)
cat<<EOF > $tmpfile
do file wild dir "$sac_dir" "$sac_fileID"
r \$file
$sac_macro
enddo
EOF

sac<<EOF
macro $tmpfile
q
EOF

rm -rf $tmpfile

#END
#!/bin/bash

wkdir=$(pwd)

dest_dir=${wkdir##/home1/03244/ktao/}

echo $dest_dir

ssh jsg19 "mkdir -p ~/$dest_dir"
rsync -avzP --include="*/" --include="*/misfit/*.pdf" --exclude="*" -m * jsg19:~/$dest_dir
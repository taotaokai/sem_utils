#!/bin/bash

tags="Z.p,P|R.p,P|Z.s,S|R.s,S|T.s,S"

for tag in ${tags//|/ }
do
  echo ====== $tag
  pdftk *$tag.pdf cat output $tag.pdf
done
#!/bin/bash

cat - | awk -F'|' '$1!~/^#/&&$4~/..Z/{printf "%s %s.%s %s %s %s %s\n", $1,$2,$3,$5,$6,$7,$8}' | sort -u | column -t

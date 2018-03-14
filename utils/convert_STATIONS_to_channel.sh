#!/bin/bash

#AZ050 DEG10 6.408650 -7.307370 0.0 0.0

#BO|WTR||BHE|34.3739|136.5748|94|48|90.000000|0.000000|STS-2.5,TSM-1|||M/S|20|2001-04-07T00:00:00.0000|2015-02-25T23:59:59.9999

while read net staloc lat lon ele dep
do

  staloc=$(echo $staloc | sed "s/\./|/")

  echo "$net|$staloc|BHE|$lat|$lon|$ele|$dep|90.0|0.0|STS-2.5,TSM-1|||M/S|20|1900-01-01T00:00:00|2050-01-01T23:59:59"
  echo "$net|$staloc|BHN|$lat|$lon|$ele|$dep|0.0|0.0|STS-2.5,TSM-1|||M/S|20|1900-01-01T00:00:00|2050-01-01T23:59:59"
  echo "$net|$staloc|BHZ|$lat|$lon|$ele|$dep|0.0|-90.0|STS-2.5,TSM-1|||M/S|20|1900-01-01T00:00:00|2050-01-01T23:59:59"

done

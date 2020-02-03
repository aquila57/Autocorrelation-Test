#!/bin/bash
if [ -z $3 ]
then
echo "Autocorrelation test with GNU Scientific Library"
echo "random number generators"
echo "Usage: $0 index lag generator"
echo "Example: $0 3 5 taus2"
echo "Example: $0 3 5 mt19937"
echo "Example: $0 3 5 ranlxd2"
echo "Where index is i, the starting offset"
echo "Where lag is m, the distance between samples"
echo "For a list of GSL generators, type:"
echo "$0 3 5 ?"
exit 1
fi
GSL_RNG_TYPE="$3"
export GSL_RNG_TYPE
gslcorr $1 $2

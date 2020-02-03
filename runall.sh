#!/bin/bash
if [ -z $2 ]
then
echo "Autocorrelation test with GNU Scientific Library"
echo "random number generators"
echo "Usage: $0 index lag"
echo "Example: $0 3 5"
echo "Where index is i, the starting offset"
echo "Where lag is m, the distance between samples"
exit 1
fi
cat allgsl.lst | allgsl.sh $1 $2 >allgsl.out

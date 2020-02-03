#!/bin/bash
if [ -z $3 ]
then
echo "Usage: $0 index lag Dieharder_RNG_number"
echo "Example: $0 3 5 053"
echo "To get Dieharder_RNG_number,"
echo "dieharder -g -l"
exit 1
fi
dieharder -g $3 -t 1000000 -o | diecorr $1 $2 $3

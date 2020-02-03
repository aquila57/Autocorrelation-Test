#!/bin/bash
if [ -z $2 ]
then
echo "Usage: $0 index lag"
echo "Example: $0 3 5"
exit 1
fi
cat alldie.lst | alldie.sh $1 $2 >alldie.out

#!/bin/bash
while [ 1 ]
do
read
if [ $REPLY = "zzz" ]
then
exit 0
fi
tstgsl.sh $1 $2 $REPLY
done

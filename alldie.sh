#!/bin/bash
while [ 1 ]
do
read
if [ $REPLY = "zzz" ]
then
exit 0
fi
tstdie.sh $1 $2 $REPLY
done

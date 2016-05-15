#!/bin/sh

while read line
do
qsub -t $line $1
done


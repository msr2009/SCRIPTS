#!/bin/bash

STARTDATE=`/bin/date`
echo "**** JOB STARTED AT $STARTDATE"

python /net/fields/vol1/home/mattrich/LIBRARIES/SCRIPTS/read-vs-mutation_quality.py -f $1 -m $2 --read_columns 1,2 --mutation_columns 3,4 --maxmut 20 -n $3

ENDDATE=`/bin/date`
echo "**** JOB ENDED AT $ENDDATE"


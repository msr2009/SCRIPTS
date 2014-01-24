#!/bin/bash

STARTDATE=`/bin/date`
echo "**** JOB STARTED AT $STARTDATE"

python /net/fields/vol1/home/mattrich/LIBRARIES/SCRIPTS/count_histogram.py -f $1 -m $2 -c 1,2,3,4 --maxmut 20

ENDDATE=`/bin/date`
echo "**** JOB ENDED AT $ENDDATE"


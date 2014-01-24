#!/bin/bash

STARTDATE=`/bin/date`
echo "**** JOB STARTED AT $STARTDATE"

echo $1
echo "mapBEDcounts.py"

python /net/fields/vol1/home/mattrich/LIBRARIES/SCRIPTS/mapBEDcounts.py -b $1

echo "annotating"

python /net/fields/vol1/home/mattrich/LIBRARIES/SCRIPTS/annotate_orfsequences.py -v $1.counts -w $2 --pos-range 1,$3 --counts-input

ENDDATE=`/bin/date`
echo "**** JOB ENDED AT $ENDDATE"


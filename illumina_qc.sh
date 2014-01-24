#!/bin/bash

STARTDATE=`/bin/date`
echo "**** JOB STARTED AT $STARTDATE"

python /net/fields/vol1/home/mattrich/LIBRARIES/SCRIPTS/splitFASTQforCross_match.py -f $1.fq
cross_match 

ENDDATE=`/bin/date`
echo "**** JOB ENDED AT $ENDDATE"


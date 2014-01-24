#!/bin/bash

STARTDATE=`/bin/date`
echo "**** JOB STARTED AT $STARTDATE"

echo $1

echo "cross_match-ing"
cross_match -discrep_lists $1 $1.qual $2 > $1.crossmatch

ENDDATE=`/bin/date`
echo "**** JOB ENDED AT $ENDDATE"


#!/bin/bash

STARTDATE=`/bin/date`
echo "**** JOB STARTED AT $STARTDATE"

perl /net/fields/vol1/home/mattrich/LIBRARIES/SCRIPTS/SolexaQA.pl $@

ENDDATE=`/bin/date`
echo "**** JOB ENDED AT $ENDDATE"


#!/bin/bash

STARTDATE=`/bin/date`
echo "**** JOB STARTED AT $STARTDATE"

perl /net/fields/vol1/home/lstarita/bin/SolexaQA.pl $@

ENDDATE=`/bin/date`
echo "**** JOB ENDED AT $ENDDATE"


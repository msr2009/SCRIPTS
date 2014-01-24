#!/bin/bash

STARTDATE=`/bin/date`
echo "**** JOB STARTED AT $STARTDATE"

cd @
wc * > wordcount

ENDDATE=`/bin/date`
echo "**** JOB ENDED AT $ENDDATE"


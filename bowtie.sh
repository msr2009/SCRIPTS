#!/bin/bash

STARTDATE=`/bin/date`
echo "**** JOB STARTED AT $STARTDATE"

bowtie2 -t -p 8 --very-fast-local --no-discordant --no-mixed --no-unal --al $1.aligned -x pMR002 -U $1 -S $1.sam

ENDDATE=`/bin/date`
echo "**** JOB ENDED AT $ENDDATE"


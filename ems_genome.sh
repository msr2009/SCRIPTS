#!/bin/bash

STARTDATE=`/bin/date`
echo "**** JOB STARTED AT $STARTDATE"

bowtie2 -t -p 8 --very-fast-local --dovetail --no-discordant --no-mixed --no-unal 4TF -1 $1 -2 $2 -S $1.sam

ENDDATE=`/bin/date`
echo "**** JOB ENDED AT $ENDDATE"


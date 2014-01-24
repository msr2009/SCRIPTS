#!/bin/bash

STARTDATE=`/bin/date`
echo "**** JOB STARTED AT $STARTDATE"

export lane=$1; e=$2; cat s_${lane}_${e}_????_qseq.txt | perl /net/fields/vol1/home/mattrich/LIBRARIES/SCRIPTS/qseq2fastq33.pl > s_${lane}_${e}.fq

ENDDATE=`/bin/date`
echo "**** JOB ENDED AT $ENDDATE"
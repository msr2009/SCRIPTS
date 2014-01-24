#!/bin/bash

STARTDATE=`/bin/date`
echo "**** JOB STARTED AT $STARTDATE"

python /net/fields/vol1/home/mattrich/LIBRARIES/SCRIPTS/annotation_pipeline/annotate_genomes.py -m $1

ENDDATE=`/bin/date`
echo "**** JOB ENDED AT $ENDDATE"
#!/bin/bash

STARTDATE=`/bin/date`
echo "**** JOB STARTED AT $STARTDATE"

python /net/fields/vol1/home/mattrich/LIBRARIES/SCRIPTS/annotation_pipeline/batch_annotate_genomes.py -f $@

ENDDATE=`/bin/date`
echo "**** JOB ENDED AT $ENDDATE"


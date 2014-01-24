#!/bin/bash

STARTDATE=`/bin/date`
echo "**** JOB STARTED AT $STARTDATE"

python /net/fields/vol1/home/mattrich/LIBRARIES/SCRIPTS/enrich-0.2.1b/enrich/enrich.py --mode run_all --config_file $1

ENDDATE=`/bin/date`
echo "**** JOB ENDED AT $ENDDATE"


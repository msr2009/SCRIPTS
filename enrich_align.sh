#!/bin/bash

STARTDATE=`/bin/date`
echo "**** JOB STARTED AT $STARTDATE"

python /net/fields/vol1/home/mattrich/LIBRARIES/SCRIPTS/enrich-0.2.1b/enrich/enrich.py --mode read_align --config_file $@

ENDDATE=`/bin/date`
echo "**** JOB ENDED AT $ENDDATE"


#!/bin/bash

STARTDATE=`/bin/date`
echo "**** JOB STARTED AT $STARTDATE"

python ~/LIBRARIES/SCRIPTS/enrich-0.2.1b/enrich/enrich.py --mode run_all --config_file $@

ENDDATE=`/bin/date`
echo "**** JOB ENDED AT $ENDDATE"


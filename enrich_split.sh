#!/bin/bash

STARTDATE=`/bin/date`
echo "**** JOB STARTED AT $STARTDATE"

enrich --mode=fastq_filter --config_file=/net/fields/vol1/home/mattrich/sequencing/enrich_stuff/input/split_barcodes_8

ENDDATE=`/bin/date`
echo "**** JOB ENDED AT $ENDDATE"


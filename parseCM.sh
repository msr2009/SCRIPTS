#!/bin/bash

STARTDATE= `/bin/date`
echo "**** JOB STARTED AT $STARTDATE"

echo $1

echo "parsing cross_match output"
python /net/fields/vol1/home/mattrich/LIBRARIES/SCRIPTS/parseCross_matchOutput.py -c $1.crossmatch -f $2 --offset $3 --cutoff $4 > $1.crossmatch.mapped
#python /net/fields/vol1/home/mattrich/LIBRARIES/SCRIPTS/parseCross_matchOutput.py -c $1.crossmatch -f $2 --offset $3 --cutoff $4 --output_bed > $1.crossmatch.bed

ENDDATE= `/bin/date`
echo "**** JOB ENDED AT $ENDDATE"


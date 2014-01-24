#!/bin/bash

STARTDATE= `/bin/date`
echo "**** JOB STARTED AT $STARTDATE"

echo $1

echo "cross_match-ing"
cross_match -discrep_lists $1 $1.qual $2 > $1.crossmatch

echo "parsing cross_match output"
python /net/fields/vol1/home/mattrich/LIBRARIES/SCRIPTS/parseCross_matchOutput.py -c $1.crossmatch -f $2 --offset $3 --cutoff $4 > $1.crossmatch.mapped
#python /net/fields/vol1/home/mattrich/LIBRARIES/SCRIPTS/parseCross_matchOutput.py -c $1.crossmatch -f $2 --offset $3 --cutoff $4 --output_bed > $1.crossmatch.bed

echo "plotting BED with R"
LENGTH=`grep -v '>' $2 | python ~/LIBRARIES/SCRIPTS/getFileLength.py -f STDIN --lowercase`

R --vanilla '--args $1.crossmatch.mapped $LENGTH' < ~/LIBRARIES/SCRIPTS/plotBED.R

ENDDATE= `/bin/date`
echo "**** JOB ENDED AT $ENDDATE"


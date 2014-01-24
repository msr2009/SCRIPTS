STARTDATE=`/bin/date`
echo "**** JOB STARTED AT $STARTDATE"

python /net/fields/vol1/home/mattrich/LIBRARIES/SCRIPTS/mutation_histogram.py -f $@

ENDDATE=`/bin/date`
echo "**** JOB ENDED AT $ENDDATE"

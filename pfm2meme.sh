#!/bin/bash
echo $1

#first, python script to turn into transfac-ish file
#python ~/LIBRARIES/SCRIPTS/convert_psm_to_transfac.py $1 > tmp.tf 
cat $1 > tmp.tf

#next, we need to use sed to turn that output into a real transfac file
sed -e '/%TRANSFAC%/{r tmp.tf' -e 'd}' ~/LIBRARIES/transfac.skel > tmp.tf1
sed 's/%NAME%/'$1'/g' tmp.tf1 > $1.transfac

#then convert to meme
transfac2meme $1.transfac > $1.meme

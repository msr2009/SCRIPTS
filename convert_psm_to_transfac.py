"""
convert_pwm_to_transfac.py

A stupidly simple script to change the formatting of YeTFasCo's pwms to a Transfac-ish format
(the same as the output of pandas_pwm.py)

Matt Rich, 07/14
"""

import sys
import pandas

df = pandas.read_csv(sys.argv[1], header=None, index_col=0, sep="\t")
df.transpose().to_csv( sys.stdout, sep="\t", index_label="P0" )


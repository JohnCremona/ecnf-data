#!/bin/bash

field_type=$1
field=$2
prec=128
verbose=0
sage_out="sage.out.$1.$2"

echo field type: ${field_type}
echo field: ${field}
echo precision: ${prec}
echo verbosity: ${verbose}
echo output file: ${sage_out}
cline="from files import recompute_real_data as rrd; rrd(${field_type}, ${field}, prec=${prec}, verbose=${verbose})"
echo command line: ${cline}
echo ${cline} | sagedev -q > ${sage_out}

# Use with parallel as follows:
#
#
# parallel -j 5 --joblog runreal.log --progress ./runreal.sh cubics :::: cubics.txt
#
# 5 parallel jobs, data in cubics/, list of fields in cubics.txt
#

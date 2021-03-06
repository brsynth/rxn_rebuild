#!/bin/bash

if [ -z "$1" ]
  then
    echo "Please specify the input filename"
    exit
fi

datadir='tests/data/retrorules'
file=$1
log_file=log.$(basename $1)

while read -r line; do
    python -m rxn_rebuild --log_file $log_file "$(echo $line | cut -f1 -d ' ')" "$(echo $line | cut -f2 -d ' ')";
done < <(python $datadir/test.py $file)

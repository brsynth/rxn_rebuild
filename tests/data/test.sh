#!/bin/bash

datadir='tests/data'

while read -r line; do
    python -m rxn_rebuild "$(echo $line | cut -f1 -d ' ')" "$(echo $line | cut -f2 -d ' ')";
done < <(python $datadir/test.py $datadir/1-rp2_pathways.csv)
#!/bin/bash -ue
outfile="sAML1_vars.csv"
files=(1_vars_df.csv 2_vars_df.csv)
cat "${files[0]}" > $outfile
for f in "${files[@]:1}"; do
    tail -n +2 "$f" >> $outfile
done

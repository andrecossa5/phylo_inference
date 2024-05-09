#!/bin/bash -ue
outfile="sAML1_dataset.csv"
files=(1_dataset_df.csv 2_dataset_df.csv)
cat "${files[0]}" > $outfile
for f in "${files[@]:1}"; do
    tail -n +2 "$f" >> $outfile
done

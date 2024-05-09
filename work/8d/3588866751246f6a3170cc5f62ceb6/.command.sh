#!/bin/bash -ue
outfile="sAML1_job.csv"
files=(1_job_df.csv 2_job_df.csv)
cat "${files[0]}" > $outfile
for f in "${files[@]:1}"; do
    tail -n +2 "$f" >> $outfile
done

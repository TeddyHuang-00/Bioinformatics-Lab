#!/usr/bin/bash

source /appsnew/source/Python-3.8.6.sh

if [[ ! -d out/samples/${1} || ! -f out/samples/${1}/Aligned.sortedByCoord.out.bam ]]; then
    echo "Input file with ID ${1} not found in './out/samples' ! Please verify your path"
    exit
else
    echo "File exists! Going on..."
fi

threadsPerTask=5

pkurun-cnlong 1 $threadsPerTask \
    "htseq-count \
-m union \
-i gene_id \
-r pos \
-s no \
-n 5 \
-p bam \
-o out/samples/${1}/filtered.bam \
-f bam \
out/samples/${1}/Aligned.sortedByCoord.out.bam \
in/genes/chrX.gtf \
1>out/samples/${1}/htseq.log
"

#!/usr/bin/bash

if [[ ! -f in/samples/${1}_chrX_1.fastq.gz || ! -e in/samples/${1}_chrX_2.fastq.gz ]]; then
    echo "Input file with ID ${1} not found in './in/samples' ! Please verify your path"
    exit
else
    echo "File exists! Going on..."
fi

if [[ ! -d out/samples/${1} ]]; then
    echo "Creating temporary directory!"
    mkdir out/samples/${1}
fi

threadsPerTask=5

pkurun-cnlong 1 $threadsPerTask \
    "./STAR \
--runThreadN ${threadsPerTask} \
--genomeDir out/genome \
--readFilesIn in/samples/${1}_chrX_1.fastq.gz in/samples/${1}_chrX_2.fastq.gz \
--readFilesCommand zcat \
--twopassMode Basic \
--outSAMattrRGline ID:${1} \
--outSAMstrandField intronMotif \
--outSAMtype BAM SortedByCoordinate \
--outSAMunmapped Within \
--outFileNamePrefix out/samples/${1}/ \
&& \
samtools index out/samples/${1}/Aligned.sortedByCoord.out.bam
"

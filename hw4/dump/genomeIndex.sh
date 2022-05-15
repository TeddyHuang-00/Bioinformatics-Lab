#!/usr/bin/bash

threadsPerTask=20

pkurun-cnlong 1 $threadsPerTask \
    "./STAR \
--runThreadN ${threadsPerTask} \
--runMode genomeGenerate \
--genomeDir out/genome \
--genomeFastaFiles in/genome/chrX.fa \
--sjdbGTFfile in/genes/chrX.gtf \
--sjdbGTFchrPrefix chr \
--sjdbGTFtagExonParentTranscript Parent \
--sjdbOverhang 100 \
--genomeSAindexNbases 12
"

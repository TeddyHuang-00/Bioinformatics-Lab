#!/usr/bin/bash

threadsPerTask=20

pkurun-cnlong 1 $threadsPerTask \
    "gatk \
--java-options '-Xmx$((threadsPerTask * 3200))M -XX:+UseParallelGC -XX:ParallelGCThreads=${threadsPerTask}' \
GenotypeGVCFs \
-R sharedResources/GRCh38.d1.vd1.fa \
-V out/combined.g.vcf \
-O out/raw.vcf"

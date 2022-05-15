#!/usr/bin/bash

threadsPerTask=20

echo "" >chr.gvcf.list
for i in {{1..22},X,Y,M}; do
    echo "out/chr${i}.g.vcf" >>chr.gvcf.list
done
pkurun-cnlong 1 $threadsPerTask \
    "gatk \
--java-options '-Xmx$((threadsPerTask * 3200))M -XX:+UseParallelGC -XX:ParallelGCThreads=${threadsPerTask}' \
CombineGVCFs \
-R sharedResources/GRCh38.d1.vd1.fa \
--variant chr.gvcf.list \
-O out/combined.g.vcf && rm chr.gvcf.list"

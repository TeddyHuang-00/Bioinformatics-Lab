#!/usr/bin/bash

threadsPerTask=4

for i in {{1..22},X,Y,M}; do
    while ((($(sq | wc -l) - 1) * threadsPerTask >= 40 || ($(sq | wc -l) - 1) >= 20)); do
        sleep 5
    done
    pkurun-cnlong 1 $threadsPerTask \
        "gatk \
--java-options '-Xmx$((threadsPerTask * 3200))M -XX:+UseParallelGC -XX:ParallelGCThreads=${threadsPerTask}' \
HaplotypeCaller \
-R sharedResources/GRCh38.d1.vd1.fa \
-L chr${i} \
--dbsnp hg38/dbsnp_138.hg38.vcf.gz \
-I out/recalibrated_chr${i}.bam \
--minimum-mapping-quality 30 \
-ERC GVCF \
-O out/chr${i}.g.vcf"
done

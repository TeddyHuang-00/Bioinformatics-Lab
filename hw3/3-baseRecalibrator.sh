#!/usr/bin/bash

threadsPerTask=1

for i in {{1..22},X,Y,M}; do
    while ((($(sq | wc -l) - 1) * threadsPerTask >= 40 || ($(sq | wc -l) - 1) >= 20)); do
        sleep 5
    done
    pkurun-cnlong 1 $threadsPerTask \
        "gatk \
--java-options '-Xmx$((threadsPerTask * 3200))M -XX:+UseParallelGC -XX:ParallelGCThreads=${threadsPerTask}' \
BaseRecalibrator \
-I out/sorted_deduplicates_reads_chr${i}.bam \
-R sharedResources/GRCh38.d1.vd1.fa \
--known-sites hg38/dbsnp_138.hg38.vcf.gz \
--known-sites hg38/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz \
--known-sites hg38/1000G_phase1.snps.high_confidence.hg38.vcf.gz \
-O out/recalculated_data_chr${i}.table"
done

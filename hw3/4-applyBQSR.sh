#!/usr/bin/bash

threadsPerTask=10

for i in {{1..22},X,Y,M}; do
    while ((($(sq | wc -l) - 1) * threadsPerTask >= 40 || ($(sq | wc -l) - 1) >= 20)); do
        sleep 5
    done
    pkurun-cnlong 1 $threadsPerTask \
        "gatk \
--java-options '-Xmx$((threadsPerTask * 3200))M -XX:+UseParallelGC -XX:ParallelGCThreads=${threadsPerTask}' \
ApplyBQSR \
-I out/sorted_deduplicates_reads_chr${i}.bam \
--bqsr-recal-file out/recalculated_data_chr${i}.table \
-R sharedResources/GRCh38.d1.vd1.fa \
-O out/recalibrated_chr${i}.bam"
done

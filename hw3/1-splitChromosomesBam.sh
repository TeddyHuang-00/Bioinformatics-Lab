#!/usr/bin/bash

threadsPerTask=1

for i in {{1..22},X,Y,M}; do
    while ((($(sq | wc -l) - 1) * threadsPerTask >= 40 || ($(sq | wc -l) - 1) >= 20)); do
        sleep 5
    done
    pkurun-cnlong 1 $threadsPerTask \
    "samtools view -b in/9086c9b5f797c7084060137f4b818525_gdc_realn.bam chr${i} > out/chr${i}.bam"
done

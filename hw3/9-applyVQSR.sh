#!/usr/bin/bash

threadsPerTask=20

pkurun-cnlong 1 $threadsPerTask \
    "gatk \
--java-options '-Xmx$((threadsPerTask * 3200))M -XX:+UseParallelGC -XX:ParallelGCThreads=${threadsPerTask}' \
ApplyVQSR \
-R sharedResources/GRCh38.d1.vd1.fa \
-V out/raw.vcf \
-O out/variant_recalibrated.vcf \
-mode SNP \
--truth-sensitivity-filter-level 99.0 \
--tranches-file out/variant_calibrated.tranches \
--recal-file out/variant_calibrated.recal"

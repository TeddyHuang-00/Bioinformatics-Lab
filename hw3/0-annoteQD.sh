#!/usr/bin/bash
# Not used in pipeline

threadsPerTask=20

pkurun-cnlong 1 $threadsPerTask \
    "gatk \
--java-options '-Xmx$((threadsPerTask * 3200))M -XX:+UseParallelGC -XX:ParallelGCThreads=${threadsPerTask}' \
VariantAnnotator \
-V out/raw.vcf \
-O out/raw_annotations.vcf \
-A QualByDepth \
-A MappingQuality \
-A MappingQualityRankSumTest \
-A ReadPosRankSumTest \
-A FisherStrand \
-A StrandOddsRatio"

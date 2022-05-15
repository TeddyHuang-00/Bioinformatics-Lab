#!/usr/bin/bash

threadsPerTask=20

pkurun-cnlong 1 $threadsPerTask \
    "gatk \
--java-options '-Xmx$((threadsPerTask * 3200))M -XX:+UseParallelGC -XX:ParallelGCThreads=${threadsPerTask}' \
VariantRecalibrator \
-R sharedResources/GRCh38.d1.vd1.fa \
-V out/raw.vcf \
--resource:hapmap,known=false,training=true,truth=true,prior=15.0 hapmap_3.3.hg38.vcf.gz \
--resource:omni,known=false,training=true,truth=false,prior=12.0 hg38/1000G_omni2.5.hg38.vcf.gz \
--resource:1000G,known=false,training=true,truth=false,prior=10.0 hg38/1000G_phase1.snps.high_confidence.hg38.vcf.gz \
--resource:dbsnp,known=true,training=false,truth=false,prior=2.0 hg38/dbsnp_138.hg38.vcf.gz \
-an QD \
-an MQ \
-an MQRankSum \
-an ReadPosRankSum \
-an FS \
-an SOR \
-mode SNP \
-O out/variant_calibrated.recal \
--tranches-file out/variant_calibrated.tranches \
--rscript-file out/variant_calibrated.plots.R"

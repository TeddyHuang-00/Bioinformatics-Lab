#!/usr/bin/bash

threadsPerTask=20

pkurun-cnlong 1 $threadsPerTask \
    "gatk \
--java-options '-Xmx$((threadsPerTask * 3200))M -XX:+UseParallelGC -XX:ParallelGCThreads=${threadsPerTask}' \
VariantFiltration \
-V out/snp.vcf \
--filter-expression \"QD<2.0||MQ<40.0||FS>60.0||SOR>3.0||MQRankSum<-12.5||ReadPosRankSum<-8.0\" \
--filter-name \"PASS\" \
-O out/filtered.vcf"

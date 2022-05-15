#!/usr/bin/bash

# # # # # # # # # # # # # # # # # # # # # # #
# 以下为在执行流程前对文件路径和程序的检查  #
# # # # # # # # # # # # # # # # # # # # # # #

# 帮助信息
helpInfo="This script is written for start processing multiple RNA-seq fastq data.\n\
Paired end sequencing result must be in 'in/samples' directory,\n\
genome reference must be in 'in/genome' directory,\n\
and annotataion data must be in 'in/genes' directory.\n\
This script relies on samtools and STAR,\n\
and they must be included in the path for the script to work properly\n\
\n\
Usage:\n\
\n\
./automate.sh\n\
\n\
Yes, this will process all data files in a loop,\n\
there is no extra operation needed!\n\
\n\
by Nan Huang (huang_nan_2019@pku.edu.cn) April 27th, 2021
"

if [[ -n $1 ]]; then
    echo -e ${helpInfo}
    exit
fi

# 程序运行计时
start=$(date +%s)

# 检查对每个样本处理的脚本是否在同一路径下
if [[ ! -f pipeline.sh ]]; then
    echo -e ${helpInfo}
    echo "This automation script relies on pipeline.sh under the same directory, please verify the status of the file"
    exit
fi

if [[ ! -d ./in ]]; then
    echo -e ${helpInfo}
    echo "All required files must be in directory named 'in'"
    exit
fi

# 获取全部样本文件 ID，供循环使用
sampleList=$(ls in/samples/ | grep ERR | cut -d "_" -f 1 | uniq)
sampleList=($sampleList)

# 样本文件计数
echo "=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+="
sampleNumber=${#sampleList[@]}
echo "${sampleNumber} sample files detected for process"
for sampleID in ${sampleList[@]}; do
    echo ${sampleID}
done
echo "=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+=+="

# 激活 Python 环境用于跑 htseq-count
if ! hash htseq-count 2>/dev/null; then
    source /appsnew/source/Python-3.8.6.sh
fi

# 设定每个任务最大线程数
threadsPerTask=20

# 检查输出文件夹是否存在
if [[ ! -d out ]]; then
    echo "Creating output directory!"
    mkdir out
fi

# 检查基因组索引文件夹是否存在
if [[ ! -d out/genome ]]; then
    echo "Creating genome index directory!"
    mkdir out/genome
fi

# 检查比对结果文件夹是否存在
if [[ ! -d out/samples ]]; then
    echo "Creating mapped result directory!"
    mkdir out/samples
fi

# # # # # # # # # # # # #
# 以下为实际流水线操作  #
# # # # # # # # # # # # #

# STAR 构建基因组索引
echo "============================================="
echo "Start building index for the reference genome"
pkurun-cnlong 1 ${threadsPerTask} \
    "STAR \
--runThreadN ${threadsPerTask} \
--runMode genomeGenerate \
--genomeDir out/genome \
--genomeFastaFiles in/genome/chrX.fa \
--sjdbGTFfile in/genes/chrX.gtf \
--sjdbGTFtagExonParentTranscript Parent \
--sjdbOverhang 100 \
--genomeSAindexNbases 12 \
&>out/genome/star_index.log
"
echo "============================================="

# # 以下为挂载到计算节点上任务的详细注释信息
# STAR                                    # 调用 STAR 程序
# --runThreadN ${threadsPerTask}          # 指定最大线程数
# --runMode genomeGenerate                # 指明此步为构建基因组索引
# --genomeDir out/genome                  # 指定基因组索引文件输出路径
# --genomeFastaFiles in/genome/chrX.fa    # 指明基因组序列文件路径
# --sjdbGTFfile in/genes/chrX.gtf         # 指明基因组注释文件路径
# --sjdbGTFtagExonParentTranscript Parent # 指明输入文件为 gff 式
# --sjdbOverhang 100                      # 指定构建连接数据库时所用两侧序列长度
# --genomeSAindexNbases 12                # 指定合适的后缀索引长度 min(14, log2(GenomeLength)/2 - 1)
# &>out/genome/star_index.log             # 输出重定向至日志文件

# 等待索引构建完毕
while ((($(sq | wc -l) - 1) > 0)); do
    sleep 5
done

# 重设每个任务最大线程数
threadsPerTask=5

# 任务计数
processedCount=0

# 循环处理每个样本
for sampleID in ${sampleList[@]}; do
    # 如果最大线程数或任务数超出上限，则等待计算资源的释放
    while ((($(sq | wc -l) - 1) * threadsPerTask >= 40 || ($(sq | wc -l) - 1) >= 20)); do
        sleep 5
    done
    echo "---------------------------------------------"
    ./pipeline.sh -s ${sampleID} -t ${threadsPerTask} -q
    processedCount=$((processedCount + 1))
    echo "Job No.${processedCount}/${sampleNumber} in total"
    echo "---------------------------------------------"
done

# 等待所有样本处理完毕
echo "Waiting for the processing of all data"
while ((($(sq | wc -l) - 1) > 0)); do
    sleep 5
done

# 成功处理所有样本
end=$(date +%s)
runtime=$((end - start))
hours=$((runtime / 3600))
minutes=$(((runtime % 3600) / 60))
seconds=$(((runtime % 3600) % 60))
echo "All samples finished! Time elapsed: ${hours}:${minutes}:${seconds} (hh:mm:ss)"

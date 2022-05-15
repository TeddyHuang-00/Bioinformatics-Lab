#!/usr/bin/bash

# # # # # # # # # # # # # #
# 以下为对输入参数的处理  #
# # # # # # # # # # # # # #

# 帮助信息
helpInfo="This script is written for processing single-sample RNA-seq fastq data\n\
Index files must be under 'out/genome' directory,\n\
and paired end sequencing result must be in 'in/samples' directory.\n\
This script relies on htseq-count, samtools and STAR,\n\
and they must be included in the path for the script to work properly\n\
\n\
Usage:\n\
\n\
./pipeline.sh [-h] [-v | -q] [-t <threadsNum>] -s <sampleID>\n\
\n\
-h | --help\n\
\t\tPrint this message\n\
-s | --sample\n\
\t\t<sampleID>: Sample ID for processing\n\
-t | --threads\n\
\t\t<threadsNum>: Max threads that can be utilized\n\
-v | --verbose\n\
\t\tVerbose output\n\
-q | --quiet\n\
\t\tPrint minimum runtime message\n\
\n\
by Nan Huang (huang_nan_2019@pku.edu.cn) April 27th, 2021
"

verboseOutput=0

# 循环处理输入参数
while [[ -n $1 ]]; do
    case "$1" in
    -h | --help)
        # 输出帮助信息
        echo -e ${helpInfo}
        exit
        ;;
    -s | --sample)
        # 指定样本 ID
        sampleID="$2"
        echo "Pending job for sample ${sampleID}"
        shift
        ;;
    -t | --threads)
        # 设定每个任务最大线程数
        threadsPerTask="$2"
        echo "Set threads limit to ${threadsPerTask}"
        shift
        ;;
    -v | --verbose)
        # 指定输出冗余信息
        verboseOutput=0
        shift
        ;;
    -q | --quiet)
        # 指定隐藏冗余信息
        verboseOutput=1
        shift
        ;;
    *)
        # 遇到位置参数则抛出错误
        echo "Unrecognized parameter name $1"
        echo -e ${helpInfo}
        exit
        ;;
    esac
    shift
done

# # # # # # # # # # # # # # # # # # # # # # #
# 以下为在执行流程前对文件路径和程序的检查  #
# # # # # # # # # # # # # # # # # # # # # # #

# 检查是否指定样本名
if [[ -z $sampleID ]]; then
    echo -e ${helpInfo}
    echo "Please set sample ID for processing with '-sample'"
    echo "FATAL ERROR: sample ID not set!"
    exit
fi

# 检查是否设置线程数
if [[ -z $threadsPerTask ]]; then
    echo "Using default thread limit of 5"
    threadsPerTask=5
fi

# 检查环境中是否有 htseq-count 程序
if hash htseq-count 2>/dev/null; then
    if [[ $verboseOutput -eq 0 ]]; then
        echo "Python environment loaded successfully!"
    fi
else
    echo -e ${helpInfo}
    echo "Please activate a Python virtual environment containing htseq-count before running this script!"
    exit
fi

# 检查环境中是否有 samtools 程序
if hash samtools 2>/dev/null; then
    if [[ $verboseOutput -eq 0 ]]; then
        echo "Samtools loaded successfully!"
    fi
else
    echo -e ${helpInfo}
    echo "Please include samtools executable file into your path before running this script!"
    exit
fi

# 检查基因组索引文件是否存在
if [[ ! -f out/genome/Genome ]]; then
    echo -e ${helpInfo}
    echo "Need genme index to run this script, plese generate it manually or using automate.sh instead!"
    exit
fi

# 检查输出文件夹是否存在，如无则创建
if [[ ! -d out/samples/${sampleID} ]]; then
    if [[ $verboseOutput -eq 0 ]]; then
        echo "Creating temporary directory!"
    fi
    mkdir out/samples/${sampleID}
fi

# # # # # # # # # # # # # # # # # # # #
# 以下为实际对样本进行比对等流程操作  #
# # # # # # # # # # # # # # # # # # # #

# 使用 STAR 将 fastq 结果比对到基因组上
# 结果以排序后的 bam 文件形式呈现
pkurun-cnlong 1 $threadsPerTask \
    "STAR \
--runThreadN ${threadsPerTask} \
--genomeDir out/genome \
--readFilesIn in/samples/${sampleID}_chrX_1.fastq.gz in/samples/${sampleID}_chrX_2.fastq.gz \
--readFilesCommand zcat \
--twopassMode Basic \
--outSAMattrRGline ID:${sampleID} \
--outSAMstrandField intronMotif \
--outSAMtype BAM SortedByCoordinate \
--outSAMunmapped Within \
--outFileNamePrefix out/samples/${sampleID}/ \
&>out/samples/${sampleID}/star_mapping.log \
&& \
samtools index out/samples/${sampleID}/Aligned.sortedByCoord.out.bam \
&>out/samples/${sampleID}/samtools.log \
&& \
htseq-count \
-m union \
-i gene_id \
-r pos \
-s no \
-n ${threadsPerTask} \
-p bam \
-o out/samples/${sampleID}/filtered.bam \
-f bam \
out/samples/${sampleID}/Aligned.sortedByCoord.out.bam \
in/genes/chrX.gtf \
&>out/samples/${sampleID}/htseq.log
"

# # 以下为挂载到计算节点上任务的详细注释信息
# STAR                                                                 # 调用 STAR 程序
# --runThreadN ${threadsPerTask}                                       # 指定最大线程数
# --genomeDir out/genome                                               # 基因组索引等相关文件路径
# --readFilesIn in/samples/${sampleID}...                              # 双端测序一个样本的两个配对的文件
# --readFilesCommand zcat                                              # 使用 zcat 来读取压缩文件格式
# --twopassMode Basic                                                  # 使用 2pass 模式（每个样本独立跑两遍）
# --outSAMattrRGline ID:${sampleID}                                    # 在结果中增加额外的样本 ID 信息
# --outSAMstrandField intronMotif                                      # 在结果中指明正反链（cufflink 软件需要）
# --outSAMtype BAM SortedByCoordinate                                  # 以 bam 格式输出文件，并且按位置排序
# --outSAMunmapped Within                                              # 在结果中包含未比对上的读段
# --outFileNamePrefix out/samples/${sampleID}/                         # 指定输出路径
# &>out/samples/${sampleID}/star_mapping.log                           # 输出重定向至日志文件
# &&                                                                   # 成功执行后继续
# samtools index out/samples/${sampleID}/Aligned.sortedByCoord.out.bam # 为生成的 bam 文件建立索引（加快运算）
# &>out/samples/${sampleID}/samtools.log                               # 输出重定向至日志文件
# &&                                                                   # 成功执行后继续
# htseq-count                                                          # 调用 htseq-count 程序
# -m union                                                             # 指定以 union 方式处理重叠的比对
# -i gene_id                                                           # 指定 gtf 文件中被认作特征 ID 的字段
# -r pos                                                               # 指明输入文件按位置排序
# -s no                                                                # 指明测序结果非链特异
# -n ${threadsPerTask}                                                 # 指定最大进程数
# -p bam                                                               # 指定输出文件格式
# -o out/samples/${sampleID}/filtered.bam                              # 指定输出文件路径
# -f bam                                                               # 指明输入文件格式
# out/samples/${sampleID}/Aligned.sortedByCoord.out.bam                # 指明输入文件路径
# in/genes/chrX.gtf                                                    # 指明注释文件路径
# &>out/samples/${sampleID}/htseq.log                                  # 输出重定向至日志文件

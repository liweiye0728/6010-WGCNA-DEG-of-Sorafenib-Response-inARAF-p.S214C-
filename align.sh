#!/bin/bash
# 创建结果目录
mkdir -p reference/hisat2_index results/alignment

# 下载并解压人类参考基因组（如果还没有）
mkdir -p reference
wget https://ftp.ensembl.org/pub/release-113/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz -P reference/
wget https://ftp.ensembl.org/pub/release-113/gff3/homo_sapiens/Homo_sapiens.GRCh38.113.gff3.gz -P reference/
gunzip reference/*.gz

# 提取转录本信息用于HISAT2索引
hisat2_extract_splice_sites.py reference/Homo_sapiens.GRCh38.109.gtf > reference/splicesites.txt
hisat2_extract_exons.py reference/Homo_sapiens.GRCh38.109.gtf > reference/exons.txt

# 构建HISAT2索引
hisat2-build -p 8 \
    --ss reference/splicesites.txt \
    --exon reference/exons.txt \
    reference/Homo_sapiens.GRCh38.dna.primary_assembly.fa \
    reference/hisat2_index/genome

# 创建一个样本数组
samples=("ARAFmu_Ctl_1" "ARAFmu_Ctl_2" "ARAFmu_sor_1" "ARAFmu_sor_2" 
         "ARAFwt_Ctl_1" "ARAFwt_Ctl_2" "ARAFwt_sor_1" "ARAFwt_sor_2")

# 比对到参考基因组
for sample in "${samples[@]}"; do
    echo "Processing sample: $sample"
    
    # 使用HISAT2进行比对
    hisat2 -p 8 \
        --dta \
        -x reference/hisat2_index/genome \
        -1 results/trimmed/${sample}_1_paired.fq.gz \
        -2 results/trimmed/${sample}_2_paired.fq.gz \
        -S results/alignment/${sample}.sam \
        --summary-file results/alignment/${sample}_alignment_summary.txt
    
    # 将SAM转换为BAM
    samtools view -bS results/alignment/${sample}.sam > results/alignment/${sample}.bam
    
    # 排序BAM文件
    samtools sort -@ 8 -o results/alignment/${sample}_sorted.bam results/alignment/${sample}.bam
    
    # 索引排序后的BAM文件
    samtools index results/alignment/${sample}_sorted.bam
    
    # 删除中间文件节省空间
    rm results/alignment/${sample}.sam results/alignment/${sample}.bam
    
    # 生成比对统计信息
    samtools flagstat results/alignment/${sample}_sorted.bam > results/alignment/${sample}_flagstat.txt
done

# 使用StringTie进行转录本组装和定量
mkdir -p results/stringtie

# 下载和准备参考注释
if [ ! -f reference/Homo_sapiens.GRCh38.109.gtf ]; then
    echo "Reference GTF file not found!"
    exit 1
fi

# 对每个样本运行StringTie
for sample in "${samples[@]}"; do
    echo "Running StringTie for sample: $sample"
    
    # 转录本组装
    stringtie results/alignment/${sample}_sorted.bam \
        -G reference/Homo_sapiens.GRCh38.109.gtf \
        -o results/stringtie/${sample}.gtf \
        -p 8 \
        -A results/stringtie/${sample}_gene_abundances.tab \
        -B \
        -e
done

# 创建样本列表文件用于StringTie合并
echo -e "sample_id\tfile_path" > results/stringtie/sample_list.txt
for sample in "${samples[@]}"; do
    echo -e "${sample}\tresults/stringtie/${sample}.gtf" >> results/stringtie/sample_list.txt
done

# 使用StringTie合并所有样本的转录本
stringtie --merge \
    -G reference/Homo_sapiens.GRCh38.109.gtf \
    -o results/stringtie/merged.gtf \
    results/stringtie/sample_list.txt

# 使用合并后的注释重新定量每个样本
for sample in "${samples[@]}"; do
    echo "Re-estimating abundance for sample: $sample with merged annotation"
    
    stringtie results/alignment/${sample}_sorted.bam \
        -G results/stringtie/merged.gtf \
        -o results/stringtie/${sample}_merged.gtf \
        -p 8 \
        -A results/stringtie/${sample}_merged_gene_abundances.tab \
        -B \
        -e
done

# 运行prepDE.py脚本准备DESeq2输入文件
# 首先下载prepDE.py脚本(如果没有)
if [ ! -f prepDE.py ]; then
    wget https://ccb.jhu.edu/software/stringtie/dl/prepDE.py
    chmod +x prepDE.py
fi

# 准备样本信息文件
echo -e "sample_id\tfile_path" > results/stringtie/sample_paths.txt
for sample in "${samples[@]}"; do
    echo -e "${sample}\tresults/stringtie/${sample}_merged.gtf" >> results/stringtie/sample_paths.txt
done

# 运行prepDE.py生成计数矩阵
python prepDE.py \
    -i results/stringtie/sample_paths.txt \
    -g results/counts/gene_count_matrix.csv \
    -t results/counts/transcript_count_matrix.csv

echo "Alignment and quantification completed successfully!"

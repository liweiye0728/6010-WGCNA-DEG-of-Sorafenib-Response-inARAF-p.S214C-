#!/bin/bash
# 创建结果目录
# 在6010目录下运行
mkdir -p results/fastqc results/trimmed results/alignment results/counts results/deseq2 results/functional

# 质量检查 - FastQC
for sample in ARAFmu_Ctl_1 ARAFmu_Ctl_2 ARAFmu_sor_1 ARAF/ARAFmu_sor_2 ARAFwt_Ctl_1 ARAFwt_Ctl_2 ARAFwt_sor_1 ARAFwt_sor_2; do
    fastqc ARAF/${sample}/${sample}_1.fq.gz ARAF/${sample}/${sample}_2.fq.gz -o results/fastqc
done

# 汇总质量报告 - MultiQC
cd results/fastqc
multiqc .
cd ../..

# 数据修剪（如果需要） - Trimmomatic
for sample in ARAFmu_Ctl_1 ARAFmu_Ctl_2 ARAFmu_sor_1 ARAFmu_sor_2 ARAFwt_Ctl_1 ARAFwt_Ctl_2 ARAFwt_sor_1 ARAFwt_sor_2; do
    trimmomatic PE -threads 8 \
    ARAF/${sample}/${sample}_1.fq.gz ARAF/${sample}/${sample}_2.fq.gz \
    results/trimmed/${sample}_1_paired.fq.gz results/trimmed/${sample}_1_unpaired.fq.gz \
    results/trimmed/${sample}_2_paired.fq.gz results/trimmed/${sample}_2_unpaired.fq.gz \
    ILLUMINACLIP:adapters.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
done

#!/bin/bash

mkdir -p reference/hisat2_index results/alignment


mkdir -p reference
wget https://ftp.ensembl.org/pub/release-113/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz -P reference/
wget https://ftp.ensembl.org/pub/release-113/gff3/homo_sapiens/Homo_sapiens.GRCh38.113.gff3.gz -P reference/
gunzip reference/*.gz


hisat2_extract_splice_sites.py reference/Homo_sapiens.GRCh38.109.gtf > reference/splicesites.txt
hisat2_extract_exons.py reference/Homo_sapiens.GRCh38.109.gtf > reference/exons.txt


hisat2-build -p 8 \
    --ss reference/splicesites.txt \
    --exon reference/exons.txt \
    reference/Homo_sapiens.GRCh38.dna.primary_assembly.fa \
    reference/hisat2_index/genome


samples=("ARAFmu_Ctl_1" "ARAFmu_Ctl_2" "ARAFmu_sor_1" "ARAFmu_sor_2" 
         "ARAFwt_Ctl_1" "ARAFwt_Ctl_2" "ARAFwt_sor_1" "ARAFwt_sor_2")


for sample in "${samples[@]}"; do
    echo "Processing sample: $sample"
    

    hisat2 -p 8 \
        --dta \
        -x reference/hisat2_index/genome \
        -1 results/trimmed/${sample}_1_paired.fq.gz \
        -2 results/trimmed/${sample}_2_paired.fq.gz \
        -S results/alignment/${sample}.sam \
        --summary-file results/alignment/${sample}_alignment_summary.txt
    

    samtools view -bS results/alignment/${sample}.sam > results/alignment/${sample}.bam
    

    samtools sort -@ 8 -o results/alignment/${sample}_sorted.bam results/alignment/${sample}.bam
    

    samtools index results/alignment/${sample}_sorted.bam
    
 
    rm results/alignment/${sample}.sam results/alignment/${sample}.bam
    

    samtools flagstat results/alignment/${sample}_sorted.bam > results/alignment/${sample}_flagstat.txt
done


mkdir -p results/stringtie


if [ ! -f reference/Homo_sapiens.GRCh38.109.gtf ]; then
    echo "Reference GTF file not found!"
    exit 1
fi


for sample in "${samples[@]}"; do
    echo "Running StringTie for sample: $sample"
    

    stringtie results/alignment/${sample}_sorted.bam \
        -G reference/Homo_sapiens.GRCh38.109.gtf \
        -o results/stringtie/${sample}.gtf \
        -p 8 \
        -A results/stringtie/${sample}_gene_abundances.tab \
        -B \
        -e
done


echo -e "sample_id\tfile_path" > results/stringtie/sample_list.txt
for sample in "${samples[@]}"; do
    echo -e "${sample}\tresults/stringtie/${sample}.gtf" >> results/stringtie/sample_list.txt
done


stringtie --merge \
    -G reference/Homo_sapiens.GRCh38.109.gtf \
    -o results/stringtie/merged.gtf \
    results/stringtie/sample_list.txt


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


if [ ! -f prepDE.py ]; then
    wget https://ccb.jhu.edu/software/stringtie/dl/prepDE.py
    chmod +x prepDE.py
fi


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

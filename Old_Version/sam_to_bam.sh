#Prepare Your Data
#Reference Genome: The genome to which you want to align your reads.
#Sequence Reads: The raw sequencing data in FASTQ format.

#Index the Reference Genome

bwa index reference.fasta

#Align the reads

bwa mem reference.fasta reads1.fastq reads2.fastq > aligned_reads.sam

bwa mem UVE_YFV_1948_UG_MR896_TVP3236.fasta /mnt/Viro_Data/Mitarbeiter/Leyla/Yellow_Fever/UVE_YFV_1948_UG_MR896_TVP3236/trimmed/YFV-Uganda_S3_L001_R1_trimmed /mnt/Viro_Data/Mitarbeiter/Leyla/Yellow_Fever/UVE_YFV_1948_UG_MR896_TVP3236/trimmed/YFV-Uganda_S3_L001_R2_trimmed > gene_mapped_aligned.sam


bwa mem UVE_YFV_1948_UG_MR896_TVP3236.fasta /mnt/Viro_Data/Mitarbeiter/Leyla/Yellow_Fever/UVE_YFV_1948_UG_MR896_TVP3236/trimmed/YFV-Uganda_S3_L001_R1_trimm > UVE_YFV_1948_UG_MR896_TVP3236.sam



#Convert SAM to BA

samtools view -bS input.sam > output.bam


#Sort and Index BAM File

samtools sort input.bam -o sorted.bam
samtools index sorted.bam

#remove duplicates:pair end
samtools rmdup A_reads.bt2.sorted.bam A_reads.bt2.sorted.noDups.bam
#picard remove duplicate
 picard MarkDuplicates  \
      I=sorted_JX949181_YFV_strain_17D.bam \
      O=marked_duplicates_sorted_JX949181_YFV_strain_17D.bam \
      M=marked_dup_metrics.txt

#Alignment Summary

samtools flagstat sorted_mapped_reads.bam

#Generate Coverage Plot
samtools depth sorted_mapped_reads.bam > coverage_samtools.txt
bedtools genomecov -ibam sorted.bam -g genome.bed -d > coverage.txt

#Obtain Consensus Sequence

bcftools mpileup -f UVE_YFV_1948_UG_MR896_TVP3236.fasta sorted_bam_UVE_YFV_1948_UG_MR896_TVP3236.bam | bcftools call -mv --ploidy 1 -Ov -o variants.vcf
#second option
bcftools mpileup -f UVE_YFV_1948_UG_MR896_TVP3236.fasta sorted_bam_UVE_YFV_1948_UG_MR896_TVP3236.bam | bcftools call -mv --ploidy 1 -Ob -o variants.bcf
bcftools mpileup -Ou -f reference.fa alignments.bam | bcftools call -mv --ploidy 1  -Oz -o calls.vcf.gz

#plot
bcftools query -f '%MyAnnotation\n' calls.bcf | my-plotting-program

# Compress the VCF file using bgzip
bgzip variants.vcf

# Index the compressed VCF file using tabix
tabix -p vcf variants.vcf.gz

bcftools consensus -f virus_reference.fasta variants.vcf -o consensus.fasta

samtools depth sorted_mapped_reads.bam > coverage.txt
gnuplot -e "set terminal png; set output 'coverage_plot.png'; plot 'coverage.txt' using 1:3 with lines title 'Coverage'"

#Visualize

gnuplot -e "plot 'coverage.txt' with lines"



##Inspect Unmapped Reads
samtools view -f 4 input_sorted.bam > unmapped_reads.fastq

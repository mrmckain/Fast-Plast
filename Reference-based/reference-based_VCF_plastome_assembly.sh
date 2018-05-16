#!/bin/bash


########### Dependencies ###########
# java1.8.0
# perl5.22.1_threads
# gcc5.4.0
# samtools1.6
# bowtie2.2.5
# trimmomatic-0.36

#You can get by using whatever Java, Perl, and GCC libraries will work with your compilation of pilon, samtools, and bowtie2

##### How to run #####

# sh reference-based_VCF_plastome_assembly.sh Reference Target_Name Paired_End_File1 Paired_End_File2 Annotation_File_Reference Reference_Position_cutoff 

#Run trimmomatic to clean raw reads. Assumes trimmomatic is in the path. 
java -Xmx20G -classpath trimmomatic-0.36.jar  org.usadellab.trimmomatic.TrimmomaticPE -threads 1 $3 $4  $2.trimmed_P1.fq $2.trimmed_U1.fq  $2.trimmed_P2.fq  $2.trimmed_U2.fq ILLUMINACLIP:/mrm/bin/Fast-Plast/bin/Trimmomatic-0.36/adapters/TruSeq3-PE.fa:1:30:10 SLIDINGWINDOW:10:20 MINLEN:40

#Run bowtie2 build for database
bowtie2-build $1 $2

#Run bowtie2
bowtie2 --very-sensitive-local -x $2 -1 $2.trimmed_P1.fq -2 $2.trimmed_P2.fq -S $2.sam

#Convert SAM to BAM
samtools view -b $2.sam -o $2.bam
#rm $2.sam

#Sort BAM	
samtools sort -m 5G -o $2_sorted.bam $2.bam
#rm $2.bam

#Index BAM	
samtools index $2_sorted.bam $2_sorted.bai

#Run mpileup to get positional info
samtools mpileup -o $2.bcf -uf $1 $2_sorted.bam 

#Create VCF file with assumed haploidy
bcftools call -c --ploidy 1 -o $2_haploid.vcf $2.bcf

#Generate plastome sequence from VCF
perl plastid_genomes_from_vcf.pl X $5 $2_haploid.vcf $6 $2


#Clean sequence with pilon
bowtie2-build $2_cp_from_vcf.fasta $2

bowtie2 --very-sensitive-local -x $2 -1 $2.trimmed_P1.fq -2 $2.trimmed_P2.fq -S $2.sam

samtools view -b $2.sam -o $2.bam

samtools sort -m 5G -o $2_sorted.bam $2.bam

samtools index $2_sorted.bam $2_sorted.bai
	
java -Xmx10G -jar /mrm/bin/pilon-1.22.jar --genome $2_cp_from_vcf.fasta --bam $2\_sorted.bam --output $2\_pilon --outdir $2\_pilon --changes 


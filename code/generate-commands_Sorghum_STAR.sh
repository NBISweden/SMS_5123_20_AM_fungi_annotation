#!/usr/bash
###################################################################
#0. Resources                                                     #
###################################################################
PROJ_ID='snic2020-5-94'
PROJECT_DIR='/proj/uppstore2017083/private/SMS_5123_20_AM_fungi_annotation/'
CODE=$PROJECT_DIR'code/'
SAMPLES=$CODE'Samples_reads.txt'
INTERMEDIATE_DIR=$PROJECT_DIR'intermediate/'
RESULTS_DIR=$PROJECT_DIR'results/'
RAW_READS=$PROJECT_DIR'data/raw_data/'
TRIMMED_READS=$PROJECT_DIR'data/raw_data/Trimmed-reads/'
STAR=$INTERMEDIATE_DIR'STAR_Sorghum/'
STAR_QC=$INTERMEDIATE_DIR'STAR_Sorghum_QC/'
GENOME_DIR=$PROJECT_DIR'data/meta_data/reference/STARIndex_Sorghum/'
GFF_FILE=$PROJECT_DIR'data/meta_data/annotation/annotation_Sorghum.gff3'
GTF_FILE=$PROJECT_DIR'data/meta_data/annotation/annotation_Sorghum.gtf'
GENOME_SIZE=$PROJECT_DIR'data/meta_data/reference/genome_Sorghum.fai'
FEATURECOUNTS_DIR=$PROJECT_DIR'intermediate/featureCounts/STAR_Sorghum/'
READ_LENGTH=151

mkdir $TRIMMED_READS $STAR $STAR_QC $FEATURECOUNTS_DIR
THREADS=20

SBATCH="#!/bin/bash -l
#SBATCH -A $PROJ_ID
#SBATCH --mail-type=all
#SBATCH --mail-user=nimarafati@gmail.com"


###################################################################
# Align the reads to genome and Index bam files                   #
###################################################################
rm -f run_align.sh
THREADS=20
while read -r sample R1 R2
do
        echo "$SBATCH
#SBATCH -p node
#SBATCH -t 2:00:00
#SBATCH -J ${sample}_align_Sor

module load bioinfo-tools star samtools
cd \$SNIC_TMP
cp $TRIMMED_READS/*$sample*/${sample}_R1_paired.fq.gz \$SNIC_TMP &
cp $TRIMMED_READS/*$sample*/${sample}_R2_paired.fq.gz \$SNIC_TMP &
wait

#Alignign the reads
mkdir -p $STAR/$sample
STAR --genomeDir $GENOME_DIR \
--sjdbGTFfile $GTF_FILE \
--readFilesIn ${sample}_R1_paired.fq.gz ${sample}_R2_paired.fq.gz \
--runThreadN  $THREADS \
--twopassMode Basic \
--outWigType bedGraph \
--outSAMtype BAM SortedByCoordinate --limitBAMsortRAM 64324509440 \
--readFilesCommand zcat \
--runDirPerm All_RWX  \
--quantMode GeneCounts \
--outFilterScoreMinOverLread 0 \
--outFileNamePrefix $sample --outSAMattrRGline ID:$sample 'SM:${sample}' 

#Simplify name of the bam and count files 
mv ${sample}Aligned.sortedByCoord.out.bam ${sample}.sort.bam

#Indexing bam files
samtools index -@ $THREADS ${sample}.sort.bam
ln -s  ${sample}.sort.bam.bai ${sample}.sort.bai
mv ${sample}* $STAR/$sample/ 
cd $CODE
sbatch Sbatch_align_${sample}.script" >Sbatch_align_${sample}_Sorghum.script
echo "sbatch Sbatch_align_${sample}_Sorghum.script" >>run_align_Sorghum.sh
done<$SAMPLES

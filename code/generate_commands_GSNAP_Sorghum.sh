#!/usr/bash
###################################################################
#0. Resources                                                     #
###################################################################
PROJ_ID='snic2020-5-24'
PROJECT_DIR='/proj/uppstore2017083/private/SMS_5123_20_AM_fungi_annotation/'
CODE=$PROJECT_DIR'code/'
SAMPLES=$CODE'Samples_reads.txt'
INTERMEDIATE_DIR=$PROJECT_DIR'intermediate/'
RESULTS_DIR=$PROJECT_DIR'results/'
RAW_READS=$PROJECT_DIR'data/raw_data/'
TRIMMED_READS=$PROJECT_DIR'data/raw_data/Trimmed-reads/'


GSNAP=$INTERMEDIATE_DIR'GSNAP_Sorghum/'
GSNAP_QC=$INTERMEDIATE_DIR'GSNAP_QC_Sorghum/'

DISAMBIGUATE=$INTERMEDIATE_DIR'Disambiguate/'
STAR_SORGHUM=$INTERMEDIATE_DIR'STAR_Sorghum/'
GENOME_DIR=$PROJECT_DIR'data/meta_data/reference/STARIndex/'
GENOME_DIR_GSNAP=$PROJECT_DIR'data/meta_data/reference/GSNAP_Sorghum/'
GENOME='GSNAP_Sorghum'
GFF_FILE=$PROJECT_DIR'data/meta_data/annotation/annotation_Sorghum.gff'
GTF_FILE=$PROJECT_DIR'data/meta_data/annotation/annotation_Sorghum.gtf'
GENOME_SIZE=$PROJECT_DIR'data/meta_data/reference/genome_Sorghum.fa.fai'
FEATURECOUNTS_DIR=$PROJECT_DIR'intermediate/featureCounts/GSNAP_Sorghum/'
READ_LENGTH=141

mkdir $GSNAP $GSNAP_QC $FEATURECOUNTS_DIR 
THREADS=20

SBATCH="#!/bin/bash -l
#SBATCH -A $PROJ_ID
#SBATCH --mail-type=END,FAIL
#SBATCH -M Snowy
#SBATCH --mail-user=nimarafati@gmail.com"


###################################################################
# Align the reads to genome and Index bam files                   #
###################################################################

##GSNAP
THREADS=20
rm -f run_align_Sorghum.sh
while read -r sample R1 R2
do
        echo "$SBATCH
#SBATCH -p node
#SBATCH -t 8:00:00
#SBATCH -J ${sample}_align_Sorghum

module load bioinfo-tools gmap-gsnap/2017-09-11 samtools zlib

cd \$SNIC_TMP
cp $TRIMMED_READS$sample/${sample}_R1_paired.fq.gz \$SNIC_TMP &
cp $TRIMMED_READS$sample/${sample}_R2_paired.fq.gz \$SNIC_TMP &
wait

mkdir -p $GSNAP/$sample
gsnap --gunzip -m 0.1 -D $GENOME_DIR_GSNAP -d $GENOME --orientation FR -B 5 -N 1 -n 30 -E 4 --nthreads=$THREADS --gmap-mode=all -A sam -J 33 -O --quiet-if-excessive \
--read-group-id=$sample --read-group-name=$sample --read-group-library=200PE --read-group-platform=Illumina \
${sample}_R1_paired.fq.gz ${sample}_R2_paired.fq.gz | samtools view -bSt $GENOME_SIZE - >${sample}.bam

#sort by coordinate
samtools sort -@ $THREADS -O bam -o ${sample}.sort.bam ${sample}.bam
#Indexing bam files
samtools index -@ $THREADS ${sample}.sort.bam
ln -s  ${sample}.sort.bam.bai ${sample}.sort.bai
mv  ${sample}.sort* $GSNAP/$sample
cd $CODE
sbatch Sbatch_align_${sample}.script" >Sbatch_align_${sample}_Sorghum.script
echo "sbatch Sbatch_align_${sample}_Sorghum.script" >>run_align_Sorghum.sh
done<$SAMPLES

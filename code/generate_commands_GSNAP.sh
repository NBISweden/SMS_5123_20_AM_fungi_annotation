#!/usr/bash
###################################################################
#0. Resources                                                     #
###################################################################
PROJ_ID='snic2020-15-100'
PROJECT_DIR='/proj/uppstore2017083/private/SMS_5123_20_AM_fungi_annotation/'
CODE=$PROJECT_DIR'code/'
SAMPLES=$CODE'Samples_reads.txt'
INTERMEDIATE_DIR=$PROJECT_DIR'intermediate/'
RESULTS_DIR=$PROJECT_DIR'results/'
RAW_READS=$PROJECT_DIR'data/raw_data/'
TRIMMED_READS=$PROJECT_DIR'data/raw_data/Trimmed-reads/'

STAR=$INTERMEDIATE_DIR'STAR/'
STAR_QC=$INTERMEDIATE_DIR'STAR_QC/'

GSNAP=$INTERMEDIATE_DIR'GSNAP/'
GSNAP_SORGHUM=$INTERMEDIATE_DIR'GSNAP_Sorghum/'
GSNAP_QC=$INTERMEDIATE_DIR'GSNAP_QC/'

DISAMBIGUATE=$INTERMEDIATE_DIR'Disambiguate_GSNAP/'
XENOFILTER=$INTERMEDIATE_DIR'XenofilteR_GSNAP/'

STRINGTIE=$INTERMEDIATE_DIR'stringtie_GSNAP/'

GENOME_DIR_GSNAP=$PROJECT_DIR'data/meta_data/reference/GSNAP_AM_Fungi/'
GENOME='GSNAP_AM_Fungi'
GFF_FILE=$PROJECT_DIR'data/meta_data/annotation/annotation.gff'
GTF_FILE=$PROJECT_DIR'data/meta_data/annotation/annotation.gtf'
GENOME_SIZE=$PROJECT_DIR'data/meta_data/reference/genome.fai'

FEATURECOUNTS_DIR=$PROJECT_DIR'intermediate/featureCounts/GSNAP_Disambiguate/'

READ_LENGTH=141


mkdir $TRIMMED_READS $GSNAP $GSNAP_QC $FEATURECOUNTS_DIR $STRINGTIE
THREADS=20

SBATCH="#!/bin/bash -l
#SBATCH -A $PROJ_ID
#SBATCH --mail-type=END,FAIL
#SBATCH -M Snowy
#SBATCH --mail-user=nimarafati@gmail.com"

###################################################################
# FastQC                                                     	  #
###################################################################
rm -rf run_fastqc.txt
THREADS=2
while read -r sample R1 R2
do
	echo "$SBATCH
#SBATCH -J ${sample}_FastQC
#SBATCH -p core -n $THREADS
#SBATCH -t 2:00:00
module load bioinfo-tools FastQC

cd \$SNIC_TMP
cp $RAW_READS/$sample/*fastq.gz  \$SNIC_TMP/ 
mkdir  ${sample}_QC
mkdir -p $INTERMEDIATE_DIR/FastQC/
fastqc -t $THREADS -o ${sample}_QC -f fastq  *fastq.gz
mv ${sample}_QC $INTERMEDIATE_DIR/FastQC/" >Sbatch_fastqc_${sample}.script
echo "sbatch Sbatch_fastqc_${sample}.script" >>run_fastqc.txt
done<$SAMPLES
#Samples_reads.txt

###################################################################
# Trimming                                                        #
###################################################################
rm -f run_trim.log
THREADS=2
while read -r sample R1 R2
do
	echo "$SBATCH
#SBATCH -J ${sample}_trim
#SBATCH -p core -n $THREADS
#SBATCH -t 12:00:00
module load bioinfo-tools trimmomatic
mkdir $TRIMMED_READS/${sample}

cp $RAW_READS/$sample/*fastq.gz \$SNIC_TMP/
cd \$SNIC_TMP
cp /home/nimar/glob_old/Contamination/adapters.fa adapters.fa 
trimmomatic PE -threads $THREADS -phred33 -trimlog ${sample}.log \
$R1 $R2 \
${sample}_R1_paired.fq ${sample}_R1_unpaired.fq \
${sample}_R2_paired.fq ${sample}_R2_unpaired.fq \
ILLUMINACLIP:adapters.fa:2:40:15:8:true SLIDINGWINDOW:4:15 LEADING:10 TRAILING:10 MINLEN:36 HEADCROP:10
#Compress paired reads
gzip ${sample}_R1_paired.fq &
gzip ${sample}_R2_paired.fq &
wait 
#Transfer only paired reads
mv ${sample}_R1_paired.fq.gz $TRIMMED_READS/${sample} &
mv ${sample}_R2_paired.fq.gz $TRIMMED_READS/${sample} &
wait
cd $CODE
sbatch Sbatch_fastqc_trimmed_${sample}.script
sbatch Sbatch_align_${sample}.script" >Sbatch_trim_${sample}.script
	echo "sbatch Sbatch_trim_${sample}.script" >>run_trim.log
done<$SAMPLES

###################################################################
# FastQC after trimming                                           #
###################################################################
rm -rf run_fastqc_trimmed.txt
THREADS=2
while read -r sample R1 R2
do
        echo "$SBATCH
#SBATCH -J ${sample}_FastQC
#SBATCH -p core -n $THREADS
#SBATCH -t 2:00:00
module load bioinfo-tools FastQC

cd \$SNIC_TMP
cp $TRIMMED_READS/$sample/*fq.gz  \$SNIC_TMP/ 
mkdir  ${sample}_QC
mkdir -p $INTERMEDIATE_DIR/FastQC_Trimmed/
fastqc -t $THREADS -o ${sample}_QC -f fastq *fq.gz
mv ${sample}_QC $INTERMEDIATE_DIR/FastQC_Trimmed/" >Sbatch_fastqc_trimmed_${sample}.script
echo "sbatch Sbatch_fastqc_trimmed_${sample}.script" >>run_fastqc_trimmed.txt
done<$SAMPLES

###################################################################
# Align the reads to genome and Index bam files                   #
###################################################################

##GSNAP
THREADS=20
rm -f run_align.sh
while read -r sample R1 R2
do
        echo "$SBATCH
#SBATCH -p node
#SBATCH -t 5:00:00
#SBATCH -J ${sample}_align

module load bioinfo-tools gmap-gsnap/2017-09-11 samtools zlib

#Alignign the reads
cd \$SNIC_TMP
cp $TRIMMED_READS$sample/${sample}_R1_paired.fq.gz \$SNIC_TMP &
cp $TRIMMED_READS$sample/${sample}_R2_paired.fq.gz \$SNIC_TMP &
wait

mkdir -p $GSNAP/$sample
gsnap --gunzip -m 0.1 -D $GENOME_DIR_GSNAP -d $GENOME --orientation FR -B 5 -N 1 -n 30 -E 4 --nthreads=$THREADS --gmap-mode=all -A sam -J 33 -O --quiet-if-excessive \
--read-group-id=$sample --read-group-name=$sample --read-group-library=200PE --read-group-platform=Illumina \
${sample}_R1_paired.fq.gz ${sample}_R2_paired.fq.gz | samtools view -bSt $GENOME_SIZE  - >${sample}.bam

#sort by coordinate
samtools sort -@ $THREADS -O bam -o ${sample}.sort.bam ${sample}.bam
#Indexing bam files
samtools index -@ $THREADS ${sample}.sort.bam
ln -s  ${sample}.sort.bam.bai ${sample}.sort.bai
mv  ${sample}.sort* $GSNAP/$sample
cd $CODE
sbatch Sbatch_Disambiguate_${sample}.script" >Sbatch_align_${sample}.script
echo "sbatch Sbatch_align_${sample}.script" >>run_align.sh
done<$SAMPLES

###################################################################
# Disambiguate                                                    #
###################################################################
rm -f run_disambiguate.sh
while read -r sample R1 R2
do
	echo "$SBATCH
#SBATCH -p core -n 4
#SBATCH -t 1:00:00
#SBATCH -J Disambg_$sample
module load bioinfo-tools samtools

#cp $GSNAP/$sample/${sample}.sort.bam \$SNIC_TMP/${sample}_AM_Fungi.sort.bam &
#cp $GSNAP_SORGHUM/$sample/${sample}.sort.bam \$SNIC_TMP/${sample}_Sorghum.sort.bam &
#wait
#cd \$SNIC_TMP
#~/git/disambiguate_conda/bin/ngs_disambiguate -s $sample -a star \
#-o ${sample} \
#${sample}_AM_Fungi.sort.bam ${sample}_Sorghum.sort.bam 
#mv $sample/  $DISAMBIGUATE/
#cd $CODE/

samtools view -@2 -f8 $DISAMBIGUATE$sample/${sample}.disambiguatedSpeciesA.bam \
| sed 's/.*NM:i://' | cut -f1 | grep -v \"A006\" \
> $DISAMBIGUATE$sample/${sample}.disambiguatedSpeciesA_nM.txt &
samtools view -@2  -f8 $DISAMBIGUATE$sample/${sample}.disambiguatedSpeciesB.bam \
| sed 's/.*NM:i://' | cut -f1 | grep -v \"A006\" \
> $DISAMBIGUATE$sample/${sample}.disambiguatedSpeciesB_nM.txt &
wait

#sbatch Sbatch_featurecounts_${sample}.script
#sbatch Sbatch_QoRTs_${sample}.script 
#sbatch Sbatch_stringtie_${sample}.script" >Sbatch_Disambiguate_${sample}.script
echo "sbatch Sbatch_Disambiguate_${sample}.script" >>run_disambiguate.sh
done<$SAMPLES


###################################################################
# XenofiltR                                                       #
###################################################################
rm -f run_xenofiltR.sh
while read -r sample R1 R2
do
        echo "$SBATCH
#SBATCH -p core -n 5
#SBATCH -t 1:00:00
#SBATCH -J XenofiltR_$sample
module load bioinfo-tools samtools R_packages
mkdir -p $XENOFILTER/$sample/

cp $GSNAP/$sample/${sample}.sort.bam \$SNIC_TMP/${sample}_AM_Fungi.sort.bam &
cp $GSNAP_SORGHUM/$sample/${sample}.sort.bam \$SNIC_TMP/${sample}_Sorghum.sort.bam &
wait

cd \$SNIC_TMP

#Sorghum
echo \"library(XenofilteR)
library(Rsamtools)
library(GenomicAlignments)
library(BiocParallel)
library(futile.logger)

bp.param <- SnowParam(workers = 5, type = 'SOCK')

samples.list <-data.frame(Sorghum = '${sample}_Sorghum.sort.bam', 
  			 AM_Fungi = '${sample}_AM_Fungi.sort.bam')
output <- c('$sample')

XenofilteR(sample.list = samples.list, destination.folder = '$XENOFILTER/$sample/',  bp.param = bp.param, output.names = output)
\" >cmd_Sorghum.Rscript

R --vanilla -q <cmd_Sorghum.Rscript

mv $XENOFILTER$sample/Filtered_bams/ \
$XENOFILTER$sample/Filtered_bams_Sorghum/

#AM_Fungi
echo \"library(XenofilteR)
library(Rsamtools)
library(GenomicAlignments)
library(BiocParallel)
library(futile.logger)

bp.param <- SnowParam(workers = 5, type = 'SOCK')

samples.list <-data.frame(AM_Fungi = '${sample}_AM_Fungi.sort.bam',
			 Sorghum = '${sample}_Sorghum.sort.bam')
output <- c('$sample')

XenofilteR(sample.list = samples.list, destination.folder = '$XENOFILTER/$sample/',  bp.param = bp.param, output.names = output)
\" >cmd_AM_Fungi.Rscript

R --vanilla -q < cmd_AM_Fungi.Rscript 

mv $XENOFILTER$sample/Filtered_bams/ \
$XENOFILTER$sample/Filtered_bams_AM_Fungi/

samtools view -f8 $XENOFILTER$sample/Filtered_bams_AM_Fungi/${sample}_Filtered.bam \
| sed 's/.*NM:i://' | cut -f1 \
> $XENOFILTER$sample/Filtered_bams_AM_Fungi/${sample}_Filtered_nM.txt &

samtools view -f8 $XENOFILTER$sample/Filtered_bams_Sorghum/${sample}_Filtered.bam \
| sed 's/.*NM:i://' | cut -f1 \
> $XENOFILTER$sample/Filtered_bams_Sorghum/${sample}_Filtered_nM.txt &

samtools view -f8 $GSNAP/$sample/${sample}.sort.bam \
| sed 's/.*NM:i://' | cut -f1 | grep -v \"A006\" \
> $GSNAP/$sample/${sample}_nM.txt &

samtools view -f8 $GSNAP_SORGHUM/$sample/${sample}.sort.bam \
| sed 's/.*NM:i://' | cut -f1 | grep -v \"A006\" \
> $GSNAP_SORGHUM/$sample/${sample}_nM.txt &
wait
" >Sbatch_XenofiltR_${sample}.script
echo "sbatch Sbatch_XenofiltR_${sample}.script" >>run_xenofiltR.sh
done<$SAMPLES



###################################################################
# featurecounts                                                   #
###################################################################
rm -f run_featureCounts.sh
while read -r sample R1 R2
do
        echo "$SBATCH
#SBATCH -p node
#SBATCH -t 0:20:00
#SBATCH -J ${sample}_FC
cd \$SNIC_TMP
cp $DISAMBIGUATE/$sample/${sample}.disambiguatedSpeciesA.bam \$SNIC_TMP
featureCounts_path=\"$FEATURECOUNTS_DIR/${sample}\"
output=\"count-s-2\"
annotation=\"$GTF_FILE\"
mkdir -p \$featureCounts_path
~/glob/Software/subread-2.0.0-source/bin/featureCounts -s 2 -t exon -g gene_id -F GTF -C -T $THREADS -p -B \
-o \$output \
-a \$annotation \
${sample}.disambiguatedSpeciesA.bam
mv \$output* \$featureCounts_path" >Sbatch_featurecounts_${sample}.script
echo "sbatch Sbatch_featurecounts_${sample}.script " >>run_featureCounts.sh
done<$SAMPLES


###################################################################
# QC by QoRTs, qualimap                                           #
###################################################################
#Runnign QoRTs
rm -f run_QoRTs.sh
THREADS=20
while read -r sample R1 R2
do
        echo "$SBATCH
#SBATCH -p node
#SBATCH -t 5:00:00
#SBATCH -J QoRTS_${sample}
module load bioinfo-tools QualiMap QoRTs

cd \$SNIC_TMP
cp $DISAMBIGUATE/$sample/${sample}.disambiguatedSpeciesA.bam \$SNIC_TMP
QoRTs QC --prefilterImproperPairs --generatePlots --stranded --generatePlots --numThreads $THREADS --maxReadLength $READ_LENGTH  --minMAPQ 20 --title ${sample} \
--addFunctions JunctionCalcs,writeJunctionSeqCounts,writeSpliceExon,writeNovelSplices,annotatedSpliceExonCounts,makeAllBrowserTracks,FPKM,writeGeneCounts,writeDESeq,writeDEXSeq,writeGeneBody,writeGenewiseGeneBody,calcDetailedGeneCounts,makeJunctionBed,writeDocs,writeGeneBodyIv --chromSizes $GENOME_SIZE ${sample}.disambiguatedSpeciesA.bam \
$GTF_FILE $GSNAP_QC${sample}_QoRTs/ " >Sbatch_QoRTs_${sample}.script
echo "sbatch Sbatch_QoRTs_${sample}.script" >>run_QoRTs.sh
done<$SAMPLES


###################################################################
# StringTie                                                       #
###################################################################
THREADS=20
rm -f run_stringtie.sh
while read sample R1 R2
do
	echo "$SBATCH
#SBATCH -p node -n 1
#SBATCH -t 1:00:00
#SBATCH -J StringTie_$sample
module load bioinfo-tools StringTie
mkdir $STRINGTIE/$sample/
cd \$SNIC_TMP
cp $DISAMBIGUATE/$sample/${sample}.disambiguatedSpeciesA.bam \$SNIC_TMP

stringtie ${sample}.disambiguatedSpeciesA.bam \
--rf \
-o ${sample}_transcript.gtf \
-l $sample \
-p $THREADS \
-A ${sample}_abundance.txt 
mv ${sample}_*  $STRINGTIE/$sample/" >Sbatch_stringtie_${sample}.script
echo "sbatch Sbatch_stringtie_${sample}.script" >>run_stringtie.sh
done<$SAMPLES



###################################################################
# bbduk                                                           #
###################################################################

#There is n  rRNA annotated in the genome

#!/usr/bash
echo "#!/bin/bash -l
#SBATCH -A snic2020-15-24
#SBATCH --mail-type=END,FAIL
#SBATCH -M Snowy
#SBATCH --mail-user=nimarafati@gmail.com
#SBATCH -p core -n 5
#SBATCH -t 5:00:00
#SBATCH -J XenofiltR_Sample_TB-2463-31A
module load bioinfo-tools samtools R_packages
mkdir -p /proj/uppstore2017083/private/SMS_5123_20_AM_fungi_annotation/intermediate/XenofilteR_GSNAP_assembly1//Sample_TB-2463-31A/" >Sbatch_XenofiltR_Sample_TB-2463-31A_4_80.script
for i in {32..80..4}
do
	echo "
cp /proj/uppstore2017083/private/SMS_5123_20_AM_fungi_annotation/intermediate/GSNAP_assembly1//Sample_TB-2463-31A/Sample_TB-2463-31A.sort.bam \$SNIC_TMP/Sample_TB-2463-31A_AM_Fungi.sort.bam &
cp /proj/uppstore2017083/private/SMS_5123_20_AM_fungi_annotation/intermediate/GSNAP_Sorghum//Sample_TB-2463-31A/Sample_TB-2463-31A.sort.bam \$SNIC_TMP/Sample_TB-2463-31A_Sorghum.sort.bam &
wait

cd \$SNIC_TMP


#Sorghum
echo \"library(XenofilteR)
library(Rsamtools)
library(GenomicAlignments)
library(BiocParallel)
library(futile.logger)

bp.param <- SnowParam(workers = 5, type = 'SOCK')

samples.list <-data.frame(Sorghum = 'Sample_TB-2463-31A_Sorghum.sort.bam', 
	                  AM_Fungi = 'Sample_TB-2463-31A_AM_Fungi.sort.bam')
output <- c('Sample_TB-2463-31A')

XenofilteR(sample.list = samples.list, destination.folder = '/proj/uppstore2017083/private/SMS_5123_20_AM_fungi_annotation/intermediate/XenofilteR_GSNAP_assembly1//Sample_TB-2463-31A/',  bp.param = bp.param, output.names = output, MM_threshold = $i)
\" >cmd_Sorghum.Rscript

R --vanilla -q <cmd_Sorghum.Rscript

mv /proj/uppstore2017083/private/SMS_5123_20_AM_fungi_annotation/intermediate/XenofilteR_GSNAP_assembly1//Sample_TB-2463-31A/Filtered_bams/ /proj/uppstore2017083/private/SMS_5123_20_AM_fungi_annotation/intermediate/XenofilteR_GSNAP_assembly1/Sample_TB-2463-31A/Filtered_bams_Sorghum_$i/

#AM_Fungi
echo \"library(XenofilteR)
library(Rsamtools)
library(GenomicAlignments)
library(BiocParallel)
library(futile.logger)

bp.param <- SnowParam(workers = 5, type = 'SOCK')

samples.list <-data.frame(AM_Fungi = 'Sample_TB-2463-31A_AM_Fungi.sort.bam',
	                  Sorghum = 'Sample_TB-2463-31A_Sorghum.sort.bam')
output <- c('Sample_TB-2463-31A')

XenofilteR(sample.list = samples.list, destination.folder = '/proj/uppstore2017083/private/SMS_5123_20_AM_fungi_annotation/intermediate/XenofilteR_GSNAP_assembly1//Sample_TB-2463-31A/', bp.param = bp.param, output.names = output, MM_threshold = $i)
\" >cmd_AM_Fungi.Rscript

R --vanilla -q < cmd_AM_Fungi.Rscript 

mv /proj/uppstore2017083/private/SMS_5123_20_AM_fungi_annotation/intermediate/XenofilteR_GSNAP_assembly1//Sample_TB-2463-31A/Filtered_bams/ /proj/uppstore2017083/private/SMS_5123_20_AM_fungi_annotation/intermediate/XenofilteR_GSNAP_assembly1/Sample_TB-2463-31A/Filtered_bams_AM_Fungi_$i/

samtools view -f8 /proj/uppstore2017083/private/SMS_5123_20_AM_fungi_annotation/intermediate/XenofilteR_GSNAP_assembly1/Sample_TB-2463-31A/Filtered_bams_AM_Fungi_$i/Sample_TB-2463-31A_Filtered.bam | sed 's/.*NM:i://' | cut -f1 > /proj/uppstore2017083/private/SMS_5123_20_AM_fungi_annotation/intermediate/XenofilteR_GSNAP_assembly1/Sample_TB-2463-31A/Filtered_bams_AM_Fungi_$i/Sample_TB-2463-31A_Filtered_${i}_nM.txt &

samtools view -f8 /proj/uppstore2017083/private/SMS_5123_20_AM_fungi_annotation/intermediate/XenofilteR_GSNAP_assembly1/Sample_TB-2463-31A/Filtered_bams_Sorghum_$i/Sample_TB-2463-31A_Filtered.bam | sed 's/.*NM:i://' | cut -f1 > /proj/uppstore2017083/private/SMS_5123_20_AM_fungi_annotation/intermediate/XenofilteR_GSNAP_assembly1/Sample_TB-2463-31A/Filtered_bams_Sorghum_$i/Sample_TB-2463-31A_Filtered_${i}_nM.txt &

wait" >>Sbatch_XenofiltR_Sample_TB-2463-31A_4_80.script
done

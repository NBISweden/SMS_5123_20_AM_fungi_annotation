#!/usr/bash
echo "cd \$SNIC_TMP" >run_qorts_xenofilter.sh

for i in {12..28..4}
do
#	echo "#!/bin/bash -l
#SBATCH -A snic2020-15-24
#SBATCH --mail-type=END,FAIL
#SBATCH -M Snowy
#SBATCH --mail-user=nimarafati@gmail.com
#SBATCH -p node
#SBATCH -t 1:00:00
#SBATCH -J QoRTS_$i
#module load bioinfo-tools QualiMap QoRTs
echo "
cp /crex/proj/uppstore2017083/nobackup/private/SMS_5123_20_AM_fungi_annotation/scratch/XenofilteR_GSNAP_assembly1/Sample_TB-2463-31A/Filtered_bams_AM_Fungi_$i/Sample_TB-2463-31A_Filtered.bam \$SNIC_TMP
QoRTs QC --prefilterImproperPairs --generatePlots --stranded --generatePlots \
--minMAPQ 20 \
--numThreads 1 \
--maxReadLength 141  \
--title  Sample_31A_$i \
--addFunctions JunctionCalcs,writeJunctionSeqCounts,writeSpliceExon,writeNovelSplices,annotatedSpliceExonCounts,makeAllBrowserTracks,FPKM,writeGeneCounts,writeDESeq,writeDEXSeq,writeGeneBody,writeGenewiseGeneBody,calcDetailedGeneCounts,makeJunctionBed,writeDocs,writeGeneBodyIv \
--chromSizes /crex/proj/uppstore2017083/private/SMS_5123_20_AM_fungi_annotation/data/meta_data/reference_1/genome.fa.fai \
Sample_TB-2463-31A_Filtered.bam \
/crex/proj/uppstore2017083/private/SMS_5123_20_AM_fungi_annotation/data/meta_data/annotation_1/annotation_funannotate.gtf \
/crex/proj/uppstore2017083/nobackup/private/SMS_5123_20_AM_fungi_annotation/scratch/XenofilteR_GSNAP_assembly1/Sample_TB-2463-31A/Filtered_bams_AM_Fungi_${i}_QoRTs/ 

#cp /crex/proj/uppstore2017083/nobackup/private/SMS_5123_20_AM_fungi_annotation/scratch/XenofilteR_GSNAP_assembly1/Sample_TB-2463-31A/Filtered_bams_Sorghum_$i/Sample_TB-2463-31A_Filtered.bam \$SNIC_TMP
#QoRTs QC --prefilterImproperPairs --generatePlots --stranded --generatePlots \
#--minMAPQ 20 \
#--numThreads 1 \
#--maxReadLength 141  \
#--title  Sample_31A_$i \
#--addFunctions JunctionCalcs,writeJunctionSeqCounts,writeSpliceExon,writeNovelSplices,annotatedSpliceExonCounts,makeAllBrowserTracks,FPKM,writeGeneCounts,writeDESeq,writeDEXSeq,writeGeneBody,writeGenewiseGeneBody,calcDetailedGeneCounts,makeJunctionBed,writeDocs,writeGeneBodyIv \
#--chromSizes /crex/proj/uppstore2017083/private/SMS_5123_20_AM_fungi_annotation/data/meta_data/reference_3n/genome_Sorghum.fa.fai \
#Sample_TB-2463-31A_Filtered.bam \
#/crex/proj/uppstore2017083/private/SMS_5123_20_AM_fungi_annotation/data/meta_data/annotation_3n/annotation_Sorghum.gtf \
#/crex/proj/uppstore2017083/nobackup/private/SMS_5123_20_AM_fungi_annotation/scratch/XenofilteR_GSNAP_assembly1/Sample_TB-2463-31A/Filtered_bams_Sorghum_${i}_QoRTs/

" >>run_qorts_xenofilter.sh
done

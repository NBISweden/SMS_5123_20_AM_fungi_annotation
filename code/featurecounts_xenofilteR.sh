#!/usr/bash
for i in {4..28..4}
do
	echo "cd \$SNIC_TMP
cp /proj/uppstore2017083/private/SMS_5123_20_AM_fungi_annotation/intermediate/XenofilteR_GSNAP_assembly1/Sample_TB-2463-31A/Filtered_bams_AM_Fungi_$i/Sample_TB-2463-31A_Filtered.bam \$SNIC_TMP
featureCounts_path='/crex/proj/uppstore2017083/nobackup/private/SMS_5123_20_AM_fungi_annotation/scratch/XenofilteR_GSNAP_assembly1/Sample_TB-2463-31A_${i}_featurecounts/'
output='count-s-2'
annotation='/proj/uppstore2017083/private/SMS_5123_20_AM_fungi_annotation/data/meta_data/annotation_1/annotation_funannotate.gtf'
mkdir -p \$featureCounts_path
~/glob/Software/subread-2.0.0-source/bin/featureCounts -s 2 -t exon -g gene_id -Q 20 --primary  -F GTF -C -T 1 -p -B -o \$output -a \$annotation Sample_TB-2463-31A_Filtered.bam 
mv \$output* \$featureCounts_path" >>run_featurecounts_xenofilter.sh
done

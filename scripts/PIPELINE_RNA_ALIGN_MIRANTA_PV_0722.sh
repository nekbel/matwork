REFERENCE_GENOME="/data/refs/Mus_musculus/UCSC/mm10/Sequence/HISAT2Index/genome"


FILE="lane1RNA1MIRANTAWT1"
hisat2 -p24 \
-x $REFERENCE_GENOME \
-1 fastq/AAANCYHHV_MIRANTA_PV_0722_22s003034-1-1_Fett_lane1RNA1MIRANTAWT1_1_sequence.txt \
-2 fastq/AAANCYHHV_MIRANTA_PV_0722_22s003034-1-1_Fett_lane1RNA1MIRANTAWT1_2_sequence.txt \
-S sam/$FILE.sam \
--rna-strandness FR \
--score-min L,0,-0.5 \
-k 3
#OK

samtools sort -@ 32 -o sam-bam/$FILE.bam sam/$FILE.sam
rm sam/$FILE.sam
samtools index sam-bam/$FILE.bam

mv sam-bam/$FILE.bam /data/www/Lavigne/tracks/miranta_bams/
mv sam-bam/$FILE.bam.bai /data/www/Lavigne/tracks/miranta_bams/

echo 'track type=bam pairEndsByName=. name=lane1RNA1MIRANTAWT1_bam description=lane1RNA1MIRANTAWT1 visibility=pack bigDataUrl=http://139.91.171.49/Lavigne/tracks/miranta_bams/lane1RNA1MIRANTAWT1.bam'

FILE="lane1RNA2MIRANTAWT2"
hisat2 -p24 \
-x $REFERENCE_GENOME \
-1 fastq/AAANCYHHV_MIRANTA_PV_0722_22s003034-1-1_Fett_lane1RNA2MIRANTAWT2_1_sequence.txt \
-2 fastq/AAANCYHHV_MIRANTA_PV_0722_22s003034-1-1_Fett_lane1RNA2MIRANTAWT2_2_sequence.txt \
-S sam/$FILE.sam \
--rna-strandness FR \
--score-min L,0,-0.5 \
-k 3

samtools sort -@ 32 -o sam-bam/$FILE.bam sam/$FILE.sam
rm sam/$FILE.sam
samtools index sam-bam/$FILE.bam

mv sam-bam/$FILE.bam /data/www/Lavigne/tracks/miranta_bams/
mv sam-bam/$FILE.bam.bai /data/www/Lavigne/tracks/miranta_bams/

echo 'track type=bam pairEndsByName=. name=lane1RNA2MIRANTAWT2_bam description=lane1RNA2MIRANTAWT2 visibility=pack bigDataUrl=http://139.91.171.49/Lavigne/tracks/miranta_bams/lane1RNA2MIRANTAWT2.bam'



FILE="lane1RNA3MIRANTAWT3"
hisat2 -p24 \
-x $REFERENCE_GENOME \
-1 fastq/AAANCYHHV_MIRANTA_PV_0722_22s003034-1-1_Fett_lane1RNA3MIRANTAWT3_1_sequence.txt \
-2 fastq/AAANCYHHV_MIRANTA_PV_0722_22s003034-1-1_Fett_lane1RNA3MIRANTAWT3_2_sequence.txt \
-S sam/$FILE.sam \
--rna-strandness FR \
--score-min L,0,-0.5 \
-k 3

samtools sort -@ 32 -o sam-bam/$FILE.bam sam/$FILE.sam
        rm sam/$FILE.sam
        samtools index sam-bam/$FILE.bam

mv sam-bam/$FILE.bam /data/www/Lavigne/tracks/miranta_bams/
mv sam-bam/$FILE.bam.bai /data/www/Lavigne/tracks/miranta_bams/

echo 'track type=bam pairEndsByName=. name=lane1RNA3MIRANTAWT3_bam description=lane1RNA3MIRANTAWT3 visibility=pack bigDataUrl=http://139.91.171.49/Lavigne/tracks/miranta_bams/lane1RNA3MIRANTAWT3.bam'



samtools merge -@ 32 /data/www/Lavigne/tracks/miranta_bams/RNA_MIRANTA_WT.merged.bam /data/www/Lavigne/tracks/miranta_bams/lane1RNA*.bam
samtools sort -@ 32 -o /data/www/Lavigne/tracks/miranta_bams/RNA_MIRANTA_WT.merged.sorted.bam /data/www/Lavigne/tracks/miranta_bams/RNA_MIRANTA_WT.merged.bam
samtools index /data/www/Lavigne/tracks/miranta_bams/RNA_MIRANTA_WT.merged.sorted.bam


echo 'track type=bam pairEndsByName=. name=RNA_MIRANTA_WT.merged_bam description=RNA_MIRANTA_WT.merged visibility=pack bigDataUrl=http://139.91.171.49/Lavigne/tracks/miranta_bams/RNA_MIRANTA_WT.merged.sorted.bam'

# KO

FILE="lane1RNA4MIRANTAKO2"
hisat2 -p24 \
-x $REFERENCE_GENOME \
-1 fastq/AAANCYHHV_MIRANTA_PV_0722_22s003034-1-1_Fett_lane1RNA4MIRANTAKO2_1_sequence.txt \
-2 fastq/AAANCYHHV_MIRANTA_PV_0722_22s003034-1-1_Fett_lane1RNA4MIRANTAKO2_2_sequence.txt \
-S sam/$FILE.sam \
--rna-strandness FR \
--score-min L,0,-0.5 \
-k 3

samtools sort -@ 32 -o sam-bam/$FILE.bam sam/$FILE.sam
        rm sam/$FILE.sam
        samtools index sam-bam/$FILE.bam

mv sam-bam/$FILE.bam /data/www/Lavigne/tracks/miranta_bams/
mv sam-bam/$FILE.bam.bai /data/www/Lavigne/tracks/miranta_bams/

echo 'track type=bam pairEndsByName=. name=lane1RNA4MIRANTAKO2_bam description=lane1RNA4MIRANTAKO2 visibility=pack bigDataUrl=http://139.91.171.49/Lavigne/tracks/miranta_bams/lane1RNA4MIRANTAKO2.bam'


FILE="lane1RNA5MIRANTAKO3"
hisat2 -p24 \
-x $REFERENCE_GENOME \
-1 fastq/AAANCYHHV_MIRANTA_PV_0722_22s003034-1-1_Fett_lane1RNA5MIRANTAKO3_1_sequence.txt \
-2 fastq/AAANCYHHV_MIRANTA_PV_0722_22s003034-1-1_Fett_lane1RNA5MIRANTAKO3_2_sequence.txt \
-S sam/$FILE.sam \
--rna-strandness FR \
--score-min L,0,-0.5 \
-k 3

samtools sort -@ 32 -o sam-bam/$FILE.bam sam/$FILE.sam
        rm sam/$FILE.sam
        samtools index sam-bam/$FILE.bam

mv sam-bam/$FILE.bam /data/www/Lavigne/tracks/miranta_bams/
mv sam-bam/$FILE.bam.bai /data/www/Lavigne/tracks/miranta_bams/

echo 'track type=bam pairEndsByName=. name=lane1RNA5MIRANTAKO3_bam description=lane1RNA5MIRANTAKO3 visibility=squish bigDataUrl=http://139.91.171.49/Lavigne/tracks/miranta_bams/lane1RNA5MIRANTAKO3.bam'



FILE="lane1RNA6MIRANTAKO4"
hisat2 -p24 \
-x $REFERENCE_GENOME \
-1 fastq/AAANCYHHV_MIRANTA_PV_0722_22s003034-1-1_Fett_lane1RNA6MIRANTAKO4_1_sequence.txt \
-2 fastq/AAANCYHHV_MIRANTA_PV_0722_22s003034-1-1_Fett_lane1RNA6MIRANTAKO4_2_sequence.txt \
-S sam/$FILE.sam \
--rna-strandness FR \
--score-min L,0,-0.5 \
-k 3

samtools sort -@ 32 -o sam-bam/$FILE.bam sam/$FILE.sam
        rm sam/$FILE.sam
        samtools index sam-bam/$FILE.bam

mv sam-bam/$FILE.bam /data/www/Lavigne/tracks/miranta_bams/
mv sam-bam/$FILE.bam.bai /data/www/Lavigne/tracks/miranta_bams/

echo 'track type=bam pairEndsByName=. name=lane1RNA6MIRANTAKO4_bam description=lane1RNA6MIRANTAKO4 visibility=pack bigDataUrl=http://139.91.171.49/Lavigne/tracks/miranta_bams/lane1RNA6MIRANTAKO4.bam'



samtools merge -@ 32 /data/www/Lavigne/tracks/miranta_bams/RNA_MIRANTA_KO.merged.bam /data/www/Lavigne/tracks/miranta_bams/*RNA*KO*.bam
samtools sort -@ 32 -o /data/www/Lavigne/tracks/miranta_bams/RNA_MIRANTA_KO.merged.sorted.bam /data/www/Lavigne/tracks/miranta_bams/RNA_MIRANTA_KO.merged.bam
samtools index /data/www/Lavigne/tracks/miranta_bams/RNA_MIRANTA_KO.merged.sorted.bam

RNA_MIRANTA_WT.merged.sorted.bam 
RNA_MIRANTA_KO.merged.sorted.bam
#/data/www/Lavigne/tracks/miranta_bams/*RNA*KO*.bam

# WT merged
samtools merge RNA_MIRANTA_WT.sorted.merged.bam ./*WT*.bam
samtools index RNA_MIRANTA_WT.sorted.merged.bam

bamCoverage \
--binSize 10 \
--normalizeUsing RPKM \
--maxFragmentLength 0 \
-b RNA_MIRANTA_WT.sorted.merged.bam \
-o RNA_MIRANTA_WT.sorted.merged.bw \
-p 24


# KO merged
samtools merge RNA_MIRANTA_KO.sorted.merged.bam ./*KO*.bam
samtools index RNA_MIRANTA_KO.sorted.merged.bam

bamCoverage \
--binSize 10 \
--normalizeUsing RPKM \
--maxFragmentLength 0 \
-b RNA_MIRANTA_KO.sorted.merged.bam \
-o RNA_MIRANTA_KO.sorted.merged.bw \
-p 24

echo 'track type=bigWig name=RNA_MIRANTA_WT.merged_bw description=RNA_MIRANTA_WT.merged_bw visibility=full maxHeightPixels=64:64:11 color=50,50,150 bigDataUrl=http://139.91.171.49/Lavigne/tracks/miranta_bams/RNA_MIRANTA_WT.merged.sorted.bw'

samtools merge RNA_MIRANTA_WT.sorted.merged.bam ./*WT*.bam ; samtools index RNA_MIRANTA_WT.sorted.merged.bam ; samtools merge RNA_MIRANTA_KO.sorted.merged.bam ./*KO*.bam ; samtools index RNA_MIRANTA_KO.sorted.merged.bam ;

bamCoverage --binSize 10 --normalizeUsing RPKM --maxFragmentLength 0 -b RNA_MIRANTA_WT.sorted.merged.bam -o RNA_MIRANTA_WT.sorted.merged.bw -p 24
bamCoverage --binSize 10 --normalizeUsing RPKM --maxFragmentLength 0 -b RNA_MIRANTA_KO.sorted.merged.bam -o RNA_MIRANTA_KO.sorted.merged.bw -p 24




##############CONVERT TO bigWig###############

#bamCoverage -b /data/www/Lavigne/tracks/miranta_bams/$FILE.bam -o /data/www/Lavigne/tracks/miranta_bams/$FILE.bw
#mv sam/$FILE.bw /data/www/Lavigne/tracks/miranta_bams/
#echo 'track type=bigWig name=RNA2MIRANTAWT2_1_bw description=RNA2MIRANTAWT2_1 visibility=pack bigDataUrl=http://139.91.171.49/Lavigne/tracks/miranta_bams/AAANCYHHV_MIRANTA_PV_0722_22s003034-1-1_Fett_lane1RNA2MIRANTAWT2_1_sequence.txt.gz.bw'

#####Process lane1RNA1MIRANTAWT1
bamCoverage \
--binSize 10 \
--normalizeUsingRPKM \
--maxFragmentLength 0 \
-b /data/www/Lavigne/tracks/miranta_bams/lane1RNA1MIRANTAWT1.bam \
-o /data/www/Lavigne/tracks/miranta_bams/lane1RNA1MIRANTAWT1.bw \
-p 24
echo 'track type=bigWig name=lane1RNA1MIRANTAWT1_bw description=lane1RNA1MIRANTAWT1_fwd visibility=full bigDataUrl=http://139.91.171.49/Lavigne/tracks/miranta_bams/lane1RNA1MIRANTAWT1.bw'

bamCoverage \
--binSize 10 \
--normalizeUsingRPKM \
--maxFragmentLength 0 \
--filterRNAstrand=forward \
-b /data/www/Lavigne/tracks/miranta_bams/lane1RNA1MIRANTAWT1.bam \
-o /data/www/Lavigne/tracks/miranta_bams/lane1RNA1MIRANTAWT1_fwd.bw \
-p 24
echo 'track type=bigWig name=lane1RNA1MIRANTAWT1_fwd_bw description=lane1RNA1MIRANTAWT1_fwd visibility=full bigDataUrl=http://139.91.171.49/Lavigne/tracks/miranta_bams/lane1RNA1MIRANTAWT1_fwd.bw'

bamCoverage \
--binSize 10 \
--normalizeUsingRPKM \
--maxFragmentLength 0 \
--filterRNAstrand=reverse \
-b /data/www/Lavigne/tracks/miranta_bams/lane1RNA1MIRANTAWT1.bam \
-o /data/www/Lavigne/tracks/miranta_bams/lane1RNA1MIRANTAWT1_rev.bw \
-p 24
echo 'track type=bigWig name=lane1RNA1MIRANTAWT1_rev_bw description=lane1RNA1MIRANTAWT1_rev visibility=full bigDataUrl=http://139.91.171.49/Lavigne/tracks/miranta_bams/lane1RNA1MIRANTAWT1_rev.bw'


#####Process lane1RNA4MIRANTAKO2
bamCoverage \
--binSize 10 \
--normalizeUsingRPKM \
--maxFragmentLength 0 \
-b /data/www/Lavigne/tracks/miranta_bams/lane1RNA4MIRANTAKO2.bam \
-o /data/www/Lavigne/tracks/miranta_bams/lane1RNA4MIRANTAKO2.bw \
-p 24
echo 'track type=bigWig name=lane1RNA4MIRANTAKO2_bw description=lane1RNA4MIRANTAKO2_bw visibility=full bigDataUrl=http://139.91.171.49/Lavigne/tracks/miranta_bams/lane1RNA4MIRANTAKO2.bw'

bamCoverage \
--binSize 10 \
--normalizeUsingRPKM \
--maxFragmentLength 0 \
--filterRNAstrand=forward \
-b /data/www/Lavigne/tracks/miranta_bams/lane1RNA4MIRANTAKO2.bam \
-o /data/www/Lavigne/tracks/miranta_bams/lane1RNA4MIRANTAKO2_fwd.bw \
-p 24
echo 'track type=bigWig name=lane1RNA4MIRANTAKO2_fwd_bw description=lane1RNA4MIRANTAKO2_fwd visibility=dense bigDataUrl=http://139.91.171.49/Lavigne/tracks/miranta_bams/lane1RNA4MIRANTAKO2_fwd.bw'

bamCoverage \
--binSize 10 \
--normalizeUsingRPKM \
--maxFragmentLength 0 \
--filterRNAstrand=reverse \
-b /data/www/Lavigne/tracks/miranta_bams/lane1RNA4MIRANTAKO2.bam \
-o /data/www/Lavigne/tracks/miranta_bams/lane1RNA4MIRANTAKO2_rev.bw \
-p 24
echo 'track type=bigWig name=lane1RNA4MIRANTAKO2_rev_bw description=lane1RNA4MIRANTAKO2_rev visibility=dense bigDataUrl=http://139.91.171.49/Lavigne/tracks/miranta_bams/lane1RNA4MIRANTAKO2_rev.bw'




#####Process RNA_MIRANTA_KO_merged
bamCoverage \
--binSize 10 \
--normalizeUsingRPKM \
--maxFragmentLength 0 \
-b /data/www/Lavigne/tracks/miranta_bams/RNA_MIRANTA_KO.merged.sorted.bam \
-o /data/www/Lavigne/tracks/miranta_bams/RNA_MIRANTA_KO.merged.sorted.bw \
-p 24
echo 'track type=bigWig name=RNA_MIRANTA_KO.merged_bw description=RNA_MIRANTA_KO.merged_bw visibility=full maxHeightPixels=64:64:11 color=0,200,0 bigDataUrl=http://139.91.171.49/Lavigne/tracks/miranta_bams/RNA_MIRANTA_KO.merged.sorted.bw'

#####Process RNA_MIRANTA_WT_merged
bamCoverage \
--binSize 10 \
--normalizeUsingRPKM \
--maxFragmentLength 0 \
-b /data/www/Lavigne/tracks/miranta_bams/RNA_MIRANTA_WT.merged.sorted.bam \
-o /data/www/Lavigne/tracks/miranta_bams/RNA_MIRANTA_WT.merged.sorted.bw \
-p 24
echo 'track type=bigWig 
name=RNA_MIRANTA_WT.merged_bw 
description=RNA_MIRANTA_WT.merged_bw 
visibility=full 
maxHeightPixels=64:64:11 
color=50,50,150 
bigDataUrl=http://139.91.171.49/Lavigne/tracks/miranta_bams/RNA_MIRANTA_WT.merged.sorted.bw'




track type=bigWig name=RNA_MIRANTA_KO.sorted.merged description=RNA_MIRANTA_KO.sorted.merged visibility=full bigDataUrl=http://139.91.171.49/Lavigne/Nektarios/tracks/miranta_PV/RNA_MIRANTA_KO.sorted.merged.bw color=0,255,0
track type=bigWig name=RNA_MIRANTA_WT.sorted.merged description=RNA_MIRANTA_WT.sorted.merged visibility=full bigDataUrl=http://139.91.171.49/Lavigne/Nektarios/tracks/miranta_PV/RNA_MIRANTA_WT.sorted.merged.bw color=255,0,0


hub MirantaHub
shortLabel Miranta Hub
longLabel RNA and ATAC runs on the Miranta files
useOneFile on
email nek.belmezos@gmail.com

genome mm10

track ComparisonKOWT
type bigWig
container multiWig
visibility full
aggregate transparentOverlay
showSubtrackColorOnUi on
maxHeightPixels 500:100:8

        track KO
        type bigWig
        bigDataUrl http://139.91.171.49/Lavigne/Nektarios/tracks/miranta_PV/RNA_MIRANTA_KO.sorted.merged.bw
        parent ComparisonKOWT
        color 70,250,150

        track WT
        type bigWig
        bigDataUrl http://139.91.171.49/Lavigne/Nektarios/tracks/miranta_PV/RNA_MIRANTA_WT.sorted.merged.bw
        parent ComparisonKOWT
        color 207,21,62
















# Hub track example
track multiWig1
type bigWig
container multiWig
shortLabel Ex. multiWig container
longLabel This multiWig overlay track graphs points from three bigWig files.
visibility full
aggregate transparentOverlay
showSubtrackColorOnUi on
maxHeightPixels 500:100:8
viewLimits 1:20
priority 1
html examplePage

                track wig1
                bigDataUrl http://hgdownload.soe.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeCshlLongRnaSeq/wgEncodeCshlLongRnaSeqA549CellLongnonpolyaMinusRawSigRep1.bigWig
                shortLabel Overlay bigWig1
                longLabel This is an example bigWig1 displaying Raw Signal from the ENCODE RNA-seq CSHL track, graphing just points as default.
                parent multiWig1
                graphTypeDefault points
                type bigWig
                color 255,0,0

                track wig2
                bigDataUrl http://hgdownload.soe.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeCshlLongRnaSeq/wgEncodeCshlLongRnaSeqA549CellLongnonpolyaPlusRawSigRep1.bigWig
                shortLabel Overlay bigWig2
                longLabel This is an example bigWig2 displaying Raw Signal from the ENCODE RNA-seq CSHL track, graphing just points as default.
                graphTypeDefault points
                parent multiWig1
                type bigWig
                color 0,255,0

                track wig3
                bigDataUrl http://hgdownload.soe.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeCshlLongRnaSeq/wgEncodeCshlLongRnaSeqAg04450CellLongnonpolyaPlusRawSigRep1.bigWig
                shortLabel Overlay bigWig3
                longLabel This is an example bigWig3 displaying Raw Signal from the ENCODE RNA-seq CSHL track, graphing just points as default. 
                graphTypeDefault points
                parent multiWig1
                type bigWig
                color 95,158,160



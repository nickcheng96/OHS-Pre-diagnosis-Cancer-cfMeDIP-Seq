#!/bin/bash
source /oicr/local/Modules/default/init/bash 
#module loadings
module load samtools/1.9
module load fastqc/0.11.8
module load bowtie2/2.3.5.1
module load picard/2.21.2
module load python/3.5
module load java/8
module load bedtools/2.27

#loading other tools
bowtie2_index="/.mounts/labs/awadallalab/private/ncheng/references/hg38/Bowtie2"
BEDOPS="/.mounts/labs/awadallalab/private/ncheng/softwares/BEDOPS/v2.4.38/bin"
trimmomatic_dir="/.mounts/labs/awadallalab/private/ncheng/softwares/trimmomatic/v0.33/trimmomatic/build/bin"
picard_dir="/.mounts/labs/awadallalab/private/ncheng/softwares/picard-tools-1.119"

#setting paths to reference files
bt2ref="hg38_F19K16_F24B22.fa"
fastq_dir="/path/to/fastq/"

#setting other directories
script_dir="/path/to/scripts"
medips_script=$script_dir/medips.R

#setting output directory
outdir="/.mounts/labs/awadallalab/private/ncheng/cfmedip_data/cptp_samples/AIX3/novaseq/novaseq_umitools_output_trimq.nontopup.hg38/"
mkdir -p $outdir

#specifying files
runfile="fileid" #replace fileid with file name
r1file=${runfile}_R1.fastq.gz
r2file=${runfile}_R2.fastq.gz
r1=$fastq_dir/$r1file
r2=$fastq_dir/$r2file
tmp=paired_trim_$r1file 
align_output=${tmp%_R1.fastq.gz}

#start#
mkdir -p $outdir
cd $outdir

#step 0 - creating directories
mkdir -p $outdir/1_trim
mkdir -p $outdir/2_fastqc
mkdir -p $outdir/3_align_qc
mkdir -p $outdir/4_preprocess
mkdir -p $outdir/5_picard
mkdir -p $outdir/6_medips_qc/window_300
mkdir -p $outdir/7_thalia_qc
mkdir -p $outdir/8_consensuscruncher
mkdir -p $outdir/8_qc_metrics


#step 1 - trimming adapters
ptr1="$outdir/1_trim/paired_trim_$r1file"
ptr2="$outdir/1_trim/paired_trim_$r2file"

#step 2 - fastqc report
fastqc $r1 -o 2_fastqc
fastqc $r2 -o 2_fastqc

#step 3 - extracting UMIs and appending to header
umi_tools extract --extract-method=regex --umi-separator="__"  --bc-pattern="(?P<umi_1>^[ACGTN]{3}[ACGN])(?P<discard_1>T)|(?P<umi_1>^[ACGTN]{3})(?P<discard_1>T)" \
                      --bc-pattern2="(?P<umi_2>^[ACGTN]{3}[ACGN])(?P<discard_2>T)|(?P<umi_2>^[ACGTN]{3})(?P<discard_2>T)" \
                       --stdin=$r1 --read2-in=$r2 --stdout=${r1}_extract --read2-out=${r2}_extract --log=${r1}.umitools.log.txt

#trimming adapters
java -jar $trimmomatic_dir/trimmomatic.jar PE \
        ${r1}_extract ${r2}_extract \
        $outdir/1_trim/paired_trim_$r1file $outdir/1_trim/unpaired_trim_$r1file $outdir/1_trim/paired_trim_$r2file $outdir/1_trim/unpaired_trim_$r2file TRAILING:3 MINLEN:30 SLIDINGWINDOW:4:7

#removing redundant/intermediate files
rm $outdir/1_trim/unpaired_trim_$r1file
rm $outdir/1_trim/unpaired_trim_$r2file
rm ${r1}_extract
rm ${r2}_extract

#aligning reads to hg38 + A. thaliana genome
bowtie2 -p 8 -x $bowtie2_index/hg38_F19K16_F24B22 \
        -1 $ptr1 \
        -2 $ptr2 \
        -S 3_align_qc/${align_output}.sam
        
        
rm $ptr1
rm $ptr2

#converting to bam file
samtools view -b -h 3_align_qc/${align_output}.sam  | samtools sort -o 4_preprocess/${align_output}.sorted.bam
rm 3_align_qc/$align_output.sam

#collapsing PCR duplicate reads according to appended header
#UMIs have variable length (6-8 base pairs) so BAM files are partitioned according to UMI lengths
samtools view 4_preprocess/${align_output}.sorted.bam >   4_preprocess/${align_output}.tmp.sam
#umis with length 6
samtools view -H  4_preprocess/${align_output}.sorted.bam >  4_preprocess/${align_output}.6.tmp.sam 
grep -P "^.*__.{6}\t"  4_preprocess/${align_output}.tmp.sam >> 4_preprocess/${align_output}.6.tmp.sam 
samtools view -Sb 4_preprocess/${align_output}.6.tmp.sam > 4_preprocess/${align_output}.6.tmp.bam 
samtools index 4_preprocess/${align_output}.6.tmp.bam 

umi_tools group -I 4_preprocess/${align_output}.6.tmp.bam  \
        --group-out=4_preprocess/${align_output}.umi_groups.tsv \
        --umi-separator="__" \
        --output-bam > 4_preprocess/${align_output}.6_deduplicated.group.tmp.bam \
        --log=group.log --paired | samtools view
umi_tools dedup -I 4_preprocess/${align_output}.6.tmp.bam --umi-separator="__" \ --paired --output-stats=deduplicated --output-stats=4_preprocess/${align_output}.6_deduplicated -S 4_preprocess/${align_output}.6.deduplicated.tmp.bam

#umis with length 7
samtools view -H  4_preprocess/${align_output}.sorted.bam >  4_preprocess/${align_output}.7.tmp.sam 
grep -P "^.*__.{7}\t"  4_preprocess/${align_output}.tmp.sam >> 4_preprocess/${align_output}.7.tmp.sam 
samtools view -Sb 4_preprocess/${align_output}.7.tmp.sam > 4_preprocess/${align_output}.7.tmp.bam 
samtools index 4_preprocess/${align_output}.7.tmp.bam 

umi_tools group -I 4_preprocess/${align_output}.7.tmp.bam  \
        --group-out=4_preprocess/${align_output}.umi_groups.tsv \
        --umi-separator="__" \
        --output-bam > 4_preprocess/${align_output}.7_deduplicated.group.tmp.bam \
        --log=group.log --paired | samtools view
umi_tools dedup -I 4_preprocess/${align_output}.7.tmp.bam --umi-separator="__" \ --paired --output-stats=deduplicated --output-stats=4_preprocess/${align_output}.7_deduplicated -S 4_preprocess/${align_output}.7.deduplicated.tmp.bam

#umis with length 8
samtools view -H  4_preprocess/${align_output}.sorted.bam >  4_preprocess/${align_output}.8.tmp.sam 
grep -P "^.*__.{8}\t"  4_preprocess/${align_output}.tmp.sam >> 4_preprocess/${align_output}.8.tmp.sam 
samtools view -Sb 4_preprocess/${align_output}.8.tmp.sam > 4_preprocess/${align_output}.8.tmp.bam 
samtools index 4_preprocess/${align_output}.8.tmp.bam 

umi_tools group -I 4_preprocess/${align_output}.8.tmp.bam  \
        --group-out=4_preprocess/${align_output}.umi_groups.tsv \
        --umi-separator="__" \
        --output-bam > 4_preprocess/${align_output}.8_deduplicated.group.tmp.bam \
        --log=group.log --paired | samtools view
umi_tools dedup -I 4_preprocess/${align_output}.8.tmp.bam --umi-separator="__" \ --paired --output-stats=deduplicated --output-stats=4_preprocess/${align_output}.8_deduplicated -S 4_preprocess/${align_output}.8.deduplicated.tmp.bam

#merge bamfiles
samtools merge $outdir/4_preprocess/${align_output}_deduplicated.bam -c -f 4_preprocess/${align_output}.6.deduplicated.tmp.bam 4_preprocess/${align_output}.7.deduplicated.tmp.bam 4_preprocess/${align_output}.8.deduplicated.tmp.bam

#filtering for reads with MAPQ >= 20 + indexing bam file
samtools view -h -q 20 $outdir/4_preprocess/${align_output}_deduplicated.bam | samtools sort -o $outdir/4_preprocess/${align_output}_deduplicated.q20.bam
samtools index $outdir/4_preprocess/${align_output}_deduplicated.q20.bam

#removing redundant files
rm 4_preprocess/*${align_output}*.sorted.bam
rm 4_preprocess/*${align_output}.*.tmp.bam 
rm 4_preprocess/*${align_output}*tmp.sam 
rm 4_preprocess/*${align_output}*umi_groups.tsv


#step 4 - preprocess to 
output_bam=$outdir/4_preprocess/${align_output}_deduplicated.q20.bam

java -jar $picard_dir/MarkDuplicates.jar \
    I=$output_bam \
    O=4_preprocess/${align_output}.sorted.dedup.bam \
    M=4_preprocess/${align_output}.sorted.dedup.metrics \
    ASSUME_SORTED=true \
    VALIDATION_STRINGENCY=SILENT \
    REMOVE_DUPLICATES=false

rm 4_preprocess/${align_output}.sorted.dedup.bam


# step 5: get some alignment metrics
echo -e "\n~ CollectMultipleMetrics ~"
java -jar $picard_dir/CollectMultipleMetrics.jar I=$output_bam R=$bt2ref O=5_picard/${align_output} 


# get Thalia stats
samtools view $output_bam | cut -f 3 | sort | uniq -c | sort -nr | sed -e 's/^ *//;s/ /\t/' | awk 'OFS="\t" {print $2,$1}' | sort -n -k1,1 > $outdir/7_thalia_qc/${align_output}_thalia.counts
total=$(samtools view $output_bam | wc -l)
unmap=$(cat $outdir/7_thalia_qc/${align_output}_thalia.counts | grep "^\*" | cut -f2); if [[ -z $unmap ]]; then unmap="0"; fi
methyl=$(cat $outdir/7_thalia_qc/${align_output}_thalia.counts | grep F19K16 | cut -f2); if [[ -z $methyl ]]; then methyl="0"; fi
unmeth=$(cat $outdir/7_thalia_qc/${align_output}_thalia.counts | grep F24B22 | cut -f2); if [[ -z $unmeth ]]; then unmeth="0"; fi
pct_thalia=$(echo "scale=3; ($methyl + $unmeth)/$total * 100" | bc -l); if [[ -z $pct_thalia ]]; then pct_thalia="0"; fi
bet_thalia=$(echo "scale=3; $methyl/($methyl + $unmeth)" | bc -l); if [[ -z $bet_thalia ]]; then bet_thalia="0"; fi
echo -e "total\tunmap\tmethyl\tunmeth\tPCT_THALIANA\tTHALIANA_BETA" > $outdir/7_thalia_qc/${align_output}_thalia_summary.txt
echo -e "$total\t$unmap\t$methyl\t$unmeth\t$pct_thalia\t$bet_thalia" >> $outdir/7_thalia_qc/${align_output}_thalia_summary.txt


#7 - medip windows
for window in 300; do
 echo -e "\n~ MEDIPS for window: $window ~"

/.mounts/labs/awadallalab/private/ncheng/softwares/R/bin/Rscript $medips_script \
    --basedir $outdir \
    --bamfile $output_bam \
    --samplename $runfile \
    --ws $window \
    --outdir $outdir/6_medips_qc/window_${window}

 # get coverage windows
   count0=$(awk '$1 == 0' 6_medips_qc/window_${window}/${runfile}_genome_count.txt | wc -l)
   count1=$(awk '$1 >= 1' 6_medips_qc/window_${window}/${runfile}_genome_count.txt | wc -l)
  count10=$(awk '$1 >= 10' 6_medips_qc/window_${window}/${runfile}_genome_count.txt | wc -l)
  count50=$(awk '$1 >= 50' 6_medips_qc/window_${window}/${runfile}_genome_count.txt | wc -l)
 count100=$(awk '$1 >= 100' 6_medips_qc/window_${window}/${runfile}_genome_count.txt | wc -l)
 echo -e "sample\tcount0\tcount1\tcount10\tcount50\tcount100" > 6_medips_qc/window_${window}/${runfile}_coverage_windows.txt
 echo -e "$NAME\t$count0\t$count1\t$count10\t$count50\t$count100" >> 6_medips_qc/window_${window}/${runfile}_coverage_windows.txt


 # convert to bed
 $BEDOPS/convert2bed -i wig < 6_medips_qc/window_${window}/${runfile}_medips_count.wig > 6_medips_qc/window_${window}/${runfile}_${window}_medips_count.bed

 echo 'window_count' > 6_medips_qc/window_${window}/medips_${runfile}_${window}_window_count.txt
 awk -F"\t" '$5>0' 6_medips_qc/window_${window}/${runfile}_${window}_medips_rpkm.bed | wc -l  >> 6_medips_qc/window_${window}/medips_${runfile}_${window}_window_count.txt

 # cobble together QC data
 echo -e "\n~ Cobbling ~"
 echo -e "samples\n$align_output" > $outdir/8_qc_metrics/${runfile}.name.txt
 paste <(cat $outdir/8_qc_metrics/${runfile}.name.txt) \
       <(cut -f2-14 6_medips_qc/window_${window}/${runfile}_enrichment_data.txt) \
       <(cut -f2-3 6_medips_qc/window_${window}/${runfile}_coverage_counts.txt) \
       <(cut -f2-6 6_medips_qc/window_${window}/${runfile}_coverage_windows.txt) \
       <(cut -f2-5 6_medips_qc/window_${window}/${runfile}_saturation_metrics.txt) \
       <(awk 'NR==7 || NR==8' 4_preprocess/${align_output}.sorted.dedup.metrics) \
       <(awk 'NR==7 || NR==8' 5_picard/${align_output}.gc_bias_metrics.txt) \
       <(awk 'NR==7 || NR==8' 5_picard/${align_output}.alignment_summary_metrics) \
       <(cut -f1 6_medips_qc/window_${window}/medips_${runfile}_${window}_window_count.txt) \
           <(cat $outdir/7_thalia_qc/${align_output}_thalia_summary.txt) >  $outdir/8_qc_metrics/medips_${runfile}_qc_metrics_${window}.txt

done
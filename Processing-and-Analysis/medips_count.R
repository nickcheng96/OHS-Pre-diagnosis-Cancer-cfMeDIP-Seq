#!/usr/bin/env Rscript
#loading libraries
library(GenomicAlignments) #
library(GenomicRanges) #
library(getopt)#
library(optparse)#
library(BSgenome.Hsapiens.UCSC.hg38)#
library(gtools)#


#genomic filters to remove encode blacklist and gap regions
filters = readRDS('/.mounts/labs/awadallalab/private/ncheng/cfmedip_data/cptp_samples/fragmentation/hg38/excluded.regions.RDS')

#setting directories for bam files and output files
#bamdir='/.mounts/labs/awadallalab/private/ncheng/cfmedip_data/cptp_samples/batchall/bam.ln/q20.100424/'
#outdir='/.mounts/labs/awadallalab/private/ncheng/cfmedip_data/cptp_samples/fragmentation/aix13.updated1/all.sample.fragmentation/'
dir.create(outdir,recursive = T) #creating output directory in case this is ran for the first time

option_list = list(
  make_option(c("-i", "--bamfile"), type="character", default=NULL, help="input bam file", metavar="character"),
  make_option(c("-b", "--bamdir"), type="character", default=NULL, help="bam directory", metavar="character"),
  make_option(c("-o", "--outdir"), type="character", default=NULL, help="output directory", metavar="character")
)
#bamfile
#bamdir
#outdir

basedir <- opt$basedir
bamfile <- opt$bamfile
outdir <- paste0(opt$outdir,'/')


### Read GAlignmentPairs ###
bamfile = file.path(bamdir, bamfile)
param = ScanBamParam(flag = scanBamFlag(isDuplicate = NA,
                                        isSecondaryAlignment = FALSE,
                                        isUnmappedQuery = FALSE,
                                        isPaired = T,
                                        isFirstMateRead = T,
                                        isSupplementaryAlignment = F),
                     mapqFilter = 20,what = c('isize','seq'))


#extracting sample id
grpid = gsub('-','_',gsub('.*AIX','AIX',gsub('_deduplicated.q20.bam','',bamfile)))

#reading bam file into R as GAlignment object
galp.paired.raw= readGAlignments(bamfile, param = param)

#removing blacklist regions and gap regions
galp.paired.raw = galp.paired.raw[-queryHits(findOverlaps(galp.paired.raw, filters)),]

#removing unpaired reads or reads with no insert size reported
galp.paired.raw = galp.paired.raw[!is.na(elementMetadata(galp.paired.raw)$isize),]

#filtering for autosomal regions
chr.select=paste0("chr",c(1:22))
galp.paired.df = data.frame(galp.paired.raw)


### qc summary for sample ###

#large fragment frequency + other qc
fragment.qc.information = data.frame(sample.id = grpid)

#estimating methylated + unmethylated spike-in proportion
galp.paired.spikein = galp.paired.raw[!galp.paired.df$seq.names %in% chr.select,] #retaining only spike-in A. thaliana chromosomes
galp.paired.mspikein = galp.paired.raw[!galp.paired.df$seq.names %in% 'F19K16',] #methylated spike-in read count
galp.paired.uspikein = galp.paired.raw[!galp.paired.df$seq.names %in% 'F24B22',] #unmethylated spike-in read count
fragment.qc.information$all.spikein.count = length(galp.paired.spikein)
fragment.qc.information$m.spikein.count  = length(galp.paired.mspikein)
fragment.qc.information$u.spikein.count = length(galp.paired.uspikein)
fragment.qc.information$mspikein.proportion = fragment.qc.information$m.spikein.count/length(galp.paired.spikein)
fragment.qc.information$total.fragment.count = length(galp.paired.raw)

#estimating fragment size distributions 
galp.paired.short.mono = galp.paired.raw[abs(elementMetadata(galp.paired.raw)$isize) <= 149 & abs(elementMetadata(galp.paired.raw)$isize) >= 100,]
galp.paired.long.mono = galp.paired.raw[abs(elementMetadata(galp.paired.raw)$isize) <= 220 & abs(elementMetadata(galp.paired.raw)$isize) >= 150,]
galp.paired.di = galp.paired.raw[abs(elementMetadata(galp.paired.raw)$isize) <= 400 & abs(elementMetadata(galp.paired.raw)$isize) > 250,]
galp.paired.tri = galp.paired.raw[abs(elementMetadata(galp.paired.raw)$isize) <= 600 & abs(elementMetadata(galp.paired.raw)$isize) > 450,]
galp.paired.long.1_10k = galp.paired.raw[abs(elementMetadata(galp.paired.raw)$isize) <= 10000 & abs(elementMetadata(galp.paired.raw)$isize) >= 1000,]
galp.paired.long.10k = galp.paired.raw[abs(elementMetadata(galp.paired.raw)$isize) > 10000,]

fragment.qc.information$mono.short.fragment.count = length(galp.paired.short.mono)
fragment.qc.information$mono.long.fragment.count = length(galp.paired.long.mono)
fragment.qc.information$di.fragment.count = length(galp.paired.di)
fragment.qc.information$tri.fragment.count = length(galp.paired.tri)
fragment.qc.information$long.1_10k.count = length(galp.paired.long.1_10k)
fragment.qc.information$long.10k.count = length(galp.paired.long.10k)

fragment.qc.information$mono.short.fragment.proportion = length(galp.paired.short.mono)/length(galp.paired.raw)
fragment.qc.information$mono.long.fragment.proportion = length(galp.paired.long.mono)/length(galp.paired.raw)
fragment.qc.information$di.fragment.proportion = length(galp.paired.di)/length(galp.paired.raw)
fragment.qc.information$tri.fragment.proportion = length(galp.paired.tri)/length(galp.paired.raw)
fragment.qc.information$long.1_10k.proportion = length(galp.paired.long.1_10k)/length(galp.paired.raw)
fragment.qc.information$long.10k.proportion = length(galp.paired.long.10k)/length(galp.paired.raw)

saveRDS(fragment.qc.information,paste0(outdir,grpid,'.fragment.qc.RDS'))


### computing counts across 300bp ### 
#reading in genomic coordinates for 300 bp nonoverlapping tiled bins
tiles.300 = readRDS('/.mounts/labs/awadallalab/private/ncheng/references/hg38/tiles/hg38.tile.300.RDS') 

#selecting for autosomal regions
galp.paired.raw =galp.paired.raw[galp.paired.df$seqnames %in% chr.select,]
galp.paired.raw = galp.paired.raw[abs(elementMetadata(galp.paired.raw)$isize) >= 100 & abs(elementMetadata(galp.paired.raw)$isize) <= 1000,] #limintg to fragments of expected size range
tiles.list = c(tiles.300 = tiles.300)

#adjusting fragment start and end positions match insert length
gr.pos = granges(galp.paired.raw[elementMetadata(galp.paired.raw)$isize > 0,],use.mcols = TRUE)
new_start = start(gr.pos)
new_end =     new_start +gr.pos$isize
ranges(gr.pos) = IRanges(new_start, new_end)


gr.neg = granges(galp.paired.raw[elementMetadata(galp.paired.raw)$isize < 0,],use.mcols = TRUE)
new_start = end(gr.neg)  + gr.neg$isize
new_end =     end(gr.neg)
ranges(gr.neg) = IRanges(new_start, new_end)
galp.paired = c(gr.pos,gr.neg)

galp.paired.raw.inserts.1000 = galp.paired
galp.paired.raw.inserts.400 = galp.paired.raw.inserts.1000[abs(elementMetadata(galp.paired.raw.inserts.1000)$isize) > 70 & abs(elementMetadata(galp.paired.raw.inserts.1000)$isize) < 400,]
galp.list = c(inserts.1000 =galp.paired.raw.inserts.1000 ) #inserts.400 =galp.paired.raw.inserts.400

#coverage across 300bp genomic bins
tiles.300 = readRDS('/.mounts/labs/awadallalab/private/ncheng/references/hg38/tiles/hg38.tile.300.RDS')
targ.grange = tiles.300
targ.df = galp.paired.raw.inserts.1000
overlaps = findOverlaps(targ.grange, targ.df,minoverlap=10)
overlaps.df = data.frame(overlaps)
tss.grange.df = data.frame(targ.grange,stringsAsFactors = F)
tss.grange.df$window = paste0(tss.grange.df[,1],':',tss.grange.df[,2],'-',tss.grange.df[,3])
olap.genefreq = data.frame(table(tss.grange.df[overlaps.df$queryHits,'window']))
colnames(olap.genefreq)[1] = 'window'
tmp = merge(tss.grange.df[,c('window','seqnames')],olap.genefreq,by = 'window',all.x = T)[,-2]
tmp = tmp[order(tmp$window),]
tmp[is.na(tmp)] = 0
colnames(tmp) = c('window',grpid)
rownames(tmp) = tmp$window
saveRDS(tmp,paste0(outdir,grpid,'.inserts.1000.all.300.q20.RDS'))


#enhancer coverage
genhancer = readRDS('/.mounts/labs/awadallalab/private/ncheng/references/enhancers/genehancer/genhancer.grange.RDS')
targ.grange = genhancer
return.df = NULL
targ.df =galp.paired.raw.inserts.1000
overlaps = findOverlaps(targ.grange, targ.df,minoverlap=10)
overlaps.df = data.frame(overlaps)
tss.grange.df = data.frame(targ.grange,stringsAsFactors = F)
tss.grange.df$window = paste0(tss.grange.df[,1],':',tss.grange.df[,2],'-',tss.grange.df[,3])
olap.genefreq = data.frame(table(tss.grange.df[overlaps.df$queryHits,'window']))
colnames(olap.genefreq)[1] = 'window'
tmp = merge(tss.grange.df[,c('window','seqnames')],olap.genefreq,by = 'window',all.x = T)[,-2]
tmp = tmp[order(tmp$window),]
tmp[is.na(tmp)] = 0
tmp = unique(tmp)
colnames(tmp) = c('window',grpid)
rownames(tmp) = tmp$window
saveRDS(tmp,paste0(outdir,grpid,'.genhancer.inserts.1000.q20.cgifilt.RDS'))

#silencer coverage
silencers = readRDS('/.mounts/labs/awadallalab/private/ncheng/references/silencers/candidate.silencers.hg38.grange.RDS') #silencers
targ.grange = silencers
return.df = NULL
targ.df = galp.paired.raw.inserts.1000
overlaps = findOverlaps(targ.grange, targ.df,minoverlap=10)
overlaps.df = data.frame(overlaps)
tss.grange.df = data.frame(targ.grange,stringsAsFactors = F)
tss.grange.df$window = paste0(tss.grange.df[,1],':',tss.grange.df[,2],'-',tss.grange.df[,3])
olap.genefreq = data.frame(table(tss.grange.df[overlaps.df$queryHits,'window']))
colnames(olap.genefreq)[1] = 'window'
tmp = merge(tss.grange.df[,c('window','seqnames')],olap.genefreq,by = 'window',all.x = T)[,-2]
tmp = tmp[order(tmp$window),]
tmp[is.na(tmp)] = 0
tmp = unique(tmp)
colnames(tmp) = c('window',grpid)
rownames(tmp) = tmp$window
saveRDS(tmp,paste0(outdir,grpid,'.silencers.inserts.1000.q20.cgifilt.RDS'))
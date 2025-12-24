#!/usr/bin/env Rscript
library(IRanges)
library(BSgenome.Hsapiens.UCSC.hg38)
library(MEDIPS)
library(optparse)
library(gtools)
MEDIPS.CpGenrichNew <-function(file=NULL, BSgenome=NULL, extend=0, shift=0, uniq=1e-3, chr.select=NULL, paired=F){
  
  ## Proof correctness....
  if(is.null(BSgenome)){stop("Must specify a BSgenome library.")}
  
  ## Read region file
  fileName=unlist(strsplit(file, "/"))[length(unlist(strsplit(file, "/")))]
  path=paste(unlist(strsplit(file, "/"))[1:(length(unlist(strsplit(file, "/"))))-1], collapse="/")
  if(path==""){path=getwd()}
  if(!fileName%in%dir(path)){stop(paste("File", fileName, " not found in", path, sep =" "))}
  
  dataset = get(ls(paste("package:", BSgenome, sep = "")))
  
  if(!paired){GRange.Reads = getGRange(fileName, path, extend, shift, chr.select, dataset, uniq)}       else{GRange.Reads = getPairedGRange(fileName, path, extend, shift, chr.select, dataset, uniq)}
  
  ## Sort chromosomes
  if(length(unique(seqlevels(GRange.Reads)))>1){chromosomes=mixedsort(unique(seqlevels(GRange.Reads)))}
  if(length(unique(seqlevels(GRange.Reads)))==1){chromosomes=unique(seqlevels(GRange.Reads))}
  
  ## Get chromosome lengths for all chromosomes within data set.
  cat(paste("Loading chromosome lengths for ",BSgenome, "...\n", sep=""))
  
  chr_lengths=as.numeric(seqlengths(dataset)[chromosomes])
  
  ranges(GRange.Reads) <- restrict(ranges(GRange.Reads),+1)
  
  ##Calculate CpG density for regions
  total=length(chromosomes)
  cat("Calculating CpG density for given regions...\n")
  
  readsChars <- unlist(getSeq(dataset, GRange.Reads, as.character=TRUE))
  
  regions.CG = sum(vcountPattern("CG",readsChars))
  regions.C  = sum(vcountPattern("C",readsChars))
  regions.G  = sum(vcountPattern("G",readsChars))
  all.genomic= sum(width(readsChars))
  
  nReads <- length(readsChars)
  
  regions.relH=as.numeric(regions.CG)/as.numeric(all.genomic)*100
  regions.GoGe=(as.numeric(regions.CG)*as.numeric(all.genomic))/(as.numeric(regions.C)*as.numeric(regions.G))
  
  CG <- DNAStringSet("CG")
  pdict0 <- PDict(CG)
  params <- new("BSParams", X = dataset, FUN = countPDict, simplify = TRUE, exclude = c("rand", "chrUn"))
  genome.CG=sum(bsapply(params, pdict = pdict0))
  params <- new("BSParams", X = dataset, FUN = alphabetFrequency, exclude = c("rand", "chrUn"), simplify=TRUE)
  alphabet=bsapply(params)
  genome.l=sum(as.numeric(alphabet))
  genome.C=as.numeric(sum(alphabet[2,]))
  genome.G=as.numeric(sum(alphabet[3,]))
  genome.relH=genome.CG/genome.l*100
  genome.GoGe=(genome.CG*genome.l)/(genome.C*genome.G);
  
  ##Calculate CpG density for reference genome
  
  enrichment.score.relH=regions.relH/genome.relH
  enrichment.score.GoGe=regions.GoGe/genome.GoGe
  
  gc()
  return(list(genome=BSgenome, regions.CG=regions.CG, regions.C=regions.C, regions.G=regions.G, regions.relH=regions.relH, regions.GoGe=regions.GoGe, genome.C=genome.C, genome.G=genome.G, genome.CG=genome.CG, genome.relH=genome.relH, genome.GoGe=genome.GoGe, enrichment.score.relH=enrichment.score.relH, enrichment.score.GoGe=enrichment.score.GoGe))
  
}
# command line options
option_list = list(
  make_option(c("-d", "--basedir"), type="character", default=NULL, help="base directory", metavar="character"),
  make_option(c("-i", "--bamfile"), type="character", default=NULL, help="input bam file", metavar="character"),
  make_option(c("-o", "--outdir"), type="character", default=NULL, help="output directory to be created", metavar="
character"),
  make_option(c("-g", "--BSgenome"), type="character", default="BSgenome.Hsapiens.UCSC.hg38", help="genome", metavar="character"),
  make_option(c("-u", "--uniq"), type="numeric", default=1e-3, help="uniq", metavar="numeric"),
  make_option(c("-e", "--extend"), type="numeric", default=0, help="extend", metavar="numeric"),
  make_option(c("-s", "--shift"), type="numeric", default=0, help="shift", metavar="numeric"),
  make_option(c("-w", "--ws"), type="numeric", default=300, help="ws", metavar="numeric"),
  make_option(c("-n", "--samplename"), type="character", default=NULL, help="ws", metavar="character"),
  make_option(c("-p", "--paired"), type="logical", default=T, help="ws", metavar="logical")
  
)

# get options
opt_parser <- OptionParser(option_list=option_list);
opt <- parse_args(opt_parser);

# set better variable names
basedir <- opt$basedir
bamfile <- opt$bamfile
outdir <- opt$outdir
BSgenome <- opt$BSgenome
uniq <- opt$uniq
extend <- opt$extend
shift <- opt$shift
ws <- opt$ws
samplename <- opt$samplename
paired = T
cigar = F
print(basedir)
# set the working dir
setwd(basedir)

# set chromosomes
chr.select=paste0("chr",c(1:22,"X","Y","M"))

# set sample name if it exists
if (exists("samplename")) {
  sample_name <- samplename
} else {
  sample_name <- basename(bamfile)
}

# print options to output
print("Running MEDIPS with the following options:")
print(opt)

# create output if does not exist
if(!dir.exists(outdir)){
  dir.create(outdir)
}

# run MEDIPS and extract counts
MeDIP.set=MEDIPS.createSet(file=bamfile, BSgenome=BSgenome, extend=extend, shift=shift, uniq=uniq, window_size=ws,chr.select=chr.select, paired = paired, simpleCigar = cigar)
#saveRDS(MeDIP.set,paste0(outdir,"/",sample_name, '_MeDIP.set.RDS'))
MEDIPS.exportWIG(Set=MeDIP.set, file=paste0(outdir, "/", samplename, "_medips_rpkm.wig"), format="rpkm", descr="")
MEDIPS.exportWIG(Set=MeDIP.set, file=paste0(outdir, "/", samplename, "_medips_count.wig"), format="count", descr="")
# coupling
CS <- MEDIPS.couplingVector(pattern="CG", refObj=MeDIP.set)

# get saturation metrics
sr <- MEDIPS.saturation(file=bamfile, BSgenome=BSgenome, uniq=uniq, extend=extend, shift=shift, window_size=ws, chr.select=chr.select)

# Coverage
cr <- MEDIPS.seqCoverage(file=bamfile, pattern="CG", BSgenome=BSgenome, extend=extend, shift=shift, uniq=uniq, chr.select=chr.select)

# CpG enrichment
er <- MEDIPS.CpGenrichNew(file=bamfile, BSgenome=BSgenome, extend=extend, shift=shift, uniq=uniq, chr.select=chr.select)

# write out counts
write.table(data.frame(MeDIP.set@genome_count), paste0(outdir, "/", samplename, "_genome_count.txt"), row.names=F, quote=F, col.names=F)

# write out saturation metrics
saturation_df <- data.frame(maxEstCorReads=sr$maxEstCor[1], maxEstCor=sr$maxEstCor[2], maxTruCorReads=sr$maxTruCor[1], maxTruCor=sr$maxTruCor[2])
rownames(saturation_df) <- sample_name
write.table(saturation_df, file=paste0(outdir, "/", samplename, "_saturation_metrics.txt"), sep="\t", row.names=T, quote=F, col.names=NA)

# write out coverage metrics
coverage_df <- data.frame(numberReadsCG=cr$numberReads, numberReadsWOCG=cr$numberReadsWO)
rownames(coverage_df) <- sample_name
write.table(coverage_df, file=paste0(outdir, "/", samplename,  "_coverage_counts.txt"), sep="\t", row.names=T, quote=F, col.names=NA)
saveRDS(cr[[1]],paste0(outdir,"/",sample_name, '_seqCoverage_sites.RDS'))
# write out enrichment metrics
enrichment_df <- t(data.frame(unlist(er)))
rownames(enrichment_df) <- sample_name
write.table(enrichment_df, file=paste0(outdir, "/", samplename,  "_enrichment_data.txt"), sep="\t", row.names=T, quote=F, col.names=NA)

# get windows
no.window <- NULL
for(i in 1:length(chr.select)){
  no.window <- c(no.window, ceiling(MeDIP.set@chr_lengths[i]/ws))
}
window.last <- NULL
for(i in 1:length(chr.select)){
  window.last <- c(window.last, sum(no.window[1:i]))
}
window.first <- window.last-no.window+1

# write out window co-ordinates
write.csv(data.frame(chr=chr.select, window.first, window.last, no.window=no.window), paste0(outdir, "/", samplename, "_MEDIPS_window_per_chr.csv"), row.names=F)
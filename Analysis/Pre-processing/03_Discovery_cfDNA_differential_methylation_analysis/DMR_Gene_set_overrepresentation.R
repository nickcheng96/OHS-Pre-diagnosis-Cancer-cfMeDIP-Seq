#####extracting enhancer target genes for GSOA####
library(GenomicRanges)
#savedir='/.mounts/labs/awadallalab/private/ncheng/manuscripts/early.cancer.risk/all.sample.dmr.cpg5/'
#figdir=paste0(savedir,'figures/')
#dir.create(figdir,recursive = T)


genhancer = readRDS('/.mounts/labs/awadallalab/private/ncheng/references/enhancers/genehancer/genhancer.grange.unreduced.RDS')
genhancer.annotation.tissue = read.table('/.mounts/labs/awadallalab/private/ncheng/references/enhancers/genehancer//GeneHancer_AnnotSV_tissues_v5.20.txt', sep='\t',stringsAsFactors = F,header=T)
genhancer.annotation.genes = read.table('/.mounts/labs/awadallalab/private/ncheng/references/enhancers/genehancer//GeneHancer_AnnotSV_gene_association_scores_v5.20.txt', sep='\t',stringsAsFactors = F,header=T)

all.features=readRDS('/.mounts/labs/awadallalab/private/ncheng/cfmedip_data/cptp_samples/thesis/figures//chapter3/breast.dmrs.cv.norm.all/median.dmr.deseq.RDS')
enhancer.promoter.gene.annotation= function(targ.windows,genhancer) {
  #convert targeted windows to grange
  tmp = data.frame(chr = gsub(':.*','',targ.windows),
                   start = gsub('-.*','',gsub('.*:','',targ.windows)),
                   end = gsub('.*-','',targ.windows),
                   window = targ.windows)
  
  colnames(tmp)[c(1:3)] = c('chrom','chromStart','chromEnd')
  
  
  targ.grange <- makeGRangesFromDataFrame(tmp,
                                          keep.extra.columns=T)
  
  genhancer.filt =data.frame(genhancer[queryHits(findOverlaps(genhancer, targ.grange))])   
  targ.enhancer.ids = genhancer.filt$GHid
  
  tissues = genhancer.annotation.tissue[genhancer.annotation.tissue$GHid %in% targ.enhancer.ids,]
  genes= genhancer.annotation.genes[genhancer.annotation.genes$GHid %in% targ.enhancer.ids,]
  window.hits.all =data.frame(targ.grange[subjectHits(findOverlaps(genhancer, targ.grange))]) 
  
  genhancer.elit = genhancer[genhancer$is_elite == 1,]
  window.hits.elit =data.frame(targ.grange[subjectHits(findOverlaps(genhancer.elit, targ.grange))]) 
  window.hits.elit$GHid = data.frame(genhancer.elit[queryHits(findOverlaps(genhancer.elit, targ.grange))])$GHid
  genhancer.annotation.genes.elite = genhancer.annotation.genes[genhancer.annotation.genes$is_elite == 1,]
  window.hits.elit = merge(window.hits.elit, unique(genhancer.annotation.genes.elite[,c('GHid','symbol')]),by='GHid')
  return(list(gene=genes,tissue=tissues,window.hits.all,window.hits.elit)) #reduce(filters.hg19)
}
sex = 'Male'


savedir='/.mounts/labs/awadallalab/private/ncheng/manuscripts/early.cancer.risk/manuscript.dmrs/'

for (sex in c('Male','Female')){
  savedir='/.mounts/labs/awadallalab/private/ncheng/manuscripts/early.cancer.risk/manuscript.dmrs/'
  
  res.df = readRDS(paste0(savedir,'discovery.',sex,'.nonage.adjusted.dmrs.inserts.1000.all.300.q20.RDS'))
  res.df = res.df[!grepl('chrX|chrY',res.df$window),]
  res.df = res.df[res.df$window %in% cpg_count[cpg_count$count>=5,'window'],]
  res.df.hyper.sig.raw = res.df[res.df$log2FoldChange > 0.25 & res.df$pvalue < 0.05,]
  res.df.hypo.sig.raw = res.df[res.df$log2FoldChange < -0.25 & res.df$pvalue < 0.05,]
  res.df.hyper.sig.raw = res.df.hyper.sig.raw[order(res.df.hyper.sig.raw$pvalue),]
  res.df.hypo.sig.raw = res.df.hypo.sig.raw[order(res.df.hypo.sig.raw$pvalue),]
  
  for (j in c(2000)) {
    savedir='~/tmp/'
    all.features.sig.hypo = res.df.hypo.sig.raw
    all.features.sig.hyper = res.df.hyper.sig.raw
    all.features.sig.hyper = all.features.sig.hyper[order(all.features.sig.hyper$pvalue),]
    all.features.sig.hypo = all.features.sig.hypo[order(all.features.sig.hypo$pvalue),]
    
    hyper.j = min(nrow(all.features.sig.hyper),j)
    hypo.j=min(nrow(all.features.sig.hypo),j)
    targ.windows.hyper =all.features.sig.hyper$window[1:hyper.j]
    targ.windows.hypo =all.features.sig.hypo$window[1:hypo.j]
    
    enhancer.genes.hyper = enhancer.promoter.gene.annotation(targ.windows.hyper,genhancer )
    #write.table(enhancer.genes.hyper[[1]]$symbol,paste0(savedir,sex,'.targ.genes.',j,'.hyper.raw.genhancer.txt'),sep='\t',row.names=F,quote=F,col.names=F)
    write.table(enhancer.genes.hyper[[1]][enhancer.genes.hyper[[1]]$is_elite == 1,'symbol'],paste0(savedir,sex,'.targ.genes.',j,'.hyper.raw.all.txt'),sep='\t',row.names=F,quote=F,col.names=F)
    
    enhancer.genes.hypo = enhancer.promoter.gene.annotation(targ.windows.hypo,genhancer)
    #write.table(enhancer.genes.hypo[[1]]$symbol,paste0(savedir,sex,'.targ.genes.',j,'.hypo.raw.genhancer.txt'),sep='\t',row.names=F,quote=F,col.names=F)
    write.table(enhancer.genes.hypo[[1]][enhancer.genes.hypo[[1]]$is_elite == 1,'symbol'],paste0(savedir,sex,'.targ.genes.',j,'.hypo.raw.all.txt'),sep='\t',row.names=F,quote=F,col.names=F)
    
    #GSOA with enhancers only/bivalent enhancer/promoter
    enhancer.genes.hyper = enhancer.promoter.gene.annotation(targ.windows.hyper,genhancer[genhancer$regulatory_element_type %in% c('Enhancer','Promoter/Enhancer'),] )
    write.table(enhancer.genes.hyper[[1]]$symbol,paste0(savedir,sex,'.targ.genes.',j,'.hyper.raw.enhancer.txt'),sep='\t',row.names=F,quote=F,col.names=F)
    write.table(enhancer.genes.hyper[[1]][enhancer.genes.hyper[[1]]$is_elite == 1,'symbol'],paste0(savedir,sex,'.targ.genes.',j,'.hyper.raw.enhancerpe.txt'),sep='\t',row.names=F,quote=F,col.names=F)
    
    enhancer.genes.hypo = enhancer.promoter.gene.annotation(targ.windows.hypo,genhancer[genhancer$regulatory_element_type %in% c('Enhancer','Promoter/Enhancer'),] )
    write.table(enhancer.genes.hypo[[1]]$symbol,paste0(savedir,sex,'.targ.genes.',j,'.hypo.raw.enhancer.txt'),sep='\t',row.names=F,quote=F,col.names=F)
    write.table(enhancer.genes.hypo[[1]][enhancer.genes.hypo[[1]]$is_elite == 1,'symbol'],paste0(savedir,sex,'.targ.genes.',j,'.hypo.raw.enhancerpe.txt'),sep='\t',row.names=F,quote=F,col.names=F)
    
    #no bivalent
    enhancer.genes.hyper = enhancer.promoter.gene.annotation(targ.windows.hyper,genhancer[genhancer$regulatory_element_type %in% c('Enhancer'),] )
    #write.table(enhancer.genes.hyper[[1]]$symbol,paste0(savedir,sex,'.targ.genes.',j,'.hyper.raw.enhancer.txt'),sep='\t',row.names=F,quote=F,col.names=F)
    write.table(enhancer.genes.hyper[[1]][enhancer.genes.hyper[[1]]$is_elite == 1,'symbol'],paste0(savedir,sex,'.targ.genes.',j,'.hyper.raw.enhancer.txt'),sep='\t',row.names=F,quote=F,col.names=F)
    
    enhancer.genes.hypo = enhancer.promoter.gene.annotation(targ.windows.hypo,genhancer[genhancer$regulatory_element_type %in% c('Enhancer'),] )
    #write.table(enhancer.genes.hypo[[1]]$symbol,paste0(savedir,sex,'.targ.genes.',j,'.hypo.raw.enhancer.txt'),sep='\t',row.names=F,quote=F,col.names=F)
    write.table(enhancer.genes.hypo[[1]][enhancer.genes.hypo[[1]]$is_elite == 1,'symbol'],paste0(savedir,sex,'.targ.genes.',j,'.hypo.raw.enhancer.txt'),sep='\t',row.names=F,quote=F,col.names=F)
    
    #only bivalent
    enhancer.genes.hyper = enhancer.promoter.gene.annotation(targ.windows.hyper,genhancer[genhancer$regulatory_element_type %in% c('Promoter/Enhancer'),] )
    #write.table(enhancer.genes.hyper[[1]]$symbol,paste0(savedir,sex,'.targ.genes.',j,'.hyper.raw.pe.txt'),sep='\t',row.names=F,quote=F,col.names=F)
    write.table(enhancer.genes.hyper[[1]][enhancer.genes.hyper[[1]]$is_elite == 1,'symbol'],paste0(savedir,sex,'.targ.genes.',j,'.hyper.raw.pe.txt'),sep='\t',row.names=F,quote=F,col.names=F)
    
    enhancer.genes.hypo = enhancer.promoter.gene.annotation(targ.windows.hypo,genhancer[genhancer$regulatory_element_type %in% c('Promoter/Enhancer'),] )
    #write.table(enhancer.genes.hypo[[1]]$symbol,paste0(savedir,sex,'.targ.genes.',j,'.hypo.raw.pe.txt'),sep='\t',row.names=F,quote=F,col.names=F)
    write.table(enhancer.genes.hypo[[1]][enhancer.genes.hypo[[1]]$is_elite == 1,'symbol'],paste0(savedir,sex,'.targ.genes.',j,'.hypo.raw.pe.txt'),sep='\t',row.names=F,quote=F,col.names=F)
    
    #GSOA with promoters only
    enhancer.genes.hyper = enhancer.promoter.gene.annotation(targ.windows.hyper,genhancer[genhancer$regulatory_element_type %in% c('Promoter','Promoter/Enhancer'),] )
    #write.table(enhancer.genes.hyper[[1]]$symbol,paste0(savedir,sex,'.targ.genes.',j,'.hyper.raw.promoter1.txt'),sep='\t',row.names=F,quote=F,col.names=F)
    write.table(enhancer.genes.hyper[[1]][enhancer.genes.hyper[[1]]$is_elite == 1,'symbol'],paste0(savedir,sex,'.targ.genes.',j,'.hyper.raw.promoterpe.txt'),sep='\t',row.names=F,quote=F,col.names=F)
    
    enhancer.genes.hypo = enhancer.promoter.gene.annotation(targ.windows.hypo,genhancer[genhancer$regulatory_element_type %in% c('Promoter','Promoter/Enhancer'),] )
    #write.table(enhancer.genes.hypo[[1]]$symbol,paste0(savedir,sex,'.targ.genes.',j,'.hypo.raw.promoter1.txt'),sep='\t',row.names=F,quote=F,col.names=F)
    write.table(enhancer.genes.hypo[[1]][enhancer.genes.hypo[[1]]$is_elite == 1,'symbol'],paste0(savedir,sex,'.targ.genes.',j,'.hypo.raw.promoterpe.txt'),sep='\t',row.names=F,quote=F,col.names=F)
    
    #GSOA with promoters only
    enhancer.genes.hyper = enhancer.promoter.gene.annotation(targ.windows.hyper,genhancer[genhancer$regulatory_element_type %in% c('Promoter'),] )
    #write.table(enhancer.genes.hyper[[1]]$symbol,paste0(savedir,sex,'.targ.genes.',j,'.hyper.raw.promoter2.txt'),sep='\t',row.names=F,quote=F,col.names=F)
    write.table(enhancer.genes.hyper[[1]][enhancer.genes.hyper[[1]]$is_elite == 1,'symbol'],paste0(savedir,sex,'.targ.genes.',j,'.hyper.raw.promoter.txt'),sep='\t',row.names=F,quote=F,col.names=F)
    
    enhancer.genes.hypo = enhancer.promoter.gene.annotation(targ.windows.hypo,genhancer[genhancer$regulatory_element_type %in% c('Promoter'),] )
    #write.table(enhancer.genes.hypo[[1]]$symbol,paste0(savedir,sex,'.targ.genes.',j,'.hypo.raw.promoter2.txt'),sep='\t',row.names=F,quote=F,col.names=F)
    write.table(enhancer.genes.hypo[[1]][enhancer.genes.hypo[[1]]$is_elite == 1,'symbol'],paste0(savedir,sex,'.targ.genes.',j,'.hypo.raw.promoter.txt'),sep='\t',row.names=F,quote=F,col.names=F)

    
  }
  
  
}


enhancer.genes.hyper = enhancer.promoter.gene.annotation(res.df.hyper.sig.fdr$window )
write.table(enhancer.genes.hyper[[1]]$symbol,'~/brca.targ.genes.all.hyper.fdr.txt',sep='\t',row.names=F,quote=F,col.names=F)
write.table(enhancer.genes.hyper[[1]][enhancer.genes.hyper[[1]]$is_elite == 1,'symbol'],'~/brca.targ.genes.elit.hyper.fdr.txt',sep='\t',row.names=F,quote=F,col.names=F)


####extracting silencer target genes for GSOA####
#silencers
gene.coords = readRDS('/.mounts/labs/awadallalab/private/ncheng/references/gene_locations/hg38/gencode.hg38.gene.coordinates.RDS')

promoter.hic1 = read.table('/.mounts/labs/awadallalab/private/ncheng/references/silencers/hic.annotation/emat.2323/TS5_CD34_promoter-other_significant_interactions.txt',sep='\t',stringsAsFactors = F,header=T)
promoter.hic2 = read.table('/.mounts/labs/awadallalab/private/ncheng/references/silencers/hic.annotation/emat.2323/TS5_GM12878_promoter-other_significant_interactions.txt',sep='\t',stringsAsFactors = F,header=T)

silencer.gene.extract.v1 = function(targ.windows, promoter.hic1){
  #convert targeted windows to grange
  tmp = data.frame(chr = gsub(':.*','',targ.windows),
                   start = gsub('-.*','',gsub('.*:','',targ.windows)),
                   end = gsub('.*-','',targ.windows),
                   window = targ.windows)
  
  colnames(tmp)[c(1:3)] = c('chrom','chromStart','chromEnd')
  
  
  targ.grange <- makeGRangesFromDataFrame(tmp,
                                          keep.extra.columns=T)
  
  #
  hic.df = promoter.hic1
  colnames(hic.df)[1:3] = c('chrom','chromStart','chromEnd')
  hic.grange <- makeGRangesFromDataFrame(hic.df[,c(1:3,4,5)],
                                         keep.extra.columns=T)
  
  
  
  hic.filt =data.frame(hic.grange[queryHits(findOverlaps(hic.grange, targ.grange))])   
  return(hic.filt)
} 


promoter.hic3 = read.table('/.mounts/labs/awadallalab/private/ncheng/references/silencers/hic.annotation/GSE86821/GSM2309025_hESC_Sample811_Lane912_washU.txt.gz',sep='\t',stringsAsFactors = F,header=F)
promoter.hic4 = read.table('/.mounts/labs/awadallalab/private/ncheng/references/silencers/hic.annotation/GSE86821/GSM2309026_hESC_Samples1096_1390_Lanes1335_1785_washU.txt.gz',sep='\t',stringsAsFactors = F,header=F)
promoter.hic5 = read.table('/.mounts/labs/awadallalab/private/ncheng/references/silencers/hic.annotation/GSE86821/GSM2309027_hNPC_Sample1145_Lane1452_washU.txt.gz',sep='\t',stringsAsFactors = F,header=F)
promoter.hic6 = read.table('/.mounts/labs/awadallalab/private/ncheng/references/silencers/hic.annotation/GSE86821/GSM2309028_hNPC_Sample1742_Lane2174_washU.txt.gz',sep='\t',stringsAsFactors = F,header=F)

silencer.gene.extract.v2 = function(targ.windows, promoter.hic3,gene.coords){
  #convert targeted windows to grange
  tmp = data.frame(chr = gsub(':.*','',targ.windows),
                   start = gsub('-.*','',gsub('.*:','',targ.windows)),
                   end = gsub('.*-','',targ.windows),
                   window = targ.windows)
  
  colnames(tmp)[c(1:3)] = c('chrom','chromStart','chromEnd')
  
  
  targ.grange <- makeGRangesFromDataFrame(tmp,
                                          keep.extra.columns=T)
  
  #
  hic.df = promoter.hic3
  hic.df[,1] = paste0('chr',hic.df[,1])
  colnames(hic.df)[1:3] = c('chrom','chromStart','chromEnd')
  hic.grange <- makeGRangesFromDataFrame(hic.df[,c(1:4)],
                                         keep.extra.columns=T)
  
  
  
  hic.filt =data.frame(hic.grange[queryHits(findOverlaps(hic.grange, targ.grange))])   
  
  gene.regions= paste0('chr',gsub(',.*','',hic.filt[,6]))
  tmp.gene.regions = data.frame(chrom = gsub(':.*','',gene.regions),
                                chromStart = gsub('-.*','',gsub('.*:','',gene.regions)),
                                chromEnd = gsub('.*-','',gene.regions))
  targ.region.grange <- makeGRangesFromDataFrame(tmp.gene.regions,
                                                 keep.extra.columns=T)
  
  gene.annotation.df=gene.coords
  colnames(gene.annotation.df)[1:3]= c('chrom','chromStart','chromEnd')
  if (sum(grepl('chr',gene.annotation.df$chrom)) == 0) {
    gene.annotation.df$chrom = paste0('chr',gene.annotation.df$chrom)
  }
  
  gene.annotation.grange <- makeGRangesFromDataFrame(gene.annotation.df,
                                                     keep.extra.columns=T)
  
  hic.targ.gene =data.frame(gene.annotation.grange[queryHits(findOverlaps(gene.annotation.grange, targ.region.grange))])   
  
  hic.targ.gene.return = unique(hic.targ.gene[,c('ensgid','gene','gene.type')])
  
  return(hic.targ.gene.return)
} 

promoter.hic7=read.table('/.mounts/labs/awadallalab/private/ncheng/references/silencers/hic.annotation/encode/ENCFF693XIL.bedpe.gz',sep='\t',stringsAsFactors = F,header=F)
silencer.gene.extract.v3 = function(targ.windows, promoter.hic3,gene.coords){
  #convert targeted windows to grange
  tmp = data.frame(chr = gsub(':.*','',targ.windows),
                   start = gsub('-.*','',gsub('.*:','',targ.windows)),
                   end = gsub('.*-','',targ.windows),
                   window = targ.windows)
  
  colnames(tmp)[c(1:3)] = c('chrom','chromStart','chromEnd')
  
  
  targ.grange <- makeGRangesFromDataFrame(tmp,
                                          keep.extra.columns=T)
  
  #
  hic.df = promoter.hic3
  colnames(hic.df)[1:3] = c('chrom','chromStart','chromEnd')
  hic.grange <- makeGRangesFromDataFrame(hic.df[,c(1:3)],
                                         keep.extra.columns=T)
  
  
  
  hic.filt =data.frame(hic.df[queryHits(findOverlaps(hic.grange, targ.grange)),])   
  
  tmp.gene.regions= hic.filt[,c(4:6)]
  colnames(tmp.gene.regions)[1:3] =  c('chrom','chromStart','chromEnd')
  targ.region.grange <- makeGRangesFromDataFrame(tmp.gene.regions,
                                                 keep.extra.columns=T)
  
  gene.annotation.df=gene.coords
  colnames(gene.annotation.df)[1:3]= c('chrom','chromStart','chromEnd')
  if (sum(grepl('chr',gene.annotation.df$chrom)) == 0) {
    gene.annotation.df$chrom = paste0('chr',gene.annotation.df$chrom)
  }
  
  gene.annotation.grange <- makeGRangesFromDataFrame(gene.annotation.df,
                                                     keep.extra.columns=T)
  
  hic.targ.gene =data.frame(gene.annotation.grange[queryHits(findOverlaps(gene.annotation.grange, targ.region.grange))])   
  
  hic.targ.gene.return = unique(hic.targ.gene[,c('ensgid','gene','gene.type')])
  
  return(hic.targ.gene.return)
} 
promoter.hic8=read.table('/.mounts/labs/awadallalab/private/ncheng/references/silencers/hic.annotation/ENCFF202FID.bed.gz',sep='\t',stringsAsFactors = F,header=F)
promoter.hic9=read.table('/.mounts/labs/awadallalab/private/ncheng/references/silencers/hic.annotation/ENCFF128MZW.bed.gz',sep='\t',stringsAsFactors = F,header=F)

silencer.gene.extract.v4 = function(targ.windows, promoter.hic3,gene.coords){
  #convert targeted windows to grange
  tmp = data.frame(chr = gsub(':.*','',targ.windows),
                   start = gsub('-.*','',gsub('.*:','',targ.windows)),
                   end = gsub('.*-','',targ.windows),
                   window = targ.windows)
  
  colnames(tmp)[c(1:3)] = c('chrom','chromStart','chromEnd')
  
  
  targ.grange <- makeGRangesFromDataFrame(tmp,
                                          keep.extra.columns=T)
  
  #
  hic.df = promoter.hic3
  colnames(hic.df)[1:3] = c('chrom','chromStart','chromEnd')
  hic.grange <- makeGRangesFromDataFrame(hic.df[,c(1:3)],
                                         keep.extra.columns=T)
  
  
  
  hic.filt =data.frame(hic.df[queryHits(findOverlaps(hic.grange, targ.grange)),])   
  
  hic.targ.gene.return = unique(hic.filt[,c('V7','V6')])
  colnames(hic.targ.gene.return) = c('ensgid','gene')
  return(hic.targ.gene.return)
} 



#savedir='/.mounts/labs/awadallalab/private/ncheng/manuscripts/early.cancer.risk/all.sample.dmr.cpg5/'
#figdir=paste0(savedir,'figures/')
#dir.create(figdir,recursive = T)

combined.silencer=NULL
for (sex in c('Male','Female')) {
  
  res.df = readRDS(paste0(savedir,'discovery.',sex,'.nonage.adjusted.dmrs.inserts.1000.all.300.q20.RDS'))
  res.df = res.df[!grepl('chrX|chrY',res.df$window),]
  res.df.hyper.sig.raw = res.df[res.df$log2FoldChange > 0.25 & res.df$pvalue < 0.05,]
  res.df.hypo.sig.raw = res.df[res.df$log2FoldChange < -0.25 & res.df$pvalue < 0.05,]
  silencer.regions.hyper = res.df.hyper.sig.raw[order(res.df.hyper.sig.raw$pvalue),][1:2000,'window']
  silencer.regions.hypo = res.df.hypo.sig.raw[order(res.df.hypo.sig.raw$pvalue),][1:2000,'window']
  
  
  
  silencer.hic1.hyper =silencer.gene.extract.v1(silencer.regions.hyper,promoter.hic1)
  silencer.hic2.hyper =silencer.gene.extract.v1(silencer.regions.hyper,promoter.hic2)
  silencer.hic3.hyper = silencer.gene.extract.v2(silencer.regions.hyper, promoter.hic3,gene.coords )
  silencer.hic4.hyper = silencer.gene.extract.v2(silencer.regions.hyper, promoter.hic4,gene.coords )
  silencer.hic5.hyper = silencer.gene.extract.v2(silencer.regions.hyper, promoter.hic5,gene.coords )
  silencer.hic6.hyper = silencer.gene.extract.v2(silencer.regions.hyper, promoter.hic6,gene.coords )
  silencer.hic7.hyper = silencer.gene.extract.v3(silencer.regions.hyper, promoter.hic7,gene.coords )
  silencer.hic8.hyper = silencer.gene.extract.v4(silencer.regions.hyper, promoter.hic8,gene.coords )
  silencer.hic9.hyper = silencer.gene.extract.v4(silencer.regions.hyper, promoter.hic9,gene.coords )
  
  
  
  silencer.hic1.hypo =silencer.gene.extract.v1(silencer.regions.hypo,promoter.hic1)
  silencer.hic2.hypo =silencer.gene.extract.v1(silencer.regions.hypo,promoter.hic2)
  silencer.hic3.hypo = silencer.gene.extract.v2(silencer.regions.hypo, promoter.hic3,gene.coords )
  silencer.hic4.hypo = silencer.gene.extract.v2(silencer.regions.hypo, promoter.hic4,gene.coords )
  silencer.hic5.hypo = silencer.gene.extract.v2(silencer.regions.hypo, promoter.hic5,gene.coords )
  silencer.hic6.hypo = silencer.gene.extract.v2(silencer.regions.hypo, promoter.hic6,gene.coords )
  silencer.hic7.hypo = silencer.gene.extract.v3(silencer.regions.hypo, promoter.hic7,gene.coords )
  silencer.hic8.hypo = silencer.gene.extract.v4(silencer.regions.hypo, promoter.hic8,gene.coords )
  silencer.hic9.hypo = silencer.gene.extract.v4(silencer.regions.hypo, promoter.hic9,gene.coords )
  
  
  
  all.hyper = c(unique(silencer.hic1.hyper$Ensembl.Gene.ID), 
                unique(silencer.hic2.hyper$Ensembl.Gene.ID),
                unique(silencer.hic3.hyper$ensgid),
                unique(silencer.hic4.hyper$ensgid),
                unique(silencer.hic5.hyper$ensgid),
                unique(silencer.hic6.hyper$ensgid),
                unique(silencer.hic7.hyper$ensgid),
                unique(silencer.hic8.hyper$ensgid),
                unique(silencer.hic9.hyper$ensgid))
  
  all.hypo = c(unique(silencer.hic1.hypo$Ensembl.Gene.ID), 
               unique(silencer.hic2.hypo$Ensembl.Gene.ID),
               unique(silencer.hic3.hypo$ensgid),
               unique(silencer.hic4.hypo$ensgid),
               unique(silencer.hic5.hypo$ensgid),
               unique(silencer.hic6.hypo$ensgid),
               unique(silencer.hic7.hyper$ensgid),
               unique(silencer.hic8.hyper$ensgid),
               unique(silencer.hic9.hyper$ensgid))
  
  all.hyper = gsub('\\..*','',all.hyper)
  all.hyper = paste(all.hyper,sep='\\|')
  all.hyper= unlist(str_split(string = all.hyper,pattern = '\\|'))
  all.hyper =data.frame(table(all.hyper))
  all.hyper$direction = 'Top 2000 Hypermethylated Regions'
  all.hypo = gsub('\\..*','',all.hypo)
  all.hypo = base::paste(all.hypo,sep='\\|')
  all.hypo= unlist(str_split(string = all.hypo,pattern = '\\|'))
  all.hypo =data.frame(table(all.hypo))
  all.hypo$direction = 'Top 2000 Hypomethylated Regions'
  colnames(all.hypo)[1] = 'ensgid'
  colnames(all.hyper)[1] = 'ensgid'
  
  combined.silencer.tmp = rbind(all.hyper,all.hypo)
  combined.silencer.tmp$group = ifelse(sex=='Female','Incident Breast Cancer','Incident Prostate Cancer')
  combined.silencer =rbind(combined.silencer.tmp,combined.silencer)
}

saveRDS(combined.silencer,paste0('silencer.dmr.genes.RDS'))



silencer.hic1.background =silencer.gene.extract.v1(paste0(promoter.hic1[,1],':',promoter.hic1[,2],'-',promoter.hic1[,3]),promoter.hic1)
silencer.hic2.background=silencer.gene.extract.v1(paste0(promoter.hic2[,1],':',promoter.hic2[,2],'-',promoter.hic2[,3]),promoter.hic2)
silencer.hic3.background = silencer.gene.extract.v2(paste0('chr',promoter.hic3[,1],':',promoter.hic3[,2],'-',promoter.hic3[,3]), promoter.hic3,gene.coords )
silencer.hic4.background = silencer.gene.extract.v2(paste0('chr',promoter.hic4[,1],':',promoter.hic4[,2],'-',promoter.hic4[,3]), promoter.hic4,gene.coords )
silencer.hic5.background = silencer.gene.extract.v2(paste0('chr',promoter.hic5[,1],':',promoter.hic5[,2],'-',promoter.hic5[,3]), promoter.hic5,gene.coords )
silencer.hic6.background = silencer.gene.extract.v2(paste0('chr',promoter.hic6[,1],':',promoter.hic6[,2],'-',promoter.hic6[,3]), promoter.hic6,gene.coords )
silencer.hic7.background = silencer.gene.extract.v3(paste0(promoter.hic7[,1],':',promoter.hic7[,2],'-',promoter.hic7[,3]), promoter.hic7,gene.coords )
silencer.hic8.background = silencer.gene.extract.v4(paste0(promoter.hic8[,1],':',promoter.hic8[,2],'-',promoter.hic8[,3]), promoter.hic8,gene.coords )
silencer.hic9.background = silencer.gene.extract.v4(paste0(promoter.hic9[,1],':',promoter.hic9[,2],'-',promoter.hic9[,3]), promoter.hic9,gene.coords )

all.genes =c(unique(promoter.hic1$Ensembl.Gene.ID), 
             unique(promoter.hic2$Ensembl.Gene.ID),
             unique(silencer.hic3.background$ensgid),
             unique(silencer.hic4.background$ensgid),
             unique(silencer.hic5.background$ensgid),
             unique(silencer.hic6.background$ensgid),
             unique(silencer.hic7.background$ensgid),
             unique(silencer.hic8.background$ensgid),
             unique(silencer.hic9.background$ensgid))
all.genes = gsub('\\..*','',all.genes)

all.genes = paste(all.genes,sep='\\|')
library(stringr)
all.genes= unlist(str_split(string = all.genes,pattern = '\\|'))
all.genes =unique(all.genes)
write.table(all.genes,paste0(savedir,'silencer.genes.background.txt'),sep='\t',quote=F,row.names=F,col.names=F)



####local GSOA####
library(fgsea)
library(msigdbr)
library(biomaRt)
library(Hmisc)
library(tools)
library(ggplot2)
library(gridExtra)
library(ggh4x)
library(ggforce)
setwd('/Users/ncheng/Desktop/ncheng/Desktop/Lab Work/Manuscripts/OHS_cancer_risk/figures/gsoa/current1/')
#enhancer.background.raw = read.table('/Users/ncheng/Desktop/ncheng/Desktop/Lab Work/Manuscripts/PhD Thesis /Chapter4/figures/discovery.cv/v5/all.dmr.genes/discovery//targ.genes.all.background.genhancer.raw.txt',sep='\t',stringsAsFactors = F,header=F)
enhancer.background=readRDS('/Users/ncheng/Desktop/ncheng/Desktop/Lab Work/Manuscripts/OHS_cancer_risk/figures/gsoa/enhancer.background.ensgid.RDS')
#converting ensgid to entrezid
hsmart <- useEnsembl(dataset = "hsapiens_gene_ensembl", biomart = "ensembl")

gsea.calc = function(targ.genes, targ.background,name) {
  mapping.targ <- getBM(
    attributes = c('ensembl_gene_id', 'entrezgene_id', 'hgnc_symbol'), 
    filters = 'ensembl_gene_id',
    values = targ.genes,
    mart = hsmart
  )
  mapping.background <- getBM(
    attributes = c('ensembl_gene_id', 'entrezgene_id', 'hgnc_symbol'), 
    filters = 'ensembl_gene_id',
    values = targ.background,
    mart = hsmart
  )
  
  
  msigdbr_df.h <- msigdbr(species = "human", category = "H")
  msigdbr_df.c2.kegg <- msigdbr(species = "human", category = "C2", subcategory = "CP:KEGG")
  msigdbr_df.c2.reactome <- msigdbr(species = "human", category = "C2", subcategory = "CP:REACTOME")
  msigdbr_df.c2.wikipath <- msigdbr(species = "human", category = "C2", subcategory = "CP:WIKIPATHWAYS")
  
  msigdbr_df.c2.cgp <- msigdbr(species = "human", category = "C2", subcategory = "CGP")
  msigdbr_df.c2.cp <- msigdbr(species = "human", category = "C2", subcategory = "CP")
  msigdbr_df.c4.cgn <- msigdbr(species = "human", category = "C4", subcategory = "CGN")
  msigdbr_df.c4.cm <- msigdbr(species = "human", category = "C4", subcategory = "CM")
  msigdbr_df.c5.bp <- msigdbr(species = "human", category = "C5", subcategory = "BP")
  msigdbr_df.c5.mf <- msigdbr(species = "human", category = "C5", subcategory = "MF")
  msigdbr_df.c5.hpo <- msigdbr(species = "human", category = "C5", subcategory = "HPO")
  msigdbr_df.c6<- msigdbr(species = "human", category = "C6")
  msigdbr_df.c7<- msigdbr(species = "human", category = "C7")
  msigdbr_df.c8<- msigdbr(species = "human", category = "C8")
  
  gene.set.list = list(H = msigdbr_df.h,
                       C2.CP = msigdbr_df.c2.cp,
                       C4.CGN = msigdbr_df.c4.cgn,
                       C4.CM = msigdbr_df.c4.cm,
                       C5.BP = msigdbr_df.c5.bp,
                       C5.MF = msigdbr_df.c5.mf,
                       C5.HPO = msigdbr_df.c5.hpo,
                       C6=msigdbr_df.c6,
                       C7=msigdbr_df.c7,
                       C8=msigdbr_df.c8
  )
  
  
  gene.set.list = list(H = msigdbr_df.h,
                       C2.KEGG = msigdbr_df.c2.kegg,
                       C2.WIKIPATH = msigdbr_df.c2.wikipath,
                       C2.REACTOME = msigdbr_df.c2.reactome,
                       # C2.CP = msigdbr_df.c2.cp,
                       #C4.CGN = msigdbr_df.c4.cgn,
                       #C4.CM = msigdbr_df.c4.cm,
                       C5.BP = msigdbr_df.c5.bp,
                       C5.MF = msigdbr_df.c5.mf
                       
                       #C5.HPO = msigdbr_df.c5.hpo,
                       #C6=msigdbr_df.c6,
                       #C7=msigdbr_df.c7,
                       #C8=msigdbr_df.c8
  )
  
  return.df = NULL
  for (gs in names(gene.set.list)){
    msigdbr_df = gene.set.list[[gs]]
    pathwaysH = split(x = msigdbr_df$entrez_gene, f = msigdbr_df$gs_name)
    foraRes <- fora(genes=mapping.targ$entrezgene_id, 
                    universe=mapping.background$entrezgene_id, 
                    pathways=pathwaysH)
    
    foraRes.sig = foraRes[foraRes$padj < 0.1,]
    if (nrow(foraRes.sig) > 1) {
      collapseres = collapsePathwaysORA(foraRes.sig, 
                                        pathways=pathwaysH, 
                                        genes=mapping.targ$entrezgene_id,
                                        universe=mapping.background$entrezgene_id, pval.threshold = 0.05)
      
      foraRes = foraRes[foraRes$pathway %in% collapseres$mainPathways,]
      
    }
    
    foraRes$marker = name
    foraRes$gene.set = gs
    return.df = rbind(return.df, foraRes)
  }
  # fixing format to work with fgsea
  #reactome
  
  
  
  return(return.df)
  
}

enhancer.id = function(enhancer.hyper.all){
  mapping.targ <- getBM(
    attributes = c('ensembl_gene_id', 'entrezgene_id', 'hgnc_symbol'), 
    filters = 'hgnc_symbol',
    values = enhancer.hyper.all,
    mart = hsmart
  )
  return(mapping.targ[,1])
}
gsea.calc = function(targ.genes, targ.background,name) {
  mapping.targ <- getBM(
    attributes = c('ensembl_gene_id', 'entrezgene_id', 'hgnc_symbol'), 
    filters = 'ensembl_gene_id',
    values = targ.genes,
    mart = hsmart
  )
  mapping.background <- getBM(
    attributes = c('ensembl_gene_id', 'entrezgene_id', 'hgnc_symbol'), 
    filters = 'ensembl_gene_id',
    values = targ.background,
    mart = hsmart
  )
  
  
  msigdbr_df.h <- msigdbr(species = "human", category = "H")
  msigdbr_df.c2.kegg <- msigdbr(species = "human", category = "C2", subcategory = "CP:KEGG")
  msigdbr_df.c2.reactome <- msigdbr(species = "human", category = "C2", subcategory = "CP:REACTOME")
  msigdbr_df.c2.wikipath <- msigdbr(species = "human", category = "C2", subcategory = "CP:WIKIPATHWAYS")
  
  msigdbr_df.c2.cgp <- msigdbr(species = "human", category = "C2", subcategory = "CGP")
  msigdbr_df.c2.cp <- msigdbr(species = "human", category = "C2", subcategory = "CP")
  msigdbr_df.c4.cgn <- msigdbr(species = "human", category = "C4", subcategory = "CGN")
  msigdbr_df.c4.cm <- msigdbr(species = "human", category = "C4", subcategory = "CM")
  msigdbr_df.c5.bp <- msigdbr(species = "human", category = "C5", subcategory = "BP")
  msigdbr_df.c5.mf <- msigdbr(species = "human", category = "C5", subcategory = "MF")
  msigdbr_df.c5.hpo <- msigdbr(species = "human", category = "C5", subcategory = "HPO")
  msigdbr_df.c6<- msigdbr(species = "human", category = "C6")
  msigdbr_df.c7<- msigdbr(species = "human", category = "C7")
  msigdbr_df.c8<- msigdbr(species = "human", category = "C8")
  
  gene.set.list = list(H = msigdbr_df.h,
                       C2.CP = msigdbr_df.c2.cp,
                       C4.CGN = msigdbr_df.c4.cgn,
                       C4.CM = msigdbr_df.c4.cm,
                       C5.BP = msigdbr_df.c5.bp,
                       C5.MF = msigdbr_df.c5.mf,
                       C5.HPO = msigdbr_df.c5.hpo,
                       C6=msigdbr_df.c6,
                       C7=msigdbr_df.c7,
                       C8=msigdbr_df.c8
  )
  
  
  gene.set.list = list(H = msigdbr_df.h,
                       C2.KEGG = msigdbr_df.c2.kegg,
                       C2.WIKIPATH = msigdbr_df.c2.wikipath,
                       C2.REACTOME = msigdbr_df.c2.reactome,
                       # C2.CP = msigdbr_df.c2.cp,
                       #C4.CGN = msigdbr_df.c4.cgn,
                       #C4.CM = msigdbr_df.c4.cm,
                       C5.BP = msigdbr_df.c5.bp,
                       C5.MF = msigdbr_df.c5.mf
                       
                       #C5.HPO = msigdbr_df.c5.hpo,
                       #C6=msigdbr_df.c6,
                       #C7=msigdbr_df.c7,
                       #C8=msigdbr_df.c8
  )
  
  return.df = NULL
  for (gs in names(gene.set.list)){
    msigdbr_df = gene.set.list[[gs]]
    pathwaysH = split(x = msigdbr_df$entrez_gene, f = msigdbr_df$gs_name)
    foraRes <- fora(genes=mapping.targ$entrezgene_id, 
                    universe=mapping.background$entrezgene_id, 
                    pathways=pathwaysH)
    
    foraRes.sig = foraRes[foraRes$padj < 0.1,]
    if (nrow(foraRes.sig) > 1) {
      collapseres = collapsePathwaysORA(foraRes.sig, 
                                        pathways=pathwaysH, 
                                        genes=mapping.targ$entrezgene_id,
                                        universe=mapping.background$entrezgene_id, pval.threshold = 0.05)
      
      foraRes = foraRes[foraRes$pathway %in% collapseres$mainPathways,]
      
    }
    
    foraRes$marker = name
    foraRes$gene.set = gs
    return.df = rbind(return.df, foraRes)
  }
  # fixing format to work with fgsea
  #reactome
  
  
  
  return(return.df)
  
}

gsea.calc.all = function(targ.genes, targ.background,name) {
  mapping.targ <- getBM(
    attributes = c('ensembl_gene_id', 'entrezgene_id', 'hgnc_symbol'), 
    filters = 'ensembl_gene_id',
    values = targ.genes,
    mart = hsmart
  )
  mapping.background <- getBM(
    attributes = c('ensembl_gene_id', 'entrezgene_id', 'hgnc_symbol'), 
    filters = 'ensembl_gene_id',
    values = targ.background,
    mart = hsmart
  )
  
  
  msigdbr_df.h <- msigdbr(species = "human", category = "H")
  msigdbr_df.c2.kegg <- msigdbr(species = "human", category = "C2", subcategory = "CP:KEGG")
  msigdbr_df.c2.reactome <- msigdbr(species = "human", category = "C2", subcategory = "CP:REACTOME")
  msigdbr_df.c2.wikipath <- msigdbr(species = "human", category = "C2", subcategory = "CP:WIKIPATHWAYS")
  
  msigdbr_df.c2.cgp <- msigdbr(species = "human", category = "C2", subcategory = "CGP")
  msigdbr_df.c2.cp <- msigdbr(species = "human", category = "C2", subcategory = "CP")
  msigdbr_df.c4.cgn <- msigdbr(species = "human", category = "C4", subcategory = "CGN")
  msigdbr_df.c4.cm <- msigdbr(species = "human", category = "C4", subcategory = "CM")
  msigdbr_df.c5.bp <- msigdbr(species = "human", category = "C5", subcategory = "BP")
  msigdbr_df.c5.mf <- msigdbr(species = "human", category = "C5", subcategory = "MF")
  msigdbr_df.c5.hpo <- msigdbr(species = "human", category = "C5", subcategory = "HPO")
  msigdbr_df.c6<- msigdbr(species = "human", category = "C6")
  msigdbr_df.c7<- msigdbr(species = "human", category = "C7")
  msigdbr_df.c8<- msigdbr(species = "human", category = "C8")
  
  gene.set.list = list(H = msigdbr_df.h,
                       C2.CP = msigdbr_df.c2.cp,
                       C4.CGN = msigdbr_df.c4.cgn,
                       C4.CM = msigdbr_df.c4.cm,
                       C5.BP = msigdbr_df.c5.bp,
                       C5.MF = msigdbr_df.c5.mf,
                       C5.HPO = msigdbr_df.c5.hpo,
                       C6=msigdbr_df.c6,
                       C7=msigdbr_df.c7,
                       C8=msigdbr_df.c8
  )
  
  
  gene.set.list = list(H = msigdbr_df.h,
                       C2.KEGG = msigdbr_df.c2.kegg,
                       #C2.WIKIPATH = msigdbr_df.c2.wikipath,
                       C2.REACTOME = msigdbr_df.c2.reactome#,
                       # C2.CP = msigdbr_df.c2.cp,
                       #C4.CGN = msigdbr_df.c4.cgn,
                       #C4.CM = msigdbr_df.c4.cm,
                       #  C5.BP = msigdbr_df.c5.bp,
                       #  C5.MF = msigdbr_df.c5.mf
                       
                       #C5.HPO = msigdbr_df.c5.hpo,
                       #C6=msigdbr_df.c6,
                       #C7=msigdbr_df.c7,
                       #C8=msigdbr_df.c8
  )
  
  return.df = NULL
  for (gs in names(gene.set.list)){
    msigdbr_df = gene.set.list[[gs]]
    pathwaysH = split(x = msigdbr_df$entrez_gene, f = msigdbr_df$gs_name)
    foraRes <- fora(genes=mapping.targ$entrezgene_id, 
                    universe=mapping.background$entrezgene_id, 
                    pathways=pathwaysH)
    
    foraRes$marker = name
    foraRes$gene.set = gs
    return.df = rbind(return.df, foraRes)
  }
  # fixing format to work with fgsea
  #reactome
  
  
  
  return(return.df)
  
}

dmr.list = list( )#,

for (Sex in c('Male','Female')) {
  dmr.list[[Sex]] = list()
  dmr.list[[Sex]][[paste0('Hypermethylated')]] = list()
  dmr.list[[Sex]][[paste0('Hypomethylated')]] = list()
  bin300.genhancer.hyper.all.raw = read.table(paste0(Sex,'.targ.genes.2000.hyper.raw.enhancer.txt'),sep='\t',stringsAsFactors = F,header=F)
  bin300.genhancer.hypo.all.raw = read.table(paste0(Sex,'.targ.genes.2000.hypo.raw.enhancer.txt'),sep='\t',stringsAsFactors = F,header=F)
  
  bin300.genhancer.hyper.all =  enhancer.id(bin300.genhancer.hyper.all.raw)
  bin300.genhancer.hypo.all =  enhancer.id(bin300.genhancer.hypo.all.raw)
  dmr.list[[Sex]][[paste0('Hypermethylated')]][['Enhancer']] = bin300.genhancer.hyper.all
  dmr.list[[Sex]][[paste0('Hypomethylated')]][['Enhancer']] = bin300.genhancer.hypo.all
  
  bin300.genhancer.hyper.all.raw = read.table(paste0(Sex,'.targ.genes.2000.hyper.raw.enhancerpe.txt'),sep='\t',stringsAsFactors = F,header=F)
  bin300.genhancer.hypo.all.raw = read.table(paste0(Sex,'.targ.genes.2000.hypo.raw.enhancerpe.txt'),sep='\t',stringsAsFactors = F,header=F)
  
  bin300.genhancer.hyper.all =  enhancer.id(bin300.genhancer.hyper.all.raw)
  bin300.genhancer.hypo.all =  enhancer.id(bin300.genhancer.hypo.all.raw)
  dmr.list[[Sex]][[paste0('Hypermethylated')]][['EnhancerBiv']] = bin300.genhancer.hyper.all
  dmr.list[[Sex]][[paste0('Hypomethylated')]][['EnhancerBiv']] = bin300.genhancer.hypo.all
  
  
  bin300.genhancer.hyper.all.raw = read.table(paste0(Sex,'.targ.genes.2000.hyper.raw.promoter.txt'),sep='\t',stringsAsFactors = F,header=F)
  bin300.genhancer.hypo.all.raw = read.table(paste0(Sex,'.targ.genes.2000.hypo.raw.promoter.txt'),sep='\t',stringsAsFactors = F,header=F)
  
  bin300.genhancer.hyper.all =  enhancer.id(bin300.genhancer.hyper.all.raw)
  bin300.genhancer.hypo.all =  enhancer.id(bin300.genhancer.hypo.all.raw)
  dmr.list[[Sex]][[paste0('Hypermethylated')]][['Promoter']] = bin300.genhancer.hyper.all
  dmr.list[[Sex]][[paste0('Hypomethylated')]][['Promoter']] = bin300.genhancer.hypo.all
  
  
  bin300.genhancer.hyper.all.raw = read.table(paste0(Sex,'.targ.genes.2000.hyper.raw.promoterpe.txt'),sep='\t',stringsAsFactors = F,header=F)
  bin300.genhancer.hypo.all.raw = read.table(paste0(Sex,'.targ.genes.2000.hypo.raw.promoterpe.txt'),sep='\t',stringsAsFactors = F,header=F)
  
  bin300.genhancer.hyper.all =  enhancer.id(bin300.genhancer.hyper.all.raw)
  bin300.genhancer.hypo.all =  enhancer.id(bin300.genhancer.hypo.all.raw)
  dmr.list[[Sex]][[paste0('Hypermethylated')]][['PromoterBiv']] = bin300.genhancer.hyper.all
  dmr.list[[Sex]][[paste0('Hypomethylated')]][['PromoterBiv']] = bin300.genhancer.hypo.all
  
  bin300.genhancer.hyper.all.raw = read.table(paste0(Sex,'.targ.genes.2000.hyper.raw.pe.txt'),sep='\t',stringsAsFactors = F,header=F)
  bin300.genhancer.hypo.all.raw = read.table(paste0(Sex,'.targ.genes.2000.hypo.raw.pe.txt'),sep='\t',stringsAsFactors = F,header=F)
  
  bin300.genhancer.hyper.all =  enhancer.id(bin300.genhancer.hyper.all.raw)
  bin300.genhancer.hypo.all =  enhancer.id(bin300.genhancer.hypo.all.raw)
  dmr.list[[Sex]][[paste0('Hypermethylated')]][['Bivalent']] = bin300.genhancer.hyper.all
  dmr.list[[Sex]][[paste0('Hypomethylated')]][['Bivalent']] = bin300.genhancer.hypo.all
  
  
}

targ.markers = names(dmr.list)
hallmark.genes.gsea.all=list()
hallmark.genes.gsea.all.background = list()

for (sex in c('Male','Female')){
  for (dir in c('Hypermethylated','Hypomethylated')){
    annotations = names(dmr.list[[sex]][[dir]])
    for (a in 'All') {#c('Enhancer','Promoter','Bivalent')
      marker = paste0(sex,'_',dir,'_',a)
      

      background =unique(enhancer.background)#)
      
      targ.dmr.genes = dmr.list[[sex]][[dir]][[a]]
      if (length(targ.dmr.genes) > 0) {
        gsea.res = gsea.calc(targ.genes=targ.dmr.genes,
                             targ.background =background,
                             name=marker)
        gsea.res$marker = marker
        hallmark.genes.gsea.all[[marker]] = gsea.res
        
        gsea.res.all = gsea.calc.all(targ.genes=targ.dmr.genes,
                                     targ.background =background,
                                     name=marker)
        gsea.res.all$marker = marker
        gsea.res.all$regulatory.element = a
        hallmark.genes.gsea.all.background[[marker]] = gsea.res.all
        
      } else {
        hallmark.genes.gsea.all[[marker]] = NULL

        hallmark.genes.gsea.all.background[[marker]] = NULL
      }
      
    }
  }
}

#####silencer GSOA####
library(fgsea)
library(msigdbr)
library(biomaRt)
library(Hmisc)
library(tools)
library(ggplot2)
library(gridExtra)
library(ggh4x)
library(ggforce)
setwd('/Users/ncheng/Desktop/ncheng/Desktop/Lab Work/Manuscripts/OHS_cancer_risk/figures/gsoa/current/')
silencer.background.raw = read.table('/Users/ncheng/Desktop/ncheng/Desktop/Lab Work/Manuscripts/OHS_cancer_risk/figures/gsoa//silencer.genes.background.txt',sep='\t',stringsAsFactors = F,header=F)





enhancer.id = function(enhancer.hyper.all){
  mapping.targ <- getBM(
    attributes = c('ensembl_gene_id', 'entrezgene_id', 'hgnc_symbol'), 
    filters = 'hgnc_symbol',
    values = enhancer.hyper.all,
    mart = hsmart
  )
  return(mapping.targ[,1])
}

silencer.dmr = readRDS('silencer.dmr.genes.RDS')#,

dmr.list.silencer = split(silencer.dmr, silencer.dmr[,c('direction','group')])
dmr.list.silencer = lapply(dmr.list.silencer, function(x) x[x$Freq >= 2,'ensgid'])


dmr.list.silencer = list(Male_hyper_silencer = c(dmr.list.silencer[[3]]),
                                    Male_hypo_silencer = c(dmr.list.silencer[[4]]),

                                    Female_hyper_silencer = c(dmr.list.silencer[[1]]),
                                    Female_hypo_silencer = c(dmr.list.silencer[[2]]))

                

targ.markers = names(dmr.list.silencer)
hallmark.genes.gsea.all=list()
hallmark.genes.gsea.all.background = list()
for (marker in targ.markers) {
  if(grepl('.*enhancer', marker) == T) {
    background =unique(enhancer.background)#)
    
  } else{
    background = unique(silencer.background.raw[,1])
    
  }
  background =unique(enhancer.background)#)
  
  
  gsea.res = gsea.calc(targ.genes=dmr.list.silencer[[marker]],
                       targ.background =background,
                       name=marker)
  gsea.res$marker = marker
  gsea.res.all$regulatory.element = 'silencer'
  
  hallmark.genes.gsea.all[[marker]] = gsea.res
  
  gsea.res.all = gsea.calc.all(targ.genes=dmr.list.silencer[[marker]],
                           targ.background =background,
                           name=marker)
  gsea.res.all$marker = marker
  gsea.res.all$regulatory.element = 'silencer'
  
  hallmark.genes.gsea.all.background[[marker]] = gsea.res.all
}



####combining into 1#####
male.hyper = unique(c(dmr.list[['Male']][['Hypermethylated']][['Enhancer']],
                      dmr.list[['Male']][['Hypermethylated']][['Promoter']],
                      dmr.list[['Male']][['Hypermethylated']][['Bivalent']]))
                      #dmr.list.silencer[[ "Top 2000 Hypermethylated Regions.Incident Prostate Cancer"]][,1]))

male.hypo = unique(c(dmr.list[['Male']][['Hypomethylated']][['Enhancer']],
                      dmr.list[['Male']][['Hypomethylated']][['Promoter']],
                      dmr.list[['Male']][['Hypomethylated']][['Bivalent']]))
                      #dmr.list.silencer[[ "Top 2000 Hypomethylated Regions.Incident Prostate Cancer"]][,1]))

female.hyper = unique(c(dmr.list[['Female']][['Hypermethylated']][['Enhancer']],
                      dmr.list[['Female']][['Hypermethylated']][['Promoter']],
                      dmr.list[['Female']][['Hypermethylated']][['Bivalent']]))
                      #dmr.list.silencer[[ "Top 2000 Hypermethylated Regions.Incident Breast Cancer"]][,1]))

male.hypo = unique(c(dmr.list[['Female']][['Hypomethylated']][['Enhancer']],
                      dmr.list[['Female']][['Hypomethylated']][['Promoter']],
                      dmr.list[['Female']][['Hypomethylated']][['Bivalent']]))
                     # dmr.list.silencer[[ "Top 2000 Hypomethylated Regions.Incident Breast Cancer"]][,1]))
dmr.list[['Male']][['Hypermethylated']][['All']] = unique(c(dmr.list[['Male']][['Hypermethylated']][['Enhancer']],
                                                              dmr.list[['Male']][['Hypermethylated']][['Promoter']],
                                                              dmr.list[['Male']][['Hypermethylated']][['Bivalent']]))
dmr.list[['Female']][['Hypermethylated']][['All']] = unique(c(dmr.list[['Female']][['Hypermethylated']][['Enhancer']],
                                                              dmr.list[['Female']][['Hypermethylated']][['Promoter']],
                                                              dmr.list[['Female']][['Hypermethylated']][['Bivalent']]))
dmr.list[['Male']][['Hypomethylated']][['All']] = unique(c(dmr.list[['Male']][['Hypomethylated']][['Enhancer']],
                                                            dmr.list[['Male']][['Hypomethylated']][['Promoter']],
                                                            dmr.list[['Male']][['Hypomethylated']][['Bivalent']]))
dmr.list[['Female']][['Hypomethylated']][['All']] = unique(c(dmr.list[['Female']][['Hypomethylated']][['Enhancer']],
                                                              dmr.list[['Female']][['Hypomethylated']][['Promoter']],
                                                              dmr.list[['Female']][['Hypomethylated']][['Bivalent']]))


#####ploting local GOSA#####
hallmark.genes.gsea.silencer.sig = lapply(hallmark.genes.gsea.all, function(x) x[x$padj < 0.1,])
combined.pathways = NULL
hallmark.genes.gsea.silencer.sig.filt = lapply(hallmark.genes.gsea.silencer.sig, function(x) x[x$gene.set %in% c('C2.KEGG','C2.REACTOME','H'),])

hallmark.genes.all = do.call('rbind',hallmark.genes.gsea.silencer.sig.filt)

targ.plot = hallmark.genes.all[order(hallmark.genes.all[,c('gene.set','pval')]),]
targ.plot = targ.plot[order(targ.plot$padj),]
targ.plot$pathway = gsub('GOBP_|GOMF_|HALLMARK_|WP_|REACTOME_|KEGG_','',targ.plot$pathway)
#targ.plot$pathway = toTitleCase(tolower(gsub('_',' ',targ.plot$pathway)))
targ.plot$pathway = gsub('_',' ',targ.plot$pathway)
targ.plot = data.frame(targ.plot,check.names=F)
write.table(targ.plot$pathway, paste0('combined.regulatory.pathways1.txt'), sep='\t',row.names=F,col.names=F,quote=F)
#write.table(targ.plot$pathway, paste0('combined.silencer.pathways.txt'), sep='\t',row.names=F,col.names=F,quote=F)


pathway.groups = read.csv('biological_pathway_groups.csv',header=T,stringsAsFactors = F)
#pathway.groups = read.csv('silencer.pathway.groups.csv',header=T,stringsAsFactors = F)
#pathway.groups = read.csv('grouped_pathways_new.csv',header=T,stringsAsFactors = F)

targ.plot.grouped= merge(targ.plot, pathway.groups,by='pathway')
db.colors =rev(c('C5.BP'='#107E7D','C5.MF'='#044B7F','C2.REACTOME'='#FF8C42','H' ='#C5283D','C2.KEGG'='#330C2F'))
targ.plot.grouped$regulatory.element =gsub('.*_','',targ.plot.grouped$marker)
exc = c('Cardiovascular Smooth Muscle','Nervous System Function','Gene Regulation & Epigenetics','Metabolism & Organelle Function')
exc = c('Cardiac / Muscle','General Signaling','Vesicle Transport','Gene Expression','Cytoskeleton / ECM','Cell Structure & Adhesion','Gene Regulation & Transcription','Other')
targ.plot.grouped = targ.plot.grouped[!targ.plot.grouped$biological.process %in% exc,]

for (reg in unique(targ.plot.grouped$regulatory.element)) {
  targ.plot.grouped.plot = unique(targ.plot.grouped[targ.plot.grouped$regulatory.element == reg,])
  targ.plot.grouped.plot$plot.group = gsub('Female_','Breast\nCancer\n',gsub('Male_','Prostate\nCancer\n',gsub('hypo_.*','Hypomethylated',gsub('hyper.*','Hypermethylated',targ.plot.grouped.plot$marker))))
  
  hallmark.genes.gsea.all.background.df= do.call('rbind',hallmark.genes.gsea.all.background)
  hallmark.genes.gsea.all.background.df.filt = hallmark.genes.gsea.all.background.df[hallmark.genes.gsea.all.background.df$regulatory.element == reg,]
  hallmark.genes.gsea.all.background.df.filt = hallmark.genes.gsea.all.background.df.filt[hallmark.genes.gsea.all.background.df.filt$pathway != 'REACTOME_MYOGENESIS',]
  geom.count.fill =function(targ.plot.grouped.plot, hallmark.genes.gsea.all.background.df.filt, pathway.groups) {
    library(tools)
    background.pathways = hallmark.genes.gsea.all.background.df.filt
    background.pathways$pathway = gsub('GOBP_|GOMF_|HALLMARK_|WP_|REACTOME_|KEGG_','',background.pathways$pathway)
    background.pathways$pathway = gsub('_',' ',background.pathways$pathway)
    filtered.background.pathways = background.pathways[background.pathways$pathway %in% targ.plot.grouped.plot$pathway,]
    filtered.background.pathways = data.frame(filtered.background.pathways[,c('pathway','padj','overlap','marker')],stringsAsFactors = F)
    filtered.background.pathways.grouped= merge(filtered.background.pathways, pathway.groups,by='pathway')
    filtered.background.pathways.grouped$regulatory.element =gsub('.*_','',filtered.background.pathways.grouped$marker)
    filtered.background.pathways.grouped$plot.group = gsub('_.*','',gsub('Female_','Breast\nCancer\n',gsub('Male_','Prostate\nCancer\n',gsub('hypo_.*','Hypomethylated',gsub('hyper.*','Hypermethylated',filtered.background.pathways.grouped$marker)))))
    filtered.background.pathways.grouped$padj = ifelse(filtered.background.pathways.grouped$padj > 0.1, 1, filtered.background.pathways.grouped$padj)
    filtered.background.pathways.grouped$overlap = ifelse(filtered.background.pathways.grouped$overlap ==0, 1, filtered.background.pathways.grouped$overlap)
    filtered.background.pathways.grouped$pathway = toTitleCase(tolower(as.character(filtered.background.pathways.grouped$pathway)))
    filtered.background.pathways.grouped = filtered.background.pathways.grouped[order(filtered.background.pathways.grouped$padj),]
    filtered.background.pathways.grouped = filtered.background.pathways.grouped[!duplicated(filtered.background.pathways.grouped[,c('pathway','marker')] ),]
    return(filtered.background.pathways.grouped)
  }
  
  gsoa.plot = geom.count.fill(targ.plot.grouped.plot, hallmark.genes.gsea.all.background.df.filt, pathway.groups)
  gsoa.plot = gsoa.plot[gsoa.plot$regulatory.element == reg,]
  gsoa.plot$pathway = gsub(' Foci Sahf','',gsoa.plot$pathway)
  gsoa.plot$pathway = gsub('.*Fcgamma Receptor Fcgr Dependent Phagocytosis.*','Fc Gamma r Mediated Phagocytosis',gsoa.plot$pathway)
  gsoa.plot = gsoa.plot[order(gsoa.plot$padj),]
  gsoa.plot = gsoa.plot[!duplicated(gsoa.plot[,c('pathway','marker')]),]
  vals = c(0, -log10(0.05),-log10(0.001))
  range = function(x){(x-min(x))/(max(x)-min(x))}
  vals = scales::rescale(vals, to = c(0, 1))
  gsoa.plot$padj =ifelse(gsoa.plot$padj <= 0.001,0.001,gsoa.plot$padj )
  gsoa.plot$plot.group = factor(as.character(gsoa.plot$plot.group), levels =c('Prostate\nCancer\nHypomethylated','Prostate\nCancer\nHypermethylated','Breast\nCancer\nHypomethylated','Breast\nCancer\nHypermethylated') )
  plot1 = ggplot(gsoa.plot[!grepl('Other',gsoa.plot$biological.process),], aes(y=pathway,x = plot.group)) +
    geom_count(aes(fill=-log10(padj),size=overlap),colour="black",pch=21,alpha=0.7)+
    
    #facet_grid2(. ~ biological.process ,scales = 'free',remove_labels=T, space = "free",ncol=1)+
    facet_col(facets = vars(biological.process), 
              scales = "free_y", 
              space = "free") +
    
    # facet_grid( biological.process ~ . ,scales = 'free_y', space='free')+
    scale_y_discrete(position = "right")+
    #scale_shape_manual(values=c(23))+
    theme_bw()+
    theme(text = element_text(size=10,face='bold'),
          axis.ticks.y = element_blank(),
          axis.ticks.x = element_blank(),
          legend.position = 'right',
          legend.title = element_blank(),
          strip.background = element_rect(fill="white")#,
          #strip.background = element_blank(),
          #strip.text.x = element_blank()
          #axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)
          
    ) +#"#e7f0fa",#"#c9e2f6",95cbee
    scale_fill_gradientn(colors = c('#FFF6F1', "#B0517D", "#482066"), #F5B0C3
                         limits = c(0,max(-log10(0.001))),
                         values = vals)+
    scale_colour_gradientn(colors = c('#FFF6F1', "#B0517D", "#482066"),
                           limits = c(0,max(-log10(0.001))),
                           values = vals)+
    
    xlab('') + ylab('Biological Pathway or Process') # + scale_fill_manual(values =db.colors )#+ ggtitle(title) +
  
  print(plot1)
  
  
  png(paste0(reg,'.combined.enhancer.gsoa.genes.png'),height=3000,width=3500,res=400)
  print(plot1)
  dev.off()  
}



#plotting multiple elements in one plot for main figure
targ.elements = c('Enhancer','Bivalent')
targ.elements = 'Bivalent'
targ.plot.grouped.plot = unique(targ.plot.grouped[targ.plot.grouped$regulatory.element %in% targ.elements,])
targ.plot.grouped.plot$plot.group = gsub('Female_','Breast\nCancer\n',gsub('Male_','Prostate\nCancer\n',gsub('hypo_.*','Hypomethylated',gsub('hyper.*','Hypermethylated',targ.plot.grouped.plot$marker))))
hallmark.genes.gsea.all.background.df= do.call('rbind',hallmark.genes.gsea.all.background)
hallmark.genes.gsea.all.background.df.filt = hallmark.genes.gsea.all.background.df[hallmark.genes.gsea.all.background.df$regulatory.element %in% targ.elements,]
hallmark.genes.gsea.all.background.df.filt = hallmark.genes.gsea.all.background.df.filt[hallmark.genes.gsea.all.background.df.filt$pathway != 'REACTOME_MYOGENESIS',]
geom.count.fill =function(targ.plot.grouped.plot, hallmark.genes.gsea.all.background.df.filt, pathway.groups) {
  library(tools)
  background.pathways = hallmark.genes.gsea.all.background.df.filt
  background.pathways$pathway = gsub('GOBP_|GOMF_|HALLMARK_|WP_|REACTOME_|KEGG_','',background.pathways$pathway)
  background.pathways$pathway = gsub('_',' ',background.pathways$pathway)
  filtered.background.pathways = background.pathways[background.pathways$pathway %in% targ.plot.grouped.plot$pathway,]
  filtered.background.pathways = data.frame(filtered.background.pathways[,c('pathway','padj','overlap','marker','regulatory.element')],stringsAsFactors = F)
  filtered.background.pathways.grouped= merge(filtered.background.pathways, pathway.groups,by='pathway')
  filtered.background.pathways.grouped$regulatory.element =gsub('.*_','',filtered.background.pathways.grouped$marker)
  filtered.background.pathways.grouped$plot.group = gsub('_.*','',gsub('Female_','Breast\nCancer\n',gsub('Male_','Prostate\nCancer\n',gsub('hypo_.*','Hypomethylated',gsub('hyper.*','Hypermethylated',filtered.background.pathways.grouped$marker)))))
  filtered.background.pathways.grouped$plot.group = paste0(filtered.background.pathways.grouped$plot.group, ' ', filtered.background.pathways.grouped$regulatory.element)
  filtered.background.pathways.grouped$padj = ifelse(filtered.background.pathways.grouped$padj > 0.1, 1, filtered.background.pathways.grouped$padj)
  filtered.background.pathways.grouped$overlap = ifelse(filtered.background.pathways.grouped$overlap ==0, 1, filtered.background.pathways.grouped$overlap)
  filtered.background.pathways.grouped$pathway = toTitleCase(tolower(as.character(filtered.background.pathways.grouped$pathway)))
  filtered.background.pathways.grouped = filtered.background.pathways.grouped[order(filtered.background.pathways.grouped$padj),]
  filtered.background.pathways.grouped = filtered.background.pathways.grouped[!duplicated(filtered.background.pathways.grouped[,c('pathway','marker')] ),]
  return(filtered.background.pathways.grouped)
}

gsoa.plot = geom.count.fill(targ.plot.grouped.plot, hallmark.genes.gsea.all.background.df.filt, pathway.groups)
gsoa.plot = gsoa.plot[gsoa.plot$regulatory.element %in% targ.elements,]
plot1 = ggplot(gsoa.plot, aes(y=pathway,x = plot.group)) +
  geom_count(aes(fill=-log10(padj),size=overlap),colour="black",pch=21,alpha=0.7)+
  
  #facet_grid2(. ~ biological.process ,scales = 'free',remove_labels=T, space = "free",ncol=1)+
  facet_col(facets = vars(biological.process), 
            scales = "free_y", 
            space = "free") +
  
  # facet_grid( biological.process ~ . ,scales = 'free_y', space='free')+
  scale_y_discrete(position = "right")+
  #scale_shape_manual(values=c(23))+
  theme_bw()+
  theme(text = element_text(size=11,face='bold'),
        axis.ticks.y = element_blank(),
        axis.ticks.x = element_blank(),
        legend.position = 'right',
        legend.title = element_blank(),
        strip.background = element_rect(fill="white")#,
        #strip.background = element_blank(),
        #strip.text.x = element_blank()
        #axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)
        
  ) +#"#e7f0fa",#"#c9e2f6",95cbee
  scale_fill_gradientn(colors = c('#FFF6F1', "#B0517D", "#482066"), #F5B0C3
                       limits = c(0,max(-log10(gsoa.plot$padj))))+
  scale_colour_gradientn(colors = c('#FFF6F1', "#B0517D", "#482066"),
                         limits = c(0,max(-log10(gsoa.plot$padj))))+
  
  xlab('') + ylab('Biological Pathway or Process') # + scale_fill_manual(values =db.colors )#+ ggtitle(title) +

print(plot1)


png(paste0('mainfig2.combined.biv.gsoa.genes.png'),height=3000,width=3200,res=400)
print(plot1)
dev.off()  

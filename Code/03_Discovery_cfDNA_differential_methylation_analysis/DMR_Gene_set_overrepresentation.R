#####Annotating DMRs with regulatory element gene targets#####
#loading packages
library(GenomicRanges)
library(stringr)
#reading in enhancer promoter annotation from genhancer
genhancer = readRDS('genhancer.grange.unreduced.RDS')
genhancer.annotation.genes = readRDS('genhancer.annotation.genes.RDS')


savedir='/save/directory/'

#annotating enhancers/promoters and gene targets among top 2000 hypermethylated and hypomethylated regions
for (sex in c('Male','Female')){
  #selecting top 2000 hypermethylated and hypomethylated regions
  res.df = readRDS(paste0('discovery.',sex,'.dmrs.RDS'))
  res.df = res.df[!grepl('chrX|chrY',res.df$window),]
  res.df = res.df[res.df$window %in% cpg_count[cpg_count$count>=5,'window'],]
  res.df.hyper.sig.raw = res.df[res.df$log2FoldChange > 0.25 & res.df$pvalue < 0.05,]
  res.df.hypo.sig.raw = res.df[res.df$log2FoldChange < -0.25 & res.df$pvalue < 0.05,]
  res.df.hyper.sig.raw = res.df.hyper.sig.raw[order(res.df.hyper.sig.raw$pvalue),]
  res.df.hypo.sig.raw = res.df.hypo.sig.raw[order(res.df.hypo.sig.raw$pvalue),]
  
  #annotating enhancers/promoters and corresponding gene targets
  all.features.sig.hypo = res.df.hypo.sig.raw
  all.features.sig.hyper = res.df.hyper.sig.raw
  all.features.sig.hyper = all.features.sig.hyper[order(all.features.sig.hyper$pvalue),]
  all.features.sig.hypo = all.features.sig.hypo[order(all.features.sig.hypo$pvalue),]
  
  hyper.j = min(nrow(all.features.sig.hyper),2000)
  hypo.j=min(nrow(all.features.sig.hypo),2000)
  targ.windows.hyper =all.features.sig.hyper$window[1:hyper.j]
  targ.windows.hypo =all.features.sig.hypo$window[1:hypo.j]
  
  enhancer.promoter.gene.annotation= function(targ.windows,genhancer) {
    #convert targeted windows to grange
    tmp = data.frame(chr = gsub(':.*','',targ.windows),
                     start = gsub('-.*','',gsub('.*:','',targ.windows)),
                     end = gsub('.*-','',targ.windows),
                     window = targ.windows)
    
    colnames(tmp)[c(1:3)] = c('chrom','chromStart','chromEnd')
    
    targ.grange <- makeGRangesFromDataFrame(tmp,  keep.extra.columns=T)
    
    genhancer.filt =data.frame(genhancer[queryHits(findOverlaps(genhancer, targ.grange))])   
    targ.enhancer.ids = genhancer.filt$GHid
    
    genes= genhancer.annotation.genes[genhancer.annotation.genes$GHid %in% targ.enhancer.ids,]
    return(genes) 
  }
  
  #Enhancers
  enhancer.genes.hyper = enhancer.promoter.gene.annotation(targ.windows.hyper,genhancer[genhancer$regulatory_element_type %in% c('Enhancer'),] )
  write.table(enhancer.genes.hyper$symbol,paste0(savedir,sex,'.targ.genes.2000.hyper.raw.enhancer.txt'),sep='\t',row.names=F,quote=F,col.names=F)
  
  enhancer.genes.hypo = enhancer.promoter.gene.annotation(targ.windows.hypo,genhancer[genhancer$regulatory_element_type %in% c('Enhancer'),] )
  write.table(enhancer.genes.hypo$symbol,paste0(savedir,sex,'.targ.genes.2000.hypo.raw.enhancer.txt'),sep='\t',row.names=F,quote=F,col.names=F)
  
  #Dual promoter/enhancers
  enhancer.genes.hyper = enhancer.promoter.gene.annotation(targ.windows.hyper,genhancer[genhancer$regulatory_element_type %in% c('Promoter/Enhancer'),] )
  write.table(enhancer.genes.hyper$symbol,paste0(savedir,sex,'.targ.genes.2000.hyper.raw.pe.txt'),sep='\t',row.names=F,quote=F,col.names=F)
  
  enhancer.genes.hypo = enhancer.promoter.gene.annotation(targ.windows.hypo,genhancer[genhancer$regulatory_element_type %in% c('Promoter/Enhancer'),] )
  write.table(enhancer.genes.hypo$symbol,paste0(savedir,sex,'.targ.genes.2000.hypo.raw.pe.txt'),sep='\t',row.names=F,quote=F,col.names=F)
  
  #Promoters
  enhancer.genes.hyper = enhancer.promoter.gene.annotation(targ.windows.hyper,genhancer[genhancer$regulatory_element_type %in% c('Promoter'),] )
  write.table(enhancer.genes.hyper$symbol,paste0(savedir,sex,'.targ.genes.2000.hyper.raw.promoter.txt'),sep='\t',row.names=F,quote=F,col.names=F)
  
  enhancer.genes.hypo = enhancer.promoter.gene.annotation(targ.windows.hypo,genhancer[genhancer$regulatory_element_type %in% c('Promoter'),] )
  write.table(enhancer.genes.hypo$symbol,paste0(savedir,sex,'.targ.genes.2000.hypo.raw.promoter.txt'),sep='\t',row.names=F,quote=F,col.names=F)
  
  
}


#annotating silencers and gene targets among top 2000 hypermethylated and hypomethylated regions
gene.coords = readRDS('gencode.hg38.gene.coordinates.RDS')

#promoter HI-C data set 1 (can be downloaded here https://www.ebi.ac.uk/biostudies/arrayexpress/studies/E-MTAB-2323)
promoter.hic1 = read.table('TS5_CD34_promoter-other_significant_interactions.txt',sep='\t',stringsAsFactors = F,header=T)
promoter.hic2 = read.table('TS5_GM12878_promoter-other_significant_interactions.txt',sep='\t',stringsAsFactors = F,header=T)

#function to extract gene targets of silencers from dataset 1
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


#promoter HI-C data set 2 (can be downloaded here https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE86821)
promoter.hic3 = read.table('GSM2309025_hESC_Sample811_Lane912_washU.txt.gz',sep='\t',stringsAsFactors = F,header=F)
promoter.hic4 = read.table('GSM2309026_hESC_Samples1096_1390_Lanes1335_1785_washU.txt.gz',sep='\t',stringsAsFactors = F,header=F)
promoter.hic5 = read.table('GSM2309027_hNPC_Sample1145_Lane1452_washU.txt.gz',sep='\t',stringsAsFactors = F,header=F)
promoter.hic6 = read.table('GSM2309028_hNPC_Sample1742_Lane2174_washU.txt.gz',sep='\t',stringsAsFactors = F,header=F)

#function to extract gene targets of silencers from dataset 2
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


#promoter HI-C data set 3 (can be downloaded here https://www.encodeproject.org/experiments/ENCSR549MGQ/)
promoter.hic7=read.table('ENCFF693XIL.bedpe.gz',sep='\t',stringsAsFactors = F,header=F)

#function to extract gene targets of silencers from dataset 3
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

#promoter HI-C data set 4 (can be downloaded here https://www.encodeproject.org/files/ENCFF202FID/ and here https://www.encodeproject.org/files/ENCFF128MZW/)
promoter.hic8=read.table('ENCFF202FID.bed.gz',sep='\t',stringsAsFactors = F,header=F)
promoter.hic9=read.table('ENCFF128MZW.bed.gz',sep='\t',stringsAsFactors = F,header=F)

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


#annotating silencer gene targets from DMRs
combined.silencer=NULL
for (sex in c('Male','Female')) {
  #selecting top 2000 hypermethylated and hypometyhlated regions
  res.df = readRDS(paste0(savedir,'discovery.',sex,'.dmrs.RDS'))
  res.df = res.df[!grepl('chrX|chrY',res.df$window),]
  res.df.hyper.sig.raw = res.df[res.df$log2FoldChange > 0.25 & res.df$pvalue < 0.05,]
  res.df.hypo.sig.raw = res.df[res.df$log2FoldChange < -0.25 & res.df$pvalue < 0.05,]
  silencer.regions.hyper = res.df.hyper.sig.raw[order(res.df.hyper.sig.raw$pvalue),][1:2000,'window']
  silencer.regions.hypo = res.df.hypo.sig.raw[order(res.df.hypo.sig.raw$pvalue),][1:2000,'window']
  
  
  #annotating silencers gene target using Hi-C data to identify promoter contacts 
  #hypermethylated regions
  silencer.hic1.hyper =silencer.gene.extract.v1(silencer.regions.hyper,promoter.hic1)
  silencer.hic2.hyper =silencer.gene.extract.v1(silencer.regions.hyper,promoter.hic2)
  silencer.hic3.hyper = silencer.gene.extract.v2(silencer.regions.hyper, promoter.hic3,gene.coords )
  silencer.hic4.hyper = silencer.gene.extract.v2(silencer.regions.hyper, promoter.hic4,gene.coords )
  silencer.hic5.hyper = silencer.gene.extract.v2(silencer.regions.hyper, promoter.hic5,gene.coords )
  silencer.hic6.hyper = silencer.gene.extract.v2(silencer.regions.hyper, promoter.hic6,gene.coords )
  silencer.hic7.hyper = silencer.gene.extract.v3(silencer.regions.hyper, promoter.hic7,gene.coords )
  silencer.hic8.hyper = silencer.gene.extract.v4(silencer.regions.hyper, promoter.hic8,gene.coords )
  silencer.hic9.hyper = silencer.gene.extract.v4(silencer.regions.hyper, promoter.hic9,gene.coords )
  
  
  #hypomethylated regions
  silencer.hic1.hypo =silencer.gene.extract.v1(silencer.regions.hypo,promoter.hic1)
  silencer.hic2.hypo =silencer.gene.extract.v1(silencer.regions.hypo,promoter.hic2)
  silencer.hic3.hypo = silencer.gene.extract.v2(silencer.regions.hypo, promoter.hic3,gene.coords )
  silencer.hic4.hypo = silencer.gene.extract.v2(silencer.regions.hypo, promoter.hic4,gene.coords )
  silencer.hic5.hypo = silencer.gene.extract.v2(silencer.regions.hypo, promoter.hic5,gene.coords )
  silencer.hic6.hypo = silencer.gene.extract.v2(silencer.regions.hypo, promoter.hic6,gene.coords )
  silencer.hic7.hypo = silencer.gene.extract.v3(silencer.regions.hypo, promoter.hic7,gene.coords )
  silencer.hic8.hypo = silencer.gene.extract.v4(silencer.regions.hypo, promoter.hic8,gene.coords )
  silencer.hic9.hypo = silencer.gene.extract.v4(silencer.regions.hypo, promoter.hic9,gene.coords )
  
  #combining hypermethylated into gene list
  all.hyper = c(unique(silencer.hic1.hyper$Ensembl.Gene.ID), 
                unique(silencer.hic2.hyper$Ensembl.Gene.ID),
                unique(silencer.hic3.hyper$ensgid),
                unique(silencer.hic4.hyper$ensgid),
                unique(silencer.hic5.hyper$ensgid),
                unique(silencer.hic6.hyper$ensgid),
                unique(silencer.hic7.hyper$ensgid),
                unique(silencer.hic8.hyper$ensgid),
                unique(silencer.hic9.hyper$ensgid))
  
  #combining hypomethylated into gene list
  all.hypo = c(unique(silencer.hic1.hypo$Ensembl.Gene.ID), 
               unique(silencer.hic2.hypo$Ensembl.Gene.ID),
               unique(silencer.hic3.hypo$ensgid),
               unique(silencer.hic4.hypo$ensgid),
               unique(silencer.hic5.hypo$ensgid),
               unique(silencer.hic6.hypo$ensgid),
               unique(silencer.hic7.hyper$ensgid),
               unique(silencer.hic8.hyper$ensgid),
               unique(silencer.hic9.hyper$ensgid))
  
  #creating dataframe of results
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


####Performing gene set overrepresentation analysis######
#loading libraries
library(fgsoa)
library(msigdbr)
library(biomaRt)
library(Hmisc)
library(tools)
library(ggplot2)
library(gridExtra)
library(ggh4x)
library(ggforce)

#gene set overrepresentation with enhancer/promoter gene targets 
enhancer.background=readRDS('enhancer.background.ensgid.RDS')

#converting ensgid to entrezid
hsmart <- useEnsembl(dataset = "hsapiens_gene_ensembl", biomart = "ensembl")

#functions for performing GSOA
enhancer.id = function(enhancer.hyper.all){
  mapping.targ <- getBM(
    attributes = c('ensembl_gene_id', 'entrezgene_id', 'hgnc_symbol'), 
    filters = 'hgnc_symbol',
    values = enhancer.hyper.all,
    mart = hsmart
  )
  return(mapping.targ[,1])
}
gsoa.calc = function(targ.genes, targ.background,name) {
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

  
  gene.set.list = list(H = msigdbr_df.h,
                       C2.KEGG = msigdbr_df.c2.kegg,
                       C2.WIKIPATH = msigdbr_df.c2.wikipath,
                       C2.REACTOME = msigdbr_df.c2.reactome)

  
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
  # fixing format to work with fgsoa
  #reactome
  
  
  
  return(return.df)
  
}
gsoa.calc.all = function(targ.genes, targ.background,name) {
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

  
  gene.set.list = list(H = msigdbr_df.h,
                       C2.KEGG = msigdbr_df.c2.kegg,
                       C2.REACTOME = msigdbr_df.c2.reactome#,

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
  # fixing format to work with fgsoa
  #reactome
  
  
  
  return(return.df)
  
}

#reading in promoter/enhancer gene target files from earlier into a list
dmr.list = list( )
for (Sex in c('Male','Female')) {
  dmr.list[[Sex]] = list()
  dmr.list[[Sex]][[paste0('Hypermethylated')]] = list()
  dmr.list[[Sex]][[paste0('Hypomethylated')]] = list()

  bin300.genhancer.hyper.all =  enhancer.id(bin300.genhancer.hyper.all.raw)
  bin300.genhancer.hypo.all =  enhancer.id(bin300.genhancer.hypo.all.raw)
  dmr.list[[Sex]][[paste0('Hypermethylated')]][['Enhancer']] = bin300.genhancer.hyper.all
  dmr.list[[Sex]][[paste0('Hypomethylated')]][['Enhancer']] = bin300.genhancer.hypo.all
  
  bin300.genhancer.hyper.all.raw = read.table(paste0(Sex,'.targ.genes.2000.hyper.raw.promoter.txt'),sep='\t',stringsAsFactors = F,header=F)
  bin300.genhancer.hypo.all.raw = read.table(paste0(Sex,'.targ.genes.2000.hypo.raw.promoter.txt'),sep='\t',stringsAsFactors = F,header=F)
  
  bin300.genhancer.hyper.all =  enhancer.id(bin300.genhancer.hyper.all.raw)
  bin300.genhancer.hypo.all =  enhancer.id(bin300.genhancer.hypo.all.raw)
  dmr.list[[Sex]][[paste0('Hypermethylated')]][['Promoter']] = bin300.genhancer.hyper.all
  dmr.list[[Sex]][[paste0('Hypomethylated')]][['Promoter']] = bin300.genhancer.hypo.all
  
  
  bin300.genhancer.hyper.all.raw = read.table(paste0(Sex,'.targ.genes.2000.hyper.raw.pe.txt'),sep='\t',stringsAsFactors = F,header=F)
  bin300.genhancer.hypo.all.raw = read.table(paste0(Sex,'.targ.genes.2000.hypo.raw.pe.txt'),sep='\t',stringsAsFactors = F,header=F)
  
  bin300.genhancer.hyper.all =  enhancer.id(bin300.genhancer.hyper.all.raw)
  bin300.genhancer.hypo.all =  enhancer.id(bin300.genhancer.hypo.all.raw)
  dmr.list[[Sex]][[paste0('Hypermethylated')]][['Bivalent']] = bin300.genhancer.hyper.all
  dmr.list[[Sex]][[paste0('Hypomethylated')]][['Bivalent']] = bin300.genhancer.hypo.all
  
  
}

targ.markers = names(dmr.list)
hallmark.genes.gsoa.all=list()
hallmark.genes.gsoa.all.background = list()

#peforming GSOA across enhancer/promoter gene targets
for (sex in c('Male','Female')){
  for (dir in c('Hypermethylated','Hypomethylated')){
    annotations = names(dmr.list[[sex]][[dir]])
    for (a in 'All') {
      marker = paste0(sex,'_',dir,'_',a)
      
      
      background =unique(enhancer.background)#)
      
      targ.dmr.genes = dmr.list[[sex]][[dir]][[a]]
      if (length(targ.dmr.genes) > 0) {
        #GSOA
        gsoa.res = gsoa.calc(targ.genes=targ.dmr.genes,
                             targ.background =background,
                             name=marker)
        gsoa.res$marker = marker
        hallmark.genes.gsoa.all[[marker]] = gsoa.res
        #GSOA including non-significant results
        gsoa.res.all = gsoa.calc.all(targ.genes=targ.dmr.genes,
                                     targ.background =background,
                                     name=marker)
        gsoa.res.all$marker = marker
        gsoa.res.all$regulatory.element = a
        hallmark.genes.gsoa.all.background[[marker]] = gsoa.res.all
        
      } else {
        hallmark.genes.gsoa.all[[marker]] = NULL
        hallmark.genes.gsoa.all.background[[marker]] = NULL
      }
      
    }
  }
}

hallmark.genes.gsoa.genhancer.sig = lapply(hallmark.genes.gsoa.all, function(x) x[x$padj < 0.1,])
combined.pathways = NULL

hallmark.genes.all = do.call('rbind',hallmark.genes.gsoa.silencer.sig.filt)
targ.plot = hallmark.genes.all[order(hallmark.genes.all[,c('gene.set','pval')]),]
targ.plot = targ.plot[order(targ.plot$padj),]
targ.plot$pathway = gsub('GOBP_|GOMF_|HALLMARK_|WP_|REACTOME_|KEGG_','',targ.plot$pathway)
targ.plot$pathway = gsub('_',' ',targ.plot$pathway)
targ.plot = data.frame(targ.plot,check.names=F)
saveRDS(targ.plot,'enriched.genhancer.genesets.RDS')



#peforming GSOA across silencer gene targets
silencer.background.raw = read.table('silencer.genes.background.txt',sep='\t',stringsAsFactors = F,header=F)
silencer.dmr = readRDS('silencer.dmr.genes.RDS')


dmr.list.silencer = split(silencer.dmr, silencer.dmr[,c('direction','group')])
#selecting silencers supported by 2 or more datasets
dmr.list.silencer = lapply(dmr.list.silencer, function(x) x[x$Freq >= 2,'ensgid'])

dmr.list.silencer = list(Male_hyper_silencer = c(dmr.list.silencer[[3]]),
                         Male_hypo_silencer = c(dmr.list.silencer[[4]]),
                         Female_hyper_silencer = c(dmr.list.silencer[[1]]),
                         Female_hypo_silencer = c(dmr.list.silencer[[2]]))


targ.markers = names(dmr.list.silencer)
hallmark.genes.gsoa.all=list()
hallmark.genes.gsoa.all.background = list()
for (marker in targ.markers) {
  background =unique(enhancer.background)
  #GSOA 
  gsoa.res = gsoa.calc(targ.genes=dmr.list.silencer[[marker]],
                       targ.background =background,
                       name=marker)
  gsoa.res$marker = marker
  gsoa.res.all$regulatory.element = 'silencer'
  
  hallmark.genes.gsoa.all[[marker]] = gsoa.res
  #GSOA including non-significant results
  gsoa.res.all = gsoa.calc.all(targ.genes=dmr.list.silencer[[marker]],
                               targ.background =background,
                               name=marker)
  gsoa.res.all$marker = marker
  gsoa.res.all$regulatory.element = 'silencer'
  
  hallmark.genes.gsoa.all.background[[marker]] = gsoa.res.all
}

hallmark.genes.gsoa.silencer.sig = lapply(hallmark.genes.gsoa.all, function(x) x[x$padj < 0.1,])
combined.pathways = NULL

hallmark.genes.all = do.call('rbind',hallmark.genes.gsoa.silencer.sig.filt)
targ.plot = hallmark.genes.all[order(hallmark.genes.all[,c('gene.set','pval')]),]
targ.plot = targ.plot[order(targ.plot$padj),]
targ.plot$pathway = gsub('GOBP_|GOMF_|HALLMARK_|WP_|REACTOME_|KEGG_','',targ.plot$pathway)
targ.plot$pathway = gsub('_',' ',targ.plot$pathway)
targ.plot = data.frame(targ.plot,check.names=F)
saveRDS(targ.plot,'enriched.silencer.genesets.RDS')

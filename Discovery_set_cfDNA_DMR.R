#loading libraries
library(DESeq2)
library(BiocParallel)
ncores = 10
register(MulticoreParam(ncores))
options(expressions = 5e5)

#reading in sample information file
combined.set = readRDS('combined.set.samples.RDS')
discovery.set = combined.set[combined.set$data.partition == 'Discovery', ]

#reading in raw couunts (deposited in zenodo)
raw.counts=  files.df = readRDS('Tiled.300bp.raw.counts.RDS')

#####setting savedir####
savedir='/savedirectory/'
figdir=paste0(savedir,'figures/')

sex.list = list('Female' = discovery.set[discovery.set$Sex == 'Female',],
                   'Male'  = discovery.set[discovery.set$Sex == 'Male',])


#performing differential methylation analysis between cases and controls using DESeq2
for (sex in names(sex.list)) {
  #target samples
  combined.sample.info.filt =sex.list[[sex]]

  #filtering raw couints for sample of interest
  combined.counts = files.df[combined.sample.info.filt$GRP_Id,!colnames(files.df) %in% 'window']

  #creating deseq object from raw counts
  dds <- DESeqDataSetFromMatrix(countData = combined.counts,
                                colData = combined.sample.info.filt,
                                design= ~  group ) 
  
  #setting factor levels
  colData(dds)$filler = factor(colData(dds)$filler,levels= c('MFiller','UFiller'))
  colData(dds)$group = factor(ifelse(colData(dds)$Cancer =='Control','Control','Cancer'),levels = c('Control','Cancer'))
  dds$condition = dds$group
  
  #estimating size factors
  dds <- estimateSizeFactors(dds) #estimating size factors for tmm normalization
  
  #defining DMR calling model
  mm = model.matrix(~ filler + group, colData(dds)) 

  #differential methylation analysis
  ddssva <- DESeq(dds,full = mm,parallel=T)
  
  #generating results table
  res.df = results(ddssva,contrast = list('groupCancer'))
  res.df$window = rownames(res.df)
  res.df = data.frame(res.df, stringsAsFactors = F)
  
  #saving results
  saveRDS(res.df,paste0(savedir,'discovery.',sex,'.dmrs.RDS'))
  
  
  
}

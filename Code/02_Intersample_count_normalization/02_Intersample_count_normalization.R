#loading libraries
library(parallel)
library(DESeq2)


#loading sample information file
sample.info.filt=readRDS('combined.set.samples.RDS')

#setting directories
savedir='/User/workdirectory/'
rawcountdir='/path/to/raw/counts/'

#reading in genomic annotation files for filtering#
blacklist = readRDS('hg38-blacklist.v2_window300.RDS')

#separating discovery and test samples
sample.list = list(Male = sample.info.filt[sample.info.filt$Sex == 'Male' & sample.info.filt$Cancer %in% c('Control','Prostate'),],
                   Female = sample.info.filt[sample.info.filt$Sex == 'Female' & sample.info.filt$Cancer %in% c('Control','Breast'),],
                   All = sample.info.filt)

library(DESeq2)
#list of raw count files
setwd(rawcountdir)


#discovery only (or not)
discovery.only = T

####loading + normalizing counts for 300bp bin matrix####
for (sex in names(sample.list))  {
  #selecting sex (breast or prostate cancer)
  targ.samples = sample.list[[sex]]
  
  #if only selecting discovery set 
  if (discovery.only == T) {
    targ.samples = targ.samples[targ.samples$data.partition %in% 'Disovery',]
  }
  
  #loading files
  files = list.files(pattern = 'AIX.*inserts.1000.all.300.q20.RDS',path=rawcountdir)
  
  #filtering only for samples of interest
  files = files[grep(paste(targ.samples$GRP_Id,collapse = '|'),files)]
  
  #reading files using mclapply
  setwd(rawcountdir)
  files.list = mclapply(files, function(x) {
    tmp = data.frame(readRDS(x)[,2],stringsAsFactors = F)
    colnames(tmp)[1] = gsub('.inserts.*','',x)
    return(tmp)
  },mc.cores = 10 )
  
  #creating dataframe to combine samples into a matrix
  windows =  readRDS(files[1])
  windows[,1] = gsub('01-','00-',windows[,1])
  return.df = data.frame(window = windows[,1])
  
  #converting list of counts to dataframe and formatting
  file.length = unlist(sapply(files.list,length))
  files.df = do.call('cbind',files.list[file.length > 0])
  files.df = cbind(return.df,files.df)
  remain.files = files[file.length == 0]
  colnames(files.df) = gsub('_Ct.*|_combined.*|_01_.*','',colnames(files.df))
  rownames(files.df) = files.df$window
  files.df = files.df[rowSums(files.df[,-1])>0, ]
  
  #saving combined raw count file 
  files.df=saveRDS(files.df, paste0('ohs.',sex,'all.rawcounts.RDS'))
  
  #normalizing raw counts using DESeQ2
  combined.filt = files.df[,-1]
  #filtering for rows with no signal
  combined.counts1 = combined.filt[rowSums(combined.filt) > 0,]
  
  #filtering out blacklist 300bp windows
  background.windows = rownames(combined.counts1)
  background.windows = background.windows[!background.windows %in% blacklist$window]
  
  dds <- DESeqDataSetFromMatrix(countData = combined.counts1[background.windows,combined.sample.info.filt$GRP_Id],
                                colData = combined.sample.info.filt,
                                design= ~  group ) #can add gender here but we're only looking at female samples currently
  
  dds <- estimateSizeFactors(dds) #estimating size factors for tmm normalization
  
  #computing normalized counts
  dds.matrix =counts(dds,normalize =T)
  
  #saving normalized counts + deseq object
  saveRDS(dds.matrix,paste0(savedir,'ohs.',sex,'all.samples.deseq.normcounts.RDS'))
  saveRDS(dds,paste0(savedir,'ohs.',sex,'all.samples.dds.RDS'))
  
  
  
  
}



####loading + normalizing counts for silencer/enhancer matrix####
marker.list = c('genhancer','silencer')

for (sex in names(sample.list))  {
  for (marker in marker.list) {
    for (inserts in c(1000)) {
      #selecting sex (breast or prostate cancer)
      targ.samples = sample.list[[sex]]
      
      #if only selecting discovery set 
      if (discovery.only == T) {
        targ.samples = targ.samples[targ.samples$data.partition %in% 'Disovery',]
      }
      
      #filtering count files for samples of interest
      files = list.files(pattern = paste0('AIX.*',marker,'.*'),path = rawcountdir)
      files = files[grepl(paste0('inserts.',inserts),files)]
      targ.samples = paste(combined.info$GRP_Id,collapse='|')
      targ.files = targ.files[grepl(targ.samples, targ.files)]
      
      #reading in individual raw counts into list
      setwd(rawcountdir)
      targ.list = lapply(targ.files, function(x){
        a = data.frame(readRDS(x)[,2])
        if (gsub('-','_',gsub('\\..*','',targ.files)) %in% sample.info$process.id) {
          colnames(a)[1] = gsub('_Ct.*|_combined.*|_01_.*','',x)
          return(a)
        }
        
      } )
      
      
      #merging list into a dataframe
      targ.df= do.call('cbind',targ.list)
      window= data.frame(window=readRDS(targ.files[1])[,1])
      combined.df = cbind(window,targ.df)
      rownames(combined.df) =combined.df$window
      
      #saving raw count matrix
      saveRDS(combined.df,paste0(savedir,'ohs.',inserts,'.',marker,'.raw.counts.RDS'))
      
      #normalizing raw counts
      
      #filtering for autosomes
      background = combined.df$window
      background = background[!grepl('chrX|chrY',background)]
      
      #adding window coordinates to rownames
      rownames(combined.df) = combined.df$window
      
      #reading raw counts into deseq
      dds <- DESeqDataSetFromMatrix(countData = combined.df[background,targ.samples$GRP_Id],
                                    colData = targ.samples,
                                    design= ~  group ) 

      #removing regions with no signal
      dds = dds[which(rowSums(combined.df[background,targ.samples$GRP_Id]) > 0),]
      dds <- estimateSizeFactors(dds) #estimating size factors for tmm normalization
      
      #computing normalized count into a dataframe
      dds.matrix =counts(dds,normalize =T)
      
      #saving fount files
      saveRDS(dds.matrix,paste0(savedir,'ohs.',inserts,'.',marker,'.norm.counts.RDS'))
      saveRDS(dds,paste0(savedir,'ohs.',inserts,'.',marker,'.dds.RDS'))
      
      
      
      
      
      
      
    }
    
    
  }
  
  
}

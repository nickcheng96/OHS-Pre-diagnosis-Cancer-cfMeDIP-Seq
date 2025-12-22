#loading libraries
library(parallel)
library(DESeq2)

#####read samples#####
targets = list('Breast' = sample.info.filt[sample.info.filt$Cancer %in% c('Breast','Control') & sample.info.filt$Sex == 'Female' ,],
               'Prostate'  = sample.info.filt[sample.info.filt$Cancer %in% c('Prostate','Control') & sample.info.filt$Sex == 'Male' ,])

wkdir='/User/workdirectory/'
figdir='/User/workdirectory/figures/qc/'
dir.create(figdir,recursive = T)
setwd(wkdir)

#####methylation insertts coverage####
cpg_count.all = readRDS('/.mounts/labs/awadallalab/private/ncheng/references/cpg_sites/cpg_site_positions/window_count/hg38_cpg_window_300_count.RDS') #number of cpg sites across 300bp regions
blacklist = readRDS('/.mounts/labs/awadallalab/private/ncheng/references/encode.blacklist/hg38-blacklist.v2_window300.RDS')

discovery.set = readRDS('/.mounts/labs/awadallalab/private/ncheng/manuscripts/early.cancer.risk/discovery.set.samples.RDS')
validation.set = readRDS('/.mounts/labs/awadallalab/private/ncheng/manuscripts/early.cancer.risk/validation.set.samples.RDS')
sample.info.filt=readRDS('/.mounts/labs/awadallalab/private/ncheng/manuscripts/early.cancer.risk/combined.set.samples.RDS')

rawcountdir='/path/to/raw/counts/'
setwd(rawcountdir)

#creating count matrix for 300bp windows
files = list.files(pattern = '.*inserts.1000.all.300.q20.RDS')
files.combined = files[grepl('combined',files)]
files.nontopup = files[!gsub('_Ct.*|_combined.*','',files) %in% gsub('_Ct.*|_combined.*','',files.combined) ]
files = c(files.combined,files.nontopup)
#file.types='inserts.400.all.300.q20.RDS'

sample.list = list(Male = discovery.set[discovery.set$Sex == 'Male' & discovery.set$Cancer %in% c('Control','Prostate'),],
                   Female = discovery.set[discovery.set$Sex == 'Female' & discovery.set$Cancer %in% c('Control','Breast'),],
                   All = discovery.set)

library(DESeq2)

#discovery only
file.types=c('inserts.1000.all.300.q20.RDS','inserts.400.all.300.q20.RDS')
for (sex in names(sample.list)[c(2,3,1)]) {
  for (f in file.types) {
    print(sex)
    print(f)
    raw=F
    targ.samples = sample.list[[sex]]
    if (raw == T){
      files = list.files(pattern = paste0('AIX.*',f))
      
      files = files[grep(paste(targ.samples$GRP_Id,collapse = '|'),files)]
      files.list = mclapply(files, function(x) {
        tmp = data.frame(readRDS(x)[,2],stringsAsFactors = F)
        colnames(tmp)[1] = gsub('.inserts.*','',x)
        return(tmp)
      },mc.cores = 10 )
      
      windows =  readRDS(files[1])
      windows[,1] = gsub('01-','00-',windows[,1])
      return.df = data.frame(window = windows[,1])
      
      file.rows = unlist(sapply(files.list,nrow))
      file.length = unlist(sapply(files.list,length))
      
      files.df = do.call('cbind',files.list[file.length > 0])
      files.df = cbind(return.df,files.df)
      remain.files = files[file.length == 0]
      colnames(files.df) = gsub('_Ct.*|_combined.*|_01_.*','',colnames(files.df))
      rownames(files.df) = files.df$window
      files.df = files.df[rowSums(files.df[,-1])>0, ]
      
      
      
    }
    
    files.df=readRDS(paste0('AIX13.discovery.',sex,'all.rawcounts.',f))
    
    
    #checking for mclapply that failed 
    savedir='/.mounts/labs/awadallalab/private/ncheng/cfmedip_data/cptp_samples/fragmentation/aix13.updated1/all.sample.fragmentation/'
    setwd(savedir)
    combined.filt = files.df[,-1]
    if (nrow(combined.filt) > 1){
      library(DESeq2)
      
      combined.counts1 = combined.filt[rowSums(combined.filt) > 0,]
      background.windows = rownames(combined.counts1)
      background.windows = background.windows[!background.windows %in% blacklist$window]
      #background.windows = background.windows[background.windows %in% cpg_count$window]
      combined.counts1 = combined.counts1[background.windows,]
      combined.sample.info.filt = targ.samples[targ.samples$GRP_Id %in% colnames(combined.counts1),]
      dds <- DESeqDataSetFromMatrix(countData = combined.counts1[,combined.sample.info.filt$GRP_Id],
                                    colData = combined.sample.info.filt,
                                    design= ~  group ) #can add gender here but we're only looking at female samples currently
      #setwd(wkdir)
      #dds = dds[which(rowMeans(combined.counts[,combined.sample.info.filt$GRP_Id]) != 1),]
      #nonzero.regions = rownames(dds)
      dds <- estimateSizeFactors(dds) #estimating size factors for tmm normalization
      #normalized counts
      dds.matrix =counts(dds,normalize =T)
      saveRDS(dds.matrix,paste0('','AIX13.allbg.discovery.',sex,'all.samples.deseq.normcounts.',f))
      saveRDS(dds,paste0('','AIX13.allbg.discovery.',sex,'all.samples.dds.',f))
      
    
      
    }
    
    
  }
  
  
}



#enhancer/silencers
marker.list = c('genhancer','silencer')#,'utr3','utr5') #,'ctcfbs','ucre')


#discoveyr set only
for (marker in marker.list) {
  for (inserts in c(1000)) {
    combined.info = discovery.set
    print(paste0('discovery.',marker))
    if (raw == T){
      files = list.files(pattern = paste0('AIX.*',marker,'.*'))
      targ.files = files[grepl(paste0('inserts.',inserts),files)]
      targ.samples = paste(combined.info$GRP_Id,collapse='|')
      targ.files = targ.files[grepl(targ.samples, targ.files)]
      
      targ.list = lapply(targ.files, function(x){
        #print(x)
        a = data.frame(readRDS(x)[,2])
        if (gsub('-','_',gsub('\\..*','',targ.files)) %in% sample.info$process.id) {
          colnames(a)[1] = gsub('_Ct.*|_combined.*|_01_.*','',x)
          return(a)
        }
        
      } )
      
      
      
      targ.df= do.call('cbind',targ.list)
      window= data.frame(window=readRDS(targ.files[1])[,1])
      combined.df = cbind(window,targ.df)
      rownames(combined.df) =combined.df$window
      saveRDS(combined.df,paste0(wkdir,'regulatory.counts/aix13.discovery.',inserts,'.',marker,'.raw.counts3.RDS'))
      
    }
    combined.df = readRDS(paste0(wkdir,'regulatory.counts/aix13.discovery.',inserts,'.',marker,'.raw.counts3.RDS'))
    print(dim(combined.df))
    #
    
    #adsf
    background = rownames(sample.matrix)
    background = background[!grepl('chrX|chrY',background)]
    sample.info = combined.info[combined.info$GRP_Id %in% colnames(combined.df),]
    medips.count_df.filt = combined.df[,sample.info$GRP_Id]
    rownames(medips.count_df.filt) = combined.df$window
    merged.df.filt = sample.info[sample.info$GRP_Id %in% colnames(medips.count_df.filt),]
    library(DESeq2)
    dds <- DESeqDataSetFromMatrix(countData = medips.count_df.filt[background,merged.df.filt$GRP_Id],
                                  colData = merged.df.filt,
                                  design= ~  seq.run ) #can add gender here but we're only looking at female samples currently
    #setwd(wkdir)
    dds = dds[which(rowSums(medips.count_df.filt[,merged.df.filt$GRP_Id]) > 0),]
    dds <- estimateSizeFactors(dds) #estimating size factors for tmm normalization
    dds.matrix =counts(dds,normalize =T)
    saveRDS(dds.matrix,paste0(wkdir,'regulatory.counts/aix13.discovery.',inserts,'.',marker,'.norm.counts.RDS'))
    saveRDS(dds,paste0(wkdir,'regulatory.counts/aix13.discovery.',inserts,'.',marker,'.dds.RDS'))
    
    
    

    
    
    
  }
  
  
}

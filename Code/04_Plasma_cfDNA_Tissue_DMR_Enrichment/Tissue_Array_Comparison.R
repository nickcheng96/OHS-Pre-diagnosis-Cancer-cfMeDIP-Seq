#loading libraries
library(ggplot2)
library(ggpubr)
library(scales)
library(GenomicRanges)
library(UpSetR)
library(Hmisc)
options(scipen=99990)
#####custom functions to be used #####
targ.filt = function(x) {
  a = x[x %in% cpg_count$window]
  a = a[!a %in% blacklist$window]
  a =a[!grepl('chrX|chrY',a)]
  #a = a[a %in% feature.filt[['brca']]]
  #a = a[!a %in% feature.filt[['blood']]]
  return(a)
}
illumina.450k.array.annotation= function(targ.windows,HM450.hg38.annotation) {
  #convert targeted windows to grange
  tmp = data.frame(chr = gsub(':.*','',targ.windows),
                   start = gsub('-.*','',gsub('.*:','',targ.windows)),
                   end = gsub('.*-','',targ.windows),
                   window = targ.windows)
  
  colnames(tmp)[c(1:3)] = c('chrom','chromStart','chromEnd')
  
  targ.grange <- makeGRangesFromDataFrame(tmp,
                                          keep.extra.columns=T)
  HM450.hg38.annotation.grange =  makeGRangesFromDataFrame(HM450.hg38.annotation,
                                                           keep.extra.columns=T)
  
  
  HM450.hg38.annotation.df =data.frame(HM450.hg38.annotation.grange[queryHits(findOverlaps(HM450.hg38.annotation.grange, targ.grange)),])   
  HM450.hg38.annotation.df$window = targ.grange[subjectHits(findOverlaps(HM450.hg38.annotation.grange, targ.grange)),]$window
  
  return(HM450.hg38.annotation.df)
}
z.score.permutation = function(cancer_blood_dmr_sig, predx.dmrs.sig.hyper,background, plot = F) {
  summary.df.all = NULL
  
  
  for (region in c('Islands','Shores','Shelves','Open Sea','All')) {
    if (region == 'All') {
      background.hm450 = background$window #[background$window %in% background$window]
      
    } else if (region == 'Islands') {
      background.hm450 = background[background$CGIposition %in% c('Island'),]$window
      
    } else if (region == 'Shores') {
      background.hm450 = background[background$CGIposition %in% c('N_Shore','S_Shore'),]$window
      
    } else if (region == 'Shelves') {
      background.hm450 = background[background$CGIposition %in% c('N_Shelf','S_Shelf'),]$window
      
    } else if (region == 'Open Sea') {
      background.hm450 = background[!background$CGIposition %in% c('Island','N_Shore','S_Shore','N_Shelf','S_Shelf'),]$window
      
    }
    
    #selecting tcga dmrs
    tcga.dmrs = background[background$probeID %in% cancer_blood_dmr_sig$probe,'window']
    #tcga.dmrs = unique(background.hm450[background.hm450 %in% tcga.dmrs])
    #selecting cfdna dmrs
    predx.hm450.dmrs = unique(background.hm450[background.hm450 %in% predx.dmrs.sig.hyper$window])
    #choosing cfdna dmrs overlapping tcga dmrs
    observed.overlap.windows = predx.hm450.dmrs[which(predx.hm450.dmrs %in% tcga.dmrs)]
    #computing size of overlap
    observed.overlap.count = length(observed.overlap.windows)
    random.selection.list = list()
    summary.df = NULL
    #length of cfdna dmrs
    n.size = length(predx.hm450.dmrs)
    #n.size=length(predx.dmrs.sig.hyper$window)
    print(n.size)
    background.windows = cpg.window.annotated[cpg.window.annotated$window %in% background.hm450 ,]
    cpg.window.annotated.sig.regions = cpg.window.annotated[cpg.window.annotated$window %in% predx.hm450.dmrs,]
    n.features = length(unique(cpg.window.annotated.sig.regions))
    search.space = unique(background.windows)
    cpg_count$window = as.character(cpg_count$window)
    search.space.window.distribution = cpg_count[cpg_count$window %in% cpg.window.annotated.sig.regions$window,]
    search.space.window.distribution = data.frame(table(search.space.window.distribution$count))
    
    cpg.count.split.background = lapply(cpg.count.split, function(x) x[x$window %in% background.windows$window,])
    
    #permutation testing with cpg distribution matching
    overall.perf.list = mclapply(1:5000, function(x) {
      set.seed(x)
      print(x)
      return.df.tmp.list = lapply(as.character(search.space.window.distribution$Var1), function(y) {
        search.space.cpg = cpg.count.split.background[[y]]
        n.features = search.space.window.distribution[search.space.window.distribution$Var1 == y,'Freq']
        return.df.tmp = search.space.cpg[sample(1:nrow(search.space.cpg),n.features),]
        
        return(as.character(return.df.tmp$window))
      } )
      return.df.random = unlist(return.df.tmp.list)
      
      
      return.df = background.windows[background.windows$window %in% return.df.random,]
      return.df$seed =x
      
      random.background = return.df$window
      tmp.overlap = unique(random.background[random.background %in% tcga.dmrs])
      tmp.df.return = data.frame(seed = x, group = 'permutation', region = region,n.overlap = length(tmp.overlap))
      
      return(tmp.df.return)
    },mc.cores=15 )
    summary.df = do.call('rbind',overall.perf.list)
    
    summary.df$zscore = (summary.df$n.overlap -mean(summary.df$n.overlap))/sd(summary.df$n.overlap)
    
    observed.df =  data.frame(seed = 0, group = 'observed', region = region,n.overlap = observed.overlap.count,
                              zscore= (observed.overlap.count -mean(summary.df$n.overlap))/sd(summary.df$n.overlap))
    
    summary.df = rbind(observed.df,summary.df)
    summary.df$sig.p = 2*pnorm(q=summary.df$zscore, lower.tail=FALSE)
    summary.df$group.size = n.size
    
    summary.df.all = rbind(summary.df.all,summary.df)
  }
  
  return(summary.df.all)
}
overlap.regions = function(cancer_blood_dmr_sig, predx.dmrs.sig.hyper,background, plot = F) {
  summary.df.all = NULL
  for (region in c('Islands','Shores','Shelves','Open Sea','All')) {
    if (region == 'All') {
      background.hm450 = background$window #[background$window %in% background$window]
      
    } else if (region == 'Islands') {
      background.hm450 = background[background$CGIposition %in% c('Island'),]$window
      
    } else if (region == 'Shores') {
      background.hm450 = background[background$CGIposition %in% c('N_Shore','S_Shore'),]$window
      
    } else if (region == 'Shelves') {
      background.hm450 = background[background$CGIposition %in% c('N_Shelf','S_Shelf'),]$window
      
    } else if (region == 'Open Sea') {
      background.hm450 = background[!background$CGIposition %in% c('Island','N_Shore','S_Shore','N_Shelf','S_Shelf'),]$window
      
    }
    
    #selecting tcga dmrs
    tcga.dmrs = background[background$probeID %in% cancer_blood_dmr_sig$probe,'window']
    #selecting cfdna dmrs
    predx.hm450.dmrs = unique(background.hm450[background.hm450 %in% predx.dmrs.sig.hyper$window])
    #choosing cfdna dmrs overlapping tcga dmrs
    observed.overlap.windows = predx.hm450.dmrs[which(predx.hm450.dmrs %in% tcga.dmrs)]
    #computing size of overlap
    observed.overlap.count = length(observed.overlap.windows)
    random.selection.list = list()
    summary.df = NULL
    if (length(observed.overlap.windows) >0) {
      summary.df = data.frame(windows = observed.overlap.windows, region = region )
      summary.df.all = rbind(summary.df.all,summary.df)
    }
    
  }
  if (plot == T) {
    overlap.regions = HM450.hg38.annotation[HM450.hg38.annotation$window %in% observed.overlap.windows,]
    
  }
  return(summary.df.all)
}
pancancer.permutation = function(predx.dmrs.sig,background,others = F){
  cancer.groups = c('prostate','breast','brain','skin','lung','colorectal','pancreatic','bladder','liver','oesophageal', 'uterus','thyroid','headneck')
  overall.results.pancanblood = NULL
  overall.results.pancannorm = NULL
  overall.results.normblood =NULL
  
  
  #all
  for (c in cancer.groups) {
    #cancer.dmr.dir='/.mounts/labs/awadallalab/private/ncheng/methylation_analysis/tissue_methylome_profiles/lm.dmp/'
    targ.cancer.dmr = readRDS(paste0(c,'.dmps.RDS')) #upload dmps
    d1=c('cfDNA Hyper, Cancer Hyper','cfDNA Hypo, Cancer Hypo')
    
    print(c)
    for (p in 0.05) {
      for (j in 0.10) {
        for (d in d1){
          if (d == 'cfDNA Hyper, Cancer Hyper') {
            cancer_blood_dmr_sig = targ.cancer.dmr[as.numeric(targ.cancer.dmr$padjust) < p & targ.cancer.dmr$diff> j & targ.cancer.dmr$comparison == 'cancer.blood',]
            cancer_normal_dmr_sig = targ.cancer.dmr[as.numeric(targ.cancer.dmr$padjust) < p & targ.cancer.dmr$diff> j & targ.cancer.dmr$comparison == 'cancer.norm',]
            normal_blood_dmr_sig = targ.cancer.dmr[as.numeric(targ.cancer.dmr$padjust) < p & targ.cancer.dmr$diff> j & targ.cancer.dmr$comparison == 'norm.blood',]
            
            predx.dmrs.sig.hyper = predx.dmrs.sig[predx.dmrs.sig$log2FoldChange > 0.25,]
            
          } else if (d == 'cfDNA Hyper, Cancer Hypo') {
            cancer_blood_dmr_sig = targ.cancer.dmr[as.numeric(targ.cancer.dmr$padjust) < p & targ.cancer.dmr$diff < -j & targ.cancer.dmr$comparison == 'cancer.blood',]
            cancer_normal_dmr_sig = targ.cancer.dmr[as.numeric(targ.cancer.dmr$padjust) < p & targ.cancer.dmr$diff < -j & targ.cancer.dmr$comparison == 'cancer.norm',]
            normal_blood_dmr_sig = targ.cancer.dmr[as.numeric(targ.cancer.dmr$padjust) < p & targ.cancer.dmr$diff <-j & targ.cancer.dmr$comparison == 'norm.blood',]
            
            predx.dmrs.sig.hyper = predx.dmrs.sig[predx.dmrs.sig$log2FoldChange > 0.25,]
            
          }  else if (d == 'cfDNA Hypo, Cancer Hypo') {
            cancer_blood_dmr_sig = targ.cancer.dmr[as.numeric(targ.cancer.dmr$padjust) < p & targ.cancer.dmr$diff < -j & targ.cancer.dmr$comparison == 'cancer.blood',]
            cancer_normal_dmr_sig = targ.cancer.dmr[as.numeric(targ.cancer.dmr$padjust) < p & targ.cancer.dmr$diff < -j & targ.cancer.dmr$comparison == 'cancer.norm',]
            normal_blood_dmr_sig = targ.cancer.dmr[as.numeric(targ.cancer.dmr$padjust) < p & targ.cancer.dmr$diff < -j & targ.cancer.dmr$comparison == 'norm.blood',]
            
            predx.dmrs.sig.hyper = predx.dmrs.sig[predx.dmrs.sig$log2FoldChange < -0.25,]
            
          }  else if (d == 'cfDNA Hypo, Cancer Hyper') {
            cancer_blood_dmr_sig = targ.cancer.dmr[as.numeric(targ.cancer.dmr$padjust) < p & targ.cancer.dmr$diff> j & targ.cancer.dmr$comparison == 'cancer.blood',]
            cancer_normal_dmr_sig = targ.cancer.dmr[as.numeric(targ.cancer.dmr$padjust) < p & targ.cancer.dmr$diff> j & targ.cancer.dmr$comparison == 'cancer.norm',]
            normal_blood_dmr_sig = targ.cancer.dmr[as.numeric(targ.cancer.dmr$padjust) < p & targ.cancer.dmr$diff> j & targ.cancer.dmr$comparison == 'norm.blood',]
            
            predx.dmrs.sig.hyper = predx.dmrs.sig[predx.dmrs.sig$log2FoldChange < -0.25,]
            
          } 
          predx.dmrs.sig.hyper = predx.dmrs.sig.hyper[order(predx.dmrs.sig.hyper$pvalue),]
          permutation.calc=T
          
          if (permutation.calc == T){
            results.df.tmp.cancerblood = z.score.permutation(cancer_blood_dmr_sig,
                                                             predx.dmrs.sig.hyper,
                                                             background) 
            results.df.tmp.cancerblood$pvalue = p
            results.df.tmp.cancerblood$abs.diff = j
            results.df.tmp.cancerblood$direction = gsub('BRCA','Tissue',d)
            results.df.tmp.cancerblood$tissue = c
            
            overall.results.pancanblood = rbind(overall.results.pancanblood,results.df.tmp.cancerblood)
            
            if (others == T){
              results.df.tmp.cancernorm = z.score.permutation(cancer_normal_dmr_sig, predx.dmrs.sig.hyper,background) 
              results.df.tmp.cancernorm$pvalue = p
              results.df.tmp.cancernorm$abs.diff = j
              results.df.tmp.cancernorm$direction = d
              results.df.tmp.cancernorm$tissue =c 
              
              overall.results.pancannorm = rbind(overall.results.pancannorm,results.df.tmp.cancernorm)
              
              results.df.tmp.normblood = z.score.permutation(normal_blood_dmr_sig, predx.dmrs.sig.hyper,background) 
              results.df.tmp.normblood$pvalue = p
              results.df.tmp.normblood$abs.diff = j
              results.df.tmp.normblood$direction = d
              results.df.tmp.normblood$tissue =c 
              
              overall.results.normblood = rbind(overall.results.normblood,results.df.tmp.normblood)
              
            }
            
            
          }
          
          
          
        }
        
      }
      
    }
    
  }
  if (others == T) {
    return(list('PanCan.Blood'=overall.results.pancanblood))#,
    
  } else {
    return(list('PanCan.Blood'=overall.results.pancanblood,
                'PanCan.Norm'=overall.results.pancannorm,
                'Norm.Blood'=overall.results.normblood))
    
  }
  
}

permutation.function = function(cancer_blood_dmr,
                                cancer_normal_dmr,
                                normal_blood_dmr ,
                                predx.dmrs.sig,
                                effect,
                                p.threshold,
                                direct,
                                HM450.hg38.annotation) { 
  overall.results.cancerpbl = NULL
  overall.results.cancernorm = NULL
  overall.results.normblood= NULL
  for (p in p.threshold) {
    for (j in effect) {
      for (d in direct){
        if (d == 'cfDNA Hyper, Cancer Hyper') {
          cancer_blood_dmr_sig = cancer_blood_dmr[as.numeric(cancer_blood_dmr$padjust) < p & cancer_blood_dmr$diff> j,]
          cancer_normal_dmr_sig = cancer_normal_dmr[as.numeric(cancer_normal_dmr$padjust) < p & cancer_normal_dmr$diff> j,]
          normal_blood_dmr_sig= normal_blood_dmr[as.numeric(normal_blood_dmr$padjust) <p & normal_blood_dmr$diff > j,]
          predx.dmrs.sig.hyper = predx.dmrs.sig[predx.dmrs.sig$log2FoldChange > 0.25,]
          
        } else if (d == 'cfDNA Hyper, Cancer Hypo') {
          cancer_blood_dmr_sig = cancer_blood_dmr[as.numeric(cancer_blood_dmr$padjust) < p & cancer_blood_dmr$diff< -j,]
          cancer_normal_dmr_sig = cancer_normal_dmr[as.numeric(cancer_normal_dmr$padjust) < p & cancer_normal_dmr$diff< -j,]
          normal_blood_dmr_sig= normal_blood_dmr[as.numeric(normal_blood_dmr$padjust) <p & normal_blood_dmr$diff < -j,]
          
          predx.dmrs.sig.hyper = predx.dmrs.sig[predx.dmrs.sig$log2FoldChange > 0.25,]
          
        }  else if (d == 'cfDNA Hypo, Cancer Hypo') {
          cancer_blood_dmr_sig = cancer_blood_dmr[as.numeric(cancer_blood_dmr$padjust) < p & cancer_blood_dmr$diff< -j,]
          cancer_normal_dmr_sig = cancer_normal_dmr[as.numeric(cancer_normal_dmr$padjust) < p & cancer_normal_dmr$diff< -j,]
          normal_blood_dmr_sig= normal_blood_dmr[as.numeric(normal_blood_dmr$padjust) <p & normal_blood_dmr$diff < -j,]
          
          predx.dmrs.sig.hyper = predx.dmrs.sig[predx.dmrs.sig$log2FoldChange < -0.25,]
          
        }  else if (d == 'cfDNA Hypo, Cancer Hyper') {
          cancer_blood_dmr_sig = cancer_blood_dmr[as.numeric(cancer_blood_dmr$padjust) < p & cancer_blood_dmr$diff> j,]
          cancer_normal_dmr_sig = cancer_normal_dmr[as.numeric(cancer_normal_dmr$padjust) < p & cancer_normal_dmr$diff> j,]
          normal_blood_dmr_sig= normal_blood_dmr[as.numeric(normal_blood_dmr$padjust) <p & normal_blood_dmr$diff > j,]
          
          predx.dmrs.sig.hyper = predx.dmrs.sig[predx.dmrs.sig$log2FoldChange < -0.25,]
          
        } 
        predx.dmrs.sig.hyper = predx.dmrs.sig.hyper[order(predx.dmrs.sig.hyper$pvalue),]
        
        results.df.tmp.cancerblood = z.score.permutation(cancer_blood_dmr_sig, 
                                                         predx.dmrs.sig.hyper,
                                                         background = HM450.hg38.annotation) 
        results.df.tmp.cancerblood.features = overlap.regions(cancer_blood_dmr_sig, 
                                                              predx.dmrs.sig.hyper,
                                                              background = HM450.hg38.annotation) 
        results.df.tmp.cancerblood$pvalue = p
        results.df.tmp.cancerblood$abs.diff = j
        results.df.tmp.cancerblood$direction = d
        results.df.tmp.cancerblood$comparison = 'Cancer.Blood'
        overall.results.cancerpbl = rbind(overall.results.cancerpbl,results.df.tmp.cancerblood)
        
        results.df.tmp.cancernorm = z.score.permutation(cancer_normal_dmr_sig, 
                                                        predx.dmrs.sig.hyper,
                                                        background = HM450.hg38.annotation) 
        results.df.tmp.cancernorm$pvalue = p
        results.df.tmp.cancernorm$abs.diff = j
        results.df.tmp.cancernorm$direction = d
        results.df.tmp.cancernorm$comparison = 'Cancer.Normal'
        overall.results.cancernorm = rbind(overall.results.cancernorm,results.df.tmp.cancernorm)
        
        
        results.df.tmp.normblood = z.score.permutation(normal_blood_dmr_sig, 
                                                       predx.dmrs.sig.hyper,
                                                       background = HM450.hg38.annotation) 
        results.df.tmp.normblood$pvalue = p
        results.df.tmp.normblood$abs.diff = j
        results.df.tmp.normblood$direction = d
        results.df.tmp.normblood$comparison = 'Normal.Blood'
        
        overall.results.normblood = rbind(overall.results.normblood,results.df.tmp.normblood)
        
        
        
      }
      
    }
    
  }
  return(list(Cancer.Blood=overall.results.cancerpbl, Cancer.Normal=overall.results.cancernorm,Normal.Blood=overall.results.normblood))
}

feature.overlaps = function(cancer_blood_dmr,
                            cancer_normal_dmr,
                            normal_blood_dmr ,
                            
                            predx.dmrs.sig,
                            effect,
                            p.threshold,
                            direct,
                            background) 
{ 
  overall.results.cancerpbl = NULL
  overall.results.cancernorm = NULL
  overall.results.normblood= NULL
  
  for (p in p.threshold) {
    for (j in effect) {
      for (d in direct){
        if (d == 'cfDNA Hyper, Cancer Hyper') {
          cancer_blood_dmr_sig = cancer_blood_dmr[as.numeric(cancer_blood_dmr$padjust) < p & cancer_blood_dmr$diff> j,]
          cancer_normal_dmr_sig = cancer_normal_dmr[as.numeric(cancer_normal_dmr$padjust) < p & cancer_normal_dmr$diff> j,]
          normal_blood_dmr_sig= normal_blood_dmr[as.numeric(normal_blood_dmr$padjust) <p & normal_blood_dmr$diff > j,]
          predx.dmrs.sig.hyper = predx.dmrs.sig[predx.dmrs.sig$log2FoldChange > 0.25,]
          
        } else if (d == 'cfDNA Hyper, Cancer Hypo') {
          cancer_blood_dmr_sig = cancer_blood_dmr[as.numeric(cancer_blood_dmr$padjust) < p & cancer_blood_dmr$diff< -j,]
          cancer_normal_dmr_sig = cancer_normal_dmr[as.numeric(cancer_normal_dmr$padjust) < p & cancer_normal_dmr$diff< -j,]
          normal_blood_dmr_sig= normal_blood_dmr[as.numeric(normal_blood_dmr$padjust) <p & normal_blood_dmr$diff < -j,]
          
          predx.dmrs.sig.hyper = predx.dmrs.sig[predx.dmrs.sig$log2FoldChange > 0.25,]
          
        }  else if (d == 'cfDNA Hypo, Cancer Hypo') {
          cancer_blood_dmr_sig = cancer_blood_dmr[as.numeric(cancer_blood_dmr$padjust) < p & cancer_blood_dmr$diff< -j,]
          cancer_normal_dmr_sig = cancer_normal_dmr[as.numeric(cancer_normal_dmr$padjust) < p & cancer_normal_dmr$diff< -j,]
          normal_blood_dmr_sig= normal_blood_dmr[as.numeric(normal_blood_dmr$padjust) <p & normal_blood_dmr$diff < -j,]
          
          predx.dmrs.sig.hyper = predx.dmrs.sig[predx.dmrs.sig$log2FoldChange < -0.25,]
          
        }  else if (d == 'cfDNA Hypo, Cancer Hyper') {
          cancer_blood_dmr_sig = cancer_blood_dmr[as.numeric(cancer_blood_dmr$padjust) < p & cancer_blood_dmr$diff> j,]
          cancer_normal_dmr_sig = cancer_normal_dmr[as.numeric(cancer_normal_dmr$padjust) < p & cancer_normal_dmr$diff> j,]
          normal_blood_dmr_sig= normal_blood_dmr[as.numeric(normal_blood_dmr$padjust) <p & normal_blood_dmr$diff > j,]
          
          predx.dmrs.sig.hyper = predx.dmrs.sig[predx.dmrs.sig$log2FoldChange < -0.25,]
          
        } 
        predx.dmrs.sig.hyper = predx.dmrs.sig.hyper[order(predx.dmrs.sig.hyper$pvalue),]
        
        
        
        results.df.tmp.cancerblood = overlap.regions(cancer_blood_dmr_sig, 
                                                     predx.dmrs.sig.hyper,
                                                     background = background) 
        results.df.tmp.cancerblood$pvalue = p
        results.df.tmp.cancerblood$abs.diff = j
        results.df.tmp.cancerblood$direction = d
        results.df.tmp.cancerblood$comparison = 'Cancer.Blood'
        overall.results.cancerpbl = rbind(overall.results.cancerpbl,results.df.tmp.cancerblood)
        
        results.df.tmp.cancernorm = overlap.regions(cancer_normal_dmr_sig, 
                                                    predx.dmrs.sig.hyper,
                                                    background) 
        results.df.tmp.cancernorm$pvalue = p
        results.df.tmp.cancernorm$abs.diff = j
        results.df.tmp.cancernorm$direction = d
        results.df.tmp.cancernorm$comparison = 'Cancer.Normal'
        overall.results.cancernorm = rbind(overall.results.cancernorm,results.df.tmp.cancernorm)
        
        
        results.df.tmp.normblood = overlap.regions(normal_blood_dmr_sig, 
                                                   predx.dmrs.sig.hyper,
                                                   background) 
        results.df.tmp.normblood$pvalue = p
        results.df.tmp.normblood$abs.diff = j
        results.df.tmp.normblood$direction = d
        results.df.tmp.normblood$comparison = 'Normal.Blood'
        
        overall.results.normblood = rbind(overall.results.normblood,results.df.tmp.normblood)
        
        
        
      }
      
    }
    
  }
  return(list(Cancer.Blood=overall.results.cancerpbl, Cancer.Normal=overall.results.cancernorm,Normal.Blood=overall.results.normblood))
  
}


permutation.plot.manuscript = function(p.threshold=0.05, effect=0.1,name,enhancer.permutation,sex) {
  p=p.threshold  
  j=effect
  dmr.overlap.list = split(enhancer.permutation,enhancer.permutation$marker)
  for (dmr.targ in names(dmr.overlap.list)){
    for (d in c('cfDNA Hyper, Cancer Hyper','cfDNA Hypo, Cancer Hypo')) {
      d.name = gsub(' ','',gsub(',','',d))
      targ.df = dmr.overlap.list[[dmr.targ]]
      targ.results = targ.df[targ.df$direction == d,]
      targ.results$region =factor(targ.results$region, levels = c('All','Islands','Shores','Shelves','Open Sea'))
      targ.results$Sig = ifelse(targ.results$sig.p < 0.05, 'Sig','Non-Sig')
      regioncol = c('Islands' ='#29335C','Shores' = '#F3A712','Shelves' = '#DB2B39','Open Sea' = '#A89B83','All' = 'black')
      plot1 = ggplot(targ.results[targ.results$group == 'permutation' ,], aes(x = region,y=zscore,col= region)) + #& targ.results$region != 'All'
        geom_boxplot(outlier.shape = NA, aes(x = region,y=zscore,color = region)) + 
        scale_color_manual(values = regioncol) +
        
        geom_point(data = targ.results[targ.results$group == 'observed' ,],mapping = aes(x = region, y = zscore,fill = Sig),shape = 23, col = 'black')+
        scale_fill_manual(values = c('Non-Sig'='grey','Border.sig' = '#2A9D8F',Sig='red')) +
        
        scale_shape_manual(values=c(23))+
        theme_bw()+
        theme(text = element_text(size=15),
              axis.ticks.y = element_blank(),
              axis.ticks.x = element_blank(),
              legend.position = 'none',
              legend.title = element_blank(),
        ) + 
        scale_y_continuous( labels = label_number(accuracy = 0.1)) +
        xlab('Region') + ylab('Overlapping DMRs\n(z-score normalized)') 
      
      plot2 = ggplot(targ.results[targ.results$group == 'observed' ,], aes(x = region,y=n.overlap,fill = region)) + #&targ.results$region != 'All'
        geom_bar(stat= 'identity')+
        scale_fill_manual(values = regioncol) +
        theme_bw()+
        theme(text = element_text(size=15),
              axis.ticks.y = element_blank(),
              axis.ticks.x = element_blank(),
              axis.title.x = element_blank(),
              axis.text.x = element_blank(),
              legend.position = 'none',
              legend.title = element_blank()) +
        xlab('Region') + ylab('Overlapping DMRs\n(Count)') 
      
      figure <- ggarrange(plot2, plot1,
                          labels = c(""),
                          ncol = 1, nrow = 2,
                          heights = c(0.3, 0.5),
                          widths = c(1.5,0.7),
                          align = 'v')
      
      png(paste0(name,'.permutation.overlap.',dmr.targ,'.',d.name,'.all','.png'), height= 1500, width= 850,res=250)
      print(figure)
      dev.off()
      
      
    }
    targ.df = dmr.overlap.list[[dmr.targ]]
    targ.results = targ.df
    
    targ.results$region =factor(targ.results$region, levels = c('All','Islands','Shores','Shelves','Open Sea'))
    targ.results$Sig = ifelse(targ.results$sig.p < 0.05, 'Sig','Non-Sig')
    
    regioncol = c('Islands' ='#29335C','Shores' = '#F3A712','Shelves' = '#DB2B39','Open Sea' = '#A89B83','All' = 'black')
    
    for (region in  c('All','Islands','Shores','Shelves','Open Sea')) {
      targ.results1 = targ.results[targ.results$region == region,]
      targ.results1$direction = gsub(', ','\n',targ.results1$direction)
      
      if (region == 'All') {
        plot1 = ggplot(targ.results1[targ.results1$group == 'permutation',], aes(x = direction,y=zscore)) +
          geom_boxplot(outlier.shape = NA, aes(x = direction,y=zscore)) + 
          scale_color_manual(values = regioncol) +
          
          geom_point(data = targ.results1[targ.results1$group == 'observed',],mapping = aes(x = direction, y = zscore,fill = Sig),shape = 23,col = 'black')+
          scale_fill_manual(values = c('Non-Sig'='grey','Border.sig' = '#2A9D8F',Sig='red')) +
          scale_shape_manual(values=c(23))+
          theme_bw()+
          theme(text = element_text(size=14),
                axis.ticks.y = element_blank(),
                axis.ticks.x = element_blank(),
                legend.position = 'none',
                legend.title = element_blank()) + 
          xlab('Region') +ylab('Overlapping DMRs\n(z-score normalized)') 
        
        
        plot2 = ggplot(targ.results[targ.results$group == 'observed' & targ.results$region != 'All' ,], aes(x = direction,y=n.overlap, fill = region)) + 
          geom_bar(stat= 'identity')+
          scale_fill_manual(values = regioncol) +
          theme_bw()+
          theme(text = element_text(size=14),
                axis.ticks.y = element_blank(),
                axis.ticks.x = element_blank(),
                axis.title.x = element_blank(),
                axis.text.x = element_blank(),
                legend.position = 'none',
                legend.title = element_blank()) +
          xlab('Region') +  ylab('Overlapping DMRs\n(Count)') 
        
      } else{
        plot1 = ggplot(targ.results1[targ.results1$group == 'permutation',], aes(x = direction,y=zscore)) +
          geom_boxplot(outlier.shape = NA, aes(x = direction,y=zscore)) + 
          scale_color_manual(values = regioncol) +
          
          geom_point(data = targ.results1[targ.results1$group == 'observed' ,],mapping = aes(x = direction, y = zscore,fill = Sig),shape = 23,col = 'black')+
          scale_fill_manual(values = c('Non-Sig'='grey','Border.sig' = '#2A9D8F',Sig='red')) +
          scale_shape_manual(values=c(23))+
          theme_bw()+
          theme(text = element_text(size=14),
                axis.ticks.y = element_blank(),
                axis.ticks.x = element_blank(),
                legend.position = 'none',
                legend.title = element_blank()) + 
          xlab('Region') +ylab('Overlapping DMRs\n(z-score normalized)') 
        
        
        plot2 = ggplot(targ.results1[targ.results1$group == 'observed' ,], aes(x = direction,y=n.overlap, fill = region)) + 
          geom_bar(stat= 'identity')+
          scale_fill_manual(values = regioncol) +
          theme_bw()+
          theme(text = element_text(size=14),
                axis.ticks.y = element_blank(),
                axis.ticks.x = element_blank(),
                axis.title.x = element_blank(),
                axis.text.x = element_blank(),
                legend.position = 'none',
                legend.title = element_blank()) +
          xlab('Region') + ylab('Overlapping DMRs\n(Count)') 
        
      }
      figure <- ggarrange(plot2, plot1,
                          labels = c(""),
                          ncol = 1, nrow = 2,
                          heights = c(0.3, 0.5),
                          widths = c(0.1,0.7),align = 'v')
      file = gsub(' ','.',paste0('permutation.overlap.dir.all',region,'.png'))
      
      
      png(paste0(name,'.permutation.overlap.',dmr.targ,'.',region,'.png'), height= 1500, width= 1400,res=250)
      print(figure)
      dev.off()
      
    }
    
    
  }
  
  
  
}

####running 450k array + cfDNA MeDIP DMR enrichment/permutation analysis####
#setting background regions for permutation analysis
cpg.window.annotated = readRDS('/.mounts/labs/awadallalab/private/ncheng/references/gene.annotations/genomic.annotations.window.300.RDS') #upload remove sex chr and blacklist regions
cpg_count = readRDS('/.mounts/labs/awadallalab/private/ncheng/references/cpg_sites/cpg_site_positions/window_count/hg38_cpg_window_300_count.RDS') #upload #number of cpg sites across 300bp regions
cpg_count$window = as.character(cpg_count$window)
cpg_count = cpg_count[cpg_count$count >= 5,] #upload version with 5+ CpGs
background=cpg_count$window
cpg.count.split = split(cpg_count,cpg_count$count)

#setting directories
savedir='/.mounts/labs/awadallalab/private/ncheng/cfmedip_data/cptp_samples/fragmentation/aix13.updated1/all.sample.fragmentation/prostate.cancer.cv/regulatory.regions.v5.ageadjusted/1/genhancer.1000/'
dmrdir='/.mounts/labs/awadallalab/private/ncheng/cfmedip_data/cptp_samples/fragmentation/aix13.updated1/all.sample.fragmentation/prostate.cancer.cv/regulatory.regions.v5.ageadjusted/1/'
figdir='/path/to/figures/'
#loading files
#450k array annotation #upload
HM450.hg38.annotation = readRDS('/.mounts/labs/awadallalab/private/ncheng/annotation_files/meth_array/hg38/HM450.hg38.annotation.window300.RDS')
HM450.hg38.annotation.base=HM450.hg38.annotation

#tissue methylation array dmrs

#loading prostate cancer tissue methylation array differentially methylated positions previously computed using a linear model
cancer.dmr.dir='/.mounts/labs/awadallalab/private/ncheng/methylation_analysis/tissue_methylome_profiles/lm.dmp/'
tcga.dmr.list = list()
for (cancer in c('Prostate','Breast')) {
  targ.cancer.dmr = readRDS(paste0(cancer.dmr.dir,cancer,'.dmps.RDS')) #upload add capitalization to cancers
  cancer_blood_dmr_sig = targ.cancer.dmr[as.numeric(targ.cancer.dmr$padjust) < 0.05 & targ.cancer.dmr$comparison == 'cancer.blood',]
  cancer_normal_dmr_sig = targ.cancer.dmr[as.numeric(targ.cancer.dmr$padjust) < 0.05 & targ.cancer.dmr$comparison == 'cancer.norm',]
  normal_blood_dmr_sig = targ.cancer.dmr[as.numeric(targ.cancer.dmr$padjust) < 0.05 &  targ.cancer.dmr$comparison == 'norm.blood',]
  
  tcga.dmr.list[[cancer]] = list('Cancer.Blood' = cancer_blood_dmr_sig,
                                       'Cancer.Normal' = cancer_normal_dmr_sig,
                                       'Normal.Blood' = normal_blood_dmr_sig)

}

#loading previously computed cfDNA MedIP-Seq DMRs 
cfdna.dmr.list=list()
for (sex in c('Male','Female')) {
  res.df = paste0(readRDS('discovery.',sex,'.dmrs'))
  res.df$window=rownames(res.df)
  res.df = res.df[!grepl('chrX|chrY',res.df$window),]
  res.df.hyper.sig.raw = res.df[res.df$log2FoldChange > 0.25 & res.df$pvalue < 0.05,]
  res.df.hypo.sig.raw = res.df[res.df$log2FoldChange < -0.25 & res.df$pvalue < 0.05,]
  res.df.hyper.sig.raw = res.df.hyper.sig.raw[order(res.df.hyper.sig.raw$pvalue),][1:2000,]
  res.df.hypo.sig.raw = res.df.hypo.sig.raw[order(res.df.hypo.sig.raw$pvalue),][1:2000,]
  ifelse(sex == 'Male','Prostate','Breast')
  cfdna.dmr.list[[cancer]] = list()
  cfdna.dmr.list[[cancer]][['Hypermethylated']] = res.df.hyper.sig.raw
  cfdna.dmr.list[[cancer]][['Hypomethylated']] = res.df.hypo.sig.raw
  
  
}

#selecting MeDIP-seq hypermethylated and hypomethylated regions that overlaps with at least one 450k methylation array CpG site
overlap.list = list()
for (cancer in c('Breast','Prostate')) {
  hyper.tmp = cfdna.dmr.list[[cancer]][['Hypermethylated']]
  hypo.tmp = cfdna.dmr.list[[cancer]][['Hypomethylated']]
  
  hyper.overlap =HM450.hg38.annotation[HM450.hg38.annotation$window %in% hyper.tmp$window,]
  hypo.overlap =HM450.hg38.annotation[HM450.hg38.annotation$window %in% hypo.tmp$window,]
  
  overlap.list[[paste0(cancer,'.Hypermethylated')]] = hyper.overlap
  overlap.list[[paste0(cancer,'.Hypomethylated')]] = hypo.overlap
}


#performing permutation tests to assess enrichment of cancer tissue methylation array DMPs and cfDNA MeDIP-Seq DMRs
p.threshold = 0.05
effect = c(0.1)
direct = c('cfDNA Hyper, Cancer Hyper','cfDNA Hypo, Cancer Hypo') #

overall.results.cancerpbl = NULL


#peforming permutation tests
combined.sample.permutation=NULL
for (marker in names(cfdna.dmr.list)) {
  
  for (r in names(cfdna.dmr.list[[marker]])) {
    sig.regions= cfdna.dmr.list[[marker]][[r]]
    sig.regions$window =rownames(sig.regions)

    #array probes overlapping with background cfDNA (for selecting background reginos later to compute null distribution for permutation test)
    background.regions.array.450k = illumina.450k.array.annotation(cpg_count$window,HM450.hg38.annotation.base)
    
    
    cancer_blood_dmr= tcga.dmr.list[[marker]][['Cancer.Blood']]
    cancer_normal_dmr= tcga.dmr.list[[marker]][['Cancer.Normal']]
    normal_blood_dmr= tcga.dmr.list[[marker]][['Normal.Blood']]
    
    direct = c('cfDNA Hyper, Cancer Hyper','cfDNA Hyper, Cancer Hypo','cfDNA Hypo, Cancer Hyper','cfDNA Hypo, Cancer Hypo')
    
    if(r == 'Hypermethylated') {
      direct =c('cfDNA Hyper, Cancer Hyper','cfDNA Hyper, Cancer Hypo')
      
    } else {
      direct= c('cfDNA Hypo, Cancer Hyper','cfDNA Hypo, Cancer Hypo')
      
    }
    
    permutation.results = permutation.function(cancer_blood_dmr,
                                                cancer_normal_dmr,
                                                normal_blood_dmr,
                                                predx.dmrs.sig=sig.regions,
                                                effect=c(0.1),
                                                p.threshold=c(0.05),
                                                direct=direct,
                                                HM450.hg38.annotation=background.regions.array.450k)
    
    
    permutation.df= do.call('rbind',permutation.results)
    permutation.df$marker = marker 
    combined.sample.permutation=rbind(combined.sample.permutation,permutation.df)
    
    
    
    
  }
  
}


#isolating specific markers overlapping between cfDNA and tissue array methylation 
combined.features.all = NULL
for (marker in names(cfdna.dmr.list)) {
  
  for (r in names(cfdna.dmr.list[[marker]])) {
    
    sig.regions= cfdna.dmr.list[[marker]][[r]]
    sig.regions$window =rownames(sig.regions)

    background.regions.array.450k = illumina.450k.array.annotation(cpg_count$window,HM450.hg38.annotation)
    
    
    cancer_blood_dmr= tcga.dmr.list[[marker]][['Cancer.Blood']]
    cancer_normal_dmr= tcga.dmr.list[[marker]][['Cancer.Normal']]
    normal_blood_dmr= tcga.dmr.list[[marker]][['Normal.Blood']]
    
    direct = c('cfDNA Hyper, Cancer Hyper','cfDNA Hyper, Cancer Hypo','cfDNA Hypo, Cancer Hyper','cfDNA Hypo, Cancer Hypo')
    
    if(r == 'Hypermethylated') {
      direct =c('cfDNA Hyper, Cancer Hyper','cfDNA Hyper, Cancer Hypo')
      
    } else {
      direct= c('cfDNA Hypo, Cancer Hyper','cfDNA Hypo, Cancer Hypo')
      
    }
    
    overlapping.markers = feature.overlaps(cancer_blood_dmr,
                                           cancer_normal_dmr,
                                           normal_blood_dmr,
                                           predx.dmrs.sig=sig.regions,
                                           effect=c(0.1),
                                           p.threshold=c(0.05),
                                           direct=direct,
                                           background=background.regions.array.450k)
    
    
    overlapping.featurse.all= do.call('rbind',overlapping.markers)
    overlapping.featurse.all$marker = marker 
    combined.features.all=rbind(combined.features.all,overlapping.featurse.all)
    
    
    
    
  }
  
}


#plotting permutation results
comparisons = c('Cancer.Blood','Cancer.Normal','Normal.Blood')
cancer = c('Prostate','Breast')
pvalue = c(0.05)
effect = c(0.10)

for (p in pvalue ) {
  for (e in effect) {
    for (comp in comparisons) {
      for (c in cancer) {
        permutation.results = combined.sample.permutation[combined.sample.permutation$marker == c & 
                                                             combined.sample.permutation$comparison == comp &
                                                             combined.sample.permutation$pvalue == p & 
                                                             combined.sample.permutation$abs.diff == e,]
        
        if (nrow(permutation.results) > 0) {

          permutation.plot.manuscript(p.threshold = 0.05,
                                      effect = effect, 
                                      name = paste0(figdir,comp,'.',c,'.',p,'.',e), 
                                      permutation.results,
                                      sex = ifelse(c == 'Prostate','Male','Female'))
          
          
        }
        
        
      }
    }
    
    
  }
}

#plotting overlapping featuers
directions = c('cfDNA Hyper, Cancer Hyper','cfDNA Hypo, Cancer Hypo')
cancers = c('Breast','Prostate')

for (m in cancers) {
  for (d in directions){
    overlapping.target.markers  =combined.features.all[combined.features.all$direction == d & 
                                                                combined.features.all$marker == m,]
    tmp1 = split(overlapping.target.markers,overlapping.target.markers$comparison)

    overlap.hyper = HM450.hg38.annotation[,c('window','CGIposition','gene')]
    overlap.hyper$"Cancer vs Normal" = ifelse(overlap.hyper$window %in% tmp1[["Cancer.Normal"]]$window,1,0)
    overlap.hyper$"Cancer vs PBL"= ifelse(overlap.hyper$window %in% tmp1[["Cancer.Blood"]]$window,1,0)
    overlap.hyper$"Normal vs PBL"= ifelse(overlap.hyper$window %in% tmp1[["Normal.Blood"]]$window,1,0)
    
    overlap.hyper$cpg_region = gsub('N_','',gsub('S_','',overlap.hyper$CGIposition))
    overlap.hyper[is.na(overlap.hyper$cpg_region),'cpg_region'] = 'null'
    overlap.hyper = unique(overlap.hyper[,c('window','cpg_region','Normal vs PBL','Cancer vs PBL',"Cancer vs Normal")])
    
    sets =c("Normal vs PBL" ,"Cancer vs PBL","Cancer vs Normal")
    sets = gsub('cancers',m,sets)
    png(paste0(figdir,m,'.tissue.blood.',d,'.overlap.png'),height = 800, width =1500,res=300)
    print(upset(overlap.hyper,nsets = 9,keep.order = TRUE, sets = sets,
                queries = list(
                  list(query = elements, 
                       params = list("cpg_region", c("Island","Shore","Shelf","null")), color = "#A89B83", active = T),
                  list(query = elements, 
                       params = list("cpg_region", c("Island","Shore","Shelf")), color = "#DB2B39", active = T),
                  list(query = elements, 
                       params = list("cpg_region", c("Island","Shore")), color = '#F3A712', active = T),
                  list(query = elements, 
                       params = list("cpg_region", c("Island")), color = '#29335C', active = T))))
    
    
    dev.off()
  }
  
  
}


#computing permutation analysis with pan-cancer tissues DMPs and cfDNA DMRs
combined.pancan.permutation = list()
for (marker in names(cfdna.dmr.list)) {
  predx.dmrs.sig = do.call('rbind',cfdna.dmr.list[[marker]])
  tmp = pancancer.permutation(predx.dmrs.sig, background=HM450.hg38.annotation)
  tmp$marker = marker
  combined.pancan.permutation[[marker]] = tmp
  
}


#plotting pancancer permutation test
for (c in c('Breast','Prostate')) {
  targ.permutation.list = combined.pancan.permutation[[c]]
  comparison.list = list('PanCan.Blood' = targ.permutation.list[['PanCan.Blood']])
  



  p=0.05     
  j=0.1
  dmr.group=c('cfDNA Hyper, Tissue Hyper','cfDNA Hypo, Tissue Hypo')
  
  dmrs = list(pancanblood = targ.permutation.list[[1]])
  
  for (dmr in names(dmrs)) {
    targ.overlap = dmrs[[dmr]]
    if (c == 'Prostate') {
      targ.results = targ.overlap[!targ.overlap$tissue %in% c('breast'),]
      
      
    } else {
      targ.results = targ.overlap
      
      
    }
    targ.results$region =factor(targ.results$region, levels = c('All','Islands','Shores','Shelves','Open Sea'))
    targ.results$Sig = ifelse(targ.results$sig.p < 0.05, 'Sig','Non-Sig')
    
    
    targ.results$direction = gsub('PRAD','Tissue',targ.results$direction)
    regioncol = c('Islands' ='#29335C','Shores' = '#F3A712','Shelves' = '#DB2B39','Open Sea' = '#A89B83','All' = 'black')
    for (d in dmr.group) {
      for(region in  c('All','Islands','Shores','Shelves','Open Sea')) { #
        if (region =='All') {
          targ.results1 = targ.results[targ.results$region != region,]
          targ.results2 = targ.results[targ.results$region == region,]
          
        } else {
          targ.results1 = targ.results[targ.results$region == region,]
          targ.results2 = targ.results[targ.results$region == region,]
          
          
        }
        
        
        targ.results1 = targ.results1[targ.results1$direction == d,]
        targ.results2 = targ.results2[targ.results2$direction == d,]
        
        
        targ.results1$direction = gsub(', ','\n',targ.results1$direction)
        targ.results1$tissue  = capitalize(targ.results1$tissue)
        targ.results2$tissue = capitalize(targ.results2$tissue)
        targ.results2[is.na(targ.results2$zscore),'zscore'] = 0
        tmp = targ.results2[targ.results2$seed ==0,]
        tmp= unique(tmp[order(tmp$zscore),])
        targ.results1$tissue  =factor( targ.results1$tissue,levels = unique(tmp$tissue))
        targ.results2$tissue = factor( targ.results2$tissue,levels = unique(tmp$tissue))
        targ.results1[is.na(targ.results1$zscore),'zscore'] =0
        
        
        plot1 = ggplot(targ.results2[targ.results2$group == 'permutation',], aes(x = tissue,y=zscore)) +
          geom_boxplot(outlier.shape = NA, aes(x = tissue,y=zscore),col = 'black') + 
          scale_color_manual(values = regioncol) +
          geom_point(data = targ.results2[targ.results2$group == 'observed',], aes(x = tissue, y = zscore,fill = Sig),shape = 23,col = 'black')+
          scale_fill_manual(values = c('Non-Sig'='grey','Border.sig' = '#2A9D8F',Sig='red')) +
          scale_shape_manual(values=c(23))+
          theme_bw()+
          coord_flip()+
          theme(text = element_text(size=8),
                #axis.ticks.y = element_blank(),
                #axis.ticks.x = element_blank(),
                legend.position = 'none',
                legend.title = element_blank(),
                axis.title.y = element_blank()) + 
          scale_x_discrete(position = "bottom")+
          scale_y_continuous(position = "left")+
          xlab('Cancer') + ylab('Overlapping Regions\n(Z-score Normalized)') 
        
        plot2 = ggplot(targ.results1[targ.results1$group == 'observed',], aes(x = tissue,y=n.overlap,fill = region)) + 
          geom_bar(stat= 'identity',position ='stack')+
          coord_flip()+
          scale_fill_manual(values = regioncol) +
          theme_bw()+
          theme(text = element_text(size=8),
                axis.ticks.y = element_blank(),
                axis.ticks.x = element_blank(),
                axis.title.y = element_blank(),
                axis.text.y = element_blank(),
                legend.position = 'none',
                legend.title = element_blank()) +
          scale_y_continuous(position = "left")+
          xlab('Cancer') + ylab('Overlapping Regions\n(Count)') 
        
        figure <- ggarrange(plot1, plot2,
                            labels = c(""),
                            ncol = 2, nrow = 1,
                            heights = c(0.5, 0.5),
                            widths = c(1.5,0.7),
                            align = 'h')
        
        file = gsub(' ','.',paste0(figdir,c,'.',dmr,'_', gsub(', ','',d),'.permutation.overlap.',region,'.raw.png'))
        png(file, height= 1250, width= 800,res=300)
        print(figure)
        dev.off()
        
        padjust.permutation = function(enhancer.permutation, vars = c('region','comparison')) {
          tmp.background = enhancer.permutation[enhancer.permutation$seed != 0,]
          tmp.background$p.adjust = 1
          tmp.observed = enhancer.permutation[enhancer.permutation$seed == 0,]
          tmp.observed.list = split(tmp.observed, tmp.observed[,vars])
          tmp.observed.list = lapply(tmp.observed.list, function(x) {
            return.df = x
            return.df$p.adjust = p.adjust(return.df$sig.p,method = 'bonferroni')
            return(return.df)
          })
          return.df = rbind(do.call('rbind',tmp.observed.list),tmp.background)
          return(return.df)
          
        }
        targ.results2 = padjust.permutation(targ.results2,c('direction'))
        targ.results2$sig.p = targ.results2$p.adjust
        targ.results2$Sig = ifelse(targ.results2$sig.p < 0.05, 'Sig','Non-Sig') #,ifelse(targ.results2$sig.p < 0.15,'Border.sig'
        
        plot1 = ggplot(targ.results2[targ.results2$group == 'permutation',], aes(x = tissue,y=zscore)) +
          geom_boxplot(outlier.shape = NA, aes(x = tissue,y=zscore),col = 'black') + 
          scale_color_manual(values = regioncol) +
          geom_point(data = targ.results2[targ.results2$group == 'observed',], aes(x = tissue, y = zscore,fill = Sig),shape = 23,col = 'black')+
          scale_fill_manual(values = c('Non-Sig'='grey','Border.sig' = '#2A9D8F',Sig='red')) +
          scale_shape_manual(values=c(23))+
          theme_bw()+
          coord_flip()+
          theme(text = element_text(size=8),
                #axis.ticks.y = element_blank(),
                #axis.ticks.x = element_blank(),
                axis.title.y = element_blank(), 
                legend.position = 'none',
                legend.title = element_blank()) + 
          scale_x_discrete(position = "bottom")+
          scale_y_continuous(position = "left")+
          xlab('Cancer') + ylab('Overlapping Regions\n(Z-score Normalized)') 
        
        plot2 = ggplot(targ.results1[targ.results1$group == 'observed',], aes(x = tissue,y=n.overlap,fill = region)) + 
          geom_bar(stat= 'identity',position ='stack')+
          coord_flip()+
          scale_fill_manual(values = regioncol) +
          theme_bw()+
          theme(text = element_text(size=8),
                axis.ticks.y = element_blank(),
                axis.ticks.x = element_blank(),
                axis.title.y = element_blank(),
                axis.text.y = element_blank(),
                legend.position = 'none',
                legend.title = element_blank()) +
          scale_y_continuous(position = "left")+
          xlab('Cancer') + ylab('Overlapping Regions\n(Count)') 
        
        figure <- ggarrange(plot1, plot2,
                            labels = c(""),
                            ncol = 2, nrow = 1,
                            heights = c(0.5, 0.5),
                            widths = c(1.5,0.7),
                            align = 'h')
        
        file = gsub(' ','.',paste0(figdir,c,'.',dmr,'_', gsub(', ','',d),'.permutation.overlap.',region,'.padj.png'))
        png(file, height= 1250, width= 800,res=300)
        print(figure)
        dev.off()
        
        
      }
      
    }
  }
  
  
  
  
  
  
}










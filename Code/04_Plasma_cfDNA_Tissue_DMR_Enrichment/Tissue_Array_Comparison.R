

cpg_count = readRDS('/.mounts/labs/awadallalab/private/ncheng/references/cpg_sites/cpg_site_positions/window_count/hg38_cpg_window_300_count.RDS') #number of cpg sites across 300bp regions
cpg_count$window = as.character(cpg_count$window)
background=cpg_count$window
cpg_count = cpg_count[cpg_count$count >= 5,] #selecting windows with at least 6 or mor CpG sites


cgi=readRDS('/.mounts/labs/awadallalab/private/ncheng/references/medip.300.references/hg38.cig.info.300.RDS')
#blood.wgbs = readRDS('/.mounts/labs/awadallalab/private/ncheng/vcfs/tissue_methylation/ihec/wgbs/wbc/bed/methylation_profile/hg38.liftover/summary_300/wbc.ihec.hg38.300.RDS')
blood.wgbs = readRDS('/.mounts/labs/awadallalab/private/ncheng/vcfs/tissue_methylation/ihec/wgbs/wbc/hg38.wbc.ihec.300.mean.RDS')
#breast.wgbs = readRDS('/.mounts/labs/awadallalab/ext_data/nxl.methylation/bed.files/hg38.liftover/300_summary/nxl.window300.mean.RDS')
blacklist = readRDS('/.mounts/labs/awadallalab/private/ncheng/references/encode.blacklist/hg38-blacklist.v2_window300.RDS')
#pbl.wgbs.information = readRDS("/.mounts/labs/awadallalab/private/ncheng/vcfs/tissue_methylation/ihec/wgbs/wbc/bed/methylation_profile/window_300_summary/ihec_wbc_sample_information.RDS")
blood.remove1 = unique(unlist(lapply(blood.wgbs[1:3], function(x) x[x$pct > 0.1 & x$cell_mean > 0.25,'window'] )))
gene.body= readRDS('/.mounts/labs/awadallalab/private/ncheng/references/gene_locations/hg38/gencode.window300.genes.RDS')
repeat.element = readRDS('/.mounts/labs/awadallalab/private/ncheng/references/repeat.elements/hg38.repeat.window.300.filter.elements.cut.RDS')

regulatory = readRDS('/.mounts/labs/awadallalab/private/ncheng/references/fantom5/hg38/hg38.window300.promoter.enhancer.RDS')
#cgi=readRDS('/.mounts/labs/awadallalab/private/ncheng/references/medip.300.references/hg38.cig.info.300.RDS')
combined.reg = readRDS('/.mounts/labs/awadallalab/private/ncheng/references/genomic.region.hg38.annot.window300.RDS')
hallmark.genes= readRDS('/.mounts/labs/awadallalab/private/ncheng/references/hallmark.genes.zapata/zapata2024.hallmark.genes.RDS')


targ.filt = function(x) {
  a = x[x %in% cpg_count$window]
  a = a[!a %in% blacklist$window]
  a = a[a %in% regulatory$window | a %in% cgi$window | a %in% repeat.element$window | a %in% gene.body$window]
  a =a[!grepl('chrX|chrY',a)]
  #a = a[a %in% feature.filt[['brca']]]
  #a = a[!a %in% feature.filt[['blood']]]
  return(a)
}

background = cpg_count$window
background.annotation = targ.filt(background)
#####prad####
library(ggpubr)
library(scales)
library(GenomicRanges)
#dds.matrix = readRDS(paste0('/.mounts/labs/awadallalab/private/ncheng/cfmedip_data/cptp_samples/fragmentation/aix13.updated1/all.sample.fragmentation/regulatory.counts/','/AIX13.Male.samples.deseq.normcounts.inserts.',1000,'.all.300.q20.RDS'))
savedir='/.mounts/labs/awadallalab/private/ncheng/cfmedip_data/cptp_samples/fragmentation/aix13.updated1/all.sample.fragmentation/prostate.cancer.cv/regulatory.regions.v5.ageadjusted/1/genhancer.1000/'
dmrdir='/.mounts/labs/awadallalab/private/ncheng/cfmedip_data/cptp_samples/fragmentation/aix13.updated1/all.sample.fragmentation/prostate.cancer.cv/regulatory.regions.v5.ageadjusted/1/'

cancer.dmr.dir = '/.mounts/labs/awadallalab/private/ncheng/vcfs/tcga_meth/cancer.pbl.rg/'
targ.cancer.dmr = readRDS(paste0(cancer.dmr.dir,'prostate','_wbc.dmrs.RDS'))

#prad.pbl.funnorm = readRDS('/.mounts/labs/awadallalab/private/ncheng/methylation_analysis/tissue_methylome_profiles/prostate_wbc/prostate_wbc_idat_filt_funnorm.RDS')
HM450.hg38.annotation = readRDS('/.mounts/labs/awadallalab/private/ncheng/annotation_files/meth_array/hg38/HM450.hg38.annotation.window300.RDS')

enhancer.promoter.array.annotation= function(targ.windows,HM450.hg38.annotation) {
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
  
  return(HM450.hg38.annotation.df) #reduce(filters.hg19)
}

#
prostate_cancer_blood_dmr =targ.cancer.dmr
prostate_cancer_blood_dmr$abs_difference = prostate_cancer_blood_dmr$cancer.mean - prostate_cancer_blood_dmr$blood.mean
prostate_cancer_blood_dmr$p.adjust = prostate_cancer_blood_dmr$cancer.blood.qval
prostate_cancer_normal_dmr = targ.cancer.dmr
prostate_cancer_normal_dmr$abs_difference = prostate_cancer_normal_dmr$cancer.mean - prostate_cancer_normal_dmr$adjnormal.mean
prostate_cancer_normal_dmr$p.adjust = prostate_cancer_normal_dmr$cancer.norm.qval
prostate_normal_blood_dmr =targ.cancer.dmr
prostate_normal_blood_dmr$abs_difference = prostate_normal_blood_dmr$adjnormal.mean - prostate_normal_blood_dmr$blood.mean
prostate_normal_blood_dmr$p.adjust = prostate_normal_blood_dmr$norm.blood.qval

setwd(figdir)
#prostate cancer vs normal
prostate_cancer_blood_dmr_sig = prostate_cancer_blood_dmr[as.numeric(prostate_cancer_blood_dmr$p.adjust) < 0.05 ,]
prostate_cancer_normal_dmr_sig = prostate_cancer_normal_dmr[as.numeric(prostate_cancer_normal_dmr$p.adjust) < 0.05,]
prostate_normal_blood_dmr_sig = prostate_normal_blood_dmr[as.numeric(prostate_normal_blood_dmr$p.adjust) < 0.05,]


tcga.dmr.dir= '/.mounts/labs/awadallalab/private/ncheng/methylation_analysis/tissue_methylome_profiles/breast_wbc/'
breast_cancer_blood_dmr = readRDS(paste0(tcga.dmr.dir,'breast_cancer_blood_dmr.RDS'))
breast_cancer_normal_dmr = readRDS(paste0(tcga.dmr.dir,'breast_cancer_normal_dmr.RDS'))
breast_normal_blood_dmr = readRDS(paste0(tcga.dmr.dir,'breast_normal_blood_dmr.RDS'))

'/.mounts/labs/awadallalab/private/ncheng/vcfs/tcga_meth/cancer.pbl.rg'

#breast cancer vs normal
breast_cancer_blood_dmr_sig = breast_cancer_blood_dmr[breast_cancer_blood_dmr$p.adjust < 0.05 ,]
breast_cancer_normal_dmr_sig = breast_cancer_normal_dmr[breast_cancer_normal_dmr$p.adjust < 0.05 ,]
#breast_blood_normal_dmr_sig = breast_normal_blood_dmr[breast_normal_blood_dmr$p.adjust < 0.05 & breast_normal_blood_dmr$abs_difference> 0.15,]

#grange sig overlap
tcga.dmr.list = list('Prostate' = list('Cancer.Blood' = prostate_cancer_blood_dmr_sig,
                                       'Cancer.Normal' = prostate_cancer_normal_dmr_sig),
                     'Breast' = list('Cancer.Blood' = breast_cancer_blood_dmr_sig,
                                     'Cancer.Normal' = breast_cancer_normal_dmr_sig))

#cfdna dmrs
cfdna.dmr.list=list()
for (sex in c('Male','Female')) {
  savedir='/.mounts/labs/awadallalab/private/ncheng/manuscripts/early.cancer.risk/all.sample.dmr.cpg5/'
  res.df = readRDS(paste0(savedir,'discovery.',sex,'.nonage.adjusted.dmrs.inserts.400.all.300.q20.RDS'))
  
  
  res.df = res.df[!grepl('chrX|chrY',res.df$window),]
  res.df = res.df[res.df$window %in% background.annotation,]
  res.df.hyper.sig.raw = res.df[res.df$log2FoldChange > 0.25 & res.df$pvalue < 0.05,]
  res.df.hypo.sig.raw = res.df[res.df$log2FoldChange < -0.25 & res.df$pvalue < 0.05,]
  res.df.hyper.sig.raw = res.df.hyper.sig.raw[order(res.df.hyper.sig.raw$pvalue),][1:2000,]
  res.df.hypo.sig.raw = res.df.hypo.sig.raw[order(res.df.hypo.sig.raw$pvalue),][1:2000,]
  cancer = ifelse(sex == 'Male','Prostate','Breast')
  cfdna.dmr.list[[cancer]] = list()
  cfdna.dmr.list[[cancer]][['Hypermethylated']] = res.df.hyper.sig.raw
  cfdna.dmr.list[[cancer]][['Hypomethylated']] = res.df.hypo.sig.raw
  
  
}

overlap.list = list()
for (sex in c('Breast','Prostate')) {
  hyper.tmp = cfdna.dmr.list[[sex]][[1]]
  hypo.tmp = cfdna.dmr.list[[sex]][[2]]
  
  hyper.overlap =HM450.hg38.annotation[HM450.hg38.annotation$window %in% hyper.tmp$window,]
  
  hypo.overlap =HM450.hg38.annotation[HM450.hg38.annotation$window %in% hypo.tmp$window,]
  overlap.list[[paste0(sex,'.Hypermethylated')]] = hyper.overlap
  overlap.list[[paste0(sex,'.Hypomethylated')]] = hypo.overlap
}


cpg.window.annotated = readRDS('/.mounts/labs/awadallalab/private/ncheng/references/gene.annotations/genomic.annotations.window.300.RDS')
cpg.window.annotated.regulatory = cpg.window.annotated[cpg.window.annotated$count >= 5,]
#cpg.window.annotated.regulatory = cpg.window.annotated.regulatory[cpg.window.annotated.regulatory$"Regulatory Element" != 'Non-Regulatory Element',]

cpg_count = readRDS('/.mounts/labs/awadallalab/private/ncheng/references/cpg_sites/cpg_site_positions/window_count/hg38_cpg_window_300_count.RDS') #number of cpg sites across 300bp regions
cpg_count$window = as.character(cpg_count$window)
background=cpg_count$window
cpg_count=cpg_count[!grepl('chrX|chrY',cpg_count$window),]

cpg_count = cpg_count[cpg_count$count >= 5,] #selecting windows with at least 6 or mor CpG sites


cpg.count.split = split(cpg_count,cpg_count$count)
#results.df.tmp.pradblood = z.score.permutation(prostate_cancer_blood_dmr_sig, predx.dmrs.sig.hyper,background,HM450.hg38.annotation)

results.df.tmp.pradblood = z.score.permutation(prostate_cancer_blood_dmr_sig, 
                                               predx.dmrs.sig.hyper,
                                               background = HM450.hg38.annotation) 
#new
z.score.permutation = function(prostate_cancer_blood_dmr_sig, predx.dmrs.sig.hyper,background, plot = F) {
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
    tcga.dmrs = background[background$probeID %in% prostate_cancer_blood_dmr_sig$probe,'window']
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
    overall.perf.list = mclapply(1:3000, function(x) {
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
#old

overlap.regions = function(prostate_cancer_blood_dmr_sig, predx.dmrs.sig.hyper,background, plot = F) {
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
    tcga.dmrs = background[background$probeID %in% prostate_cancer_blood_dmr_sig$probe,'window']
    #tcga.dmrs = unique(background.hm450[background.hm450 %in% tcga.dmrs])
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
overall.results.pradblood = NULL
overall.results.pradbrnm = NULL
overall.results.prnmblood=NULL
p.threshold =  c('0.05','0.01','0.001','0.0001')
p.threshold = c('0.05','0.01','0.001')
p.threshold = 0.05
effect = c(0.1)
direct = c('cfDNA Hyper, PRAD Hyper','cfDNA Hyper, PRAD Hypo','cfDNA Hypo, PRAD Hyper','cfDNA Hypo, PRAD Hypo')
#direct = c('cfDNA Hyper, PRAD Hyper','cfDNA Hypo, PRAD Hypo')

#breast
permutation.function(prostate_cancer_blood_dmr=cancer_blood_dmr;
                     prostate_cancer_normal_dmr=cancer_normal_dmr;
                     predx.dmrs.sig=sig.regions;
                     effect=0.1,
                     p.threshold=0.05,
                     direct=direct,
                     background=background.regions.array.450k)


permutation.function = function(prostate_cancer_blood_dmr,
                                prostate_cancer_normal_dmr,
                                #prostate_normal_blood_dmr ,
                                
                                predx.dmrs.sig,
                                effect,
                                p.threshold,
                                direct,
                                HM450.hg38.annotation) { 
  overall.results.pradblood = NULL
  overall.results.pradbrnm = NULL
  for (p in p.threshold) {
    for (j in effect) {
      for (d in direct){
        if (d == 'cfDNA Hyper, PRAD Hyper') {
          prostate_cancer_blood_dmr_sig = prostate_cancer_blood_dmr[as.numeric(prostate_cancer_blood_dmr$p.adjust) < p & prostate_cancer_blood_dmr$abs_difference> j,]
          prostate_cancer_normal_dmr_sig = prostate_cancer_normal_dmr[as.numeric(prostate_cancer_normal_dmr$p.adjust) < p & prostate_cancer_normal_dmr$abs_difference> j,]
          prostate_normal_blood_dmr_sig= prostate_normal_blood_dmr[as.numeric(prostate_normal_blood_dmr$p.adjust) <p & prostate_normal_blood_dmr$abs_difference > j,]
          predx.dmrs.sig.hyper = predx.dmrs.sig[predx.dmrs.sig$log2FoldChange > 0.25,]
          
        } else if (d == 'cfDNA Hyper, PRAD Hypo') {
          prostate_cancer_blood_dmr_sig = prostate_cancer_blood_dmr[as.numeric(prostate_cancer_blood_dmr$p.adjust) < p & prostate_cancer_blood_dmr$abs_difference< -j,]
          prostate_cancer_normal_dmr_sig = prostate_cancer_normal_dmr[as.numeric(prostate_cancer_normal_dmr$p.adjust) < p & prostate_cancer_normal_dmr$abs_difference< -j,]
          prostate_normal_blood_dmr_sig= prostate_normal_blood_dmr[as.numeric(prostate_normal_blood_dmr$p.adjust) <p & prostate_normal_blood_dmr$abs_difference < -j,]
          
          predx.dmrs.sig.hyper = predx.dmrs.sig[predx.dmrs.sig$log2FoldChange > 0.25,]
          
        }  else if (d == 'cfDNA Hypo, PRAD Hypo') {
          prostate_cancer_blood_dmr_sig = prostate_cancer_blood_dmr[as.numeric(prostate_cancer_blood_dmr$p.adjust) < p & prostate_cancer_blood_dmr$abs_difference< -j,]
          prostate_cancer_normal_dmr_sig = prostate_cancer_normal_dmr[as.numeric(prostate_cancer_normal_dmr$p.adjust) < p & prostate_cancer_normal_dmr$abs_difference< -j,]
          prostate_normal_blood_dmr_sig= prostate_normal_blood_dmr[as.numeric(prostate_normal_blood_dmr$p.adjust) <p & prostate_normal_blood_dmr$abs_difference < -j,]
          
          predx.dmrs.sig.hyper = predx.dmrs.sig[predx.dmrs.sig$log2FoldChange < -0.25,]
          
        }  else if (d == 'cfDNA Hypo, PRAD Hyper') {
          prostate_cancer_blood_dmr_sig = prostate_cancer_blood_dmr[as.numeric(prostate_cancer_blood_dmr$p.adjust) < p & prostate_cancer_blood_dmr$abs_difference> j,]
          prostate_cancer_normal_dmr_sig = prostate_cancer_normal_dmr[as.numeric(prostate_cancer_normal_dmr$p.adjust) < p & prostate_cancer_normal_dmr$abs_difference> j,]
          prostate_normal_blood_dmr_sig= prostate_normal_blood_dmr[as.numeric(prostate_normal_blood_dmr$p.adjust) <p & prostate_normal_blood_dmr$abs_difference > j,]
          
          predx.dmrs.sig.hyper = predx.dmrs.sig[predx.dmrs.sig$log2FoldChange < -0.25,]
          
        } 
        predx.dmrs.sig.hyper = predx.dmrs.sig.hyper[order(predx.dmrs.sig.hyper$pvalue),]
        
        results.df.tmp.pradblood = z.score.permutation(prostate_cancer_blood_dmr_sig, 
                                                       predx.dmrs.sig.hyper,
                                                       background = HM450.hg38.annotation) 
        results.df.tmp.pradblood.features = overlap.regions(prostate_cancer_blood_dmr_sig, 
                                                            predx.dmrs.sig.hyper,
                                                            background = HM450.hg38.annotation) 
        results.df.tmp.pradblood$pvalue = p
        results.df.tmp.pradblood$abs.diff = j
        results.df.tmp.pradblood$direction = d
        results.df.tmp.pradblood$comparison = 'PRAD.PBL'
        overall.results.pradblood = rbind(overall.results.pradblood,results.df.tmp.pradblood)
        
        results.df.tmp.pradbrnm = z.score.permutation(prostate_cancer_normal_dmr_sig, 
                                                      predx.dmrs.sig.hyper,
                                                      background = HM450.hg38.annotation) 
        results.df.tmp.pradbrnm$pvalue = p
        results.df.tmp.pradbrnm$abs.diff = j
        results.df.tmp.pradbrnm$direction = d
        results.df.tmp.pradbrnm$comparison = 'PRAD.PRNM'
        overall.results.pradbrnm = rbind(overall.results.pradbrnm,results.df.tmp.pradbrnm)
        
        
        # results.df.tmp.prnmblood = z.score.permutation(prostate_normal_blood_dmr_sig, 
        #                                                predx.dmrs.sig.hyper,
        #                                               background) 
        #results.df.tmp.prnmblood$pvalue = p
        #results.df.tmp.prnmblood$abs.diff = j
        #results.df.tmp.prnmblood$direction = d
        # results.df.tmp.prnmblood$comparison = 'PRNM.PBL'
        
        # overall.results.prnmblood = rbind(overall.results.prnmblood,results.df.tmp.prnmblood)
        
        
        
      }
      
    }
    
  }
  return(list(prad.blood=overall.results.pradblood, prad.prnm=overall.results.pradbrnm))
  #,prnm.blood = overall.results.prnmblood
  
}

feature.overlaps = function(prostate_cancer_blood_dmr,
                            prostate_cancer_normal_dmr,
                            #prostate_normal_blood_dmr ,
                            
                            predx.dmrs.sig,
                            effect,
                            p.threshold,
                            direct,
                            background) 
{ 
  overall.results.pradblood = NULL
  overall.results.pradbrnm = NULL
  for (p in p.threshold) {
    for (j in effect) {
      for (d in direct){
        if (d == 'cfDNA Hyper, PRAD Hyper') {
          prostate_cancer_blood_dmr_sig = prostate_cancer_blood_dmr[as.numeric(prostate_cancer_blood_dmr$p.adjust) < p & prostate_cancer_blood_dmr$abs_difference> j,]
          prostate_cancer_normal_dmr_sig = prostate_cancer_normal_dmr[as.numeric(prostate_cancer_normal_dmr$p.adjust) < p & prostate_cancer_normal_dmr$abs_difference> j,]
          prostate_normal_blood_dmr_sig= prostate_normal_blood_dmr[as.numeric(prostate_normal_blood_dmr$p.adjust) <p & prostate_normal_blood_dmr$abs_difference > j,]
          predx.dmrs.sig.hyper = predx.dmrs.sig[predx.dmrs.sig$log2FoldChange > 0.25,]
          
        } else if (d == 'cfDNA Hyper, PRAD Hypo') {
          prostate_cancer_blood_dmr_sig = prostate_cancer_blood_dmr[as.numeric(prostate_cancer_blood_dmr$p.adjust) < p & prostate_cancer_blood_dmr$abs_difference< -j,]
          prostate_cancer_normal_dmr_sig = prostate_cancer_normal_dmr[as.numeric(prostate_cancer_normal_dmr$p.adjust) < p & prostate_cancer_normal_dmr$abs_difference< -j,]
          prostate_normal_blood_dmr_sig= prostate_normal_blood_dmr[as.numeric(prostate_normal_blood_dmr$p.adjust) <p & prostate_normal_blood_dmr$abs_difference < -j,]
          
          predx.dmrs.sig.hyper = predx.dmrs.sig[predx.dmrs.sig$log2FoldChange > 0.25,]
          
        }  else if (d == 'cfDNA Hypo, PRAD Hypo') {
          prostate_cancer_blood_dmr_sig = prostate_cancer_blood_dmr[as.numeric(prostate_cancer_blood_dmr$p.adjust) < p & prostate_cancer_blood_dmr$abs_difference< -j,]
          prostate_cancer_normal_dmr_sig = prostate_cancer_normal_dmr[as.numeric(prostate_cancer_normal_dmr$p.adjust) < p & prostate_cancer_normal_dmr$abs_difference< -j,]
          prostate_normal_blood_dmr_sig= prostate_normal_blood_dmr[as.numeric(prostate_normal_blood_dmr$p.adjust) <p & prostate_normal_blood_dmr$abs_difference < -j,]
          
          predx.dmrs.sig.hyper = predx.dmrs.sig[predx.dmrs.sig$log2FoldChange < -0.25,]
          
        }  else if (d == 'cfDNA Hypo, PRAD Hyper') {
          prostate_cancer_blood_dmr_sig = prostate_cancer_blood_dmr[as.numeric(prostate_cancer_blood_dmr$p.adjust) < p & prostate_cancer_blood_dmr$abs_difference> j,]
          prostate_cancer_normal_dmr_sig = prostate_cancer_normal_dmr[as.numeric(prostate_cancer_normal_dmr$p.adjust) < p & prostate_cancer_normal_dmr$abs_difference> j,]
          prostate_normal_blood_dmr_sig= prostate_normal_blood_dmr[as.numeric(prostate_normal_blood_dmr$p.adjust) <p & prostate_normal_blood_dmr$abs_difference > j,]
          
          predx.dmrs.sig.hyper = predx.dmrs.sig[predx.dmrs.sig$log2FoldChange < -0.25,]
          
        } 
        predx.dmrs.sig.hyper = predx.dmrs.sig.hyper[order(predx.dmrs.sig.hyper$pvalue),]
        
        
        results.df.tmp.pradblood = overlap.regions(prostate_cancer_blood_dmr_sig, 
                                                   predx.dmrs.sig.hyper,
                                                   background = background) 
        results.df.tmp.pradblood$pvalue = p
        results.df.tmp.pradblood$abs.diff = j
        results.df.tmp.pradblood$direction = d
        results.df.tmp.pradblood$comparison = 'PRAD.PBL'
        overall.results.pradblood = rbind(overall.results.pradblood,results.df.tmp.pradblood)
        
        results.df.tmp.pradbrnm = overlap.regions(prostate_cancer_normal_dmr_sig, 
                                                  predx.dmrs.sig.hyper,
                                                  background) 
        results.df.tmp.pradbrnm$pvalue = p
        results.df.tmp.pradbrnm$abs.diff = j
        results.df.tmp.pradbrnm$direction = d
        results.df.tmp.pradbrnm$comparison = 'PRAD.PRNM'
        overall.results.pradbrnm = rbind(overall.results.pradbrnm,results.df.tmp.pradbrnm)
        
        
        #   results.df.tmp.prnmblood = overlap.regions(prostate_normal_blood_dmr_sig, 
        #                                             predx.dmrs.sig.hyper,
        #                                             background) 
        #  results.df.tmp.prnmblood$pvalue = p
        #  results.df.tmp.prnmblood$abs.diff = j
        #  results.df.tmp.prnmblood$direction = d
        #  results.df.tmp.prnmblood$comparison = 'PRNM.PBL'
        
        # overall.results.prnmblood = rbind(overall.results.prnmblood,results.df.tmp.prnmblood)
        
        
        
      }
      
    }
    
  }
  return(list(prad.blood=overall.results.pradblood, prad.prnm=overall.results.pradbrnm))
  
}




permutation.plot = function(p.threshold=0.05, effect=0.1,name,enhancer.permutation,sex) {
  p=p.threshold  
  j=effect
  d ='cfDNA Hyper, PRAD Hyper'
  n=0
  dmr.overlap.list = split(enhancer.permutation,enhancer.permutation$marker)
  for (dmr.targ in names(dmr.overlap.list)){
    for (d in c('cfDNA Hyper, PRAD Hyper','cfDNA Hypo, PRAD Hypo')) {
      d.name = gsub(' ','',gsub(',','',d))
      targ.df = dmr.overlap.list[[dmr.targ]]
      targ.results = targ.df[targ.df$direction == d,]#targ.df$pvalue == as.character(p) & targ.df$abs.diff == j & 
      targ.results$region =factor(targ.results$region, levels = c('All','Islands','Shores','Shelves','Open Sea'))
      targ.results$Sig = ifelse(targ.results$sig.p < 0.05, 'Sig','Non-Sig')
      regioncol = c('Islands' ='#29335C','Shores' = '#F3A712','Shelves' = '#DB2B39','Open Sea' = '#A89B83','All' = 'black')
      
      plot1 = ggplot(targ.results[targ.results$group == 'permutation',], aes(x = region,y=zscore,col= region)) +
        geom_boxplot(outlier.shape = NA, aes(x = region,y=zscore,color = region)) + 
        scale_color_manual(values = regioncol) +
        
        geom_point(data = targ.results[targ.results$group == 'observed',],mapping = aes(x = region, y = zscore,fill = Sig),shape = 23, col = 'black')+
        scale_fill_manual(values = c('Non-Sig'='grey',Sig='red')) +
        
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
      
      plot2 = ggplot(targ.results[targ.results$group == 'observed',], aes(x = region,y=n.overlap,fill = region)) + 
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
      
      png(paste0(name,'.permutation.overlap.',dmr.targ,'.',d.name,'.all','.png'), height= 1500, width= 1200,res=250)
      print(figure)
      dev.off()
      
      
    }
    d.name = gsub(' ','',gsub(',','',d))
    targ.df = dmr.overlap.list[[dmr.targ]]
    targ.results = targ.df#[targ.df$direction == d,]#targ.df$pvalue == as.character(p) & targ.df$abs.diff == j & 
    
    targ.results$region =factor(targ.results$region, levels = c('All','Islands','Shores','Shelves','Open Sea'))
    targ.results$Sig = ifelse(targ.results$sig.p < 0.05, 'Sig','Non-Sig')
    regioncol = c('Islands' ='#29335C','Shores' = '#F3A712','Shelves' = '#DB2B39','Open Sea' = '#A89B83','All' = 'black')
    for(region in  c('All','Islands','Shores','Shelves','Open Sea')) {
      targ.results1 = targ.results[targ.results$region == region,]
      targ.results1$direction = gsub(', ','\n',targ.results1$direction)
      if (sex == "Female") {
        targ.results1$direction = gsub('PRAD','BRCA',targ.results1$direction)
        
      }
      plot1 = ggplot(targ.results1[targ.results1$group == 'permutation',], aes(x = direction,y=zscore)) +
        geom_boxplot(outlier.shape = NA, aes(x = direction,y=zscore)) + 
        scale_color_manual(values = regioncol) +
        
        geom_point(data = targ.results1[targ.results1$group == 'observed',],mapping = aes(x = direction, y = zscore,fill = Sig),shape = 23,col = 'black')+
        scale_fill_manual(values = c('grey','red')) +
        scale_shape_manual(values=c(23))+
        theme_bw()+
        theme(text = element_text(size=14),
              axis.ticks.y = element_blank(),
              axis.ticks.x = element_blank(),
              legend.position = 'none',
              legend.title = element_blank()) + 
        xlab('Region') +ylab('Overlapping DMRs\n(z-score normalized)') 
      
      if (region == 'All') {
        plot2 = ggplot(targ.results[targ.results$group == 'observed' & targ.results$region != 'All',], aes(x = direction,y=n.overlap, fill = region)) + 
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
        plot2 = ggplot(targ.results1[targ.results1$group == 'observed',], aes(x = direction,y=n.overlap, fill = region)) + 
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
      
      
      png(paste0(name,'.permutation.overlap.',dmr.targ,'.',region,'.png'), height= 1500, width= 1200,res=250)
      print(figure)
      dev.off()
      
    }
    
    
  }
  options(scipen=99990)
  
  
  
}

permutation.plot.manuscript = function(p.threshold=0.05, effect=0.1,name,enhancer.permutation,sex) {
  p=p.threshold  
  j=effect
  d ='cfDNA Hyper, PRAD Hyper'
  n=0
  dmr.overlap.list = split(enhancer.permutation,enhancer.permutation$marker)
  for (dmr.targ in names(dmr.overlap.list)){
    for (d in c('cfDNA Hyper, PRAD Hyper','cfDNA Hypo, PRAD Hypo')) {
      d.name = gsub(' ','',gsub(',','',d))
      targ.df = dmr.overlap.list[[dmr.targ]]
      targ.results = targ.df[targ.df$direction == d,]#targ.df$pvalue == as.character(p) & targ.df$abs.diff == j & 
      #targ.results$region = gsub(' ','\n',targ.results$region)
      targ.results$region =factor(targ.results$region, levels = c('All','Islands','Shores','Shelves','Open Sea'))
      targ.results$Sig = ifelse(targ.results$sig.p < 0.05, 'Sig',ifelse(targ.results$sig.p < 0.15,'Border.sig','Non-Sig'))
      regioncol = c('Islands' ='#29335C','Shores' = '#F3A712','Shelves' = '#DB2B39','Open Sea' = '#A89B83','All' = 'black')
      plot1 = ggplot(targ.results[targ.results$group == 'permutation' & targ.results$region != 'All',], aes(x = region,y=zscore,col= region)) +
        geom_boxplot(outlier.shape = NA, aes(x = region,y=zscore,color = region)) + 
        scale_color_manual(values = regioncol) +
        
        geom_point(data = targ.results[targ.results$group == 'observed' & targ.results$region != 'All',],mapping = aes(x = region, y = zscore,fill = Sig),shape = 23, col = 'black')+
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
      
      plot2 = ggplot(targ.results[targ.results$group == 'observed' &targ.results$region != 'All',], aes(x = region,y=n.overlap,fill = region)) + 
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
      
      png(paste0(name,'.permutation.overlap.',dmr.targ,'.',d.name,'.all','.png'), height= 1500, width= 700,res=250)
      print(figure)
      dev.off()
      
      
    }
    targ.df = dmr.overlap.list[[dmr.targ]]
    targ.results = targ.df#[targ.df$direction == d,]#targ.df$pvalue == as.character(p) & targ.df$abs.diff == j & 
    #targ.results$region = gsub(' ','\n',targ.results$region)
    
    targ.results$region =factor(targ.results$region, levels = c('All','Islands','Shores','Shelves','Open Sea'))
    targ.results$Sig = ifelse(targ.results$sig.p < 0.05, 'Sig',ifelse(targ.results$sig.p < 0.15,'Border.sig','Non-Sig'))
    
    regioncol = c('Islands' ='#29335C','Shores' = '#F3A712','Shelves' = '#DB2B39','Open Sea' = '#A89B83','All' = 'black')
    
    for(region in  c('All','Islands','Shores','Shelves','Open Sea')) {
      targ.results1 = targ.results[targ.results$region == region,]
      targ.results1$direction = gsub(', ','\n',targ.results1$direction)
      #targ.results1$region = gsub(' ','\n',targ.results1$region)
      
      if (sex == "Female") {
        targ.results1$direction = gsub('PRAD','BRCA',targ.results1$direction)
        
      }
      
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
  options(scipen=99990)
  
  
  
}

combined.sample.permutation=NULL
HM450.hg38.annotation = readRDS('/.mounts/labs/awadallalab/private/ncheng/annotation_files/meth_array/hg38/HM450.hg38.annotation.window300.RDS')
HM450.hg38.annotation.base=HM450.hg38.annotation
for (marker in names(cfdna.dmr.list)) {
  
  for (r in names(cfdna.dmr.list[[marker]])) {
    
    res.df= cfdna.dmr.list[[marker]][[r]]
    #res.df=readRDS('/.mounts/labs/awadallalab/private/ncheng/cfmedip_data/cptp_samples/fragmentation/aix13.updated1/all.sample.fragmentation/prostate.cancer.cv/regulatory.regions.v5.ageadjusted/1/validation.figures/CV_Highest_Average_Coefficient/silencer/silencer.top100.features.dmrs.RDS')
    res.df$window =rownames(res.df)
    fc=0.25
    res.df = res.df[!grepl('chrX.*|chrY.*',res.df$window),]
    background.regions  = rownames(res.df)
    sig.regions = res.df[which(res.df$pvalue < 0.05 & abs(res.df$log2FoldChange) > fc & res.df$baseMean > 1),] # 
    
    sig.regions.array.450k = enhancer.promoter.array.annotation(sig.regions$window,HM450.hg38.annotation.base)
    background.regions.array.450k = enhancer.promoter.array.annotation(cpg_count$window,HM450.hg38.annotation.base)
    
    
    cancer_blood_dmr= tcga.dmr.list[[marker]][['Cancer.Blood']]
    cancer_normal_dmr= tcga.dmr.list[[marker]][['Cancer.Normal']]
    
    direct = c('cfDNA Hyper, PRAD Hyper','cfDNA Hyper, PRAD Hypo','cfDNA Hypo, PRAD Hyper','cfDNA Hypo, PRAD Hypo')
    
    if(r == 'Hypermethylated') {
      direct =c('cfDNA Hyper, PRAD Hyper','cfDNA Hyper, PRAD Hypo')
      
    } else {
      direct= c('cfDNA Hypo, PRAD Hyper','cfDNA Hypo, PRAD Hypo')
      
    }
    
    silencer.permutation = permutation.function(cancer_blood_dmr,
                                                cancer_normal_dmr,
                                                predx.dmrs.sig=sig.regions,
                                                effect=c(0.1),
                                                p.threshold=c(0.05,0.01),
                                                direct=direct,
                                                HM450.hg38.annotation=background.regions.array.450k)
    
    
    
    silencer.features = feature.overlaps(prostate_cancer_blood_dmr,
                                         prostate_cancer_normal_dmr,
                                         predx.dmrs.sig=sig.regions,
                                         effect=c(0.1),
                                         p.threshold=c(0.05,0.01),
                                         direct=direct,
                                         background=background.regions.array.450k)
    
    
    figdir='/.mounts/labs/awadallalab/private/ncheng/cfmedip_data/cptp_samples/fragmentation/aix13.updated1/all.sample.fragmentation/prostate.cancer.cv/regulatory.regions.v5.ageadjusted/1/marker.tissue/target.overlap/tcga/'
    name=paste0(figdir,marker)
    #permutation.plot(0.05,0.15,name,silencer.permutation)
    silencer.permutation.df= do.call('rbind',silencer.permutation)
    silencer.permutation.df$marker = marker 
    combined.sample.permutation=rbind(combined.sample.permutation,silencer.permutation.df)
    
    
    
    
  }
  
}

#saveRDS(combined.sample.permutation,paste0(savedir,'all.sample.tcga.permutation.RDS'))
saveRDS(combined.sample.permutation,paste0(savedir,'all.sample.tcga.permutation.fc.new25.RDS'))

#saveRDS(combined.sample.permutation,paste0(savedir,'all.sample.tcga.permutation.fc.old25.RDS'))


# observed.overlap.count

savedir='/.mounts/labs/awadallalab/private/ncheng/manuscripts/early.cancer.risk/all.sample.dmr.cpg5/'

figdir=paste0(savedir,'figures.fc0.25.nofilt/adjusted/')
dir.create(figdir,recursive = T)

comparisons = c('PRAD.PBL','PRAD.PRNM')
cancer = c('Prostate','Breast')
pvalue = c(0.05,0.01)
effect = c(0.1)


for (p in pvalue ) {
  for (e in effect) {
    for (comp in comparisons) {
      for (c in cancer) {
        enhancer.permutation = combined.sample.permutation[combined.sample.permutation$marker == c & 
                                                             combined.sample.permutation$comparison == comp &
                                                             combined.sample.permutation$pvalue == p & 
                                                             combined.sample.permutation$abs.diff == e,]
        
        if (nrow(enhancer.permutation) > 0) {
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
          
          enhancer.permutation = padjust.permutation(enhancer.permutation)
          enhancer.permutation$sig.p = enhancer.permutation$p.adjust
          permutation.plot.manuscript(p.threshold = 0.05,
                                      effect = 0.1, 
                                      name = paste0(figdir,comp,'.',c,'.',p,'.',e), 
                                      enhancer.permutation,
                                      sex = ifelse(c == 'Prostate','Male','Female'))
          
          
        }
        
        
      }
    }
    
    
  }
}


results.df.tmp = overlap.regions(prostate_cancer_blood_dmr_sig, predx.dmrs.sig.hyper,background)

if (length(results.df.tmp) > 0) {
  results.df.tmp$comparison=i
  results.df.tmp$direction = d
  results.df.tmp$pvalue = p
  results.df.tmp$abs.diff = j
  overlappingmarkers.breast = rbind(overlappingmarkers.breast,results.df.tmp)
  
}


#plotting overlapping featuers
directions = c('cfDNA Hyper, PRAD Hyper','cfDNA Hypo, PRAD Hypo')
cancers = c('Breast','Prostate')

for (m in cancers) {
  for (d in directions){
    overlappingmarkers.prostate.hyper  =combined.sample.permutation[combined.sample.permutation$direction == d & 
                                                                      combined.sample.permutation$marker == m,]
    tmp1 = split(overlappingmarkers.prostate.hyper,overlappingmarkers.prostate.hyper$comparison)
    permutation.plot(0.05, effect,name,enhancer.permutation)
    
    overlap.hyper = HM450.hg38.annotation[,c('window','CGIposition','gene')]
    overlap.hyper$"Prostate Cancer vs Adjacent Prostate Normal" = ifelse(overlap.hyper$window %in% tmp1[["PRAD.PRNM"]]$window,1,0)
    overlap.hyper$"Prostate Cancer vs PBL"= ifelse(overlap.hyper$window %in% tmp1[["PRAD.PBL"]]$window,1,0)
    overlap.hyper$cpg_region = gsub('N_','',gsub('S_','',overlap.hyper$CGIposition))
    overlap.hyper[is.na(overlap.hyper$cpg_region),'cpg_region'] = 'null'
    overlap.hyper = unique(overlap.hyper[,c('window','cpg_region','Adjacent Prostate Normal vs PBL','Prostate Cancer vs PBL')])
    png(paste0(figdir,m,'.tissue.blood.',d,'.overlap.png'),height = 800, width =1000,res=300)
    print(upset(overlap.hyper,nsets = 9,keep.order = TRUE, sets = c("Adjacent Prostate Normal vs PBL" ,"Prostate Cancer vs PBL"),
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


#####plottting pancan tcga overlap#####
pancancer.permutation = function(predx.dmrs.sig,background){
  cancer.groups = c('prostate','breast','brain','skin','lung','colorectal','pancreatic','bladder','liver','kidney','oesophageal', 'uterus','thyroid','headneck')
  overall.results.pancanblood = NULL
  overall.resultspancannorm = NULL
  overall.resultsnormblood =NULL
  
  
  #all
  for (c in cancer.groups) {
    cancer.dmr.dir = '/.mounts/labs/awadallalab/private/ncheng/vcfs/tcga_meth/cancer.pbl.rg/'
    targ.cancer.dmr = readRDS(paste0(cancer.dmr.dir,c,'_wbc.dmrs.RDS'))
    targ.cancer.dmr$cancer.blood.abs.diff = targ.cancer.dmr$cancer.mean -targ.cancer.dmr$blood.mean
    targ.cancer.dmr$cancer.norm.abs.diff = targ.cancer.dmr$cancer.mean -targ.cancer.dmr$adjnormal.mean
    targ.cancer.dmr$norm.blood.abs.diff =targ.cancer.dmr$adjnormal.mean -targ.cancer.dmr$blood.mean
    d1=c('cfDNA Hyper, PRAD Hyper','cfDNA Hypo, PRAD Hypo')
    
    print(c)
    for (p in 0.05) {
      for (j in 0.1) {
        for (d in d1){
          if (d == 'cfDNA Hyper, PRAD Hyper') {
            prostate_cancer_blood_dmr_sig = targ.cancer.dmr[targ.cancer.dmr$cancer.blood.qval < p & targ.cancer.dmr$cancer.blood.abs.diff> j,]
            prostate_cancer_normal_dmr_sig = targ.cancer.dmr[targ.cancer.dmr$cancer.norm.qval < p & targ.cancer.dmr$cancer.norm.abs.diff> j,]
            prostate_normal_blood_dmr_sig = targ.cancer.dmr[targ.cancer.dmr$norm.blood.qval < p & targ.cancer.dmr$norm.blood.abs.diff> j,]
            
            predx.dmrs.sig.hyper = predx.dmrs.sig[predx.dmrs.sig$log2FoldChange > 0.25,]
            
          } else if (d == 'cfDNA Hyper, PRAD Hypo') {
            prostate_cancer_blood_dmr_sig = targ.cancer.dmr[targ.cancer.dmr$cancer.blood.qval < p & targ.cancer.dmr$cancer.blood.abs.diff< -j,]
            prostate_cancer_normal_dmr_sig = targ.cancer.dmr[targ.cancer.dmr$cancer.norm.qval < p & targ.cancer.dmr$cancer.norm.abs.diff< -j,]
            prostate_normal_blood_dmr_sig = targ.cancer.dmr[targ.cancer.dmr$norm.blood.qval < p & targ.cancer.dmr$norm.blood.abs.diff< -j,]
            
            predx.dmrs.sig.hyper = predx.dmrs.sig[predx.dmrs.sig$log2FoldChange > 0.25,]
            
          }  else if (d == 'cfDNA Hypo, PRAD Hypo') {
            prostate_cancer_blood_dmr_sig = targ.cancer.dmr[targ.cancer.dmr$cancer.blood.qval < p & targ.cancer.dmr$cancer.blood.abs.diff< -j,]
            prostate_cancer_normal_dmr_sig = targ.cancer.dmr[targ.cancer.dmr$cancer.norm.qval < p & targ.cancer.dmr$cancer.norm.abs.diff< -j,]
            prostate_normal_blood_dmr_sig = targ.cancer.dmr[targ.cancer.dmr$norm.blood.qval < p & targ.cancer.dmr$norm.blood.abs.diff< -j,]
            
            
            predx.dmrs.sig.hyper = predx.dmrs.sig[predx.dmrs.sig$log2FoldChange < -0.25,]
            
          }  else if (d == 'cfDNA Hypo, PRAD Hyper') {
            prostate_cancer_blood_dmr_sig = targ.cancer.dmr[targ.cancer.dmr$cancer.blood.qval < p & targ.cancer.dmr$cancer.blood.abs.diff> j,]
            prostate_cancer_normal_dmr_sig = targ.cancer.dmr[targ.cancer.dmr$cancer.norm.qval < p & targ.cancer.dmr$cancer.norm.abs.diff> j,]
            prostate_normal_blood_dmr_sig = targ.cancer.dmr[targ.cancer.dmr$norm.blood.qval < p & targ.cancer.dmr$norm.blood.abs.diff> j,]
            
            predx.dmrs.sig.hyper = predx.dmrs.sig[predx.dmrs.sig$log2FoldChange < -0.25,]
            
          } 
          predx.dmrs.sig.hyper = predx.dmrs.sig.hyper[order(predx.dmrs.sig.hyper$pvalue),]
          #predx.dmrs.sig.hyper = predx.dmrs.sig.hyper[1:400,]
          permutation.calc=T
          if (permutation.calc == T){
            results.df.tmp.pradblood = z.score.permutation(prostate_cancer_blood_dmr_sig,
                                                           predx.dmrs.sig.hyper,
                                                           background) 
            results.df.tmp.pradblood$pvalue = p
            results.df.tmp.pradblood$abs.diff = j
            results.df.tmp.pradblood$direction = gsub('BRCA','Tissue',d)
            results.df.tmp.pradblood$tissue = c
            
            overall.results.pancanblood = rbind(overall.results.pancanblood,results.df.tmp.pradblood)
            
            
            results.df.tmp.pradbrnm = z.score.permutation(prostate_cancer_normal_dmr_sig, predx.dmrs.sig.hyper,background) 
            results.df.tmp.pradbrnm$pvalue = p
            results.df.tmp.pradbrnm$abs.diff = j
            results.df.tmp.pradbrnm$direction = gsub('BRCA','Tissue',d)
            results.df.tmp.pradbrnm$tissue =c 
            
            overall.resultspancannorm = rbind(overall.resultspancannorm,results.df.tmp.pradbrnm)
            
            results.df.tmp.prnmblood = z.score.permutation(prostate_normal_blood_dmr_sig, predx.dmrs.sig.hyper,background) 
            results.df.tmp.prnmblood$pvalue = p
            results.df.tmp.prnmblood$abs.diff = j
            results.df.tmp.prnmblood$direction = gsub('BRCA','Tissue',d)
            results.df.tmp.prnmblood$tissue =c 
            
            overall.resultsnormblood = rbind(overall.resultsnormblood,results.df.tmp.prnmblood)
            
          }
          
          
          
        }
        
      }
      
    }
    
  }
  return(list('PanCan.Blood'=overall.results.pancanblood,
              'PanCan.Norm'=overall.resultspancannorm,
              'Norm.Blood'=overall.resultsnormblood))
  #overlap onlyy
}

combined.pancan.permutation = list()
for (marker in names(cfdna.dmr.list)) {
  predx.dmrs.sig = do.call('rbind',cfdna.dmr.list[[marker]])
  a = pancancer.permutation(predx.dmrs.sig, background=HM450.hg38.annotation)
  a$marker = marker
  combined.pancan.permutation[[marker]] = a
  
}
saveRDS(combined.pancan.permutation,paste0(savedir,'combined.pancan.permutation.RDS'))


#
savedir='/.mounts/labs/awadallalab/private/ncheng/manuscripts/early.cancer.risk/all.sample.dmr.cpg5/'

figdir=paste0(savedir,'figures.fc0.25.nofilt/pancan/')
dir.create(figdir,recursive = T)


for (c in c('Breast','Prostate')) {
  targ.permutation.list = combined.pancan.permutation[[c]]
  comparison.list = list('pancan.blood' = targ.permutation.list[[1]],
                         'norm.blood'=targ.permutation.list[[3]])
  
  #plotting upset plot
  upset = F
  if (upset == T){
    for (p in names(comparison.list)){
      targ.features = comparison.list[[p]]
      overlappingmarkers.prostate.hyper  =targ.features#[targ.features$direction == d,]
      overlappingmarkers.prostate.hyper$cancer.targ = ifelse(overlappingmarkers.prostate.hyper$tissue == c,c,paste0('not.',c))
      tmp1 = split(overlappingmarkers.prostate.hyper,overlappingmarkers.prostate.hyper$cancer.targ)
      overlap.hyper = HM450.hg38.annotation[,c('window','CGIposition','gene')]
      overlap.hyper$"Prostate Tissue" = ifelse(overlap.hyper$window %in% tmp1[[c]]$window,1,0)
      overlap.hyper$"Non-Prostate Tissue"= ifelse(overlap.hyper$window %in% tmp1[[paste0('not.',c)]]$window,1,0)
      overlap.hyper$cpg_region = gsub('N_','',gsub('S_','',overlap.hyper$CGIposition))
      overlap.hyper[is.na(overlap.hyper$cpg_region),'cpg_region'] = 'null'
      overlap.hyper = unique(overlap.hyper[,c('window','cpg_region','Prostate Tissue','Non-Prostate Tissue')])
      plot.overlap = T
      png(paste0(figdir,'prostate.tissue.blood.',p,'.overlap.png'),height = 800, width =1000,res=300)
      print(upset(overlap.hyper,nsets = 9,keep.order = TRUE, sets = c('Prostate Tissue','Non-Prostate Tissue'),
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
      ad
      
      targ.prostate = overlappingmarkers.prostate.hyper[overlappingmarkers.prostate.hyper$tissue == 'prostate',]
      nonprostate = overlappingmarkers.prostate.hyper[overlappingmarkers.prostate.hyper$tissue != 'prostate',]
      nonprostate = nonprostate[!nonprostate$tissue %in% c('kidney','breast'),]
      
      prostate.count.overlap = NULL
      for (i in unique(targ.prostate$window)) {
        for (j in c('All','Open Sea','Islands','Shores','Shelves')) {
          nonprostate.filt = unique(nonprostate[nonprostate$windows %in% i & nonprostate$region %in% j,])
          targ.prostate.filt = targ.prostate[targ.prostate$windows %in% i & targ.prostate$region %in% j,]
          if (nrow(targ.prostate.filt) > 0) {
            targ.prostate.filt$non.prostate.count = nrow(nonprostate.filt)
            
            prostate.count.overlap = rbind(prostate.count.overlap,targ.prostate.filt)
            #
            
          }
          
        }
      }
      
      png(paste0(figdir,'pancan.',p,,'.',c,'.overlap.count.png'),height = 800, width = 1500,res=300)
      plot1 = ggplot(prostate.count.overlap[prostate.count.overlap$region != 'All',], aes(x = non.prostate.count,fill =region)) +  #Diagnosis_Time
        geom_bar()+
        theme_bw()+
        scale_x_continuous(breaks = seq(0,20,2))+
        theme(text = element_text(size=10),
              axis.ticks.y = element_blank(),
              axis.ticks.x = element_blank(),
              legend.position = 'none',
              legend.title = element_blank()) +
        xlab('') + ylab('')+ scale_fill_manual(values =c(Islands='#29335C','Open Sea'="#A89B83", 'Shelves'="#DB2B39",Shores= '#F3A712')) 
      # geom_hline(yintercept=0) +
      # geom_vline(xintercept=0)# + geom_smooth(method = "lm")
      print(plot1)
      dev.off()
      
      
    }
    
    
  }
  
  #plotting hyper-hypo overlap all#
  library(Hmisc)
  options(scipen=99990)
  p=0.05     
  j=0.15
  dmr.group=c('cfDNA Hyper, Tissue Hyper','cfDNA Hypo, Tissue Hypo')
  
  dmrs = list(pancanblood = targ.permutation.list[[1]], pancannorm=targ.permutation.list[[2]], normblood = targ.permutation.list[[3]])
  
  for (dmr in names(dmrs)) {
    targ.overlap = dmrs[[dmr]]
    if (c == 'Prostate') {
      targ.results = targ.overlap[!targ.overlap$tissue %in% c('kidney','breast'),]
      
      
    } else {
      targ.results = targ.overlap#[!targ.overlap$tissue %in% c('kidney','breast'),]
      
      
    }
    targ.results$region =factor(targ.results$region, levels = c('All','Islands','Shores','Shelves','Open Sea'))
    #targ.results$Sig = ifelse(targ.results$sig.p < 0.01, 'Sig','Non-Sig')
    targ.results$Sig = ifelse(targ.results$sig.p < 0.05, 'Sig',ifelse(targ.results$sig.p < 0.15,'Border.sig','Non-Sig'))
    
    
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
          theme(text = element_text(size=11),
                #axis.ticks.y = element_blank(),
                #axis.ticks.x = element_blank(),
                legend.position = 'none',
                legend.title = element_blank()) + 
          scale_x_discrete(position = "bottom")+
          scale_y_continuous(position = "right")+
          xlab('Cancer') + ylab('Overlapping Regions\n(Z-score Normalized)') 
        
        plot2 = ggplot(targ.results1[targ.results1$group == 'observed',], aes(x = tissue,y=n.overlap,fill = region)) + 
          geom_bar(stat= 'identity',position ='stack')+
          coord_flip()+
          scale_fill_manual(values = regioncol) +
          theme_bw()+
          theme(text = element_text(size=11),
                axis.ticks.y = element_blank(),
                axis.ticks.x = element_blank(),
                axis.title.y = element_blank(),
                axis.text.y = element_blank(),
                legend.position = 'none',
                legend.title = element_blank()) +
          scale_y_continuous(position = "right")+
          xlab('Cancer') + ylab('Overlapping Regions\n(Count)') 
        
        figure <- ggarrange(plot1, plot2,
                            labels = c(""),
                            ncol = 2, nrow = 1,
                            heights = c(0.5, 0.5),
                            widths = c(1.5,0.7),
                            align = 'h')
        
        file = gsub(' ','.',paste0(figdir,c,'.',dmr,'_', gsub(', ','',d),'.permutation.overlap.',region,'.raw.png'))
        png(file, height= 1500, width= 1500,res=300)
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
        targ.results2$Sig = ifelse(targ.results2$sig.p < 0.05, 'Sig',ifelse(targ.results2$sig.p < 0.15,'Border.sig','Non-Sig'))
        
        plot1 = ggplot(targ.results2[targ.results2$group == 'permutation',], aes(x = tissue,y=zscore)) +
          geom_boxplot(outlier.shape = NA, aes(x = tissue,y=zscore),col = 'black') + 
          scale_color_manual(values = regioncol) +
          geom_point(data = targ.results2[targ.results2$group == 'observed',], aes(x = tissue, y = zscore,fill = Sig),shape = 23,col = 'black')+
          scale_fill_manual(values = c('Non-Sig'='grey','Border.sig' = '#2A9D8F',Sig='red')) +
          scale_shape_manual(values=c(23))+
          theme_bw()+
          coord_flip()+
          theme(text = element_text(size=11),
                #axis.ticks.y = element_blank(),
                #axis.ticks.x = element_blank(),
                legend.position = 'none',
                legend.title = element_blank()) + 
          scale_x_discrete(position = "bottom")+
          scale_y_continuous(position = "right")+
          xlab('Cancer') + ylab('Overlapping Regions\n(Z-score Normalized)') 
        
        plot2 = ggplot(targ.results1[targ.results1$group == 'observed',], aes(x = tissue,y=n.overlap,fill = region)) + 
          geom_bar(stat= 'identity',position ='stack')+
          coord_flip()+
          scale_fill_manual(values = regioncol) +
          theme_bw()+
          theme(text = element_text(size=11),
                axis.ticks.y = element_blank(),
                axis.ticks.x = element_blank(),
                axis.title.y = element_blank(),
                axis.text.y = element_blank(),
                legend.position = 'none',
                legend.title = element_blank()) +
          scale_y_continuous(position = "right")+
          xlab('Cancer') + ylab('Overlapping Regions\n(Count)') 
        
        figure <- ggarrange(plot1, plot2,
                            labels = c(""),
                            ncol = 2, nrow = 1,
                            heights = c(0.5, 0.5),
                            widths = c(1.5,0.7),
                            align = 'h')
        
        file = gsub(' ','.',paste0(figdir,c,'.',dmr,'_', gsub(', ','',d),'.permutation.overlap.',region,'.padj.png'))
        png(file, height= 1500, width= 1500,res=300)
        print(figure)
        dev.off()
        
        
      }
      
    }
  }
  
  
  
  
  
  
}

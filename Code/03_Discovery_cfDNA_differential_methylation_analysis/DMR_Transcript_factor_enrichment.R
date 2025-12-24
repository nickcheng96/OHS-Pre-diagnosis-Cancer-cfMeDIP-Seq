####tfbs enrichment#####
library("LOLA")
library(rGREAT)
library(Hmisc)
#figdir='/.mounts/labs/awadallalab/private/ncheng/cfmedip_data/cptp_samples/cancer_analysis/feature_analysis/BRCA/CV.figures/hyper.dmrs.vst/'
#targ.markers = readRDS(paste0(figdir,'freq.50.filt.lc.reg.RDS'))
cpg.window.annotated = readRDS('/.mounts/labs/awadallalab/private/ncheng/references/gene.annotations/genomic.annotations.window.300.RDS')
cpg.window.annotated.regulatory = cpg.window.annotated[cpg.window.annotated$count >= 4,]
cpg.window.annotated.regulatory = cpg.window.annotated.regulatory[cpg.window.annotated.regulatory$"Regulatory Element" != 'Non-Regulatory Element',]

#savedir='/.mounts/labs/awadallalab/private/ncheng/manuscripts/early.cancer.risk/all.sample.dmr.cpg5/'
#figdir=paste0(savedir,'figures/')
#dir.create(figdir,recursive = T)


permutation.calc.cpgbalance = function(windows, background, vars = c('CpG Region','Repeat Element','Regulatory Element','Gene Location')) {
  n.vars =length(vars)
  complete.sankey = NULL
  background.windows = cpg.window.annotated[cpg.window.annotated$window %in% background$window,]
  cpg.window.annotated.sig.regions = cpg.window.annotated[cpg.window.annotated$window %in% windows,]
  n.features = length(unique(cpg.window.annotated.sig.regions$window))
  search.space = unique(background.windows$window)
  search.space.window.distribution = cpg_count.all[cpg_count.all$window %in% windows,]
  search.space.window.distribution = data.frame(table(search.space.window.distribution$count))
  sampled.background.list = mclapply(1:3000, function(x) {
    set.seed(x)
    print(x)
    return.df.tmp.list = lapply(as.character(search.space.window.distribution$Var1), function(y) {
      search.space.cpg = cpg.count.split[[y]]
      n.features = search.space.window.distribution[search.space.window.distribution$Var1 == y,'Freq']
      return.df.tmp = search.space.cpg[sample(1:nrow(search.space.cpg),n.features),]
      return(as.character(return.df.tmp$window))
    } )
    return.df.random = unlist(return.df.tmp.list)
    return.df = background.windows[background.windows$window %in% return.df.random,]
    return.df$seed =x
    return(return.df)
  },mc.cores=2 )
  
  permutation.list = list()
  for (v in vars) {
    print(v)
    sub.vars = unique(background.windows[,v])
    for (s in sub.vars) {
      print(s)
      
      targ.var = cpg.window.annotated.sig.regions[cpg.window.annotated.sig.regions[,v] %in% s,]
      targ.df.region = data.frame(seed=0,
                                  permutation='observed',
                                  n.overlap = length(unique(targ.var$window)))
      permutations.list = mclapply(sampled.background.list, function(x) {
        permutated.window.overlap =  x[x[,v] %in% s,]
        return.df.tmp1 =  data.frame(seed=x$seed[1],
                                     permutation='permutation',
                                     n.overlap = length(unique(permutated.window.overlap$window)))
        return(return.df.tmp1)
        
      },mc.cores = 2 )
      
      permutation.calc = do.call('rbind',permutations.list)
      permutation.calc = rbind(permutation.calc, targ.df.region)
      permutation.calc$zscore = (permutation.calc$n.overlap - mean(permutation.calc$n.overlap))/sd(permutation.calc$n.overlap)
      
      permutation.calc$pvalue = pnorm(q=permutation.calc$zscore, 
                                      mean = mean(permutation.calc$zscore),
                                      sd = sd(permutation.calc$zscore), lower.tail=FALSE)
      permutation.calc$variable = v
      permutation.calc$variable.group = s
      permutation.calc$group.size = length(windows)
      permutation.list[[length(permutation.list)+1]] = permutation.calc
      
    }
    
  }
  
  
  permutation.return= do.call('rbind',permutation.list)
  return(permutation.return)
}


for (sex in c('Male','Female')) {
  res.df = readRDS(paste0(savedir,'discovery.',sex,'.nonage.adjusted.dmrs.inserts.1000.all.300.q20.RDS'))
  
  res.df = res.df[!grepl('chrX|chrY',res.df$window),]
  res.df.hyper.sig.raw = res.df[res.df$log2FoldChange > 0.25 & res.df$pvalue < 0.05,]
  res.df.hypo.sig.raw = res.df[res.df$log2FoldChange < -0.25 & res.df$pvalue < 0.05,]
  res.df.hyper.sig.raw = res.df.hyper.sig.raw[order(res.df.hyper.sig.raw$pvalue),]
  res.df.hypo.sig.raw = res.df.hypo.sig.raw[order(res.df.hypo.sig.raw$pvalue),]
  bed.setup = function(sig.dmrs) {
    windows = sig.dmrs
    chr =gsub(':.*','',windows)
    start =gsub('.*:','',gsub('-.*','',windows))
    end =gsub('.*-','',windows)
    return.df = data.frame(chr = chr, start = as.numeric(start),end = as.numeric(end), window = windows)
    
  }
  background.regions = bed.setup(cpg.window.annotated[cpg.window.annotated$window %in% res.df$window,'window'])#[res.df.hyper.sig.raw$window %in% ]
  hyper.regions = bed.setup(res.df.hyper.sig.raw$window[1:2000])#[res.df.hyper.sig.raw$window %in% ]
  hypo.regions = bed.setup(res.df.hypo.sig.raw$window[1:2000])#[res.df.hyper.sig.raw$window %in% ]
  
  tfbs.dir=paste0(savedir,'tfbs/')
  dir.create(tfbs.dir)
  write.table(hyper.regions,paste0(tfbs.dir,sex,'.hyper2000.bed'),col.names=F,row.names=F,sep='\t',quote=F)
  write.table(hypo.regions,paste0(tfbs.dir,sex,'.hypo2000.bed'),col.names=F,row.names=F,sep='\t',quote=F)
  
}

#rgreat
res$window = paste0(res$seqnames,':',res$start,'-',res$end)

res$gene = factor(as.character(res$gene), levels = unique(res$gene))
res.proximal = res[abs(res$distTSS) < 2000,]

setwd(figdir)
write.table(res.proximal$gene,'../proximal.tss.filt.genes.txt', sep = '\t', quote = F,row.names = F,col.names = F)


#res.proximal
res.proximal = ddply(res.proximal[,c('gene','count','porportion.dmrs')],c('gene'),numcolwise(max))
pdf(paste0(figdir,paste0(names(go.list)[go],'.sig.binomfdr.proximalpromoters.pdf')),height = 6 ,width = 10)
plot2 = ggplot(res.proximal[res.proximal$porportion.dmrs > 0.75,], aes(x = gene, y=porportion.dmrs)) +  #Diagnosis_Time
  geom_bar(stat = 'identity') +
  theme_bw()+
  theme(text = element_text(size=10),
        axis.ticks.y = element_blank(),
        axis.ticks.x = element_blank(),
        legend.position = 'bottom',
        axis.text.x = element_text(angle = 70),
        legend.title = element_blank()) + xlab('Gene') + ylab('Significant DMR calls') 
print(plot2)
dev.off()

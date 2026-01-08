#loading in libraries 
library(scales)
library(ggh4x)
library(ggpubr)
library(ggplot2)
library(parallel)
library(plyr)
library(ggforce)

#permutating calculation function
permutation.calc.cpgbalance = function(windows, vars = c('CpG Region','Repeat Element','Regulatory Element','Gene Location')) {
  n.vars =length(vars)
  complete.sankey = NULL
  background.windows = cpg.window.annotated$window
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

#reading in CpG density across windows
cpg_count.all = readRDS('genomic.annotations.window.300.RDS')
cpg.count.split = split(cpg_count.all,cpg_count.all$count)

#setting save directories
savedir='/direcotry/to/save/'
figdir=paste0(savedir,'figures/')
dir.create(figdir,recursive = T)

#assessing enrichment of genomic elements across genomic elements
combined.permutation = NULL
for (sex in c('Female','Male')){
  #reading in results from differential methylation analysis
  res.df = readRDS(paste0(savedir,'discovery.',sex,'.dmrs.RDS'))
  res.df = res.df[!grepl('chrX|chrY',res.df$window),]
  
  #defining top 2000 hypermethylated and hypomethylated regions
  res.df.hyper.sig.raw = res.df[res.df$log2FoldChange > 0.25 & res.df$pvalue < 0.05,]
  res.df.hypo.sig.raw = res.df[res.df$log2FoldChange < -0.25 & res.df$pvalue < 0.05,]
  res.df.hyper.sig.raw = res.df.hyper.sig.raw[order(res.df.hyper.sig.raw$pvalue),]
  res.df.hypo.sig.raw = res.df.hypo.sig.raw[order(res.df.hypo.sig.raw$pvalue),]
  
  #permutationt esting to assess enrichment for genomic elements taking into account CpG density
  all.hyper.permutation.cpgbalance =  permutation.calc.cpgbalance(windows=res.df.hyper.sig.raw$window[1:2000])
  all.hypo.permutation.cpgbalance =  permutation.calc.cpgbalance(windows=res.df.hypo.sig.raw$window[1:2000])
  
  permutation.list = list('Hypermethylated' = all.hyper.permutation.cpgbalance,
                          'Hypomethylated' = all.hypo.permutation.cpgbalance)
    
  #annotating results
  for (p in names(permutation.list)) {
    targ.results = permutation.list[[p]]
    targ.results$direction = p
    targ.results$Sig = ifelse(targ.results$pvalue < 0.01,'Sig','Non-Sig')
    targ.results = targ.results[targ.results$variable %in% c('Regulatory Element','Repeat Element'),]
    targ.results = targ.results[!targ.results$variable.group %in% c('Non-Regulatory Element'),]
    targ.results$comparison =ifelse(sex == 'Female','Breast Cancer','Prostate Cancer')
    targ.results$comparison = factor(as.character(targ.results$comparison),levels = c('Prostate Cancer','Breast Cancer'))
    
    combined.permutation = rbind(combined.permutation,targ.results)
  }
  
  
}

#plotting permutation analysis results 
#defining colours
cpg.cols = c("CpG Islands" ='#29335C',"CpG Shelves"= "#DB2B39","CpG Shores" = '#F3A712',"Open Sea" = "#A89B83")
repeat.cols = c('DNA Repeat'='#F4A0BD','LINE'='#0E9974','LTR'='#435A9D','Non-Repeat Region'='#E5E7E6','Other'='#FFD166','Satellite'='#E98A15','SINE'='#9FBE6E','STR'='#815D94')
reg.cols = c('Enhancer'='#FF8360','Non-Regulatory Element'='#E5E7E6','Promoter'='#7DCE82','Promoter/Enhancer'='#C42847','Silencer'='#2E86AB')
gene.cols = c('3\' UTR'='#593C8F','5\' UTR'='#DB5461','Exon'='#8EF9F3','Intergenic'='#136F63')
base.colors = c(cpg.cols,repeat.cols,reg.cols,gene.cols)

#creating permutation enrichment plot
plot1 = ggplot(combined.permutation[combined.permutation$permutation == 'permutation',], aes(x = variable.group,y=zscore,col= variable.group)) +
  geom_boxplot(outlier.shape = NA, aes(x = variable.group,y=zscore,color = variable.group)) + 
  scale_color_manual(values = base.colors) +
  geom_point(data = combined.permutation[combined.permutation$permutation == 'observed',],mapping = aes(x = variable.group, y = zscore,fill = Sig),shape = 23, col = 'black')+
  scale_fill_manual(values = c('Non-Sig'='grey',Sig='red')) +
  facet_grid2(variable  ~ comparison+direction,  scales = 'free_y', remove_labels=T, space = 'free_y')+
  scale_shape_manual(values=c(23))+
  theme_bw()+
  theme(text = element_text(size=10,face='bold'),
        axis.ticks.y = element_blank(),
        axis.ticks.x = element_blank(),
        legend.position = 'none',
        legend.title = element_blank(),
        panel.spacing.y = unit(0, "lines"),
        strip.background = element_rect(fill="white")

        
  ) + coord_flip()+
  
  xlab('Region') + ylab('DMR Enrichment Z-Score') 

png(paste0(figdir,'combined.permutation','.png'), height= 1200, width= 1850,res=250)
print(plot1)
dev.off()


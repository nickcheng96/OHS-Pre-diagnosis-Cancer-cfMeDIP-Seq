####permutation test for enrichment of these regions####

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
permutation.calc = function(windows, background.regions, vars = c('CpG Region','Repeat Element','Regulatory Element','Gene Location')) {
  n.vars =length(vars)
  complete.sankey = NULL
  background.windows = cpg.window.annotated[cpg.window.annotated$window %in% background.regions,]
  cpg.window.annotated.sig.regions = cpg.window.annotated[cpg.window.annotated$window %in% windows,]
  n.features = length(unique(cpg.window.annotated.sig.regions$window))
  search.space = unique(background.windows$window)
  
  sampled.background.list = mclapply(1:2500, function(x) {
    set.seed(x)
    print(x)
    return.df.tmp = search.space[sample(1:length(search.space),n.features)]
    return.df = background.windows[background.windows$window %in% return.df.tmp,]
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

background = cpg_count.all[cpg_count.all$window %in% rownames(res.df),]
cpg.count.split = split(cpg_count.all,cpg_count.all$count)
#savedir='/.mounts/labs/awadallalab/private/ncheng/manuscripts/early.cancer.risk/all.sample.dmr.cpg5/'
figdir=paste0(savedir,'figures/')
dir.create(figdir,recursive = T)
sex = 'Male'
for (sex in c('Female','Male')){
  
  res.df = readRDS(paste0(savedir,'discovery.',sex,'.nonage.adjusted.dmrs.inserts.1000.all.300.q20.RDS'))
  res.df = res.df[!grepl('chrX|chrY',res.df$window),]
  res.df.hyper.sig.raw = res.df[res.df$log2FoldChange > 0.25 & res.df$pvalue < 0.05,]
  res.df.hypo.sig.raw = res.df[res.df$log2FoldChange < -0.25 & res.df$pvalue < 0.05,]
  res.df.hyper.sig.raw = res.df.hyper.sig.raw[order(res.df.hyper.sig.raw$pvalue),]
  res.df.hypo.sig.raw = res.df.hypo.sig.raw[order(res.df.hypo.sig.raw$pvalue),]
  
  
  all.hyper.permutation.cpgbalance =  permutation.calc.cpgbalance(windows=res.df.hyper.sig.raw$window[1:2000] , background)
  all.hypo.permutation.cpgbalance =  permutation.calc.cpgbalance(windows=res.df.hypo.sig.raw$window[1:2000] , background)
  
  #all.hyper.permutation.all =  permutation.calc(windows=res.df.hyper.sig.raw$window[1:2000] , background.regions=background$window)
  #all.hypo.permutation.all =  permutation.calc(windows=res.df.hypo.sig.raw$window[1:2000] , background.regions=background$window)
  permutation.list = list('all.hyper.raw' = all.hyper.permutation.cpgbalance,
                          'all.hypo.raw' = all.hypo.permutation.cpgbalance)
  saveRDS(permutation.list, paste0(savedir,sex,'.functional.element.permutation.list.RDS'))
  #'top500.hyper' = top500.hyper.permutation,
  # 'top500.hypo' = top500.hypo.permutation)
  
}

#permutation test performance plotting
library(scales)
library(ggh4x)
library(ggpubr)
sex = 'Male'
for (sex in c('Female','Male')){
  permutation.list = readRDS(paste0(savedir,sex,'.functional.element.permutation.list.RDS'))
  
  for (p in names(permutation.list)) {
    targ.results = permutation.list[[p]]
    targ.results$Sig = ifelse(targ.results$pvalue < 0.01,'Sig','Non-Sig')
    targ.results = targ.results[targ.results$variable %in% c('Regulatory Element','Repeat Element'),]
    targ.results = targ.results[!targ.results$variable.group %in% c('Non-Regulatory Element'),]
    plot1 = ggplot(targ.results[targ.results$permutation == 'permutation',], aes(x = variable.group,y=zscore,col= variable.group)) +
      geom_boxplot(outlier.shape = NA, aes(x = variable.group,y=zscore,color = variable.group)) + 
      scale_color_manual(values = base.colors) +
      
      geom_point(data = targ.results[targ.results$permutation == 'observed',],mapping = aes(x = variable.group, y = zscore,fill = Sig),shape = 23, col = 'black')+
      scale_fill_manual(values = c('Non-Sig'='grey',Sig='red')) +
      facet_grid2(variable ~ .,scales = 'free',remove_labels=T)+
      scale_shape_manual(values=c(23))+
      theme_bw()+
      theme(axis.text = element_text(size=11),
            axis.ticks.y = element_blank(),
            axis.ticks.x = element_blank(),
            legend.position = 'none',
            legend.title = element_blank(),
            strip.background = element_blank(),
            strip.text = element_blank()
            
      ) + coord_flip()+
      scale_y_continuous( labels = label_number(accuracy = 0.1)) +
      xlab('Region') + ylab('Overlapping DMRs\n(z-score normalized)') 
    
    plot2 = ggplot(targ.results[targ.results$permutation == 'observed',], aes(x = variable.group,y=n.overlap,fill = variable.group)) + 
      geom_bar(stat= 'identity')+
      scale_fill_manual(values = base.colors) +
      theme_bw()+
      facet_grid2(variable ~ .,scales = 'free',remove_labels=T)+
      theme(axis.text = element_text(size=10),
            axis.ticks.y = element_blank(),
            #axis.ticks.x = element_blank(),
            axis.title.y = element_blank(),
            axis.text.y = element_blank(),
            strip.background = element_blank(),
            strip.text = element_blank(),
            legend.position = 'none') +
      coord_flip()+
      xlab('Region') + ylab('Overlapping DMRs\n(Count)') 
    
    figure <- ggarrange(plot1,plot2,
                        labels = c(""),
                        ncol = 2, nrow = 1,
                        widths = c(0.8, 0.4),
                        heights = c(1,1),
                        align = 'h')
    
    png(paste0(figdir,sex,'.',p,'.permutation','.png'), height =1400,width= 3000,res=400)
    print(figure)
    dev.off()
  }
}

combined.permutation = NULL
for (sex in c('Male','Female')) {
  permutation.list = readRDS(paste0(savedir,sex,'.functional.element.permutation.list.RDS'))
  
  for (p in names(permutation.list)) {
    
    targ.results = permutation.list[[p]]
    targ.results$direction = p
    targ.results$Sig = ifelse(targ.results$pvalue < 0.01,'Sig','Non-Sig')
    targ.results = targ.results[targ.results$variable %in% c('Regulatory Element','Repeat Element'),]
    targ.results = targ.results[!targ.results$variable.group %in% c('Non-Regulatory Element'),]
    targ.results$comparison =ifelse(sex == 'Female','Incident Breast Cancer','Incident Prostate Cancer')
    combined.permutation = rbind(combined.permutation,targ.results)
  }
  
}

cpregion.col = c("CpG Island" ='#29335C',"CpG Shelf"= "#DB2B39","CpG Shore" = '#F3A712',"Open Sea" = "#A89B83")
repeat.col = c('SINE'='#9FBE6E','LINE'='#0E9974','STR'='#815D94','LTR'='#435A9D','DNA Repeat'='#F4A0BD','Satellite'='#E98A15','MER'='#FFD166','Non-repetitive Region'='#E5E7E6')
regulatory.col = c('Promoter'='#7DCE82','Enhancer'='#FF8360','Silencer'='#2E86AB','Non-Regulatory' ='#E5E7E6')
coding.col = c('Protein Coding Gene'='#593C8F','lncRNA'='#8EF9F3','miRNA'='#DB5461','Other Coding Region' = '#136F63','Non-Coding Region' ='#E5E7E6')

cpg.cols = c("CpG Islands" ='#29335C',"CpG Shelves"= "#DB2B39","CpG Shores" = '#F3A712',"Open Sea" = "#A89B83")
repeat.cols = c('DNA Repeat'='#F4A0BD','LINE'='#0E9974','LTR'='#435A9D','Non-Repeat Region'='#E5E7E6','Other'='#FFD166','Satellite'='#E98A15','SINE'='#9FBE6E','STR'='#815D94')
reg.cols = c('Enhancer'='#FF8360','Non-Regulatory Element'='#E5E7E6','Promoter'='#7DCE82','Promoter/Enhancer'='#C42847','Silencer'='#2E86AB')
gene.cols = c('3\' UTR'='#593C8F','5\' UTR'='#DB5461','Exon'='#8EF9F3','Intergenic'='#136F63')
base.colors = c(cpg.cols,repeat.cols,reg.cols,gene.cols)

targ.results=combined.permutation
targ.results$direction = ifelse(targ.results$direction == 'all.hyper.raw','Top 2000 Hypermethylated','Top 2000 Hypomethylated')
targ.results$direction= gsub('Top 2000 ','',targ.results$direction)
targ.results$comparison = gsub('Incident ','',targ.results$comparison)
targ.results$comparison = factor(as.character(targ.results$comparison),levels = c('Prostate Cancer','Breast Cancer'))
library(ggforce)
plot1 = ggplot(targ.results[targ.results$permutation == 'permutation',], aes(x = variable.group,y=zscore,col= variable.group)) +
  geom_boxplot(outlier.shape = NA, aes(x = variable.group,y=zscore,color = variable.group)) + 
  scale_color_manual(values = base.colors) +
  
  geom_point(data = targ.results[targ.results$permutation == 'observed',],mapping = aes(x = variable.group, y = zscore,fill = Sig),shape = 23, col = 'black')+
  scale_fill_manual(values = c('Non-Sig'='grey',Sig='red')) +
  facet_grid2(variable  ~ comparison+direction,  scales = 'free_y', remove_labels=T, space = 'free_y')+
  
  
  scale_shape_manual(values=c(23))+
  theme_bw()+
  theme(text = element_text(size=10,face='bold'),
        # strip.text.x = element_text(size=10,face='bold'),
        axis.ticks.y = element_blank(),
        axis.ticks.x = element_blank(),
        legend.position = 'none',
        legend.title = element_blank(),
        panel.spacing.y = unit(0, "lines"),
        strip.background = element_rect(fill="white")
        
        # strip.background.y = element_blank(),
        # strip.text.y = element_blank()
        
  ) + coord_flip()+
  
  # scale_y_continuous( labels = label_number(accuracy = 0.1)) +
  xlab('Region') + ylab('DMR Enrichment Z-Score') 
png(paste0(figdir,'combnied.all.permutation','.png'), height= 1000, width= 4500,res=400)
print(plot1)
dev.off()

plot2 = ggplot(targ.results[targ.results$permutation == 'observed',], aes(x = variable.group,y=n.overlap,fill = variable.group)) + 
  geom_bar(stat= 'identity')+
  scale_fill_manual(values = base.colors) +
  theme_bw()+
  facet_grid2(variable ~ .,scales = 'free',remove_labels=T)+
  theme(axis.text = element_text(size=11),
        axis.ticks.y = element_blank(),
        #axis.ticks.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        strip.background = element_blank(),
        strip.text = element_blank(),
        legend.position = 'none') +
  coord_flip()+
  xlab('Region') + ylab('Overlapping DMRs\n(Count)') 

figure <- ggarrange(plot1,plot2,
                    labels = c(""),
                    ncol = 2, nrow = 1,
                    widths = c(0.8, 0.4),
                    heights = c(1,1),
                    align = 'h')

png(paste0(figdir,'combined.permutation','.png'), height= 1200, width= 1850,res=250)
print(plot1)
dev.off()


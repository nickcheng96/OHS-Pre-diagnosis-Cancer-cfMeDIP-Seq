###enrichment for tissue specific markers results####
library(ggplot2)
library(reshape2)
savedir="/.mounts/labs/awadallalab/private/ncheng/cfmedip_data/cptp_samples/fragmentation/aix13.updated1/motif.analysis/tissue.specificity/"
figdir="/directory/to/plot/"
setwd(savedir)

#previously identified cell type specific markers from WGBS data
#upload
onehot.list = list('cell.type'='cell.type.onehot.average.RDS',
                   'tissue.type'='tissue.type.onehot.average.RDS',
                   'cell.tissue.type'='cell.tissue.type.onehot.average.RDS')



for (o in names(onehot.list)) {
  targ.files = list.files(pattern=paste0(o,'.cutoff.*','enrichment.summary.RDS'))
  targ.list = lapply(targ.files, readRDS)
  combined.df = do.call('rbind',targ.list)
  combined.df.sig = combined.df[combined.df$seed == 0,]
  for (cutoff.filter in unique(combined.df.sig$onehot.cutoff)) {
    for (comp in unique(combined.df.sig$group.comparison)) {
      targ.permutations = combined.df.sig[combined.df.sig$onehot.cutoff == cutoff.filter & combined.df.sig$group.comparison == comp,]
      targ.permutations$group.comparison = gsub('\\.','\n',targ.permutations$group.comparison)
      padjust.permutation.nobackground = function(enhancer.permutation, vars = c('region','comparison')) {
        
        tmp.observed = enhancer.permutation[enhancer.permutation$seed == 0,]
        tmp.observed.list = split(tmp.observed, tmp.observed[,vars])
        tmp.observed.list = tmp.observed.list
        tmp.observed.list = lapply(tmp.observed.list, function(x) {
          return.df = x
          return.df$p.adjust = p.adjust(return.df$sig.p,method = 'bonferroni')
          return(return.df)
        })
        return.df = rbind(do.call('rbind',tmp.observed.list))
        return(return.df)
        
      }
      targ.permutations = padjust.permutation.nobackground(targ.permutations,vars = 'shared.min')
      targ.permutations$sig.p = targ.permutations$p.adjust
      targ.permutations = targ.permutations[targ.permutations$shared.min <= 10,]
      targ.permutations$sig.p = ifelse(targ.permutations$sig.p > 0.05,1,targ.permutations$sig.p)
      plot1 = ggplot(targ.permutations, aes(x=shared.min,y = tissue)) +
        geom_tile(aes(fill=-log10(sig.p)))+
        scale_y_discrete(position = "right")+
        theme_bw()+
        theme(text = element_text(size=11,face='bold'),
              axis.ticks.y = element_blank(),
              axis.ticks.x = element_blank(),
              legend.position = 'right',
              legend.title = element_blank(),
              strip.background = element_rect(fill="white")
              
        ) +
        scale_fill_gradientn(colors = c('#FFF6F1', "#B0517D", "#482066"), #F5B0C3
                             limits = c(0,max(-log10(targ.permutations$sig.p))))+
        scale_colour_gradientn(colors = c('#FFF6F1', "#B0517D", "#482066"),
                               limits = c(0,max(-log10(targ.permutations$sig.p))))+
        xlab('')
      
      
      png(paste0(figdir,comp,'.',o,'.',cutoff.filter,'.enrichment.tile.png'),height=3000,width=3200,res=400)
      print(plot1)
      dev.off()  
      
    }
    
  }

  
}


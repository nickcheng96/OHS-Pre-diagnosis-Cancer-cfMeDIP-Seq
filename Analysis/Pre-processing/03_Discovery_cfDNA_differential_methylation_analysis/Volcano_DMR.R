####volcano plot of bin 300 DMRs####
female.dmr.df =readRDS(paste0(savedir,'discovery.Female.nonage.adjusted.dmrs.inserts.1000.all.300.q20.RDS'))
male.dmr.df =readRDS(paste0(savedir,'discovery.Male.nonage.adjusted.dmrs.inserts.1000.all.300.q20.RDS'))
female.dmr.df$Comparison = c('Incident Breast Cancer vs Control')
male.dmr.df$Comparison = c('Incident Prostate Cancer vs Control')

dmr.list=list(Male = male.dmr.df,
              Female = female.dmr.df)

combined.dmr.df = do.call('rbind',dmr.list)

#sex separate
for (sex in names(dmr.list)) {
  res.df= dmr.list[[sex]]
  res.df = data.frame(res.df)
  res.df = res.df[res.df$window %in% cpg_count$window,]
  res.df$p.adjust = p.adjust(res.df$pvalue,method = 'fdr')
  res.df$sig =ifelse(res.df$p.adjust < 0.05 & res.df$log2FoldChange > 0.25,'Significantly Hypermethylated',
                     ifelse(res.df$p.adjust < 0.05 & res.df$log2FoldChange < -0.25,'Significantly Hypermethylated',
                            'Non-Significant'))
  
  res.df$sig =ifelse(res.df$pvalue < 0.05 & res.df$log2FoldChange > 0.25,'Significantly Hypermethylated',
                     ifelse(res.df$pvalue < 0.05 & res.df$log2FoldChange < -0.25,'Significantly Hypomethylated',
                            'Non-Significant'))
  
  res.df$sig.overall = ifelse(res.df$p.adjust < 0.05 & abs(res.df$log2FoldChange) > 0.25,
                              'p < 0.05',
                              ifelse(res.df$pvalue < 0.05 & abs(res.df$log2FoldChange) > 0.25,'p < 0.05','p > 0.05'))
  
  plot1 = ggplot(res.df, aes(x = log2FoldChange,y=-log10(pvalue),col= sig.overall)) +
    geom_point(alpha= 0.7) + 
    scale_color_manual(values = c('p < 0.05' = '#2B4570','FDR < 0.05'= '#982649','p > 0.05'='#C0BCB5')) +
    #scale_fill_manual(values = c('Non-Sig'='grey',Sig='red')) +
    theme_bw()+
    theme(text = element_text(size=10,face='bold'),
          axis.ticks.y = element_blank(),
          axis.ticks.x = element_blank(),
          legend.position = 'none',
          legend.title = element_blank(),
    )+
    geom_hline(yintercept = -log10(0.05), linetype ='dotted')+
    
    xlab('log2(Methylation Fold-Change)') + ylab('-log10(p-value)')
  
  if (max(abs(res.df$log2FoldChange),na.rm=T) < 2) {
    plot1 = plot1+scale_x_continuous(limits= c(-2,2))
  } else {
    plot1 = plot1+scale_x_continuous(limits= c(-4,4))
    
  }
  
    
  png(paste0(figdir,sex,'.dmr.raw.logfvolcano.png'),height =1400,width= 800,res=400)
  print(plot1)
  dev.off()

  
}

#sex together
combined.volcano.df = NULL
for (sex in names(dmr.list)) {
  res.df= dmr.list[[sex]]
  res.df = data.frame(res.df)
  res.df = res.df[res.df$window %in% cpg_count$window,]
  res.df$p.adjust = p.adjust(res.df$pvalue,method = 'fdr')
  res.df = res.df[order(res.df$pvalue),]
  
  res.df.top2000.hyper = res.df[res.df$pvalue < 0.05 & res.df$log2FoldChange > 0.25, ][1:2000,]
  res.df.top2000.hypo = res.df[res.df$pvalue < 0.05 & res.df$log2FoldChange  < -0.25, ][1:2000,]
  

  res.df$sig =ifelse(res.df$pvalue < 0.05 & res.df$log2FoldChange > 0.25,'Significantly Hypermethylated',
                     ifelse(res.df$pvalue < 0.05 & res.df$log2FoldChange < -0.25,'Significantly Hypomethylated',
                            'Non-Significant'))
  
  res.df$sig.overall = ifelse(res.df$p.adjust < 0.05 & abs(res.df$log2FoldChange) > 0.25,
                              'p < 0.05',
                              ifelse(res.df$pvalue < 0.05 & abs(res.df$log2FoldChange) > 0.25,'p < 0.05','p > 0.05'))
  
  #top 2000
  res.df$sig = ifelse(res.df$window %in% res.df.top2000.hyper$window,'Significantly Hypermethylated',
                      ifelse(res.df$window %in% res.df.top2000.hypo$window,'Significantly Hypomethylated','Non-Significant'))
  res.df$sig.overall = ifelse(res.df$window %in% c(res.df.top2000.hyper$window,res.df.top2000.hypo$window),'p < 0.05','p > 0.05')
  

  combined.volcano.df= rbind(combined.volcano.df,res.df)
  
  
  
  
}


plot.df = combined.volcano.df[!grepl('chrX|chrY',combined.volcano.df$window),]#[sample(1:nrow(combined.volcano.df),50000),]
plot.df$Comparison = factor(as.character(gsub('Incident','',plot.df$Comparison)),levels=c(' Prostate Cancer vs Control',' Breast Cancer vs Control'))
plot1 = ggplot(plot.df, aes(x = log2FoldChange,y=-log10(pvalue),col= sig)) +
  geom_point(alpha= 0.7) + 
  scale_color_manual(values = c('Significantly Hypomethylated' = '#2B4570','Significantly Hypermethylated'= '#982649','Non-Significant'='#C0BCB5')) +
  #scale_fill_manual(values = c('Non-Sig'='grey',Sig='red')) +
  theme_bw() +
  facet_grid(. ~ Comparison)+#,scales = 'free',remove_labels=T
  scale_shape_manual(values=c(23))+
  theme_bw()+
  theme(text = element_text(size=10,face='bold'),
        axis.ticks.y = element_blank(),
        axis.ticks.x = element_blank(),
        legend.position = 'none',
        legend.title = element_blank(),
        strip.background = element_rect(fill="white")
  ) +
  scale_x_continuous(limits= c(-2,2))+
  scale_y_continuous(limits= c(0,7))+
  
 # scale_x_continuous(limits= c(-3,3))+
  
  geom_hline(yintercept = -log10(0.05), linetype ='dotted')+
  
  
  xlab('log2(cfDNA Methylation Fold-Change)') + ylab('-log10(p-value)')



png(paste0(figdir,'combined.dmr.raw.logfvolcano.png'),height =1200,width= 1600,res=400)
print(plot1)
dev.off()


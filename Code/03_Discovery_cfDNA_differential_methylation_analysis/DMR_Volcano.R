#loading in DMR table output from Discovery_set_cfDNA_DMR.R
female.dmr.df =readRDS(paste0('discovery.Female.dmrs.RDS'))
male.dmr.df =readRDS(paste0('discovery.Male.dmrs.RDS'))

female.dmr.df$Comparison = c('Incident Breast Cancer vs Control')
male.dmr.df$Comparison = c('Incident Prostate Cancer vs Control')

dmr.list=list(Male = male.dmr.df,
              Female = female.dmr.df)

#defining top 2000 regions
combined.volcano.df = NULL
for (sex in names(dmr.list)) {
  #loading dmr table
  res.df= dmr.list[[sex]]
  res.df = data.frame(res.df)
  res.df = res.df[!grepl('chrX|chrY',res.df$window),]
  res.df = res.df[order(res.df$pvalue),]
  
  #defining top 2000 hypomethylated and hypermethylated regions
  res.df.top2000.hyper = res.df[res.df$pvalue < 0.05 & res.df$log2FoldChange > 0.25, ][1:2000,]
  res.df.top2000.hypo = res.df[res.df$pvalue < 0.05 & res.df$log2FoldChange  < -0.25, ][1:2000,]

  res.df$sig = ifelse(res.df$window %in% res.df.top2000.hyper$window,'Top Hypermethylated',
                      ifelse(res.df$window %in% res.df.top2000.hypo$window,'Top Hypomethylated','Non-Significant'))

  
  combined.volcano.df= rbind(combined.volcano.df,res.df)

}

#plotting volcano plot 
plot.df$Comparison = factor(as.character(plot.df$Comparison),levels=c('Prostate Cancer vs Control',' Breast Cancer vs Control'))
plot1 = ggplot(plot.df, aes(x = log2FoldChange,y=-log10(pvalue),col= sig)) +
  geom_point(alpha= 0.7) + 
  scale_color_manual(values = c('Significantly Hypomethylated' = '#2B4570','Significantly Hypermethylated'= '#982649','Non-Significant'='#C0BCB5')) +
  theme_bw() +
  facet_grid(. ~ Comparison)+
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
  
  geom_hline(yintercept = -log10(0.05), linetype ='dotted')+
  
  
  xlab('log2(cfDNA Methylation Fold-Change)') + ylab('-log10(p-value)')



png('combined.dmr.raw.logfvolcano.png',height =1200,width= 1600,res=400)
print(plot1)
dev.off()

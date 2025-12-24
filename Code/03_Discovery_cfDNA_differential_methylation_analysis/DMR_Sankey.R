cpg.window.annotated = readRDS('/.mounts/labs/awadallalab/private/ncheng/references/gene.annotations/genomic.annotations.window.300.RDS')
cpg.window.annotated.regulatory = cpg.window.annotated[cpg.window.annotated$count >= 4,]
cpg.window.annotated.regulatory = cpg.window.annotated.regulatory[cpg.window.annotated.regulatory$"Regulatory Element" != 'Non-Regulatory Element',]


sankey.format = function(windows, vars = c('CpG Region','Repeat Element','Regulatory Element','Gene Location')) {
  n.vars =length(vars)
  complete.sankey = NULL
  cpg.window.annotated.hypo.sig = cpg.window.annotated[cpg.window.annotated$window %in% windows,]
  
  nodes = unique(cpg.window.annotated.hypo.sig[,vars[1]])
  for (i in 1:n.vars){
    if (i != n.vars) {
      lv.tmp = cpg.window.annotated.hypo.sig[,c(vars[i],vars[i+1])]
      colnames(lv.tmp) = c('node','next_node')
      lv.tmp$x = vars[i]
      lv.tmp$next_x = vars[i+1]
      
    } else {
      lv.tmp = data.frame(node = cpg.window.annotated.hypo.sig[,c(vars[i])], next_node = NA)
      colnames(lv.tmp) = c('node','next_node')
      lv.tmp$x = vars[i]
      lv.tmp$next_x = NA
      
    }
    
    
    
    complete.sankey = rbind(complete.sankey,lv.tmp)
    nodes = c(nodes, unique(cpg.window.annotated.hypo.sig[,vars[i]]))
  }
  
  complete.sankey$x = factor(complete.sankey$x,levels= vars)
  complete.sankey$next_x = factor(complete.sankey$next_x,levels =vars)
  
  nodes = data.frame(name = nodes,check.names=F,stringsAsFactors = F)
  return.list = list(complete.sankey,nodes)
  
  return(complete.sankey)
}
sankey.format = function(windows, vars = c('CpG Region','Regulatory Element','Repeat Element','Gene Location')) {
  n.vars =length(vars)
  complete.sankey = NULL
  cpg.window.annotated.hypo.sig = cpg.window.annotated[cpg.window.annotated$window %in% windows,]
  
  nodes = unique(cpg.window.annotated.hypo.sig[,vars[1]])
  for (i in 1:n.vars){
    if (i != n.vars) {
      lv.tmp = cpg.window.annotated.hypo.sig[,c(vars[i],vars[i+1])]
      colnames(lv.tmp) = c('node.tmp','next_node.tmp')
      lv.tmp$x = vars[i]
      lv.tmp$next_x = vars[i+1]
      
    } else {
      lv.tmp = data.frame(node = cpg.window.annotated.hypo.sig[,c(vars[i])], next_node = NA)
      colnames(lv.tmp) = c('node.tmp','next_node.tmp')
      lv.tmp$x = vars[i]
      lv.tmp$next_x = NA
      
    }
    
    
    
    complete.sankey = rbind(complete.sankey,lv.tmp)
    nodes = c(nodes, unique(cpg.window.annotated.hypo.sig[,vars[i]]))
  }
  
  #node percentages
  complete.sankey.list= split(complete.sankey, complete.sankey$x)
  complete.sankey.list.percentage = lapply(complete.sankey.list, function(x) {
    node.percentage = data.frame(round(100*(table(x$node.tmp)/sum(table(x$node.tmp))),digits=1),stringsAsFactors = F)
    colnames(node.percentage) = c('node.tmp','node.pct')
    
    next_node.percentage = data.frame(round(100*(table(x$next_node.tmp)/sum(table(x$next_node.tmp))),digits=1),stringsAsFactors = F)
    if(nrow(next_node.percentage) > 0 ){
      colnames(next_node.percentage) = c('next_node.tmp','next_node.pct')
      x.tmp = merge(x, node.percentage,by= 'node.tmp')
      x.tmp = merge(x.tmp, next_node.percentage,by= 'next_node.tmp')
      x.tmp$node = paste0(x.tmp$node.tmp,' (',x.tmp$node.pct,'%)')
      x.tmp$next_node = paste0(x.tmp$next_node.tmp,' (',x.tmp$next_node.pct,'%)')
    } else {
      x.tmp = merge(x, node.percentage,by= 'node.tmp')
      x.tmp$node = paste0(x.tmp$node.tmp,' (',x.tmp$node.pct,'%)')
      x.tmp$next_node = x.tmp$next_node.tmp
      
    }
    
    
    return(x.tmp[,c('node','next_node','x','next_x')])
    
  } )
  complete.sankey = do.call('rbind',complete.sankey.list.percentage)
  complete.sankey$x = factor(complete.sankey$x,levels= vars)
  complete.sankey$next_x = factor(complete.sankey$next_x,levels =vars)
  
  nodes = data.frame(name = nodes,check.names=F,stringsAsFactors = F)
  return.list = list(complete.sankey,nodes)
  
  return(complete.sankey)
}
#savedir='/.mounts/labs/awadallalab/private/ncheng/manuscripts/early.cancer.risk/all.sample.dmr.cpg5/'
#figdir=paste0(savedir,'figures/')
dir.create(figdir,recursive = T)
sex = 'Female'
for (sex in c('Male','Female')) {
  res.df = readRDS(paste0(savedir,'discovery.',sex,'.nonage.adjusted.dmrs.inserts.1000.all.300.q20.RDS'))
  
  res.df = res.df[!grepl('chrX|chrY',res.df$window),]
  res.df.hyper.sig.raw = res.df[res.df$log2FoldChange > 0.25 & res.df$pvalue < 0.05,]
  res.df.hypo.sig.raw = res.df[res.df$log2FoldChange < -0.25 & res.df$pvalue < 0.05,]
  res.df.hyper.sig.raw = res.df.hyper.sig.raw[order(res.df.hyper.sig.raw$pvalue),]
  res.df.hypo.sig.raw = res.df.hypo.sig.raw[order(res.df.hypo.sig.raw$pvalue),]
  
  sig.hypo.raw.sankey.df = sankey.format(res.df.hypo.sig.raw$window[1:2000] ,  vars = c('CpG Region','Regulatory Element','Repeat Element'))
  sig.hyper.raw.sankey.df = sankey.format(res.df.hyper.sig.raw$window[1:2000] ,  vars = c('CpG Region','Regulatory Element','Repeat Element'))
  
  saveRDS(sig.hypo.raw.sankey.df,paste0(savedir,sex,'.hypo.sankey.RDS'))
  
  saveRDS(sig.hyper.raw.sankey.df,paste0(savedir,sex,'.hyper.sankey.RDS'))
  
}
# Make the Network
library(ggplot2)
library(ggsankey)
setwd(savedir)

cpregion.col = c("CpG Island" ='#29335C',"CpG Shelf"= "#DB2B39","CpG Shore" = '#F3A712',"Open Sea" = "#A89B83")
repeat.col = c('SINE'='#9FBE6E','LINE'='#0E9974','STR'='#815D94','LTR'='#435A9D','DNA Repeat'='#F4A0BD','Satellite'='#E98A15','MER'='#FFD166','Non-repetitive Region'='#E5E7E6')
regulatory.col = c('Promoter'='#7DCE82','Enhancer'='#FF8360','Silencer'='#2E86AB','Non-Regulatory' ='#E5E7E6')
coding.col = c('Protein Coding Gene'='#593C8F','lncRNA'='#8EF9F3','miRNA'='#DB5461','Other Coding Region' = '#136F63','Non-Coding Region' ='#E5E7E6')

cpg.cols = c("CpG Islands" ='#29335C',"CpG Shelves"= "#DB2B39","CpG Shores" = '#F3A712',"Open Sea" = "#A89B83")
repeat.cols = c('DNA Repeat'='#F4A0BD','LINE'='#0E9974','LTR'='#435A9D','Non-Repeat Region'='#E5E7E6','Other'='#FFD166','Satellite'='#E98A15','SINE'='#9FBE6E','STR'='#815D94')
reg.cols = c('Enhancer'='#FF8360','Non-Regulatory Element'='#E5E7E6','Promoter'='#7DCE82','Promoter/Enhancer'='#C42847','Silencer'='#2E86AB')
gene.cols = c('3\' UTR'='#593C8F','5\' UTR'='#DB5461','Exon'='#8EF9F3','Intergenic'='#136F63')
base.colors = c(cpg.cols,repeat.cols,reg.cols,gene.cols)


####sankey plot local####
library(ggplot2)
library(ggsankeyfier)
library(plyr)
library(stringr)
library(ggh4x)
setwd('/Users/ncheng/Desktop/ncheng/Desktop/Lab Work/Manuscripts/OHS_cancer_risk/figures/sankey/')
setwd('/Users/ncheng/Desktop/ncheng/Desktop/Lab Work/Manuscripts/OHS_cancer_risk/figures/sankey/final/')
combined.snankey.all = NULL
for (sex in c('Male','Female')) {
  sig.hyper.raw.sankey.df = readRDS(paste0(sex,'.hyper.sankey.RDS'))
  sig.hypo.raw.sankey.df = readRDS(paste0(sex,'.hypo.sankey.RDS'))
  
  
  sankey.ggsankfier.reformat = function(sig.hyper.raw.sankey.df) {
    return.df=sig.hyper.raw.sankey.df
    return.df.list = split(return.df,return.df$x)
    return.df.list = lapply(return.df.list, function(x) {
      return.tmp = x
      return.tmp$weight = 1/nrow(return.tmp)
      return(return.tmp)
    } )
    return.df = do.call('rbind',return.df.list)
    return.df = ddply(return.df,colnames(return.df)[1:4],numcolwise(sum))
    j = 0
    sankey.df = NULL
    x=return.df
    for (i in 1:nrow(x)) {
      sankey.tmp = data.frame(proportion=x[i,'weight'],
                              connector = c('from','to'),
                              node = c(x[i,1],x[i,2]),
                              stage = c(x[i,3],x[i,4]),
                              edge_id=c(i,i))
      sankey.df = rbind(sankey.df,sankey.tmp)
    }
    sankey.df=sankey.df[!is.na(sankey.df$stage),]
    return(sankey.df)
  }
  hyper.sankey = sankey.ggsankfier.reformat(sig.hyper.raw.sankey.df)
  hypo.sankey = sankey.ggsankfier.reformat(sig.hypo.raw.sankey.df)
  
  
  
  cpregion.col = c("CpG Island" ='#29335C',"CpG Shelf"= "#DB2B39","CpG Shore" = '#F3A712',"Open Sea" = "#A89B83")
  repeat.col = c('SINE'='#9FBE6E','LINE'='#0E9974','STR'='#815D94','LTR'='#435A9D','DNA Repeat'='#F4A0BD','Satellite'='#E98A15','MER'='#FFD166','Non-repetitive Region'='#E5E7E6')
  regulatory.col = c('Promoter'='#7DCE82','Enhancer'='#FF8360','Silencer'='#2E86AB','Non-Regulatory' ='#E5E7E6')
  coding.col = c('Protein Coding Gene'='#593C8F','lncRNA'='#8EF9F3','miRNA'='#DB5461','Other Coding Region' = '#136F63','Non-Coding Region' ='#E5E7E6')
  
  cpg.cols = c("CpG Islands" ='#29335C',"CpG Shelves"= "#DB2B39","CpG Shores" = '#F3A712',"Open Sea" = "#A89B83")
  repeat.cols = c('DNA Repeat'='#F4A0BD','LINE'='#0E9974','LTR'='#435A9D','Non-Repeat Region'='#B0C0BC','Other'='#FFD166','Satellite'='#E98A15','SINE'='#9FBE6E','STR'='#815D94')
  reg.cols = c('Enhancer'='#FF8360','Non-Regulatory Element'='#B0C0BC','Promoter'='#7DCE82','Promoter/Enhancer'='#C42847','Silencer'='#2E86AB')
  gene.cols = c('3\' UTR'='#593C8F','5\' UTR'='#DB5461','Exon'='#8EF9F3','Intergenic'='#136F63')
  base.colors = c(cpg.cols,repeat.cols,reg.cols,gene.cols)
  
  
  #renaming color vectors
  #cpg values
  recolor.labels = function(sig.hypo.raw.sankey.df) {
    unique.vars = unique(as.character(sig.hypo.raw.sankey.df$node),as.character(sig.hypo.raw.sankey.df$next_node))
    unique.vars = unique(gsub(' \\(.*','',gsub('\n.*','',unique.vars)))
    cpg.cols = c("CpG Islands" ='#29335C',"CpG Shelves"= "#DB2B39","CpG Shores" = '#F3A712',"Open Sea" = "#A89B83")
    cpg.cols = cpg.cols[names(cpg.cols) %in% unique.vars]
    cpg.cols = cpg.cols[order(names(cpg.cols))]
    
    cpg.names = unique(sig.hypo.raw.sankey.df[sig.hypo.raw.sankey.df$x == 'CpG Region','node'])
    cpg.names = cpg.names[order(cpg.names)] 
    names(cpg.cols) = cpg.names
    #repeats
    repeat.cols = c('DNA Repeat'='#F4A0BD','LINE'='#0E9974','LTR'='#435A9D','Non-Repeat Region'='#B0C0BC','Other'='#FFD166','Satellite'='#E98A15','SINE'='#9FBE6E','STR'='#815D94')
    repeat.cols = repeat.cols[names(repeat.cols) %in% unique.vars]
    repeat.cols = repeat.cols[order(names(repeat.cols))]
    repeat.names =  unique(sig.hypo.raw.sankey.df[sig.hypo.raw.sankey.df$x == 'Repeat Element','node'])
    repeat.names = repeat.names[order(repeat.names)]
    names(repeat.cols) = repeat.names
    #reg elements
    reg.cols = c('Enhancer'='#FF8360','Non-Regulatory Element'='#B0C0BC','Promoter'='#7DCE82','Promoter/Enhancer'='#C42847','Silencer'='#2E86AB')
    reg.cols = reg.cols[names(reg.cols) %in% unique.vars]
    reg.cols = reg.cols[order(names(reg.cols))]
    
    reg.names =  unique(sig.hypo.raw.sankey.df[sig.hypo.raw.sankey.df$x == 'Regulatory Element','node'])
    reg.names = reg.names[order(reg.names)]
    names(reg.cols) = reg.names
    #Gene Location
    gene.cols = c('3\' UTR'='#593C8F','5\' UTR'='#DB5461','Exon'='#8EF9F3','Intergenic'='#136F63')
    gene.cols = gene.cols[names(gene.cols) %in% unique.vars]
    gene.cols = gene.cols[order(names(gene.cols))]
    
    gene.names =  unique(sig.hypo.raw.sankey.df[sig.hypo.raw.sankey.df$x == 'Gene Location','node'])
    gene.names = gene.names[order(gene.names)]
    names(gene.cols) = gene.names
    
    targ.colors = c(cpg.cols,repeat.cols,reg.cols,gene.cols)
    return(targ.colors)
  }
  
  targ.colors.hyper = recolor.labels(sig.hyper.raw.sankey.df)
  pos <- position_sankey(v_space = 0, align = "justify")
  pos_text <- position_sankey(v_space = 0.02, align = "bottom", nudge_x = 0.1,direction='forward')
  #pos_text <- position_sankey(v_space = 'auto', align = "bottom", nudge_x = 0.1,direction='forward')
  
  hyper.sankey$node = gsub(' \\(','\n\\(',hyper.sankey$node)
  names(targ.colors.hyper) = gsub(' \\(','\n\\(',names(targ.colors.hyper))
  hyper.sankey$direction = 'Top 2000 Hypermethylated'
  
  targ.colors.hypo = recolor.labels(sig.hypo.raw.sankey.df)
  pos <- position_sankey(v_space = 0, align = "justify")
  pos_text <- position_sankey(v_space = 0.02, align = "bottom", nudge_x = 0.1,direction='forward')
  #pos_text <- position_sankey(v_space = 'auto', align = "bottom", nudge_x = 0.1,direction='forward')
  
  hypo.sankey$node = gsub(' \\(','\n\\(',hypo.sankey$node)
  names(targ.colors.hypo) = gsub(' \\(','\n\\(',names(targ.colors.hypo))
  hypo.sankey$direction = 'Top 2000 Hypomethylated'
  hyper.sankey$edge_id1=hyper.sankey$edge_id
  hypo.sankey$edge_id1 = hypo.sankey$edge_id + max(hyper.sankey$edge_id)
  
  combined.snankey.tmp = rbind(hyper.sankey,hypo.sankey)
  combined.snankey.tmp$node = gsub('\n.*','',combined.snankey.tmp$node)
  combined.snankey.tmp$sex = ifelse(sex == 'Male','Incident Prostate Cancer','Incident Breast Cancer')
  combined.snankey.all = rbind(combined.snankey.all,combined.snankey.tmp)
  
  
  plot1= ggplot(hyper.sankey,
                aes(x = stage, y = proportion, group = node,
                    connector = connector, edge_id = edge_id, fill =node)) +
    geom_sankeyedge(v_space = 0) +
    geom_sankeynode(v_space = 0) +
    theme_bw() +
    
    theme(text = element_text(size=10,face='bold'),
          axis.ticks.y = element_blank(),
          axis.ticks.x = element_blank(),
          axis.title.x=element_blank(),
          legend.position = 'none',
          legend.title = element_blank(),
    )+
    scale_x_discrete(expand=c(0,0,0,0.7))+
    scale_y_continuous(limits=c(0,1))+
    scale_fill_manual(values=targ.colors.hyper)+
    
    scale_color_manual(values=targ.colors.hyper)+
    ylab('Proportion of Regions')+
    facet_grid2(direction ~ .,scales = 'free')
  print(plot1)
  png(paste0(sex,'.hyper.blank.sankey.png'),height=1000,width=1500,res=400)
  print(plot1)
  dev.off()
  
  png(paste0(sex,'.hyper.labeled.sankey.png'),height=1000,width=1500,res=300)
  
  plot1=plot1+geom_label(aes(label = str_wrap(node, 20),col=node),
                         color='white', position = pos_text, stat = "sankeynode",
                         hjust = 0,cex=2,alpha=0.7)
  print(plot1)
  dev.off()
  
  targ.colors.hypo = recolor.labels(sig.hypo.raw.sankey.df)
  pos <- position_sankey(v_space = 0, align = "justify")
  pos_text <- position_sankey(v_space = 0.02, align = "bottom", nudge_x = 0.1,direction='forward')
  #pos_text <- position_sankey(v_space = 'auto', align = "bottom", nudge_x = 0.1,direction='forward')
  
  hypo.sankey$node = gsub(' \\(','\n\\(',hypo.sankey$node)
  names(targ.colors.hypo) = gsub(' \\(','\n\\(',names(targ.colors.hypo))
  
  plot1= ggplot(hypo.sankey,
                aes(x = stage, y = proportion, group = node,
                    connector = connector, edge_id = edge_id, fill =node)) +
    geom_sankeyedge(v_space = 0) +
    geom_sankeynode(v_space = 0) +
    theme_bw() +
    
    theme(text = element_text(size=10,face='bold'),
          axis.ticks.y = element_blank(),
          axis.ticks.x = element_blank(),
          axis.title.x=element_blank(),
          legend.position = 'none',
          legend.title = element_blank(),
    )+
    scale_x_discrete(expand=c(0,0,0,0.7))+
    scale_fill_manual(values=targ.colors.hypo)+
    
    scale_color_manual(values=targ.colors.hypo)+
    ylab('Proportion of Regions')
  print(plot1)
  png(paste0(sex,'.hypo.blank.sankey.png'),height=1000,width=1500,res=400)
  print(plot1)
  dev.off()
  
  png(paste0(sex,'.hypo.labeled.sankey.png'),height=1000,width=1500,res=300)
  
  print(plot1+geom_label(aes(label = str_wrap(node, 20),col=node),
                         color='white', position = pos_text, stat = "sankeynode",
                         hjust = 0,cex=2,alpha=0.7))
  dev.off()
  
  
}

levels = c('Open Sea','CpG Islands','CpG Shores','CpG Shelves',
           'Non-Regulatory Element','Silencer','Enhancer','Promoter/Enhancer','Promoter',
           'Non-Repeat Region','SINE','LINE','LTR','DNA Repeat','STR','Satellite','Other')
combined.snankey.all$sex = factor(as.character(combined.snankey.all$sex), levels= c("Incident Prostate Cancer","Incident Breast Cancer"))
combined.snankey.all$node = factor(as.character(combined.snankey.all$node),levels =rev(levels))
plot1= ggplot(combined.snankey.all,
              aes(x = stage, y = proportion, node=node, group = node,
                  connector = connector, edge_id = edge_id1, fill =node)) +
  geom_sankeyedge(position= position_sankey(order = "as_is")) +  #v_space = 0,+
  geom_sankeynode(position= position_sankey(order = "as_is")) +
  theme_bw() +
  facet_grid(direction ~ sex)+
  theme(text = element_text(size=10,face='bold'),
        axis.ticks.y = element_blank(),
        axis.ticks.x = element_blank(),
        #axis.title.x=element_blank(),
        legend.position = 'none',
        legend.title = element_blank(),
        strip.background = element_rect(fill="white")
  )+
  scale_x_discrete(expand=c(0,0,0,0.7))+
  scale_y_continuous(limits=c(0,1.01))+
  scale_fill_manual(values=base.colors)+#base.colors
  xlab('Genomic Annotation')#+ facet_grid2(direction ~ .)
print(plot1)

png(paste0('combined.blank.sankey.png'),height=1800,width=3200,res=500)
print(plot1)
dev.off()
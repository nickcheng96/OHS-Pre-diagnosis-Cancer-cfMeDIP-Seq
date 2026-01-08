#creating BED files for top 2000 hypermethylated and hypomethylated regions
for (sex in c('Male','Female')) {
  res.df = readRDS(paste0(savedir,'discovery.',sex,'.dmrs.RDS'))
  
  #selecting top 2000 hypermethyled and hypomethylated regions
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
  hyper.regions = bed.setup(res.df.hyper.sig.raw$window[1:2000])
  hypo.regions = bed.setup(res.df.hypo.sig.raw$window[1:2000])
  

  write.table(hyper.regions,paste0(sex,'.hyper2000.bed'),col.names=F,row.names=F,sep='\t',quote=F)
  write.table(hypo.regions,paste0(sex,'.hypo2000.bed'),col.names=F,row.names=F,sep='\t',quote=F)
  
}


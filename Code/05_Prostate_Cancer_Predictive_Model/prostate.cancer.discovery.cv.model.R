
options(expressions = 5e5)
library(DESeq2)
library("BiocParallel")
ncores = 10
register(MulticoreParam(ncores))

#wi
#####functions#####
auc_calc = function(prediction_table,labels = c('Control','Cancer')) {
  tmp = prediction_table
  tmp = tmp[order(-tmp$methylation_score),]
  #labels = c(tmp$reported)
  pred = prediction(predictions = c(tmp$methylation_score) ,labels =  tmp$reported, labels)
  perf_AUC=performance(pred,"auc") #Calculate the AUC value
  AUC=perf_AUC@y.values[[1]]
  return(AUC)
}

pca.plot.fillerlabel.project =function(deseq.matrix,predx.dmrs.sig.hyper,combined.info,dx.time,name,cv.label = T){
  library(ggfortify)
  library(gridExtra)
  library(ggpubr)
  library(ggplot2)
  library(RColorBrewer)
  
  diagnosis_time_grouping = function(diagnosis_time) {
    tmp = ifelse(diagnosis_time > 2920, '8-10', diagnosis_time)
    tmp = ifelse(diagnosis_time <= 2920, '6-8', tmp)
    tmp = ifelse(diagnosis_time <= 2190, '4-6', tmp)
    tmp = ifelse(diagnosis_time <= 1460, '2-4', tmp)
    tmp = ifelse(diagnosis_time <= 730, '0-2', tmp)
    tmp[is.na(diagnosis_time)] = 'Control'
    diagnosis_time_groups = c('Control','0-2','2-4','4-6','6-8','8-10')
    tmp = factor(tmp, levels = rev(diagnosis_time_groups))
    return(tmp)
  }
  
  if (cv.label == T){
    tmp.info = combined.info[combined.info$GRP_Id  %in% colnames(deseq.matrix),]
    tmp.info$Diagnosis_Time = diagnosis_time_grouping(tmp.info[,dx.time])
    predx.dmrs.sig.hyper.filt = predx.dmrs.sig.hyper[predx.dmrs.sig.hyper %in% rownames(deseq.matrix)]
    pca.obj = prcomp(data.frame(t(deseq.matrix[predx.dmrs.sig.hyper.filt,tmp.info$GRP_Id]),check.names=F))
    
    # project new data onto the PCA space
    
    return.df = data.frame(pca.obj$x,check.names=F)
    pca.plot = return.df
    pca.plot$GRP_Id = rownames(pca.plot)
    pca.plot = merge(pca.plot, tmp.info,by='GRP_Id')
    
    pc1 = paste0('PC1 (',100*round(summary(pca.obj)$importance[2,1],digits = 3),'%)')
    pc2 = paste0('PC2 (',100*round(summary(pca.obj)$importance[2,2],digits = 3),'%)')
    pc3 = paste0('PC3 (',100*round(summary(pca.obj)$importance[2,3],digits = 3),'%)')
    diagnosis_time_colors1 = c('#7A797C',"#048BA8",'#AAF683','#FFD97D','#FF9B85','#C8553D')
    names(diagnosis_time_colors1) = c('Control','0-2','2-4','4-6','6-8','8-10')
    plot1 = ggplot(pca.plot, aes(x = PC1, y = PC2, col = Diagnosis_Time,shape = filler))+ geom_point(size = 3,alpha = 0.8) + 
      scale_color_manual(values = diagnosis_time_colors1)+#+ ggtitle(title) +
      theme_bw()+
      scale_shape_manual(values = c("MFiller" = 19,"UFiller" = 17))+
      theme(text = element_text(size=12),
            axis.text=element_text(size=12, face = "bold"),
            axis.title=element_text(size=14,face="bold"),
            legend.position = "none")+ guides(colour = guide_legend(override.aes = list(size=4),nrow=7)) + xlab(pc1) + ylab(pc2)+ ggtitle(paste0( ' Features: ',length(predx.dmrs.sig.hyper)))
    
    plot2 = ggplot(pca.plot, aes(x = PC1, y = PC3, col = Diagnosis_Time,shape = filler))+ geom_point(size = 3,alpha = 0.8) + 
      scale_color_manual(values = diagnosis_time_colors1)+#+ ggtitle(title) +
      theme_bw()+
      scale_shape_manual(values = c("MFiller" = 19,"UFiller" = 17))+
      theme(text = element_text(size=12),
            axis.text=element_text(size=12, face = "bold"),
            axis.title=element_text(size=14,face="bold"),
            legend.position = "none")+ guides(colour = guide_legend(override.aes = list(size=4),nrow=7)) + xlab(pc1) + ylab(pc3)+ ggtitle(paste0( ' Features: ',length(predx.dmrs.sig.hyper)))
    
    
    plot3 = ggplot(pca.plot, aes(x = PC2, y = PC3, col = Diagnosis_Time,shape = filler))+ geom_point(size = 3,alpha = 0.8) + 
      scale_color_manual(values = diagnosis_time_colors1)+#+ ggtitle(title) +
      theme_bw()+
      scale_shape_manual(values = c("MFiller" = 19,"UFiller" = 17))+
      theme(text = element_text(size=12),
            axis.text=element_text(size=12, face = "bold"),
            axis.title=element_text(size=14,face="bold"),
            legend.position = "none")+ guides(colour = guide_legend(override.aes = list(size=4),nrow=7)) + xlab(pc2) + ylab(pc3)+ ggtitle(paste0( ' Features: ',length(predx.dmrs.sig.hyper)))
    
    title=paste0(name,'.dx.time.pca.png')
    # title=paste0(figdir,title)
    png(title,height = 1000, width = 3000,res=200)
    figure <- ggarrange(plot1, plot2,plot3,
                        labels = c(""),
                        ncol = 3, nrow = 1
                        # heights = c(0.3, 0.5),
                        # widths = c(1.5,0.7)
    )
    print(figure)
    dev.off()
    plot4 = ggplot(pca.plot, aes(x = PC1, y = PC2, col = Diagnosis_Time,shape = filler))+ geom_point(size = 3,alpha = 0.8) + 
      scale_color_manual(values = diagnosis_time_colors1)+#+ ggtitle(title) +
      theme_bw()+
      scale_shape_manual(values = c("Methylated" = 19,"Unmethylated" = 17))+
      theme(text = element_text(size=12),
            axis.text=element_text(size=12, face = "bold"),
            axis.title=element_text(size=14,face="bold"),
            legend.position = "right")+ guides(colour = guide_legend(override.aes = list(size=4),nrow=7)) + xlab(pc2) + ylab(pc3)+ ggtitle(paste0( ' Features: ',length(predx.dmrs.sig.hyper)))
    png(paste0(name,'.dxtime.pca.labs.png'),height = 1000, width = 3000,res=200)
    print(plot4)
    dev.off()
    
    
    n <- length(unique(tmp.info$seq.date))
    dates = unique((tmp.info$seq.date))
    dates = as.character(dates[order(dates)])
    pca.plot$seq.date=factor(as.character(pca.plot$seq.date),levels = dates)
    qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
    col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
    
    plot1 = ggplot(pca.plot, aes(x = PC1, y = PC2, col = seq.date,shape = filler))+ geom_point(size = 3,alpha = 0.8) + 
      scale_color_manual(values = col_vector)+#+ ggtitle(title) +
      theme_bw()+
      scale_shape_manual(values = c("MFiller" = 19,"UFiller" = 17))+
      theme(text = element_text(size=12),
            axis.text=element_text(size=12, face = "bold"),
            axis.title=element_text(size=14,face="bold"),
            legend.position = "none")+ guides(colour = guide_legend(override.aes = list(size=4),nrow=7)) + xlab(pc1) + ylab(pc3)+ ggtitle(paste0( ' Features: ',length(predx.dmrs.sig.hyper)))
    
    plot2 = ggplot(pca.plot, aes(x = PC2, y = PC3, col = seq.date,shape = filler))+ geom_point(size = 3,alpha = 0.8) + 
      scale_color_manual(values = col_vector)+#+ ggtitle(title) +
      theme_bw()+
      scale_shape_manual(values = c("MFiller" = 19,"UFiller" = 17))+
      theme(text = element_text(size=12),
            axis.text=element_text(size=12, face = "bold"),
            axis.title=element_text(size=14,face="bold"),
            legend.position = "none")+ guides(colour = guide_legend(override.aes = list(size=4),nrow=7)) + xlab(pc1) + ylab(pc3)+ ggtitle(paste0( ' Features: ',length(predx.dmrs.sig.hyper)))
    
    
    plot3 = ggplot(pca.plot, aes(x = PC1, y = PC3, col = seq.date,shape = filler))+ geom_point(size = 3,alpha = 0.8) + 
      scale_color_manual(values = col_vector)+#+ ggtitle(title) +
      theme_bw()+
      scale_shape_manual(values = c("MFiller" = 19,"UFiller" = 17))+
      theme(text = element_text(size=12),
            axis.text=element_text(size=12, face = "bold"),
            axis.title=element_text(size=14,face="bold"),
            legend.position = "none")+ guides(colour = guide_legend(override.aes = list(size=4),nrow=7)) + xlab(pc2) + ylab(pc3)+ ggtitle(paste0( ' Features: ',length(predx.dmrs.sig.hyper)))
    
    title=paste0(name,'.batch.pca.png')
    # title=paste0(figdir,title)
    png(title,height = 1000, width = 3000,res=200)
    figure <- ggarrange(plot1, plot2,plot3,
                        labels = c(""),
                        ncol = 3, nrow = 1
                        # heights = c(0.3, 0.5),
                        # widths = c(1.5,0.7)
    )
    print(figure)
    dev.off()
    
    plot3 = ggplot(pca.plot, aes(x = PC1, y = PC2, col = seq.date,shape = filler))+ geom_point(size = 3,alpha = 0.8) + 
      scale_color_manual(values = col_vector)+#+ ggtitle(title) +
      theme_bw()+
      scale_shape_manual(values = c("MFiller" = 19,"UFiller" = 17))+
      theme(text = element_text(size=12),
            axis.text=element_text(size=12, face = "bold"),
            axis.title=element_text(size=14,face="bold"),
            legend.position = "right")+ guides(colour = guide_legend(override.aes = list(size=4),nrow=7)) + xlab(pc2) + ylab(pc3)+ ggtitle(paste0( ' Features: ',length(predx.dmrs.sig.hyper)))
    png(paste0(name,'.batch.pca.labs.png'),height = 1000, width = 3000,res=200)
    print(plot3)
    dev.off()
    
    
    n <- length(unique(tmp.info$filler))
    qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
    col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
    
    plot1 = ggplot(pca.plot, aes(x = PC1, y = PC2, col = filler,shape = filler))+ geom_point(size = 3,alpha = 0.8) + 
      scale_color_manual(values = col_vector)+#+ ggtitle(title) +
      theme_bw()+
      scale_shape_manual(values = c("MFiller" = 19,"UFiller" = 17))+
      theme(text = element_text(size=12),
            axis.text=element_text(size=12, face = "bold"),
            axis.title=element_text(size=14,face="bold"),
            legend.position = "none")+ guides(colour = guide_legend(override.aes = list(size=4),nrow=7)) + xlab(pc1) + ylab(pc2)+ ggtitle(paste0( ' Features: ',length(predx.dmrs.sig.hyper)))
    
    plot2 = ggplot(pca.plot, aes(x = PC1, y = PC3, col = filler,shape = filler))+ geom_point(size = 3,alpha = 0.8) + 
      scale_color_manual(values = col_vector)+#+ ggtitle(title) +
      theme_bw()+
      scale_shape_manual(values = c("MFiller" = 19,"UFiller" = 17))+
      theme(text = element_text(size=12),
            axis.text=element_text(size=12, face = "bold"),
            axis.title=element_text(size=14,face="bold"),
            legend.position = "none")+ guides(colour = guide_legend(override.aes = list(size=4),nrow=7)) + xlab(pc1) + ylab(pc3)+ ggtitle(paste0( ' Features: ',length(predx.dmrs.sig.hyper)))
    
    
    plot3 = ggplot(pca.plot, aes(x = PC2, y = PC3, col = filler,shape = filler))+ geom_point(size = 3,alpha = 0.8) + 
      scale_color_manual(values = col_vector)+#+ ggtitle(title) +
      theme_bw()+
      scale_shape_manual(values = c("MFiller" = 19,"UFiller" = 17))+
      theme(text = element_text(size=12),
            axis.text=element_text(size=12, face = "bold"),
            axis.title=element_text(size=14,face="bold"),
            legend.position = "none")+ guides(colour = guide_legend(override.aes = list(size=4),nrow=7)) + xlab(pc2) + ylab(pc3)+ ggtitle(paste0( ' Features: ',length(predx.dmrs.sig.hyper)))
    
    title=paste0(name,'.filler.pca.png')
    # title=paste0(figdir,title)
    png(title,height = 1000, width = 3000,res=200)
    figure <- ggarrange(plot1, plot2,plot3,
                        labels = c(""),
                        ncol = 3, nrow = 1
                        # heights = c(0.3, 0.5),
                        # widths = c(1.5,0.7)
    )
    print(figure)
    dev.off()
    
    plot3 = ggplot(pca.plot, aes(x = PC1, y = PC2, col = filler,shape = filler))+ geom_point(size = 3,alpha = 0.8) + 
      scale_color_manual(values = col_vector)+#+ ggtitle(title) +
      theme_bw()+
      scale_shape_manual(values = c("MFiller" = 19,"UFiller" = 17))+
      theme(text = element_text(size=12),
            axis.text=element_text(size=12, face = "bold"),
            axis.title=element_text(size=14,face="bold"),
            legend.position = "right")+ guides(colour = guide_legend(override.aes = list(size=4),nrow=7)) + xlab(pc2) + ylab(pc3)+ ggtitle(paste0( ' Features: ',length(predx.dmrs.sig.hyper)))
    png(paste0(name,'.filler.pca.labs.png'),height = 1000, width = 3000,res=200)
    print(plot3)
    dev.off()
    cancer.colors = c('Control' = '#BFCBC2','Pancreas' ='#09814A','Breast'='#f564a9','Prostate'= '#0094C6')
    
    plot1 = ggplot(pca.plot, aes(x = PC1, y = PC2, col = Cancer,shape = Sex))+ geom_point(size = 3,alpha = 0.8) + 
      scale_color_manual(values = cancer.colors)+#+ ggtitle(title) +
      theme_bw()+
      scale_shape_manual(values = c("Male" = 19,"Female" = 17))+
      theme(text = element_text(size=12),
            axis.text=element_text(size=12, face = "bold"),
            axis.title=element_text(size=14,face="bold"),
            legend.position = "none")+ guides(colour = guide_legend(override.aes = list(size=4),nrow=7)) + xlab(pc1) + ylab(pc2)+ ggtitle(paste0( ' Features: ',length(predx.dmrs.sig.hyper)))
    
    plot2 = ggplot(pca.plot, aes(x = PC1, y = PC3, col = Cancer,shape = Sex))+ geom_point(size = 3,alpha = 0.8) + 
      scale_color_manual(values = cancer.colors)+#+ ggtitle(title) +
      theme_bw()+
      scale_shape_manual(values = c("Male" = 19,"Female" = 17))+
      theme(text = element_text(size=12),
            axis.text=element_text(size=12, face = "bold"),
            axis.title=element_text(size=14,face="bold"),
            legend.position = "none")+ guides(colour = guide_legend(override.aes = list(size=4),nrow=7)) + xlab(pc1) + ylab(pc3)+ ggtitle(paste0( ' Features: ',length(predx.dmrs.sig.hyper)))
    
    
    plot3 = ggplot(pca.plot, aes(x = PC2, y = PC3, col = Cancer,shape = Sex))+ geom_point(size = 3,alpha = 0.8) + 
      scale_color_manual(values = cancer.colors)+#+ ggtitle(title) +
      theme_bw()+
      scale_shape_manual(values = c("Male" = 19,"Female" = 17))+
      theme(text = element_text(size=12),
            axis.text=element_text(size=12, face = "bold"),
            axis.title=element_text(size=14,face="bold"),
            legend.position = "none")+ guides(colour = guide_legend(override.aes = list(size=4),nrow=7)) + xlab(pc2) + ylab(pc3)+ ggtitle(paste0( ' Features: ',length(predx.dmrs.sig.hyper)))
    
    title=paste0(name,'.filler.pca.png')
    # title=paste0(figdir,title)
    png(title,height = 1000, width = 3000,res=200)
    figure <- ggarrange(plot1, plot2,plot3,
                        labels = c(""),
                        ncol = 3, nrow = 1
                        # heights = c(0.3, 0.5),
                        # widths = c(1.5,0.7)
    )
    print(figure)
    dev.off()
    
    plot3 = ggplot(pca.plot, aes(x = PC1, y = PC2, col = Cancer,shape = Sex))+ geom_point(size = 3,alpha = 0.8) + 
      scale_color_manual(values = cancer.colors)+#+ ggtitle(title) +
      theme_bw()+
      scale_shape_manual(values = c("Male" = 19,"Female" = 17))+
      theme(text = element_text(size=12),
            axis.text=element_text(size=12, face = "bold"),
            axis.title=element_text(size=14,face="bold"),
            legend.position = "right")+ guides(colour = guide_legend(override.aes = list(size=4),nrow=7)) + xlab(pc2) + ylab(pc3)+ ggtitle(paste0( ' Features: ',length(predx.dmrs.sig.hyper)))
    png(paste0(name,'.group.pca.labs.png'),height = 1000, width = 3000,res=200)
    print(plot3)
    dev.off()
    
  } 
  
  return(return.df)
}

pca.plot.traintest =function(deseq.matrix,predx.dmrs.sig.hyper.all,combined.info,dx.time,name,groups = F,labs = F,plot.pc=T){
  library(ggfortify)
  library(gridExtra)
  library(ggpubr)
  
  
  diagnosis_time_grouping = function(diagnosis_time) {
    tmp = ifelse(diagnosis_time > 2920, '8-10', diagnosis_time)
    tmp = ifelse(diagnosis_time <= 2920, '6-8', tmp)
    tmp = ifelse(diagnosis_time <= 2190, '4-6', tmp)
    tmp = ifelse(diagnosis_time <= 1460, '2-4', tmp)
    tmp = ifelse(diagnosis_time <= 730, '0-2', tmp)
    tmp[is.na(diagnosis_time)] = 'Control'
    diagnosis_time_groups = c('Control','0-2','2-4','4-6','6-8','8-10')
    tmp = factor(tmp, levels = rev(diagnosis_time_groups))
    return(tmp)
  }
  
  
  
  predx.dmrs.sig.hyper = predx.dmrs.sig.hyper.all[!is.na(predx.dmrs.sig.hyper.all)]
  tmp.info = combined.info[combined.info$GRP_Id  %in% colnames(deseq.matrix),]
  tmp.info$Diagnosis_Time = diagnosis_time_grouping(tmp.info[,dx.time])
  tmp.info$train.test = ifelse(tmp.info$GRP_Id %in%train.set$GRP_Id,'Train','Test')
  tmp.info.train = tmp.info[tmp.info$GRP_Id %in% train.set$GRP_Id,]
  tmp.info.test = tmp.info[tmp.info$GRP_Id %in% test.set$GRP_Id,]
  
  pca.obj = prcomp(t(deseq.matrix[predx.dmrs.sig.hyper,tmp.info.train$GRP_Id]))
  tmp.scale= scale(t(deseq.matrix[predx.dmrs.sig.hyper,tmp.info.test$GRP_Id]), pca.obj$center, pca.obj$scale) %*% pca.obj$rotation 
  pc1 = paste0('PC1 (',100*round(summary(pca.obj)$importance[2,1],digits = 3),'%)')
  pc2 = paste0('PC2 (',100*round(summary(pca.obj)$importance[2,2],digits = 3),'%)')
  pc3 = paste0('PC3 (',100*round(summary(pca.obj)$importance[2,3],digits = 3),'%)')
  diagnosis_time_colors1 = c('#7A797C',"#048BA8",'#AAF683','#FFD97D','#FF9B85','#C8553D')
  names(diagnosis_time_colors1) = c('Control','0-2','2-4','4-6','6-8','8-10')
  pc.matrix = data.frame(pca.obj$x)
  pc.matrix$GRP_Id = rownames(pc.matrix)
  
  pc.matrix1 = data.frame(tmp.scale)
  pc.matrix1$GRP_Id = rownames(pc.matrix1)
  
  if (plot.pc == T){
    pca.matrix = merge(tmp.info,rbind(pc.matrix,pc.matrix1),  by = 'GRP_Id')
    #pc.matrix$Batch = factor(as.character(pc.matrix$Batch), levels 
    
    plot1 = ggplot(pca.matrix,aes(x = PC1, y = PC2, col = Diagnosis_Time,shape = train.test)) + geom_point(alpha = 0.75,size=3) +
      scale_color_manual(values = diagnosis_time_colors1)+#+ ggtitle(title) +
      scale_shape_manual(values = c(Train = 16,Test = 15))+
      theme_bw()+
      theme(text = element_text(size=12),
            axis.text=element_text(size=12, face = "bold"),
            axis.title=element_text(size=14,face="bold"),
            legend.position = "none")+ guides(colour = guide_legend(override.aes = list(size=4),nrow=4)) + xlab(pc1) + ylab(pc2)
    plot2 = ggplot(pca.matrix,aes(x = PC1, y = PC3, col = Diagnosis_Time,shape = train.test)) + geom_point(alpha = 0.75,size=3) +
      scale_color_manual(values = diagnosis_time_colors1)+#+ ggtitle(title) +
      scale_shape_manual(values = c(Train = 16,Test = 15))+
      theme_bw()+
      theme(text = element_text(size=12),
            axis.text=element_text(size=12, face = "bold"),
            axis.title=element_text(size=14,face="bold"),
            legend.position = "none")+ guides(colour = guide_legend(override.aes = list(size=4),nrow=4)) + xlab(pc1) + ylab(pc3)
    
    plot3 = ggplot(pca.matrix,aes(x = PC2, y = PC3, col = Diagnosis_Time,shape = train.test)) + geom_point(alpha = 0.75,size=3) +
      scale_color_manual(values = diagnosis_time_colors1)+#+ ggtitle(title) +
      scale_shape_manual(values = c(Train = 16,Test = 15))+
      theme_bw()+
      theme(text = element_text(size=12),
            axis.text=element_text(size=12, face = "bold"),
            axis.title=element_text(size=14,face="bold"),
            legend.position = "none")+ guides(colour = guide_legend(override.aes = list(size=4),nrow=4)) + xlab(pc2) + ylab(pc3)
    
    
    title=paste0(name,'.dx.time.pca.combined.png')
    # title=paste0(figdir,title)
    png(title,height = 1000, width = 3000,res=200)
    figure <- ggarrange(plot1, plot2,plot3,
                        labels = c(""),
                        ncol = 3, nrow = 1
                        # heights = c(0.3, 0.5),
                        # widths = c(1.5,0.7)
    )
    print(figure)
    dev.off()
    
    
    plot1 = ggplot(pca.matrix[pca.matrix$train.test == 'Test',],aes(x = PC1, y = PC2, col = Diagnosis_Time,shape = train.test)) + geom_point(alpha = 0.75,size=3) +
      scale_color_manual(values = diagnosis_time_colors1)+#+ ggtitle(title) +
      scale_shape_manual(values = c(Train = 16,Test = 15))+
      theme_bw()+
      theme(text = element_text(size=12),
            axis.text=element_text(size=12, face = "bold"),
            axis.title=element_text(size=14,face="bold"),
            legend.position = "none")+ guides(colour = guide_legend(override.aes = list(size=4),nrow=4)) + xlab(pc1) + ylab(pc2)
    plot2 = ggplot(pca.matrix[pca.matrix$train.test == 'Test',],aes(x = PC1, y = PC3, col = Diagnosis_Time,shape = train.test)) + geom_point(alpha = 0.75,size=3) +
      scale_color_manual(values = diagnosis_time_colors1)+#+ ggtitle(title) +
      scale_shape_manual(values = c(Train = 16,Test = 15))+
      theme_bw()+
      theme(text = element_text(size=12),
            axis.text=element_text(size=12, face = "bold"),
            axis.title=element_text(size=14,face="bold"),
            legend.position = "none")+ guides(colour = guide_legend(override.aes = list(size=4),nrow=4)) + xlab(pc1) + ylab(pc3)
    
    plot3 = ggplot(pca.matrix[pca.matrix$train.test == 'Test',],aes(x = PC2, y = PC3, col = Diagnosis_Time,shape = train.test)) + geom_point(alpha = 0.75,size=3) +
      scale_color_manual(values = diagnosis_time_colors1)+#+ ggtitle(title) +
      scale_shape_manual(values = c(Train = 16,Test = 15))+
      theme_bw()+
      theme(text = element_text(size=12),
            axis.text=element_text(size=12, face = "bold"),
            axis.title=element_text(size=14,face="bold"),
            legend.position = "none")+ guides(colour = guide_legend(override.aes = list(size=4),nrow=4)) + xlab(pc2) + ylab(pc3)
    
    
    title=paste0(name,'.dx.time.pca.test.png')
    # title=paste0(figdir,title)
    png(title,height = 1000, width = 3000,res=200)
    figure <- ggarrange(plot1, plot2,plot3,
                        labels = c(""),
                        ncol = 3, nrow = 1
                        # heights = c(0.3, 0.5),
                        # widths = c(1.5,0.7)
    )
    print(figure)
    dev.off()
    
    cancer.colors = c('Control' = '#BFCBC2','Pancreas' ='#09814A','Breast'='#f564a9','Prostate'= '#0094C6')
    
    plot1 = ggplot(pca.matrix, aes(x = PC1, y = PC2, col = Cancer,shape = train.test))+ geom_point(size = 3,alpha = 0.8) + 
      scale_color_manual(values = cancer.colors)+#+ ggtitle(title) +
      theme_bw()+
      scale_shape_manual(values = c(Train = 16,Test = 15))+
      theme(text = element_text(size=12),
            axis.text=element_text(size=12, face = "bold"),
            axis.title=element_text(size=14,face="bold"),
            legend.position = "none")+ guides(colour = guide_legend(override.aes = list(size=4),nrow=7)) + xlab(pc1) + ylab(pc2)+ ggtitle(paste0( ' Features: ',length(predx.dmrs.sig.hyper)))
    
    plot2 = ggplot(pca.matrix, aes(x = PC1, y = PC3, col = Cancer,shape = train.test))+ geom_point(size = 3,alpha = 0.8) + 
      scale_color_manual(values = cancer.colors)+#+ ggtitle(title) +
      theme_bw()+
      scale_shape_manual(values = c(Train = 16,Test = 15))+
      theme(text = element_text(size=12),
            axis.text=element_text(size=12, face = "bold"),
            axis.title=element_text(size=14,face="bold"),
            legend.position = "none")+ guides(colour = guide_legend(override.aes = list(size=4),nrow=7)) + xlab(pc1) + ylab(pc3)+ ggtitle(paste0( ' Features: ',length(predx.dmrs.sig.hyper)))
    
    
    plot3 = ggplot(pca.matrix, aes(x = PC2, y = PC3, col = Cancer,shape = train.test))+ geom_point(size = 3,alpha = 0.8) + 
      scale_color_manual(values = cancer.colors)+#+ ggtitle(title) +
      theme_bw()+
      scale_shape_manual(values = c(Train = 16,Test = 15))+
      theme(text = element_text(size=12),
            axis.text=element_text(size=12, face = "bold"),
            axis.title=element_text(size=14,face="bold"),
            legend.position = "none")+ guides(colour = guide_legend(override.aes = list(size=4),nrow=7)) + xlab(pc2) + ylab(pc3)+ ggtitle(paste0( ' Features: ',length(predx.dmrs.sig.hyper)))
    
    title=paste0(name,'.filler.pca.combined.png')
    # title=paste0(figdir,title)
    png(title,height = 1000, width = 3000,res=200)
    figure <- ggarrange(plot1, plot2,plot3,
                        labels = c(""),
                        ncol = 3, nrow = 1
                        # heights = c(0.3, 0.5),
                        # widths = c(1.5,0.7)
    )
    print(figure)
    dev.off()
    
    
    
    
    plot1 = ggplot(pca.matrix[pca.matrix$train.test == 'Test',], aes(x = PC1, y = PC2, col = Cancer,shape = train.test))+ geom_point(size = 3,alpha = 0.8) + 
      scale_color_manual(values = cancer.colors)+#+ ggtitle(title) +
      theme_bw()+
      scale_shape_manual(values = c(Train = 16,Test = 15))+
      theme(text = element_text(size=12),
            axis.text=element_text(size=12, face = "bold"),
            axis.title=element_text(size=14,face="bold"),
            legend.position = "none")+ guides(colour = guide_legend(override.aes = list(size=4),nrow=7)) + xlab(pc1) + ylab(pc2)+ ggtitle(paste0( ' Features: ',length(predx.dmrs.sig.hyper)))
    
    plot2 = ggplot(pca.matrix[pca.matrix$train.test == 'Test',], aes(x = PC1, y = PC3, col = Cancer,shape = train.test))+ geom_point(size = 3,alpha = 0.8) + 
      scale_color_manual(values = cancer.colors)+#+ ggtitle(title) +
      theme_bw()+
      scale_shape_manual(values = c(Train = 16,Test = 15))+
      theme(text = element_text(size=12),
            axis.text=element_text(size=12, face = "bold"),
            axis.title=element_text(size=14,face="bold"),
            legend.position = "none")+ guides(colour = guide_legend(override.aes = list(size=4),nrow=7)) + xlab(pc1) + ylab(pc3)+ ggtitle(paste0( ' Features: ',length(predx.dmrs.sig.hyper)))
    
    
    plot3 = ggplot(pca.matrix[pca.matrix$train.test == 'Test',], aes(x = PC2, y = PC3, col = Cancer,shape = train.test))+ geom_point(size = 3,alpha = 0.8) + 
      scale_color_manual(values = cancer.colors)+#+ ggtitle(title) +
      theme_bw()+
      scale_shape_manual(values = c(Train = 16,Test = 15))+
      theme(text = element_text(size=12),
            axis.text=element_text(size=12, face = "bold"),
            axis.title=element_text(size=14,face="bold"),
            legend.position = "none")+ guides(colour = guide_legend(override.aes = list(size=4),nrow=7)) + xlab(pc2) + ylab(pc3)+ ggtitle(paste0( ' Features: ',length(predx.dmrs.sig.hyper)))
    
    title=paste0(name,'.filler.pca.test.png')
    # title=paste0(figdir,title)
    png(title,height = 1000, width = 3000,res=200)
    figure <- ggarrange(plot1, plot2,plot3,
                        labels = c(""),
                        ncol = 3, nrow = 1
                        # heights = c(0.3, 0.5),
                        # widths = c(1.5,0.7)
    )
    print(figure)
    dev.off()
    
    
    
    pca.matrix$grp = gsub('_.*','',pca.matrix$GRP_Id)
    group.cols= c('#92140C','#111D4A')
    plot1 = ggplot(pca.matrix[pca.matrix$train.test == 'Test',], aes(x = PC1, y = PC2, col = grp,shape = train.test))+ geom_point(size = 3,alpha = 0.8) + 
      scale_color_manual(values = group.cols)+#+ ggtitle(title) +
      theme_bw()+
      scale_shape_manual(values = c(Train = 16,Test = 15))+
      theme(text = element_text(size=12),
            axis.text=element_text(size=12, face = "bold"),
            axis.title=element_text(size=14,face="bold"),
            legend.position = "none")+ guides(colour = guide_legend(override.aes = list(size=4),nrow=7)) + xlab(pc1) + ylab(pc2)+ ggtitle(paste0( ' Features: ',length(predx.dmrs.sig.hyper)))
    
    plot2 = ggplot(pca.matrix[pca.matrix$train.test == 'Test',], aes(x = PC1, y = PC3, col = grp,shape = train.test))+ geom_point(size = 3,alpha = 0.8) + 
      scale_color_manual(values = group.cols)+#+ ggtitle(title) +
      theme_bw()+
      scale_shape_manual(values = c(Train = 16,Test = 15))+
      theme(text = element_text(size=12),
            axis.text=element_text(size=12, face = "bold"),
            axis.title=element_text(size=14,face="bold"),
            legend.position = "none")+ guides(colour = guide_legend(override.aes = list(size=4),nrow=7)) + xlab(pc1) + ylab(pc3)+ ggtitle(paste0( ' Features: ',length(predx.dmrs.sig.hyper)))
    
    
    plot3 = ggplot(pca.matrix[pca.matrix$train.test == 'Test',], aes(x = PC2, y = PC3, col = grp,shape = train.test))+ geom_point(size = 3,alpha = 0.8) + 
      scale_color_manual(values = group.cols)+#+ ggtitle(title) +
      theme_bw()+
      scale_shape_manual(values = c(Train = 16,Test = 15))+
      theme(text = element_text(size=12),
            axis.text=element_text(size=12, face = "bold"),
            axis.title=element_text(size=14,face="bold"),
            legend.position = "none")+ guides(colour = guide_legend(override.aes = list(size=4),nrow=7)) + xlab(pc2) + ylab(pc3)+ ggtitle(paste0( ' Features: ',length(predx.dmrs.sig.hyper)))
    
    title=paste0(name,'.aix.pca.test.png')
    # title=paste0(figdir,title)
    png(title,height = 1000, width = 3000,res=200)
    figure <- ggarrange(plot1, plot2,plot3,
                        labels = c(""),
                        ncol = 3, nrow = 1
                        # heights = c(0.3, 0.5),
                        # widths = c(1.5,0.7)
    )
    print(figure)
    dev.off()
    pca.matrix$grp = gsub('_.*','',pca.matrix$GRP_Id)
    group.cols= c('#92140C','#111D4A')
    plot1 = ggplot(pca.matrix, aes(x = PC1, y = PC2, col = grp,shape = train.test))+ geom_point(size = 3,alpha = 0.8) + 
      scale_color_manual(values = group.cols)+#+ ggtitle(title) +
      theme_bw()+
      scale_shape_manual(values = c(Train = 16,Test = 15))+
      theme(text = element_text(size=12),
            axis.text=element_text(size=12, face = "bold"),
            axis.title=element_text(size=14,face="bold"),
            legend.position = "none")+ guides(colour = guide_legend(override.aes = list(size=4),nrow=7)) + xlab(pc1) + ylab(pc2)+ ggtitle(paste0( ' Features: ',length(predx.dmrs.sig.hyper)))
    
    plot2 = ggplot(pca.matrix, aes(x = PC1, y = PC3, col = grp,shape = train.test))+ geom_point(size = 3,alpha = 0.8) + 
      scale_color_manual(values = group.cols)+#+ ggtitle(title) +
      theme_bw()+
      scale_shape_manual(values = c(Train = 16,Test = 15))+
      theme(text = element_text(size=12),
            axis.text=element_text(size=12, face = "bold"),
            axis.title=element_text(size=14,face="bold"),
            legend.position = "none")+ guides(colour = guide_legend(override.aes = list(size=4),nrow=7)) + xlab(pc1) + ylab(pc3)+ ggtitle(paste0( ' Features: ',length(predx.dmrs.sig.hyper)))
    
    
    plot3 = ggplot(pca.matrix, aes(x = PC2, y = PC3, col = grp,shape = train.test))+ geom_point(size = 3,alpha = 0.8) + 
      scale_color_manual(values = group.cols)+#+ ggtitle(title) +
      theme_bw()+
      scale_shape_manual(values = c(Train = 16,Test = 15))+
      theme(text = element_text(size=12),
            axis.text=element_text(size=12, face = "bold"),
            axis.title=element_text(size=14,face="bold"),
            legend.position = "none")+ guides(colour = guide_legend(override.aes = list(size=4),nrow=7)) + xlab(pc2) + ylab(pc3)+ ggtitle(paste0( ' Features: ',length(predx.dmrs.sig.hyper)))
    
    title=paste0(name,'.aix.pca.combined.png')
    # title=paste0(figdir,title)
    png(title,height = 1000, width = 3000,res=200)
    figure <- ggarrange(plot1, plot2,plot3,
                        labels = c(""),
                        ncol = 3, nrow = 1
                        # heights = c(0.3, 0.5),
                        # widths = c(1.5,0.7)
    )
    print(figure)
    dev.off()
    
    
    plot3 = ggplot(pca.matrix, aes(x = PC2, y = PC3, col = grp,shape = train.test))+ geom_point(size = 3,alpha = 0.8) + 
      scale_color_manual(values = group.cols)+#+ ggtitle(title) +
      theme_bw()+
      scale_shape_manual(values = c(Train = 16,Test = 15))+
      theme(text = element_text(size=12),
            axis.text=element_text(size=12, face = "bold"),
            axis.title=element_text(size=14,face="bold"),
            legend.position = "right")+ guides(colour = guide_legend(override.aes = list(size=4),nrow=7)) + xlab(pc2) + ylab(pc3)+ ggtitle(paste0( ' Features: ',length(predx.dmrs.sig.hyper)))
    
    title=paste0(name,'.aix.pca.combined.labs.png')
    # title=paste0(figdir,title)
    png(title,height = 1000, width = 1500,res=200)
    
    print(plot3)
    dev.off()
    
  }
  return.df = rbind(pc.matrix,pc.matrix1)
  
  return(list(loading = return.df,variance=summary(pca.obj)$importance,weighting=pca.obj$rotation))
  
}


#coxph model
all_cause_mortality_CAN<-structure(list(x = c(0, 2.5, 7, 12, 17, 22, 27, 32, 37, 42, 47, 
                                              52, 57, 62, 67, 72, 77, 82, 87, 90), qx = c(0.0047, 2e-04, 1e-04, 
                                                                                          1e-04, 5e-04, 8e-04, 0.001, 0.0011, 0.0013, 0.0015, 0.0023, 0.0036, 
                                                                                          0.0056, 0.0088, 0.0134, 0.0216, 0.0346, 0.061, 0.1089, 0.2163
                                              ), lx = c(1e+05, 97650, 97552.35, 97503.573825, 97454.8220380875, 
                                                        97211.1849829923, 96822.3402430603, 96338.228541845, 95808.3682848649, 
                                                        95185.6138910133, 94471.7217868307, 93385.2969862821, 91704.361640529, 
                                                        89136.6395145942, 85214.6273759521, 79505.2473417633, 70918.6806288528, 
                                                        58649.7488800613, 40761.5754716426, 18566.8976273332), dx = c(470, 
                                                                                                                      19.53, 9.755235, 9.7503573825, 48.7274110190437, 77.7689479863938, 
                                                                                                                      96.8223402430603, 105.97205139603, 124.550878770324, 142.77842083652, 
                                                                                                                      217.284960109711, 336.187069150616, 513.544425186962, 784.402427728429, 
                                                                                                                      1141.87600683776, 1717.31334258209, 2453.78634975831, 3577.63468168374, 
                                                                                                                      4438.93556886188, 4016.01995679217), qx.1 = c(0.004, 1e-04, 1e-04, 
                                                                                                                                                                    1e-04, 2e-04, 3e-04, 4e-04, 5e-04, 6e-04, 9e-04, 0.0015, 0.0023, 
                                                                                                                                                                    0.0037, 0.0058, 0.0087, 0.0142, 0.0232, 0.0416, 0.0768, 0.179
                                                                                                                      ), lx.1 = c(1e+05, 98000, 97951, 97902.0245, 97853.07348775, 
                                                                                                                                  97755.2204142622, 97608.5875836408, 97413.3704084736, 97169.8369824524, 
                                                                                                                                  96878.327471505, 96442.3749978832, 95719.0571853991, 94618.288027767, 
                                                                                                                                  92867.8496992533, 90174.682057975, 86252.0833884531, 80128.1854678729, 
                                                                                                                                  70833.3159535996, 56099.9862352509, 34557.5915209146), dx.1 = c(400, 
                                                                                                                                                                                                  9.8, 9.7951, 9.79020245, 19.57061469755, 29.3265661242787, 39.0434350334563, 
                                                                                                                                                                                                  48.7066852042368, 58.3019021894714, 87.1904947243545, 144.663562496825, 
                                                                                                                                                                                                  220.153831526418, 350.087665702738, 538.633528255669, 784.519733904382, 
                                                                                                                                                                                                  1224.77958411603, 1858.97390285465, 2946.66594366975, 4308.47894286727, 
                                                                                                                                                                                                  6185.80888224371)), class = "data.frame", row.names = c(50436L, 
                                                                                                                                                                                                                                                          50442L, 50448L, 50454L, 50460L, 50466L, 50472L, 50478L, 50484L, 
                                                                                                                                                                                                                                                          50490L, 50496L, 50502L, 50508L, 50514L, 50520L, 50526L, 50532L, 
                                                                                                                                                                                                                                                          50538L, 50544L, 50550L))

age_incidence<-structure(list(Age.group = c("0 to 04", "05 to 09", "10 to 14", 
                                            "15 to 19", "20 to 24", "25 to 29", "30 to 34", "35 to 39", "40 to 44", 
                                            "45 to 49", "50 to 54", "55 to 59", "60 to 64", "65 to 69", "70 to 74", 
                                            "75 to 79", "80 to 84", "85 to 89", "90+"), PanB = c(0, 0, 0, 
                                                                                                 0, 0, 0, 0.8, 1.3, 1.6, 5, 9.7, 15.1, 27.1, 40.9, 56.5, 71.4, 
                                                                                                 75, 77, 56.2), PanM = c(0, 0, 0, 0, 0, 0.5, 0.5, 1.1, 2.2, 5.9, 
                                                                                                                         9.8, 19.2, 31.1, 46.7, 66.7, 79.2, 82, 85.1, 64.4), PanF = c(0, 
                                                                                                                                                                                      0, 0, 0.6, 0, 0, 1, 1, 1.1, 4.2, 10.2, 11.1, 23.3, 35.4, 47, 
                                                                                                                                                                                      65.8, 71, 71.6, 55.7), LungB = c(0, 0, 0.3, 0, 0.5, 0.5, 1.3, 
                                                                                                                                                                                                                       1.9, 4.7, 10, 32.1, 71.6, 132.3, 208.3, 302.4, 387.2, 373.1, 
                                                                                                                                                                                                                       302.4, 188), LungM = c(0, 0, 0, 0, 0.5, 0.5, 1, 2.2, 3.9, 9.6, 
                                                                                                                                                                                                                                              28.8, 65.8, 134.1, 217, 314.4, 423.3, 439.9, 384.8, 300.7), LungF = c(0, 
                                                                                                                                                                                                                                                                                                                    0, 0, 0, 0.6, 0.5, 1.5, 2.1, 5.4, 10.5, 35.4, 77.3, 130, 200.8, 
                                                                                                                                                                                                                                                                                                                    292.3, 357.4, 322, 248.3, 139.3), BreastB = c(0, 0, 0, 0, 0.8, 
                                                                                                                                                                                                                                                                                                                                                                  4.8, 13.6, 28, 55, 84, 101, 115, 144, 175, 208, 188, 196.7, 212, 
                                                                                                                                                                                                                                                                                                                                                                  175), BreastM = c(0, 0, 0, 0, 0, 0, 0, 0.5, 0.6, 1.1, 0.5, 1, 
                                                                                                                                                                                                                                                                                                                                                                                    2.8, 4.1, 5.5, 6.8, 6, 3.4, 0), BreastF = c(0, 0, 0, 0.6, 1.1, 
                                                                                                                                                                                                                                                                                                                                                                                                                                9.3, 26.8, 55, 109, 166, 200, 227, 280, 335, 393, 346, 345.6, 
                                                                                                                                                                                                                                                                                                                                                                                                                                349, 250.7), ProsB = c(0, 0, 0, 0, 0, 0, 0, 0, 2.7, 9.8, 40.6, 
                                                                                                                                                                                                                                                                                                                                                                                                                                                       93, 161.6, 265.1, 296.3, 278.6, 228.5, 167.4, 110.2), ProsM = c(0, 
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                       0, 0, 0, 0, 0, 0, 0, 5.6, 19.2, 81.6, 187.6, 330.4, 547.9, 618.7, 
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                       602.1, 517.9, 422.2, 365.1)), class = "data.frame", row.names = c(NA, 
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                         -19L))




aml_inc <- function(gender, x){
  if(gender==1){
    tmp = splinefun(x=c(seq(0,90,5)), y=c(cumsum(age_incidence$ProsM/100000)*5), method="mono")(x)
    
  }
  else if (gender == 0) {
    tmp = splinefun(x=c(seq(0,90,5)), y=c(cumsum(age_incidence$BreastF/100000)*5), method="mono")(x)
    
  }
  return(tmp)
}

all_surv <- function(gender, age1, age2){
  if(gender==1) {
    s <- all_cause_mortality_CAN$lx
    
  } else {
    s <- all_cause_mortality_CAN$lx.1
    
  }
  f <- function(x) exp(splinefun(all_cause_mortality_CAN$x, log(s), method="mono")(x))
  return(f(age2) / f(age1))
  
  
}

aml_inc_cr <- Vectorize(function(gender, age1, age2) sum(diff(aml_inc(gender, seq(age1,age2,1) ))*all_surv(gender, age1, seq(age1,age2-1,1)) ), c("gender","age1","age2"))



weightsf.females<-function (matched_samples){
  matched_samples_temp<- matched_samples[matched_samples$censorship_time > 0 | matched_samples$group %in% c('control','Control'),]#[na.omit(match(matched_samples$GRP_Id ,ids)),]
  latest.linkage=  as.Date('2019-01-01',format= '%Y-%m-%d')
  matched_samples_temp$gender= 1
  control.incidence = matched_samples_temp[matched_samples_temp$group %in% c('control','Control'),]
  
  expected_rate_breast_cr <- mean(aml_inc_cr(matched_samples_temp$gender, matched_samples_temp$SDC_AGE_CALC, matched_samples_temp$SDC_AGE_CALC+pmax(1, matched_samples_temp$censorship_time/365))[matched_samples_temp$group %in% c('control','Control')])
  
  n_total_breast <- sum(!matched_samples_temp $group %in% c("Control",'control'))/expected_rate_breast_cr
  n_total_breast
  weights <- rep(1, nrow(matched_samples_temp))
  weights[matched_samples_temp$group %in% c("Control",'control')] <- n_total_breast/sum(matched_samples_temp$group %in% c("Control",'control'))
  
  return(weights)
}

coxph.calling.all = function(fragment.df.ratio.raw.t,merged.df.filt,filler = T,ncores=5) {
  library(survival)
  library(parallel)
  #fragment.df.ratio.raw.t #matrix of methylation counts per region (rows) across samples (columns)
  #merged.df.filt #sample information with covariates
  #ncores #number of cores for parallel computing
  combined.df.cox = cbind(merged.df.filt,fragment.df.ratio.raw.t) #combining count matrix wiith sample information table
  combined.df.cox = combined.df.cox[combined.df.cox$censorship_time > 0 ,]
  combined.df.cox$event = ifelse(combined.df.cox$group %in% c('control','Control'),0,1)
  
  windows = colnames(fragment.df.ratio.raw.t) #list of regions or cpg sites
  
  combined.df.cox$TOTAL = combined.df.cox$total/1000000
  #running coxph per region or cpg site using lapply ran in parallel
  return.list= mclapply(windows, function(z) {
    targ.df= combined.df.cox[,c('group','SDC_AGE_CALC','TOTAL','event','censorship_time','THALIANA_BETA','filler',z)] #
    targ.df$methyl = targ.df[,z]
    targ.df$filler = factor(targ.df$filler, levels= c('UFiller','MFiller'))
    if (sum(is.na(targ.df$methyl)) == 0){
      if (filler == T){
        test = summary(coxph(Surv(as.numeric(targ.df$censorship_time), targ.df$event) ~ TOTAL+ THALIANA_BETA+ filler  + methyl,targ.df)) #
        
      } else {
        test = summary(coxph(Surv(as.numeric(targ.df$censorship_time), targ.df$event) ~ TOTAL +THALIANA_BETA + methyl,targ.df,weights = female.weights)) #
        
      }
      return.df = data.frame(window =z, pvalue = test$coefficients[nrow(test$coefficients),6], HR =  test$coefficients[nrow(test$coefficients),2],HR.SE = test$coefficients[nrow(test$coefficients),3])
      return(return.df)
    }
    
    
  },mc.cores = ncores)
  res.df=do.call('rbind',return.list)
  return(res.df)
}

predictive.models = function(targ.matrix1,merged.df,train.set,test.set,targ.features,feats,feature.type,stat.test,no_cores=5, standardize = F) {
  library(caret)
  library(glmnet)
  library(ROCR)
  auc_calc = function(prediction_table,labels = c('Control','Cancer')) {
    tmp = prediction_table
    tmp = tmp[order(-tmp$methylation_score),]
    #labels = c(tmp$reported)
    pred = prediction(predictions = c(tmp$methylation_score) ,labels =  tmp$reported, labels)
    perf_AUC=performance(pred,"auc") #Calculate the AUC value
    AUC=perf_AUC@y.values[[1]]
    return(AUC)
  }
  
  prediction.setup.lasso= function(lambda.min,lasso.fit,test_set_matrix,features,model = 'logreg') {
    lasso.predict.prob = predict(lasso.fit, s='lambda.min', newx=as.matrix(test_set_matrix[,c(features)]), type="response") # check help(predict.glmnet)
    predictions = ifelse(lasso.predict.prob > 0.5,'Cancer','Control')
    prediction_table = data.frame(GRP_Id=test_set_matrix$GRP_Id, 'predictions' = predictions[,1], reported = test_set_matrix$group, methylation_score = lasso.predict.prob[,1])
    prediction_table = prediction_table[order(-prediction_table$methylation_score),]
    prediction_table$auroc = auc_calc(prediction_table,c('Control','Cancer'))
    prediction_table$model = 'logreg'
    #train_performance = getTrainPerf(lasso.fit)
    #prediction_table$TrainROC = train_performance$TrainROC
    #prediction_table$TrainSens = train_performance$TrainSens
    #prediction_table$TrainSpec = train_performance$TrainSpec
    
    
    return(prediction_table)
  }
  prediction.setup.rf= function(rf_model,test_set_matrix,features,model = 'RF') {
    
    predictions = predict(rf_model, test_set_matrix)
    predictions_prob = predict(rf_model, test_set_matrix[,features], type = 'prob')
    prediction_table = data.frame(GRP_Id= rownames(test_set_matrix), predictions = predictions, reported = test_set_matrix$group, methylation_score = predictions_prob[,colnames(predictions_prob) != 'Control'], model = model)
    prediction_table = prediction_table[order(-prediction_table$methylation_score),]
    prediction_table$auroc = auc_calc(prediction_table,c('Control','Cancer'))
    prediction_table$model = 'rf'
    #train_performance = getTrainPerf(rf_model)
    #prediction_table$TrainROC = train_performance$TrainROC
    # prediction_table$TrainSens = train_performance$TrainSens
    # prediction_table$TrainSpec = train_performance$TrainSpec
    
    return(prediction_table)
  }
  prediction.setup.rfsurv= function(rf.surv,targ.matrix.test,model = 'RFSurv') {
    
    surv.test.set = targ.matrix.test
    surv.test.set$event =ifelse(surv.test.set$group == 'Control',0,1)
    surv.test.set$censorship_time = test.set[match(surv.test.set$GRP_Id,test.set$GRP_Id),'censorship_time']
    
    o.pred <- predict(object = rf.surv, surv.test.set)
    
    predictions = ifelse(o.pred$predicted >median(o.pred$predicted ),'Cancer','Control')
    prediction_table = data.frame(GRP_Id=surv.test.set$GRP_Id, 
                                  'predictions' = predictions, reported = targ.matrix.test$group, 
                                  methylation_score = o.pred$predicted)
    prediction_table = prediction_table[order(-prediction_table$methylation_score),]
    prediction_table$auroc = auc_calc(prediction_table,c('Control','Cancer'))
    prediction_table$model = 'RFSurv'
    #train_performance = getTrainPerf(lasso.fit)
    #prediction_table$TrainROC = train_performance$TrainROC
    #prediction_table$TrainSens = train_performance$TrainSens
    #prediction_table$TrainSpec = train_performance$TrainSpec
    
    
    return(prediction_table)
  }
  
  results.df.all = NULL
  feature.weights.all = NULL
  require(doMC)
  library(doParallel)
  cl <- makeCluster(ncores)
  registerDoParallel(cl)
  for (fl in feats){
    start= Sys.time()
    results.df = NULL
    f = min(fl, length(targ.features))
    
    #all features as features
    targ.features1 = targ.features[1:f]
    targ.matrix = targ.matrix1[,targ.features[1:f]]
    
    targ.matrix$GRP_Id = rownames(targ.matrix)
    if (standardize == T){
      targ.matrix[,targ.features1] = apply(targ.matrix[,targ.features1],2,scale)
      rownames(targ.matrix) = targ.matrix$GRP_Id
    }
    targ.matrix=  merge(targ.matrix,merged.df,by='GRP_Id')
    
    
    targ.matrix.train = targ.matrix[targ.matrix$GRP_Id %in% train.set$GRP_Id,]
    targ.matrix.test = targ.matrix[targ.matrix$GRP_Id %in% test.set$GRP_Id,]
    
    rownames(targ.matrix.train)= targ.matrix.train$GRP_Id
    rownames(targ.matrix.test) = targ.matrix.test$GRP_Id
    
    targ.matrix.train$group = factor(ifelse(targ.matrix.train$group == 'Control','Control','Cancer'),levels = c('Control','Cancer'))
    targ.matrix.test$group = factor(ifelse(targ.matrix.test$group == 'Control','Control','Cancer'),levels=  c('Control','Cancer'))
    
    
    #logreg new
    ##
    alpha_grid <- seq(0, 1, by = 0.1)
    
    # Initialize variables to store the results
    cv_errors <- numeric(length(alpha_grid))
    best_lambdas <- numeric(length(alpha_grid))
    
    # Perform grid search over alpha
    for (i in seq_along(alpha_grid)) {
      alpha <- alpha_grid[i]
      #targ.matrix.train will be your train set with features for each column and sample acroww each row
      cv_fit <- cv.glmnet(as.matrix(targ.matrix.train[,c(targ.features[1:f])]), targ.matrix.train$group, family="binomial", type.measure="class",  nfolds = 10,parallel=F,alpha=alpha)
      cv_errors[i] <- min(cv_fit$cvm)  # cvm is the mean cross-validated error
      best_lambdas[i] <- cv_fit$lambda.min
    }
    
    # Find the alpha with the smallest cross-validated error
    optimal_index <- which.min(cv_errors)
    optimal_alpha <- alpha_grid[optimal_index]
    optimal_lambda <- best_lambdas[optimal_index]
    
    
    # Refit the model using the optimal alpha and lambda
    best_fit <- glmnet(as.matrix(targ.matrix.train[,c(targ.features[1:f])]), y=as.numeric(targ.matrix.train$group), alpha = optimal_alpha, lambda = optimal_lambda, family = "binomial", type.measure = "class")
    #test_set_matrix will be your test set matrix (features per column, samples per row) formated the same way as the train set matrix 
    lasso.predict.prob = predict(best_fit, s='lambda.min', newx=as.matrix(targ.matrix.test[,c(targ.features[1:f])]), type="response") 
    
    prediction_table.logreg.predx.predictors.norelh =prediction.setup.lasso(optimal_lambda,best_fit,targ.matrix.test,c(targ.features[1:f]),model = 'logreg')
    prediction_table.logreg.predx.predictors.norelh$model = c('logreg.new')
    results.df = rbind(results.df,prediction_table.logreg.predx.predictors.norelh)
    feature.coef =coef(best_fit,optimal_lambda)
    lasso.feature.weights =data.frame(PC =rownames(feature.coef), coef= feature.coef[1:length(feature.coef)])
    lasso.feature.weights$model = 'logreg.new'
    
    
    #logreg old
    
    prediction.setup.lasso.prev= function(lambda.min,lasso.fit,test_set_matrix,features,model = 'logreg') {
      lasso.predict.prob = predict(lasso.fit, s='lambda.min', newx=as.matrix(test_set_matrix[,c(features)]), type="response") # check help(predict.glmnet)
      predictions = ifelse(lasso.predict.prob > 0.5,'Cancer','Control')
      prediction_table = data.frame(GRP_Id=test_set_matrix$GRP_Id, 'predictions' = predictions[,1], reported = test_set_matrix$group, methylation_score = lasso.predict.prob[,1])
      prediction_table = prediction_table[order(-prediction_table$methylation_score),]
      prediction_table$auroc = auc_calc(prediction_table,c('Control','Cancer'))
      prediction_table$model = 'logreg.old'
      #train_performance = getTrainPerf(lasso.fit)
      #prediction_table$TrainROC = train_performance$TrainROC
      #prediction_table$TrainSens = train_performance$TrainSens
      #prediction_table$TrainSpec = train_performance$TrainSpec
      
      
      return(prediction_table)
    }
    #glmnet
    ##
    MSEs = NULL
    for (i in 1:10){
      lasso.cv = cv.glmnet(as.matrix(targ.matrix.train[,c(targ.features[1:f])]), targ.matrix.train$group, family="binomial", type.measure="class",  nfolds = 10)
      MSEs <- cbind(MSEs, lasso.cv$cvm)
    }
    rownames(MSEs) <- lasso.cv$lambda
    lambda.min <- as.numeric(names(which.min(rowMeans(MSEs))))
    
    lasso.fit = glmnet(as.matrix(targ.matrix.train[,c(targ.features[1:f])]), 
                       targ.matrix.train$group, 
                       family="binomial",
                       nfolds=10,
                       parallel=F,
                       standardize=F,
                       lambda= lambda.min)
    
    #lasso.fit = cv.glmnet(as.matrix(targ.matrix.train[,c(targ.features)]), targ.matrix.train$group, family="binomial")
    prediction_table.logreg.old =prediction.setup.lasso.prev(lambda.min,lasso.fit,targ.matrix.test,c(targ.features[1:f]),model = 'logreg')
    
    prediction_table.logreg.old$model = c('logreg.old')
    results.df = rbind(results.df,prediction_table.logreg.old)
    feature.coef =coef(lasso.fit,lambda.min)
    lasso.old.feature.weights =data.frame(PC =rownames(feature.coef), coef= feature.coef[1:length(feature.coef)])
    lasso.old.feature.weights$model = 'logreg.old'
    
    
    #rf
    start1=Sys.time()
    
    control = trainControl(method = 'repeatedcv', number = 10, repeats = 10, search = 'random', classProbs = TRUE, summaryFunction = twoClassSummary)
    tunegrid <- expand.grid(.mtry=seq(10,30,10))
    metric = 'Accuracy'
    
    ntrees = c(1000)
    ntree_list = list()
    acc = c()
    targ.matrix.train$group = factor(targ.matrix.train$group, levels = c('Control','Cancer'))
    rf_model = train(group ~ ., data = targ.matrix.train[,c('group',targ.features[1:f])], method = 'rf', metric = metric, tuneGrid = tunegrid, trControl = control, ntrees = ntrees)
    
    
    rf.performance = prediction.setup.rf(rf_model,targ.matrix.test,targ.features[1:f],model='rf')
    rf.performance$model ='RF'
    rf.feature.weights=varImp(rf_model)$importance
    rf.feature.weights=data.frame(PC = rownames(rf.feature.weights), 
                                  coef = rf.feature.weights[,1],
                                  model = 'RF')
    
    results.df = rbind(results.df,rf.performance)
    end1=Sys.time()
    
    
    #random forest survival
    library(randomForestSRC)
    
    surv.train.set = targ.matrix.train
    surv.train.set$event =ifelse(surv.train.set$group == 'Control',0,1)
    surv.train.set$censorship_time = train.set[match(surv.train.set$GRP_Id,train.set$GRP_Id),'censorship_time']
    
    
    o <- tune(Surv(censorship_time, event) ~ ., data = surv.train.set[,c('censorship_time','event',targ.features[1:f])] )
    rf.surv = rfsrc(Surv(censorship_time, event) ~ ., 
                    data = surv.train.set[,c('censorship_time','event',targ.features[1:f])] ,
                    ntree=500,
                    mtry = o$optimal[2],
                    nodesize=o$optimal[1])  
    
    
    
    rf.surv.predictions = prediction.setup.rfsurv(rf.surv,targ.matrix.test)
    rf.surv.predictions$model = 'RFSurv'
    results.df = rbind(results.df,rf.surv.predictions)
    
    rfsurv.feature.weights <- vimp(rf.surv)$importance
    rfsurv.feature.weights=data.frame(PC = names(rfsurv.feature.weights), 
                                      coef = rfsurv.feature.weights,
                                      model = 'RFSurv')
    
    #ncvsurv
    library(ncvreg)
    library(survival)
    train.fu = train.set[order(match(train.set$GRP_Id, rownames(targ.matrix.train))),'censorship_time']
    train.event = train.set[order(match(train.set$GRP_Id, rownames(targ.matrix.train))),'Cancer']
    
    tmp.surv = data.frame(fu.time = train.fu, event= ifelse(train.event == 'Control',0,1))
    
    ncvsurv.weights = NULL
    for (p in c('MCP','SCAD','lasso')) {
      print(paste0(p,'.',f))
      alpha_grid <- seq(0.1, 1, by = 0.1)
      
      # Initialize variables to store the results
      cv_errors <- numeric(length(alpha_grid))
      best_lambdas <- numeric(length(alpha_grid))
      
      # Perform grid search over alpha
      for (i in seq_along(alpha_grid)) {
        alpha <- alpha_grid[i]
        # Perform cross-validation
        cv_fit <- cv.ncvsurv(as.matrix(targ.matrix.train[,targ.features[1:f]]),tmp.surv, penalty = p, alpha = alpha)
        cv_errors[i] <- min(cv_fit$cve)  # cve is the cross-validated error
        best_lambdas[i] <- cv_fit$lambda.min
      }
      
      # Find the alpha with the smallest cross-validated error
      optimal_index <- which.min(cv_errors)
      optimal_alpha <- alpha_grid[optimal_index]
      optimal_lambda <- best_lambdas[optimal_index]
      
      # Print the optimal alpha and lambda
      cat("Optimal alpha:", optimal_alpha, "\n")
      cat("Optimal lambda:", optimal_lambda, "\n")
      
      # Refit the model using the optimal alpha and lambda
      best_fit <- ncvsurv(as.matrix(targ.matrix.train[,targ.features[1:f]]),tmp.surv, penalty = p, alpha = optimal_alpha, lambda = optimal_lambda)
      prediction.setup.ncvsurv= function(lasso.fit,test_set_matrix,features,model = p) {
        lasso.predict.prob = predict(lasso.fit,as.matrix(test_set_matrix[,c(features)])) # check help(predict.glmnet)
        predictions = ifelse(lasso.predict.prob > median(lasso.predict.prob),'Cancer','Control')
        prediction_table = data.frame(GRP_Id=test_set_matrix$GRP_Id, 'predictions' = predictions, reported = test_set_matrix$group, methylation_score = lasso.predict.prob)
        prediction_table = prediction_table[order(-prediction_table$methylation_score),]
        prediction_table$auroc = auc_calc(prediction_table,c('Control','Cancer'))
        prediction_table$model = paste0('NCVSurv.',model)
        #train_performance = getTrainPerf(lasso.fit)
        #prediction_table$TrainROC = train_performance$TrainROC
        #prediction_table$TrainSens = train_performance$TrainSens
        #prediction_table$TrainSpec = train_performance$TrainSpec
        
        
        return(prediction_table)
      }
      ncvsurv.prediction = prediction.setup.ncvsurv(best_fit, targ.matrix.test,targ.features[1:f],model = p)
      results.df = rbind(results.df,ncvsurv.prediction)
      
      ncvsurv.feature.weights <- coef(best_fit)
      ncvsurv.feature.weights=data.frame(PC = names(ncvsurv.feature.weights), 
                                         coef = ncvsurv.feature.weights,
                                         model =paste0('NCVSurv.',p))
      ncvsurv.weights = rbind(ncvsurv.weights,ncvsurv.feature.weights)
      
    }
    
    
    #annotating scores
    results.df$fold = fold
    results.df$seed = seedno
    results.df$feature = feature.type
    results.df$test = stat.test
    results.df$n.features = fl
    results.df.all = rbind(results.df.all,results.df)
    
    #combining feature weights
    feature.weights = rbind(lasso.feature.weights,lasso.old.feature.weights, rf.feature.weights,rfsurv.feature.weights,ncvsurv.weights)
    feature.weights$fold = fold
    feature.weights$seed = seedno
    feature.weights$feature = feature.type
    feature.weights$test = stat.test
    feature.weights$n.features = fl
    feature.weights.all = rbind(feature.weights.all,feature.weights)
    end = Sys.time()
    
    print(start-end)
  }
  
  return.list = list(predictions = results.df.all,
                     weights = feature.weights.all)
  
  return(return.list)
}


predictive.models = function(targ.matrix1,
                             merged.df,
                             train.set,
                             test.set,
                             targ.features,
                             feats,
                             feature.type,
                             stat.test,no_cores=5, standardize = F) {
  library(caret)
  library(glmnet)
  library(ROCR)
  auc_calc = function(prediction_table,labels = c('Control','Cancer')) {
    tmp = prediction_table
    tmp = tmp[order(-tmp$methylation_score),]
    #labels = c(tmp$reported)
    pred = prediction(predictions = c(tmp$methylation_score) ,labels =  tmp$reported, labels)
    perf_AUC=performance(pred,"auc") #Calculate the AUC value
    AUC=perf_AUC@y.values[[1]]
    return(AUC)
  }
  
  prediction.setup.lasso= function(lambda.min,lasso.fit,test_set_matrix,features,model = 'logreg') {
    targ.matrix.test.model = model.matrix(~. -1,targ.matrix.test[,features] )
    
    lasso.predict.prob = predict(lasso.fit, s='lambda.min', newx=targ.matrix.test.model, type="response") # check help(predict.glmnet)
    predictions = ifelse(lasso.predict.prob > 0.5,'Cancer','Control')
    prediction_table = data.frame(GRP_Id=test_set_matrix$GRP_Id, 'predictions' = predictions[,1], reported = test_set_matrix$group, methylation_score = lasso.predict.prob[,1])
    prediction_table = prediction_table[order(-prediction_table$methylation_score),]
    prediction_table$auroc = auc_calc(prediction_table,c('Control','Cancer'))
    prediction_table$model = 'logreg'
    #train_performance = getTrainPerf(lasso.fit)
    #prediction_table$TrainROC = train_performance$TrainROC
    #prediction_table$TrainSens = train_performance$TrainSens
    #prediction_table$TrainSpec = train_performance$TrainSpec
    
    
    return(prediction_table)
  }
  prediction.setup.rf= function(rf_model,test_set_matrix,features,model = 'RF') {
    
    predictions = predict(rf_model, test_set_matrix)
    #targ.matrix.test.model = model.matrix( ~. -1,test_set_matrix[,c(targ.features[1:f])])
    
    predictions_prob = predict(rf_model,test_set_matrix[,features], type = 'prob')
    prediction_table = data.frame(GRP_Id= rownames(test_set_matrix), predictions = predictions, reported = test_set_matrix$group, methylation_score = predictions_prob[,colnames(predictions_prob) != 'Control'], model = model)
    prediction_table = prediction_table[order(-prediction_table$methylation_score),]
    prediction_table$auroc = auc_calc(prediction_table,c('Control','Cancer'))
    prediction_table$model = 'rf'
    #train_performance = getTrainPerf(rf_model)
    #prediction_table$TrainROC = train_performance$TrainROC
    # prediction_table$TrainSens = train_performance$TrainSens
    # prediction_table$TrainSpec = train_performance$TrainSpec
    
    return(prediction_table)
  }
  prediction.setup.rfsurv= function(rf.surv,targ.matrix.test,model = 'RFSurv') {
    
    surv.test.set = targ.matrix.test
    surv.test.set$event =ifelse(surv.test.set$group == 'Control',0,1)
    surv.test.set$censorship_time = test.set[match(surv.test.set$GRP_Id,test.set$GRP_Id),'censorship_time']
    
    o.pred <- predict(object = rf.surv, surv.test.set)
    
    predictions = ifelse(o.pred$predicted >median(o.pred$predicted ),'Cancer','Control')
    prediction_table = data.frame(GRP_Id=surv.test.set$GRP_Id, 
                                  'predictions' = predictions, reported = targ.matrix.test$group, 
                                  methylation_score = o.pred$predicted)
    prediction_table = prediction_table[order(-prediction_table$methylation_score),]
    prediction_table$auroc = auc_calc(prediction_table,c('Control','Cancer'))
    prediction_table$model = 'RFSurv'
    #train_performance = getTrainPerf(lasso.fit)
    #prediction_table$TrainROC = train_performance$TrainROC
    #prediction_table$TrainSens = train_performance$TrainSens
    #prediction_table$TrainSpec = train_performance$TrainSpec
    
    
    return(prediction_table)
  }
  
  results.df.all = NULL
  feature.weights.all = NULL
  require(doMC)
  library(doParallel)
  cl <- makeCluster(ncores)
  registerDoParallel(cl)
  for (fl in feats){
    start= Sys.time()
    results.df = NULL
    f = min(fl, length(targ.features))
    
    #all features as features
    targ.features1 = targ.features[1:f]
    targ.matrix = targ.matrix1[,targ.features[1:f]]
    
    targ.matrix$GRP_Id = rownames(targ.matrix)
    if (standardize == T){
      targ.matrix[,targ.features1] = apply(targ.matrix[,targ.features1],2,scale)
      rownames(targ.matrix) = targ.matrix$GRP_Id
    }
    targ.matrix=  merge(targ.matrix,merged.df,by='GRP_Id')
    
    
    targ.matrix.train = targ.matrix[targ.matrix$GRP_Id %in% train.set$GRP_Id,]
    targ.matrix.test = targ.matrix[targ.matrix$GRP_Id %in% test.set$GRP_Id,]
    
    rownames(targ.matrix.train)= targ.matrix.train$GRP_Id
    rownames(targ.matrix.test) = targ.matrix.test$GRP_Id
    
    targ.matrix.train$group = factor(ifelse(targ.matrix.train$group == 'Control','Control','Cancer'),levels = c('Control','Cancer'))
    targ.matrix.test$group = factor(ifelse(targ.matrix.test$group == 'Control','Control','Cancer'),levels=  c('Control','Cancer'))
    
    
    targ.matrix.train.model = model.matrix(~. -1,targ.matrix.train[,c(targ.features[1:f])] )
    targ.matrix.test.model = model.matrix(~. -1,targ.matrix.test[,c(targ.features[1:f])] )
    
    #logreg new
    ##
    alpha_grid <- seq(0, 1, by = 0.1)
    
    # Initialize variables to store the results
    cv_errors <- numeric(length(alpha_grid))
    best_lambdas <- numeric(length(alpha_grid))
    
    # Perform grid search over alpha
    for (i in seq_along(alpha_grid)) {
      alpha <- alpha_grid[i]
      #targ.matrix.train will be your train set with features for each column and sample acroww each row
      #cv_fit <- cv.glmnet(as.matrix(targ.matrix.train[,c(targ.features[1:f])]), targ.matrix.train$group, family="binomial", type.measure="class",  nfolds = 10,parallel=F,alpha=alpha)
      cv_fit <- cv.glmnet(targ.matrix.train.model, targ.matrix.train$group, family="binomial", type.measure="class",  nfolds = 10,parallel=F,alpha=alpha)
      
      cv_errors[i] <- min(cv_fit$cvm)  # cvm is the mean cross-validated error
      best_lambdas[i] <- cv_fit$lambda.min
    }
    
    # Find the alpha with the smallest cross-validated error
    optimal_index <- which.min(cv_errors)
    optimal_alpha <- alpha_grid[optimal_index]
    optimal_lambda <- best_lambdas[optimal_index]
    
    
    # Refit the model using the optimal alpha and lambda
    best_fit <- glmnet(targ.matrix.train.model, y=as.numeric(targ.matrix.train$group), alpha = optimal_alpha, lambda = optimal_lambda, family = "binomial", type.measure = "class")
    #test_set_matrix will be your test set matrix (features per column, samples per row) formated the same way as the train set matrix 
    lasso.predict.prob = predict(best_fit, s='lambda.min', newx=targ.matrix.test.model, type="response") 
    
    prediction_table.logreg.predx.predictors.norelh =prediction.setup.lasso(optimal_lambda,best_fit,targ.matrix.test,targ.features[1:f],model = 'logreg')
    prediction_table.logreg.predx.predictors.norelh$model = c('logreg.new')
    results.df = rbind(results.df,prediction_table.logreg.predx.predictors.norelh)
    feature.coef =coef(best_fit,optimal_lambda)
    lasso.feature.weights =data.frame(PC =rownames(feature.coef), coef= feature.coef[1:length(feature.coef)])
    lasso.feature.weights$model = 'logreg.new'
    
    
    #logreg old
    
    prediction.setup.lasso.prev= function(lambda.min,lasso.fit,test_set_matrix,features,model = 'logreg') {
      test_set_matrix.model= model.matrix(~. -1,targ.matrix.test[,features] )
      lasso.predict.prob = predict(lasso.fit, s=lambda.min, newx=test_set_matrix.model, type="response") # check help(predict.glmnet)
      predictions = ifelse(lasso.predict.prob > 0.5,'Cancer','Control')
      prediction_table = data.frame(GRP_Id=test_set_matrix$GRP_Id, 'predictions' = predictions[,1], reported = test_set_matrix$group, methylation_score = lasso.predict.prob[,1])
      prediction_table = prediction_table[order(-prediction_table$methylation_score),]
      prediction_table$auroc = auc_calc(prediction_table,c('Control','Cancer'))
      prediction_table$model = 'logreg.old'
      #train_performance = getTrainPerf(lasso.fit)
      #prediction_table$TrainROC = train_performance$TrainROC
      #prediction_table$TrainSens = train_performance$TrainSens
      #prediction_table$TrainSpec = train_performance$TrainSpec
      
      
      return(prediction_table)
    }
    #glmnet
    ##
    MSEs = NULL
    for (i in 1:10){
      lasso.cv = cv.glmnet(targ.matrix.train.model, targ.matrix.train$group, family="binomial", type.measure="class",  nfolds = 10)
      MSEs <- cbind(MSEs, lasso.cv$cvm)
    }
    rownames(MSEs) <- lasso.cv$lambda
    lambda.min <- as.numeric(names(which.min(rowMeans(MSEs))))
    
    lasso.fit = glmnet(as.matrix(targ.matrix.train.model), 
                       as.numeric(targ.matrix.train$group), 
                       family="binomial",
                       nfolds=10,
                       lambda= lambda.min,
                       type.measure = "class")
    
    #lasso.fit = cv.glmnet(targ.matrix.train.model, targ.matrix.train$group, family="binomial")
    lambda.min = lasso.fit$lambda.min
    prediction_table.logreg.old =prediction.setup.lasso.prev(lasso.fit$lambda.min,lasso.fit,targ.matrix.test,c(targ.features[1:f]),model = 'logreg')
    
    prediction_table.logreg.old$model = c('logreg.old')
    results.df = rbind(results.df,prediction_table.logreg.old)
    feature.coef =coef(best_fit,optimal_lambda)
    lasso.old.feature.weights =data.frame(PC =rownames(feature.coef), coef= feature.coef[1:length(feature.coef)])
    lasso.old.feature.weights$model = 'logreg.old'
    
    #
    prediction.setup.lasso.cv= function(lasso.fit,lambda.min,test_set_matrix,features,model = 'logreg') {
      test_set_matrix.model= model.matrix(~. -1,targ.matrix.test[,features] )
      
      lasso.predict.prob = predict(lasso.fit, s=lambda.min, newx=test_set_matrix.model, type="response") # check help(predict.glmnet)
      predictions = ifelse(lasso.predict.prob > 0.5,'Cancer','Control')
      prediction_table = data.frame(GRP_Id=test_set_matrix$GRP_Id, 'predictions' = predictions[,1], reported = test_set_matrix$group, methylation_score = lasso.predict.prob[,1])
      prediction_table = prediction_table[order(-prediction_table$methylation_score),]
      prediction_table$auroc = auc_calc(prediction_table,c('Control','Cancer'))
      prediction_table$model = 'logreg'
      #train_performance = getTrainPerf(lasso.fit)
      #prediction_table$TrainROC = train_performance$TrainROC
      #prediction_table$TrainSens = train_performance$TrainSens
      #prediction_table$TrainSpec = train_performance$TrainSpec
      
      
      return(prediction_table)
    }
    
    MSEs = NULL
    for (i in 1:10){
      print(i)
      lasso.cv = cv.glmnet(targ.matrix.train.model, targ.matrix.train$group, family="binomial", type.measure="class",  nfolds = 5,parallel=T)
      MSEs <- cbind(MSEs, lasso.cv$cvm)
    }
    rownames(MSEs) <- lasso.cv$lambda
    lambda.min <- as.numeric(names(which.min(rowMeans(MSEs))))
    
    lasso.fit = cv.glmnet(targ.matrix.train.model, 
                          targ.matrix.train$group, 
                          family="binomial",
                          nfolds=10,
                          parallel=T)
    
    feature.coef =coef(lasso.fit,lambda.min) 
    lasso.feature.weights =data.frame(PC =rownames(feature.coef), coef= feature.coef[1:length(feature.coef)])
    lasso.feature.weights$model = 'logreg'
    # lambda.minUsing cross validation glmnet
    lambda.min = lasso.fit$lambda.min
    prediction_table.logreg.predx.predictors.norelh =prediction.setup.lasso.cv(lasso.fit,lambda.min,targ.matrix.test,c(targ.features[1:f]),model = 'logreg')
    prediction_table.logreg.predx.predictors.norelh$model = c('logreg.cv')
    results.df = rbind(results.df,prediction_table.logreg.predx.predictors.norelh)
    
    
    #rf
    start1=Sys.time()
    
    control = trainControl(method = 'repeatedcv', number = 10, repeats = 10, search = 'random', classProbs = TRUE, summaryFunction = twoClassSummary)
    tunegrid <- expand.grid(.mtry=seq(10,30,10))
    metric = 'Accuracy'
    
    ntrees = c(1000)
    ntree_list = list()
    acc = c()
    
    targ.matrix.train$group =factor(targ.matrix.train$group, levels = c('Control','Cancer'))
    targ.matrix.train.model = model.matrix(group ~. -1,targ.matrix.train[,c('group',targ.features[1:f])] )
    targ.matrix.train.model$group = targ.matrix.train$group
    targ.matrix.test.model = model.matrix( ~. -1,targ.matrix.test[,c(targ.features[1:f])])
    #targ.matrix.test.model$group = targ.matrix.test.model$group
    
    #targ.matrix.train$group = factor(targ.matrix.train$group, levels = c('Control','Cancer'))
    rf_model = train(group ~ . -1, data = targ.matrix.train[,c('group',targ.features[1:f])] , method = 'rf', metric = metric, tuneGrid = tunegrid, trControl = control, ntrees = ntrees)
    
    
    rf.performance = prediction.setup.rf(rf_model,targ.matrix.test,targ.features[1:f],model='rf')
    rf.performance$model ='RF'
    rf.feature.weights=varImp(rf_model)$importance
    rf.feature.weights=data.frame(PC = rownames(rf.feature.weights), 
                                  coef = rf.feature.weights[,1],
                                  model = 'RF')
    
    results.df = rbind(results.df,rf.performance)
    end1=Sys.time()
    
    
    
    #annotating scores
    results.df$fold = fold
    results.df$seed = seedno
    results.df$feature = feature.type
    results.df$test = stat.test
    results.df$n.features = fl
    results.df.all = rbind(results.df.all,results.df)
    
    #combining feature weights
    feature.weights = rbind(lasso.feature.weights,lasso.old.feature.weights, rf.feature.weights)
    feature.weights$fold = fold
    feature.weights$seed = seedno
    feature.weights$feature = feature.type
    feature.weights$test = stat.test
    feature.weights$n.features = fl
    feature.weights.all = rbind(feature.weights.all,feature.weights)
    end = Sys.time()
    
    print(start-end)
  }
  
  return.list = list(predictions = results.df.all,
                     weights = feature.weights.all)
  
  return(return.list)
}


options(expressions = 5e5) 


predictive.models = function(targ.matrix1, #modified to allow covariates
                             merged.df,
                             train.set,
                             test.set,
                             targ.features,
                             feats,
                             feature.type,
                             stat.test,no_cores=5, standardize = F) {
  library(caret)
  library(glmnet)
  library(ROCR)
  auc_calc = function(prediction_table,labels = c('Control','Cancer')) {
    tmp = prediction_table
    tmp = tmp[order(-tmp$methylation_score),]
    #labels = c(tmp$reported)
    pred = prediction(predictions = c(tmp$methylation_score) ,labels =  tmp$reported, labels)
    perf_AUC=performance(pred,"auc") #Calculate the AUC value
    AUC=perf_AUC@y.values[[1]]
    return(AUC)
  }
  
  prediction.setup.lasso= function(lambda.min,lasso.fit,test_set_matrix,features,model = 'logreg') {
    targ.matrix.test.model = model.matrix(~. -1,targ.matrix.test[,features] )
    
    lasso.predict.prob = predict(lasso.fit, s='lambda.min', newx=targ.matrix.test.model, type="response") # check help(predict.glmnet)
    predictions = ifelse(lasso.predict.prob > 0.5,'Cancer','Control')
    prediction_table = data.frame(GRP_Id=test_set_matrix$GRP_Id, 'predictions' = predictions[,1], reported = test_set_matrix$group, methylation_score = lasso.predict.prob[,1])
    prediction_table = prediction_table[order(-prediction_table$methylation_score),]
    prediction_table$auroc = auc_calc(prediction_table,c('Control','Cancer'))
    prediction_table$model = 'logreg'
    #train_performance = getTrainPerf(lasso.fit)
    #prediction_table$TrainROC = train_performance$TrainROC
    #prediction_table$TrainSens = train_performance$TrainSens
    #prediction_table$TrainSpec = train_performance$TrainSpec
    
    
    return(prediction_table)
  }
  prediction.setup.rf= function(rf_model,test_set_matrix,features,model = 'RF') {
    
    predictions = predict(rf_model, test_set_matrix)
    #targ.matrix.test.model = model.matrix( ~. -1,test_set_matrix[,c(targ.features[1:f])])
    
    predictions_prob = predict(rf_model,test_set_matrix[,features], type = 'prob')
    prediction_table = data.frame(GRP_Id= rownames(test_set_matrix), predictions = predictions, reported = test_set_matrix$group, methylation_score = predictions_prob[,colnames(predictions_prob) != 'Control'], model = model)
    prediction_table = prediction_table[order(-prediction_table$methylation_score),]
    prediction_table$auroc = auc_calc(prediction_table,c('Control','Cancer'))
    prediction_table$model = 'rf'
    #train_performance = getTrainPerf(rf_model)
    #prediction_table$TrainROC = train_performance$TrainROC
    # prediction_table$TrainSens = train_performance$TrainSens
    # prediction_table$TrainSpec = train_performance$TrainSpec
    
    return(prediction_table)
  }
  prediction.setup.rfsurv= function(rf.surv,targ.matrix.test,model = 'RFSurv') {
    
    surv.test.set = targ.matrix.test
    surv.test.set$event =ifelse(surv.test.set$group == 'Control',0,1)
    surv.test.set$censorship_time = test.set[match(surv.test.set$GRP_Id,test.set$GRP_Id),'censorship_time']
    
    o.pred <- predict(object = rf.surv, surv.test.set)
    
    predictions = ifelse(o.pred$predicted >median(o.pred$predicted ),'Cancer','Control')
    prediction_table = data.frame(GRP_Id=surv.test.set$GRP_Id, 
                                  'predictions' = predictions, reported = targ.matrix.test$group, 
                                  methylation_score = o.pred$predicted)
    prediction_table = prediction_table[order(-prediction_table$methylation_score),]
    prediction_table$auroc = auc_calc(prediction_table,c('Control','Cancer'))
    prediction_table$model = 'RFSurv'
    #train_performance = getTrainPerf(lasso.fit)
    #prediction_table$TrainROC = train_performance$TrainROC
    #prediction_table$TrainSens = train_performance$TrainSens
    #prediction_table$TrainSpec = train_performance$TrainSpec
    
    
    return(prediction_table)
  }
  
  results.df.all = NULL
  feature.weights.all = NULL
  require(doMC)
  library(doParallel)
  cl <- makeCluster(ncores)
  registerDoParallel(cl)
  for (fl in feats){
    start= Sys.time()
    results.df = NULL
    f = min(fl, length(targ.features))
    
    #all features as features
    targ.features1 = targ.features[1:f]
    targ.matrix = targ.matrix1[,targ.features[1:f]]
    
    targ.matrix$GRP_Id = rownames(targ.matrix)
    if (standardize == T){
      targ.matrix[,targ.features1] = apply(targ.matrix[,targ.features1],2,scale)
      rownames(targ.matrix) = targ.matrix$GRP_Id
    }
    targ.matrix=  merge(targ.matrix,merged.df,by='GRP_Id')
    
    
    targ.matrix.train = targ.matrix[targ.matrix$GRP_Id %in% train.set$GRP_Id,]
    targ.matrix.test = targ.matrix[targ.matrix$GRP_Id %in% test.set$GRP_Id,]
    
    rownames(targ.matrix.train)= targ.matrix.train$GRP_Id
    rownames(targ.matrix.test) = targ.matrix.test$GRP_Id
    
    targ.matrix.train$group = factor(ifelse(targ.matrix.train$group == 'Control','Control','Cancer'),levels = c('Control','Cancer'))
    targ.matrix.test$group = factor(ifelse(targ.matrix.test$group == 'Control','Control','Cancer'),levels=  c('Control','Cancer'))
    
    
    targ.matrix.train.model = model.matrix(~. -1,targ.matrix.train[,c(targ.features[1:f])] )
    targ.matrix.test.model = model.matrix(~. -1,targ.matrix.test[,c(targ.features[1:f])] )
    
    
    
    #logreg old
    
    prediction.setup.lasso.prev= function(lambda.min,lasso.fit,test_set_matrix,features,model = 'logreg') {
      test_set_matrix.model= model.matrix(~. -1,targ.matrix.test[,features] )
      lasso.predict.prob = predict(lasso.fit, s=lambda.min, newx=test_set_matrix.model, type="response") # check help(predict.glmnet)
      predictions = ifelse(lasso.predict.prob > 0.5,'Cancer','Control')
      prediction_table = data.frame(GRP_Id=test_set_matrix$GRP_Id, 'predictions' = predictions[,1], reported = test_set_matrix$group, methylation_score = lasso.predict.prob[,1])
      prediction_table = prediction_table[order(-prediction_table$methylation_score),]
      prediction_table$auroc = auc_calc(prediction_table,c('Control','Cancer'))
      prediction_table$model = 'logreg.old'
      #train_performance = getTrainPerf(lasso.fit)
      #prediction_table$TrainROC = train_performance$TrainROC
      #prediction_table$TrainSens = train_performance$TrainSens
      #prediction_table$TrainSpec = train_performance$TrainSpec
      
      
      return(prediction_table)
    }
    #glmnet
    ##
    MSEs = NULL
    for (i in 1:10){
      lasso.cv = cv.glmnet(targ.matrix.train.model, targ.matrix.train$group, family="binomial", type.measure="class",  nfolds = 10)
      MSEs <- cbind(MSEs, lasso.cv$cvm)
    }
    rownames(MSEs) <- lasso.cv$lambda
    lambda.min <- as.numeric(names(which.min(rowMeans(MSEs))))
    
    lasso.fit = glmnet(as.matrix(targ.matrix.train.model), 
                       as.numeric(targ.matrix.train$group), 
                       family="binomial",
                       nfolds=10,
                       lambda= lambda.min,
                       type.measure = "class")
    
    #lasso.fit = cv.glmnet(targ.matrix.train.model, targ.matrix.train$group, family="binomial")
    lambda.min = lasso.fit$lambda.min
    prediction_table.logreg.old =prediction.setup.lasso.prev(lasso.fit$lambda.min,lasso.fit,targ.matrix.test,c(targ.features[1:f]),model = 'logreg')
    
    prediction_table.logreg.old$model = c('logreg.old')
    results.df = rbind(results.df,prediction_table.logreg.old)
    #feature.coef =coef(lasso.fit,optimal_lambda)
    feature.coef =coef(lasso.fit,lambda.min)
    
    
    lasso.old.feature.weights =data.frame(PC =rownames(feature.coef), coef= feature.coef[1:length(feature.coef)])
    lasso.old.feature.weights$model = 'logreg.old'
    
    
    #annotating scores
    results.df$fold = fold
    results.df$seed = 75
    results.df$feature = feature.type
    results.df$test = stat.test
    results.df$n.features = fl
    results.df.all = rbind(results.df.all,results.df)
    
    #combining feature weights
    feature.weights = rbind(lasso.old.feature.weights)
    feature.weights$fold = fold
    feature.weights$seed = 75
    feature.weights$feature = feature.type
    feature.weights$test = stat.test
    feature.weights$n.features = fl
    feature.weights.all = rbind(feature.weights.all,feature.weights)
    end = Sys.time()
    
    print(start-end)
  }
  
  return.list = list(predictions = results.df.all,
                     weights = feature.weights.all)
  
  return(return.list)
}


predictive.models = function(targ.matrix1,
                             merged.df,
                             train.set,
                             test.set,
                             targ.features,
                             feats,
                             feature.type,
                             stat.test,no_cores=5, standardize = F) {
  library(caret)
  library(glmnet)
  library(ROCR)
  auc_calc = function(prediction_table,labels = c('Control','Cancer')) {
    tmp = prediction_table
    tmp = tmp[order(-tmp$methylation_score),]
    #labels = c(tmp$reported)
    pred = prediction(predictions = c(tmp$methylation_score) ,labels =  tmp$reported, labels)
    perf_AUC=performance(pred,"auc") #Calculate the AUC value
    AUC=perf_AUC@y.values[[1]]
    return(AUC)
  }
  
  prediction.setup.lasso= function(lambda.min,lasso.fit,test_set_matrix,features,model = 'logreg') {
    targ.matrix.test.model = model.matrix(~. -1,targ.matrix.test[,features] )
    
    lasso.predict.prob = predict(lasso.fit, s='lambda.min', newx=targ.matrix.test.model, type="response") # check help(predict.glmnet)
    predictions = ifelse(lasso.predict.prob > 0.5,'Cancer','Control')
    prediction_table = data.frame(GRP_Id=test_set_matrix$GRP_Id, 'predictions' = predictions[,1], reported = test_set_matrix$group, methylation_score = lasso.predict.prob[,1])
    prediction_table = prediction_table[order(-prediction_table$methylation_score),]
    prediction_table$auroc = auc_calc(prediction_table,c('Control','Cancer'))
    prediction_table$model = 'logreg'
    #train_performance = getTrainPerf(lasso.fit)
    #prediction_table$TrainROC = train_performance$TrainROC
    #prediction_table$TrainSens = train_performance$TrainSens
    #prediction_table$TrainSpec = train_performance$TrainSpec
    
    
    return(prediction_table)
  }
  prediction.setup.rf= function(rf_model,test_set_matrix,features,model = 'RF') {
    
    predictions = predict(rf_model, test_set_matrix)
    #targ.matrix.test.model = model.matrix( ~. -1,test_set_matrix[,c(targ.features[1:f])])
    
    predictions_prob = predict(rf_model,test_set_matrix[,features], type = 'prob')
    prediction_table = data.frame(GRP_Id= rownames(test_set_matrix), predictions = predictions, reported = test_set_matrix$group, methylation_score = predictions_prob[,colnames(predictions_prob) != 'Control'], model = model)
    prediction_table = prediction_table[order(-prediction_table$methylation_score),]
    prediction_table$auroc = auc_calc(prediction_table,c('Control','Cancer'))
    prediction_table$model = 'rf'
    #train_performance = getTrainPerf(rf_model)
    #prediction_table$TrainROC = train_performance$TrainROC
    # prediction_table$TrainSens = train_performance$TrainSens
    # prediction_table$TrainSpec = train_performance$TrainSpec
    
    return(prediction_table)
  }
  prediction.setup.rfsurv= function(rf.surv,targ.matrix.test,model = 'RFSurv') {
    
    surv.test.set = targ.matrix.test
    surv.test.set$event =ifelse(surv.test.set$group == 'Control',0,1)
    surv.test.set$censorship_time = test.set[match(surv.test.set$GRP_Id,test.set$GRP_Id),'censorship_time']
    
    o.pred <- predict(object = rf.surv, surv.test.set)
    
    predictions = ifelse(o.pred$predicted >median(o.pred$predicted ),'Cancer','Control')
    prediction_table = data.frame(GRP_Id=surv.test.set$GRP_Id, 
                                  'predictions' = predictions, reported = targ.matrix.test$group, 
                                  methylation_score = o.pred$predicted)
    prediction_table = prediction_table[order(-prediction_table$methylation_score),]
    prediction_table$auroc = auc_calc(prediction_table,c('Control','Cancer'))
    prediction_table$model = 'RFSurv'
    #train_performance = getTrainPerf(lasso.fit)
    #prediction_table$TrainROC = train_performance$TrainROC
    #prediction_table$TrainSens = train_performance$TrainSens
    #prediction_table$TrainSpec = train_performance$TrainSpec
    
    
    return(prediction_table)
  }
  
  results.df.all = NULL
  feature.weights.all = NULL
  require(doMC)
  library(doParallel)
  cl <- makeCluster(ncores)
  registerDoParallel(cl)
  for (fl in feats){
    start= Sys.time()
    results.df = NULL
    f = min(fl, length(targ.features))
    
    #all features as features
    targ.features1 = targ.features[1:f]
    targ.matrix = targ.matrix1[,targ.features[1:f]]
    
    targ.matrix$GRP_Id = rownames(targ.matrix)
    if (standardize == T){
      targ.matrix[,targ.features1] = apply(targ.matrix[,targ.features1],2,scale)
      rownames(targ.matrix) = targ.matrix$GRP_Id
    }
    targ.matrix=  merge(targ.matrix,merged.df,by='GRP_Id')
    
    
    targ.matrix.train = targ.matrix[targ.matrix$GRP_Id %in% train.set$GRP_Id,]
    targ.matrix.test = targ.matrix[targ.matrix$GRP_Id %in% test.set$GRP_Id,]
    
    rownames(targ.matrix.train)= targ.matrix.train$GRP_Id
    rownames(targ.matrix.test) = targ.matrix.test$GRP_Id
    
    targ.matrix.train$group = factor(ifelse(targ.matrix.train$group == 'Control','Control','Cancer'),levels = c('Control','Cancer'))
    targ.matrix.test$group = factor(ifelse(targ.matrix.test$group == 'Control','Control','Cancer'),levels=  c('Control','Cancer'))
    
    
    targ.matrix.train.model = model.matrix(~. -1,targ.matrix.train[,c(targ.features[1:f])] )
    targ.matrix.test.model = model.matrix(~. -1,targ.matrix.test[,c(targ.features[1:f])] )
    
    #logreg new
    ##
    alpha_grid <- seq(0, 1, by = 0.1)
    
    # Initialize variables to store the results
    cv_errors <- numeric(length(alpha_grid))
    best_lambdas <- numeric(length(alpha_grid))
    
    # Perform grid search over alpha
    for (i in seq_along(alpha_grid)) {
      alpha <- alpha_grid[i]
      #targ.matrix.train will be your train set with features for each column and sample acroww each row
      #cv_fit <- cv.glmnet(as.matrix(targ.matrix.train[,c(targ.features[1:f])]), targ.matrix.train$group, family="binomial", type.measure="class",  nfolds = 10,parallel=F,alpha=alpha)
      cv_fit <- cv.glmnet(targ.matrix.train.model, targ.matrix.train$group, family="binomial", type.measure="class",  nfolds = 10,parallel=F,alpha=alpha)
      
      cv_errors[i] <- min(cv_fit$cvm)  # cvm is the mean cross-validated error
      best_lambdas[i] <- cv_fit$lambda.min
    }
    
    # Find the alpha with the smallest cross-validated error
    optimal_index <- which.min(cv_errors)
    optimal_alpha <- alpha_grid[optimal_index]
    optimal_lambda <- best_lambdas[optimal_index]
    
    
    # Refit the model using the optimal alpha and lambda
    best_fit <- glmnet(targ.matrix.train.model, y=as.numeric(targ.matrix.train$group), alpha = optimal_alpha, lambda = optimal_lambda, family = "binomial", type.measure = "class")
    #test_set_matrix will be your test set matrix (features per column, samples per row) formated the same way as the train set matrix 
    lasso.predict.prob = predict(best_fit, s='lambda.min', newx=targ.matrix.test.model, type="response") 
    
    prediction_table.logreg.predx.predictors.norelh =prediction.setup.lasso(optimal_lambda,best_fit,targ.matrix.test,targ.features[1:f],model = 'logreg')
    prediction_table.logreg.predx.predictors.norelh$model = c('logreg.new')
    results.df = rbind(results.df,prediction_table.logreg.predx.predictors.norelh)
    feature.coef =coef(best_fit,optimal_lambda)
    lasso.feature.weights =data.frame(PC =rownames(feature.coef), coef= feature.coef[1:length(feature.coef)])
    lasso.feature.weights$model = 'logreg.new'
    
    
    #logreg old
    
    prediction.setup.lasso.prev= function(lambda.min,lasso.fit,test_set_matrix,features,model = 'logreg') {
      test_set_matrix.model= model.matrix(~. -1,targ.matrix.test[,features] )
      lasso.predict.prob = predict(lasso.fit, s=lambda.min, newx=test_set_matrix.model, type="response") # check help(predict.glmnet)
      predictions = ifelse(lasso.predict.prob > 0.5,'Cancer','Control')
      prediction_table = data.frame(GRP_Id=test_set_matrix$GRP_Id, 'predictions' = predictions[,1], reported = test_set_matrix$group, methylation_score = lasso.predict.prob[,1])
      prediction_table = prediction_table[order(-prediction_table$methylation_score),]
      prediction_table$auroc = auc_calc(prediction_table,c('Control','Cancer'))
      prediction_table$model = 'logreg.old'
      #train_performance = getTrainPerf(lasso.fit)
      #prediction_table$TrainROC = train_performance$TrainROC
      #prediction_table$TrainSens = train_performance$TrainSens
      #prediction_table$TrainSpec = train_performance$TrainSpec
      
      
      return(prediction_table)
    }
    #glmnet
    ##
    MSEs = NULL
    for (i in 1:10){
      lasso.cv = cv.glmnet(targ.matrix.train.model, targ.matrix.train$group, family="binomial", type.measure="class",  nfolds = 10)
      MSEs <- cbind(MSEs, lasso.cv$cvm)
    }
    rownames(MSEs) <- lasso.cv$lambda
    lambda.min <- as.numeric(names(which.min(rowMeans(MSEs))))
    
    lasso.fit = glmnet(as.matrix(targ.matrix.train.model), 
                       as.numeric(targ.matrix.train$group), 
                       family="binomial",
                       nfolds=10,
                       lambda= lambda.min,
                       type.measure = "class")
    
    #lasso.fit = cv.glmnet(targ.matrix.train.model, targ.matrix.train$group, family="binomial")
    #lambda.min = lasso.fit$lambda.min
    prediction_table.logreg.old =prediction.setup.lasso.prev(lambda.min,lasso.fit,targ.matrix.test,c(targ.features[1:f]),model = 'logreg')
    
    prediction_table.logreg.old$model = c('logreg.old')
    results.df = rbind(results.df,prediction_table.logreg.old)
    feature.coef =coef(lasso.fit,lambda.min)
    lasso.old.feature.weights =data.frame(PC =rownames(feature.coef), coef= feature.coef[1:length(feature.coef)])
    lasso.old.feature.weights$model = 'logreg.old'
    
    #
    prediction.setup.lasso.cv= function(lasso.fit,lambda.min,test_set_matrix,features,model = 'logreg') {
      test_set_matrix.model= model.matrix(~. -1,targ.matrix.test[,features] )
      
      lasso.predict.prob = predict(lasso.fit, s=lambda.min, newx=test_set_matrix.model, type="response") # check help(predict.glmnet)
      predictions = ifelse(lasso.predict.prob > 0.5,'Cancer','Control')
      prediction_table = data.frame(GRP_Id=test_set_matrix$GRP_Id, 'predictions' = predictions[,1], reported = test_set_matrix$group, methylation_score = lasso.predict.prob[,1])
      prediction_table = prediction_table[order(-prediction_table$methylation_score),]
      prediction_table$auroc = auc_calc(prediction_table,c('Control','Cancer'))
      prediction_table$model = 'logreg'
      #train_performance = getTrainPerf(lasso.fit)
      #prediction_table$TrainROC = train_performance$TrainROC
      #prediction_table$TrainSens = train_performance$TrainSens
      #prediction_table$TrainSpec = train_performance$TrainSpec
      
      
      return(prediction_table)
    }
    
    MSEs = NULL
    for (i in 1:10){
      print(i)
      lasso.cv = cv.glmnet(targ.matrix.train.model, targ.matrix.train$group, family="binomial", type.measure="class",  nfolds = 5,parallel=T)
      MSEs <- cbind(MSEs, lasso.cv$cvm)
    }
    rownames(MSEs) <- lasso.cv$lambda
    lambda.min <- as.numeric(names(which.min(rowMeans(MSEs))))
    
    lasso.fit = cv.glmnet(targ.matrix.train.model, 
                          targ.matrix.train$group, 
                          family="binomial",
                          nfolds=10,
                          parallel=T)
    
    feature.coef =coef(lasso.fit,lambda.min) 
    lasso.feature.weights =data.frame(PC =rownames(feature.coef), coef= feature.coef[1:length(feature.coef)])
    lasso.feature.weights$model = 'logreg'
    # lambda.minUsing cross validation glmnet
    lambda.min = lasso.fit$lambda.min
    prediction_table.logreg.predx.predictors.norelh =prediction.setup.lasso.cv(lasso.fit,lambda.min,targ.matrix.test,c(targ.features[1:f]),model = 'logreg')
    prediction_table.logreg.predx.predictors.norelh$model = c('logreg.cv')
    results.df = rbind(results.df,prediction_table.logreg.predx.predictors.norelh)
    
    
    #rf
    start1=Sys.time()
    
    control = trainControl(method = 'repeatedcv', number = 10, repeats = 10, search = 'random', classProbs = TRUE, summaryFunction = twoClassSummary)
    tunegrid <- expand.grid(.mtry=seq(10,30,10))
    metric = 'Accuracy'
    
    ntrees = c(1000)
    ntree_list = list()
    acc = c()
    
    targ.matrix.train$group =factor(targ.matrix.train$group, levels = c('Control','Cancer'))
    targ.matrix.train.model = model.matrix(group ~. -1,targ.matrix.train[,c('group',targ.features[1:f])] )
    targ.matrix.train.model$group = targ.matrix.train$group
    targ.matrix.test.model = model.matrix( ~. -1,targ.matrix.test[,c(targ.features[1:f])])
    #targ.matrix.test.model$group = targ.matrix.test.model$group
    
    #targ.matrix.train$group = factor(targ.matrix.train$group, levels = c('Control','Cancer'))
    rf_model = train(group ~ . -1, data = targ.matrix.train[,c('group',targ.features[1:f])] , method = 'rf', metric = metric, tuneGrid = tunegrid, trControl = control, ntrees = ntrees)
    
    
    rf.performance = prediction.setup.rf(rf_model,targ.matrix.test,targ.features[1:f],model='rf')
    rf.performance$model ='RF'
    rf.feature.weights=varImp(rf_model)$importance
    rf.feature.weights=data.frame(PC = rownames(rf.feature.weights), 
                                  coef = rf.feature.weights[,1],
                                  model = 'RF')
    
    results.df = rbind(results.df,rf.performance)
    end1=Sys.time()
    
    
    
    #annotating scores
    results.df$fold = fold
    results.df$seed = seedno
    results.df$feature = feature.type
    results.df$test = stat.test
    results.df$n.features = fl
    results.df.all = rbind(results.df.all,results.df)
    
    #combining feature weights
    feature.weights = rbind(lasso.feature.weights,lasso.old.feature.weights, rf.feature.weights)
    feature.weights$fold = fold
    feature.weights$seed = seedno
    feature.weights$feature = feature.type
    feature.weights$test = stat.test
    feature.weights$n.features = fl
    feature.weights.all = rbind(feature.weights.all,feature.weights)
    end = Sys.time()
    
    print(start-end)
  }
  
  return.list = list(predictions = results.df.all,
                     weights = feature.weights.all)
  
  return(return.list)
}




####reading in annotations####

cpg_count = readRDS('/.mounts/labs/awadallalab/private/ncheng/references/cpg_sites/cpg_site_positions/window_count/hg38_cpg_window_300_count.RDS') #number of cpg sites across 300bp regions
background = as.character(cpg_count$window)

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

colnames(hallmark.genes) = c('Reference','Hallmark','gene')


compile.annotation.groups = function(x) {
  
  genomic.region = list(CpG_Island = cgi[which(cgi$cpg_region == 'cpg_islands'),'window'],
                        CpG_Shore = cgi[which(cgi$cpg_region == 'cpg_shores'),'window'],
                        CpG_Shelf = cgi[which(cgi$cpg_region == 'cpg_shelves'),'window']
  )
  
  repeat.groups = unique(repeat.element$Repeat)
  for (r in repeat.groups) {
    targ.windows = repeat.element[repeat.element$Repeat == r,'window']
    genomic.region[[paste0('Repeat_',r)]] = targ.windows
  }
  
  gene.types = unique(gene.body$gene.type)
  for (g in gene.types) {
    targ.windows = gene.body[gene.body$gene.type == g,'window']
    genomic.region[[paste0('Gene.Type_',g)]] = targ.windows
  }
  
  genomic.region[['Regulatory_Enhancer']] = regulatory[regulatory$Enhancer == 'Enhancer','window']
  genomic.region[['Regulatory_Promoter']] = regulatory[regulatory$Promoter == 'Promoter','window']
  
  hallmarks = unique(hallmark.genes$Hallmark)
  for (h in hallmarks) {
    targ.genes= hallmark.genes[hallmark.genes$Hallmark == h,]
    targ.windows = gene.body[gene.body$genes %in% targ.genes$gene,'window']
    genomic.region[[paste0('Hallmark_',h)]] = targ.windows
  }
  repeat.annotation = repeat.element
  repeat.annotation$Repeat.simp = ifelse(repeat.annotation$Repeat == 'Non-repetitive Region','Non-repetitive Region', 
                                         ifelse(repeat.annotation$Repeat %in% c('Alu','AluJ','AluS','AluY','AmnSINE1','AmnSINE2','MIR','SVA','MER'),'SINE',
                                                ifelse( repeat.annotation$Repeat %in% c('ALR/Alpha'),'Satellite',
                                                        ifelse(repeat.annotation$Repeat %in% c('LTR','MLT','MST','SATR','THE1'),'LTR',
                                                               ifelse(repeat.annotation$Repeat %in% c('L1','L2','L3'),'LINE',
                                                                      ifelse(repeat.annotation$Repeat %in% c('G-rich','A-rich','GA-rich'),'STR', 
                                                                             ifelse(repeat.annotation$Repeat %in% c('Tigger'),'DNA Repeat','Other')))))))
  
  
  
  simp.repeats = unique(repeat.annotation$Repeat.simp)
  for (r in simp.repeats) {
    targ.windows = repeat.annotation[repeat.annotation$Repeat.simp == r,'window']
    genomic.region[[paste0('Repeat.group_',r)]] = targ.windows
  }
  return(genomic.region)
  
  
}

targ.gene.annotation = function(targ.windows, gene.body,hallmark.genes) {
  targ.genes= gene.body[gene.body$window %in% targ.windows,c('window','genes')]
  colnames(targ.genes)[2] = 'gene'
  
  overlap.hallmarks = hallmark.genes[hallmark.genes$gene %in% targ.genes$gene,]
  targ.windows.df = data.frame(gene = targ.windows)
  return.df = merge(targ.genes, overlap.hallmarks[,c('Hallmark','gene')],by='gene')
  return.df = merge(targ.windows.df, return.df,by='window',all.x = T)
  return(return.df)
  
}
autosome_filt = function(medips.count_df) {
  tmp =gsub(':.*','',medips.count_df)
  return_df = medips.count_df[!tmp %in% c('chrY','chrX','chrM')]
  return(return_df)
}

targ.filt = function(x) {
  a = x[x %in% cpg_count$window]
  a = a[!a %in% blacklist$window]
  a = autosome_filt(a)
  #a = a[a %in% feature.filt[['brca']]]
  #a = a[!a %in% feature.filt[['blood']]]
  return(a)
}
background.autosome = targ.filt(background)

targ.filt = function(x) {
  a = x[x %in% cpg_count$window]
  a = a[!a %in% blacklist$window]
  a = a[a %in% regulatory$window | a %in% cgi$window | a %in% repeat.element$window | a %in% gene.body$window]
  a = autosome_filt(a)
  #a = a[a %in% feature.filt[['brca']]]
  #a = a[!a %in% feature.filt[['blood']]]
  return(a)
}
background.annotation = targ.filt(background)


targ.filt = function(x) {
  a = x[x %in% cpg_count$window]
  a = a[!a %in% blacklist$window]
  a = a[a %in% regulatory$window | a %in% cgi$window | a %in% repeat.element$window | a %in% gene.body$window]
  a = a[!a %in% blood.remove1]
  a = autosome_filt(a)
  #a = a[a %in% feature.filt[['brca']]]
  #a = a[!a %in% feature.filt[['blood']]]
  return(a)
}
background.pbl = targ.filt(background)

filters = list(
  'background' = background,
  'background.autosome' = background.autosome,
  'background.annotation' = background.annotation,
  'background.pbl' = background.pbl)
#

###read in files####
library(DESeq2)
library(caret)
library(ROCR)
library(glmnet)
qc.info =readRDS('/.mounts/labs/awadallalab/private/ncheng/cfmedip_data/cptp_samples/participant_data/complete.data/combined.sequencing.qc.RDS')
qx.info = readRDS('/.mounts/labs/awadallalab/private/ncheng/cfmedip_data/cptp_samples/participant_data/complete.data/combined.aix13.clin.qx.data.RDS')

sample.info = merge(qc.info,qx.info,by='GRP_Id')
sample.info$group =ifelse(sample.info$Cancer =='Control','Control','Cancer')

outliers=readRDS(paste0('/.mounts/labs/awadallalab/private/ncheng/cfmedip_data/cptp_samples/hg38.counts//aix13.brca.exclude.RDS'))
male.outliers = readRDS('/.mounts/labs/awadallalab/private/ncheng/cfmedip_data/cptp_samples/fragmentation/aix13.updated1/pancan.matrix/male.outliers.RDS')
fraglength.outliers = readRDS('/.mounts/labs/awadallalab/private/ncheng/cfmedip_data/cptp_samples/fragmentation/aix13.updated1/all.sample.fragmentation/fraglength//fraglength.outliers.RDS')
sample.info.filt.pretime = sample.info
sample.info.filt = sample.info[!sample.info$GRP_Id %in% c(outliers,male.outliers ,fraglength.outliers),]



sample.info.filt = sample.info.filt[which(sample.info.filt$diff_in_days > 0 | sample.info.filt$Diagnosis_Time == 'Control')  ,]
sample.info.filt = sample.info.filt[-which(sample.info.filt$F1_DIS_CANCER_EVER == 1 & sample.info.filt$Cancer %in% c('control','Control')),]
sample.info.filt = sample.info.filt[which(sample.info.filt$enrichment.score.relH > 2.5 & 
                                            sample.info.filt$THALIANA_BETA > 0.95 & 
                                            sample.info.filt$total.fragment.count > 5000000&
                                            sample.info.filt$total > 5000000 & 
                                            sample.info.filt$maxTruCor > 0.6 & 
                                            #sample.info.filt$numberReadsWOCG/(sample.info.filt$numberReadsWOCG+sample.info.filt$numberReadsCG) > 0.2 &
                                            sample.info.filt$long.10k.proportion < 0.001),]
#sample.info.filt = sample.info.filt[!sample.info.filt$seq.date %in% c('22-06-15','22-04-29','22-05-31'),]
sample.info.filt = sample.info.filt[!as.character(sample.info.filt$seq.date) %in% c('22-06-15','22-04-29','22-05-31'),]

sample.info.filt = sample.info.filt[order(sample.info.filt$Diagnosis_Time),]
sample.info.filt = sample.info.filt[!duplicated(sample.info.filt$GRP_Id),]


sample.info.filt = sample.info.filt[which(sample.info.filt$diff_in_days < 365*5 | sample.info.filt$Cancer =='Control'),]

opt.samples =readRDS('/.mounts/labs/awadallalab/private/ncheng/cfmedip_data/cptp_samples/fragmentation/aix13.updated1/methylation.zero.cv10/cv.top.median.methscores.RDS')
opt.samples = readRDS('/.mounts/labs/awadallalab/private/ncheng/cfmedip_data/cptp_samples/fragmentation/brca.discovery.RDS')
opt.samples$GRP_Id = gsub('_combined.*','',opt.samples$GRP_Id)



#113 brca, 100 control
targets = list('Breast' = sample.info.filt[sample.info.filt$Sex =='Female' & sample.info.filt$Cancer %in% c('Breast','Control'),],#[sample.info.filt$GRP_Id %in% opt.samples$GRP_Id,],
               'Prostate'  = sample.info.filt[sample.info.filt$Cancer %in% c('Prostate','Control') & sample.info.filt$Sex == 'Male' & grepl('AIX_',sample.info.filt$GRP_Id) == T ,],
               'Pancreas'= sample.info.filt[sample.info.filt$Cancer %in% c('Pancreas','Control') ,],
               'PanCan' = sample.info.filt[sample.info.filt$Cancer %in% c('Breast','Prostate','Pancreas','Control'),])

discovery.set = readRDS('/.mounts/labs/awadallalab/private/ncheng/manuscripts/early.cancer.risk/discovery.set2.samples.RDS')

targets = targets[1:2]
targets = lapply(targets, function(x) x[x$GRP_Id %in% discovery.set$GRP_Id,] )

####
train_test_partition_cv.singlesex = function(year2_diagnosis_samples, splits = 10, seed = 10) { 
  set.seed(seed)
  year2_diagnosis_samples$diagnosis_time_group =year2_diagnosis_samples$Cancer
  tmp_list = split(year2_diagnosis_samples, year2_diagnosis_samples[,c('Cancer','filler')]) #split into list according to diagnosis time group
  tmp_list =tmp_list[sapply(tmp_list,nrow) > 0]
  tmp_list_train = lapply(tmp_list, function(x) {
    yourData<-x[sample(nrow(x)),]
    #Create 10 equally size folds
    folds <- cut(seq(1,nrow(yourData)),breaks=splits,labels=FALSE)
    yourData$fold = folds
    return(yourData)
  })
  return_df = do.call('rbind',tmp_list_train)
  return(return_df)
}
train_test_partition_cv.validation = function(year2_diagnosis_samples, splits = 10, seed = 10) { 
  set.seed(seed)
  year2_diagnosis_samples$diagnosis_time_group =year2_diagnosis_samples$Cancer
  tmp_list = split(year2_diagnosis_samples, year2_diagnosis_samples[,c('Cancer','filler','Diagnosis_Time','Gleason.score')]) #split into list according to diagnosis time group
  tmp_list =tmp_list[sapply(tmp_list,nrow) > 0]
  tmp_list_train = lapply(tmp_list, function(x) {
    yourData<-x[sample(nrow(x)),]
    #Create 10 equally size folds
    folds <- cut(seq(1,nrow(yourData)),breaks=splits,labels=FALSE)
    yourData$fold = folds
    return(yourData)
  })
  return_df = do.call('rbind',tmp_list_train)
  return(return_df)
}
targ.samples = targets[[1]]
targ.samples =targ.samples


####coverage across regulatory elements####
library(DESeq2)
#matrix/proportions comparison#
wkdir='/.mounts/labs/awadallalab/private/ncheng/cfmedip_data/cptp_samples/fragmentation/aix13.updated1/all.sample.fragmentation/regulatory.counts/'
dmrdir='/.mounts/labs/awadallalab/private/ncheng/cfmedip_data/cptp_samples/fragmentation/aix13.updated1/all.sample.fragmentation/breast.cancer.cv/'

setwd(wkdir)
marker.list = c('genhancer','silencer') #,'utr3','utr5','ctcfbs') #ctcfbs #paste0(wkdir,'aix13.combined.',inserts,'.',marker,'.norm.counts.RDS')
#marker.list='ctcfbs'

#dmr calling w/o age
#fold=fold
for (fold in foldno){
  for (marker in marker.list[1]) {
    print(paste0(fold,'.',marker))
    for (inserts in c(1000)) {
      savedir=paste0(dmrdir,'regulatory.regions.v10.newsplit5/',marker,'.',inserts,'.ntotalTB.all.v2/')
      dir.create(savedir,recursive = T)
      for (i in names(targets)) {
        if (file.exists(paste0(savedir,i,'.',seedno,'.',fold,'.dmr.RDS')) == F){
          
          start = Sys.time()
          #dds = readRDS(paste0(wkdir,'aix13.combined.',inserts,'.',marker,'.dds.RDS'))
          dds = readRDS(paste0(wkdir,'/aix13.all.dv.1000.genhancer.dds.RDS'))
          #dds.matrix = readRDS(paste0(wkdir,'aix13.combined.',inserts,'.',marker,'.norm.counts.RDS'))
          
          colData(dds)$total_reads= colData(dds)$total/1000000
          
          
          #
          targ.samples = targets[[i]]
          targ.samples = targ.samples[targ.samples$GRP_Id %in% colnames(dds),]
          if (i %in% c('Breast','Prostate')) {
            sample.split = train_test_partition_cv.singlesex(targ.samples, splits = 10, seed = seedno)
          } else if (i %in% c('Pancreas')) {
            sample.split = train_test_partition_cv.multisex(targ.samples, splits = 10, seed = seedno)
            
          } else if (i == 'PanCan') {
            sample.split = train_test_partition_cv.multisex(targ.samples, splits = 10, seed = seedno)
            
          }
          
          train.set = sample.split[sample.split$fold != fold,]
          test.set = sample.split[sample.split$fold == fold,]
          merged.df.filt = train.set
          
          windows = rownames(dds)
          windows = windows[!grepl('chrX|chrY',windows)]
          ##breast cancer vs controls##
          merged.df.filt.targ = train.set
          colData(dds)$total.fragment.count.m = colData(dds)$total.fragment.count/1000000
          colData(dds)$filler = factor(colData(dds)$filler,levels= c('MFiller','UFiller'))
          colData(dds)$group = factor(ifelse(colData(dds)$Cancer =='Control','Control','Cancer'),levels = c('Control','Cancer'))
          
          dds.filt = dds[windows,merged.df.filt.targ$GRP_Id]
          if (i %in% c('Breast','Prostate')) {
            mm = model.matrix(~ filler + group, colData(dds.filt)) # + filler:nested.batch +
            
          } else {
            mm = model.matrix(~ total.fragment.count.m +THALIANA_BETA +filler + Sex + group, colData(dds.filt)) # + filler:nested.batch +
            
          }
          
          dds.filt$condition = dds.filt$group
          
          
          ddssva <- DESeq(dds.filt,full = mm,parallel=T) #differential methylation analysis
          res.df = results(ddssva,contrast = list('groupCancer')) #generating results table
          res.df$feature =paste0(marker,'.',inserts)
          res.df$cancer= i
          res.df$window = rownames(res.df)
          res.df$seed = seedno
          res.df$fold = fold
          end = Sys.time()
          saveRDS(res.df,paste0(savedir,i,'.',seedno,'.',fold,'.dmr.RDS'))
          saveRDS(merged.df.filt.targ,paste0(savedir,i,'.',seedno,'.',fold,'.samplesplit.RDS'))
          print(start- end)
          
        }
        
        #end
      }
      
    }}
  
  #dmr calling w age
  
}


#ml modelling
for (fold in foldno){
  for (marker in marker.list) {
    #age adjusted
    for (inserts in c(1000)) {
      dds.matrix = readRDS(paste0(wkdir,'aix13.combined.',inserts,'.',marker,'.norm.counts.RDS'))
      
      #savedir=paste0(dmrdir,'regulatory.regions.v10.ageadjusted/',marker,'.',inserts,'/')
      savedir=paste0(dmrdir,'regulatory.regions.v10.newsplit5/',marker,'.',inserts,'.ntotal/')
      
      
      
      
      for (i in names(targets)[2]) {
        if (file.exists(paste0(savedir,i,'.',seedno,'.',fold,'.dmr.RDS')) == T) {
          
          #
          #start
          # dds.matrix = readRDS('/.mounts/labs/awadallalab/private/ncheng/cfmedip_data/cptp_samples/fragmentation/aix13.updated1/all.sample.fragmentation/AIX13.all.samples.deseq.normcounts.inserts.400.all.300.q20.RDS')
          res.df = data.frame(readRDS(paste0(savedir,i,'.',seedno,'.',fold,'.dmr.RDS')),check.names=F)
          targ.samples = targets[[i]]
          targ.samples = targ.samples[targ.samples$GRP_Id %in% colnames(dds.matrix),]
          targ.samples$group = ifelse(targ.samples$Cancer == 'Control','Control','Cancer')
          train.set = readRDS(paste0(savedir,i,'.',seedno,'.',fold,'.samplesplit.RDS'))
          test.set = targ.samples[!targ.samples$GRP_Id %in% train.set$GRP_Id,]
          
          
          #
          for (dir in c('abs','hyper')){
            predir=paste0(savedir,'predictions.bm1.',dir,'/')
            dir.create(predir, recursive = T)
            
            if(file.exists(paste0(predir,i,'.',seedno,'.',fold,'.predictions.RDS')) == F) {
              
              if (dir == 'hyper'){
                res.df.sig= res.df[which(res.df$log2FoldChange > 0.25 & res.df$pvalue < 0.05 & res.df$baseMean  > 1),]#
                
              } else {
                res.df.sig= res.df[which(abs(res.df$log2FoldChange) > 0.25 & res.df$pvalue < 0.05 & res.df$baseMean  > 1),]#
                
              }
              res.df.sig = res.df.sig[order(res.df.sig$pvalue),]
              #res.df.sig = res.df.sig[rownames(res.df.sig) %in% pbl.regions$window,]
              
              train.set$group = ifelse(train.set$Cancer=='Control','Control','Cancer')
              test.set$group = ifelse(test.set$Cancer=='Control','Control','Cancer')
              
              
              if (nrow(res.df.sig) >1 ){
                
                log2.filt= data.frame(t(log2(dds.matrix[res.df.sig$window,unique(targ.samples$GRP_Id)]+1)),check.names=F) 
                #norm.filt = data.frame(t(dds.matrix[res.df.sig$window,unique(targ.samples$GRP_Id)]),check.names=F)
                
                #log2.std.filt= data.frame(t(apply(log2.filt,1,scale)),check.names=F)
                #rownames(log2.std.filt)= rownames(log2.filt)
                #colnames(log2.std.filt) = colnames(log2.filt)
                
                #raw features
                matrix.list= list(log2.std = log2.filt)
                
                complete.res.base = NULL
                complete.weight.base= NULL
                for (mat in names(matrix.list)) {
                  print(mat)
                  targ.matrix1 = matrix.list[[mat]]
                  #adding filler as feature
                  #a = dds.matrix[,unique(targ.samples$GRP_Id)]
                  targ.features = res.df.sig$window
                  feature.sizes= c(seq(25,200,25),seq(200,500,50))#,1000,min(5000,length(targ.features)))
                  feature.sizes = feature.sizes[feature.sizes <= length(targ.features) ]
                  
                  for (fl in feature.sizes) {
                    res.df.all = predictive.models(targ.matrix1,
                                                   unique(targ.samples[,c('GRP_Id','group','Cancer')]),
                                                   train.set,
                                                   test.set,
                                                   targ.features = c(targ.features[1:fl]),
                                                   feats = feature.sizes,
                                                   feature.type = paste0(marker),
                                                   stat.test = paste0('deseq.',mat))
                    perf.df = res.df.all[[1]] 
                    perf.df$mat = mat
                    perf.df$comparison = i
                    perf.df$direction = dir
                    perf.df$insert.size = inserts
                    
                    perf.df = res.df.all[[2]] 
                    perf.df$mat = mat
                    perf.df$comparison = i
                    perf.df$direction = dir
                    
                    perf.df$insert.size = inserts
                    
                    complete.res.base= rbind(complete.res.base,perf.df)
                    
                    complete.weight.base= rbind(complete.weight.base,perf.df)
                    
                  }
                  
                  
                  
                  
                }
                
                saveRDS(complete.weight.base,paste0(predir,i,'.',seedno,'.',fold,'.feature.weights.RDS'))
                saveRDS(complete.res.base,paste0(predir,i,'.',seedno,'.',fold,'.predictions.RDS'))
                
                
                
              }
              
              
              
              #
              
            }
            
          }
          
          
          
          
          
          
          
          
        }
        #end
      }
      
    }
    
  }
  
  
  
}


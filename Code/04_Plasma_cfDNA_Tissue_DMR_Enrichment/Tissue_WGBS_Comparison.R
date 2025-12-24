options(expressions = 5e5)
library(DESeq2)
library("BiocParallel")
ncores = 10
register(MulticoreParam(ncores))

library(reshape2)

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

#####background region annotation####
cpg_count = readRDS('/.mounts/labs/awadallalab/private/ncheng/references/cpg_sites/cpg_site_positions/window_count/hg38_cpg_window_300_count.RDS') #number of cpg sites across 300bp regions
background = as.character(cpg_count$window)

cpg_count = cpg_count[cpg_count$count >= 4,] #selecting windows with at least 6 or mor CpG sites

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
      lasso.predict.prob = predict(lasso.fit, s=lambda.min, newx=as.matrix(test_set_matrix[,c(features)]), type="response") # check help(predict.glmnet)
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
    feature.coef =coef(best_fit,optimal_lambda)
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



options(expressions = 5e5) 

coxph.calling.all = function(fragment.df.ratio.raw.t,merged.df.filt,filler = T,ncores=10) {
  library(survival)
  combined.df.cox = cbind(merged.df.filt,fragment.df.ratio.raw.t)
  # tmp$censorship_time = ifelse(is.na(tmp$diff_in_days) == F, tmp$diff_in_days,as.Date('01/01/2019',format = '%m/%d/%Y')-as.Date(tmp$collection_date,format='%Y-%m-%d'))
  combined.df.cox = combined.df.cox[combined.df.cox$censorship_time > 0 ,]
  combined.df.cox$event = ifelse(combined.df.cox$group %in% c('control','Control'),0,1)
  windows = colnames(fragment.df.ratio.raw.t)
  female.weights = weightsf.females(combined.df.cox)
  combined.df.cox$TOTAL = combined.df.cox$total/1000000
  return.list= mclapply(windows, function(z) {
    targ.df= combined.df.cox[,c('group','SDC_AGE_CALC','TOTAL','event','censorship_time','THALIANA_BETA','filler',z)] #
    targ.df$methyl = targ.df[,z]
    targ.df$filler = factor(targ.df$filler, levels= c('UFiller','MFiller'))
    if (sum(is.na(targ.df$methyl)) == 0){
      if (filler == T){
        test = summary(coxph(Surv(as.numeric(targ.df$censorship_time), targ.df$event) ~ TOTAL+ THALIANA_BETA+ filler  + methyl,targ.df,weights = female.weights)) #
        
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
###read in files####

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

sample.info.filt.pretimefilt = sample.info.filt

prad.samples.all = sample.info.filt.pretimefilt[sample.info.filt.pretimefilt$Sex == 'Male' & sample.info.filt.pretimefilt$Cancer %in% c('Control','Prostate') & grepl('AIX_',sample.info.filt.pretimefilt$GRP_Id) == T,]
brca.samples.all = sample.info.filt.pretimefilt[sample.info.filt.pretimefilt$Sex == 'Female' & sample.info.filt.pretimefilt$Cancer %in% c('Control','Breast') ,]

combined.all.samples = rbind(prad.samples.all,brca.samples.all )
#saveRDS(combined.all.samples,'/.mounts/labs/awadallalab/private/ncheng/manuscripts/early.cancer.risk/study.sample.info.RDS')

sample.info.filt = sample.info.filt[which(sample.info.filt$diff_in_days < 365*5 | sample.info.filt$Cancer =='Control'),]

brca.opt.samples = readRDS('/.mounts/labs/awadallalab/private/ncheng/cfmedip_data/cptp_samples/fragmentation/brca.discovery.RDS')
brca.opt.samples$GRP_Id = gsub('_combined.*','',brca.opt.samples$GRP_Id)

validation = readRDS('/.mounts/labs/awadallalab/private/ncheng/cfmedip_data/cptp_samples/fragmentation/aix13.updated1/all.sample.fragmentation/qc.filt.validation.samples.updated.RDS')
discovery.prad = validation[validation$fold !=1 & validation$Sex == 'Male',]
discovery.set = sample.info.filt[sample.info.filt$GRP_Id %in% c(discovery.prad$GRP_Id,brca.opt.samples$GRP_Id ),]
discovery.set = discovery.set[discovery.set$Cancer %in% c('Control','Breast','Prostate'),]
#saveRDS(discovery.set,'/.mounts/labs/awadallalab/private/ncheng/manuscripts/early.cancer.risk/discovery.set.samples.RDS' )
validation.set = combined.all.samples[!combined.all.samples$GRP_Id %in% discovery.set$GRP_Id,]
#saveRDS(validation.set,'/.mounts/labs/awadallalab/private/ncheng/manuscripts/early.cancer.risk/validation.set.samples.RDS' )




###identifying tissue specifici regions####
dataset1 = '/.mounts/labs/awadallalab/private/ncheng/vcfs/tissue_methylation/tissue/public.pantissue/GSE233417/bed/'

dataset2 = '/.mounts/labs/awadallalab/private/ncheng/vcfs/tissue_methylation/tissue/public.pantissue/GSE186458/hg38/bed/'

#GSE233417.matrix = readRDS('/.mounts/labs/awadallalab/private/ncheng/vcfs/tissue_methylation/tissue/public.pantissue/GSE233417/bed/rds.files/GSE233417.matrix.RDS')
#GSE186458.matrix = readRDS('/.mounts/labs/awadallalab/private/ncheng/vcfs/tissue_methylation/tissue/public.pantissue/GSE186458/hg38/bed/GSE186458.matrix.RDS')
GSE186458.matrix = apply(GSE186458.matrix,2,function(x) ifelse(x == -1,NA,x))

#visualizing distribution of CpG sites in 10 randomly slected samples

#plotting methylation distribution
GSE186458.matrix$window = rownames(GSE186458.matrix)

GSE186458.matrix.melt = melt(GSE186458.matrix[,c('window',sample(colnames(GSE186458.matrix),30))],'window')
#GSE186458.matrix.melt = melt(GSE186458.matrix,'window')

plot1 = ggplot(GSE186458.matrix.melt, aes(x = value, col = variable))+ 
  geom_density()+
  #scale_color_manual(values = cancer.colors)+#+ ggtitle(title) +
  theme_bw()+
  theme(text = element_text(size=10),
        axis.text=element_text(size=10, face = "bold"),
        axis.title=element_text(size=10,face="bold"),
        legend.position = "none")+ guides(colour = guide_legend(override.aes = list(size=4),nrow=7)) +
  xlab('Methylation Level') + ylab('Proportion of CpG Sites')



savedir="/.mounts/labs/awadallalab/private/ncheng/cfmedip_data/cptp_samples/fragmentation/aix13.updated1/motif.analysis/tissue.specificity/"
dir.create(savedir,recursive = T)
setwd(savedir)
title='GSE186458.10tissue.methyl.distribution.png'
png(title,height = 1000, width = 3000,res=200)
print(plot1)
dev.off()


###looking at sample distribution among brca and prad dmrs####
library(ggplot2)
library(ggh4x)
library(GenomicRanges)
wgbs.cpg.windows = GSE186458.matrix$window
wgbs.cpg.windows.df =data.frame(chr=gsub(':.*','',wgbs.cpg.windows),
                                start= gsub('.*:','',gsub('-.*','',wgbs.cpg.windows)),
                                end = gsub('.*-','',wgbs.cpg.windows),
                                window = wgbs.cpg.windows)

wgbs.grange = makeGRangesFromDataFrame(wgbs.cpg.windows.df,
                                       keep.extra.columns=T,
                                       
                                       seqnames.field=c("chr"),
                                       start.field="start",
                                       end.field=c("end"),
                                       strand.field="strand")

dmrdir='/.mounts/labs/awadallalab/private/ncheng/cfmedip_data/cptp_samples/fragmentation/aix13.updated1/all.sample.fragmentation/breast.cancer.cv/regulatory.regions.v10/validation.breast.test//new.split2/'
res.df = readRDS(paste0(savedir,'discovery.allbg.',sex,'.All.1000.AutoChr.before.dmrs.ntotalTB.RDS'))

#dmrdir='/.mounts/labs/awadallalab/private/ncheng/manuscripts/early.cancer.risk/all.sample.dmr.cpg5/'

female.dmr.df =readRDS(paste0(dmrdir,'discovery.Female.nonage.adjusted.dmrs.inserts.400.all.300.q20.RDS'))
male.dmr.df =readRDS(paste0(dmrdir,'discovery.Male.nonage.adjusted.dmrs.inserts.400.all.300.q20.RDS'))
female.dmr.df$Comparison = c('Incident Breast Cancer vs Control')
male.dmr.df$Comparison = c('Incident Prostate Cancer vs Control')

dmr.list=list(Male = male.dmr.df,
              Female = female.dmr.df)




dmr.bed.overlap = function(wgbs.grange,sig.dmrs.hyper) {
  
  sig.dmrs.windows = sig.dmrs.hyper$window
  sig.dmrs.windows.df = data.frame(chr=gsub(':.*','',sig.dmrs.windows),
                                   start= gsub('.*:','',gsub('-.*','',sig.dmrs.windows)),
                                   end = gsub('.*-','',sig.dmrs.windows),
                                   window = sig.dmrs.windows)
  
  dmr.grange = makeGRangesFromDataFrame(sig.dmrs.windows.df,
                                        keep.extra.columns=T,
                                        
                                        seqnames.field=c("chr"),
                                        start.field="start",
                                        end.field=c("end"),
                                        strand.field="strand")
  
  
  
  overlaps = wgbs.grange[queryHits(findOverlaps(wgbs.grange,dmr.grange )),]
  overlaps.df = data.frame(overlaps)[,c(1:3,6)]
  return(overlaps.df)
  
  
  
}

dmr.overlap.list = list()

#computing CpG sites overlapping with dmrs
for (d in names(dmr.list)) {
  dmr.df=dmr.list[[d]]
  sig.dmrs = dmr.df[!grepl('chrX|chrY',dmr.df$window),]
  sig.dmrs = sig.dmrs[order(sig.dmrs$pvalue),]
  sig.dmrs.hyper = sig.dmrs[sig.dmrs$log2FoldChange > 0.25, ][1:2000,]
  sig.dmrs.hypo = sig.dmrs[sig.dmrs$log2FoldChange < -0.25, ][1:2000,]
  
  sig.dmrs.windows.hyper =sig.dmrs.hyper$window
  sig.dmrs.windows.hypo =sig.dmrs.hyper$window
  
  hypermethylated.cpgs = dmr.bed.overlap(wgbs.grange,sig.dmrs.hyper) 
  hypomethylated.cpgs = dmr.bed.overlap(wgbs.grange,sig.dmrs.hypo) 
  
  dmr.overlap.list[[paste0('Hypermethylated.',d)]] = hypermethylated.cpgs
  dmr.overlap.list[[paste0('Hypomethylated.',d)]] = hypomethylated.cpgs
}

#plotting pheatmap
for (d in names(dmr.overlap.list)) {
  library(pheatmap)
  targ.dmr.df = dmr.overlap.list[[d]]
  targ.tissue.cpgs = GSE186458.matrix[targ.dmr.df$window,]
  colnames(targ.tissue.cpgs) = gsub('.*_','',gsub('-Z0.*','',colnames(targ.tissue.cpgs)))
  col.names = colnames(targ.tissue.cpgs)
  col.names = col.names[!col.names == 'window']
  targ.tissue.cpgs[targ.tissue.cpgs == -1] <- 0
  ph.short =(pheatmap(targ.tissue.cpgs[,col.names],
                      cluster_rows = T,
                      cluster_cols = T,
                      clustering_distance_rows = "euclidean",
                      clustering_distance_cols = "euclidean",
                      border_color=NA,
                      show_rownames = F,
                      show_colnames = T, #
                      fontsize = 4
  ))
  dev.off()
  name = paste0(figdir,d)
  png(paste0(name,'.pheatmap.png'), height = 1500, width = 3000,res=300)
  print(ph.short)
  dev.off()
  
  
}


####labelling tissue specific regions####

#for each region average cpgs within each cell type 
savedir=paste0(dmrdir,'tissue.specificity.wgbs/')
savedir="/.mounts/labs/awadallalab/private/ncheng/cfmedip_data/cptp_samples/fragmentation/aix13.updated1/motif.analysis/tissue.specificity/"

dir.create(savedir,recursive = T)
setwd(savedir)
full.tissue.types = gsub('.*_','',gsub('-Z0.*','',colnames(GSE186458.matrix)))
organ.tissue.types = gsub('-.*','',full.tissue.types)

full.tissue.average = data.frame(window = rownames(GSE186458.matrix))
#GSE186458.matrix = GSE186458.matrix[,!colnames(GSE186458.matrix) == 'window'][,1:10]


for (t in unique(full.tissue.types)) {
  print(t)
  start= Sys.time()
  targ.matrix = GSE186458.matrix[,full.tissue.types == t]
  # targ.matrix = apply(targ.matrix,1,function(x) ifelse())
  if (sum(full.tissue.types == t) > 1) {
    full.tissue.average[,t] = rowMeans(GSE186458.matrix[,full.tissue.types == t],na.rm = T)
    
  } else {
    full.tissue.average[,t] = GSE186458.matrix[,full.tissue.types == t]
    
  }
  end= Sys.time()
  print(start-end)
}
rownames(full.tissue.average) = full.tissue.average$window

saveRDS(full.tissue.average,paste0(savedir,'tissue.celltype.average.RDS'))
#shapiro test for normality
#full.tissue.average = readRDS(paste0(savedir,'tissue.celltype.average.RDS'))

GSE186458.matrix.onehot = apply(GSE186458.matrix,1,function(x) ifelse(x > 0.8,1,0))


saveRDS(GSE186458.matrix.onehot,paste0(savedir,'GSE186458.onehot.tissue.average.matrix.RDS'))

####averaging one hot tissue#####
cell.types=data.frame(ids = colnames(GSE186458.matrix))
cell.types$tissue.cell.type = gsub('.*_','',gsub('-Z.*','',cell.types$ids))
cell.types$tissue = gsub('-.*','',cell.types$tissue.cell.type )
cell.types$cell.type = gsub('.*-','',cell.types$tissue.cell.type )

single.tissue.celltypes = data.frame(table(cell.types$cell.type,cell.types$tissue))
single.tissue.celltypes= single.tissue.celltypes[single.tissue.celltypes$Freq > 0,]
cell.type.freq= data.frame(single.tissue.celltypes$Var1)

#average one hot + methylation level per cell type
tissue.type.average.methyl = data.frame(window = rownames(GSE186458.matrix))
cell.type.average.onehot = data.frame(window=colnames(GSE186458.matrix.onehot))

single.tissue.celltypes = data.frame(table(cell.types$tissue,cell.types$cell.type))
single.tissue.celltypes= single.tissue.celltypes[single.tissue.celltypes$Freq > 0,]
tissue.type.freq= data.frame(single.tissue.celltypes$Var1)


for (ct in unique(tissue.type.freq[,1])) {
  print(ct)
  start= Sys.time()
  targ.cell.type.samples = cell.types[cell.types$tissue == ct,]
  targ.onehot.matrix = GSE186458.matrix.onehot[targ.cell.type.samples$ids,]
  targ.norm.matrix = GSE186458.matrix[,targ.cell.type.samples$ids]
  # targ.matrix = apply(targ.matrix,1,function(x) ifelse())
  if (nrow(targ.cell.type.samples) > 1) {
    cell.type.average.methyl[,ct] = rowMeans(targ.norm.matrix,na.rm = T)
    cell.type.average.onehot[,ct] = colMeans(targ.onehot.matrix,na.rm = T)
    
  } else {
    cell.type.average.methyl[,ct] = targ.norm.matrix
    cell.type.average.onehot[,ct] =targ.onehot.matrix
    
  }
  end= Sys.time()
  print(start-end)
  
}

saveRDS(cell.type.average.methyl,'cell.type.methylation.average.RDS')
saveRDS(cell.type.average.onehot,'cell.type.onehot.average.RDS')

#average one hot + methylation level per tissue  type
#average one hot + methylation level per cell type
tissue.type.average.methyl = data.frame(window = rownames(GSE186458.matrix))
tissue.type.average.onehot = data.frame(window=colnames(GSE186458.matrix.onehot))

single.tissue.celltypes = data.frame(table(cell.types$tissue,cell.types$cell.type))
single.tissue.celltypes= single.tissue.celltypes[single.tissue.celltypes$Freq > 0,]
tissue.type.freq= data.frame(single.tissue.celltypes$Var1)

for (ct in unique(tissue.type.freq[,1])) {
  print(ct)
  start= Sys.time()
  targ.cell.type.samples = cell.types[cell.types$tissue == ct,]
  targ.onehot.matrix = GSE186458.matrix.onehot[targ.cell.type.samples$ids,]
  targ.norm.matrix = GSE186458.matrix[,targ.cell.type.samples$ids]
  # targ.matrix = apply(targ.matrix,1,function(x) ifelse())
  if (nrow(targ.cell.type.samples) > 1) {
    tissue.type.average.methyl[,ct] = rowMeans(targ.norm.matrix,na.rm = T)
    tissue.type.average.onehot[,ct] = colMeans(targ.onehot.matrix,na.rm = T)
    
  } else {
    tissue.type.average.methyl[,ct] = targ.norm.matrix
    tissue.type.average.onehot[,ct] =targ.onehot.matrix
    
  }
  end= Sys.time()
  print(start-end)
  
}

saveRDS(tissue.type.average.methyl,'tissue.type.methylation.average.RDS')
saveRDS(tissue.type.average.onehot,'tissue.type.onehot.average.RDS')


#average one hot + methylation level perc tissue cell type
cell.tissue.type.average.methyl = data.frame(window = rownames(GSE186458.matrix))
cell.tissue.type.average.onehot = data.frame(window=colnames(GSE186458.matrix.onehot))

single.tissue.celltypes = data.frame(table(cell.types$tissue,cell.types$cell.type))
single.tissue.celltypes= single.tissue.celltypes[single.tissue.celltypes$Freq > 0,]
cell.tissue.tissue.type.freq= data.frame(table(cell.types$tissue.cell.type))



for (ct in unique(cell.tissue.tissue.type.freq[,1])) {
  print(ct)
  start= Sys.time()
  targ.cell.type.samples = cell.types[cell.types$tissue.cell.type == ct,]
  targ.onehot.matrix = GSE186458.matrix.onehot[targ.cell.type.samples$ids,]
  targ.norm.matrix = GSE186458.matrix[,targ.cell.type.samples$ids]
  # targ.matrix = apply(targ.matrix,1,function(x) ifelse())
  if (nrow(targ.cell.type.samples) > 1) {
    cell.tissue.type.average.methyl[,ct] = rowMeans(targ.norm.matrix,na.rm = T)
    cell.tissue.type.average.onehot[,ct] = colMeans(targ.onehot.matrix,na.rm = T)
    
  } else {
    cell.tissue.type.average.methyl[,ct] = targ.norm.matrix
    cell.tissue.type.average.onehot[,ct] =targ.onehot.matrix
    
  }
  end= Sys.time()
  print(start-end)
  
}
saveRDS(cell.tissue.type.average.methyl,'cell.tissue.type.methylation.average.RDS')
saveRDS(cell.tissue.type.average.onehot,'cell.tissue.type.onehot.average.RDS')


####total one hot permutation####
library(reshape2)
library(ggplot2)
library(GenomicRanges)
savedir="/.mounts/labs/awadallalab/private/ncheng/cfmedip_data/cptp_samples/fragmentation/aix13.updated1/motif.analysis/tissue.specificity/"
figdir=paste0(savedir,'figures/')
dir.create(figdir,recursive = T)
setwd(savedir)

dmrdir='/.mounts/labs/awadallalab/private/ncheng/manuscripts/early.cancer.risk/all.sample.dmr.cpg5/'

female.dmr.df =readRDS(paste0(dmrdir,'discovery.Female.nonage.adjusted.dmrs.inserts.400.all.300.q20.RDS'))
male.dmr.df =readRDS(paste0(dmrdir,'discovery.Male.nonage.adjusted.dmrs.inserts.400.all.300.q20.RDS'))
female.dmr.df$Comparison = c('Incident Breast Cancer vs Control')
male.dmr.df$Comparison = c('Incident Prostate Cancer vs Control')

dmr.list=list(Male = male.dmr.df,
              Female = female.dmr.df)
dmr.sig.windows =list()
for (d in names(dmr.list)) {
  dmr.df=dmr.list[[d]]
  dmr.df$window = rownames(dmr.df)
  sig.dmrs = dmr.df[!grepl('chrX|chrY',dmr.df$window),]
  sig.dmrs = sig.dmrs[order(sig.dmrs$pvalue),]
  sig.dmrs.hyper = sig.dmrs[sig.dmrs$log2FoldChange > 0.25 , ][1:2000,]#& sig.dmrs$baseMean > 1
  sig.dmrs.hypo = sig.dmrs[sig.dmrs$log2FoldChange < -0.25, ][1:2000,]# & sig.dmrs$baseMean > 1 
  
  sig.dmrs.windows.hyper =sig.dmrs.hyper$window
  sig.dmrs.windows.hypo =sig.dmrs.hyper$window
  
  dmr.sig.windows[[paste0('Hypermethylated.',d)]] = sig.dmrs.hyper
  dmr.sig.windows[[paste0('Hypomethylated.',d)]] = sig.dmrs.hypo
}

cell.type.average.onehot =readRDS('cell.type.onehot.average.RDS')
tissue.type.average.onehot = readRDS('tissue.type.onehot.average.RDS')
cell.tissue.type.average.onehot = readRDS('cell.tissue.type.onehot.average.RDS')

onehot.list = list('cell.type'=cell.type.average.onehot)#,
'tissue.type'=tissue.type.average.onehot,
'cell.tissue.type'=cell.tissue.type.average.onehot)

methyl.list = list('cell.type'=cell.type.average.methyl,
                   'tissue.type'=tissue.type.average.methyl,
                   'cell.tissue.type'=cell.tissue.type.average.methyl)


#tmp = rowSums(GSE186458.matrix[,-1],na.rm=T)
#GSE186458.matrix.filtered = GSE186458.matrix[tmp > 0,]


#distribution of one hot
onehot.list = list('cell.type'='cell.type.onehot.average.RDS',
                   'tissue.type'='tissue.type.onehot.average.RDS',
                   'cell.tissue.type'='cell.tissue.type.onehot.average.RDS')
#performing permutation analysis
for (o in names(onehot.list)[seedno]) {
  
  targ.df = readRDS(onehot.list[[o]])
  max.groups = ncol(targ.df)-1
  upper.limit = 0.75*max.groups
  targ.df = targ.df[which(rowSums(targ.df[,-1]) < upper.limit),]
  
  rownames(targ.df) = targ.df$window
  print(o)
  onehot.nonzero.freq.all = NULL
  enrichment.df.overall.all = NULL
  
  onehot.nonzero.freq = data.frame(window=targ.df[,1])
  background.windows= targ.df$window
  base.grange.df =data.frame(chr=gsub(':.*','',background.windows),
                             start= gsub('.*:','',gsub('-.*','',background.windows)),
                             end = gsub('.*-','',background.windows),
                             window = background.windows)
  base.grange.df$start= as.numeric(base.grange.df$start)
  base.grange.df$end = as.numeric(base.grange.df$end)
  
  start = Sys.time()
  targ.df.melt = melt(targ.df,'window')
  end = Sys.time()
  print(start-end)
  plot1 = ggplot(targ.df.melt, aes(x = value, col = variable))+ 
    geom_density()+
    #scale_color_manual(values = cancer.colors)+#+ ggtitle(title) +
    theme_bw()+
    theme(text = element_text(size=10),
          axis.text=element_text(size=10, face = "bold"),
          axis.title=element_text(size=10,face="bold"),
          legend.position = "none")+ guides(colour = guide_legend(override.aes = list(size=4),nrow=7)) +
    xlab('One Hot Averaged Level') + ylab('Proportion of CpG Sites')
  
  
  title=paste0(figdir,o,'.onehot.mean.distribution.png')
  png(title,height = 1000, width = 3000,res=200)
  print(plot1)
  dev.off()
  
  #looking at frequency > 0 across group type
  onehot.nonzero.freq[,paste0(o,'.nonzero.groups')] =rowSums(targ.df[,-1] > 0)
  onehot.nonzero.freq[,paste0(o,'.total.groups')] = ncol(targ.df[,-1])
  
  #plotting distribution of nonzero groups per region
  nonzero.freq.df= data.frame(table(onehot.nonzero.freq[,paste0(o,'.nonzero.groups')]))
  
  
  plot1 = ggplot(nonzero.freq.df, aes(x = Var1, y =Freq))+ 
    geom_bar(stat='identity')+
    #scale_color_manual(values = cancer.colors)+#+ ggtitle(title) +
    theme_bw()+
    theme(text = element_text(size=10),
          axis.text=element_text(size=10, face = "bold"),
          axis.title=element_text(size=10,face="bold"),
          legend.position = "none")+ guides(colour = guide_legend(override.aes = list(size=4),nrow=7)) +
    xlab('Number of Groups Methylated') + ylab('Frequency')
  
  
  title=paste0(figdir,o,'.onehot.nonzero.distribution.png')
  png(title,height = 1000, width = 3000,res=200)
  print(plot1)
  dev.off()
  
  #0.25 cutoff freq
  onehot.nonzero.freq[,paste0(o,'.0.25min.groups')] =rowSums(targ.df[,-1] > 0.25)
  nonzero.freq.df= data.frame(table(onehot.nonzero.freq[,paste0(o,'.0.25min.groups')]))
  
  
  plot1 = ggplot(nonzero.freq.df, aes(x = Var1, y =Freq))+ 
    geom_bar(stat='identity')+
    #scale_color_manual(values = cancer.colors)+#+ ggtitle(title) +
    theme_bw()+
    theme(text = element_text(size=10),
          axis.text=element_text(size=10, face = "bold"),
          axis.title=element_text(size=10,face="bold"),
          legend.position = "none")+ guides(colour = guide_legend(override.aes = list(size=4),nrow=7)) +
    xlab('Number of Groups Methylated') + ylab('Frequency')
  
  
  title=paste0(figdir,o,'.onehot.0.25min.distribution.png')
  png(title,height = 1000, width = 3000,res=200)
  print(plot1)
  dev.off()
  
  #0.6 cutoff freq
  onehot.nonzero.freq[,paste0(o,'.0.6min.groups')] =rowSums(targ.df[,-1] > 0.6)
  nonzero.freq.df= data.frame(table(onehot.nonzero.freq[,paste0(o,'.0.6min.groups')]))
  
  
  plot1 = ggplot(nonzero.freq.df, aes(x = Var1, y =Freq))+ 
    geom_bar(stat='identity')+
    #scale_color_manual(values = cancer.colors)+#+ ggtitle(title) +
    theme_bw()+
    theme(text = element_text(size=10),
          axis.text=element_text(size=10, face = "bold"),
          axis.title=element_text(size=10,face="bold"),
          legend.position = "none")+ guides(colour = guide_legend(override.aes = list(size=4),nrow=7)) +
    xlab('Number of Groups Methylated') + ylab('Frequency')
  
  
  title=paste0(figdir,o,'.onehot.0.6min.distribution.png')
  png(title,height = 1000, width = 3000,res=200)
  print(plot1)
  dev.off()
  plot.stuff=T
  enrichment.df.overall = NULL
  #onehot.nonzero.freq.all = rbind(onehot.nonzero.freq,onehot.nonzero.freq.all)
  #0.6 cutoff freq
  onehot.nonzero.freq[,paste0(o,'.0.8min.groups')] =rowSums(targ.df[,-1] > 0.8)
  nonzero.freq.df= data.frame(table(onehot.nonzero.freq[,paste0(o,'.0.8min.groups')]))
  
  
  plot1 = ggplot(nonzero.freq.df, aes(x = Var1, y =Freq))+ 
    geom_bar(stat='identity')+
    #scale_color_manual(values = cancer.colors)+#+ ggtitle(title) +
    theme_bw()+
    theme(text = element_text(size=10),
          axis.text=element_text(size=10, face = "bold"),
          axis.title=element_text(size=10,face="bold"),
          legend.position = "none")+ guides(colour = guide_legend(override.aes = list(size=4),nrow=7)) +
    xlab('Number of Groups Methylated') + ylab('Frequency')
  
  
  title=paste0(figdir,o,'.onehot.0.8min.distribution.png')
  png(title,height = 1000, width = 3000,res=200)
  print(plot1)
  dev.off()
  plot.stuff=T
  enrichment.df.overall = NULL
  #onehot.nonzero.freq.all = rbind(onehot.nonzero.freq,onehot.nonzero.freq.all)
  saveRDS(onehot.nonzero.freq,paste0(o,'.onehot.tissue.specific.count.summary.RDS'))
  
  #plotting tissue average
  if (plot.stuff == F) {
    #tissue clustering among regions with different cutoffs
    cutoff.limit = c(1:10)
    targ.cutoffs= c(paste0(o,'.0.6min.groups'),paste0(o,'.nonzero.groups'),paste0(o,'.0.25min.groups',paste0(o,'.0.8min.groups'))) #,
    for (c in cutoff.limit[foldno]) {
      for (t in targ.cutoffs) {
        
        filtered.positions = onehot.nonzero.freq[which(onehot.nonzero.freq[,t] <= c & onehot.nonzero.freq[,t] > 0) ,'window']
        filtered.grange.df= base.grange.df[base.grange.df$window %in% filtered.positions,]
        
        filtered.grange = makeGRangesFromDataFrame(filtered.grange.df,
                                                   keep.extra.columns=T,
                                                   
                                                   seqnames.field=c("chr"),
                                                   start.field="start",
                                                   end.field=c("end"),
                                                   strand.field="strand")
        
        
        
        for (d in names(dmr.sig.windows)) {
          library(pheatmap)
          targ.dmr.df = dmr.sig.windows[[d]]
          enriched.dmrs = dmr.bed.overlap(filtered.grange,targ.dmr.df) 
          
          targ.df.plot =targ.df
          rownames(targ.df) = targ.df$window
          targ.tissue.cpgs = targ.df[enriched.dmrs$window,-1]
          colnames(targ.tissue.cpgs) = gsub('.*_','',gsub('-Z0.*','',colnames(targ.tissue.cpgs)))
          col.names = colnames(targ.tissue.cpgs)
          col.names = col.names[!col.names == 'window']
          targ.tissue.cpgs[targ.tissue.cpgs == -1] <- 0
          targ.tissue.cpgs = targ.tissue.cpgs[rowSums(targ.tissue.cpgs) >0 & apply(targ.tissue.cpgs,1,sd) >0,]
          pheatmap.plot.df=targ.tissue.cpgs[,col.names]
          pheatmap.plot.df = pheatmap.plot.df[,colSums(pheatmap.plot.df) >0,]
          ph.short =(pheatmap(pheatmap.plot.df,
                              cluster_rows = T,
                              cluster_cols = T,
                              clustering_distance_rows = "correlation",
                              clustering_distance_cols = "correlation",
                              border_color=NA,
                              show_rownames = F,
                              show_colnames = T, #
                              fontsize = 4
          ))
          dev.off()
          name = paste0(figdir,d,'.',o,'.',c,'.',t)
          png(paste0(name,'.pheatmap.png'), height = 1500, width = 3000,res=300)
          print(ph.short)
          dev.off()
          
          ##enrichment for tissue specific regions among cfdna dmrs
          #randomly selecting background regions
          c.cutoff = ifelse(grepl('25',c) == T, 0.25,ifelse(grepl('0.6',c)==T,0.6,0))
          dmr.window.regions = targ.dmr.df
          background.window.regions = male.dmr.df[!grepl('chrX|chrY',male.dmr.df$window),]
          dmr.cpgs =  dmr.bed.overlap(filtered.grange,dmr.window.regions) 
          start=Sys.time()
          targ.filtered.df = targ.df.plot[ filtered.positions,]#long
          targ.tissues = colnames(targ.filtered.df)
          targ.tissues = targ.tissues[!targ.tissues == 'window']
          targ.filtered.df.onehot = targ.filtered.df[,targ.tissues]
          targ.filtered.df.onehot[targ.filtered.df.onehot > c.cutoff] <- 1
          end = Sys.time()
          
          targ.enrichment.onehot= data.frame(specific.regions=colSums(targ.filtered.df.onehot[unique(dmr.cpgs$window),]))
          targ.enrichment.onehot$group ='observed'
          targ.enrichment.onehot$seed = 0
          targ.enrichment.onehot$tissue =rownames(targ.enrichment.onehot)
          start= Sys.time()
          random.background.list = mclapply(1:5000,function(x) {
            set.seed(x)
            random.enriched.dmrs = dmr.bed.overlap(filtered.grange,background.window.regions[sample(background.window.regions$window,nrow(dmr.window.regions)),]) 
            
            targ.enrichment.onehot.tmp= data.frame(specific.regions = colSums(targ.filtered.df.onehot[random.enriched.dmrs$window,]))
            targ.enrichment.onehot.tmp$group ='expected'
            targ.enrichment.onehot.tmp$seed = x
            targ.enrichment.onehot.tmp$tissue =rownames(targ.enrichment.onehot.tmp)
            
            return(targ.enrichment.onehot.tmp)
            
          },mc.cores=ncores )
          random.background.df= do.call('rbind',random.background.list)
          enrichment.df = rbind(targ.enrichment.onehot,random.background.df)
          enrichment.df$shared.min =c
          enrichment.df$onehot.cutoff=t
          enrichment.df$group.comparison = d
          enrichment.df$onehot.level = o
          enrichment.df.split = split(enrichment.df,enrichment.df$tissue)
          enrichment.df.split= lapply(enrichment.df.split, function(x) {
            summary.df = x
            summary.df$zscore = (summary.df$specific.regions -mean(summary.df$specific.regions))/sd(summary.df$specific.regions)
            
            summary.df$sig.p = 2*pnorm(q=summary.df$zscore, lower.tail=FALSE)
            return(summary.df)
          } )
          enrichment.df.return= do.call('rbind',enrichment.df.split)
          end= Sys.time()
          print(start-end)
          #enrichment for tissue specific markers among dmrs
          enrichment.df.overall = rbind(enrichment.df.overall,enrichment.df.return)
          
        }
        
        
      }
      
      saveRDS(enrichment.df.overall,paste0(o,'.cutofffoldno.onehot.enrichment.summary.RDS'))
      
      
      
      
      
    }
    
  }
  
  
}

#extracting tissue specific regions/annotationg
for (o in names(onehot.list)[-1]) {
  
  targ.df = readRDS(onehot.list[[o]])
  max.groups = ncol(targ.df)-1
  upper.limit = 0.75*max.groups
  targ.df = targ.df[which(rowSums(targ.df[,-1]) < upper.limit),]
  
  rownames(targ.df) = targ.df$window
  print(o)
  onehot.nonzero.freq.all = NULL
  enrichment.df.overall.all = NULL
  
  onehot.nonzero.freq = data.frame(window=targ.df[,1])
  background.windows= targ.df$window
  base.grange.df =data.frame(chr=gsub(':.*','',background.windows),
                             start= gsub('.*:','',gsub('-.*','',background.windows)),
                             end = gsub('.*-','',background.windows),
                             window = background.windows)
  base.grange.df$start= as.numeric(base.grange.df$start)
  base.grange.df$end = as.numeric(base.grange.df$end)
  
  start = Sys.time()
  targ.df.melt = melt(targ.df,'window')
  end = Sys.time()
  print(start-end)
  
  
  #looking at frequency > 0 across group type
  onehot.nonzero.freq[,paste0(o,'.nonzero.groups')] =rowSums(targ.df[,-1] > 0)
  onehot.nonzero.freq[,paste0(o,'.total.groups')] = ncol(targ.df[,-1])
  
  #plotting distribution of nonzero groups per region
  nonzero.freq.df= data.frame(table(onehot.nonzero.freq[,paste0(o,'.nonzero.groups')]))
  
  
  #0.25 cutoff freq
  onehot.nonzero.freq[,paste0(o,'.0.25min.groups')] =rowSums(targ.df[,-1] > 0.25)
  nonzero.freq.df= data.frame(table(onehot.nonzero.freq[,paste0(o,'.0.25min.groups')]))
  
  
  
  #0.6 cutoff freq
  onehot.nonzero.freq[,paste0(o,'.0.6min.groups')] =rowSums(targ.df[,-1] > 0.6)
  nonzero.freq.df= data.frame(table(onehot.nonzero.freq[,paste0(o,'.0.6min.groups')]))
  
  
  #0.8 cutoff freq
  onehot.nonzero.freq[,paste0(o,'.0.8min.groups')] =rowSums(targ.df[,-1] > 0.8)
  nonzero.freq.df= data.frame(table(onehot.nonzero.freq[,paste0(o,'.0.8min.groups')]))
  
  
  plot.stuff=T
  enrichment.df.overall = NULL
  #onehot.nonzero.freq.all = rbind(onehot.nonzero.freq,onehot.nonzero.freq.all)
  #saveRDS(onehot.nonzero.freq,paste0(o,'.onehot.tissue.specific.count.summary.RDS'))
  
  #plotting tissue average
  #tissue clustering among regions with different cutoffs
  targ.cutoffs= c(paste0(o,'.0.6min.groups'),paste0(o,'.nonzero.groups'),paste0(o,'.0.25min.groups'),paste0(o,'.0.8min.groups')) #,
  
  c=10
  filter = onehot.nonzero.freq[which(onehot.nonzero.freq[,targ.cutoffs[1]] <= c | 
                                       onehot.nonzero.freq[,targ.cutoffs[2]] <= c |
                                       onehot.nonzero.freq[,targ.cutoffs[3]] <= c|
                                       onehot.nonzero.freq[,targ.cutoffs[4]] <= c) ,]
  
  filter = filter[which(filter[,targ.cutoffs[1]]> 0 | 
                          filter[,targ.cutoffs[2]] > 0|
                          filter[,targ.cutoffs[3]] > 0|
                          filter[,targ.cutoffs[4]] > 0) ,]     
  
  filtered.list = list()
  tissue.df.summary = data.frame(window= filter$window)
  
  for (t in targ.cutoffs) {
    if (grepl(0.6,t) == T) {
      onehot.threshold = 0.6
    } else if (grepl(0.25,t) == T) {
      onehot.threshold = 0.25
      
    } else if (grepl('zero',t) == T) {
      onehot.threshold = 0
      
    } else if (grepl('0.8',t) == T) {
      onehot.threshold = 0.8
      
    }
    filter.targ = filter#[filter[,t] < c,]
    tissue.df= targ.df[which(targ.df$window %in% filter.targ$window),]
    tissue.df.summary[,paste0(onehot.threshold,'.','onehotcutoff')] = 0
    for (tis in colnames(tissue.df)[-1]) {
      tissue.df.summary[,paste0(onehot.threshold,'.','onehotcutoff')] = ifelse(tissue.df[,tis] > onehot.threshold, paste0(tissue.df.summary[,paste0(onehot.threshold,'.','onehotcutoff')],'_',tis),tissue.df.summary[,paste0(onehot.threshold,'.','onehotcutoff')]) 
    }
    tissue.df.summary[,paste0(onehot.threshold,'.','onehotcutoff')] = gsub('0_','',tissue.df.summary[,paste0(onehot.threshold,'.','onehotcutoff')])
    
    non.zero = tissue.df.summary[tissue.df.summary[,paste0(onehot.threshold,'.','onehotcutoff')] != 0,'window']
    targ.df = readRDS(onehot.list[[o]])
    
    cor_mat <- cor(targ.df[,-1], method = "spearman",use="complete.obs") #non.zero
    library(pheatmap)
    ph.short =(pheatmap(cor_mat,
                        cluster_rows = T,
                        cluster_cols = T,
                        clustering_distance_rows = "correlation",
                        clustering_distance_cols = "correlation",
                        border_color=NA,
                        show_rownames = T,
                        show_colnames = T, #
                        fontsize = 4
    ))
    dev.off()
    name = paste0(figdir,d)
    png(paste0(figdir,o,'.',t,'.pairwisecor.pheatmap.png'), height = 2000, width = 2000,res=300)
    #png(paste0(figdir,'all.background.',o,'.pairwisecor.pheatmap.png'), height = 2000, width = 2000,res=300)
    
    print(ph.short)
    dev.off()
    
    
    
    
    
    
  }
  
  
  saveRDS(tissue.df.summary, paste0(savedir,o,'.onehot.tissue.markers.RDS'))
  
  
  
}



#distribution of averaged methylation level
for (m in names(methyl.list)) {
  targ.df = methyl.list[[m]]
  targ.df.melt = melt(targ.df[sample(1:nrow(targ.df),10000),],'window')
  
  plot1 = ggplot(targ.df.melt, aes(x = value, col = variable))+ 
    geom_density()+
    #scale_color_manual(values = cancer.colors)+#+ ggtitle(title) +
    theme_bw()+
    theme(text = element_text(size=10),
          axis.text=element_text(size=10, face = "bold"),
          axis.title=element_text(size=10,face="bold"),
          legend.position = "none")+ guides(colour = guide_legend(override.aes = list(size=4),nrow=7)) +
    xlab('One Hot Averaged Level') + ylab('Proportion of CpG Sites')
  
  
  
  
  title=paste0(figdir,o,'.onehot.mean.distribution.png')
  png(title,height = 1000, width = 3000,res=200)
  print(plot1)
  dev.off()
}

#average methylation level per cell type

#one hot - compute distrubionted of methylated regions (overall, by organ, by cell type)
GSE186458.matrix.onehot.tmp = GSE186458.matrix.onehot[,1:100000]


#one hot - plot frequency of tissue per region
#one hot select 


GSE186458.matrix.onehot.average = apply(GSE186458.matrix,1,function(x) ifelse(x > 0.8,1,0))


for (t in )
  shapiro.bakground= lapply(1:100,  function(x) shapiro.test(data.frame(t(full.tissue.average[x,-1]))[,1])$p.value)
sample(2:nrow(full.tissue.average),100)


#single sample t test/stat test


tissue.mean = data.frame(mean = rowMeans(full.tissue.average[,-1]))

plot1 = ggplot(GSE186458.matrix.melt, aes(x = value, col = variable))+ 
  geom_density()+
  #scale_color_manual(values = cancer.colors)+#+ ggtitle(title) +
  theme_bw()+
  theme(text = element_text(size=10),
        axis.text=element_text(size=10, face = "bold"),
        axis.title=element_text(size=10,face="bold"),
        legend.position = "none")+ guides(colour = guide_legend(override.aes = list(size=4),nrow=7)) +
  xlab('Methylation Level') + ylab('Proportion of CpG Sites')



savedir="/.mounts/labs/awadallalab/private/ncheng/cfmedip_data/cptp_samples/fragmentation/aix13.updated1/motif.analysis/tissue.specificity/"
dir.create(savedir,recursive = T)
setwd(savedir)
title='GSE186458.10tissue.methyl.distribution.png'
png(title,height = 1000, width = 3000,res=200)
print(plot1)
dev.off()

#one sample t test


full.tissue.average.onehot = full.tissue.average[,-1]

full.tissue.average.onehot[full.tissue.average.onehot == '-1'] <- NA
full.tissue.average.onehot = apply(full.tissue.average.onehot,2, function(x) ifelse(x < 0.8,0,1))

full.tissue.average.adjusted = apply(full.tissue.average.onehot, 2, function(x) {
  tmp = as.numeric(x)*(1/length(x))
  return(tmp/sum(tmp,na.rm=T))
} )




###enrichment for tissue specific markers results####
library(ggplot2)
library(pheatmap)
library(reshape2)
savedir="/.mounts/labs/awadallalab/private/ncheng/cfmedip_data/cptp_samples/fragmentation/aix13.updated1/motif.analysis/tissue.specificity/"
dir.create(savedir,recursive = T)
setwd(savedir)
figdir=paste0(savedir,'figures/onehot.enrichment.summary/')
dir.create(figdir,recursive = T)
#distribution of one hot
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
        
        # facet_grid( biological.process ~ . ,scales = 'free_y', space='free')+
        scale_y_discrete(position = "right")+
        #scale_shape_manual(values=c(23))+
        theme_bw()+
        theme(text = element_text(size=11,face='bold'),
              axis.ticks.y = element_blank(),
              axis.ticks.x = element_blank(),
              legend.position = 'right',
              legend.title = element_blank(),
              strip.background = element_rect(fill="white")#,
              #strip.background = element_blank(),
              #strip.text.x = element_blank()
              #axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)
              
        ) +#"#e7f0fa",#"#c9e2f6",95cbee
        scale_fill_gradientn(colors = c('#FFF6F1', "#B0517D", "#482066"), #F5B0C3
                             limits = c(0,max(-log10(targ.permutations$sig.p))))+
        scale_colour_gradientn(colors = c('#FFF6F1', "#B0517D", "#482066"),
                               limits = c(0,max(-log10(targ.permutations$sig.p))))+
        xlab('') #+ ylab('Biological Pathway or Process') # + scale_fill_manual(values =db.colors )#+ ggtitle(title) +
      
      #print(plot1)
      
      
      png(paste0(figdir,comp,'.',o,'.',cutoff.filter,'.enrichment.tile.png'),height=3000,width=3200,res=400)
      print(plot1)
      dev.off()  
      
    }
    
  }
  #plotting geomtile of cell type significance
  
}


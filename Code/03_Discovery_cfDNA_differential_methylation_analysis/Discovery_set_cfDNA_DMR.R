
options(expressions = 5e5)
library(DESeq2)
library("BiocParallel")
ncores = 10
register(MulticoreParam(ncores))


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


cpg_count = readRDS('/.mounts/labs/awadallalab/private/ncheng/references/cpg_sites/cpg_site_positions/window_count/hg38_cpg_window_300_count.RDS') #number of cpg sites across 300bp regions
cpg_count$window = as.character(cpg_count$window)
background=cpg_count$window
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
saveRDS(combined.all.samples,'/.mounts/labs/awadallalab/private/ncheng/manuscripts/early.cancer.risk/study.sample.info.RDS')

sample.info.filt = sample.info.filt[which(sample.info.filt$diff_in_days < 365*5 | sample.info.filt$Cancer =='Control'),]

brca.opt.samples = readRDS('/.mounts/labs/awadallalab/private/ncheng/cfmedip_data/cptp_samples/fragmentation/brca.discovery.RDS')
brca.opt.samples$GRP_Id = gsub('_combined.*','',brca.opt.samples$GRP_Id)d

validation = readRDS('/.mounts/labs/awadallalab/private/ncheng/cfmedip_data/cptp_samples/fragmentation/aix13.updated1/all.sample.fragmentation/qc.filt.validation.samples.updated.RDS')
discovery.prad = validation[validation$fold !=1 & validation$Sex == 'Male',]
discovery.set = sample.info.filt[sample.info.filt$GRP_Id %in% c(discovery.prad$GRP_Id,brca.opt.samples$GRP_Id ),]
discovery.set = discovery.set[discovery.set$Cancer %in% c('Control','Breast','Prostate'),]
#saveRDS(discovery.set,'/.mounts/labs/awadallalab/private/ncheng/manuscripts/early.cancer.risk/discovery.set.samples.RDS' )
validation.set = combined.all.samples[!combined.all.samples$GRP_Id %in% discovery.set$GRP_Id,]
#saveRDS(validation.set,'/.mounts/labs/awadallalab/private/ncheng/manuscripts/early.cancer.risk/validation.set.samples.RDS' )
discovery.set = readRDS('/.mounts/labs/awadallalab/private/ncheng/manuscripts/early.cancer.risk/discovery.set3.samples.RDS')
validation.set = readRDS('/.mounts/labs/awadallalab/private/ncheng/manuscripts/early.cancer.risk/validation.set3.samples.RDS')



combined.set = rbind(discovery.set,validation.set)

combined.set$Event = ifelse(combined.set$Cancer == 'Control',0,1)

combined.set$Alcohol.Frequency = ifelse(combined.set$ALC_CUR_FREQ == '-7' | is.na(combined.set$ALC_CUR_FREQ) == T,'Never',
                                        ifelse(combined.set$ALC_CUR_FREQ == '0','Former',
                                               ifelse(combined.set$ALC_CUR_FREQ == '1','< 1 per month',
                                                      ifelse(combined.set$ALC_CUR_FREQ == '2','~1 per month',
                                                             ifelse(combined.set$ALC_CUR_FREQ == '3','2-3 per month',
                                                                    ifelse(combined.set$ALC_CUR_FREQ == '4','1 per week',
                                                                           ifelse(combined.set$ALC_CUR_FREQ == '5','2-3 per week',
                                                                                  ifelse(combined.set$ALC_CUR_FREQ == '6','4-5 per week',
                                                                                         ifelse(combined.set$ALC_CUR_FREQ == '7','6-7 per week','Not Reported')))))))))
combined.set$Alcohol.Frequency = factor(combined.set$Alcohol.Frequency, levels = c('Never','Former','< 1 per month','~1 per month','2-3 per month','1 per week','2-3 per week','4-5 per week','6-7 per week'))

combined.set$Alch.con.group = ifelse(combined.set$Alcohol.Frequency %in% c('Never'),'Never',
                                     ifelse(combined.set$Alcohol.Frequency %in% c('Former'),'Former',
                                            ifelse(combined.set$Alcohol.Frequency %in% c('< 1 per month','~1 per month','2-3 per month'),'Infrequent',
                                                   ifelse(combined.set$Alcohol.Frequency %in% c('1 per week','2-3 per week'),'Moderate',
                                                          ifelse(combined.set$Alcohol.Frequency %in% c('4-5 per week','6-7 per week'),'Frequent','Other')))))
combined.set$Alch.con.group  = factor(combined.set$Alch.con.group , levels = c('Never','Former','Infrequent','Moderate','Frequent'))

combined.set$Smoking.Frequency = ifelse(combined.set$SMK_CIG_STATUS == 0, 'Never',
                                        ifelse(combined.set$SMK_CIG_STATUS == 1, 'Former', 
                                               ifelse(combined.set$SMK_CIG_STATUS == 2, 'Occasional',
                                                      ifelse(combined.set$SMK_CIG_STATUS == 3, 'Daily','Other'))))
combined.set$Smoking.Frequency = factor(as.character(combined.set$Smoking.Frequency), levels=c('Never','Former','Occasional','Daily'))
combined.set$Sex = ifelse(combined.set$SDC_GENDER == 1,'Male','Female')
combined.set$Event= ifelse(combined.set$Cancer == 'Control',0,1)
combined.set$censorship_time = abs(combined.set$censorship_time)
combined.set = combined.set[!is.na(combined.set$Alch.con.group),]
combined.set = combined.set[combined.set$censorship_time < 365*12,]
combined.set$Family.history.prostate =ifelse(combined.set$DIS_CANCER_F_PROSTATE == 1 | 
                                               combined.set$DIS_CANCER_SIB_PROSTATE == 1 |
                                               combined.set$DIS_CANCER_CHILD_PROSTATE == 1, 'Family History','No Reported Family History')
combined.set[is.na(combined.set$Family.history.prostate),'Family.history.prostate'] = 'No Reported Family History'
combined.set$Family.history.prostate  = factor(as.character(combined.set$Family.history.prostate),levels = c('No Reported Family History','Family History'))
#combined.set = merge(combined.set,ohs.qx.pulled[,c('GRP_Id',colnames(ohs.qx.pulled)[grepl('.*BMI.*',colnames(ohs.qx.pulled))])] )
combined.set$age_group = ifelse(combined.set$SDC_AGE_CALC >= 30 & combined.set$SDC_AGE_CALC < 40, '30-40',
                                ifelse(combined.set$SDC_AGE_CALC >= 40 & combined.set$SDC_AGE_CALC < 50, '40-50',
                                       ifelse(combined.set$SDC_AGE_CALC >= 50 & combined.set$SDC_AGE_CALC < 60, '50-60',
                                              ifelse(combined.set$SDC_AGE_CALC >= 60 & combined.set$SDC_AGE_CALC < 70, '60-70',
                                                     ifelse(combined.set$SDC_AGE_CALC >= 70 & combined.set$SDC_AGE_CALC < 80, '70-80','80+')))))
combined.set$age_group =factor(as.character(combined.set$age_group), levels = c('30-40','40-50','50-60','60-70','70-80','80+'))
combined.set$data.partition = ifelse(combined.set$GRP_Id %in% discovery.set$GRP_Id,'Discovery','Validation')
#saveRDS(combined.set,'/.mounts/labs/awadallalab/private/ncheng/manuscripts/early.cancer.risk/combined.set.samples.RDS')
#####setting savedir####
library(ggplot2)
library(ggh4x)
savedir='/.mounts/labs/awadallalab/private/ncheng/manuscripts/early.cancer.risk/all.sample.dmr.cpg5/'
savedir='/.mounts/labs/awadallalab/private/ncheng/manuscripts/early.cancer.risk/manuscript.dmrs/'

figdir=paste0(savedir,'figures/')
dir.create(figdir,recursive = T)
###merginging male and female####
matrix.dir='/.mounts/labs/awadallalab/private/ncheng/cfmedip_data/cptp_samples/fragmentation/aix13.updated1/all.sample.fragmentation/'

discovery.set = readRDS('/.mounts/labs/awadallalab/private/ncheng/manuscripts/early.cancer.risk/discovery.set.samples.RDS')
cpg_count.all = readRDS('/.mounts/labs/awadallalab/private/ncheng/references/cpg_sites/cpg_site_positions/window_count/hg38_cpg_window_300_count.RDS') #number of cpg sites across 300bp regions
cpg_count = cpg_count.all[cpg_count.all$count >= 5,] #selecting windows with at least 6 or mor CpG sites
blacklist = readRDS('/.mounts/labs/awadallalab/private/ncheng/references/encode.blacklist/hg38-blacklist.v2_window300.RDS')

f = 'inserts.1000.all.300.q20.RDS'
files.df.male = readRDS(paste0(matrix.dir,'AIX13.','Male','.rawcounts.',f))
files.df.female = readRDS(paste0(matrix.dir,'AIX13.','Female','.rawcounts.',f))

files.df.male = files.df.male[,colnames(files.df.male) %in% discovery.set$GRP_Id,]
files.df.female = files.df.female[,colnames(files.df.female) %in% discovery.set$GRP_Id,]

files.df.combined = merge(files.df.male, files.df.female,by='window',all=T)
files.df.combined[is.na(files.df.combined)] = 0
files.df.combined = files.df.combined[files.df.combined$window %in% cpg_count$window,]
saveRDS(files.df.combined,paste0(matrix.dir,'AIX13.','All','.rawcounts.',f))
#
files.df.combined = readRDS(paste0(matrix.dir,'AIX13.','All','.rawcounts.',f))
combined.sample.info.filt = discovery.set[discovery.set$GRP_Id %in% colnames(files.df.combined),]

background.windows = files.df.combined$window
background.windows = background.windows[!background.windows %in% blacklist$window]

rownames(files.df.combined)= files.df.combined$window
combined.counts1 = files.df.combined[background.windows,combined.sample.info.filt$GRP_Id]

dds <- DESeqDataSetFromMatrix(countData = files.df.combined[,discovery.set$GRP_Id],
                              colData = combined.sample.info.filt,
                              design= ~  group ) #can add gender here but we're only looking at female samples currently
colData(dds)$total.fragment.count.m = colData(dds)$total.fragment.count/1000000
colData(dds)$filler = factor(colData(dds)$filler,levels= c('MFiller','UFiller'))
colData(dds)$group = factor(ifelse(colData(dds)$Cancer =='Control','Control','Cancer'),levels = c('Control','Cancer'))

dds <- estimateSizeFactors(dds) 

saveRDS(dds,paste0(matrix.dir,'AIX13.','All','.dds.',f))
dds.matrix =counts(dds,normalize = T)
saveRDS(dds.matrix,paste0(matrix.dir,'AIX13.','All','.norm.matrix.',f))


#dds = dds[which(rowMeans(combined.counts[,combined.sample.info.filt$GRP_Id]) != 1),]
#nonzero.regions = rownames(dds)
####discovery set only#####
library(DESeq2)
library(caret)
library(ROCR)
library(glmnet)
discovery.set = readRDS('/.mounts/labs/awadallalab/private/ncheng/manuscripts/early.cancer.risk/discovery.set3.samples.RDS')

#dmr calling across 300bp regions#
cpg_count.all = readRDS('/.mounts/labs/awadallalab/private/ncheng/references/cpg_sites/cpg_site_positions/window_count/hg38_cpg_window_300_count.RDS') #number of cpg sites across 300bp regions
cpg_count = cpg_count.all[cpg_count.all$count >= 5,] #selecting windows with at least 6 or mor CpG sites
blacklist = readRDS('/.mounts/labs/awadallalab/private/ncheng/references/encode.blacklist/hg38-blacklist.v2_window300.RDS')


galpdir='/.mounts/labs/awadallalab/private/ncheng/cfmedip_data/cptp_samples/fragmentation/aix13.updated1/all.sample.fragmentation/'
setwd(galpdir)

f='inserts.1000.all.300.q20.RDS'
files = list.files(pattern = paste0('AIX.*',f))
files.combined = files[grepl('combined',files)]
files.nontopup = files[!gsub('_Ct.*|_combined.*','',files) %in% gsub('_Ct.*|_combined.*','',files.combined) ]
files = c(files.combined,files.nontopup)
file.types=c('inserts.400.all.300.q20.RDS','inserts.1000.all.300.q20.RDS')[2]
sample.list = list('Breast' = discovery.set[discovery.set$Sex == 'Female',],
                   'Prostate'  = discovery.set[discovery.set$Sex == 'Male',])
sex.list = sample.list
names(sex.list) = c('Female','Male')
savedir='/.mounts/labs/awadallalab/private/ncheng/manuscripts/early.cancer.risk/manuscript.dmrs/'
dir.create(savedir,recursive = T)

for (sex in names(sex.list)[2]){
  for (f in file.types) {
    if (file.exists(paste0(savedir,'discovery.',sex,'.nonage.adjusted.dmrs.',f)) == F) {
      print(paste0(sex,'.',f))
      matrix.dir='/.mounts/labs/awadallalab/private/ncheng/cfmedip_data/cptp_samples/fragmentation/aix13.updated1/all.sample.fragmentation/'
      files.df = readRDS(paste0(matrix.dir,'AIX13.',sex,'.rawcounts.',f))
      targ.samples =sex.list[[sex]]
      #checking for mclapply that failed 
      setwd(savedir)
      combined.filt = files.df[,-1]
      if (nrow(combined.filt) > 1){
        library(DESeq2)
        
        combined.counts1 = combined.filt[rowSums(combined.filt) > 0,]
        background.windows = rownames(combined.counts1)
        background.windows = background.windows[!background.windows %in% blacklist$window]
        background.windows = background.windows[background.windows %in% cpg_count$window]
        combined.sample.info.filt = targ.samples[targ.samples$GRP_Id %in% colnames(combined.counts1),]
        combined.counts1 = combined.counts1[background.windows,combined.sample.info.filt$GRP_Id]
        
        dds <- DESeqDataSetFromMatrix(countData = combined.counts1,
                                      colData = combined.sample.info.filt,
                                      design= ~  group ) #can add gender here but we're only looking at female samples currently
        colData(dds)$total.fragment.count.m = colData(dds)$total.fragment.count/1000000
        colData(dds)$filler = factor(colData(dds)$filler,levels= c('MFiller','UFiller'))
        colData(dds)$group = factor(ifelse(colData(dds)$Cancer =='Control','Control','Cancer'),levels = c('Control','Cancer'))
        
        #dds = dds[which(rowMeans(combined.counts[,combined.sample.info.filt$GRP_Id]) != 1),]
        #nonzero.regions = rownames(dds)
        dds <- estimateSizeFactors(dds) #estimating size factors for tmm normalization
        #normalized counts
        dds.matrix =counts(dds,normalize =T)
        saveRDS(dds.matrix,paste0(savedir,'discovery.',sex,'.samples.deseq.normcounts.',f))
        saveRDS(dds,paste0(savedir,'discovery.',sex,'.samples.dds.',f))
        
        #age adjusted
        mm = model.matrix(~ filler + group, colData(dds)) # + filler:nested.batch +
        
        dds$condition = dds$group
        
        
        ddssva <- DESeq(dds,full = mm,parallel=T) #differential methylation analysis
        res.df = results(ddssva,contrast = list('groupCancer')) #generating results table
        res.df$feature =f
        res.df$cancer= sex
        res.df$window = rownames(res.df)
        res.df = data.frame(res.df, stringsAsFactors = F)
        
        saveRDS(res.df,paste0(savedir,'discovery.',sex,'.nonage.adjusted.dmrs.',f))
        
        
      }
      
    }
    
    
  }
  
}

###coxph###
names(sex.list) = c('Female','Male')

for (sex in names(sex.list)){
  for (f in file.types) {
    if (file.exists(paste0(savedir,'discovery.',sex,'.nonage.adjusted.dmrs.',f)) == F) {
      print(paste0(sex,'.',f))
      matrix.dir='/.mounts/labs/awadallalab/private/ncheng/cfmedip_data/cptp_samples/fragmentation/aix13.updated1/all.sample.fragmentation/'
      files.df = readRDS(paste0(matrix.dir,'AIX13.',sex,'.rawcounts.',f))
      targ.samples =sex.list[[sex]]
      #checking for mclapply that failed 
      savedir='/.mounts/labs/awadallalab/private/ncheng/manuscripts/early.cancer.risk/all.sample.dmr/'
      dir.create(savedir)
      setwd(savedir)
      combined.filt = files.df[,-1]
      if (nrow(combined.filt) > 1){
        library(DESeq2)
        
        combined.counts1 = combined.filt[rowSums(combined.filt) > 0,]
        background.windows = rownames(combined.counts1)
        background.windows = background.windows[!background.windows %in% blacklist$window]
        background.windows = background.windows[background.windows %in% cpg_count$window]
        combined.sample.info.filt = targ.samples[targ.samples$GRP_Id %in% colnames(combined.counts1),]
        combined.counts1 = combined.counts1[background.windows,combined.sample.info.filt$GRP_Id]
        a= coxph.calling.all(combined.counts1,combined.sample.info.filt,filler = T,ncores=10)
        saveRDS(a,paste0(savedir,'discovery.',sex,'.samples.deseq.normcounts.',f))
        
        #saveRDS(dds.matrix,paste0(savedir,'discovery.',sex,'.samples.deseq.normcounts.',f))
        #saveRDS(dds,paste0(savedir,'discovery.',sex,'.samples.dds.',f))
        
        
      }
      
    }
    
    
  }
  
}
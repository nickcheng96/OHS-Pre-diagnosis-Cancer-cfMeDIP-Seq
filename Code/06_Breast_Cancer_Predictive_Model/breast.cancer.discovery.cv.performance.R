####loading packages#####
#!/usr/bin/env Rscript
#libraries to load (you will need to install DESeq2 and sva)
library(plyr)
library(ROCR)
library(edgeR)
library(rstatix)
library(DESeq2)
library(sva)
library(matrixTests)
library(car)
library(gridExtra)
library(caret)
library(ROCR)
library(gridExtra)
library(caret)
library(ROCR)
library(pROC)
library(stringr)
library(parallel)
library(survival)
library(ncvreg)
library(survminer)
library(RColorBrewer)
auc_calc = function(prediction_table,labels = c('Control','Cancer')) {
  tmp = prediction_table
  tmp = tmp[order(-tmp$methylation_score),]
  #labels = c(tmp$reported)
  pred = prediction(predictions = c(tmp$methylation_score) ,labels =  tmp$reported, labels)
  perf_AUC=performance(pred,"auc") #Calculate the AUC value
  AUC=perf_AUC@y.values[[1]]
  return(AUC)
}
tpr.fpr.calc = function(x){
  tmp1 = x
  tmp1$f = as.integer(ifelse(tmp1$reported == 'Control', 0, 1))
  tmp1$f.str = tmp1$reported
  
  tmp1 = tmp1[order(-tmp1$methylation_score),]
  case_no = nrow(tmp1[tmp1$reported!='Control',])
  control_no = nrow(tmp1[tmp1$reported=='Control',])
  auc.df =data.frame(matrix(nrow = 0, ncol=3))
  for(l in 1:nrow(tmp1)) {
    x = tmp1[1:l,]
    case_cum = nrow(x[x$reported!='Control',])
    control_cum = nrow(x[x$reported=='Control',])
    tpr = case_cum/case_no
    fpr = control_cum/control_no
    return.tmp = data.frame(l,tpr,fpr)
    auc.df = rbind(auc.df, return.tmp)
  }
  base = data.frame(l = 0, tpr = 0, fpr = 0)
  auc.df = rbind(base,auc.df)
  return(auc.df)
  #auc.df = auc.df[auc.df$fpr > 0,]
}

tpr.fpr.calc = function(x){
  tmp1 = x
  tmp1$f = as.integer(ifelse(tmp1$reported == 'Control', 0, 1))
  tmp1$f.str = tmp1$reported
  
  tmp1 = tmp1[order(-tmp1$methylation_score),]
  case_no = nrow(tmp1[tmp1$reported!='Control',])
  control_no = nrow(tmp1[tmp1$reported=='Control',])
  auc.df =data.frame(matrix(nrow = 0, ncol=3))
  for(l in 1:nrow(tmp1)) {
    x = tmp1[1:l,]
    case_cum = nrow(x[x$reported!='Control',])
    control_cum = nrow(x[x$reported=='Control',])
    tpr = case_cum/case_no
    fpr = control_cum/control_no
    return.tmp = data.frame(l,tpr,fpr)
    auc.df = rbind(auc.df, return.tmp)
  }
  base = data.frame(l = 0, tpr = 0, fpr = 0)
  auc.df = rbind(base,auc.df)
  return(auc.df)
  #auc.df = auc.df[auc.df$fpr > 0,]
}

group.similarity = function(targ.pcs, sample.information,groups = c('group','filler')) {
  library(cluster)
  library(kBET)
  library(stats)
  library(scPOP)
  #LISI
  sample.information.filt = sample.information[sample.information$GRP_Id %in% rownames(targ.pcs),]
  
  lisi.index.tmp = lisi(targ.pcs[sample.information.filt$GRP_Id,],sample.information.filt,groups)
  lisi.index.tmp = apply(lisi.index.tmp,2, mean)
  
  #Silhouette width
  library(plyr)
  
  sl.calc = function(targ.pcs, grouping) {
    sw.tmp.group.tmp = silhouette(as.numeric(as.factor(grouping)),dist(targ.pcs[sample.information.filt$GRP_Id,]))
    a = ddply(data.frame(sw.tmp.group.tmp[,-2]),'cluster', numcolwise(mean))
    return(mean(a[,2]))
  }
  
  k.means.fun = function(targ.pcs, sample.information,l) {
    kmeans.clusters = kmeans(targ.pcs, length(unique(sample.information[,l])), iter.max = 10, nstart = 2,
                             algorithm = c("Hartigan-Wong", "Lloyd", "Forgy",
                                           "MacQueen"), trace=FALSE)
    kmeans.clusters = data.frame(cluster= kmeans.clusters$cluster)
    kmeans.clusters$GRP_Id =rownames(kmeans.clusters)
    kmeans.clusters = merge(kmeans.clusters, sample.information,by='GRP_Id')
    return.df = ari(kmeans.clusters$cluster, kmeans.clusters[,l])
    return(return.df)
  }
  
  sl.calc.res= NULL
  kbet.res = NULL
  ari.res = NULL
  for (g in groups){
    sw.tmp.group = sl.calc(targ.pcs, sample.information.filt[,g])
    sl.calc.res[g] = sw.tmp.group
    kbet.group = kBET(targ.pcs, sample.information.filt[,g])$summary$kBET.observed[1]
    kbet.res[g] = kbet.group
    ari.tmp.group = k.means.fun(targ.pcs,sample.information.filt ,g)
    ari.res[g] = ari.tmp.group
  }
  
  
  #ARI
  return.df=NULL
  for (g in groups ) {
    metrics.results.tmp = data.frame(variable = g, 
                                     LISI = lisi.index.tmp[g],
                                     ASW = sl.calc.res[g],
                                     kBET = kbet.res[g],
                                     ARI = ari.res[g])
    return.df = rbind(return.df,metrics.results.tmp)
  }
  
  return(return.df)
  
  
  
  
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

coxph.calling.all = function(fragment.df.ratio.raw.t,merged.df.filt,filler = T) {
  tmp = cbind(merged.df.filt,fragment.df.ratio.raw.t)
  # tmp$censorship_time = ifelse(is.na(tmp$diff_in_days) == F, tmp$diff_in_days,as.Date('01/01/2019',format = '%m/%d/%Y')-as.Date(tmp$collection_date,format='%Y-%m-%d'))
  tmp = tmp[tmp$censorship_time > 0 ,]
  tmp$event = ifelse(tmp$group %in% c('control','Control'),0,1)
  windows = colnames(fragment.df.ratio.raw.t)
  female.weights = weightsf.females(tmp)
  tmp$TOTAL = tmp$total/1000000
  return.list= mclapply(windows, function(i) {
    targ.df= tmp[,c('group','SDC_AGE_CALC','TOTAL','event','censorship_time','sequence.run','filler',i)] #
    targ.df$methyl = targ.df[,i]
    if (sum(is.na(targ.df$methyl)) == 0){
      if (filler == T){
        test = summary(coxph(Surv(as.numeric(targ.df$censorship_time), targ.df$event) ~ TOTAL+ filler  + methyl,targ.df,weights = female.weights)) #
        
      } else {
        test = summary(coxph(Surv(as.numeric(targ.df$censorship_time), targ.df$event) ~ TOTAL  + methyl,targ.df,weights = female.weights)) #
        
      }
      return.df = data.frame(window =i, pvalue = test$coefficients[nrow(test$coefficients),6], HR =  test$coefficients[nrow(test$coefficients),2],HR.SE = test$coefficients[nrow(test$coefficients),3])
      return(return.df)
    }
    
    
  },mc.cores = detectCores()/2)
  res.df=do.call('rbind',return.list)
  return(res.df)
}
wcx.calling.all = function(fragment.df.ratio.raw.t,merged.df.filt) {
  tmp = cbind(merged.df.filt,fragment.df.ratio.raw.t)
  tmp = tmp[tmp$censorship_time > 0,]
  tmp$event = ifelse(tmp$group %in% c('control','Control'),0,1)
  windows = colnames(fragment.df.ratio.raw.t)
  female.weights = weightsf.females(tmp)
  tmp$TOTAL = tmp$total/1000000
  
  wcx.test = mclapply(windows, function(x) {
    if (is.numeric(tmp[,x]) == T & var(tmp[,x]) > 0) {
      targ.df = tmp[,c('GRP_Id',x,'group')]
      colnames(targ.df)[2] = 'targ'
      res <- wilcox.test(targ ~ group, data = targ.df,
                         exact = FALSE,conf.int=TRUE)
      return.df = data.frame(window = x, pvalue = res$p.value,est = res$estimate)
      return(return.df)
    }
    
  },mc.cores = detectCores()/2)
  wcx.res = do.call('rbind',wcx.test)
  
  return(wcx.res)
}


pca.plot.fillerlabel =function(deseq.matrix,predx.dmrs.sig.hyper,combined.info,dx.time,name,cv.label = T){
  library(ggfortify)
  library(RColorBrewer)
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
  
  if (cv.label == T & length(predx.dmrs.sig.hyper) >= 3){
    tmp.info = combined.info[combined.info$GRP_Id  %in% colnames(deseq.matrix),]
    tmp.info$Diagnosis_Time = diagnosis_time_grouping(tmp.info[,dx.time])
    predx.dmrs.sig.hypo.filt = predx.dmrs.sig.hyper[predx.dmrs.sig.hyper %in% rownames(deseq.matrix)]
    pca.obj = prcomp(data.frame(t(deseq.matrix[predx.dmrs.sig.hypo.filt,tmp.info$GRP_Id]),check.names=F))
    return.df = data.frame(pca.obj$x,check.names=F)
    return.df = cbind(tmp.info,return.df)
    pc1 = paste0('PC1 (',100*round(summary(pca.obj)$importance[2,1],digits = 3),'%)')
    pc2 = paste0('PC2 (',100*round(summary(pca.obj)$importance[2,2],digits = 3),'%)')
    pc3 = paste0('PC3 (',100*round(summary(pca.obj)$importance[2,3],digits = 3),'%)')
    
    diagnosis_time_colors1 = c('#7A797C',"#048BA8",'#AAF683','#FFD97D','#FF9B85','#C8553D')
    names(diagnosis_time_colors1) = c('Control','0-2','2-4','4-6','6-8','8-10')
    plot1 = ggplot(return.df,aes(col = Diagnosis_Time,shape = filler, x= PC1, y= PC2))+ geom_point(size = 3,alpha = 0.8) + 
      scale_color_manual(values = diagnosis_time_colors1)+ #+ ggtitle(title) +
      theme_bw()+
      scale_shape_manual(values = c("Methylated" = 19,"Unmethylated" = 17))+
      theme(text = element_text(size=12),
            axis.text=element_text(size=12, face = "bold"),
            axis.title=element_text(size=14,face="bold"),
            legend.position = "none")+ guides(colour = guide_legend(override.aes = list(size=4),nrow=7)) + xlab(pc1) + ylab(pc2)+ ggtitle(paste0( ' Features: ',length(predx.dmrs.sig.hyper)))
    
    plot2 = ggplot(return.df,aes(col = Diagnosis_Time,shape = filler, x= PC1, y= PC3))+ geom_point(size = 3,alpha = 0.8) + 
      scale_color_manual(values = diagnosis_time_colors1)+#+ ggtitle(title) +
      theme_bw()+ 
      scale_shape_manual(values = c("Methylated" = 19,"Unmethylated" = 17))+
      theme(text = element_text(size=12),
            axis.text=element_text(size=12, face = "bold"),
            axis.title=element_text(size=14,face="bold"),
            legend.position = "none")+ guides(colour = guide_legend(override.aes = list(size=4),nrow=7)) + xlab(pc1) + ylab(pc3)+ ggtitle(paste0( ' Features: ',length(predx.dmrs.sig.hyper)))
    
    
    plot3=ggplot(return.df,aes(col = Diagnosis_Time,shape = filler, x= PC2, y= PC3))+ geom_point(size = 3,alpha = 0.8) + 
      scale_color_manual(values = diagnosis_time_colors1)+#+ ggtitle(title) +
      theme_bw()+
      scale_shape_manual(values = c("Methylated" = 19,"Unmethylated" = 17))+
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
    plot4 = autoplot(pca.obj, data = tmp.info, col = 'Diagnosis_Time',shape = 'filler',size = 3,alpha = 0.8,x=2,y=3) + 
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
    
    
    n <- length(unique(tmp.info$sequence.run))
    qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
    col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
    
    plot1 = ggplot(return.df,aes(shape = filler,col = sequence.run, x= PC1, y= PC2))+ geom_point(size = 3,alpha = 0.8) + 
      scale_color_manual(values = col_vector)+#      theme_bw()+ stat_ellipse(geom = "polygon", aes(fill = after_scale(alpha(colour, 0))))+
      scale_shape_manual(values = c("Methylated" = 19,"Unmethylated" = 17))+
      theme(text = element_text(size=12),
            axis.text=element_text(size=12, face = "bold"),
            axis.title=element_text(size=14,face="bold"),
            legend.position = "none")+ guides(colour = guide_legend(override.aes = list(size=4),nrow=7)) + xlab(pc1) + ylab(pc2)+ ggtitle(paste0( ' Features: ',length(predx.dmrs.sig.hyper)))
    
    plot2 =ggplot(return.df,aes(shape = filler,col = sequence.run, x= PC1, y= PC3))+ geom_point(size = 3,alpha = 0.8) + 
      scale_color_manual(values = col_vector)+#  
      theme_bw()+ stat_ellipse(geom = "polygon", aes(fill = after_scale(alpha(colour, 0))))+
      scale_shape_manual(values = c("Methylated" = 19,"Unmethylated" = 17))+
      theme(text = element_text(size=12),
            axis.text=element_text(size=12, face = "bold"),
            axis.title=element_text(size=14,face="bold"),
            legend.position = "none")+ guides(colour = guide_legend(override.aes = list(size=4),nrow=7)) + xlab(pc1) + ylab(pc3)+ ggtitle(paste0( ' Features: ',length(predx.dmrs.sig.hyper)))
    
    
    plot3 =ggplot(return.df,aes(shape = filler,col = sequence.run, x= PC2, y= PC3))+ geom_point(size = 3,alpha = 0.8) + 
      scale_color_manual(values = col_vector)+#  
      theme_bw()+ stat_ellipse(geom = "polygon", aes(fill = after_scale(alpha(colour, 0))))+
      scale_shape_manual(values = c("Methylated" = 19,"Unmethylated" = 17))+
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
    
    plot3 = autoplot(pca.obj, data = tmp.info, col = 'sequence.run',shape = 'filler',size = 3,alpha = 0.8,x=2,y=3) + 
      scale_color_manual(values = col_vector)+#+ ggtitle(title) +
      theme_bw()+
      scale_shape_manual(values = c("Methylated" = 19,"Unmethylated" = 17))+
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
    
    plot1 = ggplot(return.df,aes(col = filler,shape = filler, x= PC1, y= PC2))+ geom_point(size = 3,alpha = 0.8) + 
      scale_color_manual(values = col_vector)+#  
      theme_bw()+ stat_ellipse(geom = "polygon", aes(fill = after_scale(alpha(colour, 0))))+
      scale_shape_manual(values = c("Methylated" = 19,"Unmethylated" = 17))+
      theme(text = element_text(size=12),
            axis.text=element_text(size=12, face = "bold"),
            axis.title=element_text(size=14,face="bold"),
            legend.position = "none")+ guides(colour = guide_legend(override.aes = list(size=4),nrow=7)) + xlab(pc1) + ylab(pc2)+ ggtitle(paste0( ' Features: ',length(predx.dmrs.sig.hyper)))
    
    plot2 =ggplot(return.df,aes(col = filler,shape = filler, x= PC1, y= PC3))+ geom_point(size = 3,alpha = 0.8) + 
      scale_color_manual(values = col_vector)+#        theme_bw()+ stat_ellipse(geom = "polygon", aes(fill = after_scale(alpha(colour, 0))))+
      scale_shape_manual(values = c("Methylated" = 19,"Unmethylated" = 17))+
      theme(text = element_text(size=12),
            axis.text=element_text(size=12, face = "bold"),
            axis.title=element_text(size=14,face="bold"),
            legend.position = "none")+ guides(colour = guide_legend(override.aes = list(size=4),nrow=7)) + xlab(pc1) + ylab(pc3)+ ggtitle(paste0( ' Features: ',length(predx.dmrs.sig.hyper)))
    
    
    plot3 = ggplot(return.df,aes(col = filler,shape = filler, x= PC2, y= PC3))+ geom_point(size = 3,alpha = 0.8) + 
      scale_color_manual(values = col_vector)+#        theme_bw()+ stat_ellipse(geom = "polygon", aes(fill = after_scale(alpha(colour, 0))))+
      scale_shape_manual(values = c("Methylated" = 19,"Unmethylated" = 17))+
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
    
    plot3 = autoplot(pca.obj, data = tmp.info, col = 'filler',shape = 'filler',size = 3,alpha = 0.8,x=2,y=3) + 
      scale_color_manual(values = col_vector)+#+ ggtitle(title) +
      theme_bw()+
      scale_shape_manual(values = c("Methylated" = 19,"Unmethylated" = 17))+
      theme(text = element_text(size=12),
            axis.text=element_text(size=12, face = "bold"),
            axis.title=element_text(size=14,face="bold"),
            legend.position = "right")+ guides(colour = guide_legend(override.aes = list(size=4),nrow=7)) + xlab(pc2) + ylab(pc3)+ ggtitle(paste0( ' Features: ',length(predx.dmrs.sig.hyper)))
    png(paste0(name,'.filler.pca.labs.png'),height = 1000, width = 3000,res=200)
    print(plot3)
    dev.off()
    
  } else {
    tmp.info = combined.info[combined.info$GRP_Id  %in% colnames(deseq.matrix),]
    tmp.info$Diagnosis_Time = diagnosis_time_grouping(tmp.info[,dx.time])
    predx.dmrs.sig.hypo.filt = predx.dmrs.sig.hyper[predx.dmrs.sig.hyper %in% rownames(deseq.matrix)]
    pca.obj = prcomp(data.frame(t(deseq.matrix[predx.dmrs.sig.hypo.filt,tmp.info$GRP_Id]),check.names=F),scale=T)
    return.df = data.frame(pca.obj$x,check.names=F)
    pc1 = paste0('PC1 (',100*round(summary(pca.obj)$importance[2,1],digits = 3),'%)')
    pc2 = paste0('PC2 (',100*round(summary(pca.obj)$importance[2,2],digits = 3),'%)')
    diagnosis_time_colors1 = c('#7A797C',"#048BA8",'#AAF683','#FFD97D','#FF9B85','#C8553D')
    names(diagnosis_time_colors1) = c('Control','0-2','2-4','4-6','6-8','8-10')
    plot1 = ggplot(return.df,aes(col = Diagnosis_Time,shape = filler, x= PC1, y= PC2))+ geom_point(size = 3,alpha = 0.8) + 
      scale_color_manual(values = diagnosis_time_colors1)+#+ ggtitle(title) +
      theme_bw()+ stat_ellipse(geom = "polygon", aes(fill = after_scale(alpha(colour, 0))))+
      scale_shape_manual(values = c("Methylated" = 19,"Unmethylated" = 17))+
      theme(text = element_text(size=12),
            axis.text=element_text(size=12, face = "bold"),
            axis.title=element_text(size=14,face="bold"),
            legend.position = "none")+ guides(colour = guide_legend(override.aes = list(size=4),nrow=7)) + xlab(pc1) + ylab(pc2)+ ggtitle(paste0( ' Features: ',length(predx.dmrs.sig.hyper)))
    
    title=paste0(name,'.dx.time.pca.png')
    # title=paste0(figdir,title)
    # title=paste0(figdir,title)
    png(title,height = 1000, width = 1000,res=200)
    print(plot1)
    dev.off()
    
    
    
    n <- length(unique(tmp.info$sequence.run))
    qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
    col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
    
    plot1 = ggplot(return.df,aes(col = sequence.run,shape = filler, x= PC1, y= PC2))+ geom_point(size = 3,alpha = 0.8) + 
      scale_color_manual(values = col_vector)+#+ ggtitle(title) +
      theme_bw()+
      scale_shape_manual(values = c("Methylated" = 19,"Unmethylated" = 17))+
      theme(text = element_text(size=12),
            axis.text=element_text(size=12, face = "bold"),
            axis.title=element_text(size=14,face="bold"),
            legend.position = "none")+ guides(colour = guide_legend(override.aes = list(size=4),nrow=7)) + xlab(pc1) + ylab(pc2)+ ggtitle(paste0( ' Features: ',length(predx.dmrs.sig.hyper)))
    
    
    
    title=paste0(name,'.batch.pca.png')
    # title=paste0(figdir,title)
    # title=paste0(figdir,title)
    png(title,height = 1000, width = 1000,res=200)
    print(plot1)
    dev.off()
    
    
    
    n <- length(unique(tmp.info$filler))
    qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
    col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
    
    plot1 =  ggplot(return.df,aes(col = filler,shape = filler, x= PC1, y= PC2))+ geom_point(size = 3,alpha = 0.8) + 
      scale_color_manual(values = col_vector)+#+ ggtitle(title) +
      theme_bw()+
      scale_shape_manual(values = c("Methylated" = 19,"Unmethylated" = 17))+
      theme(text = element_text(size=12),
            axis.text=element_text(size=12, face = "bold"),
            axis.title=element_text(size=14,face="bold"),
            legend.position = "none")+ guides(colour = guide_legend(override.aes = list(size=4),nrow=7)) + xlab(pc1) + ylab(pc2)+ ggtitle(paste0( ' Features: ',length(predx.dmrs.sig.hyper)))
    
    
    title=paste0(name,'.filler.pca.png')
    # title=paste0(figdir,title)
    # title=paste0(figdir,title)
    png(title,height = 1000, width = 1000,res=200)
    print(plot1)
    dev.off()
    
    
    
    
  }
  
  return(return.df)
}
pca.plot.fillerlabel.project =function(deseq.matrix,predx.dmrs.sig.hyper,combined.info,dx.time,name,cv.label = T,filler){
  library(ggfortify)
  library(gridExtra)
  library(ggpubr)
  library(ggplot2)
  
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
    predx.dmrs.sig.hypo.filt = predx.dmrs.sig.hyper[predx.dmrs.sig.hyper %in% rownames(deseq.matrix)]
    pca.obj = prcomp(data.frame(t(deseq.matrix[predx.dmrs.sig.hypo.filt,tmp.info[tmp.info$GRP_Id %in% train.set$GRP_Id,'GRP_Id']]),check.names=F))
    
    # project new data onto the PCA space
    pca.obj1 = scale(t(deseq.matrix[predx.dmrs.sig.hypo.filt,tmp.info[tmp.info$GRP_Id %in% test.set$GRP_Id,'GRP_Id']]), pca.obj$center, pca.obj$scale) %*% pca.obj$rotation 
    
    return.df = data.frame(pca.obj$x,check.names=F)
    return.df1 = data.frame(pca.obj1,check.names=F)
    pca.plot = rbind(return.df,return.df1)
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
      scale_shape_manual(values = c("Methylated" = 19,"Unmethylated" = 17))+
      theme(text = element_text(size=12),
            axis.text=element_text(size=12, face = "bold"),
            axis.title=element_text(size=14,face="bold"),
            legend.position = "none")+ guides(colour = guide_legend(override.aes = list(size=4),nrow=7)) + xlab(pc1) + ylab(pc2)+ ggtitle(paste0( ' Features: ',length(predx.dmrs.sig.hyper)))
    
    plot2 = ggplot(pca.plot, aes(x = PC1, y = PC3, col = Diagnosis_Time,shape = filler))+ geom_point(size = 3,alpha = 0.8) + 
      scale_color_manual(values = diagnosis_time_colors1)+#+ ggtitle(title) +
      theme_bw()+
      scale_shape_manual(values = c("Methylated" = 19,"Unmethylated" = 17))+
      theme(text = element_text(size=12),
            axis.text=element_text(size=12, face = "bold"),
            axis.title=element_text(size=14,face="bold"),
            legend.position = "none")+ guides(colour = guide_legend(override.aes = list(size=4),nrow=7)) + xlab(pc1) + ylab(pc3)+ ggtitle(paste0( ' Features: ',length(predx.dmrs.sig.hyper)))
    
    
    plot3 = ggplot(pca.plot, aes(x = PC2, y = PC3, col = Diagnosis_Time,shape = filler))+ geom_point(size = 3,alpha = 0.8) + 
      scale_color_manual(values = diagnosis_time_colors1)+#+ ggtitle(title) +
      theme_bw()+
      scale_shape_manual(values = c("Methylated" = 19,"Unmethylated" = 17))+
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
    
    
    n <- length(unique(tmp.info$sequence.run))
    qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
    col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
    
    plot1 = ggplot(pca.plot, aes(x = PC1, y = PC2, col = sequence.run,shape = filler))+ geom_point(size = 3,alpha = 0.8) + 
      scale_color_manual(values = col_vector)+#+ ggtitle(title) +
      theme_bw()+
      scale_shape_manual(values = c("Methylated" = 19,"Unmethylated" = 17))+
      theme(text = element_text(size=12),
            axis.text=element_text(size=12, face = "bold"),
            axis.title=element_text(size=14,face="bold"),
            legend.position = "none")+ guides(colour = guide_legend(override.aes = list(size=4),nrow=7)) + xlab(pc1) + ylab(pc3)+ ggtitle(paste0( ' Features: ',length(predx.dmrs.sig.hyper)))
    
    plot2 = ggplot(pca.plot, aes(x = PC2, y = PC3, col = sequence.run,shape = filler))+ geom_point(size = 3,alpha = 0.8) + 
      scale_color_manual(values = col_vector)+#+ ggtitle(title) +
      theme_bw()+
      scale_shape_manual(values = c("Methylated" = 19,"Unmethylated" = 17))+
      theme(text = element_text(size=12),
            axis.text=element_text(size=12, face = "bold"),
            axis.title=element_text(size=14,face="bold"),
            legend.position = "none")+ guides(colour = guide_legend(override.aes = list(size=4),nrow=7)) + xlab(pc1) + ylab(pc3)+ ggtitle(paste0( ' Features: ',length(predx.dmrs.sig.hyper)))
    
    
    plot3 = ggplot(pca.plot, aes(x = PC1, y = PC3, col = sequence.run,shape = filler))+ geom_point(size = 3,alpha = 0.8) + 
      scale_color_manual(values = col_vector)+#+ ggtitle(title) +
      theme_bw()+
      scale_shape_manual(values = c("Methylated" = 19,"Unmethylated" = 17))+
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
    
    plot3 = ggplot(pca.plot, aes(x = PC1, y = PC2, col = sequence.run,shape = filler))+ geom_point(size = 3,alpha = 0.8) + 
      scale_color_manual(values = col_vector)+#+ ggtitle(title) +
      theme_bw()+
      scale_shape_manual(values = c("Methylated" = 19,"Unmethylated" = 17))+
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
      scale_shape_manual(values = c("Methylated" = 19,"Unmethylated" = 17))+
      theme(text = element_text(size=12),
            axis.text=element_text(size=12, face = "bold"),
            axis.title=element_text(size=14,face="bold"),
            legend.position = "none")+ guides(colour = guide_legend(override.aes = list(size=4),nrow=7)) + xlab(pc1) + ylab(pc2)+ ggtitle(paste0( ' Features: ',length(predx.dmrs.sig.hyper)))
    
    plot2 = ggplot(pca.plot, aes(x = PC1, y = PC3, col = filler,shape = filler))+ geom_point(size = 3,alpha = 0.8) + 
      scale_color_manual(values = col_vector)+#+ ggtitle(title) +
      theme_bw()+
      scale_shape_manual(values = c("Methylated" = 19,"Unmethylated" = 17))+
      theme(text = element_text(size=12),
            axis.text=element_text(size=12, face = "bold"),
            axis.title=element_text(size=14,face="bold"),
            legend.position = "none")+ guides(colour = guide_legend(override.aes = list(size=4),nrow=7)) + xlab(pc1) + ylab(pc3)+ ggtitle(paste0( ' Features: ',length(predx.dmrs.sig.hyper)))
    
    
    plot3 = ggplot(pca.plot, aes(x = PC2, y = PC3, col = filler,shape = filler))+ geom_point(size = 3,alpha = 0.8) + 
      scale_color_manual(values = col_vector)+#+ ggtitle(title) +
      theme_bw()+
      scale_shape_manual(values = c("Methylated" = 19,"Unmethylated" = 17))+
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
      scale_shape_manual(values = c("Methylated" = 19,"Unmethylated" = 17))+
      theme(text = element_text(size=12),
            axis.text=element_text(size=12, face = "bold"),
            axis.title=element_text(size=14,face="bold"),
            legend.position = "right")+ guides(colour = guide_legend(override.aes = list(size=4),nrow=7)) + xlab(pc2) + ylab(pc3)+ ggtitle(paste0( ' Features: ',length(predx.dmrs.sig.hyper)))
    png(paste0(name,'.filler.pca.labs.png'),height = 1000, width = 3000,res=200)
    print(plot3)
    dev.off()
    
  } 
  
  return(return.df)
}


########sample splitting######

merged.df = readRDS('/.mounts/labs/awadallalab/private/ncheng/cfmedip_data/cptp_samples/AIX3/metadata/AIX13.sample.information.RDS')
merged.df = merged.df[order(merged.df$diff_in_days),]
merged.df = merged.df[!duplicated(merged.df$GRP_Id),]
merged.df[is.na(merged.df$topup.order),'topup.order'] = 0
merged.df$topup.order = ifelse(grepl('combined',merged.df$GRP_Id) == T,0,merged.df$topup.order)

topup.combined = merged.df[grepl('combined',merged.df$GRP_Id) == T,]
merged.df.filt = merged.df[!merged.df$GRP_Id.sample %in% topup.combined$GRP_Id.sample,]
merged.df.filt = rbind(merged.df.filt,topup.combined)


outliers1=readRDS('/.mounts/labs/awadallalab/private/ncheng/cfmedip_data/cptp_samples/hg38.counts/outliers.pca.background.RDS')
#outliers = readRDS('/.mounts/labs/awadallalab/private/ncheng/cfmedip_data/cptp_samples/hg38.counts/outliers.pca.dmrs.RDS')
outliers2=readRDS('/.mounts/labs/awadallalab/private/ncheng/cfmedip_data/cptp_samples/hg38.counts/outliers.conditions.RDS')
outliers=readRDS(paste0('/.mounts/labs/awadallalab/private/ncheng/cfmedip_data/cptp_samples/hg38.counts//aix13.brca.exclude.RDS'))

merged.df.filt = merged.df.filt[which(merged.df.filt$diff_in_days > 0 | merged.df.filt$Diagnosis_Time == 'Control')  ,]
merged.df.filt = merged.df.filt[!merged.df.filt$GRP_Id.sample %in% c(outliers),] #outliers

merged.df.filt = merged.df.filt[-which(merged.df.filt$F1_DIS_CANCER_EVER == 1 & merged.df.filt$group %in% c('control','Control')),]
merged.df.filt = merged.df.filt[which(merged.df.filt$enrichment.score.relH > 2.4 & merged.df.filt$THALIANA_BETA > 0.95 & merged.df.filt$total > 5000000 & merged.df.filt$maxTruCor > 0.6 & merged.df.filt$numberReadsWOCG/(merged.df.filt$numberReadsWOCG+merged.df.filt$numberReadsCG) > 0.15),]

merged.df.all = merged.df.filt[merged.df.filt$Cancer %in% c('Control','Breast'),]
merged.df.all = merged.df.all#[merged.df.all$Sex == 'Female',]
merged.df.all = merged.df.all[!merged.df.all$GRP_Id %in% 'AIX_0030',]
#merged.df.all = merged.df.all[merged.df.all$Cancer == 'Control' | merged.df.all$diff_in_days < 365*5,]

weightsf.males<-function (matched_samples){
  matched_samples_temp<- matched_samples[matched_samples$censorship_time > 0 | matched_samples$group %in% c('control','Control'),]#[na.omit(match(matched_samples$GRP_Id ,ids)),]
  latest.linkage=  as.Date('2019-01-01',format= '%Y-%m-%d')
  matched_samples_temp$gender= 1
  control.incidence = matched_samples_temp[matched_samples_temp$group %in% c('control','Control'),]
  
  expected_rate_breast_cr <- mean(aml_inc_cr(matched_samples_temp$gender, matched_samples_temp$SDC_AGE_CALC, matched_samples_temp$SDC_AGE_CALC+pmax(1, matched_samples_temp$censorship_time/365))[matched_samples_temp$group %in% c('control','Control')])
  
  n_total_prostate <- sum(!matched_samples_temp $group %in% c("Control",'control'))/expected_rate_breast_cr
  n_total_prostate
  weights <- rep(1, nrow(matched_samples_temp))
  weights[matched_samples_temp$group %in% c("Control",'control')] <- n_total_prostate/sum(matched_samples_temp$group %in% c("Control",'control'))
  
  return(weights)
}


predir=paste0(savedir,'/predictions/')
figdir=paste0(savedir,'figures/')
dir.create(predir,recursive = T)
dir.create(figdir,recursive = T)

discovery.set

####results analysis####
library(cutpointr)
library(dplyr)
library(plyr)
library(ggplot2)
library(ggh4x)
library(plyr)
library(ROCR)


#####methylation dmr performance#####
savedir='/.mounts/labs/awadallalab/private/ncheng/cfmedip_data/cptp_samples/fragmentation/aix13.updated1/methylation.insert.select/predictions/'
savedir='/.mounts/labs/awadallalab/private/ncheng/cfmedip_data/cptp_samples/fragmentation/aix13.updated1/methylation.insert.select/predictions.updated/'
savedir='/.mounts/labs/awadallalab/private/ncheng/cfmedip_data/cptp_samples/fragmentation/aix13.updated1/methylation.insert.select1/predictions1/'
savedir='/.mounts/labs/awadallalab/private/ncheng/cfmedip_data/cptp_samples/fragmentation/aix13.updated1/methylation.insert.validation.split/predictions/'
savedir='/.mounts/labs/awadallalab/private/ncheng/cfmedip_data/cptp_samples/fragmentation/aix13.updated1/methylation.insert.validation.split1/predictions/'
savedir='/.mounts/labs/awadallalab/private/ncheng/cfmedip_data/cptp_samples/fragmentation/aix13.updated1/methylation.insert.check/predictions/'
savedir='/.mounts/labs/awadallalab/private/ncheng/cfmedip_data/cptp_samples/fragmentation/aix13.updated1/methylation.insert.validation.split.8020/predictions'
#savedir='//.mounts/labs/awadallalab/private/ncheng/cfmedip_data/cptp_samples/fragmentation/aix13.updated1/methylation.insert.validation.80.male.controls/predictions/'
savedir='/.mounts/labs/awadallalab/private/ncheng/cfmedip_data/cptp_samples/fragmentation/thesis/brca.cv10/preupdate.opt/predictions/'
savedir='/.mounts/labs/awadallalab/private/ncheng/cfmedip_data/cptp_samples/fragmentation/aix13.updated1/methylation.insert.select1/predictions1/'
#savedir='//.mounts/labs/awadallalab/private/ncheng/cfmedip_data/cptp_samples/thesis/figures/chapter3/breast.dmrs.cv.norm.all/predictions.hyper/'
#savedir='/.mounts/labs/awadallalab/private/ncheng/cfmedip_data/cptp_samples/thesis/figures//chapter3/breast.dmrs.cv.norm.brcaonly//predictions.hyper/'
savedir='/.mounts/labs/awadallalab/private/ncheng/cfmedip_data/cptp_samples/fragmentation/aix13.updated1/all.sample.fragmentation/breast.cancer.cv/'
savedir='/.mounts/labs/awadallalab/private/ncheng/cfmedip_data/cptp_samples/brca/regulatory.300/predictions.basefilt.hyper/'


####updated####
library(DESeq2)
#matrix/proportions comparison#
#wkdir='/.mounts/labs/awadallalab/private/ncheng/cfmedip_data/cptp_samples/fragmentation/aix13.updated1/all.sample.fragmentation/regulatory.counts/'
marker.list = c('genhancer','silencer','utr3','utr5','TSS.2','Bin300')
marker.list= c('DNA\ Repeat','Alu','L1','LINE','LTR','genhancer','silencer')
marker.list = paste0(marker.list,'.1000')
#marker.list = c('Bin300.400','Bin300.1000')
marker.list='Bin300.combined.matrix.400'
marker.list = c('Bin300.combined.matrix.400','genhancer.1000','silencer.1000','TSS.2.400')
#marker.list = c('genhancer','silencer','TSS.2','Bin300')#,'ctcfbs') #ctcfbs #paste0(wkdir,'aix13.combined.',inserts,'.',marker,'.norm.counts.RDS')
validation.samples = readRDS('/.mounts/labs/awadallalab/private/ncheng/cfmedip_data/cptp_samples/fragmentation/aix13.updated1/all.sample.fragmentation/qc.filt.validation.samples.updated.RDS')
sample.info = readRDS('/.mounts/labs/awadallalab/private/ncheng/manuscripts/early.cancer.risk/discovery.set3.samples.RDS')

#validation.samples = readRDS('/.mounts/labs/awadallalab/private/ncheng/cfmedip_data/cptp_samples/fragmentation/aix13.updated1/all.sample.fragmentation/qc.filt.validation.samples.updated.thesis.RDS')
#wkdir='/.mounts/labs/awadallalab/private/ncheng/cfmedip_data/cptp_samples/fragmentation/aix13.updated1/all.sample.fragmentation/prostate.cancer.cv/regulatory.regions.v5/'
pno=1
#wkdir='/.mounts/labs/awadallalab/private/ncheng/cfmedip_data/cptp_samples/fragmentation/aix13.updated1/all.sample.fragmentation/breast.cancer.cv/regulatory.regions.v10.ageadjusted/'
#wkdir='/.mounts/labs/awadallalab/private/ncheng/cfmedip_data/cptp_samples/fragmentation/aix13.updated1/all.sample.fragmentation/breast.cancer.cv/regulatory.regions.v10/'
#wkdir='/.mounts/labs/awadallalab/private/ncheng/cfmedip_data/cptp_samples/fragmentation/aix13.updated1/all.sample.fragmentation/breast.cancer.cv/regulatory.regions.v6/'
wkdir='/.mounts/labs/awadallalab/private/ncheng/cfmedip_data/cptp_samples/fragmentation/aix13.updated1/all.sample.fragmentation/breast.cancer.cv/regulatory.regions.v10/'
wkdir='/.mounts/labs/awadallalab/private/ncheng/cfmedip_data/cptp_samples/fragmentation/aix13.updated1/all.sample.fragmentation/breast.cancer.cv/regulatory.regions.v10.newsplit7/'
discovery.set = readRDS('/.mounts/labs/awadallalab/private/ncheng/manuscripts/early.cancer.risk/discovery.set3.samples.RDS')
marker.list = c('silencer.1000','genhancer.1000')[1]
marker.list = c('Bin300.combined.matrix.1000')
marker.list = c('Bin300.combined.matrix.ntotalTB.1000')#,'Bin300.combined.matrix.ntotal.1000')#
marker.list = c('Bin300.combined.matrix.ntotalTB.1000','Bin300.combined.matrix.femaleall.ntotalTB.1000')
marker.list= c('genhancer.1000.ntotalTB.femalenorm.v2','genhancer.1000.ntotal.femalenorm.v2','genhancer.1000.ntotal.allnorm.v2','genhancer.1000.ntotalTB.allnorm.v2')
marker.list = c('silencer.1000.ntotalTB.allnorm.v')
#marker.list = c('silencer.1000.ntotalTB.allnorm.v3','genhancer.1000.ntotalTB.allnorm.v3')#'genhancer.1000.ntotalTB.femalenorm.v3',

marker.list = c('genhancer.1000.ntotalTB.allnorm.v1')#'genhancer.1000.ntotalTB.femalenorm.v3',
#marker.list = c('Bin300.combined.matrix.femaleall.ntotalTB.1000','Bin300.combined.matrix.ntotalTB.1000')
#dir.list = c('base.abs','base.annotfilt1.hyper','base.annotfilt.hyper','base.hyper','bm1.abs','bm1.hyper','hyper')
for (pno in 1){
  
  median.perf.list = list()
  mean.perf.list = list()
  for (m in marker.list){
    print(m)
    # validation.train= validation.samples[validation.samples$fold != pno,]
    # validation.test =  validation.samples[validation.samples$fold == pno,]
    directions = c('covar.bm1.abs','covar.base.abs')#,'abs','hypo')
    directions = c('test.hyper')#,'abs','hypo')
    directions = c('bm1.final.abs')#,'abs','hypo')
    directions='base.annotfilt1.hyper'
    directions = c('base.abs','base.annotfilt1.hyper','base.annotfilt.hyper','base.hyper','bm1.abs','bm1.hyper','hyper')
    directions='base.basefilt1.hyper'
    directions = paste0(rep(c('abs','hyper'),3),rep(c('.None','.filler','.filler.total'),2))
    directions = c('abs.filler','abs.total','abs.filler.total','hyper.filler','hyper.total','hyper.filler.total')
    directions = c('abs.filler','abs.total','abs.filler.total','hyper.filler','hyper.total','hyper.filler.total')
    directions = list.files(pattern='predictionsbm1.*|predictionsft.*',path =paste0(wkdir,m,'/'))
    directions = list.files(pattern='predictionsft.*0.45|predictionsft.*0|predictionsft.*0.1',path =paste0(wkdir,m,'/'))
    #directions = list.files(pattern='predictionsft..*|predictionsftcpg5..*',path =paste0(wkdir,m,'/'))
    directions = list.files(pattern='predictionsu6.all.*|predictionsu6.reg.*',path =paste0(wkdir,m,'/'))
    directions = list.files(pattern='predictions.bm1.*',path =paste0(wkdir,m,'/'))
    
    
    if (length(directions) > 0) {
      pred.df.all = NULL
      for (dir in directions) {
        print(dir)
        savedir=paste0(wkdir,'/',m,'/',dir,'/')
        if(file.exists(savedir) == T) { 
          setwd(savedir)
          pred.files = list.files(pattern='.*predictions.RDS')
          if (length(pred.files) > 1) {
            pred.list =mclapply(pred.files,function(x) {
              tmp = readRDS(x)
              tmp$feature = m
              if (sum(grepl('filler',dir)) == 1) {
                tmp = tmp[tmp$covariate == 'filler',]
              } else if (sum(grepl('None',dir)) == 1) {
                tmp = tmp[tmp$covariate == 'None',]
                
              }
              if (sum(grepl('covariate',colnames(tmp))) == 0) {tmp$covariate= 'None'}
              
              if (sum(grepl('fc',colnames(tmp))) == 0) {tmp$fc= 0.25}
              tmp$dir = dir
              tmp=tmp[tmp$GRP_Id %in% discovery.set$GRP_Id,]
              return(tmp)
              
              
            },mc.cores=5 )
            #pred.df.all=do.call('rbind',pred.list)[,c('GRP_Id','reported','methylation_score','model','feature','test','n.features','filter','fc')]
            pred.df.all.tmp=do.call('rbind',pred.list)
            pred.df.all.tmp$direction = dir
            feature.size.count = data.frame(table(pred.df.all.tmp$n.features),stringsAsFactors = F)
            n.base = max(feature.size.count$Freq)
            #feature.sizes= c(seq(25,200,25),seq(200,1000,100))
            #low.freq = feature.size.count[feature.size.count$Freq < n.base*0.15  | !feature.size.count$Var1 %in%  feature.sizes,]
            
            #pred.df.all.tmp$n.features = ifelse(pred.df.all.tmp$n.features %in% feature.sizes,  pred.df.all.tmp$n.features, 1100)
            pred.df.all.tmp$GRP_Id = gsub('_combined.*','',pred.df.all.tmp$GRP_Id)
            pred.df.all =rbind(pred.df.all, unique(pred.df.all.tmp))
            
          }
          
        }
        
      }
      #pred.df.all$feature = pred.df.all$covariate#m
      pred.df.all$freq = 1
      #pred.df.all = pred.df.all[pred.df.all$seed %in% c(1:10),]
      #pred.df.all$covariate = 'None'
      pred.df.all.freq =ddply(pred.df.all[,c('GRP_Id','reported','model','feature','test','n.features','mat','direction','comparison','freq','covariate','dir','fc')],
                              c('GRP_Id','reported','model','feature','test','n.features','mat','direction','comparison','covariate','dir','fc'),
                              numcolwise(sum))
      
      
      pred.df.collapse.mean = ddply(pred.df.all[,c('GRP_Id','reported','methylation_score','model','feature','test','n.features','mat','direction','comparison','covariate','dir','fc')],
                                    c('GRP_Id','reported','model','feature','test','n.features','mat','direction','comparison','covariate','dir','fc'),
                                    numcolwise(mean))
      pred.df.collapse.mean =merge(pred.df.collapse.mean,pred.df.all.freq,by=c('GRP_Id','reported','model','feature','test','n.features','mat','direction','comparison','covariate','dir','fc'))
      # pred.df.collapse.mean = pred.df.collapse.mean[pred.df.collapse.mean$freq >= 25,]
      if(max(pred.df.collapse.mean$freq) > 100) {
        print('exceed')
        stop()
      }
      #pred.df.collapse.median$validation = ifelse(pred.df.collapse.median$GRP_Id %in% validation.train$GRP_Id,'Discovery','Validation')
      mean.perf.list[[length(mean.perf.list)+1]] = pred.df.collapse.mean
      
    }
    
  }
  
  figdir = paste0(wkdir,'/cv.figures.all.cancer/')
  figdir = paste0(wkdir,'/cv.figures.all.cancer.covar/') #manuscript
  figdir = paste0(wkdir,'/cv.figures.prad.final/')
  figdir = paste0(wkdir,'/cv.figures.brca.300assess1/')
  dir.create(figdir,recursive = T)
  perfdir= paste0(wkdir,'/')
  discovery.female = discovery.set[discovery.set$Sex =='Female',]
  mean.perf.df = do.call('rbind',mean.perf.list)
  perf.list = list('mean.cv'= mean.perf.df)
  
  for (p in names(perf.list)) {
    library(ggh4x)
    targ.perf = perf.list[[p]]
    targ.perf = targ.perf#[targ.perf$GRP_Id %in% discovery.female$GRP_Id,]#[targ.perf$filler == 'All',] #[!targ.perf$feature %in% c('fraglengths.v1.tiles.100kb','fraglengths.v2.tiles.100kb','fraglengths.v3.tiles.100kb'),]
    #targ.perf$feature = 'Methylation'
    c('model','feature','test','n.features','mat','direction','comparison')
    targ.perf.list = split(targ.perf, targ.perf[,c('model','feature','test','n.features','mat','direction','comparison','covariate','fc')])
    targ.perf.list = targ.perf.list[sapply(targ.perf.list,nrow)>0 ]
    targ.perf.list =lapply(targ.perf.list, function(x) {
      tmp = x#[x$GRP_Id %in% sample.splits$GRP_Id,]
      tmp$All = auc_calc(tmp,c('Control','Cancer'))
      concordance_calc = function(tmp,sample.info){
        x = merge(tmp,sample.info)
        x$event=ifelse(x$reported =='Cancer',1,0)
        x$reported.surv = ifelse(x$reported == 'Cancer',1,0)
        library(survcomp)
        #male.weights= weightsf.females(x)
        
        ci= concordance.index(x$methylation_score, x$'censorship_time', surv.event = x$event, comppairs=10, na.rm = FALSE)#weights
        return(ci$c.index)
      }
      tmp$Concordance = concordance_calc(tmp,sample.info)
      
      # tmp1 = tmp[tmp$GRP_Id %in% aix1.samples$GRP_Id,]
      # tmp$UFiller =  auc_calc(tmp1,c('Control','Cancer'))
      # tmp1 = tmp[tmp$GRP_Id %in% aix3.samples$GRP_Id,]
      # tmp$MFiller =  auc_calc(tmp1,c('Control','Cancer'))
      return(tmp)
    } )
    
    targ.perf.all=do.call('rbind',targ.perf.list)
    targ.perf.all = targ.perf.all[order(-targ.perf.all$All),]
    rownames(targ.perf.all) = c(1:nrow(targ.perf.all))
    a1=targ.perf.all[targ.perf.all$GRP_Id == 'AIX_0011',]
    a=targ.perf.all[targ.perf.all$GRP_Id == 'AIX_0010',]
    #,'Methylated','Unmethylated'
    #
    
    afdf
    
    
    #saveRDS(targ.perf.all,paste0(perfdir,p,'.overall.performance.RDS'))
    #saveRDS(targ.perf.all,paste0(perfdir,p,'.select.overall.performance.RDS'))
    #saveRDS(targ.perf.all,paste0(perfdir,p,'.brca.300.performance.newsplit.RDS'))
    #saveRDS(targ.perf.all,paste0(perfdir,p,'.genhancer.performance.newsplit.hr.75.RDS'))
    saveRDS(targ.perf.all,paste0(perfdir,p,'.genhancer.top.performance.RDS'))
    #saveRDS(targ.perf.all,paste0(perfdir,p,'.prad.silencer.performance.RDS'))
    #saveRDS(targ.perf.all,paste0(perfdir,p,'.topprad.silencer.performance.RDS'))
    
    #targ.perf.all = readRDS(paste0(perfdir,p,'.overall.performance.RDS'))
    print(paste0(perfdir,p,'.overall.performance.RDS'))
    for (comp in c('Breast','Prostate') ) {
      for (dir in directions) {#'base.abs','bm1.abs''background.annotation','background.autosome','background.pbl'
        # targ.perf.all$test = ifelse(targ.perf.all$test %in% c('cph','coxph'),'Cox Regression','Wilcoxon Test')
        targ.perf.all.filt = targ.perf.all[targ.perf.all$direction == dir & targ.perf.all$comparison == comp ,]
        if (nrow(targ.perf.all.filt) > 0) {
          targ.perf.all.melt = reshape2::melt(targ.perf.all.filt,colnames(targ.perf.all.filt)[!colnames(targ.perf.all.filt) %in% c('All','UFiller','MFiller')])
          
          plot1 = ggplot(targ.perf.all.melt,aes(x = n.features, y = value,color = model)) + geom_point() +
            #scale_color_manual(values = group_col)+#+ ggtitle(title) +
            # scale_color_manual(values=c('rf' ='#247BA0','logreg'='#B8336A'))+
            theme_bw()+
            #geom_line()+
            
            scale_y_continuous(limits = c(0,1))+
            facet_grid2(covariate ~ mat)+
            theme(text = element_text(size=12),
                  axis.text=element_text(size=12, face = "bold"),
                  axis.title=element_text(size=14,face="bold"),
                  legend.position = "right",
                  axis.text.x=element_text(angle=90,vjust=0.5,hjust=1))+
            guides(colour = guide_legend(override.aes = list(size=4),nrow=4)) +
            xlab('Feature Type') + ylab('CV AUROC')
          
          
          png(paste0(figdir,p,'.performance.auroc.',dir,'.',comp,'.methyl.png'),height=2000,width = 3000,res = 200)
          print(plot1)
          dev.off()
          plot1 = ggplot(targ.perf.all.melt,aes(x = n.features, y = Concordance,color = model)) + geom_point() +
            #scale_color_manual(values = group_col)+#+ ggtitle(title) +
            #scale_color_manual(values=c('rf' ='#247BA0','logreg'='#B8336A','entropy' ='#3DDC97'))+
            theme_bw()+
            #geom_line()+
            scale_y_continuous(limits = c(0,1))+
            facet_grid2(covariate ~ mat)+
            theme(text = element_text(size=12),
                  axis.text=element_text(size=12, face = "bold"),
                  axis.title=element_text(size=14,face="bold"),
                  legend.position = "right",
                  axis.text.x=element_text(angle=90,vjust=0.5,hjust=1))+
            guides(colour = guide_legend(override.aes = list(size=4),nrow=4)) +
            xlab('Feature Type') + ylab('CV Concordance Index')
          
          png(paste0(figdir,p,'.performance.ci.',dir,'.',comp,'.methyl.png'),height=2000,width = 3000,res = 200)
          print(plot1)
          dev.off()
          
          
        }#[targ.perf.all$validation == 'Validation',]#[targ.perf.all$filter == filter,]
        
        
      }
      
      
    }
    
    #other covars
    oc = F
    if (oc == T) {
      for (comp in c('Breast','Prostate')[1] ) {
        for (dir in directions) {#'background.annotation','background.autosome','background.pbl'
          # targ.perf.all$test = ifelse(targ.perf.all$test %in% c('cph','coxph'),'Cox Regression','Wilcoxon Test')
          targ.perf.all.filt = targ.perf.all[targ.perf.all$direction == dir & targ.perf.all$comparison == comp,]
          if (nrow(targ.perf.all.filt) > 0) {
            targ.perf.all.melt = reshape2::melt(targ.perf.all.filt,colnames(targ.perf.all.filt)[!colnames(targ.perf.all.filt) %in% c('All','UFiller','MFiller')])
            sadf
            plot1 = ggplot(targ.perf.all.melt,aes(x = n.features, y = value,color = model)) + geom_point() +
              #scale_color_manual(values = group_col)+#+ ggtitle(title) +
              # scale_color_manual(values=c('rf' ='#247BA0','logreg'='#B8336A'))+
              theme_bw()+
              #geom_line()+
              
              scale_y_continuous(limits = c(0,1))+
              facet_grid2(covariate ~ mat)+
              theme(text = element_text(size=12),
                    axis.text=element_text(size=12, face = "bold"),
                    axis.title=element_text(size=14,face="bold"),
                    legend.position = "right",
                    axis.text.x=element_text(angle=90,vjust=0.5,hjust=1))+
              guides(colour = guide_legend(override.aes = list(size=4),nrow=4)) +
              xlab('Feature Type') + ylab('CV AUROC')
            
            
            png(paste0(figdir,p,'.performance.auroc.',dir,'.',comp,'.covaradjust.methyl.png'),height=2000,width = 3000,res = 200)
            print(plot1)
            dev.off()
            plot1 = ggplot(targ.perf.all.melt,aes(x = n.features, y = Concordance,color = model)) + geom_point() +
              #scale_color_manual(values = group_col)+#+ ggtitle(title) +
              #scale_color_manual(values=c('rf' ='#247BA0','logreg'='#B8336A','entropy' ='#3DDC97'))+
              theme_bw()+
              #geom_line()+
              scale_y_continuous(limits = c(0,1))+
              facet_grid2(covariate ~ mat)+
              theme(text = element_text(size=12),
                    axis.text=element_text(size=12, face = "bold"),
                    axis.title=element_text(size=14,face="bold"),
                    legend.position = "right",
                    axis.text.x=element_text(angle=90,vjust=0.5,hjust=1))+
              guides(colour = guide_legend(override.aes = list(size=4),nrow=4)) +
              xlab('Feature Type') + ylab('CV Concordance Index')
            
            png(paste0(figdir,p,'.performance.ci.',dir,'.',comp,'.covaradjust.methyl.png'),height=2000,width = 3000,res = 200)
            print(plot1)
            dev.off()
            
            
          }#[targ.perf.all$validation == 'Validation',]#[targ.perf.all$filter == filter,]
          
          
        }
        
        
      }
      
    }
  }
  
  
}


####selecting top performing parameter#####
library(cutpointr)
wkdir='/.mounts/labs/awadallalab/private/ncheng/cfmedip_data/cptp_samples/fragmentation/aix13.updated1/all.sample.fragmentation/breast.cancer.cv/regulatory.regions.v10.newsplit7/'
figdir = paste0(wkdir,'/cv.figures.all.cancer/')
figdir = paste0(wkdir,'/cv.figures.300/')
figdir = paste0(wkdir,'/cv.figures.genhancer.final.90/')

dir.create(figdir,recursive = T)
perfdir= paste0(wkdir,'/')
#overall.perf = readRDS(paste0(perfdir,'mean.cv.select.overall.performance.RDS'))
#overall.perf = readRDS(paste0(perfdir,'.brca.300.performance.newsplit.RDS'))
overall.perf = readRDS(paste0(perfdir,'mean.cv.genhancer.performance.newsplit.75.RDS'))
overall.perf = readRDS(paste0(perfdir,'mean.cv.prad.silencer.performance.RDS'))
overall.perf = readRDS(paste0(perfdir,'mean.cv.topprad.silencer.performance.RDS'))
overall.perf = readRDS(paste0(perfdir,'mean.cv.genhancer.top.performance.RDS'))


#overall.perf= readRDS(paste0(perfdir,'mean.cv.genhancer.performance.newsplit.hr.75.RDS'))
overall.perf = overall.perf[as.numeric(overall.perf$freq) >= 1 & overall.perf$n.features == 90 ,] ##& overall.perf$covariate == 'None'

overall.perf.list = split(overall.perf,overall.perf$comparison)
overall.perf.list.top = lapply(overall.perf.list, function(x) {
  targ.id = x$GRP_Id[1]
  return.df = x[x$GRP_Id %in% targ.id,]
  return.df = return.df[order(-return.df$All),]
  rownames(return.df) = c(1:nrow(return.df))
  return(return.df)
  
} )
lapply(overall.perf.list.top,head)

parameters = c('model','feature','test','n.features', 'mat','comparison','direction')#covariate

#top breast cancer parameters
top.perf = overall.perf.list.top[[1]]#[1,]
top.perf = top.perf[top.perf$n.features == 90,]
overall.perf.targ.breast = overall.perf
for (p in parameters) {
  print(top.perf[,p])
  overall.perf.targ.breast = overall.perf.targ.breast[overall.perf.targ.breast[,p] == top.perf[,p],]
  
}
overall.perf.targ.breast = overall.perf[overall.perf$All == max(overall.perf$All),]
score.cutoff = cutpointr(overall.perf.targ.breast$methylation_score,overall.perf.targ.breast$reported, method = maximize_metric, metric = youden)
score.cutoff.breast = score.cutoff$optimal_cutpoint
#saveRDS(score.cutoff.breast,paste0(wkdir,'breast.optimal.scorecutoff.RDS'))
#top prostate cancer parameters
top.perf = overall.perf.list.top[['Prostate']]#[1,]
top.perf = top.perf[top.perf$n.features == 100,][1,]
overall.perf.targ.prostate = overall.perf.list[['Prostate']]
for (p in parameters) {
  overall.perf.targ.prostate = overall.perf.targ.prostate[overall.perf.targ.prostate[,p] == top.perf[,p],]
  
}
score.cutoff = cutpointr(overall.perf.targ.prostate$methylation_score,overall.perf.targ.prostate$reported, method = maximize_metric, metric = youden)
score.cutoff.prostate = score.cutoff$optimal_cutpoint
#saveRDS(score.cutoff.prostate,paste0(wkdir,'prostate.optimal.scorecutoff.RDS'))
###breast cancer plotting####
library(cutpointr)

figdir = paste0(wkdir,'/cv.figures.300/')
figdir = paste0(wkdir,'/cv.figures.genhancer.base.discovery.final.90.v3/')

dir.create(figdir,recursive = T)

pred.df.targ = overall.perf.targ.breast
name = paste0(figdir,'brca.top')
merged.df.all = discovery.set
dx.all = F
score.cutoff=score.cutoff.breast
cutpoint.use=T
sample.info = discovery.set

brca.auc.plot = function(pred.df.targ,name,merged.df.all,figdir,dx.all=F,score.cutoff=0.5,cutpoint.use=F) {
  library(cenROC)
  library(cutpointr)
  library(plyr)
  library(ggplot2)
  library(ggh4x)
  library(pROC)
  library(ROCR)
  library(caret)
  tpr.fpr.calc = function(x){
    tmp1 = x
    tmp1$f = as.integer(ifelse(tmp1$reported == 'Control', 0, 1))
    tmp1$f.str = tmp1$reported
    
    tmp1 = tmp1[order(-tmp1$methylation_score),]
    case_no = nrow(tmp1[tmp1$reported!='Control',])
    control_no = nrow(tmp1[tmp1$reported=='Control',])
    auc.df =data.frame(matrix(nrow = 0, ncol=3))
    for(l in 1:nrow(tmp1)) {
      x = tmp1[1:l,]
      case_cum = nrow(x[x$reported!='Control',])
      control_cum = nrow(x[x$reported=='Control',])
      tpr = case_cum/case_no
      fpr = control_cum/control_no
      return.tmp = data.frame(l,tpr,fpr)
      auc.df = rbind(auc.df, return.tmp)
    }
    base = data.frame(l = 0, tpr = 0, fpr = 0)
    auc.df = rbind(base,auc.df)
    return(auc.df)
    #auc.df = auc.df[auc.df$fpr > 0,]
  }
  summary.calc = function(pred.df.targ.collapse,merged.df.all.tmp) {
    if (length(unique(pred.df.targ.collapse$reported)) > 1){
      auc.plot.all = tpr.fpr.calc(pred.df.targ.collapse)
      auc.plot.all$auc = auc_calc(pred.df.targ.collapse,c('Control','Cancer'))
      roc.a= pROC::roc(data=pred.df.targ.collapse,response='reported',predictor='methylation_score')
      a.confint = ci.auc(roc.a, conf.level=0.95, method=c( "bootstrap"), boot.n = 2000)
      auc.plot.all$auc.lower = a.confint[1]
      auc.plot.all$auc.upper = a.confint[3]
      concordance_calc = function(tmp,sample.info,group='reported'){
        x = merge(tmp,sample.info)
        x$event=ifelse(x[,group] =='Cancer',1,0)
        x$reported.surv = ifelse(x[,group] == 'Cancer',1,0)
        library(survcomp)
        male.weights= weightsf.females(x)
        
        ci= concordance.index(x$methylation_score, x$'censorship_time', surv.event = x$event, comppairs=10, na.rm = FALSE,weights = male.weights)#
        return(c(ci$c.index,ci$lower,ci$upper))
      }
      
      
      auc.plot.all$ci = concordance_calc(pred.df.targ.collapse, merged.df.all.tmp,'reported')[1]
      auc.plot.all$ci.lower = concordance_calc(pred.df.targ.collapse, merged.df.all.tmp,'reported')[2]
      auc.plot.all$ci.upper = concordance_calc(pred.df.targ.collapse, merged.df.all.tmp,'reported')[3]
      
      spec.threshold=data.frame(ci.se(roc.a,c(0.95,0.9)))
      
      auc.plot.all$se.spec.95 = spec.threshold[1,2]
      auc.plot.all$se.spec.95.lower = spec.threshold[1,1]
      auc.plot.all$se.spec.95.upper = spec.threshold[1,3]
      
      
      auc.plot.all$se.spec.90 = spec.threshold[2,2]
      auc.plot.all$se.spec.90.lower = spec.threshold[2,1]
      auc.plot.all$se.spec.90.upper = spec.threshold[2,3]
      
      #new
      jcutpoint.stats = coords(roc.a,cp.youden, 'threshold', ret=c("sensitivity","specificity","ppv","npv"))
      f1cutpoint.stats = coords(roc.a,cp.F1_score, 'threshold', ret=c("sensitivity","specificity","ppv","npv"))
      
      jcutpoint.stats.ci = ci.coords(roc.a,cp.youden, 'threshold', ret=c("sensitivity","specificity","ppv","npv"),conf.level=0.95)
      f1cutpoint.stats.ci = ci.coords(roc.a,cp.F1_score, 'threshold', ret=c("sensitivity","specificity","ppv","npv"),conf.level=0.95)
      
      
      auc.plot.all$jcutpoint.sens = jcutpoint.stats$sensitivity
      auc.plot.all$jcutpoint.sens.lower = jcutpoint.stats.ci$sensitivity[1]
      auc.plot.all$jcutpoint.sens.upper = jcutpoint.stats.ci$sensitivity[3]
      
      auc.plot.all$jcutpoint.spec = jcutpoint.stats$specificity
      auc.plot.all$jcutpoint.spec.lower = jcutpoint.stats.ci$specificity[1]
      auc.plot.all$jcutpoint.spec.upper = jcutpoint.stats.ci$specificity[3]
      
      auc.plot.all$jcutpoint.ppv = jcutpoint.stats$ppv
      auc.plot.all$jcutpoint.ppv.lower = jcutpoint.stats.ci$ppv[1]
      auc.plot.all$jcutpoint.ppv.upper = jcutpoint.stats.ci$ppv[3]
      
      auc.plot.all$jcutpoint.npv = jcutpoint.stats$npv
      auc.plot.all$jcutpoint.npv.lower = jcutpoint.stats.ci$npv[1]
      auc.plot.all$jcutpoint.npv.upper = jcutpoint.stats.ci$npv[3]
      
      
      auc.plot.all$f1cutpoint.sens = f1cutpoint.stats$sensitivity
      auc.plot.all$f1cutpoint.sens.lower = f1cutpoint.stats.ci$sensitivity[1]
      auc.plot.all$f1cutpoint.sens.upper = f1cutpoint.stats.ci$sensitivity[3]
      
      auc.plot.all$f1cutpoint.spec = f1cutpoint.stats$specificity
      auc.plot.all$f1cutpoint.spec.lower = f1cutpoint.stats.ci$specificity[1]
      auc.plot.all$f1cutpoint.spec.upper = f1cutpoint.stats.ci$specificity[3]
      
      auc.plot.all$f1cutpoint.ppv = f1cutpoint.stats$ppv
      auc.plot.all$f1cutpoint.ppv.lower = f1cutpoint.stats.ci$ppv[1]
      auc.plot.all$f1cutpoint.ppv.upper = f1cutpoint.stats.ci$ppv[3]
      
      auc.plot.all$f1cutpoint.npv = f1cutpoint.stats$npv
      auc.plot.all$f1cutpoint.npv.lower = f1cutpoint.stats.ci$npv[1]
      auc.plot.all$f1cutpoint.npv.upper = f1cutpoint.stats.ci$npv[3]
      
      auc.plot.all$jcutpoint.value = cp.youden
      auc.plot.all$f1cutpoint.value = cp.F1_score
      
      #npv ppv weighting
      mean.perf.df.targ.tmp.merged = pred.df.targ.collapse
      mean.perf.df.targ.tmp.merged$Event=ifelse(mean.perf.df.targ.tmp.merged$reported == 'Control',0,1)
      mean.perf.df.targ.tmp.merged$group = mean.perf.df.targ.tmp.merged$reported
      mean.perf.df.targ.tmp.merged = mean.perf.df.targ.tmp.merged[!is.na(mean.perf.df.targ.tmp.merged$methylation_score),]
      mean.perf.df.targ.tmp.merged$Risk.group = ifelse(mean.perf.df.targ.tmp.merged$methylation_score > cp.youden,'High Predicted Risk','Low Predicted Risk')
      mean.perf.df.targ.tmp.merged$Risk.group = factor(as.character(mean.perf.df.targ.tmp.merged$Risk.group ),levels = c('Low Predicted Risk','High Predicted Risk'))
      mean.perf.df.targ.tmp.merged$Event= ifelse(mean.perf.df.targ.tmp.merged$reported == 'Control',0,1 )
      subject1 = mean.perf.df.targ.tmp.merged
      female.ohs.qx.weighting = function(subject,combined.info.all) {
        #subject = merge(mean.perf.df.targ.tmp.merged,sample.info[,c('GRP_Id','age_group','Family.history.breast','Alch.con.group','Smoking.Frequency')],by='GRP_Id')
        combined.info.all.male = combined.info.all[combined.info.all$Sex == 'Female',]
        combined.info.all.male = combined.info.all.male[combined.info.all.male$Cancer %in% c('Control','Breast'),]
        combined.info.all.male$Smoking.Frequency = as.character(combined.info.all.male$Smoking.Frequency)
        combined.info.all.male[is.na(combined.info.all.male$Smoking.Frequency),'Smoking.Frequency'] = 'Never'
        #bmi.short.groups = unique(combined.info.all.male$bmi.group.short)
        #bmi.long.groups = unique(combined.info.all.male$bmi.group.long)
        age.groups =  unique(combined.info.all.male$age_group)
        fh.groups = unique(combined.info.all.male$Family.history.breast)
        alc.groups = unique(combined.info.all.male$Alch.con.group)
        smk.groups = unique(combined.info.all.male$Smoking.Frequency)
        control.samples = subject[subject$reported == 'Control',]
        cancer.samples = subject[subject$reported != 'Control',]
        
        control.samples = combined.info.all[combined.info.all$GRP_Id %in% control.samples$GRP_Id,]
        cancer.samples = combined.info.all[combined.info.all$GRP_Id %in% cancer.samples$GRP_Id,]
        
        ohs.pop.samples = combined.info.all[combined.info.all$Cohort == 'EpiCan', ]
        control.weight.samples = NULL
        for (age in age.groups ) {
          for (fh in fh.groups ){
            for (alc in alc.groups) {
              for (smk in smk.groups ){
                group.combination = paste(age,fh,alc,smk,sep='.')
                targ.control.cohort = control.samples[which(control.samples$age_group == age &
                                                              control.samples$Family.history.breast == fh &
                                                              control.samples$Alch.con.group == alc  &
                                                              control.samples$Smoking.Frequency == smk ) ,]
                if (nrow(targ.control.cohort)>0){
                  cohort.freq =  ohs.pop.samples[which(ohs.pop.samples$age_group == age &
                                                         ohs.pop.samples$Family.history.breast == fh &
                                                         ohs.pop.samples$Alch.con.group == alc &
                                                         ohs.pop.samples$Smoking.Frequency == smk ),]
                  
                  targ.control.cohort$weight= nrow(cohort.freq)/nrow(targ.control.cohort)
                  
                  control.weight.samples = rbind(control.weight.samples,targ.control.cohort)
                }
                
                
              }
            }
          }
        }
        
        if (nrow(control.weight.samples) < nrow(control.samples)) {
          exc.samples = control.samples[!control.samples$GRP_Id %in% control.weight.samples$GRP_Id,]
          for (r in 1:nrow(exc.samples)) {
            print(r)
            targ.sample = exc.samples[r,]
            #fh, alc, smk, bmi
            non.na.vars = c('age_group','ALC_CUR_FREQ','SMK_CIG_STATUS','DIS_CANCER_F_EVER')
            non.na.vars = non.na.vars[which(!is.na(targ.sample[,non.na.vars]))]
            
            cohort.freq =  ohs.pop.samples
            for (v in non.na.vars) {
              cohort.freq =  cohort.freq[which(cohort.freq[,v] == targ.sample[,v]) ,]
              targ.freq = control.samples[which(control.samples[,v] == targ.sample[,v])  ,]
            }
            targ.sample$weight = nrow(cohort.freq)/nrow(targ.freq)
            control.weight.samples = rbind(control.weight.samples,targ.sample)
          }
        }
        cancer.samples$weight=1
        return.df = rbind(control.weight.samples,cancer.samples)
        return.df = return.df[order(match(return.df$GRP_Id,subject$GRP_Id)),]
        return(return.df$weight)
        
      }
      
      x = merge(mean.perf.df.targ.tmp.merged,sample.info)
      x$event=ifelse(x[,'reported'] =='Cancer',1,0)
      x$reported.surv = ifelse(x[,'reported'] == 'Cancer',1,0)
      library(survcomp)
      x = x[order(match(x$GRP_Id,mean.perf.df.targ.tmp.merged$GRP_Id )),]
      male.weights= weightsf.females(x)
      mean.perf.df.targ.tmp.merged = x[,c(colnames(mean.perf.df.targ.tmp.merged))]
      
      mean.perf.df.targ.tmp.merged$weighting = male.weights
      
      #mean.perf.df.targ.tmp.merged$weighting = female.ohs.qx.weighting(mean.perf.df.targ.tmp.merged,merged.df.all.tmp[merged.df.all.tmp$GRP_Id %in% mean.perf.df.targ.tmp.merged$GRP_Id,])
      
      
      #normalized weighting
      weighted.df = NULL
      for (j in 1:nrow(mean.perf.df.targ.tmp.merged)) {
        weighted.df.tmp = mean.perf.df.targ.tmp.merged[rep(j,mean.perf.df.targ.tmp.merged[j,'weighting']),]
        weighted.df=rbind(weighted.df,weighted.df.tmp)
      }
      roc.a.weighted= pROC::roc(data=weighted.df,response='reported',predictor='methylation_score')
      
      jcutpoint.stats.weighted = coords(roc.a.weighted,cp.youden, 'threshold', ret=c("sensitivity","specificity","ppv","npv"))
      
      jcutpoint.stats.ci.weighted = ci.coords(roc.a.weighted,cp.youden, 'threshold', ret=c("sensitivity","specificity","ppv","npv"),conf.level=0.95)
      
      
      
      auc.plot.all$jcutpoint.sens.weighted = jcutpoint.stats.weighted$sensitivity
      auc.plot.all$jcutpoint.sens.lower.weighted = jcutpoint.stats.ci.weighted$sensitivity[1]
      auc.plot.all$jcutpoint.sens.upper.weighted = jcutpoint.stats.ci.weighted$sensitivity[3]
      
      auc.plot.all$jcutpoint.spec.weighted = jcutpoint.stats.weighted$specificity
      auc.plot.all$jcutpoint.spec.lower.weighted = jcutpoint.stats.ci.weighted$specificity[1]
      auc.plot.all$jcutpoint.spec.upper.weighted = jcutpoint.stats.ci.weighted$specificity[3]
      
      auc.plot.all$jcutpoint.ppv.weighted = jcutpoint.stats.weighted$ppv
      auc.plot.all$jcutpoint.ppv.lower.weighted = jcutpoint.stats.ci.weighted$ppv[1]
      auc.plot.all$jcutpoint.ppv.upper.weighted = jcutpoint.stats.ci.weighted$ppv[3]
      
      auc.plot.all$jcutpoint.npv.weighted = jcutpoint.stats.weighted$npv
      auc.plot.all$jcutpoint.npv.lower.weighted = jcutpoint.stats.ci.weighted$npv[1]
      auc.plot.all$jcutpoint.npv.upper.weighted = jcutpoint.stats.ci.weighted$npv[3]
      
      
      
      #
      auc.plot.all$count = length(unique(merged.df.all.tmp[merged.df.all.tmp$GRP_Id %in% pred.df.targ.collapse$GRP_Id & merged.df.all.tmp$group == 'Cancer','GRP_Id']))
      #time dependent auroc
      time.roc.calc = function(tmp){
        x = merge(tmp,sample.info)
        x$event=ifelse(x[,'reported'] =='Cancer',1,0)
        x$reported.surv = ifelse(x[,'reported'] == 'Cancer',1,0)
        library(survcomp)
        male.weights= weightsf.females(x)
        return.roct = NULL
        targ.times = seq(30,3000,100)
        targ.times = targ.times[targ.times < max(x$censorship_time)]
        for (time in targ.times) {
          roct = cenROC(Y=x$censorship_time, 
                        M=x$methylation_score, 
                        censor=x$event, 
                        t=time,
                        #  h = male.weights,
                        bw = "NR", method = "emp",
                        ktype = "normal", ktype1 = "normal", B = 100, alpha = 0.05, plot = "TRUE")
          tmp.df = data.frame(time = time, 
                              auct = roct$AUC[1],
                              auct.lower = roct$AUC[3],
                              auct.upper = roct$AUC[4])
          return.roct =rbind(return.roct,tmp.df)
        }
        return(return.roct)
      }
      troc = time.roc.calc(pred.df.targ.collapse)
      return.list = list(roc.a, auc.plot.all,troc)
      return(return.list)
    } else {
      troc=data.frame('time'=0 ,AUC=0, LCL=0, UCL=0)
      auc.plot.all = data.frame(l=0,
                                count = 0,
                                tpr=0,
                                fpr=0,
                                auc=0,
                                auc.lower=0,
                                auc.upper=0,
                                ci=0,
                                ci.lower=0,
                                ci.upper=0,
                                se.spec.95=0,
                                se.spec.95.lower=0,
                                se.spec.95.upper=0,
                                se.spec.90=0,
                                se.spec.90.lower=0,
                                se.spec.90.upper=0,
                                
                                jcutpoint.sens = 0,
                                jcutpoint.sens.lower = 0,
                                jcutpoint.sens.upper = 0,
                                
                                jcutpoint.spec = 0,
                                jcutpoint.spec.lower =0,
                                jcutpoint.spec.upper = 0,
                                
                                jcutpoint.ppv =0,
                                jcutpoint.ppv.lower = 0,
                                jcutpoint.ppv.upper = 0,
                                
                                jcutpoint.npv = 0,
                                jcutpoint.npv.lower = 0,
                                jcutpoint.npv.upper = 0,
                                
                                
                                f1cutpoint.sens = 0,
                                f1cutpoint.sens.lower = 0,
                                f1cutpoint.sens.upper = 0,
                                
                                f1cutpoint.spec =0,
                                f1cutpoint.spec.lower = 0,
                                f1cutpoint.spec.upper = 0,
                                
                                f1cutpoint.ppv = 0,
                                f1cutpoint.ppv.lower = 0,
                                f1cutpoint.ppv.upper = 0,
                                
                                f1cutpoint.npv = 0,
                                f1cutpoint.npv.lower = 0,
                                f1cutpoint.npv.upper = 0,
                                
                                jcutpoint.value = cp.youden,
                                f1cutpoint.value = cp.F1_score,
                                
                                jcutpoint.sens.weighted = 0,
                                jcutpoint.sens.lower.weighted = 0,
                                jcutpoint.sens.upper.weighted = 0,
                                
                                jcutpoint.spec.weighted = 0,
                                jcutpoint.spec.lower.weighted =0,
                                jcutpoint.spec.upper.weighted = 0,
                                
                                jcutpoint.ppv.weighted =0,
                                jcutpoint.ppv.lower.weighted = 0,
                                jcutpoint.ppv.upper.weighted = 0,
                                
                                jcutpoint.npv.weighted = 0,
                                jcutpoint.npv.lower.weighted = 0,
                                jcutpoint.npv.upper.weighted = 0
                                
                                
                                
                                
      )
      
      return(list(NULL,auc.plot.all,troc))
    }
    
    
  }
  
  summary.calc.dxtime = function(pred.df.targ.collapse,merged.df.all.tmp) {
    if (length(unique(pred.df.targ.collapse$reported)) > 1){
      auc.plot.all = tpr.fpr.calc(pred.df.targ.collapse)
      auc.plot.all$auc = auc_calc(pred.df.targ.collapse,c('Control','Cancer'))
      roc.a= pROC::roc(data=pred.df.targ.collapse,response='reported',predictor='methylation_score')
      a.confint = ci.auc(roc.a, conf.level=0.95, method=c( "bootstrap"), boot.n = 2000)
      auc.plot.all$auc.lower = a.confint[1]
      auc.plot.all$auc.upper = a.confint[3]
      concordance_calc = function(tmp,sample.info,group='reported'){
        x = merge(tmp,sample.info)
        x$event=ifelse(x[,group] =='Cancer',1,0)
        x$reported.surv = ifelse(x[,group] == 'Cancer',1,0)
        library(survcomp)
        male.weights= weightsf.females(x)
        
        ci= concordance.index(x$methylation_score, x$'censorship_time', surv.event = x$event, comppairs=10, na.rm = FALSE,weights = male.weights)#
        return(c(ci$c.index,ci$lower,ci$upper))
      }
      
      
      auc.plot.all$ci = concordance_calc(pred.df.targ.collapse, merged.df.all.tmp,'reported')[1]
      auc.plot.all$ci.lower = concordance_calc(pred.df.targ.collapse, merged.df.all.tmp,'reported')[2]
      auc.plot.all$ci.upper = concordance_calc(pred.df.targ.collapse, merged.df.all.tmp,'reported')[3]
      
      spec.threshold=data.frame(ci.se(roc.a,c(0.95,0.9)))
      
      auc.plot.all$se.spec.95 = spec.threshold[1,2]
      auc.plot.all$se.spec.95.lower = spec.threshold[1,1]
      auc.plot.all$se.spec.95.upper = spec.threshold[1,3]
      
      
      auc.plot.all$se.spec.90 = spec.threshold[2,2]
      auc.plot.all$se.spec.90.lower = spec.threshold[2,1]
      auc.plot.all$se.spec.90.upper = spec.threshold[2,3]
      
      #new
      jcutpoint.stats = coords(roc.a,cp.youden, 'threshold', ret=c("sensitivity","specificity","ppv","npv"))
      f1cutpoint.stats = coords(roc.a,cp.F1_score, 'threshold', ret=c("sensitivity","specificity","ppv","npv"))
      
      jcutpoint.stats.ci = ci.coords(roc.a,cp.youden, 'threshold', ret=c("sensitivity","specificity","ppv","npv"),conf.level=0.95)
      f1cutpoint.stats.ci = ci.coords(roc.a,cp.F1_score, 'threshold', ret=c("sensitivity","specificity","ppv","npv"),conf.level=0.95)
      
      
      auc.plot.all$jcutpoint.sens = jcutpoint.stats$sensitivity
      auc.plot.all$jcutpoint.sens.lower = jcutpoint.stats.ci$sensitivity[1]
      auc.plot.all$jcutpoint.sens.upper = jcutpoint.stats.ci$sensitivity[3]
      
      auc.plot.all$jcutpoint.spec = jcutpoint.stats$specificity
      auc.plot.all$jcutpoint.spec.lower = jcutpoint.stats.ci$specificity[1]
      auc.plot.all$jcutpoint.spec.upper = jcutpoint.stats.ci$specificity[3]
      
      auc.plot.all$jcutpoint.ppv = jcutpoint.stats$ppv
      auc.plot.all$jcutpoint.ppv.lower = jcutpoint.stats.ci$ppv[1]
      auc.plot.all$jcutpoint.ppv.upper = jcutpoint.stats.ci$ppv[3]
      
      auc.plot.all$jcutpoint.npv = jcutpoint.stats$npv
      auc.plot.all$jcutpoint.npv.lower = jcutpoint.stats.ci$npv[1]
      auc.plot.all$jcutpoint.npv.upper = jcutpoint.stats.ci$npv[3]
      
      
      auc.plot.all$f1cutpoint.sens = f1cutpoint.stats$sensitivity
      auc.plot.all$f1cutpoint.sens.lower = f1cutpoint.stats.ci$sensitivity[1]
      auc.plot.all$f1cutpoint.sens.upper = f1cutpoint.stats.ci$sensitivity[3]
      
      auc.plot.all$f1cutpoint.spec = f1cutpoint.stats$specificity
      auc.plot.all$f1cutpoint.spec.lower = f1cutpoint.stats.ci$specificity[1]
      auc.plot.all$f1cutpoint.spec.upper = f1cutpoint.stats.ci$specificity[3]
      
      auc.plot.all$f1cutpoint.ppv = f1cutpoint.stats$ppv
      auc.plot.all$f1cutpoint.ppv.lower = f1cutpoint.stats.ci$ppv[1]
      auc.plot.all$f1cutpoint.ppv.upper = f1cutpoint.stats.ci$ppv[3]
      
      auc.plot.all$f1cutpoint.npv = f1cutpoint.stats$npv
      auc.plot.all$f1cutpoint.npv.lower = f1cutpoint.stats.ci$npv[1]
      auc.plot.all$f1cutpoint.npv.upper = f1cutpoint.stats.ci$npv[3]
      
      
      
      auc.plot.all$jcutpoint.value = cp.youden
      auc.plot.all$f1cutpoint.value = cp.F1_score
      
      #npv ppv weighting
      mean.perf.df.targ.tmp.merged = pred.df.targ.collapse
      mean.perf.df.targ.tmp.merged$Event=ifelse(mean.perf.df.targ.tmp.merged$reported == 'Control',0,1)
      mean.perf.df.targ.tmp.merged$group = mean.perf.df.targ.tmp.merged$reported
      mean.perf.df.targ.tmp.merged = mean.perf.df.targ.tmp.merged[!is.na(mean.perf.df.targ.tmp.merged$methylation_score),]
      mean.perf.df.targ.tmp.merged$Risk.group = ifelse(mean.perf.df.targ.tmp.merged$methylation_score >= cp.youden,'High Predicted Risk','Low Predicted Risk')
      mean.perf.df.targ.tmp.merged$Risk.group = factor(as.character(mean.perf.df.targ.tmp.merged$Risk.group ),levels = c('Low Predicted Risk','High Predicted Risk'))
      mean.perf.df.targ.tmp.merged$Event= ifelse(mean.perf.df.targ.tmp.merged$reported == 'Control',0,1 )
      subject1 = mean.perf.df.targ.tmp.merged
      female.ohs.qx.weighting = function(subject1,combined.info.all) {
        subject = merge(subject1, sample.info.filt.pretime[,c('GRP_Id','Smoking.Frequency','age_group','Family.history.breast','Alch.con.group')],by='GRP_Id')
        #subject = merge(mean.perf.df.targ.tmp.merged,sample.info[,c('GRP_Id','age_group','Family.history.breast','Alch.con.group','Smoking.Frequency')],by='GRP_Id')
        combined.info.all.male = combined.info.all[combined.info.all$Sex == 'Male',]
        combined.info.all.male$Smoking.Frequency = as.character(combined.info.all.male$Smoking.Frequency)
        combined.info.all.male[is.na(combined.info.all.male$Smoking.Frequency),'Smoking.Frequency'] = 'Never'
        #bmi.short.groups = unique(combined.info.all.male$bmi.group.short)
        #bmi.long.groups = unique(combined.info.all.male$bmi.group.long)
        age.groups =  unique(combined.info.all.male$age_group)
        fh.groups = unique(combined.info.all.male$Family.history.breast)
        alc.groups = unique(combined.info.all.male$Alch.con.group)
        smk.groups = unique(combined.info.all.male$Smoking.Frequency)
        control.samples = subject[subject$reported == 'Control',]
        cancer.samples = subject[subject$reported != 'Control',]
        
        ohs.pop.samples = combined.info.all[combined.info.all$Cohort == 'EpiCan', ]
        control.weight.samples = NULL
        for (age in age.groups ) {
          for (fh in fh.groups ){
            for (alc in alc.groups) {
              for (smk in smk.groups ){
                group.combination = paste(age,fh,alc,smk,sep='.')
                targ.control.cohort = control.samples[control.samples$age_group == age &
                                                        control.samples$Family.history.breast == fh &
                                                        control.samples$Alch.con.group == alc &
                                                        control.samples$Smoking.Frequency == smk ,]
                if (nrow(targ.control.cohort)>0){
                  cohort.freq =  ohs.pop.samples[ohs.pop.samples$age_group == age &
                                                   ohs.pop.samples$Family.history.breast == fh &
                                                   ohs.pop.samples$Alch.con.group == alc &
                                                   ohs.pop.samples$Smoking.Frequency == smk ,]
                  
                  targ.control.cohort$weight= nrow(cohort.freq)/nrow(targ.control.cohort)
                  
                  control.weight.samples = rbind(control.weight.samples,targ.control.cohort)
                }
                
                
              }
            }
          }
        }
        cancer.samples$weight=1
        return.df = rbind(control.weight.samples,cancer.samples)
        return.df = return.df[match(subject1$GRP_Id,return.df$GRP_Id),]
        return(return.df$weight)
        
      }
      
      mean.perf.df.targ.tmp.merged$weighting = female.ohs.qx.weighting(mean.perf.df.targ.tmp.merged, combined.info.all)
      
      
      #normalized weighting
      weighted.df = NULL
      for (j in 1:nrow(mean.perf.df.targ.tmp.merged)) {
        weighted.df.tmp = mean.perf.df.targ.tmp.merged[rep(j,mean.perf.df.targ.tmp.merged[j,'weighting']),]
        weighted.df=rbind(weighted.df,weighted.df.tmp)
      }
      roc.a.weighted= pROC::roc(data=weighted.df,response='reported',predictor='methylation_score')
      
      jcutpoint.stats.weighted = coords(roc.a.weighted,cp.youden, 'threshold', ret=c("sensitivity","specificity","ppv","npv"))
      
      jcutpoint.stats.ci.weighted = ci.coords(roc.a.weighted,cp.youden, 'threshold', ret=c("sensitivity","specificity","ppv","npv"),conf.level=0.95)
      
      
      
      
      auc.plot.all$jcutpoint.sens.weighted = jcutpoint.stats.weighted$sensitivity
      auc.plot.all$jcutpoint.sens.lower.weighted = jcutpoint.stats.ci.weighted$sensitivity[1]
      auc.plot.all$jcutpoint.sens.upper.weighted = jcutpoint.stats.ci.weighted$sensitivity[3]
      
      auc.plot.all$jcutpoint.spec.weighted = jcutpoint.stats.weighted$specificity
      auc.plot.all$jcutpoint.spec.lower.weighted = jcutpoint.stats.ci.weighted$specificity[1]
      auc.plot.all$jcutpoint.spec.upper.weighted = jcutpoint.stats.ci.weighted$specificity[3]
      
      auc.plot.all$jcutpoint.ppv.weighted = jcutpoint.stats.weighted$ppv
      auc.plot.all$jcutpoint.ppv.lower.weighted = jcutpoint.stats.ci.weighted$ppv[1]
      auc.plot.all$jcutpoint.ppv.upper.weighted = jcutpoint.stats.ci.weighted$ppv[3]
      
      auc.plot.all$jcutpoint.npv.weighted = jcutpoint.stats.weighted$npv
      auc.plot.all$jcutpoint.npv.lower.weighted = jcutpoint.stats.ci.weighted$npv[1]
      auc.plot.all$jcutpoint.npv.upper.weighted = jcutpoint.stats.ci.weighted$npv[3]
      
      
      
      #
      auc.plot.all$count = length(unique(merged.df.all.tmp[merged.df.all.tmp$GRP_Id %in% pred.df.targ.collapse$GRP_Id & merged.df.all.tmp$group == 'Cancer','GRP_Id']))
      #time dependent auroc
      return.list = list(roc.a, auc.plot.all)
      return(return.list)
    } else {
      troc=data.frame('time'=0 ,AUC=0, LCL=0, UCL=0)
      auc.plot.all = data.frame(l=0,
                                count = 0,
                                tpr=0,
                                fpr=0,
                                auc=0,
                                auc.lower=0,
                                auc.upper=0,
                                ci=0,
                                ci.lower=0,
                                ci.upper=0,
                                se.spec.95=0,
                                se.spec.95.lower=0,
                                se.spec.95.upper=0,
                                se.spec.90=0,
                                se.spec.90.lower=0,
                                se.spec.90.upper=0,
                                
                                jcutpoint.sens.weighted = 0,
                                jcutpoint.sens.lower.weighted = 0,
                                jcutpoint.sens.upper.weighted = 0,
                                
                                jcutpoint.spec.weighted = 0,
                                jcutpoint.spec.lower.weighted =0,
                                jcutpoint.spec.upper.weighted = 0,
                                
                                jcutpoint.ppv.weighted =0,
                                jcutpoint.ppv.lower.weighted = 0,
                                jcutpoint.ppv.upper.weighted = 0,
                                
                                jcutpoint.npv.weighted = 0,
                                jcutpoint.npv.lower.weighted = 0,
                                jcutpoint.npv.upper.weighted = 0
                                
                                
                                
      )
      
      return(list(NULL,auc.plot.all))
    }
    
    
  }
  
  
  pred.df.targ.collapse = ddply(pred.df.targ[,c('GRP_Id','reported','methylation_score')],c('GRP_Id','reported'),numcolwise(median))
  pred.df.targ.collapse= merge(pred.df.targ.collapse, merged.df.all[,c('GRP_Id','ResearchId')],by='GRP_Id')
  auc.plot.all = tpr.fpr.calc(pred.df.targ.collapse)
  auc.plot.all$gleason = 'All'
  auc.plot.all$auc = auc_calc(pred.df.targ.collapse,labels= c('Control','Cancer'))
  combined.auroc = NULL
  auc.t.combined = NULL
  
  
  for (c in c('All')) {
    print(c)
    if (c == 'All'){
      merged.df.all.tmp = merged.df.all
    } else {
      merged.df.all.tmp = merged.df.all[merged.df.all$filler == c,]
    }
    library(cutpointr)
    
    
    pred.df.targ.collapse.all=pred.df.targ.collapse
    if (cutpoint.use == F) {
      cp.youden= cutpointr(pred.df.targ.collapse.all$methylation_score,pred.df.targ.collapse.all$reported, method = maximize_metric, metric = youden)$optimal_cutpoint
      cp.F1_score= cutpointr(pred.df.targ.collapse.all$methylation_score,pred.df.targ.collapse.all$reported, method = maximize_metric, metric = F1_score)$optimal_cutpoint
      print(cp.youden)
    } else {
      cp.youden= score.cutoff
      cp.F1_score= cutpointr(pred.df.targ.collapse.all$methylation_score,pred.df.targ.collapse.all$reported, method = maximize_metric, metric = F1_score)$optimal_cutpoint
      print(cp.youden)
    }
    #plotting dx time grouped by 2 years
    dx12 = merged.df.all.tmp[merged.df.all.tmp$Diagnosis_Time %in% c('Control','0-1','1-2','0-2'),]
    dx34 = merged.df.all.tmp[merged.df.all.tmp$Diagnosis_Time %in% c('Control','2-3','3-4','2-4'),]
    dx45 = merged.df.all.tmp[merged.df.all.tmp$Diagnosis_Time %in% c('Control','4-5','5+','4-6'),]
    dx56 = merged.df.all.tmp[merged.df.all.tmp$Diagnosis_Time %in% c('Control','6-8','8-10'),]
    #dx45 = merged.df.all.tmp[merged.df.all.tmp$Diagnosis_Time %in% c('Control','4-5','5+','4-6'),]
    
    
    auc.plot.all = tpr.fpr.calc(pred.df.targ.collapse)
    if (dx.all == T){
      auc.plot.all$dx.time = 'All'
      auc.plot.all$auc = auc_calc(pred.df.targ.collapse,c('Control','Cancer'))
      auc.plot.all.dx12 = tpr.fpr.calc(pred.df.targ.collapse[pred.df.targ.collapse$GRP_Id %in% dx12$GRP_Id,])
      auc.plot.all.dx12$dx.time = '0-2'
      auc.plot.all.dx12$auc = auc_calc(pred.df.targ.collapse[pred.df.targ.collapse$GRP_Id %in% dx12$GRP_Id,],c('Control','Cancer'))
      
      print('dx tmie 2 years')
      auc.plot.all.dx34 = tpr.fpr.calc(pred.df.targ.collapse[pred.df.targ.collapse$GRP_Id %in% dx34$GRP_Id,])
      auc.plot.all.dx34$dx.time = '2-4'
      auc.plot.all.dx34$auc = auc_calc(pred.df.targ.collapse[pred.df.targ.collapse$GRP_Id %in% dx34$GRP_Id,],c('Control','Cancer'))
      
      
      auc.plot.all.dx45 = tpr.fpr.calc(pred.df.targ.collapse[pred.df.targ.collapse$GRP_Id %in% dx45$GRP_Id,])
      auc.plot.all.dx45$dx.time = '4-6'
      auc.plot.all.dx45$auc = auc_calc(pred.df.targ.collapse[pred.df.targ.collapse$GRP_Id %in% dx45$GRP_Id,],c('Control','Cancer'))
      
      auc.plot.all.dx56 = tpr.fpr.calc(pred.df.targ.collapse[pred.df.targ.collapse$GRP_Id %in% dx56$GRP_Id,])
      auc.plot.all.dx56$dx.time = '6-10'
      auc.plot.all.dx56$auc = auc_calc(pred.df.targ.collapse[pred.df.targ.collapse$GRP_Id %in% dx56$GRP_Id,],c('Control','Cancer'))
      
      
      auc.res1 = rbind(auc.plot.all[1,c('auc','dx.time')],
                       auc.plot.all.dx12[1,c('auc','dx.time')],
                       auc.plot.all.dx34[1,c('auc','dx.time')],
                       auc.plot.all.dx45[1,c('auc','dx.time')],
                       auc.plot.all.dx56[1,c('auc','dx.time')])
      colnames(auc.res1) = c('auc','subgroup')
      diagnosis_time_colors1 = c('#7A797C',"#048BA8",'#AAF683','#FFD97D','#C8553D')
      names(diagnosis_time_colors1) = c('All','0-2','2-4','4-6','6-10')
      combined.auc.plot.all = rbind(auc.plot.all,auc.plot.all.dx12,auc.plot.all.dx34,auc.plot.all.dx45,auc.plot.all.dx56)
      plot1 = ggplot(combined.auc.plot.all,aes(x = fpr, y = tpr, col = dx.time)) + geom_line(linewidth=1) +
        scale_color_manual(values = c(diagnosis_time_colors1))+ #+ ggtitle(title) +
        theme_bw()+ 
        theme(text = element_text(size=8),
              axis.text=element_text(size=8),
              axis.title=element_text(size=14,face="bold"),
              legend.position = "none")+ guides(colour = guide_legend(override.aes = list(size=4),nrow=4)) + xlab('False Positive Rate') + ylab('Sensitivity')
      
      png(paste0(name,'.filler.',c,'.auroc.dxtime2.png'),height = 1000,width=1000,res = 300)
      print(plot1)
      dev.off()
      
      plot1 = ggplot(combined.auc.plot.all,aes(x = fpr, y = tpr, col = dx.time)) + geom_line(linewidth=1) +
        scale_color_manual(values = c(diagnosis_time_colors1))+ #+ ggtitle(title) +
        theme_bw()+ 
        theme(text = element_text(size=8),
              axis.text=element_text(size=8, face = "bold"),
              axis.title=element_text(size=14,face="bold"),
              legend.position = "right")+ guides(colour = guide_legend(override.aes = list(size=4),ncol=1)) + xlab('False Positive Rate') + ylab('Sensitivity')
      
      png(paste0(name,'.filler.',c,'.auroc.dxtim2e.labs.png'),height = 1000,width=1000,res = 300)
      print(plot1)
      dev.off()
      
      
      
    } else {
      auc.res1 = NULL
    }
    
    library(pROC)
    #dx time 1 year
    print('dx tmie 1 years')
    diagnosis_time_grouping = function(diagnosis_time) {
      tmp = ifelse(diagnosis_time > 3285, '9+', diagnosis_time)
      tmp = ifelse(diagnosis_time <= 3285, '8-9', tmp)
      tmp = ifelse(diagnosis_time <= 2920, '7-8', tmp)
      tmp = ifelse(diagnosis_time <= 2555, '6-7', tmp)
      tmp = ifelse(diagnosis_time <= 2190, '5-6', tmp)
      tmp = ifelse(diagnosis_time <= 1825, '4-5', tmp)
      tmp = ifelse(diagnosis_time <= 1460, '3-4', tmp)
      tmp = ifelse(diagnosis_time <= 1095, '2-3', tmp)
      tmp = ifelse(diagnosis_time <= 730, '1-2', tmp)
      tmp = ifelse(diagnosis_time <= 365, '0-1', tmp)
      tmp[is.na(diagnosis_time)] = 'Control'
      diagnosis_time_groups = c('Control', '0-1','1-2','2-3','3-4','4-5','5-6','6-7','7-8','8-9','9+')
      tmp = factor(tmp, levels = diagnosis_time_groups) 
      return(tmp)
    }
    if (c == 'All'){
      merged.df.all.tmp = merged.df.all
    } else {
      merged.df.all.tmp = merged.df.all[merged.df.all$filler == c,]
    }
    
    #plotting by dx time
    merged.df.all.tmp$Diagnosis_Time = diagnosis_time_grouping(merged.df.all.tmp$diff_in_days)
    dx1 = merged.df.all.tmp[merged.df.all.tmp$Diagnosis_Time %in% c('Control','0-1'),]
    dx2 = merged.df.all.tmp[merged.df.all.tmp$Diagnosis_Time %in% c('Control','1-2'),]
    dx3 = merged.df.all.tmp[merged.df.all.tmp$Diagnosis_Time %in% c('Control','2-3'),]
    dx4 = merged.df.all.tmp[merged.df.all.tmp$Diagnosis_Time %in% c('Control','3-4'),]
    dx5 = merged.df.all.tmp[merged.df.all.tmp$Diagnosis_Time %in% c('Control','4-5'),]
    
    
    
    pred.df.targ.collapse$reported = factor(pred.df.targ.collapse$reported,levels = c('Control','Cancer'))
    auc.plot.all = summary.calc(pred.df.targ.collapse,merged.df.all.tmp)
    auc.plot.all[[2]]$dx.time = 'All'
    auc.plot.all[[2]]$var = 'All'
    
    auc.plot.all[[2]]$var.group = 'Time to Diagnosis'
    
    auc.plot.all.dx1 = summary.calc(pred.df.targ.collapse[pred.df.targ.collapse$GRP_Id %in% dx1$GRP_Id,],merged.df.all.tmp)
    auc.plot.all.dx1[[2]]$dx.time = '0-1'
    auc.plot.all.dx1[[2]]$var = '0-1'
    
    auc.plot.all.dx1[[2]]$var.group = 'Time to Diagnosis'
    
    
    auc.plot.all.dx2 = summary.calc(pred.df.targ.collapse[pred.df.targ.collapse$GRP_Id %in% dx2$GRP_Id,],merged.df.all.tmp)
    auc.plot.all.dx2[[2]]$dx.time = '1-2'
    auc.plot.all.dx2[[2]]$var = '1-2'
    
    auc.plot.all.dx2[[2]]$var.group = 'Time to Diagnosis'
    
    auc.plot.all.dx3 = summary.calc(pred.df.targ.collapse[pred.df.targ.collapse$GRP_Id %in% dx3$GRP_Id,],merged.df.all.tmp)
    auc.plot.all.dx3[[2]]$dx.time = '2-3'
    auc.plot.all.dx3[[2]]$var = '2-3'
    
    auc.plot.all.dx3[[2]]$var.group = 'Time to Diagnosis'
    
    auc.plot.all.dx4 = summary.calc(pred.df.targ.collapse[pred.df.targ.collapse$GRP_Id %in% dx4$GRP_Id,],merged.df.all.tmp)
    auc.plot.all.dx4[[2]]$dx.time = '3-4'
    auc.plot.all.dx4[[2]]$var = '3-4'
    auc.plot.all.dx4[[2]]$var.group = 'Time to Diagnosis'
    
    auc.plot.all.dx5 = summary.calc(pred.df.targ.collapse[pred.df.targ.collapse$GRP_Id %in% dx5$GRP_Id,],merged.df.all.tmp)
    auc.plot.all.dx5[[2]]$dx.time = '4-5'
    auc.plot.all.dx5[[2]]$var = '4-5'
    
    auc.plot.all.dx5[[2]]$var.group = 'Time to Diagnosis'
    
    if (dx.all == T){
      dx6 = merged.df.all.tmp[merged.df.all.tmp$Diagnosis_Time %in% c('Control','5-6'),]
      dx7 = merged.df.all.tmp[merged.df.all.tmp$Diagnosis_Time %in% c('Control','6-7','7-8','8-9','9+'),]
      dx8 = merged.df.all.tmp[merged.df.all.tmp$Diagnosis_Time %in% c('Control','7-8'),]
      dx9 = merged.df.all.tmp[merged.df.all.tmp$Diagnosis_Time %in% c('Control','8-9','9+'),]
      
      auc.plot.all.dx6 = summary.calc(pred.df.targ.collapse[pred.df.targ.collapse$GRP_Id %in% dx6$GRP_Id,],merged.df.all.tmp)
      auc.plot.all.dx6[[2]]$dx.time = '5-6'
      auc.plot.all.dx6[[2]]$var = '5-6'
      
      auc.plot.all.dx6[[2]]$var.group = 'Time to Diagnosis'
      
      auc.plot.all.dx7 = summary.calc(pred.df.targ.collapse[pred.df.targ.collapse$GRP_Id %in% dx7$GRP_Id,],merged.df.all.tmp)
      auc.plot.all.dx7[[2]]$dx.time = '6-7'
      auc.plot.all.dx7[[2]]$var = '6-7'
      
      auc.plot.all.dx7[[2]]$var.group = 'Time to Diagnosis'
      
      
      auc.plot.all.dx8 = summary.calc(pred.df.targ.collapse[pred.df.targ.collapse$GRP_Id %in% dx8$GRP_Id,],merged.df.all.tmp)
      auc.plot.all.dx8[[2]]$dx.time = '7-8'
      auc.plot.all.dx8[[2]]$var = '7-8'
      
      auc.plot.all.dx8[[2]]$var.group = 'Time to Diagnosis'
      
      
      auc.plot.all.dx9 = summary.calc(pred.df.targ.collapse[pred.df.targ.collapse$GRP_Id %in% dx9$GRP_Id,],merged.df.all.tmp)
      auc.plot.all.dx9[[2]]$dx.time = '8-10'
      auc.plot.all.dx9[[2]]$var = '8-10'
      
      auc.plot.all.dx9[[2]]$var.group = 'Time to Diagnosis'
      
      #
      
      auc.res2 = rbind(auc.plot.all[[2]][1,],
                       auc.plot.all.dx1[[2]][1,],
                       auc.plot.all.dx2[[2]][1,],
                       auc.plot.all.dx3[[2]][1,],
                       auc.plot.all.dx4[[2]][1,],
                       auc.plot.all.dx5[[2]][1,],
                       auc.plot.all.dx6[[2]][1,],
                       auc.plot.all.dx7[[2]][1,],
                       auc.plot.all.dx8[[2]][1,],
                       auc.plot.all.dx9[[2]][1,])
      auc.res2$title = 'Diagnosis Time'
      combined.auroc = rbind(combined.auroc,auc.res2)
      
      #annotating auc t
      auc.plot.all[[3]]$var = 'All'
      auc.plot.all.dx1[[3]]$var = '0-1'
      auc.plot.all.dx2[[3]]$var = '1-2'
      auc.plot.all.dx3[[3]]$var = '2-3'
      auc.plot.all.dx4[[3]]$var = '3-4'
      auc.plot.all.dx5[[3]]$var = '4-5'
      auc.plot.all.dx6[[3]]$var = '5-6'
      auc.plot.all.dx7[[3]]$var = '6+'
      auc.plot.all.dx8[[3]]$var = '7-8'
      auc.plot.all.dx9[[3]]$var = '8-10'
      
      
      auc.t =rbind( auc.plot.all[[3]],
                    auc.plot.all.dx1[[3]],
                    auc.plot.all.dx2[[3]],
                    auc.plot.all.dx3[[3]],
                    auc.plot.all.dx4[[3]],
                    auc.plot.all.dx5[[3]],
                    auc.plot.all.dx6[[3]],
                    auc.plot.all.dx7[[3]],
                    auc.plot.all.dx8[[3]],
                    auc.plot.all.dx9[[3]]
      )
      auc.t$title = 'Diagnosis Time'
      auc.t.combined = rbind(auc.t.combined,auc.t)
      
      #annotating others
      
      diagnosis_time_colors1 = c('#7A797C',"#048BA8",'#60D394','#AAF683','#FFD97D','#FF9B85','#C8553D','#F46197','#C3C3E6','#442B48')
      names(diagnosis_time_colors1) = c('Control','0-1','1-2','2-3','3-4','4-5','5-6','6-7','7-8','8-10')
      combined.auc.plot.all = rbind(auc.plot.all[[2]],
                                    auc.plot.all.dx1[[2]],
                                    auc.plot.all.dx2[[2]],
                                    auc.plot.all.dx3[[2]],
                                    auc.plot.all.dx4[[2]],
                                    auc.plot.all.dx5[[2]],
                                    auc.plot.all.dx6[[2]],
                                    auc.plot.all.dx7[[2]],
                                    auc.plot.all.dx8[[2]],
                                    auc.plot.all.dx9[[2]])
      combined.auc.plot.all$title = 'Diagnosis Time'
      
      plot1 = ggplot(combined.auc.plot.all,aes(x = fpr, y = tpr, col = dx.time)) + geom_line(linewidth=1) +
        scale_color_manual(values = c(diagnosis_time_colors1))+ #+ ggtitle(title) +
        theme_bw()+ 
        #facet_grid2(.~title)+
        theme(text = element_text(size=8),
              axis.text=element_text(size=8, face = "bold"),
              axis.title=element_text(size=8,face="bold"),
              legend.position = "none")+ guides(colour = guide_legend(override.aes = list(size=4),nrow=4)) + xlab('False Positive Rate') + ylab('Sensitivity')
      
      png(paste0(name,'.filler.',c,'.auroc.dxtime1.png'),height = 1100,width=1100,res = 300)
      print(plot1)
      dev.off()
      
      
    } else {
      
      auc.res2 = rbind(auc.plot.all[[2]][1,],
                       auc.plot.all.dx1[[2]][1,],
                       auc.plot.all.dx2[[2]][1,],
                       auc.plot.all.dx3[[2]][1,],
                       auc.plot.all.dx4[[2]][1,],
                       auc.plot.all.dx5[[2]][1,])
      auc.res2$title = 'Diagnosis Time'
      combined.auroc = rbind(combined.auroc,auc.res2)
      
      #annotating auc t
      auc.plot.all[[3]]$var = 'All'
      auc.plot.all.dx1[[3]]$var = '0-1'
      auc.plot.all.dx2[[3]]$var = '1-2'
      auc.plot.all.dx3[[3]]$var = '2-3'
      auc.plot.all.dx4[[3]]$var = '3-4'
      auc.plot.all.dx5[[3]]$var = '4-5'
      
      auc.t =rbind( auc.plot.all[[3]],
                    auc.plot.all.dx1[[3]],
                    auc.plot.all.dx2[[3]],
                    auc.plot.all.dx3[[3]],
                    auc.plot.all.dx4[[3]],
                    auc.plot.all.dx5[[3]]
      )
      auc.t$title = 'Diagnosis Time'
      auc.t.combined = rbind(auc.t.combined,auc.t)
      
      #annotating others
      diagnosis_time_colors1 = c('#7A797C',"#048BA8",'#60D394','#AAF683','#FFD97D','#FF9B85','#C8553D')
      names(diagnosis_time_colors1) = c('Control','0-1','1-2','2-3','3-4','4-5','5+')
      combined.auc.plot.all = rbind(auc.plot.all[[2]],
                                    auc.plot.all.dx1[[2]],
                                    auc.plot.all.dx2[[2]],
                                    auc.plot.all.dx3[[2]],
                                    auc.plot.all.dx4[[2]],
                                    auc.plot.all.dx5[[2]])
      combined.auc.plot.all$title = 'Diagnosis Time'
      plot1 = ggplot(combined.auc.plot.all,aes(x = fpr, y = tpr, col = dx.time)) + geom_line(linewidth=1) +
        scale_color_manual(values = c(diagnosis_time_colors1))+ #+ ggtitle(title) +
        theme_bw()+ 
        # facet_grid2(.~title)+
        theme(text = element_text(size=8),
              axis.text=element_text(size=8, face = "bold"),
              axis.title=element_text(size=8,face="bold"),
              legend.position = "none")+ guides(colour = guide_legend(override.aes = list(size=4),nrow=4)) + xlab('False Positive Rate') + ylab('Sensitivity')
      
      png(paste0(name,'.filler.',c,'.auroc.dxtime1.png'),height = 1100,width=1100,res = 300)
      print(plot1)
      dev.off()
    }
    
    if (dx.all == T){
      
      dx6 = merged.df.all.tmp[merged.df.all.tmp$Diagnosis_Time %in% c('Control','5-6','6-7','7-8','8-9','9+'),]
      
      
      auc.plot.all.dx6 = summary.calc(pred.df.targ.collapse[pred.df.targ.collapse$GRP_Id %in% dx6$GRP_Id,],merged.df.all.tmp)
      auc.plot.all.dx6[[2]]$dx.time = '5+'
      auc.plot.all.dx6[[2]]$var = '5+'
      
      auc.plot.all.dx6[[2]]$var.group = 'Time to Diagnosis'
      
      
      auc.res2 = rbind(auc.plot.all[[2]][1,],
                       auc.plot.all.dx1[[2]][1,],
                       auc.plot.all.dx2[[2]][1,],
                       auc.plot.all.dx3[[2]][1,],
                       auc.plot.all.dx4[[2]][1,],
                       auc.plot.all.dx5[[2]][1,],
                       auc.plot.all.dx6[[2]][1,])
      auc.res2$title = 'Diagnosis Time'
      combined.auroc = rbind(combined.auroc,auc.res2)
      
      #annotating auc t
      auc.plot.all[[3]]$var = 'All'
      auc.plot.all.dx1[[3]]$var = '0-1'
      auc.plot.all.dx2[[3]]$var = '1-2'
      auc.plot.all.dx3[[3]]$var = '2-3'
      auc.plot.all.dx4[[3]]$var = '3-4'
      auc.plot.all.dx5[[3]]$var = '4-5'
      auc.plot.all.dx6[[3]]$var = '5+'
      
      
      
      auc.t =rbind( auc.plot.all[[3]],
                    auc.plot.all.dx1[[3]],
                    auc.plot.all.dx2[[3]],
                    auc.plot.all.dx3[[3]],
                    auc.plot.all.dx4[[3]],
                    auc.plot.all.dx5[[3]],
                    auc.plot.all.dx6[[3]]
      )
      auc.t$title = 'Diagnosis Time'
      auc.t.combined = rbind(auc.t.combined,auc.t)
      
      #annotating others
      
      diagnosis_time_colors1 = c('#7A797C',"#048BA8",'#60D394','#AAF683','#FFD97D','#FF9B85','#C8553D')
      names(diagnosis_time_colors1) = c('Control','0-1','1-2','2-3','3-4','4-5','5+')
      combined.auc.plot.all = rbind(auc.plot.all[[2]],
                                    auc.plot.all.dx1[[2]],
                                    auc.plot.all.dx2[[2]],
                                    auc.plot.all.dx3[[2]],
                                    auc.plot.all.dx4[[2]],
                                    auc.plot.all.dx5[[2]],
                                    auc.plot.all.dx6[[2]])
      combined.auc.plot.all$title = 'Diagnosis Time'
      
      plot1 = ggplot(combined.auc.plot.all,aes(x = fpr, y = tpr, col = dx.time)) + geom_line(linewidth=1) +
        scale_color_manual(values = c(diagnosis_time_colors1))+ #+ ggtitle(title) +
        theme_bw()+ 
        #facet_grid2(.~title)+
        theme(text = element_text(size=8),
              axis.text=element_text(size=8, face = "bold"),
              axis.title=element_text(size=8,face="bold"),
              legend.position = "none")+ guides(colour = guide_legend(override.aes = list(size=4),nrow=4)) + xlab('False Positive Rate') + ylab('Sensitivity')
      
      png(paste0(name,'.filler.',c,'.auroc.dxtime1.png'),height = 1100,width=1100,res = 300)
      print(plot1)
      dev.off()
      
      
    } 
    
    
    
    
    pred.df.targ.collapse.annotated = merge(pred.df.targ.collapse, merged.df.all.tmp[,c('GRP_Id','Diagnosis_Time')],by='GRP_Id')
    
    if (dx.all == F){
      my_comparisons = list(c('Control','0-1'),
                            c('Control','1-2'),
                            c('Control','2-3'),
                            c('Control','3-4'),
                            c('Control','4-5'))
      options(scipen=2)
      pred.df.targ.collapse.annotated$Diagnosis.time.ks = factor(as.character(pred.df.targ.collapse.annotated$Diagnosis_Time),levels = c('0-1','1-2','2-3','3-4','4-5','Control'))
      print(kruskal.test(methylation_score ~ Diagnosis.time.ks,  data = pred.df.targ.collapse.annotated)  )
      
      
      plot1 = ggplot(pred.df.targ.collapse.annotated,aes(x = Diagnosis_Time, y = methylation_score, col = Diagnosis_Time)) + geom_boxplot(outlier.shape=NA)+geom_jitter(width=0.1)+
        scale_color_manual(values = c(diagnosis_time_colors1))+ #+ ggtitle(title) +
        theme_bw()+ 
        scale_y_continuous(limits=c(0,1.45),breaks = seq(0,1,0.25))+
        theme(text = element_text(size=8),
              axis.text=element_text(size=8, face = "bold"),
              axis.title=element_text(size=8,face="bold"),
              legend.position = "none")+ guides(colour = guide_legend(override.aes = list(size=4),nrow=4)) + xlab('Time to Diagnosis (Years)') + ylab('Methylation Score') +
        stat_compare_means(comparisons = my_comparisons,label.y=c(1,1.1,1.2,1.3,1.4,1.5),size=3)
      
      
      png(paste0(name,'.dx.methscore.',c,'.png'),height = 1100,width=1100,res = 300)
      print(plot1)
      dev.off()
      
    } else {
      my_comparisons = list(c('Control','0-1'),
                            c('Control','1-2'),
                            c('Control','2-3'),
                            c('Control','3-4'),
                            c('Control','4-5'),
                            c('Control','5-6'),
                            c('Control','6-7'),
                            c('Control','7-8'))
      diagnosis_time_colors1 = c('grey',"#048BA8",'#60D394','#AAF683','#FFD97D','#FF9B85','#C8553D')#,'#F46197','#C3C3E6','#442B48'
      names(diagnosis_time_colors1) = c('Control','0-1','1-2','2-3','3-4','4-5','5+')#,'6-7','7-8','8-10')
      
      options(scipen=2)
      pred.df.targ.collapse.annotated$Diagnosis.time.ks = factor(as.character(pred.df.targ.collapse.annotated$Diagnosis_Time),levels = c('0-1','1-2','2-3','3-4','4-5','5-6','6-7','7-8','8-10','Control'))
      print(kruskal.test(methylation_score ~ Diagnosis.time.ks,  data = pred.df.targ.collapse.annotated)  )
      
      plot1 = ggplot(pred.df.targ.collapse.annotated,aes(x = Diagnosis_Time, y = methylation_score, col = Diagnosis_Time)) + geom_boxplot(outlier.shape=NA)+geom_jitter(width=0.1)+
        scale_color_manual(values = c(diagnosis_time_colors1))+ #+ ggtitle(title) +
        theme_bw()+ 
        scale_y_continuous(limits=c(0,1.8),breaks = seq(0,1.85,0.25))+
        theme(text = element_text(size=8),
              axis.text=element_text(size=8, face = "bold"),
              axis.title=element_text(size=8,face="bold"),
              legend.position = "none")+ guides(colour = guide_legend(override.aes = list(size=4),nrow=4)) + xlab('Time to Diagnosis (Years)') + ylab('Methylation Score') +
        stat_compare_means(comparisons = my_comparisons,label.y=seq(1,1.85,0.1),size=3)
      
      
      png(paste0(name,'.dx.methscore.',c,'.png'),height = 1100,width=1100,res = 300)
      print(plot1)
      dev.off()
      
    }
    
    
    ####plotting stage####
    print('stage')
    stage1 = merged.df.all.tmp[merged.df.all.tmp$Stage %in% ('I') | merged.df.all.tmp$group =='Control',]
    stage2 = merged.df.all.tmp[merged.df.all.tmp$Stage %in% c('II')| merged.df.all.tmp$group =='Control',]
    stage34 = merged.df.all.tmp[merged.df.all.tmp$Stage %in% c('III','IV')| merged.df.all.tmp$group =='Control',]
    auc.plot.all = summary.calc(pred.df.targ.collapse,merged.df.all.tmp)
    auc.plot.all[[2]]$stage = 'All'
    auc.plot.all[[2]]$var = 'All'
    auc.plot.all[[2]]$var.group = 'Stage'
    
    stagereported = merged.df.all.tmp[merged.df.all.tmp$Stage %in% c('I','II','III','IV')| merged.df.all.tmp$group =='Control',]
    auc.plot.all.sr = summary.calc(pred.df.targ.collapse[pred.df.targ.collapse$GRP_Id %in% stagereported$GRP_Id,],merged.df.all.tmp)
    auc.plot.all.sr[[2]]$stage = 'I-IV'
    auc.plot.all.sr[[2]]$var.group = 'Stage'
    auc.plot.all.sr[[2]]$var ='I-IV'
    
    
    
    auc.plot.all.s1 = summary.calc(pred.df.targ.collapse[pred.df.targ.collapse$GRP_Id %in% stage1$GRP_Id,],merged.df.all.tmp)
    auc.plot.all.s1[[2]]$stage = 'I'
    auc.plot.all.s1[[2]]$var.group = 'Stage'
    auc.plot.all.s1[[2]]$var ='I'
    
    auc.plot.all.s2 = summary.calc(pred.df.targ.collapse[pred.df.targ.collapse$GRP_Id %in% stage2$GRP_Id,],merged.df.all.tmp)
    auc.plot.all.s2[[2]]$stage = 'II'
    auc.plot.all.s2[[2]]$var =  'II'
    
    auc.plot.all.s2[[2]]$var.group = 'Stage'
    
    
    auc.plot.all.s34 = summary.calc(pred.df.targ.collapse[pred.df.targ.collapse$GRP_Id %in% stage34$GRP_Id,],merged.df.all.tmp)
    auc.plot.all.s34[[2]]$stage = 'III/IV'
    auc.plot.all.s34[[2]]$var = 'III/IV'
    
    auc.plot.all.s34[[2]]$var.group = 'Stage'
    
    #annotating auc t
    auc.plot.all[[3]]$var = 'All'
    auc.plot.all.s1[[3]]$var = 'I'
    auc.plot.all.s2[[3]]$var = 'II'
    auc.plot.all.s34[[3]]$var = 'III/IV'
    
    auc.t =rbind(# auc.plot.all[[3]],
      auc.plot.all.s1[[3]],
      auc.plot.all.s2[[3]],
      auc.plot.all.s34[[3]]
    )
    auc.t$title = 'Stage'
    auc.t.combined = rbind(auc.t.combined,auc.t)
    
    #auc t
    
    
    
    #
    auc.res1 = rbind(
      auc.plot.all.s1[[2]][1,],
      auc.plot.all.s2[[2]][1,],
      auc.plot.all.s34[[2]][1,])
    auc.res1$title = 'Stage'
    
    combined.auroc = combined.auroc[,colnames(combined.auroc) %in% colnames(auc.res1)]
    
    combined.auroc = rbind(combined.auroc,auc.res1[,colnames(auc.res1) %in% colnames(combined.auroc)])
    
    stage_colors = c('All' = "#000004FF", 'I' = "#3B0F70FF",'II' = "#8C2981FF",'III/IV' = "#DE4968FF")
    
    grade_colors = c('1' = '#FEF6C9','2' = '#D4DFC7', '3' = '#96C0B7', '9' ='#878E88')
    
    
    plot1 = ggplot(auc.t[auc.t$time >= 365*1,],aes(x = time/365, y = AUC, col = var)) + geom_line(linewidth=1) +
      scale_color_manual(values = c(stage_colors))+ #+ ggtitle(title) +
      scale_y_continuous(limits=c(0,1))+
      theme_bw()+ 
      #facet_grid2(.~title)+
      theme(text = element_text(size=8),
            axis.text=element_text(size=8, face = "bold"),
            axis.title=element_text(size=8,face="bold"),
            legend.position = "none")+ guides(colour = guide_legend(override.aes = list(size=4),nrow=4)) +
      xlab('Time (Years)') + ylab('Time-Dependent CV AUROC')
    
    png(paste0(name,'.filler.',c,'.timeauroc.stage.png'),height = 1100,width=1100,res = 300)
    print(plot1)
    dev.off()
    
    combined.auc.plot.all = rbind(
      auc.plot.all.s1[[2]],
      auc.plot.all.s2[[2]],
      auc.plot.all.s34[[2]])
    plot1 = ggplot(combined.auc.plot.all,aes(x = fpr, y = tpr, col = stage)) + geom_line(linewidth=1) +
      scale_color_manual(values = c(stage_colors))+ #+ ggtitle(title) +
      theme_bw()+ 
      theme(text = element_text(size=8),
            axis.text=element_text(size=8, face = "bold"),
            axis.title=element_text(size=8,face="bold"),
            legend.position = "none")+ guides(colour = guide_legend(override.aes = list(size=4),nrow=4)) + xlab('False Positive Rate') + ylab('Sensitivity')
    
    png(paste0(name,'.filler.',c,'.auroc.stage.png'),height = 1100,width=1100,res = 300)
    print(plot1)
    dev.off()
    
    plot1 = ggplot(combined.auc.plot.all,aes(x = fpr, y = tpr, col = stage)) + geom_line(linewidth=1) +
      scale_color_manual(values = c(stage_colors))+ #+ ggtitle(title) +
      theme_bw()+ 
      theme(text = element_text(size=8),
            axis.text=element_text(size=8, face = "bold"),
            axis.title=element_text(size=8,face="bold"),
            legend.position = "right")+ guides(colour = guide_legend(override.aes = list(size=4),nrow=4)) + xlab('False Positive Rate') + ylab('Sensitivity')
    
    png(paste0(name,'.filler.',c,'.auroc.stage.labs.png'),height = 1000,width=1000,res = 300)
    print(plot1)
    dev.off()
    
    stage_colors.sample = c('Not Reported' = '#66462C','Control'= '#7A797C','0' = "#000004FF", 'I' = "#3B0F70FF",'II' = "#8C2981FF",'III' = "#DE4968FF",'IV' = "#FE9F6DFF")
    
    
    
    plot1 = ggplot(auc.res1,aes(x = stage, y = auc, col = stage)) +
      geom_errorbar(aes(ymin=auc.lower, ymax=auc.upper), width=.2,
                    position=position_dodge(0.05))+
      geom_point()+
      scale_y_continuous(limits=c(0,1))+
      scale_color_manual(values = c(stage_colors))+ #+ ggtitle(title) +
      theme_bw()+ 
      #facet_grid(.~title)+
      
      theme(text = element_text(size=8),
            axis.text=element_text(size=8, face = "bold"),
            axis.title=element_text(size=8,face="bold"),
            legend.position = "none")+ 
      guides(colour = guide_legend(override.aes = list(size=4),nrow=4)) +
      xlab('Stage at Diagnosis') + ylab('CV AUROC (95% CI)')
    
    png(paste0(name,'.filler.',c,'.auroc.ci.stage.png'),height = 1100,width=1100,res = 300)
    print(plot1)
    dev.off()
    plot1 = ggplot(auc.res1,aes(x = stage, y = ci, col = stage)) +
      geom_errorbar(aes(ymin=ci.lower, ymax=ci.upper), width=.2,
                    position=position_dodge(0.05))+
      geom_point()+
      scale_y_continuous(limits=c(0,1))+
      scale_color_manual(values = c(stage_colors))+ #+ ggtitle(title) +
      theme_bw()+ 
      #facet_grid(.~title)+
      theme(text = element_text(size=8),
            axis.text=element_text(size=8, face = "bold"),
            axis.title=element_text(size=8,face="bold"),
            legend.position = "none")+ 
      guides(colour = guide_legend(override.aes = list(size=4),nrow=4)) +
      xlab('Stage at Diagnosis') + ylab('CV Concordance Index (95% CI)')
    
    png(paste0(name,'.filler.',c,'.ci.ci.stage.png'),height = 1100,width=1100,res = 300)
    print(plot1)
    dev.off()
    plot1 = ggplot(auc.res1,aes(x = stage, y =se.spec.95 , col = stage)) +
      geom_errorbar(aes(ymin=se.spec.95.lower, ymax=se.spec.95.upper), width=.2,
                    position=position_dodge(0.05))+
      geom_point()+
      scale_y_continuous(limits=c(0,1))+
      scale_color_manual(values = c(stage_colors))+ #+ ggtitle(title) +
      theme_bw()+ 
      #facet_grid(.~title)+
      theme(text = element_text(size=8),
            axis.text=element_text(size=8, face = "bold"),
            axis.title=element_text(size=8,face="bold"),
            legend.position = "none")+ 
      guides(colour = guide_legend(override.aes = list(size=4),nrow=4)) +
      xlab('Stage at Diagnosis') + ylab('Sensitivity at 95% Specificity (95% CI)')
    
    png(paste0(name,'.filler.',c,'.sens.spec95.ci.stage.png'),height = 1100,width=1100,res = 300)
    print(plot1)
    dev.off()
    
    plot1 = ggplot(auc.res1,aes(x = stage, y = se.spec.90, col = stage)) +
      geom_errorbar(aes(ymin=se.spec.90.lower, ymax=se.spec.90.upper), width=.2,
                    position=position_dodge(0.05))+
      geom_point()+
      scale_y_continuous(limits=c(0,1))+
      scale_color_manual(values = c(stage_colors))+ #+ ggtitle(title) +
      theme_bw()+ 
      #facet_grid(.~title)+
      theme(text = element_text(size=8),
            axis.text=element_text(size=8, face = "bold"),
            axis.title=element_text(size=8,face="bold"),
            legend.position = "none")+ 
      guides(colour = guide_legend(override.aes = list(size=4),nrow=4)) +
      xlab('Stage at Diagnosis')+ ylab('Sensitivity at 90% Specificity (95% CI)')
    
    png(paste0(name,'.filler.',c,'.sens.spec90.ci.stage.png'),height = 1100,width=1100,res = 300)
    print(plot1)
    dev.off()
    
    plot1 = ggplot(combined.auc.plot.all,aes(x = fpr, y = tpr, col = stage)) + geom_line(linewidth=1) +
      scale_color_manual(values = c(stage_colors.sample))+ #+ ggtitle(title) +
      theme_bw()+ 
      theme(text = element_text(size=8),
            axis.text=element_text(size=8, face = "bold"),
            axis.title=element_text(size=8,face="bold"),
            legend.position = "right")+ guides(colour = guide_legend(override.aes = list(size=4),ncol=1)) + xlab('False Positive Rate') + ylab('Sensitivity')
    
    png(paste0(name,'.filler.',c,'.auroc.stage.labs.png'),height = 1100,width=1100,res = 300)
    print(plot1)
    dev.off()
    
    
    #boxplot of score
    pred.df.targ.collapse.annotated = merge(pred.df.targ.collapse, merged.df.all.tmp[,c('GRP_Id','Diagnosis_Time','filler','Stage','GRADE_CD')],by='GRP_Id')
    stage_colors.score = c('Control' = "grey", 'I' = "#3B0F70FF",'II' = "#8C2981FF",'III' = "#DE4968FF",'IV' = '#FE9F6DFF','NR' = '#66462C')
    
    pred.df.targ.collapse.annotated$Stage = ifelse(pred.df.targ.collapse.annotated$Stage %in% c('Not Reported','Not Reported/Unknown'),'NR',pred.df.targ.collapse.annotated$Stage)
    pred.df.targ.collapse.annotated$Stage = ifelse(pred.df.targ.collapse.annotated$reported == 'Control','Control',pred.df.targ.collapse.annotated$Stage)
    pred.df.targ.collapse.annotated$Stage = factor(pred.df.targ.collapse.annotated$Stage,levels = c('Control','I','II','III','IV','NR'))
    
    my_comparisons= list()
    if (wilcox.test(methylation_score ~ Stage, pred.df.targ.collapse.annotated[pred.df.targ.collapse.annotated$Stage %in% c('Control','I'),])$p.value < 0.05 ) {
      my_comparisons[[length(my_comparisons)+1]] = c('Control','I')
    }
    
    if (wilcox.test(methylation_score ~ Stage, pred.df.targ.collapse.annotated[pred.df.targ.collapse.annotated$Stage %in% c('Control','II'),])$p.value < 0.05 ) {
      my_comparisons[[length(my_comparisons)+1]] = c('Control','II')
    }
    
    if (wilcox.test(methylation_score ~ Stage, pred.df.targ.collapse.annotated[pred.df.targ.collapse.annotated$Stage %in% c('Control','III'),])$p.value < 0.05 ) {
      my_comparisons[[length(my_comparisons)+1]] = c('Control','III')
    }
    
    if (wilcox.test(methylation_score ~ Stage, pred.df.targ.collapse.annotated[pred.df.targ.collapse.annotated$Stage %in% c('Control','IV'),])$p.value < 0.05 ) {
      my_comparisons[[length(my_comparisons)+1]] = c('Control','IV')
    }
    
    
    
    
    print(kruskal.test(methylation_score ~ Stage,  data = pred.df.targ.collapse.annotated[as.character(pred.df.targ.collapse.annotated$Stage) %in%  c('Control','I','II','III','IV'),])  )
    plot1 = ggplot(pred.df.targ.collapse.annotated,aes(x = Stage, y = methylation_score, col = Stage)) +
      geom_boxplot(outlier.shape=NA)+geom_jitter(width=0.1)+
      scale_color_manual(values = c(stage_colors.score))+ #+ ggtitle(title) +
      theme_bw()+ 
      scale_y_continuous(limits=c(0,1+0.1*length(my_comparisons)))+
      
      theme(text = element_text(size=8),
            axis.text=element_text(size=7, face = "bold"),
            axis.title=element_text(size=8,face="bold"),
            legend.position = "none")+ 
      guides(colour = guide_legend(override.aes = list(size=4),nrow=4)) +
      xlab('Stage at Diagnosis') +
      ylab('Methylation Score')+
      stat_compare_means(comparisons = my_comparisons,label.y=seq(1,1+length(my_comparisons)*0.1,0.1),size=3)
    
    
    png(paste0(name,'.stage.methscore.',c,'.png'),height = 1100,width=1000,res = 300)
    print(plot1)
    dev.off()
    #copy
    plot1 = ggplot(pred.df.targ.collapse.annotated,aes(x = Stage, y = methylation_score, col = Stage)) +
      geom_boxplot(outlier.shape=NA)+geom_jitter(width=0.1)+
      scale_color_manual(values = c(stage_colors.sample))+ #+ ggtitle(title) +
      theme_bw()+ 
      scale_y_continuous(limits=c(0,1))+
      
      theme(text = element_text(size=8),
            axis.text=element_text(size=8, face = "bold"),
            axis.title=element_text(size=8,face="bold"),
            legend.position = "none")+ guides(colour = guide_legend(override.aes = list(size=4),nrow=4)) + xlab('Time to Diagnosis') + ylab('Methylation Score')
    
    png(paste0(name,'.stage.methscore.',c,'.png'),height = 900,width=1100,res = 300)
    print(plot1)
    dev.off()
    #copy
    
    #####morphology#####
    {
      print('morphology')
      brca.morphology = c('85003' = 'Infiltrating ductal carcinoma',
                          '85203' = 'Lobular NOS',
                          '85223'='Infiltrating ductal and lobular carcinoma',
                          '85233' = 'Infiltrating ductal mixed with other',
                          '85433' = 'Paget disease + infiltrating ductal carcinoma',
                          '85073' = 'Micropapillary carcinoma',
                          '84803' = 'Mucinous adenocarcinoma',
                          '82603' = 'Papillary carcinoma NOS',
                          '82113' = 'Tubular adenocarcinoma',
                          '80103' = 'Carcinoma NOS',
                          'Control' = 'Control'
      )
      
      pred.df.targ.collapse.annotated = merge(pred.df.targ.collapse, merged.df.all.tmp[,c('GRP_Id','Cancer','Diagnosis_Time','Sex','SDC_GENDER')],by='GRP_Id')
      pred.df.targ.collapse.annotated = merge(pred.df.targ.collapse.annotated, merged.df.all.tmp[,c('GRP_Id','CURR_MORPH_CD')])
      pred.df.targ.collapse.annotated$CURR_MORPH_CD = ifelse(pred.df.targ.collapse.annotated$reported =='Control','Control',pred.df.targ.collapse.annotated$CURR_MORPH_CD)
      pred.df.targ.collapse.annotated$morphology = brca.morphology[pred.df.targ.collapse.annotated$CURR_MORPH_CD]
      
      morph1 = pred.df.targ.collapse.annotated[pred.df.targ.collapse.annotated$morphology %in% ('Infiltrating ductal and lobular carcinoma') | pred.df.targ.collapse.annotated$reported =='Control',]
      morph2 = pred.df.targ.collapse.annotated[pred.df.targ.collapse.annotated$morphology %in% c('Infiltrating ductal carcinoma')| pred.df.targ.collapse.annotated$reported =='Control',]
      morph3 = pred.df.targ.collapse.annotated[pred.df.targ.collapse.annotated$morphology %in% c('Lobular NOS')| pred.df.targ.collapse.annotated$reported =='Control',]
      
      
      auc.plot.all.s1 = summary.calc(pred.df.targ.collapse[pred.df.targ.collapse$GRP_Id %in% morph1$GRP_Id,],merged.df.all.tmp)
      auc.plot.all.s1[[2]]$morph = 'Ductal+Lobular'
      auc.plot.all.s1[[2]]$var.group = 'Morphology'
      auc.plot.all.s1[[2]]$var = 'Ductal+Lobular'
      
      auc.plot.all.s2 = summary.calc(pred.df.targ.collapse[pred.df.targ.collapse$GRP_Id %in% morph2$GRP_Id,],merged.df.all.tmp)
      auc.plot.all.s2[[2]]$morph = 'Ductal'
      auc.plot.all.s2[[2]]$var =  'Ductal'
      
      auc.plot.all.s2[[2]]$var.group = 'Morphology'
      
      
      auc.plot.all.s34 = summary.calc(pred.df.targ.collapse[pred.df.targ.collapse$GRP_Id %in% morph3$GRP_Id,],merged.df.all.tmp)
      auc.plot.all.s34[[2]]$morph = 'Lobular'
      auc.plot.all.s34[[2]]$var = 'Lobular'
      
      auc.plot.all.s34[[2]]$var.group = 'Morphology'
      
      
      #annotating auc t
      auc.plot.all.s1[[3]]$var = 'Ductal+Lobular'
      auc.plot.all.s2[[3]]$var ='Ductal'
      auc.plot.all.s34[[3]]$var = 'Lobular'
      auc.t =rbind( 
        auc.plot.all.s1[[3]],
        auc.plot.all.s2[[3]],
        auc.plot.all.s34[[3]]
      )
      auc.t$title = 'Morphology'
      auc.t.combined = rbind(auc.t.combined,auc.t)
      
      #auc t
      
      #
      auc.res1 = rbind(
        auc.plot.all.s1[[2]][1,],
        auc.plot.all.s2[[2]][1,],
        auc.plot.all.s34[[2]][1,])
      auc.res1$title = 'morph'
      
      combined.auroc = combined.auroc[,colnames(combined.auroc) %in% colnames(auc.res1)]
      
      combined.auroc = rbind(combined.auroc,auc.res1[,colnames(auc.res1) %in% colnames(combined.auroc)])
      
      morph_colors = c('Ductal+Lobular' = "#AFA2FF",'Ductal' = "#EF8275",'Lobular' = "#8A4F7D",'All' = 'black')
      
      grade_colors = c('1' = '#FEF6C9','2' = '#D4DFC7', '3' = '#96C0B7', '9' ='#878E88')
      
      
      plot1 = ggplot(auc.t[auc.t$time >= 365*1,],aes(x = time/365, y = AUC, col = var)) + geom_line(linewidth=1) +
        scale_color_manual(values = c(morph_colors))+ #+ ggtitle(title) +
        scale_y_continuous(limits=c(0,1))+
        theme_bw()+ 
        #facet_grid2(.~title)+
        theme(text = element_text(size=8),
              axis.text=element_text(size=8, face = "bold"),
              axis.title=element_text(size=8,face="bold"),
              legend.position = "none")+ guides(colour = guide_legend(override.aes = list(size=4),nrow=4)) +
        xlab('Time (Years)') + ylab('Time-Dependent CV AUROC')
      
      png(paste0(name,'.filler.',c,'.timeauroc.morph.png'),height = 1100,width=1100,res = 300)
      print(plot1)
      dev.off()
      
      combined.auc.plot.all = rbind(
        auc.plot.all.s1[[2]],
        auc.plot.all.s2[[2]],
        auc.plot.all.s34[[2]])
      plot1 = ggplot(combined.auc.plot.all,aes(x = fpr, y = tpr, col = morph)) + geom_line(linewidth=1) +
        scale_color_manual(values = c(morph_colors))+ #+ ggtitle(title) +
        theme_bw()+ 
        theme(text = element_text(size=8),
              axis.text=element_text(size=8, face = "bold"),
              axis.title=element_text(size=8,face="bold"),
              legend.position = "none")+ guides(colour = guide_legend(override.aes = list(size=4),nrow=4)) + xlab('False Positive Rate') + ylab('Sensitivity')
      
      png(paste0(name,'.filler.',c,'.auroc.morph.png'),height = 1100,width=1100,res = 300)
      print(plot1)
      dev.off()
      
      plot1 = ggplot(combined.auc.plot.all,aes(x = fpr, y = tpr, col = morph)) + geom_line(linewidth=1) +
        scale_color_manual(values = c(morph_colors))+ #+ ggtitle(title) +
        theme_bw()+ 
        theme(text = element_text(size=8),
              axis.text=element_text(size=8, face = "bold"),
              axis.title=element_text(size=8,face="bold"),
              legend.position = "right")+ guides(colour = guide_legend(override.aes = list(size=4),nrow=4)) + xlab('False Positive Rate') + ylab('Sensitivity')
      
      png(paste0(name,'.filler.',c,'.auroc.morph.labs.png'),height = 1000,width=1000,res = 300)
      print(plot1)
      dev.off()
      
      morph_colors.sample = c('Not Reported' = '#66462C','Control'= '#7A797C','0' = "#000004FF", 'I' = "#3B0F70FF",'II' = "#8C2981FF",'III' = "#DE4968FF",'IV' = "#FE9F6DFF")
      morph_colors = c('Ductal+Lobular' = "#AFA2FF",'Ductal' = "#EF8275",'Lobular' = "#8A4F7D")
      
      
      
      plot1 = ggplot(auc.res1,aes(x = morph, y = auc, col = morph)) +
        geom_errorbar(aes(ymin=auc.lower, ymax=auc.upper), width=.2,
                      position=position_dodge(0.05))+
        geom_point()+
        scale_y_continuous(limits=c(0,1))+
        scale_color_manual(values = c(morph_colors))+ #+ ggtitle(title) +
        theme_bw()+ 
        #facet_grid(.~title)+
        
        theme(text = element_text(size=8),
              axis.text=element_text(size=8, face = "bold"),
              axis.title=element_text(size=8,face="bold"),
              legend.position = "none")+ 
        guides(colour = guide_legend(override.aes = list(size=4),nrow=4)) +
        xlab('morph at Diagnosis') + ylab('CV AUROC (95% CI)')
      
      png(paste0(name,'.filler.',c,'.auroc.ci.morph.png'),height = 1100,width=1100,res = 300)
      print(plot1)
      dev.off()
      plot1 = ggplot(auc.res1,aes(x = morph, y = ci, col = morph)) +
        geom_errorbar(aes(ymin=ci.lower, ymax=ci.upper), width=.2,
                      position=position_dodge(0.05))+
        geom_point()+
        scale_y_continuous(limits=c(0,1))+
        scale_color_manual(values = c(morph_colors))+ #+ ggtitle(title) +
        theme_bw()+ 
        #facet_grid(.~title)+
        theme(text = element_text(size=8),
              axis.text=element_text(size=8, face = "bold"),
              axis.title=element_text(size=8,face="bold"),
              legend.position = "none")+ 
        guides(colour = guide_legend(override.aes = list(size=4),nrow=4)) +
        xlab('morph at Diagnosis') + ylab('CV Concordance Index (95% CI)')
      
      png(paste0(name,'.filler.',c,'.ci.ci.morph.png'),height = 1100,width=1100,res = 300)
      print(plot1)
      dev.off()
      plot1 = ggplot(auc.res1,aes(x = morph, y =se.spec.95 , col = morph)) +
        geom_errorbar(aes(ymin=se.spec.95.lower, ymax=se.spec.95.upper), width=.2,
                      position=position_dodge(0.05))+
        geom_point()+
        scale_y_continuous(limits=c(0,1))+
        scale_color_manual(values = c(morph_colors))+ #+ ggtitle(title) +
        theme_bw()+ 
        #facet_grid(.~title)+
        theme(text = element_text(size=8),
              axis.text=element_text(size=8, face = "bold"),
              axis.title=element_text(size=8,face="bold"),
              legend.position = "none")+ 
        guides(colour = guide_legend(override.aes = list(size=4),nrow=4)) +
        xlab('morph at Diagnosis') + ylab('Sensitivity at 95% Specificity (95% CI)')
      
      png(paste0(name,'.filler.',c,'.sens.spec95.ci.morph.png'),height = 1100,width=1100,res = 300)
      print(plot1)
      dev.off()
      
      plot1 = ggplot(auc.res1,aes(x = morph, y = se.spec.90, col = morph)) +
        geom_errorbar(aes(ymin=se.spec.90.lower, ymax=se.spec.90.upper), width=.2,
                      position=position_dodge(0.05))+
        geom_point()+
        scale_y_continuous(limits=c(0,1))+
        scale_color_manual(values = c(morph_colors))+ #+ ggtitle(title) +
        theme_bw()+ 
        #facet_grid(.~title)+
        theme(text = element_text(size=8),
              axis.text=element_text(size=8, face = "bold"),
              axis.title=element_text(size=8,face="bold"),
              legend.position = "none")+ 
        guides(colour = guide_legend(override.aes = list(size=4),nrow=4)) +
        xlab('morph at Diagnosis')+ ylab('Sensitivity at 90% Specificity (95% CI)')
      
      png(paste0(name,'.filler.',c,'.sens.spec90.ci.morph.png'),height = 1100,width=1100,res = 300)
      print(plot1)
      dev.off()
      
      plot1 = ggplot(combined.auc.plot.all,aes(x = fpr, y = tpr, col = morph)) + geom_line(linewidth=1) +
        scale_color_manual(values = c(morph_colors.sample))+ #+ ggtitle(title) +
        theme_bw()+ 
        theme(text = element_text(size=8),
              axis.text=element_text(size=8, face = "bold"),
              axis.title=element_text(size=8,face="bold"),
              legend.position = "right")+ guides(colour = guide_legend(override.aes = list(size=4),ncol=1)) + xlab('False Positive Rate') + ylab('Sensitivity')
      
      png(paste0(name,'.filler.',c,'.auroc.morph.labs.png'),height = 1100,width=1100,res = 300)
      print(plot1)
      dev.off()
      
      
      #boxplot of score
      morph_colors.score = c('Control' = "grey", 'I' = "#3B0F70FF",'II' = "#8C2981FF",'III' = "#DE4968FF",'IV' = '#FE9F6DFF','NR' = '#66462C')
      morph_colors.score = c('Control' = "grey",'Ductal+Lobular' = "#AFA2FF",'Ductal' = "#EF8275",'Lobular' = "#8A4F7D")
      
      
      my_comparisons= list()
      pred.df.targ.collapse.annotated$morph= ifelse(pred.df.targ.collapse.annotated$morphology == 'Lobular NOS','Lobular',
                                                    ifelse(pred.df.targ.collapse.annotated$morphology == 'Infiltrating ductal carcinoma','Ductal',
                                                           ifelse(pred.df.targ.collapse.annotated$morphology == 'Infiltrating ductal and lobular carcinoma','Ductal+Lobular',pred.df.targ.collapse.annotated$morphology)))
      if (wilcox.test(methylation_score ~ morph, pred.df.targ.collapse.annotated[pred.df.targ.collapse.annotated$morph %in% c('Control','Ductal+Lobular'),])$p.value < 0.05 ) {
        my_comparisons[[length(my_comparisons)+1]] = c('Control','Ductal+Lobular')
      }
      
      if (wilcox.test(methylation_score ~ morph, pred.df.targ.collapse.annotated[pred.df.targ.collapse.annotated$morph %in% c('Control','Ductal'),])$p.value < 0.05 ) {
        my_comparisons[[length(my_comparisons)+1]] = c('Control','Ductal')
      }
      
      if (wilcox.test(methylation_score ~ morph, pred.df.targ.collapse.annotated[pred.df.targ.collapse.annotated$morph %in% c('Control','Lobular'),])$p.value < 0.05 ) {
        my_comparisons[[length(my_comparisons)+1]] = c('Control','Lobular')
      }
      
      morph_colors.score = c('Control' = "grey",'Ductal+Lobular' = "#AFA2FF",'Ductal' = "#EF8275",'Lobular' = "#8A4F7D")
      
      
      print(kruskal.test(methylation_score ~ morph,  data = pred.df.targ.collapse.annotated[as.character(pred.df.targ.collapse.annotated$morph) %in%  names(morph_colors.score),])  )
      plot1 = ggplot(pred.df.targ.collapse.annotated,aes(x = morph, y = methylation_score, col = morph)) +
        geom_boxplot(outlier.shape=NA)+geom_jitter(width=0.1)+
        scale_color_manual(values = c(morph_colors.score))+ #+ ggtitle(title) +
        theme_bw()+ 
        scale_y_continuous(limits=c(0,1+0.1*length(my_comparisons)))+
        
        theme(text = element_text(size=8),
              axis.text=element_text(size=7, face = "bold"),
              axis.title=element_text(size=8,face="bold"),
              legend.position = "none")+ 
        guides(colour = guide_legend(override.aes = list(size=4),nrow=4)) +
        xlab('morph at Diagnosis') +
        ylab('Methylation Score')+
        stat_compare_means(comparisons = my_comparisons,label.y=seq(1,1+length(my_comparisons)*0.1,0.1),size=3)
      
      
      png(paste0(name,'.morph.methscore.',c,'.png'),height = 1100,width=1000,res = 300)
      print(plot1)
      dev.off()
      
    }
    
    
    
    
    #plotting grade
    print('grade')
    
    grade1 = merged.df.all.tmp[merged.df.all.tmp$GRADE %in% ('1') | merged.df.all.tmp$group =='Control',]
    grade2 = merged.df.all.tmp[merged.df.all.tmp$GRADE %in% c('2')| merged.df.all.tmp$group =='Control',]
    grade3 = merged.df.all.tmp[merged.df.all.tmp$GRADE %in% c('3')| merged.df.all.tmp$group =='Control',]
    grade9 = merged.df.all.tmp[merged.df.all.tmp$GRADE %in% c('9')| merged.df.all.tmp$group =='Control',]
    
    
    auc.plot.all.g1 = summary.calc(pred.df.targ.collapse[pred.df.targ.collapse$GRP_Id %in% grade1$GRP_Id,],merged.df.all.tmp)
    auc.plot.all.g1[[2]]$grade = '1'
    auc.plot.all.g1[[2]]$var.group = 'Grade'
    auc.plot.all.g1[[2]]$var = '1'
    
    
    auc.plot.all$auc = auc_calc(pred.df.targ.collapse,c('Control','Cancer'))
    auc.plot.all.g2 = summary.calc(pred.df.targ.collapse[pred.df.targ.collapse$GRP_Id %in% grade2$GRP_Id,],merged.df.all.tmp)
    auc.plot.all.g2[[2]]$grade = '2'
    auc.plot.all.g2[[2]]$var.group = 'Grade'
    auc.plot.all.g2[[2]]$var = '2'
    
    
    
    auc.plot.all$auc = auc_calc(pred.df.targ.collapse,c('Control','Cancer'))
    auc.plot.all.g3 = summary.calc(pred.df.targ.collapse[pred.df.targ.collapse$GRP_Id %in% grade3$GRP_Id,],merged.df.all.tmp)
    auc.plot.all.g3[[2]]$grade = '3'
    auc.plot.all.g3[[2]]$var.group = 'Grade'
    auc.plot.all.g3[[2]]$var = '3'
    
    auc.plot.all$auc = auc_calc(pred.df.targ.collapse,c('Control','Cancer'))
    auc.plot.all.g9 = summary.calc(pred.df.targ.collapse[pred.df.targ.collapse$GRP_Id %in% grade9$GRP_Id,],merged.df.all.tmp)
    auc.plot.all.g9[[2]]$grade = '9'
    auc.plot.all.g9[[2]]$var.group = 'Grade'
    auc.plot.all.g9[[2]]$var = '9'
    
    
    auc.plot.all.g1[[3]]$var = '1'
    auc.plot.all.g2[[3]]$var = '2'
    auc.plot.all.g3[[3]]$var = '3'
    auc.plot.all.g9[[3]]$var = '9'
    
    
    auc.t =rbind( auc.plot.all.g1[[3]],
                  auc.plot.all.g2[[3]],
                  auc.plot.all.g3[[3]],
                  auc.plot.all.g9[[3]]
                  
    )
    auc.t$title = 'Grade'
    auc.t.combined = rbind(auc.t.combined,auc.t)
    
    
    auc.res1 = rbind(
      auc.plot.all.g1[[2]][1,],
      auc.plot.all.g2[[2]][1,],
      auc.plot.all.g3[[2]][1,],
      auc.plot.all.g9[[2]][1,])
    
    auc.res1$title = 'morph'
    combined.auroc = combined.auroc[,colnames(combined.auroc) %in% colnames(auc.res1)]
    
    auc.res1 = auc.res1[,colnames(auc.res1) %in% colnames(combined.auroc)]
    combined.auroc = rbind(combined.auroc, auc.res1)
    grade_colors = c('Control'= 'grey','1' = '#82A6B1','2' = '#35605A', '3' = '#2A324B','9'='#66462C')
    grade_colors.sample = grade_colors
    plot1 = ggplot(auc.t[auc.t$time >= 365*1,],aes(x = time/365, y = AUC, col = var)) + 
      geom_line(linewidth=1) +
      scale_color_manual(values = c(grade_colors))+ #+ ggtitle(title) +
      scale_y_continuous(limits=c(0,1))+
      theme_bw()+ 
      #facet_grid2(.~title)+
      theme(text = element_text(size=8),
            axis.text=element_text(size=8, face = "bold"),
            axis.title=element_text(size=8,face="bold"),
            legend.position = "none")+ guides(colour = guide_legend(override.aes = list(size=4),nrow=4)) +
      xlab('Time (Years)') + ylab('Time-Dependent CV AUROC')
    
    png(paste0(name,'.filler.',c,'.timeauroc.grade.png'),height = 1100,width=1100,res = 300)
    print(plot1)
    dev.off()
    
    
    combined.auc.plot.all = rbind(auc.plot.all.g1[[2]],auc.plot.all.g2[[2]],auc.plot.all.g3[[2]],auc.plot.all.g9[[2]])
    plot1 = ggplot(combined.auc.plot.all,aes(x = fpr, y = tpr, col = grade)) + geom_line(linewidth=1) +
      scale_color_manual(values = c(grade_colors))+ #+ ggtitle(title) +
      theme_bw()+ 
      theme(text = element_text(size=8),
            axis.text=element_text(size=8, face = "bold"),
            axis.title=element_text(size=8,face="bold"),
            legend.position = "none")+ guides(colour = guide_legend(override.aes = list(size=4),nrow=4)) + xlab('False Positive Rate') + ylab('Sensitivity')
    
    png(paste0(name,'.filler.',c,'.auroc.grade.png'),height = 1000,width=1000,res = 300)
    print(plot1)
    dev.off()
    
    plot1 = ggplot(combined.auc.plot.all,aes(x = fpr, y = tpr, col = grade)) + geom_line(linewidth=1) +
      scale_color_manual(values = c(grade_colors))+ #+ ggtitle(title) +
      theme_bw()+ 
      theme(text = element_text(size=8),
            axis.text=element_text(size=8, face = "bold"),
            axis.title=element_text(size=8,face="bold"),
            legend.position = "right")+ guides(colour = guide_legend(override.aes = list(size=4),nrow=4)) + xlab('False Positive Rate') + ylab('Sensitivity')
    
    png(paste0(name,'.filler.',c,'.auroc.grade.labs.png'),height = 1000,width=1000,res = 300)
    print(plot1)
    dev.off()
    merged.df.all.tmp = data.frame(merged.df.all.tmp,check.names=F)
    pred.df.targ.collapse.annotated= merge(pred.df.targ.collapse,merged.df.all.tmp[,c('GRP_Id','GRADE_CD')],by='GRP_Id')
    pred.df.targ.collapse.annotated$GRADE = ifelse(pred.df.targ.collapse.annotated$reported == 'Control','Control',pred.df.targ.collapse.annotated$GRADE_CD)
    pred.df.targ.collapse.annotated$GRADE = ifelse(pred.df.targ.collapse.annotated$GRADE == '9','NR',pred.df.targ.collapse.annotated$GRADE)
    pred.df.targ.collapse.annotated[is.na(pred.df.targ.collapse.annotated$GRADE),'GRADE'] ='NR'
    my_comparisons= list()
    if (wilcox.test(methylation_score ~ reported, pred.df.targ.collapse.annotated[pred.df.targ.collapse.annotated$GRADE %in% c('Control','1'),])$p.value < 0.05 ) {
      my_comparisons[[length(my_comparisons)+1]] = c('Control','1')
    }
    
    if (wilcox.test(methylation_score ~ reported, pred.df.targ.collapse.annotated[pred.df.targ.collapse.annotated$GRADE %in% c('Control','2'),])$p.value < 0.05 ) {
      my_comparisons[[length(my_comparisons)+1]] = c('Control','2')
    }
    
    if (wilcox.test(methylation_score ~ reported, pred.df.targ.collapse.annotated[pred.df.targ.collapse.annotated$GRADE %in% c('Control','3'),])$p.value < 0.05 ) {
      my_comparisons[[length(my_comparisons)+1]] = c('Control','3')
    }
    
    if (wilcox.test(methylation_score ~ reported, pred.df.targ.collapse.annotated[pred.df.targ.collapse.annotated$GRADE %in% c('Control','NR'),])$p.value < 0.05 ) {
      my_comparisons[[length(my_comparisons)+1]] = c('Control','NR')
    }
    
    
    
    print(kruskal.test(methylation_score ~ GRADE,  data = pred.df.targ.collapse.annotated[as.character(pred.df.targ.collapse.annotated$GRADE) %in%  c('Control','1','2','3'),])  )
    grade_colors.sample = c('Control'= 'grey','1' = '#82A6B1','2' = '#35605A', '3' = '#2A324B','NR'='#66462C')
    
    
    plot1 = ggplot(pred.df.targ.collapse.annotated,aes(x = GRADE, y = methylation_score, col = GRADE)) + geom_boxplot(outlier.shape=NA)+geom_jitter(width=0.1)+
      scale_color_manual(values = c(grade_colors.sample))+ #+ ggtitle(title) +
      theme_bw()+ 
      scale_y_continuous(limits=c(0,1+0.1*length(my_comparisons)))+
      
      theme(text = element_text(size=8),
            axis.text=element_text(size=8, face = "bold"),
            axis.title=element_text(size=8,face="bold"),
            legend.position = "none")+ guides(colour = guide_legend(override.aes = list(size=4),nrow=4)) + 
      xlab('Grade at Diagnosis') +
      ylab('Methylation Score')+
      stat_compare_means(comparisons = my_comparisons,label.y=seq(1,1+length(my_comparisons)*0.1,0.1),size=3)
    
    png(paste0(name,'.grade.methscore.',c,'.png'),height = 1100,width=1000,res = 300)
    print(plot1)
    dev.off()
    
    
    
    #plotting age##
    
    library(ggh4x)
    age_grouping2 = function(age) {
      tmp = ifelse(age >= 35 & age < 45, '35-45', age)
      tmp = ifelse(age >= 45 & age < 55, '45-55', tmp)
      tmp = ifelse(age >= 55 & age < 65, '55-65', tmp)
      tmp = ifelse(age >= 65 & age < 75, '65-75', tmp)
      tmp = ifelse(age >= 75, '75+', tmp)
      
      age_groups = c('35-45','45-55','55-65','65-75','75+')
      tmp = factor(tmp, levels = age_groups) 
    }
    
    merged.df.all.tmp$age_group2 = age_grouping2(merged.df.all.tmp$SDC_AGE_CALC)
    
    age1 = merged.df.all.tmp[merged.df.all.tmp$age_group2 %in% ('35-45') | merged.df.all.tmp$Cancer =='Control',]
    age2 = merged.df.all.tmp[merged.df.all.tmp$age_group2 %in% c('45-55') | merged.df.all.tmp$Cancer =='Control',]
    age3 = merged.df.all.tmp[merged.df.all.tmp$age_group2 %in% c('55-65') | merged.df.all.tmp$Cancer =='Control',]
    age4 = merged.df.all.tmp[merged.df.all.tmp$age_group2 %in% c('65-75') | merged.df.all.tmp$Cancer =='Control',]
    
    auc.plot.all = summary.calc(pred.df.targ.collapse,merged.df.all.tmp)
    auc.plot.all[[2]]$var = 'All'
    auc.plot.all[[2]]$var.group = 'Age'
    
    auc.plot.all.s1 = summary.calc(pred.df.targ.collapse[pred.df.targ.collapse$GRP_Id %in% age1$GRP_Id,],merged.df.all.tmp)
    auc.plot.all.s1[[2]]$var = '35-45'
    auc.plot.all.s1[[2]]$var.group = 'Age'
    
    auc.plot.all.s2 = summary.calc(pred.df.targ.collapse[pred.df.targ.collapse$GRP_Id %in% age2$GRP_Id,],merged.df.all.tmp)
    auc.plot.all.s2[[2]]$var = '45-55'
    auc.plot.all.s2[[2]]$var.group = 'Age'
    
    
    auc.plot.all.s3 = summary.calc(pred.df.targ.collapse[pred.df.targ.collapse$GRP_Id %in% age3$GRP_Id,],merged.df.all.tmp)
    auc.plot.all.s3[[2]]$var = '55-65'
    auc.plot.all.s3[[2]]$var.group = 'Age'
    
    
    auc.plot.all.s4 = summary.calc(pred.df.targ.collapse[pred.df.targ.collapse$GRP_Id %in% age4$GRP_Id,],merged.df.all.tmp)
    auc.plot.all.s4[[2]]$var = '65-75'
    auc.plot.all.s4[[2]]$var.group = 'Age'
    
    
    auc.res1 = rbind(
      auc.plot.all.s1[[2]][1,],
      auc.plot.all.s2[[2]][1,],
      auc.plot.all.s3[[2]][1,],
      auc.plot.all.s4[[2]][1,])
    auc.res1$title = 'Age'
    combined.auroc = rbind(combined.auroc,auc.res1[,colnames(auc.res1) %in% colnames(combined.auroc)])
    
    age_colors = c('All' = "#000004FF",'35-45' = '#C5C392','45-55'= '#FFB20F','55-65'='#FF4B3E', '65-75'='#972D07')
    
    
    #annotating auc t
    auc.plot.all[[3]]$var = 'All'
    auc.plot.all.s1[[3]]$var = '35-45'
    auc.plot.all.s2[[3]]$var = '45-55'
    auc.plot.all.s3[[3]]$var = '55-65'
    auc.plot.all.s4[[3]]$var = '65-75'
    
    
    auc.t =rbind( auc.plot.all[[3]],
                  auc.plot.all.s1[[3]],
                  auc.plot.all.s2[[3]],
                  auc.plot.all.s3[[3]],
                  auc.plot.all.s4[[3]]
                  
    )
    auc.t$title = 'Age'
    auc.t.combined = rbind(auc.t.combined,auc.t)
    
    #auc t
    
    plot1 = ggplot(auc.t[auc.t$time >= 365*1,],aes(x = time/365, y = AUC, col = var)) + geom_line(linewidth=1) +
      scale_color_manual(values = c(age_colors))+ #+ ggtitle(title) +
      scale_y_continuous(limits=c(0,1))+
      theme_bw()+ 
      facet_grid2(.~title)+
      theme(text = element_text(size=8),
            axis.text=element_text(size=8, face = "bold"),
            axis.title=element_text(size=8,face="bold"),
            legend.position = "none")+ guides(colour = guide_legend(override.aes = list(size=4),nrow=4)) +
      xlab('Time (Years)') + ylab('Time-Dependent CV AUROC')
    
    png(paste0(name,'.filler.',c,'.timeauroc.age.png'),height = 1100,width=1100,res = 300)
    print(plot1)
    dev.off()
    
    
    combined.auc.plot.all = rbind(
      auc.plot.all.s1[[2]],
      auc.plot.all.s2[[2]],
      auc.plot.all.s3[[2]],
      auc.plot.all.s4[[2]])
    plot1 = ggplot(combined.auc.plot.all,aes(x = fpr, y = tpr, col = var)) + geom_line(linewidth=1) +
      scale_color_manual(values = c(age_colors))+ #+ ggtitle(title) +
      theme_bw()+ 
      theme(text = element_text(size=8),
            axis.text=element_text(size=8, face = "bold"),
            axis.title=element_text(size=8,face="bold"),
            legend.position = "none")+ guides(colour = guide_legend(override.aes = list(size=4),nrow=4)) + xlab('False Positive Rate') + ylab('Sensitivity')
    
    png(paste0(name,'.filler.',c,'.auroc.age.png'),height = 1000,width=1000,res = 300)
    print(plot1)
    dev.off()
    
    plot1 = ggplot(combined.auc.plot.all,aes(x = fpr, y = tpr, col = var)) + geom_line(linewidth=1) +
      scale_color_manual(values = c(age_colors))+ #+ ggtitle(title) +
      theme_bw()+ 
      theme(text = element_text(size=8),
            axis.text=element_text(size=8, face = "bold"),
            axis.title=element_text(size=8,face="bold"),
            legend.position = "right")+ guides(colour = guide_legend(override.aes = list(size=4),nrow=4)) + xlab('False Positive Rate') + ylab('Sensitivity')
    
    png(paste0(name,'.filler.',c,'.auroc.age.labs.png'),height = 1000,width=1000,res = 300)
    print(plot1)
    dev.off()
    
    
    
    
    
    plot1 = ggplot(auc.res1,aes(x = var, y = auc, col = var)) +
      geom_errorbar(aes(ymin=auc.lower, ymax=auc.upper), width=.2,
                    position=position_dodge(0.05))+
      geom_point()+
      scale_y_continuous(limits=c(0,1))+
      scale_color_manual(values = c(age_colors))+ #+ ggtitle(title) +
      theme_bw()+ 
      #facet_grid(.~title)+
      
      theme(text = element_text(size=8),
            axis.text=element_text(size=8, face = "bold"),
            axis.title=element_text(size=8,face="bold"),
            legend.position = "none")+ 
      guides(colour = guide_legend(override.aes = list(size=4),nrow=4)) +
      xlab('Age') + ylab('CV AUROC (95% CI)')
    
    png(paste0(name,'.filler.',c,'.auroc.ci.age.png'),height = 1100,width=1100,res = 300)
    print(plot1)
    dev.off()
    plot1 = ggplot(auc.res1,aes(x = var, y = ci, col = var)) +
      geom_errorbar(aes(ymin=ci.lower, ymax=ci.upper), width=.2,
                    position=position_dodge(0.05))+
      geom_point()+
      scale_y_continuous(limits=c(0,1))+
      scale_color_manual(values = c(age_colors))+ #+ ggtitle(title) +
      theme_bw()+ 
      #facet_grid(.~title)+
      theme(text = element_text(size=8),
            axis.text=element_text(size=8, face = "bold"),
            axis.title=element_text(size=8,face="bold"),
            legend.position = "none")+ 
      guides(colour = guide_legend(override.aes = list(size=4),nrow=4)) +
      xlab('age at Diagnosis') + ylab('CV Concordance Index (95% CI)')
    
    png(paste0(name,'.filler.',c,'.ci.ci.age.png'),height = 1100,width=1100,res = 300)
    print(plot1)
    dev.off()
    plot1 = ggplot(auc.res1,aes(x = var, y =se.spec.95 , col = var)) +
      geom_errorbar(aes(ymin=se.spec.95.lower, ymax=se.spec.95.upper), width=.2,
                    position=position_dodge(0.05))+
      geom_point()+
      scale_y_continuous(limits=c(0,1))+
      scale_color_manual(values = c(age_colors))+ #+ ggtitle(title) +
      theme_bw()+ 
      #facet_grid(.~title)+
      theme(text = element_text(size=8),
            axis.text=element_text(size=8, face = "bold"),
            axis.title=element_text(size=8,face="bold"),
            legend.position = "none")+ 
      guides(colour = guide_legend(override.aes = list(size=4),nrow=4)) +
      xlab('Age') + ylab('Sensitivity at 95% Specificity (95% CI)')
    
    png(paste0(name,'.filler.',c,'.sens.spec95.ci.age.png'),height = 1100,width=1100,res = 300)
    print(plot1)
    dev.off()
    
    plot1 = ggplot(auc.res1,aes(x = var, y = se.spec.90, col = var)) +
      geom_errorbar(aes(ymin=se.spec.90.lower, ymax=se.spec.90.upper), width=.2,
                    position=position_dodge(0.05))+
      geom_point()+
      scale_y_continuous(limits=c(0,1))+
      scale_color_manual(values = c(age_colors))+ #+ ggtitle(title) +
      theme_bw()+ 
      #facet_grid(.~title)+
      theme(text = element_text(size=8),
            axis.text=element_text(size=8, face = "bold"),
            axis.title=element_text(size=8,face="bold"),
            legend.position = "none")+ 
      guides(colour = guide_legend(override.aes = list(size=4),nrow=4)) +
      xlab('Age')+ ylab('Sensitivity at 90% Specificity (95% CI)')
    
    png(paste0(name,'.filler.',c,'.sens.spec90.ci.age.png'),height = 1100,width=1100,res = 300)
    print(plot1)
    dev.off()
    
    plot1 = ggplot(combined.auc.plot.all,aes(x = fpr, y = tpr, col = var)) + geom_line(linewidth=1) +
      scale_color_manual(values = c(age_colors))+ #+ ggtitle(title) +
      theme_bw()+ 
      theme(text = element_text(size=8),
            axis.text=element_text(size=8, face = "bold"),
            axis.title=element_text(size=8,face="bold"),
            legend.position = "right")+ guides(colour = guide_legend(override.aes = list(size=4),ncol=1)) + 
      xlab('False Positive Rate') + ylab('Sensitivity')
    
    png(paste0(name,'.filler.',c,'.auroc.age.labs.png'),height = 1100,width=1100,res = 300)
    print(plot1)
    dev.off()
    
    pred.df.targ.collapse.annotated = merge(pred.df.targ.collapse, merged.df.all.tmp[,c('GRP_Id','Cancer','Diagnosis_Time','filler','age_group2','GRADE_CD')],by='GRP_Id')
    
    #boxplot of score
    plot1 = ggplot(pred.df.targ.collapse.annotated[pred.df.targ.collapse.annotated$Cancer %in% c('Control','Breast'),],aes(x = age_group2, y = methylation_score, col = age_group2)) +
      geom_boxplot(outlier.shape=NA)+geom_jitter(width=0.1)+
      scale_color_manual(values = c(age_colors))+ #+ ggtitle(title) +
      theme_bw()+ 
      scale_y_continuous(limits=c(0,1))+
      
      facet_grid2(. ~ Cancer)+
      
      theme(text = element_text(size=8),
            axis.text=element_text(size=8, face = "bold"),
            axis.title=element_text(size=8,face="bold"),
            legend.position = "none")+ 
      guides(colour = guide_legend(override.aes = list(size=4),nrow=4)) +
      xlab('Age') + ylab('Methylation Score')
    
    png(paste0(name,'.age.methscore.',c,'.png'),height = 900,width=1100,res = 300)
    print(plot1)
    dev.off()
    
    
    
    
    
    #plotting dxage##
    
    library(ggh4x)
    age_grouping3 = function(age) {
      tmp = ifelse(age >= 30 & age < 50, '30-50', age)
      # tmp = ifelse(age >= 40 & age < 50, '40-50', tmp)
      tmp = ifelse(age >= 50 & age < 60, '50-60', tmp)
      tmp = ifelse(age >= 60 & age < 70, '60-70', tmp)
      tmp = ifelse(age >= 70, '70-80', tmp)
      
      age_groups = c('30-50','50-60','60-70','70-80')
      tmp = factor(tmp, levels = age_groups) 
    }
    #targline 
    merged.df.all.tmp$dx.age = ifelse(merged.df.all.tmp$Cancer != 'Control',merged.df.all.tmp$SDC_AGE_CALC + merged.df.all.tmp$diff_in_days/365,merged.df.all.tmp$SDC_AGE_CALC)
    merged.df.all.tmp$age_group3 = age_grouping3(merged.df.all.tmp$dx.age)
    
    age1 = merged.df.all.tmp[merged.df.all.tmp$age_group3 %in% ('30-50') | merged.df.all.tmp$Cancer =='Control',]
    age2 = merged.df.all.tmp[merged.df.all.tmp$age_group3 %in% c('50-60') | merged.df.all.tmp$Cancer =='Control',]
    age3 = merged.df.all.tmp[merged.df.all.tmp$age_group3 %in% c('60-70') | merged.df.all.tmp$Cancer =='Control',]
    age4 = merged.df.all.tmp[merged.df.all.tmp$age_group3 %in% c('70-80') | merged.df.all.tmp$Cancer =='Control',]
    
    
    
    auc.plot.all.s1 = summary.calc(pred.df.targ.collapse[pred.df.targ.collapse$GRP_Id %in% age1$GRP_Id,],merged.df.all.tmp)
    auc.plot.all.s1[[2]]$var = '30-50'
    auc.plot.all.s1[[2]]$var.group = 'DxAge'
    
    auc.plot.all.s2 = summary.calc(pred.df.targ.collapse[pred.df.targ.collapse$GRP_Id %in% age2$GRP_Id,],merged.df.all.tmp)
    auc.plot.all.s2[[2]]$var = '50-60'
    auc.plot.all.s2[[2]]$var.group = 'DxAge'
    
    
    auc.plot.all.s3 = summary.calc(pred.df.targ.collapse[pred.df.targ.collapse$GRP_Id %in% age3$GRP_Id,],merged.df.all.tmp)
    auc.plot.all.s3[[2]]$var = '60-70'
    auc.plot.all.s3[[2]]$var.group = 'DxAge'
    
    
    auc.plot.all.s4 = summary.calc(pred.df.targ.collapse[pred.df.targ.collapse$GRP_Id %in% age4$GRP_Id,],merged.df.all.tmp)
    auc.plot.all.s4[[2]]$var = '70-80'
    auc.plot.all.s4[[2]]$var.group = 'DxAge'
    
    
    auc.res1 = rbind(
      auc.plot.all.s1[[2]][1,],
      auc.plot.all.s2[[2]][1,],
      auc.plot.all.s3[[2]][1,],
      auc.plot.all.s4[[2]][1,])
    auc.res1$title = 'DxAge'
    combined.auroc = rbind(combined.auroc,auc.res1[,colnames(auc.res1) %in% colnames(combined.auroc)])
    
    dxage_colors = c('All' = "#000004FF",'30-50' = '#C5C392','50-60'= '#FFB20F','60-70'='#FF4B3E', '70-80'='#972D07')
    
    
    #annotating auc t
    auc.plot.all.s1[[3]]$var = '30-50'
    auc.plot.all.s2[[3]]$var = '50-60'
    auc.plot.all.s3[[3]]$var = '60-70'
    auc.plot.all.s4[[3]]$var = '70-80'
    
    
    auc.t =rbind( 
      auc.plot.all.s1[[3]],
      auc.plot.all.s2[[3]],
      auc.plot.all.s3[[3]],
      auc.plot.all.s4[[3]]
      
    )
    auc.t$title = 'DxAge'
    auc.t.combined = rbind(auc.t.combined,auc.t)
    
    #auc t
    
    plot1 = ggplot(auc.t[auc.t$time >= 365*1,],aes(x = time/365, y = AUC, col = var)) + geom_line(linewidth=1) +
      scale_color_manual(values = c(dxage_colors))+ #+ ggtitle(title) +
      scale_y_continuous(limits=c(0,1))+
      theme_bw()+ 
      #facet_grid2(.~title)+
      theme(text = element_text(size=8),
            axis.text=element_text(size=8, face = "bold"),
            axis.title=element_text(size=8,face="bold"),
            legend.position = "none")+ guides(colour = guide_legend(override.aes = list(size=4),nrow=4)) +
      xlab('Time (Years)') + ylab('Time-Dependent CV AUROC')
    
    png(paste0(name,'.filler.',c,'.timeauroc.dxage.png'),height = 1100,width=1100,res = 300)
    print(plot1)
    dev.off()
    
    
    combined.auc.plot.all = rbind(
      auc.plot.all.s1[[2]],
      auc.plot.all.s2[[2]],
      auc.plot.all.s3[[2]],
      auc.plot.all.s4[[2]])
    plot1 = ggplot(combined.auc.plot.all,aes(x = fpr, y = tpr, col = var)) + geom_line(linewidth=1) +
      scale_color_manual(values = c(dxage_colors))+ #+ ggtitle(title) +
      theme_bw()+ 
      theme(text = element_text(size=8),
            axis.text=element_text(size=8, face = "bold"),
            axis.title=element_text(size=8,face="bold"),
            legend.position = "none")+ guides(colour = guide_legend(override.aes = list(size=4),nrow=4)) + xlab('False Positive Rate') + ylab('Sensitivity')
    
    png(paste0(name,'.filler.',c,'.auroc.dxage.png'),height = 1100,width=1100,res = 300)
    print(plot1)
    dev.off()
    
    plot1 = ggplot(combined.auc.plot.all,aes(x = fpr, y = tpr, col = var)) + geom_line(linewidth=1) +
      scale_color_manual(values = c(dxage_colors))+ #+ ggtitle(title) +
      theme_bw()+ 
      theme(text = element_text(size=8),
            axis.text=element_text(size=8, face = "bold"),
            axis.title=element_text(size=8,face="bold"),
            legend.position = "right")+ guides(colour = guide_legend(override.aes = list(size=4),nrow=4)) + xlab('False Positive Rate') + ylab('Sensitivity')
    
    png(paste0(name,'.filler.',c,'.auroc.dxage.labs.png'),height = 1100,width=1100,res = 300)
    print(plot1)
    dev.off()
    
    
    
    
    
    plot1 = ggplot(auc.res1,aes(x = var, y = auc, col = var)) +
      geom_errorbar(aes(ymin=auc.lower, ymax=auc.upper), width=.2,
                    position=position_dodge(0.05))+
      geom_point()+
      scale_y_continuous(limits=c(0,1))+
      scale_color_manual(values = c(dxage_colors))+ #+ ggtitle(title) +
      theme_bw()+ 
      #facet_grid(.~title)+
      
      theme(text = element_text(size=8),
            axis.text=element_text(size=8, face = "bold"),
            axis.title=element_text(size=8,face="bold"),
            legend.position = "none")+ 
      guides(colour = guide_legend(override.aes = list(size=4),nrow=4)) +
      xlab('Age') + ylab('CV AUROC (95% CI)')
    
    png(paste0(name,'.filler.',c,'.auroc.ci.dxage.png'),height = 1100,width=1100,res = 300)
    print(plot1)
    dev.off()
    plot1 = ggplot(auc.res1,aes(x = var, y = ci, col = var)) +
      geom_errorbar(aes(ymin=ci.lower, ymax=ci.upper), width=.2,
                    position=position_dodge(0.05))+
      geom_point()+
      scale_y_continuous(limits=c(0,1))+
      scale_color_manual(values = c(dxage_colors))+ #+ ggtitle(title) +
      theme_bw()+ 
      #facet_grid(.~title)+
      theme(text = element_text(size=8),
            axis.text=element_text(size=8, face = "bold"),
            axis.title=element_text(size=8,face="bold"),
            legend.position = "none")+ 
      guides(colour = guide_legend(override.aes = list(size=4),nrow=4)) +
      xlab('Diagnosis Age') + ylab('CV Concordance Index (95% CI)')
    
    png(paste0(name,'.filler.',c,'.ci.ci.dxage.png'),height = 1100,width=1100,res = 300)
    print(plot1)
    dev.off()
    plot1 = ggplot(auc.res1,aes(x = var, y =se.spec.95 , col = var)) +
      geom_errorbar(aes(ymin=se.spec.95.lower, ymax=se.spec.95.upper), width=.2,
                    position=position_dodge(0.05))+
      geom_point()+
      scale_y_continuous(limits=c(0,1))+
      scale_color_manual(values = c(dxage_colors))+ #+ ggtitle(title) +
      theme_bw()+ 
      #facet_grid(.~title)+
      theme(text = element_text(size=8),
            axis.text=element_text(size=8, face = "bold"),
            axis.title=element_text(size=8,face="bold"),
            legend.position = "none")+ 
      guides(colour = guide_legend(override.aes = list(size=4),nrow=4)) +
      xlab('Diagnosis Age') + ylab('Sensitivity at 95% Specificity (95% CI)')
    
    png(paste0(name,'.filler.',c,'.sens.spec95.ci.dxage.png'),height = 1100,width=1100,res = 300)
    print(plot1)
    dev.off()
    
    plot1 = ggplot(auc.res1,aes(x = var, y = se.spec.90, col = var)) +
      geom_errorbar(aes(ymin=se.spec.90.lower, ymax=se.spec.90.upper), width=.2,
                    position=position_dodge(0.05))+
      geom_point()+
      scale_y_continuous(limits=c(0,1))+
      scale_color_manual(values = c(dxage_colors))+ #+ ggtitle(title) +
      theme_bw()+ 
      #facet_grid(.~title)+
      theme(text = element_text(size=8),
            axis.text=element_text(size=8, face = "bold"),
            axis.title=element_text(size=8,face="bold"),
            legend.position = "none")+ 
      guides(colour = guide_legend(override.aes = list(size=4),nrow=4)) +
      xlab('Diagnosis Age')+ ylab('Sensitivity at 90% Specificity (95% CI)')
    
    png(paste0(name,'.filler.',c,'.sens.spec90.ci.dxage.png'),height = 1100,width=1100,res = 300)
    print(plot1)
    dev.off()
    
    plot1 = ggplot(combined.auc.plot.all,aes(x = fpr, y = tpr, col = var)) + geom_line(linewidth=1) +
      scale_color_manual(values = c(dxage_colors))+ #+ ggtitle(title) +
      theme_bw()+ 
      theme(text = element_text(size=8),
            axis.text=element_text(size=8, face = "bold"),
            axis.title=element_text(size=8,face="bold"),
            legend.position = "right")+ guides(colour = guide_legend(override.aes = list(size=4),ncol=1)) + 
      xlab('False Positive Rate') + ylab('Sensitivity')
    
    png(paste0(name,'.filler.',c,'.auroc.dxage.labs.png'),height = 1100,width=1100,res = 300)
    print(plot1)
    dev.off()
    
    pred.df.targ.collapse.annotated = merge(pred.df.targ.collapse, merged.df.all.tmp[,c('GRP_Id','Cancer','Diagnosis_Time','filler','age_group3','GRADE_CD')],by='GRP_Id')
    
    #boxplot of score
    plot1 = ggplot(pred.df.targ.collapse.annotated[pred.df.targ.collapse.annotated$Cancer %in% c('Control','Breast'),],aes(x = age_group3, y = methylation_score, col = age_group3)) +
      geom_boxplot(outlier.shape=NA)+geom_jitter(width=0.1)+
      scale_color_manual(values = c(dxage_colors))+ #+ ggtitle(title) +
      theme_bw()+ 
      scale_y_continuous(limits=c(0,1))+
      
      facet_grid2(. ~ Cancer)+
      
      theme(text = element_text(size=8),
            axis.text=element_text(size=8, face = "bold"),
            axis.title=element_text(size=8,face="bold"),
            legend.position = "none")+ 
      guides(colour = guide_legend(override.aes = list(size=4),nrow=4)) +
      xlab('Diagnosis Age') + ylab('Methylation Score')
    
    png(paste0(name,'.dxage.methscore.',c,'.png'),height = 900,width=1100,res = 300)
    print(plot1)
    dev.off()
    
    
    #last mmg
    print('mmg ever')
    
    last.mmg = read.csv('/.mounts/labs/awadallalab/private/ncheng/cfmedip_data/cptp_samples/participant_data/ohs.sample.qx/combined.baseline.additional.var.csv',header=T,stringsAsFactors = F)
    last.mmg$mammogram = ifelse(last.mmg$HS_MMG_LAST == '-7','Never Had','Unknown')
    last.mmg$mammogram = ifelse(last.mmg$HS_MMG_LAST == '1','< 6 months',last.mmg$mammogram)
    last.mmg$mammogram = ifelse(last.mmg$HS_MMG_LAST == '2','0.5-1 year',last.mmg$mammogram)
    last.mmg$mammogram = ifelse(last.mmg$HS_MMG_LAST == '3','1-2 years',last.mmg$mammogram)
    last.mmg$mammogram = ifelse(last.mmg$HS_MMG_LAST == '4','2-3 years',last.mmg$mammogram)
    last.mmg$mammogram = ifelse(last.mmg$HS_MMG_LAST == '5','3+ years',last.mmg$mammogram)
    
    
    library(ggh4x)
    
    
    merged.df.all.tmp = merge(merged.df.all.tmp,last.mmg[,c('ResearchId','mammogram')],all.x=T)
    age0 = merged.df.all.tmp[merged.df.all.tmp$mammogram %in% ('Never Had') | merged.df.all.tmp$group =='Control',]
    age1 = merged.df.all.tmp[merged.df.all.tmp$mammogram %in% ('< 6 months') | merged.df.all.tmp$group =='Control',]
    age2 = merged.df.all.tmp[merged.df.all.tmp$mammogram %in% c('0.5-1 year')| merged.df.all.tmp$group =='Control',]
    age3 = merged.df.all.tmp[merged.df.all.tmp$mammogram %in% c('1-2 years')| merged.df.all.tmp$group =='Control',]
    age4 = merged.df.all.tmp[merged.df.all.tmp$mammogram %in% c('2-3 years','3+ years')| merged.df.all.tmp$group =='Control',]
    
    
    auc.plot.all = summary.calc(pred.df.targ.collapse,merged.df.all.tmp)
    auc.plot.all[[2]]$var = 'All'
    auc.plot.all[[2]]$var.group = 'Last Mammogram'
    
    auc.plot.all.s0 = summary.calc(pred.df.targ.collapse[pred.df.targ.collapse$GRP_Id %in% age0$GRP_Id,],merged.df.all.tmp)
    auc.plot.all.s0[[2]]$var = 'Never'
    auc.plot.all.s0[[2]]$var.group =  'Last Mammogram'
    
    auc.plot.all.s1 = summary.calc(pred.df.targ.collapse[pred.df.targ.collapse$GRP_Id %in% age1$GRP_Id,],merged.df.all.tmp)
    auc.plot.all.s1[[2]]$var = '< 0.5'
    auc.plot.all.s1[[2]]$var.group =  'Last Mammogram'
    
    auc.plot.all.s2 = summary.calc(pred.df.targ.collapse[pred.df.targ.collapse$GRP_Id %in% age2$GRP_Id,],merged.df.all.tmp)
    auc.plot.all.s2[[2]]$var = '0.5-1'
    auc.plot.all.s2[[2]]$var.group =  'Last Mammogram'
    
    
    auc.plot.all.s3 = summary.calc(pred.df.targ.collapse[pred.df.targ.collapse$GRP_Id %in% age3$GRP_Id,],merged.df.all.tmp)
    auc.plot.all.s3[[2]]$var = '1-2'
    auc.plot.all.s3[[2]]$var.group =  'Last Mammogram'
    
    
    auc.plot.all.s4 = summary.calc(pred.df.targ.collapse[pred.df.targ.collapse$GRP_Id %in% age4$GRP_Id,],merged.df.all.tmp)
    auc.plot.all.s4[[2]]$var = '2+'
    auc.plot.all.s4[[2]]$var.group =  'Last Mammogram'
    
    
    
    auc.res1 = rbind(
      auc.plot.all.s0[[2]][1,],
      auc.plot.all.s1[[2]][1,],
      auc.plot.all.s2[[2]][1,],
      auc.plot.all.s3[[2]][1,],
      auc.plot.all.s4[[2]][1,])
    auc.res1$title = 'Last Mammogram'
    #combined.auroc = combined.auroc[,colnames(combined.auroc) %in% colnames(auc.res1)]
    combined.auroc = rbind(combined.auroc,auc.res1[,colnames(auc.res1) %in% colnames(combined.auroc)])
    
    mmg_colors = c('All' = "#000004FF",
                   'Never' = '#93A3B1',
                   '< 0.5'='#78BC61',
                   '0.5-1'='#C0C781',
                   '1-2'='#C59B76',
                   '2+'='#AD343E')
    
    #annotating auc t
    auc.plot.all[[3]]$var = 'All'
    auc.plot.all.s0[[3]]$var = 'Never'
    auc.plot.all.s1[[3]]$var = '< 0.5'
    auc.plot.all.s2[[3]]$var = '0.5-1'
    auc.plot.all.s3[[3]]$var = '1-2'
    auc.plot.all.s4[[3]]$var = '2+'
    
    
    auc.t =rbind(auc.plot.all.s0[[3]],
                 auc.plot.all.s1[[3]],
                 auc.plot.all.s2[[3]],
                 auc.plot.all.s3[[3]],
                 auc.plot.all.s4[[3]])
    
    
    auc.t$title = 'Last Mammogram'
    auc.t.combined = rbind(auc.t.combined,auc.t)
    
    #auc t
    
    plot1 = ggplot(auc.t[auc.t$time >= 365*1,],aes(x = time/365, y = AUC, col = var)) + geom_line(linewidth=1) +
      scale_color_manual(values = c(mmg_colors))+ #+ ggtitle(title) +
      scale_y_continuous(limits=c(0,1))+
      theme_bw()+ 
      #facet_grid2(.~title)+
      theme(text = element_text(size=8),
            axis.text=element_text(size=8, face = "bold"),
            axis.title=element_text(size=8,face="bold"),
            legend.position = "none")+ guides(colour = guide_legend(override.aes = list(size=4),nrow=4)) +
      xlab('Time (Years)') + ylab('Time-Dependent CV AUROC')
    
    png(paste0(name,'.filler.',c,'.timeauroc.mmg.png'),height = 1100,width=1100,res = 300)
    print(plot1)
    dev.off()
    
    
    combined.auc.plot.all = rbind(
      auc.plot.all.s0[[2]],
      auc.plot.all.s1[[2]],
      auc.plot.all.s2[[2]],
      auc.plot.all.s3[[2]],
      auc.plot.all.s4[[2]])
    plot1 = ggplot(combined.auc.plot.all,aes(x = fpr, y = tpr, col = var)) + geom_line(linewidth=1) +
      scale_color_manual(values = c(mmg_colors))+ #+ ggtitle(title) +
      theme_bw()+ 
      theme(text = element_text(size=8),
            axis.text=element_text(size=8, face = "bold"),
            axis.title=element_text(size=8,face="bold"),
            legend.position = "none")+ guides(colour = guide_legend(override.aes = list(size=4),nrow=4)) + xlab('False Positive Rate') + ylab('Sensitivity')
    
    png(paste0(name,'.filler.',c,'.auroc.mmg.png'),height = 1000,width=1000,res = 300)
    print(plot1)
    dev.off()
    
    plot1 = ggplot(combined.auc.plot.all,aes(x = fpr, y = tpr, col = var)) + geom_line(linewidth=1) +
      scale_color_manual(values = c(mmg_colors))+ #+ ggtitle(title) +
      theme_bw()+ 
      theme(text = element_text(size=8),
            axis.text=element_text(size=8, face = "bold"),
            axis.title=element_text(size=8,face="bold"),
            legend.position = "right")+ guides(colour = guide_legend(override.aes = list(size=4),nrow=4)) + xlab('False Positive Rate') + ylab('Sensitivity')
    
    png(paste0(name,'.filler.',c,'.auroc.mmg.labs.png'),height = 1000,width=1000,res = 300)
    print(plot1)
    dev.off()
    pred.df.targ.collapse.annotated = merge(pred.df.targ.collapse, last.mmg[,c('ResearchId','mammogram')],all.x=T)
    plot1 = ggplot(pred.df.targ.collapse.annotated[!is.na(pred.df.targ.collapse.annotated$mammogram),],aes(x = mammogram, y = methylation_score, col = mammogram)) +
      geom_boxplot(outlier.shape=NA)+geom_jitter(width=0.1)+
      scale_color_manual(values = c(mmg_colors))+ #+ ggtitle(title) +
      theme_bw()+ 
      scale_y_continuous(limits=c(0,1))+
      
      theme(text = element_text(size=8),
            axis.text=element_text(size=8, face = "bold"),
            axis.title=element_text(size=8,face="bold"),
            legend.position = "none")+ guides(colour = guide_legend(override.aes = list(size=4),nrow=4)) + xlab('Last Mammogram (years)') + ylab('Methylation Score')
    
    png(paste0(name,'.mmg.methscore.',c,'.png'),height = 900,width=1000,res = 300)
    print(plot1)
    dev.off()
    
    
    
    plot1 = ggplot(auc.res1,aes(x = var, y = auc, col = var)) +
      geom_errorbar(aes(ymin=auc.lower, ymax=auc.upper), width=.2,
                    position=position_dodge(0.05))+
      geom_point()+
      scale_y_continuous(limits=c(0,1))+
      scale_color_manual(values = c(mmg_colors))+ #+ ggtitle(title) +
      theme_bw()+ 
      #facet_grid(.~title)+
      
      theme(text = element_text(size=8),
            axis.text=element_text(size=8, face = "bold"),
            axis.title=element_text(size=8,face="bold"),
            legend.position = "none")+ 
      guides(colour = guide_legend(override.aes = list(size=4),nrow=4)) +
      xlab('Last Mammogram (years)') + ylab('CV AUROC (95% CI)')
    
    png(paste0(name,'.filler.',c,'.auroc.ci.mmg.png'),height = 1100,width=1100,res = 300)
    print(plot1)
    dev.off()
    plot1 = ggplot(auc.res1,aes(x = var, y = ci, col = var)) +
      geom_errorbar(aes(ymin=ci.lower, ymax=ci.upper), width=.2,
                    position=position_dodge(0.05))+
      geom_point()+
      scale_y_continuous(limits=c(0,1))+
      scale_color_manual(values = c(mmg_colors))+ #+ ggtitle(title) +
      theme_bw()+ 
      #facet_grid(.~title)+
      theme(text = element_text(size=8),
            axis.text=element_text(size=8, face = "bold"),
            axis.title=element_text(size=8,face="bold"),
            legend.position = "none")+ 
      guides(colour = guide_legend(override.aes = list(size=4),nrow=4)) +
      xlab('Last Mammogram (years)') + ylab('CV Concordance Index (95% CI)')
    
    png(paste0(name,'.filler.',c,'.ci.ci.mmg.png'),height = 1100,width=1100,res = 300)
    print(plot1)
    dev.off()
    plot1 = ggplot(auc.res1,aes(x = var, y =se.spec.95 , col = var)) +
      geom_errorbar(aes(ymin=se.spec.95.lower, ymax=se.spec.95.upper), width=.2,
                    position=position_dodge(0.05))+
      geom_point()+
      scale_y_continuous(limits=c(0,1))+
      scale_color_manual(values = c(mmg_colors))+ #+ ggtitle(title) +
      theme_bw()+ 
      #facet_grid(.~title)+
      theme(text = element_text(size=8),
            axis.text=element_text(size=8, face = "bold"),
            axis.title=element_text(size=8,face="bold"),
            legend.position = "none")+ 
      guides(colour = guide_legend(override.aes = list(size=4),nrow=4)) +
      xlab('Last Mammogram (years)') + ylab('Sensitivity at 95% Specificity (95% CI)')
    
    png(paste0(name,'.filler.',c,'.sens.spec95.ci.mmg.png'),height = 1100,width=1100,res = 300)
    print(plot1)
    dev.off()
    
    plot1 = ggplot(auc.res1,aes(x = var, y = se.spec.90, col = var)) +
      geom_errorbar(aes(ymin=se.spec.90.lower, ymax=se.spec.90.upper), width=.2,
                    position=position_dodge(0.05))+
      geom_point()+
      scale_y_continuous(limits=c(0,1))+
      scale_color_manual(values = c(mmg_colors))+ #+ ggtitle(title) +
      theme_bw()+ 
      #facet_grid(.~title)+
      theme(text = element_text(size=8),
            axis.text=element_text(size=8, face = "bold"),
            axis.title=element_text(size=8,face="bold"),
            legend.position = "none")+ 
      guides(colour = guide_legend(override.aes = list(size=4),nrow=4)) +
      xlab('Last Mammogram (years)')+ ylab('Sensitivity at 90% Specificity (95% CI)')
    
    png(paste0(name,'.filler.',c,'.sens.spec90.ci.mmg.png'),height = 1100,width=1100,res = 300)
    print(plot1)
    dev.off()
    
    plot1 = ggplot(combined.auc.plot.all,aes(x = fpr, y = tpr, col = var)) + geom_line(linewidth=1) +
      scale_color_manual(values = c(mmg_colors))+ #+ ggtitle(title) +
      theme_bw()+ 
      theme(text = element_text(size=8),
            axis.text=element_text(size=8, face = "bold"),
            axis.title=element_text(size=8,face="bold"),
            legend.position = "right")+ guides(colour = guide_legend(override.aes = list(size=4),ncol=1)) + 
      xlab('False Positive Rate') + ylab('Sensitivity')
    
    png(paste0(name,'.filler.',c,'.auroc.mmg.labs.png'),height = 1100,width=1100,res = 300)
    print(plot1)
    dev.off()
    
    
    #boxplot of score
    plot1 = ggplot(pred.df.targ.collapse.annotated,aes(x = mammogram, y = methylation_score, col = mammogram)) +
      geom_boxplot(outlier.shape=NA)+geom_jitter(width=0.1)+
      scale_color_manual(values = c(mmg_colors))+ #+ ggtitle(title) +
      theme_bw()+ 
      scale_y_continuous(limits=c(0,1))+
      
      facet_grid2(. ~ reported)+
      
      theme(text = element_text(size=8),
            axis.text=element_text(size=8, face = "bold"),
            axis.title=element_text(size=8,face="bold"),
            legend.position = "none")+ 
      guides(colour = guide_legend(override.aes = list(size=4),nrow=4)) +
      xlab('Last Mammogram (years)') + ylab('Methylation Score')
    
    png(paste0(name,'.mmg.methscore.',c,'.png'),height = 900,width=2000,res = 300)
    print(plot1)
    dev.off()
    
    #hr type
    hr_information=readRDS('/.mounts/labs/awadallalab/private/ncheng/cfmedip_data/cptp_samples/pathology_records/breast_hr_information.RDS')
    hr_information$HR_status = ifelse(hr_information$ER %in% c('positive','low positive') | hr_information$PR %in% c('positive','low positive'), 'HR+','Not Reported/Tested')
    hr_information$HR_status = ifelse(hr_information$ER %in% c('negative') & hr_information$PR %in% c('negative'), 'HR-',hr_information$HR_status)
    hr_information$HER2_status = ifelse(hr_information$her2_status %in% c('Equivocal','Overexpressed/Amplified'),'HER2+',
                                        ifelse(hr_information$her2_status %in% 'Negative','HER2-', 'Not Reported/Tested'))
    
    hr_information$HR.HER2 = paste0(hr_information$HR_status,'\n',hr_information$HER2_status)
    
    merged.df.all.tmp = merge(merged.df.all.tmp,hr_information[,c('GRP_Id','HR.HER2')],all.x=T)
    merged.df.all.tmp$HR.HER2 = ifelse(merged.df.all.tmp$group == 'Control','Control',merged.df.all.tmp$HR.HER2)
    merged.df.all.tmp$HR.HER2 = ifelse(merged.df.all.tmp$HR.HER2 %in% c('Control','HR+\nHER2-','HR+\nHER2+','HR-\nHER2-','HR-\nHER2+'),merged.df.all.tmp$HR.HER2,'NR')
    
    age0 = merged.df.all.tmp[merged.df.all.tmp$HR.HER2 %in% ('HR+\nHER2-') | merged.df.all.tmp$group =='Control',]
    age1 = merged.df.all.tmp[merged.df.all.tmp$HR.HER2 %in% ('HR+\nHER2+') | merged.df.all.tmp$group =='Control',]
    age2 = merged.df.all.tmp[merged.df.all.tmp$HR.HER2 %in% c('HR-\nHER2-')| merged.df.all.tmp$group =='Control',]
    age3 = merged.df.all.tmp[merged.df.all.tmp$HR.HER2 %in% c('HR-\nHER2+')| merged.df.all.tmp$group =='Control',]
    age4 = merged.df.all.tmp[merged.df.all.tmp$HR.HER2 %in% c('NR')| merged.df.all.tmp$group =='Control',]
    
    
    auc.plot.all.s0 = summary.calc(pred.df.targ.collapse[pred.df.targ.collapse$GRP_Id %in% age0$GRP_Id,],merged.df.all.tmp)
    auc.plot.all.s0[[2]]$var = 'HR+\nHER2-'
    auc.plot.all.s0[[2]]$var.group =  'Subtype'
    
    auc.plot.all.s1 = summary.calc(pred.df.targ.collapse[pred.df.targ.collapse$GRP_Id %in% age1$GRP_Id,],merged.df.all.tmp)
    auc.plot.all.s1[[2]]$var = 'HR+\nHER2+'
    auc.plot.all.s1[[2]]$var.group =  'Subtype'
    
    auc.plot.all.s2 = summary.calc(pred.df.targ.collapse[pred.df.targ.collapse$GRP_Id %in% age2$GRP_Id,],merged.df.all.tmp)
    auc.plot.all.s2[[2]]$var = 'HR-\nHER2-'
    auc.plot.all.s2[[2]]$var.group =  'Subtype'
    
    
    #auc.plot.all.s3 = summary.calc(pred.df.targ.collapse[pred.df.targ.collapse$GRP_Id %in% age3$GRP_Id,],merged.df.all.tmp)
    #auc.plot.all.s3[[2]]$var = 'HR-\nHER2+'
    #auc.plot.all.s3[[2]]$var.group = 'Subtype'
    
    auc.plot.all.s4 = summary.calc(pred.df.targ.collapse[pred.df.targ.collapse$GRP_Id %in% age4$GRP_Id,],merged.df.all.tmp)
    auc.plot.all.s4[[2]]$var = 'NR'
    auc.plot.all.s4[[2]]$var.group =  'Subtype'
    #annotating auc t
    auc.plot.all = summary.calc(pred.df.targ.collapse,merged.df.all.tmp)
    auc.plot.all[[2]]$stage = 'All'
    auc.plot.all[[2]]$var = 'All'
    auc.plot.all[[2]]$var.group = 'Subtype'
    
    auc.plot.all[[3]]$var = 'All'
    
    
    auc.res1 = rbind(
      auc.plot.all.s0[[2]][1,],
      auc.plot.all.s1[[2]][1,],
      auc.plot.all.s2[[2]][1,],
      auc.plot.all.s4[[2]][1,])
    auc.res1$title = 'Subtype'
    #combined.auroc = combined.auroc[,colnames(combined.auroc) %in% colnames(auc.res1)]
    combined.auroc = rbind(combined.auroc,auc.res1[,colnames(auc.res1) %in% colnames(combined.auroc)])
    
    subtype_colors =  c('HR+\nHER2+' = '#8FD5A6','HR+\nHER2-' = '#BF1363',"HR-\nHER2+" = '#E09F3E','HR-\nHER2-'='#545E75','NR' ='#66462C','Control' = 'black')
    combined.auc.plot.all = rbind(
      auc.plot.all.s0[[2]],
      auc.plot.all.s1[[2]],
      auc.plot.all.s2[[2]],
      auc.plot.all.s4[[2]])#,
    # auc.plot.all.s3[[2]])
    
    
    #annotating auc t
    auc.plot.all[[3]]$var = 'All'
    auc.plot.all.s0[[3]]$var = 'HR+\nHER2-'
    auc.plot.all.s1[[3]]$var = 'HR+\nHER2+'
    auc.plot.all.s2[[3]]$var = 'HR-\nHER2-'
    auc.plot.all.s4[[3]]$var = 'NR'
    
    
    
    auc.t =rbind( auc.plot.all[[3]],
                  auc.plot.all.s0[[3]],
                  
                  auc.plot.all.s1[[3]],
                  auc.plot.all.s2[[3]],
                  auc.plot.all.s4[[3]]
                  
                  
                  
    )
    auc.t$title = 'Subtype'
    auc.t.combined = rbind(auc.t.combined,auc.t)
    
    #auc t
    subtype_colors.groups =  c('HR+\nHER2+' = '#8FD5A6','HR+\nHER2-' = '#BF1363',"HR-\nHER2+" = '#E09F3E','HR-\nHER2-'='#545E75','NR' ='#66462C','All' = 'black')
    
    plot1 = ggplot(auc.t[auc.t$time >= 365*1,],aes(x = time/365, y = AUC, col = var)) + geom_line(linewidth=1) +
      scale_color_manual(values = c(subtype_colors.groups))+ #+ ggtitle(title) +
      scale_y_continuous(limits=c(0,1))+
      theme_bw()+ 
      #facet_grid2(.~title)+
      theme(text = element_text(size=8),
            axis.text=element_text(size=8, face = "bold"),
            axis.title=element_text(size=8,face="bold"),
            legend.position = "none")+ guides(colour = guide_legend(override.aes = list(size=4),nrow=4)) +
      xlab('Time (Years)') + ylab('Time-Dependent CV AUROC')
    
    png(paste0(name,'.filler.',c,'.timeauroc.subtype.png'),height = 1100,width=1100,res = 300)
    print(plot1)
    dev.off()
    
    
    
    plot1 = ggplot(combined.auc.plot.all,aes(x = fpr, y = tpr, col = var)) + geom_line(linewidth=1) +
      scale_color_manual(values = c(subtype_colors))+ #+ ggtitle(title) +
      theme_bw()+ 
      theme(text = element_text(size=8),
            axis.text=element_text(size=8, face = "bold"),
            axis.title=element_text(size=8,face="bold"),
            legend.position = "none")+ guides(colour = guide_legend(override.aes = list(size=4),nrow=4)) + xlab('False Positive Rate') + ylab('Sensitivity')
    
    png(paste0(name,'.filler.',c,'.auroc.subtype.png'),height = 1000,width=1000,res = 300)
    print(plot1)
    dev.off()
    
    plot1 = ggplot(combined.auc.plot.all,aes(x = fpr, y = tpr, col = var)) + geom_line(linewidth=1) +
      scale_color_manual(values = c(subtype_colors))+ #+ ggtitle(title) +
      theme_bw()+ 
      theme(text = element_text(size=8),
            axis.text=element_text(size=8, face = "bold"),
            axis.title=element_text(size=8,face="bold"),
            legend.position = "right")+ guides(colour = guide_legend(override.aes = list(size=4),nrow=4)) + xlab('False Positive Rate') + ylab('Sensitivity')
    
    png(paste0(name,'.filler.',c,'.auroc.subtype.labs.png'),height = 1000,width=1000,res = 300)
    print(plot1)
    dev.off()
    pred.df.targ.collapse.annotated = merge(pred.df.targ.collapse, hr_information[,c('GRP_Id','HR.HER2')],all.x=T)
    pred.df.targ.collapse.annotated$HR.HER2 = ifelse(pred.df.targ.collapse.annotated$reported == 'Control','Control',pred.df.targ.collapse.annotated$HR.HER2)
    pred.df.targ.collapse.annotated$HR.HER2 = ifelse(pred.df.targ.collapse.annotated$HR.HER2 %in% c('Control','HR+\nHER2-','HR+\nHER2+','HR-\nHER2-','HR-\nHER2+'),pred.df.targ.collapse.annotated$HR.HER2,'NR')
    
    plot1 = ggplot(pred.df.targ.collapse.annotated[pred.df.targ.collapse.annotated$HR.HER2 %in% c(names(subtype_colors),'Control'),],aes(x = HR.HER2, y = methylation_score, col = HR.HER2)) +
      geom_boxplot(outlier.shape=NA)+geom_jitter(width=0.1)+
      scale_color_manual(values = c(subtype_colors))+ #+ ggtitle(title) +
      theme_bw()+ 
      scale_y_continuous(limits=c(0,1))+
      
      theme(text = element_text(size=8),
            axis.text=element_text(size=8, face = "bold"),
            axis.title=element_text(size=8,face="bold"),
            legend.position = "none")+ guides(colour = guide_legend(override.aes = list(size=4),nrow=4)) +
      xlab('Subtype') + ylab('Methylation Score')
    
    png(paste0(name,'.subtype.methscore.',c,'.png'),height = 900,width=1100,res = 300)
    print(plot1)
    dev.off()
    
    
    
    plot1 = ggplot(auc.res1,aes(x = var, y = auc, col = var)) +
      geom_errorbar(aes(ymin=auc.lower, ymax=auc.upper), width=.2,
                    position=position_dodge(0.05))+
      geom_point()+
      scale_y_continuous(limits=c(0,1))+
      scale_color_manual(values = c(subtype_colors))+ #+ ggtitle(title) +
      theme_bw()+ 
      #facet_grid(.~title)+
      
      theme(text = element_text(size=8),
            axis.text=element_text(size=8, face = "bold"),
            axis.title=element_text(size=8,face="bold"),
            legend.position = "none")+ 
      guides(colour = guide_legend(override.aes = list(size=4),nrow=4)) +
      xlab('Subtype') + ylab('CV AUROC (95% CI)')
    
    png(paste0(name,'.filler.',c,'.auroc.ci.subtype.png'),height = 1100,width=1500,res = 300)
    print(plot1)
    dev.off()
    plot1 = ggplot(auc.res1,aes(x = var, y = ci, col = var)) +
      geom_errorbar(aes(ymin=ci.lower, ymax=ci.upper), width=.2,
                    position=position_dodge(0.05))+
      geom_point()+
      scale_y_continuous(limits=c(0,1))+
      scale_color_manual(values = c(subtype_colors))+ #+ ggtitle(title) +
      theme_bw()+ 
      #facet_grid(.~title)+
      theme(text = element_text(size=8),
            axis.text=element_text(size=8, face = "bold"),
            axis.title=element_text(size=8,face="bold"),
            legend.position = "none")+ 
      guides(colour = guide_legend(override.aes = list(size=4),nrow=4)) +
      xlab('Subtype') + ylab('CV Concordance Index (95% CI)')
    
    png(paste0(name,'.filler.',c,'.ci.ci.subtype.png'),height = 1100,width=1500,res = 300)
    print(plot1)
    dev.off()
    plot1 = ggplot(auc.res1,aes(x = var, y =se.spec.95 , col = var)) +
      geom_errorbar(aes(ymin=se.spec.95.lower, ymax=se.spec.95.upper), width=.2,
                    position=position_dodge(0.05))+
      geom_point()+
      scale_y_continuous(limits=c(0,1))+
      scale_color_manual(values = c(subtype_colors))+ #+ ggtitle(title) +
      theme_bw()+ 
      #facet_grid(.~title)+
      theme(text = element_text(size=8),
            axis.text=element_text(size=8, face = "bold"),
            axis.title=element_text(size=8,face="bold"),
            legend.position = "none")+ 
      guides(colour = guide_legend(override.aes = list(size=4),nrow=4)) +
      xlab('Subtype') + ylab('Sensitivity at 95% Specificity (95% CI)')
    
    png(paste0(name,'.filler.',c,'.sens.spec95.ci.subtype.png'),height = 1100,width=1500,res = 300)
    print(plot1)
    dev.off()
    
    plot1 = ggplot(auc.res1,aes(x = var, y = se.spec.90, col = var)) +
      geom_errorbar(aes(ymin=se.spec.90.lower, ymax=se.spec.90.upper), width=.2,
                    position=position_dodge(0.05))+
      geom_point()+
      scale_y_continuous(limits=c(0,1))+
      scale_color_manual(values = c(subtype_colors))+ #+ ggtitle(title) +
      theme_bw()+ 
      #facet_grid(.~title)+
      theme(text = element_text(size=8),
            axis.text=element_text(size=8, face = "bold"),
            axis.title=element_text(size=8,face="bold"),
            legend.position = "none")+ 
      guides(colour = guide_legend(override.aes = list(size=4),nrow=4)) +
      xlab('Subtype')+ ylab('Sensitivity at 90% Specificity (95% CI)')
    
    png(paste0(name,'.filler.',c,'.sens.spec90.ci.subtype.png'),height = 1100,width=1500,res = 300)
    print(plot1)
    dev.off()
    
    plot1 = ggplot(combined.auc.plot.all,aes(x = fpr, y = tpr, col = var)) + geom_line(linewidth=1) +
      scale_color_manual(values = c(subtype_colors))+ #+ ggtitle(title) +
      theme_bw()+ 
      theme(text = element_text(size=8),
            axis.text=element_text(size=8, face = "bold"),
            axis.title=element_text(size=8,face="bold"),
            legend.position = "right")+ guides(colour = guide_legend(override.aes = list(size=4),ncol=1)) + 
      xlab('False Positive Rate') + ylab('Sensitivity')
    
    png(paste0(name,'.filler.',c,'.auroc.subtype.labs.png'),height = 1100,width=1500,res = 300)
    print(plot1)
    dev.off()
    
    
    #boxplot of score
    plot1 = ggplot(pred.df.targ.collapse.annotated[pred.df.targ.collapse.annotated$HR.HER2 %in% c(names(subtype_colors),'Control'),],aes(x = HR.HER2, y = methylation_score, col = HR.HER2)) +
      geom_boxplot(outlier.shape=NA)+geom_jitter(width=0.1)+
      scale_color_manual(values = c(subtype_colors))+ #+ ggtitle(title) +
      theme_bw()+ 
      scale_y_continuous(limits=c(0,1))+
      
      theme(text = element_text(size=8),
            axis.text=element_text(size=8, face = "bold"),
            axis.title=element_text(size=8,face="bold"),
            legend.position = "none")+ 
      guides(colour = guide_legend(override.aes = list(size=4),nrow=4)) +
      xlab('Subtype') + ylab('Methylation Score')
    
    png(paste0(name,'.hr.methscore.',c,'.png'),height = 1100,width=1500,res = 300)
    print(plot1)
    dev.off()
    #
    
    targ.plots.spec = unique(combined.auroc[combined.auroc$var == 'All',c('jcutpoint.spec.weighted',
                                                                          'jcutpoint.spec.lower.weighted',
                                                                          'jcutpoint.spec.upper.weighted')])
    colnames(targ.plots.spec) = c('value','lower','upper')
    targ.plots.spec$stat = 'Specificity'
    targ.plots.npv = unique(combined.auroc[combined.auroc$var == 'All',c('jcutpoint.npv.weighted',
                                                                         'jcutpoint.npv.lower.weighted',
                                                                         'jcutpoint.npv.upper.weighted')])
    colnames(targ.plots.npv) = c('value','lower','upper')
    targ.plots.npv$stat = 'NPV'
    
    targ.plots.sens = unique(combined.auroc[combined.auroc$var == 'All',c('jcutpoint.sens.weighted',
                                                                          'jcutpoint.sens.lower.weighted',
                                                                          'jcutpoint.sens.upper.weighted')])
    colnames(targ.plots.sens) = c('value','lower','upper')
    targ.plots.sens$stat = 'Sensitivity'
    
    targ.plots.ppv = unique(combined.auroc[combined.auroc$var == 'All',c('jcutpoint.ppv.weighted',
                                                                         'jcutpoint.ppv.lower.weighted',
                                                                         'jcutpoint.ppv.upper.weighted')])
    colnames(targ.plots.ppv) = c('value','lower','upper')
    targ.plots.ppv$stat = 'PPV'
    
    combined.targ = rbind(targ.plots.spec,targ.plots.npv,targ.plots.sens, targ.plots.ppv)
    
    
    plot1 = ggplot(combined.targ[!duplicated(combined.targ$stat),],aes(x = stat, y = value)) +
      geom_errorbar(aes(ymin=lower, ymax=upper), width=.2,
                    position=position_dodge(0.05))+
      geom_point()+
      scale_y_continuous(limits=c(0,1))+
      #scale_color_manual(values = c(dxage_colors))+ #+ ggtitle(title) +
      theme_bw()+ 
      #facet_grid(.~title)+
      theme(text = element_text(size=8),
            axis.text=element_text(size=8, face = "bold"),
            axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
            axis.title.y=element_text(size=8,face="bold"),
            axis.title.x=element_blank(),
            legend.position = "none")+ 
      guides(colour = guide_legend(override.aes = list(size=4),nrow=4)) +
      xlab('')+ ylab('Value')
    
    png(paste0(name,'.filler.',c,'.npv.spec.youdenoverall.png'),height = 1100,width=800,res = 300)
    print(plot1)
    dev.off()
    
    png(paste0(name,'.filler.',c,'.npv.spec.youdenoverall.short.png'),height = 600,width=800,res = 300)
    print(plot1)
    dev.off()
    #specificty/npv plot youdens index cutoff
    targ.plots.spec = unique(combined.auroc[combined.auroc$var == 'All',c('f1cutpoint.spec',
                                                                          'f1cutpoint.spec.lower',
                                                                          'f1cutpoint.spec.upper')])
    colnames(targ.plots.spec) = c('value','lower','upper')
    targ.plots.spec$stat = 'Specificity'
    targ.plots.npv = unique(combined.auroc[combined.auroc$var == 'All',c('f1cutpoint.npv',
                                                                         'f1cutpoint.npv.lower',
                                                                         'f1cutpoint.npv.upper')])
    colnames(targ.plots.npv) = c('value','lower','upper')
    targ.plots.npv$stat = 'NPV'
    
    combined.targ = rbind(targ.plots.spec,targ.plots.npv)
    
    
    plot1 = ggplot(combined.targ,aes(x = stat, y = value)) +
      geom_errorbar(aes(ymin=lower, ymax=upper), width=.2,
                    position=position_dodge(0.05))+
      geom_point()+
      scale_y_continuous(limits=c(0,1))+
      #scale_color_manual(values = c(dxage_colors))+ #+ ggtitle(title) +
      theme_bw()+ 
      #facet_grid(.~title)+
      theme(text = element_text(size=8),
            axis.text=element_text(size=8, face = "bold"),
            axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
            axis.title.y=element_text(size=8,face="bold"),
            axis.title.x=element_blank(),
            legend.position = "none")+ 
      guides(colour = guide_legend(override.aes = list(size=4),nrow=4)) +
      xlab('')+ ylab('Value')
    
    png(paste0(name,'.filler.',c,'.npv.spec.f1overall.png'),height = 1100,width=600,res = 300)
    print(plot1)
    dev.off()
    
    png(paste0(name,'.filler.',c,'.npv.spec.f1overall.short.png'),height = 600,width=800,res = 300)
    print(plot1)
    dev.off()
    
    
    #km curves.all#
    library(survminer)
    pca.coords =merge(pred.df.targ.collapse.annotated,merged.df.all.tmp) # merge(pred.df.targ.collapse.annotated, merged.df.all.tmp
    pca.coords =merge(pred.df.targ.collapse,merged.df.all.tmp) # merge(pred.df.targ.collapse.annotated, merged.df.all.tmp
    
    pca.coords$median.g = factor(ifelse(pca.coords$methylation_score > score.cutoff,'Higher Test Median','Lower Test Median'),levels=c('Lower Test Median','Higher Test Median'))
    
    pca.coords$Cancer.g= ifelse(pca.coords$reported == 'Control',0,1)
    female.weights= weightsf.females(pca.coords)
    cph = summary(coxph(Surv(censorship_time, Cancer.g) ~ median.g, data =  pca.coords,weights=female.weights))$logtest[3]
    hr = round(summary(coxph(Surv(censorship_time, Cancer.g) ~ median.g, data =  pca.coords,weights=female.weights))$coefficients[2],digits=3)
    median.g = pca.coords$median.g
    a=coxph(formula = Surv(censorship_time, Cancer.g) ~ median.g, data =  pca.coords,weights=female.weights)
    b = survfit(a,newdata=data.frame(365*seq(1:5)))
    test =T
    if (test == T){      
      score.cutoff.val = score.cutoff
      pca.coords$median.g = factor(ifelse(pca.coords$methylation_score > score.cutoff.val,'Higher Test Median','Lower Test Median'),levels=c('Lower Test Median','Higher Test Median'))
      
      pca.coords$Cancer.g= ifelse(pca.coords$reported == 'Control',0,1)
      female.weights= weightsf.females(pca.coords)
      cph = summary(coxph(Surv(censorship_time, Cancer.g) ~ median.g, data =  pca.coords,weights=female.weights))$logtest[3]
      hr = round(summary(coxph(Surv(censorship_time, Cancer.g) ~ median.g, data =  pca.coords,weights=female.weights))$coefficients[2],digits=3)
      median.g = pca.coords$median.g
      if (dx.all == T) {
        max.time = 10
      } else {
        max.time=5
      }
      a=coxph(formula = Surv(censorship_time, Cancer.g) ~ median.g, data =  pca.coords,weights=female.weights)
      b = survfit(a,newdata=data.frame(365*seq(1:max.time)))
      plot4 = ggsurvplot(
        fit = survfit(Surv(censorship_time/365, Cancer.g) ~ median.g, data =  pca.coords,weights = female.weights), 
        risk.table = F,
        #cumevents = TRUE,
        palette = c('#2a9d8f','#e9c46a'),
        legend.labs=c("Higher Test Median","Lower Test Median"),
        xlim = c(0,10), xlab = "Time (Years)", ylab = c("Survival probability"), break.x.by = 1,
        ylim = c(0.85,1),
        alpha = 0.6,
        ggtheme = theme_bw(),
        font.tickslab = c(8,'bold'), 
        font.main = c(8, "bold",'black'),
        font.x = c(8, "bold",'black'),
        font.y = c(8, "bold",'black'),
        font.legend = c(8)) + 
        ggtitle(paste(hr,'.',formatC(cph, format = "e", digits = 2))) #+ theme(legend.position = 'none')
      
      png(paste0(name,".scores.surv_median.all.png"),height = 1100, width = 1100,res=300)
      print(plot4)
      dev.off()
      
      plot4 = ggsurvplot(
        fit = survfit(Surv(censorship_time/365, Cancer.g) ~ median.g, data =  pca.coords,weights = female.weights), 
        risk.table = F,
        #cumevents = TRUE,
        palette = c('#2a9d8f','#e9c46a'),
        legend.labs=c("Higher Test Median","Lower Test Median"),
        xlim = c(0,10), xlab = "Time (Years)", ylab = c("Survival probability"), break.x.by = 1,
        ylim = c(0,1),
        ggtheme = theme_bw(),
        
        font.tickslab = c(8,'bold'), 
        font.main = c(8, "bold",'black'),
        font.x = c(8, "bold",'black'),
        font.y = c(8, "bold",'black'),
        alpha = 0.6,
        font.legend = c(8), risk.table.fontsize = 8, risk.table.col = "black", size = 1) + ggtitle(paste(hr,'.',formatC(cph, format = "e", digits = 2)))   #+ theme(legend.position = 'none')
      
      png(paste0(name,".scores.surv_median.all.full.png"),height = 1100, width = 1100,res=300)
      print(plot4)
      dev.off()
      
      pca.coords$quartile <- dplyr::ntile(pca.coords$methylation_score, 4) 
      quartile = pca.coords$quartile 
      cph = summary(coxph(Surv(censorship_time, Cancer.g) ~ quartile, data =  pca.coords,weights = female.weights))$logtest[3]
      hr = round(summary(coxph(Surv(censorship_time, Cancer.g) ~ quartile, data =  pca.coords,weights=female.weights))$coefficients[2],digits=3)
      survfit.obj = survfit(Surv(censorship_time/365, Cancer.g) ~ quartile, data =  pca.coords,weights = female.weights)
      quart.col = c('#2a9d8f','#e9c46a','#EA526F','#586994')
      names(quart.col) = c("1","2","3","4")
      quart.col.plot = quart.col[names(quart.col) %in%  pca.coords$quartile ]
      plot4 = ggsurvplot(
        fit = survfit.obj, 
        risk.table = F,
        #cumevents = TRUE,
        palette = quart.col.plot,
        legend.labs=names(quart.col.plot),
        xlim = c(0,10), xlab = "Time (Years)", ylab = c("Survival probability"), break.x.by = 1,
        ylim=c(0.8,1),
        ggtheme = theme_bw(),
        
        font.tickslab = c(8,'bold'), 
        font.main = c(8, "bold",'black'),
        font.x = c(8, "bold",'black'),
        font.y = c(8, "bold",'black'),
        alpha = 0.6,
        font.legend = c(8), risk.table.fontsize = 5, risk.table.col = "black", size = 1)  + 
        ggtitle(paste0(hr,'.',formatC(cph, format = "e", digits = 2)) )  #+ the #+ theme(legend.position = 'none')
      
      
      png(paste0(name,".scores.surv_quart.all.cut.png"),height = 1100, width = 1100,res=300)
      print(plot4)
      dev.off()
      
      pca.coords$quartile <- ifelse(pca.coords$methylation_score >=0 & pca.coords$methylation_score < 0.25, c(1),
                                    ifelse(pca.coords$methylation_score >=0.25 & pca.coords$methylation_score < 0.5, c(2),
                                           ifelse(pca.coords$methylation_score >=0.5 & pca.coords$methylation_score < 0.75, c(3),
                                                  ifelse(pca.coords$methylation_score >=0.75 & pca.coords$methylation_score < 1, c(4),'Other'))))
      quartile = pca.coords$quartile 
      cph = summary(coxph(Surv(censorship_time, Cancer.g) ~ quartile, data =  pca.coords,weights = female.weights))$logtest[3]
      hr = round(summary(coxph(Surv(censorship_time, Cancer.g) ~ quartile, data =  pca.coords,weights=female.weights))$coefficients[2],digits=3)
      survfit.obj = survfit(Surv(censorship_time/365, Cancer.g) ~ quartile, data =  pca.coords,weights = female.weights)
      quart.col = c('#2a9d8f','#e9c46a','#EA526F','#586994')
      names(quart.col) = c("1","2","3","4")
      quart.col.plot = quart.col[names(quart.col) %in%  pca.coords$quartile ]
      plot4 = ggsurvplot(
        fit = survfit.obj, 
        risk.table = F,
        #cumevents = TRUE,
        palette = quart.col.plot,
        legend.labs=names(quart.col.plot),
        xlim = c(0,10), xlab = "Time (Years)", ylab = c("Survival probability"), break.x.by = 1,
        ylim=c(0.85,1),
        ggtheme = theme_bw(),
        
        font.tickslab = c(8,'bold'), 
        font.main = c(8, "bold",'black'),
        font.x = c(8, "bold",'black'),
        font.y = c(8, "bold",'black'),
        alpha = 0.6,
        font.legend = c(8), risk.table.fontsize = 5, risk.table.col = "black", size = 1)  + 
        ggtitle(paste0(hr,'.',formatC(cph, format = "e", digits = 2)) )  #+ the #+ theme(legend.position = 'none')
      
      
      png(paste0(name,".scores.surv_quart.breaks.cut.png"),height = 1100, width = 1100,res=300)
      print(plot4)
      dev.off()
      
      
      
      plot4 = ggsurvplot(
        fit = survfit.obj, 
        risk.table = F,
        #cumevents = TRUE,
        palette = quart.col.plot,
        legend.labs=names(quart.col.plot),
        xlim = c(0,10), xlab = "Time (Years)", ylab = c("Survival probability"), break.x.by = 1,
        ylim=c(0,1),
        ggtheme = theme_bw(),
        
        font.tickslab = c(8,'bold'), 
        font.main = c(8, "bold",'black'),
        font.x = c(8, "bold",'black'),
        font.y = c(8, "bold",'black'),
        alpha = 0.6,
        font.legend = c(8), risk.table.fontsize = 5, risk.table.col = "black", size = 1)  + 
        ggtitle(paste0(hr,'.',formatC(cph, format = "e", digits = 2)) )  #+ the #+ theme(legend.position = 'none')
      
      png(paste0(name,".scores.surv_quart.all.png"),height = 1100, width = 1100,res=300)
      print(plot4)
      dev.off() 
      
      
      km=survfit(formula = Surv(censorship_time, Cancer.g) ~ median.g, data =  pca.coords,weights=female.weights)
      tmp = summary(km,times=c(1:3650))
      km.summary.weighted.test = data.frame(time = tmp$time, surv.prob =tmp$surv, surv.prob.l =tmp$lower,surv.prob.h = tmp$upper,risk = gsub('risk=','',tmp$strata))
      colnames(km.summary.weighted.test)[2:4]= c('surv','lower','upper')
      #combined.validation.surv = rbind(km.summary.weighted.test,surv.prob.validation[,-3])
      plot1= ggplot(km.summary.weighted.test) +
        geom_line(aes(x = time/365, y = 100-surv*100, col = risk),size = 1,alpha= 0.7) + 
        
        scale_fill_manual(values = c('#B9F3EC','#FDF2C4'))+
        scale_color_manual(values = c('#2a9d8f','#e9c46a'))+
        scale_x_continuous(breaks = seq(0,max.time,1),limits = c(0,max.time))+
        theme_bw()+
        theme(text = element_text(size=8),
              axis.text=element_text(size=8, face = "bold"),
              axis.title=element_text(size=8,face="bold"),
              legend.position = "none")+
        #scale_y_continuous(limits = c(0,1), breaks = seq(0,1,0.2)) +
        ylab('Breast Cancer Probability (%)') + xlab('Time (Years)') 
      png(paste0(name,".scores.weighted.risk.all.all.png"),height = 900, width = 1100,res=300)
      print(plot1)
      dev.off()
      
      for (dxage in c('30-50','50-60','60-70','70+')) {
        pca.coords.filt =  pca.coords[as.character(pca.coords$age_group3) %in% c(dxage) | pca.coords$Cancer == 'Control',]
        female.weights.filt = female.weights[as.character(pca.coords$age_group3) %in% c(dxage)| pca.coords$Cancer == 'Control']
        
        
        cph = summary(coxph(Surv(censorship_time, Cancer.g) ~ median.g, data = pca.coords.filt,weights=female.weights.filt))$logtest[3]
        hr = round(summary(coxph(Surv(censorship_time, Cancer.g) ~ median.g, data =  pca.coords.filt,weights=female.weights.filt))$coefficients[2],digits=3)
        median.g = pca.coords.filt$median.g
        a=coxph(formula = Surv(censorship_time, Cancer.g) ~ median.g, data =  pca.coords.filt,weights=female.weights.filt)
        b = survfit(a,newdata=data.frame(365*seq(1:5)))
        
        if (length(unique(pca.coords.filt$median.g)) > 1) {
          plot4 = ggsurvplot(
            fit = survfit(Surv(censorship_time/365, Cancer.g) ~ median.g, data =  pca.coords.filt,weights = female.weights.filt), 
            risk.table = F,
            #cumevents = TRUE,
            palette = c('#2a9d8f','#e9c46a'),
            legend.labs=c("Higher Test Median","Lower Test Median"),
            xlim = c(0,10), xlab = "Time (Years)", ylab = c("Survival probability"), break.x.by = 1,
            ylim = c(0.9,1),
            ggtheme = theme_bw(),
            
            font.tickslab = c(8,'bold'), 
            font.main = c(8, "bold",'black'),
            font.x = c(8, "bold",'black'),
            font.y = c(8, "bold",'black'),
            alpha = 0.6,
            font.legend = c(8), risk.table.fontsize = 5, risk.table.col = "black", size = 1) + 
            ggtitle(paste(hr,'.',formatC(cph, format = "e", digits = 2))) #+ theme(legend.position = 'none')
          
          png(paste0(name,".scores.surv_median.dxage.",dxage,".cut.png"),height = 1100, width = 1100,res=300)
          print(plot4)
          dev.off()
          
          plot4 = ggsurvplot(
            fit = survfit(Surv(censorship_time/365, Cancer.g) ~ median.g, data =  pca.coords.filt,weights = female.weights.filt), 
            risk.table = F,
            #cumevents = TRUE,
            palette = c('#2a9d8f','#e9c46a'),
            legend.labs=c("Higher Test Median","Lower Test Median"),
            xlim = c(0,10), xlab = "Time (Years)", ylab = c("Survival probability"), break.x.by = 1,
            ylim = c(0,1),
            ggtheme = theme_bw(),
            
            font.tickslab = c(8,'bold'), 
            font.main = c(8, "bold",'black'),
            font.x = c(8, "bold",'black'),
            font.y = c(8, "bold",'black'),
            alpha = 0.6,
            font.legend = c(8), risk.table.fontsize = 8, risk.table.col = "black", size = 1) + ggtitle(paste(hr,'.',formatC(cph, format = "e", digits = 2)))   #+ theme(legend.position = 'none')
          
          png(paste0(name,".scores.surv_median.dxage.",dxage,".all.png"),height = 1100, width = 1100,res=300)
          print(plot4)
          dev.off()
          
          
          pca.coords.filt$quartile <- dplyr::ntile(pca.coords.filt$methylation_score, 4) 
          quartile = pca.coords.filt$quartile 
          cph = summary(coxph(Surv(censorship_time, Cancer.g) ~ quartile, data =  pca.coords.filt,weights = female.weights.filt))$logtest[3]
          hr = round(summary(coxph(Surv(censorship_time, Cancer.g) ~ quartile, data =  pca.coords.filt,weights=female.weights.filt))$coefficients[2],digits=3)
          survfit.obj = survfit(Surv(censorship_time/365, Cancer.g) ~ quartile, data =  pca.coords.filt,weights = female.weights.filt)
          quart.col = c('#2a9d8f','#e9c46a','#EA526F','#586994')
          names(quart.col) = c("1","2","3","4")
          quart.col.plot = quart.col[names(quart.col) %in%  pca.coords.filt$quartile ]
          plot4 = ggsurvplot(
            fit = survfit.obj, 
            risk.table = F,
            #cumevents = TRUE,
            palette = quart.col.plot,
            legend.labs=names(quart.col.plot),
            xlim = c(0,10), xlab = "Time (Years)", ylab = c("Survival probability"), break.x.by = 1,
            ylim=c(0.9,1),
            ggtheme = theme_bw(),
            
            font.tickslab = c(8,'bold'), 
            font.main = c(8, "bold",'black'),
            font.x = c(8, "bold",'black'),
            font.y = c(8, "bold",'black'),  alpha = 0.6,
            font.legend = c(8), risk.table.fontsize = 5, risk.table.col = "black", size = 1)  + 
            ggtitle(paste0(hr,'.',formatC(cph, format = "e", digits = 2)) )  #+ the #+ theme(legend.position = 'none')
          
          
          png(paste0(name,".scores.surv_quart.dxage.",dxage,".png"),height = 1100, width = 1100,res=300)
          print(plot4)
          dev.off()
          
          
          plot4 = ggsurvplot(
            fit = survfit.obj, 
            risk.table = F,
            #cumevents = TRUE,
            palette = quart.col.plot,
            legend.labs=names(quart.col.plot),
            xlim = c(0,10), xlab = "Time (Years)", ylab = c("Survival probability"), break.x.by = 1,
            ylim=c(0,1),
            ggtheme = theme_bw(),
            
            font.tickslab = c(8,'bold'), 
            font.main = c(8, "bold",'black'),
            font.x = c(8, "bold",'black'),
            font.y = c(8, "bold",'black'), alpha = 0.6,
            font.legend = c(8), risk.table.fontsize = 5, risk.table.col = "black", size = 1)  + 
            ggtitle(paste0(hr,'.',formatC(cph, format = "e", digits = 2)) )  #+ the #+ theme(legend.position = 'none')
          
          png(paste0(name,".scores.surv_quart.dxage.",dxage,".full.png"),height = 1100, width = 1100,res=300)
          print(plot4)
          dev.off()
          
          pca.coords$quartile <- ifelse(pca.coords$methylation_score >=0 & pca.coords$methylation_score < 0.25, c(1),
                                        ifelse(pca.coords$methylation_score >=0.25 & pca.coords$methylation_score < 0.5, c(2),
                                               ifelse(pca.coords$methylation_score >=0.5 & pca.coords$methylation_score < 0.75, c(3),
                                                      ifelse(pca.coords$methylation_score >=0.75 & pca.coords$methylation_score < 1, c(4),'Other'))))
          quartile = pca.coords.filt$quartile 
          cph = summary(coxph(Surv(censorship_time, Cancer.g) ~ quartile, data =  pca.coords.filt,weights = female.weights.filt))$logtest[3]
          hr = round(summary(coxph(Surv(censorship_time, Cancer.g) ~ quartile, data =  pca.coords.filt,weights=female.weights.filt))$coefficients[2],digits=3)
          survfit.obj = survfit(Surv(censorship_time/365, Cancer.g) ~ quartile, data =  pca.coords.filt,weights = female.weights.filt)
          quart.col = c('#2a9d8f','#e9c46a','#EA526F','#586994')
          names(quart.col) = c("1","2","3","4")
          quart.col.plot = quart.col[names(quart.col) %in%  pca.coords.filt$quartile ]
          plot4 = ggsurvplot(
            fit = survfit.obj, 
            risk.table = F,
            #cumevents = TRUE,
            palette = quart.col.plot,
            legend.labs=names(quart.col.plot),
            xlim = c(0,10), xlab = "Time (Years)", ylab = c("Survival probability"), break.x.by = 1,
            ylim=c(0.9,1),
            ggtheme = theme_bw(),
            
            font.tickslab = c(8,'bold'), 
            font.main = c(8, "bold",'black'),
            font.x = c(8, "bold",'black'),
            font.y = c(8, "bold",'black'),
            alpha = 0.6,
            font.legend = c(8), risk.table.fontsize = 5, risk.table.col = "black", size = 1)  + 
            ggtitle(paste0(hr,'.',formatC(cph, format = "e", digits = 2)) )  #+ the #+ theme(legend.position = 'none')
          
          
          png(paste0(name,".scores.surv_quart.fixed.dxage.",dxage,".png"),height = 1100, width = 1100,res=300)
          print(plot4)
          dev.off()
          
          
          km=survfit(formula = Surv(censorship_time, Cancer.g) ~ median.g, data =  pca.coords.filt,weights=female.weights.filt)
          tmp = summary(km,times=c(1:3650))
          km.summary.weighted.test = data.frame(time = tmp$time, surv.prob =tmp$surv, surv.prob.l =tmp$lower,surv.prob.h = tmp$upper,risk = gsub('risk=','',tmp$strata))
          colnames(km.summary.weighted.test)[2:4]= c('surv','lower','upper')
          #combined.validation.surv = rbind(km.summary.weighted.test,surv.prob.validation[,-3])
          plot1= ggplot(km.summary.weighted.test) +
            geom_line(aes(x = time/365, y = 100-surv*100, col = risk),size = 1,alpha= 0.7) + 
            
            scale_fill_manual(values = c('#B9F3EC','#FDF2C4'))+
            scale_color_manual(values = c('#2a9d8f','#e9c46a'))+
            scale_x_continuous(breaks = seq(0,max.time,1),limits = c(0,max.time))+
            theme_bw()+
            theme(text = element_text(size=8),
                  axis.text=element_text(size=8, face = "bold"),
                  axis.title=element_text(size=8,face="bold"),
                  legend.position = "none")+
            #scale_y_continuous(limits = c(0,1), breaks = seq(0,1,0.2)) +
            ylab('Breast Cancer Probability (%)') + xlab('Time (Years)') 
          png(paste0(name,".scores.weighted.risk.dxage.",dxage,".all.png"),height = 900, width = 1100,res=300)
          print(plot1)
          dev.off()
          
          
        }
        
        
        
        
      }
      
      for (mmg.groups in c('< 6 months','0.5-1 year','1-2 years','2+ years','Never Had')) {
        if (mmg.groups == c('2+ years')) {
          pca.coords =merge(pred.df.targ.collapse.annotated, merged.df.all.tmp[merged.df.all.tmp$mammogram %in% c('2-3 years','3+ years') | merged.df.all.tmp$Cancer == 'Control',]) # merge(pred.df.targ.collapse.annotated, merged.df.all.tmp
          
        } else {
          pca.coords =merge(pred.df.targ.collapse.annotated,merged.df.all.tmp[merged.df.all.tmp$mammogram %in% (mmg.groups) | merged.df.all.tmp$Cancer == 'Control',]) # merge(pred.df.targ.collapse.annotated, merged.df.all.tmp
          
        }
        
        pca.coords$median.g = ifelse(pca.coords$methylation_score > score.cutoff,'Higher Test Median','Lower Test Median')
        pca.coords$Cancer.g= ifelse(pca.coords$reported == 'Control',0,1)
        female.weights= weightsf.females(pca.coords)
        female.weights.filt =female.weights
        cph = summary(coxph(Surv(censorship_time, Cancer.g) ~ median.g, data =  pca.coords,weights=female.weights))$logtest[3]
        hr = round(summary(coxph(Surv(censorship_time, Cancer.g) ~ median.g, data =  pca.coords,weights=female.weights))$coefficients[2],digits=3)
        if (length(unique(pca.coords$median.g)) > 1) {
          plot4 = ggsurvplot(
            fit = survfit(Surv(censorship_time/365, Cancer.g) ~ median.g, data =  pca.coords,weights = female.weights), 
            risk.table = F,
            #cumevents = TRUE,
            palette = c('#2a9d8f','#e9c46a'),
            legend.labs=c("Higher Test Median","Lower Test Median"),
            xlim = c(0,10), xlab = "Time (Years)", ylab = c("Survival probability"), break.x.by = 1,
            ylim = c(0.9,1),
            ggtheme = theme_bw(),
            
            font.tickslab = c(8,'bold'), 
            font.main = c(8, "bold",'black'),
            font.x = c(8, "bold",'black'),
            font.y = c(8, "bold",'black'),
            alpha = 0.6,
            font.legend = c(8), risk.table.fontsize = 5, risk.table.col = "black", size = 1) + 
            ggtitle(paste(hr,'.',formatC(cph, format = "e", digits = 2))) #+ theme(legend.position = 'none')
          
          png(paste0(name,".scores.surv_median.mmg.",mmg.groups,".cut.png"),height = 1100, width = 1100,res=300)
          print(plot4)
          dev.off()
          
          plot4 = ggsurvplot(
            fit = survfit(Surv(censorship_time/365, Cancer.g) ~ median.g, data =  pca.coords,weights = female.weights), 
            risk.table = F,
            #cumevents = TRUE,
            palette = c('#2a9d8f','#e9c46a'),
            legend.labs=c("Higher Test Median","Lower Test Median"),
            xlim = c(0,10), xlab = "Time (Years)", ylab = c("Survival probability"), break.x.by = 1,
            ylim = c(0,1),
            ggtheme = theme_bw(),
            
            font.tickslab = c(8,'bold'), 
            font.main = c(8, "bold",'black'),
            font.x = c(8, "bold",'black'),
            font.y = c(8, "bold",'black'),
            alpha = 0.6,
            font.legend = c(8), risk.table.fontsize = 8, risk.table.col = "black", size = 1) + ggtitle(paste(hr,'.',formatC(cph, format = "e", digits = 2)))   #+ theme(legend.position = 'none')
          
          png(paste0(name,".scores.surv_median.mmg.",mmg.groups,".all.png"),height = 1100, width = 1100,res=300)
          print(plot4)
          dev.off()
          
          
          pca.coords$quartile <- dplyr::ntile(pca.coords$methylation_score, 4) 
          quartile = pca.coords$quartile 
          cph = summary(coxph(Surv(censorship_time, Cancer.g) ~ quartile, data =  pca.coords,weights = female.weights.filt))$logtest[3]
          hr = round(summary(coxph(Surv(censorship_time, Cancer.g) ~ quartile, data =  pca.coords,weights=female.weights.filt))$coefficients[2],digits=3)
          survfit.obj = survfit(Surv(censorship_time/365, Cancer.g) ~ quartile, data =  pca.coords,weights = female.weights.filt)
          quart.col = c('#2a9d8f','#e9c46a','#EA526F','#586994')
          names(quart.col) = c("1","2","3","4")
          quart.col.plot = quart.col[names(quart.col) %in%  pca.coords.filt$quartile ]
          plot4 = ggsurvplot(
            fit = survfit.obj, 
            risk.table = F,
            #cumevents = TRUE,
            palette = quart.col.plot,
            legend.labs=names(quart.col.plot),
            xlim = c(0,10), xlab = "Time (Years)", ylab = c("Survival probability"), break.x.by = 1,
            ylim=c(0.9,1),
            ggtheme = theme_bw(),
            
            font.tickslab = c(8,'bold'), 
            font.main = c(8, "bold",'black'),
            font.x = c(8, "bold",'black'),
            font.y = c(8, "bold",'black'),  alpha = 0.6,
            font.legend = c(8), risk.table.fontsize = 5, risk.table.col = "black", size = 1)  + 
            ggtitle(paste0(hr,'.',formatC(cph, format = "e", digits = 2)) )  #+ the #+ theme(legend.position = 'none')
          
          
          png(paste0(name,".scores.surv_quart.mmg.",mmg.groups,".png"),height = 1100, width = 1100,res=300)
          print(plot4)
          dev.off()
          
          
          plot4 = ggsurvplot(
            fit = survfit.obj, 
            risk.table = F,
            #cumevents = TRUE,
            palette = quart.col.plot,
            legend.labs=names(quart.col.plot),
            xlim = c(0,10), xlab = "Time (Years)", ylab = c("Survival probability"), break.x.by = 1,
            ylim=c(0,1),
            ggtheme = theme_bw(),
            
            font.tickslab = c(8,'bold'), 
            font.main = c(8, "bold",'black'),
            font.x = c(8, "bold",'black'),
            font.y = c(8, "bold",'black'), alpha = 0.6,
            font.legend = c(8), risk.table.fontsize = 5, risk.table.col = "black", size = 1)  + 
            ggtitle(paste0(hr,'.',formatC(cph, format = "e", digits = 2)) )  #+ the #+ theme(legend.position = 'none')
          
          png(paste0(name,".scores.surv_quart.mmg.",mmg.groups,".full.png"),height = 1100, width = 1100,res=300)
          print(plot4)
          dev.off()
          
          pca.coords$quartile <- ifelse(pca.coords$methylation_score >=0 & pca.coords$methylation_score < 0.25, c(1),
                                        ifelse(pca.coords$methylation_score >=0.25 & pca.coords$methylation_score < 0.5, c(2),
                                               ifelse(pca.coords$methylation_score >=0.5 & pca.coords$methylation_score < 0.75, c(3),
                                                      ifelse(pca.coords$methylation_score >=0.75 & pca.coords$methylation_score < 1, c(4),'Other'))))
          quartile = pca.coords$quartile 
          cph = summary(coxph(Surv(censorship_time, Cancer.g) ~ quartile, data =  pca.coords,weights = female.weights.filt))$logtest[3]
          hr = round(summary(coxph(Surv(censorship_time, Cancer.g) ~ quartile, data =  pca.coords,weights=female.weights.filt))$coefficients[2],digits=3)
          survfit.obj = survfit(Surv(censorship_time/365, Cancer.g) ~ quartile, data =  pca.coords,weights = female.weights.filt)
          quart.col = c('#2a9d8f','#e9c46a','#EA526F','#586994')
          names(quart.col) = c("1","2","3","4")
          quart.col.plot = quart.col[names(quart.col) %in%  pca.coords$quartile ]
          plot4 = ggsurvplot(
            fit = survfit.obj, 
            risk.table = F,
            #cumevents = TRUE,
            palette = quart.col.plot,
            legend.labs=names(quart.col.plot),
            xlim = c(0,10), xlab = "Time (Years)", ylab = c("Survival probability"), break.x.by = 1,
            ylim=c(0.9,1),
            ggtheme = theme_bw(),
            
            font.tickslab = c(8,'bold'), 
            font.main = c(8, "bold",'black'),
            font.x = c(8, "bold",'black'),
            font.y = c(8, "bold",'black'),
            alpha = 0.6,
            font.legend = c(8), risk.table.fontsize = 5, risk.table.col = "black", size = 1)  + 
            ggtitle(paste0(hr,'.',formatC(cph, format = "e", digits = 2)) )  #+ the #+ theme(legend.position = 'none')
          
          
          png(paste0(name,".scores.surv_quart.fixed.mmg.",mmg.groups,".png"),height = 1100, width = 1100,res=300)
          print(plot4)
          dev.off()
          
          
          km=survfit(formula = Surv(censorship_time, Cancer.g) ~ median.g, data =  pca.coords,weights=female.weights.filt)
          tmp = summary(km,times=c(1:3650))
          km.summary.weighted.test = data.frame(time = tmp$time, surv.prob =tmp$surv, surv.prob.l =tmp$lower,surv.prob.h = tmp$upper,risk = gsub('risk=','',tmp$strata))
          colnames(km.summary.weighted.test)[2:4]= c('surv','lower','upper')
          #combined.validation.surv = rbind(km.summary.weighted.test,surv.prob.validation[,-3])
          plot1= ggplot(km.summary.weighted.test) +
            geom_line(aes(x = time/365, y = 100-surv*100, col = risk),size = 1,alpha= 0.7) + 
            
            scale_fill_manual(values = c('#B9F3EC','#FDF2C4'))+
            scale_color_manual(values = c('#2a9d8f','#e9c46a'))+
            scale_x_continuous(breaks = seq(0,max.time,1),limits = c(0,max.time))+
            theme_bw()+
            theme(text = element_text(size=8),
                  axis.text=element_text(size=8, face = "bold"),
                  axis.title=element_text(size=8,face="bold"),
                  legend.position = "none")+
            #scale_y_continuous(limits = c(0,1), breaks = seq(0,1,0.2)) +
            ylab('Breast Cancer Probability (%)') + xlab('Time (Years)') 
          png(paste0(name,".scores.weighted.risk.mmg.",mmg.groups,".all.png"),height = 900, width = 1100,res=300)
          print(plot1)
          dev.off()
        }
      }
      
      for (subtypes in c('HR+\nHER2-','HR-\nHER2-','HR+\nHER2+','HR-\nHER2+','NR')) {
        subtypes.name = gsub('\n','',subtypes)
        pca.coords =merge(pred.df.targ.collapse.annotated,merged.df.all.tmp[merged.df.all.tmp[,'HR.HER2'] %in% c(subtypes,'Control'),]) # merge(pred.df.targ.collapse.annotated, merged.df.all.tmp
        print(table(pca.coords$HR.HER2))
        pca.coords$median.g = ifelse(pca.coords$methylation_score > score.cutoff,'Higher Test Median','Lower Test Median')
        pca.coords$Cancer.g= ifelse(pca.coords$reported == 'Control',0,1)
        female.weights= weightsf.females(pca.coords)
        female.weights.filt =female.weights
        
        cph = summary(coxph(Surv(censorship_time, Cancer.g) ~ median.g, data =  pca.coords,weights=female.weights))$logtest[3]
        hr = round(summary(coxph(Surv(censorship_time, Cancer.g) ~ median.g, data =  pca.coords,weights=female.weights))$coefficients[2],digits=3)
        if (length(unique(pca.coords$median.g)) > 1) {
          plot4 = ggsurvplot(
            fit = survfit(Surv(censorship_time/365, Cancer.g) ~ median.g, data =  pca.coords,weights = female.weights.filt), 
            risk.table = F,
            #cumevents = TRUE,
            palette = c('#2a9d8f','#e9c46a'),
            legend.labs=c("Higher Test Median","Lower Test Median"),
            xlim = c(0,10), xlab = "Time (Years)", ylab = c("Survival probability"), break.x.by = 1,
            ylim = c(0.9,1),
            ggtheme = theme_bw(),
            
            font.tickslab = c(8,'bold'), 
            font.main = c(8, "bold",'black'),
            font.x = c(8, "bold",'black'),
            font.y = c(8, "bold",'black'),
            alpha = 0.6,
            font.legend = c(8), risk.table.fontsize = 5, risk.table.col = "black", size = 1) + 
            ggtitle(paste(hr,'.',formatC(cph, format = "e", digits = 2))) #+ theme(legend.position = 'none')
          
          png(paste0(name,".scores.surv_median.subtype.",subtypes.name,".cut.png"),height = 1100, width = 1100,res=300)
          print(plot4)
          dev.off()
          
          plot4 = ggsurvplot(
            fit = survfit(Surv(censorship_time/365, Cancer.g) ~ median.g, data =  pca.coords,weights = female.weights.filt), 
            risk.table = F,
            #cumevents = TRUE,
            palette = c('#2a9d8f','#e9c46a'),
            legend.labs=c("Higher Test Median","Lower Test Median"),
            xlim = c(0,10), xlab = "Time (Years)", ylab = c("Survival probability"), break.x.by = 1,
            ylim = c(0,1),
            ggtheme = theme_bw(),
            
            font.tickslab = c(8,'bold'), 
            font.main = c(8, "bold",'black'),
            font.x = c(8, "bold",'black'),
            font.y = c(8, "bold",'black'),
            alpha = 0.6,
            font.legend = c(8), risk.table.fontsize = 8, risk.table.col = "black", size = 1) + ggtitle(paste(hr,'.',formatC(cph, format = "e", digits = 2)))   #+ theme(legend.position = 'none')
          
          png(paste0(name,".scores.surv_median.subtype.",subtypes.name,".all.png"),height = 1100, width = 1100,res=300)
          print(plot4)
          dev.off()
          
          pca.coords.filt =pca.coords
          pca.coords.filt$quartile <- dplyr::ntile(pca.coords.filt$methylation_score, 4) 
          quartile = pca.coords.filt$quartile 
          cph = summary(coxph(Surv(censorship_time, Cancer.g) ~ quartile, data =  pca.coords.filt,weights = female.weights.filt))$logtest[3]
          hr = round(summary(coxph(Surv(censorship_time, Cancer.g) ~ quartile, data =  pca.coords.filt,weights=female.weights.filt))$coefficients[2],digits=3)
          survfit.obj = survfit(Surv(censorship_time/365, Cancer.g) ~ quartile, data =  pca.coords.filt,weights = female.weights.filt)
          quart.col = c('#2a9d8f','#e9c46a','#EA526F','#586994')
          names(quart.col) = c("1","2","3","4")
          quart.col.plot = quart.col[names(quart.col) %in%  pca.coords.filt$quartile ]
          plot4 = ggsurvplot(
            fit = survfit.obj, 
            risk.table = F,
            #cumevents = TRUE,
            palette = quart.col.plot,
            legend.labs=names(quart.col.plot),
            xlim = c(0,10), xlab = "Time (Years)", ylab = c("Survival probability"), break.x.by = 1,
            ylim=c(0.9,1),
            ggtheme = theme_bw(),
            
            font.tickslab = c(8,'bold'), 
            font.main = c(8, "bold",'black'),
            font.x = c(8, "bold",'black'),
            font.y = c(8, "bold",'black'),  alpha = 0.6,
            font.legend = c(8), risk.table.fontsize = 5, risk.table.col = "black", size = 1)  + 
            ggtitle(paste0(hr,'.',formatC(cph, format = "e", digits = 2)) )  #+ the #+ theme(legend.position = 'none')
          
          
          png(paste0(name,".scores.surv_quart.subtype.",subtypes.name,".png"),height = 1100, width = 1100,res=300)
          print(plot4)
          dev.off()
          
          
          plot4 = ggsurvplot(
            fit = survfit.obj, 
            risk.table = F,
            #cumevents = TRUE,
            palette = quart.col.plot,
            legend.labs=names(quart.col.plot),
            xlim = c(0,10), xlab = "Time (Years)", ylab = c("Survival probability"), break.x.by = 1,
            ylim=c(0,1),
            ggtheme = theme_bw(),
            
            font.tickslab = c(8,'bold'), 
            font.main = c(8, "bold",'black'),
            font.x = c(8, "bold",'black'),
            font.y = c(8, "bold",'black'), alpha = 0.6,
            font.legend = c(8), risk.table.fontsize = 5, risk.table.col = "black", size = 1)  + 
            ggtitle(paste0(hr,'.',formatC(cph, format = "e", digits = 2)) )  #+ the #+ theme(legend.position = 'none')
          
          png(paste0(name,".scores.surv_quart.subtype.",subtypes.name,".full.png"),height = 1100, width = 1100,res=300)
          print(plot4)
          dev.off()
          
          pca.coords$quartile <- ifelse(pca.coords$methylation_score >=0 & pca.coords$methylation_score < 0.25, c(1),
                                        ifelse(pca.coords$methylation_score >=0.25 & pca.coords$methylation_score < 0.5, c(2),
                                               ifelse(pca.coords$methylation_score >=0.5 & pca.coords$methylation_score < 0.75, c(3),
                                                      ifelse(pca.coords$methylation_score >=0.75 & pca.coords$methylation_score < 1, c(4),'Other'))))
          quartile = pca.coords.filt$quartile 
          cph = summary(coxph(Surv(censorship_time, Cancer.g) ~ quartile, data =  pca.coords.filt,weights = female.weights.filt))$logtest[3]
          hr = round(summary(coxph(Surv(censorship_time, Cancer.g) ~ quartile, data =  pca.coords.filt,weights=female.weights.filt))$coefficients[2],digits=3)
          survfit.obj = survfit(Surv(censorship_time/365, Cancer.g) ~ quartile, data =  pca.coords.filt,weights = female.weights.filt)
          quart.col = c('#2a9d8f','#e9c46a','#EA526F','#586994')
          names(quart.col) = c("1","2","3","4")
          quart.col.plot = quart.col[names(quart.col) %in%  pca.coords.filt$quartile ]
          plot4 = ggsurvplot(
            fit = survfit.obj, 
            risk.table = F,
            #cumevents = TRUE,
            palette = quart.col.plot,
            legend.labs=names(quart.col.plot),
            xlim = c(0,10), xlab = "Time (Years)", ylab = c("Survival probability"), break.x.by = 1,
            ylim=c(0.9,1),
            ggtheme = theme_bw(),
            
            font.tickslab = c(8,'bold'), 
            font.main = c(8, "bold",'black'),
            font.x = c(8, "bold",'black'),
            font.y = c(8, "bold",'black'),
            alpha = 0.6,
            font.legend = c(8), risk.table.fontsize = 5, risk.table.col = "black", size = 1)  + 
            ggtitle(paste0(hr,'.',formatC(cph, format = "e", digits = 2)) )  #+ the #+ theme(legend.position = 'none')
          
          
          png(paste0(name,".scores.surv_quart.fixed.subtype.",subtypes.name,".png"),height = 1100, width = 1100,res=300)
          print(plot4)
          dev.off()
          
          
          km=survfit(formula = Surv(censorship_time, Cancer.g) ~ median.g, data =  pca.coords.filt,weights=female.weights.filt)
          tmp = summary(km,times=c(1:3650))
          km.summary.weighted.test = data.frame(time = tmp$time, surv.prob =tmp$surv, surv.prob.l =tmp$lower,surv.prob.h = tmp$upper,risk = gsub('risk=','',tmp$strata))
          colnames(km.summary.weighted.test)[2:4]= c('surv','lower','upper')
          #combined.validation.surv = rbind(km.summary.weighted.test,surv.prob.validation[,-3])
          plot1= ggplot(km.summary.weighted.test) +
            geom_line(aes(x = time/365, y = 100-surv*100, col = risk),size = 1,alpha= 0.7) + 
            
            scale_fill_manual(values = c('#B9F3EC','#FDF2C4'))+
            scale_color_manual(values = c('#2a9d8f','#e9c46a'))+
            scale_x_continuous(breaks = seq(0,max.time,1),limits = c(0,max.time))+
            theme_bw()+
            theme(text = element_text(size=8),
                  axis.text=element_text(size=8, face = "bold"),
                  axis.title=element_text(size=8,face="bold"),
                  legend.position = "none")+
            #scale_y_continuous(limits = c(0,1), breaks = seq(0,1,0.2)) +
            ylab('Breast Cancer Probability (%)') + xlab('Time (Years)') 
          png(paste0(name,".scores.weighted.risk.subtype.",subtypes.name,".all.png"),height = 900, width = 1100,res=300)
          print(plot1)
          dev.off()
        }
        
      }
    }
    
    #subtyping
    
    
  }
  
  
  
  return(combined.auroc)
  
}

combined.auroc.breast = combined.auroc
#savedir='/.mounts/labs/awadallalab/private/ncheng/cfmedip_data/cptp_samples/fragmentation/aix13.updated1/all.sample.fragmentation/breast.cancer.cv/regulatory.regions.v10.newsplit7//cv.figures.genhancer.base.discovery.final.90.v2/'
saveRDS(combined.auroc,paste0(figdir,'discovery.brca.auroc.RDS'))
saveRDS(pred.df.targ,paste0(figdir,'discovery.brca.scores.RDS'))

#savedir='/.mounts/labs/awadallalab/private/ncheng/cfmedip_data/cptp_samples/fragmentation/aix13.updated1/all.sample.fragmentation/breast.cancer.cv/regulatory.regions.v10/validation.breast.test/'
saveRDS(combined.auroc,paste0(figdir,'validation.brca.auroc.RDS'))
saveRDS(pred.df.targ,paste0(figdir,'validation.brca.scores.RDS'))
integrated.calc = function(combined.predictions.df.male.filt) {
  
  
  auroc.summary=split(combined.predictions.df.male.filt,combined.predictions.df.male.filt$feature.type1)
  auroc.summary=sapply(auroc.summary, function(x) auc_calc(x,labels= c('Control','Cancer')))
  tmp = combined.predictions.df.male.filt#[combined.predictions.df.male.filt$feature.type1 %in% names(auroc.summary[auroc.summary>0.55]),c('GRP_Id','methylation_score','reported')]
  
  tmp = ddply(tmp, c('GRP_Id','reported'),numcolwise(mean))
  tmp$feature.type1 = 'Integrated'
  return.df = rbind( combined.predictions.df.male.filt,tmp)
  
  return(return.df)
}

train.opt =readRDS('/.mounts/labs/awadallalab/private/ncheng/cfmedip_data/cptp_samples/fragmentation/aix13.updated1/methylation.zero.cv10/cv.top.median.methscores.RDS')

####ploting auroc summary####
variables = c('ci','auroc')
#combined.auroc = readRDS(paste0(savedir,'discovery.brca.auroc.RDS'))

stage_colors.score = c('I' = "#3B0F70FF",'II' = "#8C2981FF",'III/IV' = "#DE4968FF",'NR' = '#66462C')
morph_colors.score = c('Ductal+Lobular' = "#AFA2FF",'Ductal' = "#EF8275",'Lobular' = "#8A4F7D")
morph_colors.score = c('Duct\n+Lob' = "#AFA2FF",'Duct' = "#EF8275",'Lob' = "#8A4F7D")

grade_colors.sample = c('1' = '#82A6B1','2' = '#35605A', '3' = '#2A324B','NR'='#66462C')

age_colors = c('35-45' = '#C5C392','45-55'= '#FFB20F','55-65'='#FF4B3E', '65-75'='#972D07')
dxage_colors = c('30-50' = '#C5C392','50-60'= '#FFB20F','60-70'='#FF4B3E', '70-80'='#972D07')


mmg_colors = c('Never' = '#93A3B1',
               '< 0.5'='#78BC61',
               '0.5-1'='#C0C781',
               '1-2'='#C59B76',
               '2+'='#AD343E')
subtype_colors =  c('HR+\nHER2+' = '#8FD5A6','HR+\nHER2-' = '#BF1363',"HR-\nHER2+" = '#E09F3E','HR-\nHER2-'='#545E75','NR' ='#66462C')
combined.colors = c(subtype_colors,stage_colors.score,morph_colors.score,grade_colors.sample,dxage_colors,mmg_colors,age_colors)
#combined.colors = unique(combined.colors)
#boxplot of score
targ.vars = c('Stage','DxAge','Subtype','Grade','Morphology','Last Mammogram')
combined.auroc$var.names=ifelse(combined.auroc$var.group %in% 'DxAge','Age at Diagnosis',combined.auroc$var.group)
combined.auroc$var.names=ifelse(combined.auroc$var.group %in% 'Last Mammogram','Last Mammogram\nBefore Baseline (Years)',combined.auroc$var.names)

combined.auroc = combined.auroc[!duplicated(combined.auroc$var),]
combined.auroc$var = gsub('\\+Lob','\n+Lob',gsub('al','',gsub('ular','',combined.auroc$var)))

plot1 = ggplot(combined.auroc[combined.auroc$var.group %in% targ.vars,],aes(x = var, y = auc, fill = var)) +
  geom_bar(stat='identity')+
  geom_errorbar(aes(ymin=auc.lower, ymax=auc.upper), width=.2,
                position=position_dodge(.9)) +
  scale_fill_manual(values = c(combined.colors))+ #+ ggtitle(title) +
  theme_bw()+ 
  scale_y_continuous(limits=c(0,1))+
  
  theme(text = element_text(size=8),
        axis.text=element_text(size=8, face = "bold"),
        axis.title=element_text(size=8,face="bold"),
        legend.position = "none",
        strip.text=element_text(size=8, face = "bold"),
        strip.background = element_rect(fill = "white"))+ 
  facet_grid2(. ~ var.names,scales = 'free_x') +
  guides(colour = guide_legend(override.aes = list(size=4),nrow=4)) +
  xlab('') + ylab('AUROC (95% CI)')

png(paste0(name,'all.combined.auroc.png'),height = 600,width=520*length(targ.vars),res = 300)
print(plot1)
dev.off()


plot1 = ggplot(combined.auroc[combined.auroc$var.group %in% targ.vars,],aes(x = var, y = ci, fill = var)) +
  geom_bar(stat='identity')+
  geom_errorbar(aes(ymin=ci.lower, ymax=ci.upper), width=.2,
                position=position_dodge(.9)) +
  scale_fill_manual(values = c(combined.colors))+ #+ ggtitle(title) +
  theme_bw()+ 
  scale_y_continuous(limits=c(0,1))+
  
  theme(text = element_text(size=8),
        axis.text=element_text(size=8, face = "bold"),
        axis.title=element_text(size=8,face="bold"),
        legend.position = "none",
        strip.text=element_text(size=8, face = "bold"),
        strip.background = element_rect(fill = "white"))+ 
  facet_grid2(. ~ var.names,scales = 'free_x') +
  guides(colour = guide_legend(override.aes = list(size=4),nrow=4)) +
  xlab('') + ylab('Concordance Index (95% CI)')

png(paste0(name,'all.combined.ci.png'),height = 600,width=520*length(targ.vars),res = 300)
print(plot1)
dev.off()



####plotting auc for brca#####
#top.perf.df = readRDS('/.mounts/labs/awadallalab/private/ncheng/cfmedip_data/cptp_samples/fragmentation/aix13.updated1/methylation.insert.select1/predictions1//cv.top.median.methscores.RDS')
top.perf.df = targ.perf.all[targ.perf.all$feature == 'TB' & 
                              targ.perf.all$test =='deseq.log2.std' &
                              targ.perf.all$direction == 'bm1.hyper' & 
                              targ.perf.all$model == 'RF' & 
                              targ.perf.all$n.features == 25,]
top.perf.df$GRP_Id= gsub('_combined.*','',top.perf.df$GRP_Id)
pred.df.targ = top.perf.df
#pred.df.targ=readRDS(paste0(savedir,'top.perf.discovery.cv.RDS'))
pred.df.targ.collapse = ddply(pred.df.targ[,c('GRP_Id','reported','methylation_score')],c('GRP_Id','reported'),numcolwise(mean))
last.mmg = read.csv('/.mounts/labs/awadallalab/private/ncheng/cfmedip_data/cptp_samples/participant_data/ohs_baseline_mmg_vars.csv',header=T,stringsAsFactors = F)
discovery.set = readRDS('/.mounts/labs/awadallalab/private/ncheng/manuscripts/early.cancer.risk/discovery.set3.samples.RDS')
brca.morphology = c('85003' = 'Infiltrating ductal carcinoma',
                    '85203' = 'Lobular NOS',
                    '85223'='Infiltrating ductal and lobular carcinoma',
                    '85233' = 'Infiltrating ductal mixed with other',
                    '85433' = 'Paget disease + infiltrating ductal carcinoma',
                    '85073' = 'Micropapillary carcinoma',
                    '84803' = 'Mucinous adenocarcinoma',
                    '82603' = 'Papillary carcinoma NOS',
                    '82113' = 'Tubular adenocarcinoma',
                    '80103' = 'Carcinoma NOS'
)
sample.info.filt.pretime = discovery.set
targ.auc =pred.df.targ.collapse

my_comparisons = list(c('Control\n(MFiller)','Control\n(UFiller)'),
                      c('Breast\n(MFiller)','Breast\n(UFiller)'),
                      c('Breast\n(MFiller)','Control\n(MFiller)'),
                      c('Breast\n(UFiller)','Control\n(UFiller)'))
cancer.colors = c('Control' = '#BFCBC2','Pancreas' ='#09814A','Breast'='#f564a9','Prostate'= '#0094C6')
mean.perf.df.targ.tmp.annotated = merge(targ.auc,sample.info.filt.pretime,by='GRP_Id')
mean.perf.df.targ.tmp.annotated$filler.cancer = ifelse(mean.perf.df.targ.tmp.annotated$filler == 'MFiller',
                                                       paste0(mean.perf.df.targ.tmp.annotated$Cancer,'\n(MFiller)'),
                                                       paste0(mean.perf.df.targ.tmp.annotated$Cancer,'\n(UFiller)'))
mean.perf.df.targ.tmp.annotated = mean.perf.df.targ.tmp.annotated[mean.perf.df.targ.tmp.annotated$Cancer %in% c('Control','Breast') & mean.perf.df.targ.tmp.annotated$Sex == 'Female',]

plot1 = ggplot(mean.perf.df.targ.tmp.annotated,aes(x = filler.cancer, y = methylation_score, col = Cancer)) +
  geom_boxplot() + geom_jitter(width = 0.1)+
  scale_color_manual(values = c(cancer.colors))+
  scale_y_continuous(limits = c(0,1.5),breaks = seq(0,1,0.25))+#+ ggtitle(title) +
  theme_bw()+ 
  #facet_grid2(.~ filler)+
  theme(text = element_text(size=8),
        axis.text=element_text(size=8),
        axis.title=element_text(size=8,face="bold"),
        legend.position = "none")+ guides(colour = guide_legend(override.aes = list(size=4),nrow=4)) +
  xlab('Group') + ylab('Classification Score') +  
  stat_compare_means(comparisons = my_comparisons,label.y=c(1,1.1,1.2,1.3),size=3)

#,method='wilcox.test')
#+ # Add pairwise comparisons p-valu  
png(paste0(figdir,'methscore.filler.brca.png'),height = 1100,width=1000,res = 300)
print(plot1)
dev.off()



merged.df.all.filt = sample.info.filt.pretime#merged.df.all[merged.df.all$ResearchId %in% last.mmg.6m$ResearchId| merged.df.all$Cancer =='Control',]
brca.auroc = brca.auc.plot(pred.df.targ = pred.df.targ.collapse[pred.df.targ.collapse$GRP_Id %in% discovery.set$GRP_Id,],
                           name=paste0(figdir,'methylation.auroc'),
                           merged.df.all = discovery.set,
                           figdir,
                           dx.all = F,
                           score.cutoff=median(pred.df.targ.collapse$methylation_score))
pred.df.targ = pred.df.targ.collapse[pred.df.targ.collapse$GRP_Id %in% merged.df.all.filt$GRP_Id,];
name=paste0(figdir,'methylation.auroc');
merged.df.all = merged.df.all.filt;
figdir;
dx.all = F;
score.cutoff=median(pred.df.targ.collapse$methylation_score)





####computing HR breast#####
ohs.qx.pulled = readRDS('/.mounts/labs/awadallalab/private/ncheng/cfmedip_data/cptp_samples/participant_data/complete.qx.pulled.samples.RDS')
ohs.qx.pulled.multimatch = readRDS('/.mounts/labs/awadallalab/private/ncheng/cfmedip_data/cptp_samples/participant_data/complete.qx.multimatched.samples.RDS')
#
discovery.set = readRDS('/.mounts/labs/awadallalab/private/ncheng/manuscripts/early.cancer.risk/discovery.set3.samples.RDS')
validation.set = readRDS('/.mounts/labs/awadallalab/private/ncheng/manuscripts/early.cancer.risk/validation.set3.samples.RDS')
sample.info =  rbind(discovery.set,validation.set)
combined.info = rbind(discovery.set,validation.set)
combined.info.filt = combined.info[,c('GRP_Id','SDC_GENDER','SDC_AGE_CALC','ALC_CUR_FREQ','SMK_CIG_STATUS','censorship_time','Cancer','DIS_CANCER_M_BREAST','DIS_CANCER_SIB_BREAST','DIS_CANCER_CHILD_BREAST','DIS_CANCER_F_EVER','DIS_CANCER_M_EVER','DIS_CANCER_SIB_EVER','DIS_CANCER_CHILD_EVER')]
combined.info.filt$Cohort = 'OHS'
combined.info.filt[is.na(combined.info.filt$ALC_CUR_FREQ),'ALC_CUR_FREQ'] = '-7'
combined.info.filt = merge(combined.info.filt,ohs.qx.pulled[,c('GRP_Id',colnames(ohs.qx.pulled)[grepl('.*BMI.*',colnames(ohs.qx.pulled))])] )
combined.info.filt= unique(combined.info.filt)
#cohort wide
canpath.master= read.csv('/.mounts/labs/awadallalab/public/EpiCan/CanPath_MasterFile_May2020_With_BMI_Cardio_PRS_IDsmatchs.txt',header = T, stringsAsFactors = F)
canpath.master = canpath.master[canpath.master$SDC_GENDER == 2,]
canpath.master.control = canpath.master[grepl('OHS.*',canpath.master$ADM_STUDY_DATASET) == T,c('EpiCan_ID','ALC_CUR_FREQ','SDC_GENDER','ADM_QX_COMPLETION','SDC_AGE_CALC','SMK_CIG_STATUS','DIS_CANCER_M_BREAST','DIS_CANCER_SIB_BREAST','DIS_CANCER_CHILD_BREAST','BMI','DIS_CANCER_F_EVER','DIS_CANCER_M_EVER','DIS_CANCER_SIB_EVER','DIS_CANCER_CHILD_EVER')]
incident.cases = read.table("/.mounts/labs/awadallalab/public/EpiCan/OHS.incident.cancers.txt", sep ='\t',quote = '"', header = T,stringsAsFactors = F)
incident.cases = incident.cases[grepl('.*Breast*',incident.cases$CURR_TOPOG_DESC),]
epican.ids = read.csv("/.mounts/labs/awadallalab/public/EpiCan/EpiCAN_genotype_crosswalk_nondup.csv",quote = '"', header = T,stringsAsFactors = F)
epican.ids.incident = epican.ids[epican.ids[,1] %in% incident.cases$participantid,]

canpath.master.control$Cancer = ifelse(canpath.master.control$EpiCan_ID %in% epican.ids.incident$EpiCan_ID,'Cancer','Control')

canpath.master.control$censorship_time = as.Date('2020-01-01',format = '%Y-%m-%d') - as.Date(canpath.master.control$ADM_QX_COMPLETION,format = '%Y-%m-%d') 
canpath.master.control$GRP_Id =canpath.master.control$EpiCan_ID
canpath.master.control = canpath.master.control[,c('GRP_Id','SDC_GENDER','ALC_CUR_FREQ','SDC_AGE_CALC','SMK_CIG_STATUS','censorship_time','Cancer','DIS_CANCER_M_BREAST','DIS_CANCER_SIB_BREAST','DIS_CANCER_CHILD_BREAST','BMI','DIS_CANCER_F_EVER','DIS_CANCER_M_EVER','DIS_CANCER_SIB_EVER','DIS_CANCER_CHILD_EVER')]
canpath.master.control$Cohort = 'EpiCan'
canpath.master.control$Event = ifelse(canpath.master.control$Cancer == 'Control',0,1)
#ohs.master = ohs.master[ohs.master$]
#

#matching epican to canpath#
combined.info = sample.info
#colnames(combined.info) =  gsub('MD_','',colnames(combined.info) )
overlapping.vars = colnames(combined.info)
overlapping.vars = overlapping.vars[overlapping.vars %in% colnames(ohs.master)]
#a = merge(ohs.master, sample.info, by= overlapping.vars ) 
#exc = sample.info[!sample.info$GRP_Id %in% a$GRP_Id,]
qx.pull.match=F
if (qx.pull.match == T){
  ohs.qx.pulled = NULL
  ohs.qx.pulled.multimatch = NULL
  for (j in 1:nrow(combined.info)) {
    print(j/nrow(combined.info))
    targ.sample = sample.info[j,]
    targ.sample.filt =targ.sample[,overlapping.vars]
    non.na.vars = colnames(targ.sample.filt)[which(is.na(targ.sample.filt) == F)]
    na.vars =  colnames(targ.sample.filt)[which(is.na(targ.sample.filt) == T)]
    var.index = 1
    
    var.tmp = non.na.vars[1]
    exc.vars= c('DIN')
    ohs.master.filtered = ohs.master[ohs.master[,var.tmp] %in% targ.sample[,var.tmp],]
    while (nrow(ohs.master.filtered) > 1  & var.index < length(non.na.vars) ) {
      var.index = var.index+1
      
      var.tmp = non.na.vars[var.index]
      if (nrow(ohs.master.filtered[ohs.master.filtered[,var.tmp] %in% targ.sample[,var.tmp],]) > 0 ){
        ohs.master.filtered = ohs.master.filtered[ohs.master.filtered[,var.tmp] %in% targ.sample[,var.tmp],]
        
      }
      
    }
    if (nrow(ohs.master.filtered) ==1 ) {
      ohs.master.filtered$GRP_Id = targ.sample$GRP_Id
      ohs.qx.pulled = rbind(ohs.qx.pulled,ohs.master.filtered)
      
    } else {
      for (i in na.vars) {
        ohs.master.filtered[,i] = ifelse(is.na(ohs.master.filtered[,i]) == T | ohs.master.filtered[,i] == '',NA,ohs.master.filtered[,i])
        
        if (nrow(ohs.master.filtered[is.na(ohs.master.filtered[,i]) == T,])  > 0) {
          # ohs.master.filtered1=ohs.master.filtered
          
          ohs.master.filtered = ohs.master.filtered[is.na(ohs.master.filtered[,i]) == T,]
          
        }
      }
      if (nrow(ohs.master.filtered) == 1) {
        ohs.master.filtered$GRP_Id = targ.sample$GRP_Id
        ohs.qx.pulled = rbind(ohs.qx.pulled,ohs.master.filtered)
        
      } else {
        multiple.matches.remain = ohs.master.filtered
        multiple.matches.remain$GRP_Id =  targ.sample$GRP_Id
        ohs.qx.pulled.multimatch = rbind(ohs.qx.pulled.multimatch,multiple.matches.remain)
      }
      
    }
    
  }
  #setwd('/.mounts/labs/awadallalab/private/ncheng/cfmedip_data/cptp_samples/participant_data')
  #saveRDS(ohs.qx.pulled,'complete.qx.pulled.samples.RDS')
  #saveRDS(ohs.qx.pulled.multimatch,'complete.qx.multimatched.samples.RDS')
  #
  
}
combined.info.filt$Event = ifelse(combined.info.filt$Cancer == 'Control',0,1)
combined.info.all =rbind(combined.info.filt[,colnames(combined.info.filt) %in% colnames(canpath.master.control)],canpath.master.control)

combined.info.all$Alcohol.Frequency = ifelse(combined.info.all$ALC_CUR_FREQ == '-7' | is.na(combined.info.all$ALC_CUR_FREQ) == T,'Never',
                                             ifelse(combined.info.all$ALC_CUR_FREQ == '0','Former',
                                                    ifelse(combined.info.all$ALC_CUR_FREQ == '1','< 1 per month',
                                                           ifelse(combined.info.all$ALC_CUR_FREQ == '2','~1 per month',
                                                                  ifelse(combined.info.all$ALC_CUR_FREQ == '3','2-3 per month',
                                                                         ifelse(combined.info.all$ALC_CUR_FREQ == '4','1 per week',
                                                                                ifelse(combined.info.all$ALC_CUR_FREQ == '5','2-3 per week',
                                                                                       ifelse(combined.info.all$ALC_CUR_FREQ == '6','4-5 per week',
                                                                                              ifelse(combined.info.all$ALC_CUR_FREQ == '7','6-7 per week','Not Reported')))))))))
combined.info.all$Alcohol.Frequency = factor(combined.info.all$Alcohol.Frequency, levels = c('Never','Former','< 1 per month','~1 per month','2-3 per month','1 per week','2-3 per week','4-5 per week','6-7 per week'))

combined.info.all$Alch.con.group = ifelse(combined.info.all$Alcohol.Frequency %in% c('Never'),'Never',
                                          ifelse(combined.info.all$Alcohol.Frequency %in% c('Former'),'Former',
                                                 ifelse(combined.info.all$Alcohol.Frequency %in% c('< 1 per month','~1 per month','2-3 per month'),'Infrequent',
                                                        ifelse(combined.info.all$Alcohol.Frequency %in% c('1 per week','2-3 per week'),'Moderate',
                                                               ifelse(combined.info.all$Alcohol.Frequency %in% c('4-5 per week','6-7 per week'),'Frequent','Other')))))
combined.info.all$Alch.con.group  = factor(combined.info.all$Alch.con.group , levels = c('Never','Former','Infrequent','Moderate','Frequent'))

combined.info.all$Smoking.Frequency = ifelse(combined.info.all$SMK_CIG_STATUS == 0, 'Never',
                                             ifelse(combined.info.all$SMK_CIG_STATUS == 1, 'Former', 
                                                    ifelse(combined.info.all$SMK_CIG_STATUS == 2, 'Occasional',
                                                           ifelse(combined.info.all$SMK_CIG_STATUS == 3, 'Daily','Other'))))
combined.info.all$Smoking.Frequency = factor(as.character(combined.info.all$Smoking.Frequency), levels=c('Never','Former','Occasional','Daily'))
combined.info.all$Sex = ifelse(combined.info.all$SDC_GENDER == 1,'Male','Female')
combined.info.all$Event= ifelse(combined.info.all$Cancer == 'Control',0,1)
combined.info.all$censorship_time = abs(combined.info.all$censorship_time)
combined.info.all = combined.info.all[!is.na(combined.info.all$Alch.con.group),]
combined.info.all = combined.info.all[combined.info.all$censorship_time < 365*12,]
combined.info.all$Family.history.breast =ifelse(combined.info.all$DIS_CANCER_M_BREAST == 1 | 
                                                  combined.info.all$DIS_CANCER_SIB_BREAST == 1 |
                                                  combined.info.all$DIS_CANCER_CHILD_BREAST == 1, 'Family History','No Reported Family History')
combined.info.all$Family.history.all =ifelse(combined.info.all$DIS_CANCER_M_EVER == 1 | 
                                               combined.info.all$DIS_CANCER_F_EVER == 1 | 
                                               combined.info.all$DIS_CANCER_SIB_EVER == 1 |
                                               combined.info.all$DIS_CANCER_CHILD_EVER == 1, 'Family History','No Reported Family History')
combined.info.all[is.na(combined.info.all$Family.history.breast),'Family.history.breast'] = 'No Reported Family History'
combined.info.all$Family.history.breast  = factor(as.character(combined.info.all$Family.history.breast),levels = c('No Reported Family History','Family History'))
combined.info.all$Family.history.all  = factor(as.character(combined.info.all$Family.history.all),levels = c('No Reported Family History','Family History'))

#combined.info.all = merge(combined.info.all,ohs.qx.pulled[,c('GRP_Id',colnames(ohs.qx.pulled)[grepl('.*BMI.*',colnames(ohs.qx.pulled))])] )
combined.info.all$age_group = ifelse(combined.info.all$SDC_AGE_CALC >= 30 & combined.info.all$SDC_AGE_CALC < 40, '30-40',
                                     ifelse(combined.info.all$SDC_AGE_CALC >= 40 & combined.info.all$SDC_AGE_CALC < 50, '40-50',
                                            ifelse(combined.info.all$SDC_AGE_CALC >= 50 & combined.info.all$SDC_AGE_CALC < 60, '50-60',
                                                   ifelse(combined.info.all$SDC_AGE_CALC >= 60 & combined.info.all$SDC_AGE_CALC < 70, '60-70',
                                                          ifelse(combined.info.all$SDC_AGE_CALC >= 70 & combined.info.all$SDC_AGE_CALC < 80, '70-80','80+')))))
combined.info.all$age_group =factor(as.character(combined.info.all$age_group), levels = c('30-40','40-50','50-60','60-70','70-80','80+'))
combined.info.all$bmi.group.short = ifelse(combined.info.all$BMI <= 18.5, 'Underweight',
                                           ifelse(combined.info.all$BMI > 18.5 & combined.info.all$BMI <= 25,'Normal Weight',
                                                  ifelse(combined.info.all$BMI > 25 & combined.info.all$BMI <= 30,'Overweight',
                                                         ifelse(combined.info.all$BMI > 30,'Obese','Other'))))
combined.info.all$bmi.group.long = ifelse(combined.info.all$BMI <= 18.5, 'Underweight',
                                          ifelse(combined.info.all$BMI > 18.5 & combined.info.all$BMI <= 25,'Normal Weight',
                                                 ifelse(combined.info.all$BMI > 25 & combined.info.all$BMI <= 30,'Overweight',
                                                        ifelse(combined.info.all$BMI > 30 & combined.info.all$BMI <= 35,'Class 1 Obesity',
                                                               ifelse(combined.info.all$BMI > 35 & combined.info.all$BMI <= 40, 'Class 2 Obesity', 'Class 3 Obesity')))))
combined.info.all$bmi.group.short = factor(combined.info.all$bmi.group.short, levels = c('Normal Weight','Underweight','Overweight','Obese'))
combined.info.all$bmi.group.long = factor(combined.info.all$bmi.group.long, levels =c('Normal Weight','Underweight','Overweight','Class 1 Obesity','Class 2 Obesity','Class 3 Obesity'))
#saveRDS(combined.info.all,paste0('/.mounts/labs/awadallalab/private/ncheng/cfmedip_data/cptp_samples/participant_data/combined.ohs.full.qx.RDS'))
#
library(survival)
library(ggfortify)
library(survminer)
#survival plots on cancer dx with grouped alcohol
km_trt_fit <- survfit(Surv(censorship_time/365, Event) ~ Alch.con.group , data=combined.info.all)
km_trt_fit.male <- survfit(Surv(censorship_time/365, Event) ~ Alch.con.group , data=combined.info.all[combined.info.all$Sex == 'Male',])
km_trt_fit.female <- survfit(Surv(censorship_time/365, Event) ~ Alch.con.group , data=combined.info.all[combined.info.all$Sex == 'Female',])

sex.goup.col = c("Female Cancer" ='#EF476F',"Female Control" = '#FAC7D3', "Male Cancer"  ='#26547C', "Male Control"='#B1CEE7')
alc.group.col =c( Never = 'Grey',Former="#63A375",Infrequent="#EDC79B",Moderate="#D57A66",Frequent="#713E5A")
diagnosis_time_colors = c("white","#440154FF","#414487FF","#2A788EFF","#22A884FF","#7AD151FF","#FDE725FF")
cancer.col=c(Breast='#A40E4C',Control='#AF8D86', Prostate='#4C86A8', Pancreatic='#ACC3A6')
sex.goup.col = c( "Female Cancer" ='#EF476F',"Female Control" = '#FAC7D3',    "Male Cancer"  ='#26547C', "Male Control"='#B1CEE7')
sex.col = c(Female = '#FAC7D3', Male='#B1CEE7')


plot1 = ggsurvplot(fit=km_trt_fit,conf.int = F,ylim = c(0.9,1),censor=F,legend = 'right', 
                   palette = unname(alc.group.col),
                   ggtheme = theme_bw(),
                   legend.title = "Alcohol Frequency",
                   legend.labs = c("Never (n=29,415)", "Former (n=17,703)","1 or less drinks per month (n=57,520)","1-3 drinks per week (n=88,691)","4+ drinks per week (n=62,381)"),
                   xlab='Time to Cancer Event (Years)',
                   ylab='Cancer-Free Fraction'
) 

plot2 = ggsurvplot(fit=km_trt_fit.male, conf.int = F,ylim = c(0.9,1),censor=F,legend = 'bottom',
                   palette = unname(alc.group.col),
                   ggtheme = theme_bw(),
                   legend.title = "Alcohol Frequency",
                   legend.labs = c("Never\n(n=11,021)", "Former\n(n=6,685)","Infrequent\n(n=31,682)","Moderate\n(n=37,131)","Frequent\n(n=30,769)"),
                   xlab='Time to Cancer Event',
                   ylab='Cancer-Free Fraction') 

plot3 = ggsurvplot(fit=km_trt_fit.female, conf.int = F,ylim = c(0.9,1),censor=F,legend = 'bottom',
                   palette = unname(alc.group.col),
                   ggtheme = theme_bw(),
                   legend.title = "Alcohol Frequency",
                   legend.labs = c("Never\n(n=18,397)", "Former\n(n=11,040)","Infrequent\n(n=75,600)","Moderate\n(n=51,621)","Frequent\n(n=31,668)"),
                   xlab='Time to Cancer Event',
                   ylab='Cancer-Free Fraction')  

plot.list = list(plot1,plot2,plot3)
setwd(figdir)
#png(paste0('alc.only.cph.png'),height= 350, width= 600,type = 'cairo')
pdf(paste0('alc.only.cph.png'),height= 350, width= 600,type = 'cairo')

print(plot1)
dev.off()

png(paste0('alc.only.cph.female.png'),height= 350, width= 600,type = 'cairo')
print(plot3)
dev.off()


#survival plots on cancer dx with grouped smoking

age.group.col =c('30-40'='#CBF7ED','40-50'='#7B2D26','50-60'='#D7C9AA','60-70'='#0B7A75','70-80'='#19535F')

km_trt_fit.male <- survfit(Surv(censorship_time/365, Event) ~ age_group , data=combined.info.all[combined.info.all$Sex == 'Female',])

plot2 = ggsurvplot(fit=km_trt_fit.female, conf.int = F,ylim = c(0.9,1),censor=F,legend = 'bottom',
                   palette = unname(age.group.col),
                   ggtheme = theme_bw(),
                   legend.title = "Age Group",
                   legend.labs = c("30-40\n(n=11,771)", "40-50\n(n=15,024)","50-60\n(n=18,361)","60-70\n(n=17,211)","70-80\n(n=3670)"),
                   xlab='Time to Cancer Event',
                   ylab='Cancer-Free Fraction') 


#plot.list = list(plot1,plot2,plot3)
setwd(figdir)

png(paste0('age.only.cph.female.png'),height= 350, width= 600,type = 'cairo')
print(plot2)
dev.off()

###bmi short


bmi.short.group.col =c('Normal Weight' = 'grey','Underweight'='#C0F0B9','Overweight'='#5FB49C','Obese'='#414288')

km_trt_fit.female.bmi.short <- survfit(Surv(censorship_time/365, Event) ~ bmi.group.short , data=combined.info.all[combined.info.all$Sex == 'Female',])

plot2 = ggsurvplot(fit=km_trt_fit.female.bmi.short, conf.int = F,ylim = c(0.9,1),censor=F,legend = 'bottom',
                   palette = unname(bmi.short.group.col),
                   ggtheme = theme_bw(),
                   legend.title = "BMI",
                   legend.labs = c("Normal Weight\n(n=47,018)", "Underweight\n(n=1,512)","Overweight\n(n=49,375)","Obese\n(n=33,705)"),
                   xlab='Time to Cancer Event',
                   ylab='Cancer-Free Fraction') 


#plot.list = list(plot1,plot2,plot3)
setwd(figdir)

png(paste0('bmi.short.only.cph.female.png'),height= 350, width= 600,type = 'cairo')
print(plot2)
dev.off()


bmi.long.group.col =c('Normal Weight' = 'grey','Underweight'='#C0F0B9','Overweight'='#5FB49C','Class 1 Obesity'='#414288','Class 2 Obesity'='#682D63','Class 3 Obesity'='#0A1128')

km_trt_fit.female.bmi.long <- survfit(Surv(censorship_time/365, Event) ~ bmi.group.long , data=combined.info.all[combined.info.all$Sex == 'Female',])

plot2 = ggsurvplot(fit=km_trt_fit.female.bmi.long, conf.int = F,ylim = c(0.9,1),censor=F,legend = 'bottom',
                   palette = unname(bmi.long.group.col),
                   ggtheme = theme_bw(),
                   legend.title = "BMI",
                   legend.labs = c("Normal Weight\n(n=47,018)", "Underweight\n(n=1,512)","Overweight\n(n=49,375)","Class 1 Obesity\n(n=21,482)","Class 2 Obesity\n(n=7,640)","Class 3 Obesity\n(n=4,584)"),
                   xlab='Time to Cancer Event',
                   ylab='Cancer-Free Fraction') 


#plot.list = list(plot1,plot2,plot3)
setwd(figdir)

png(paste0('bmi.long.only.cph.female.png'),height= 350, width= 600,type = 'cairo')
print(plot2)
dev.off()


#survival plots on cancer dx with grouped smoking
km_trt_fit.female <- survfit(Surv(censorship_time/365, Event) ~ Smoking.Frequency , data=combined.info.all[combined.info.all$Sex == 'Female',])

smk.group.col =c( Never = 'grey',Former="#C99DA3",Occasional="#996888",Daily="#5E4955")


plot2 = ggsurvplot(fit=km_trt_fit.female, conf.int = F,ylim = c(0.9,1),censor=F,legend = 'bottom',
                   palette = unname(smk.group.col),
                   ggtheme = theme_bw(),
                   legend.title = "Smoking Frequency",
                   legend.labs = c("Never\n(n=31,885)", "Former\n(n=21,360)","Occasional\n(n=1,903)","Daily\n(n=6,064)"),
                   xlab='Time to Cancer Event',
                   ylab='Cancer-Free Fraction') 


#plot.list = list(plot1,plot2,plot3)
setwd(figdir)

png(paste0('smk.only.cph.female.png'),height= 350, width= 600,type = 'cairo')
print(plot2)
dev.off()


#survival plots on cancer dx with grouped family history
km_trt_fit.female <- survfit(Surv(censorship_time/365, Event) ~ Family.history.breast , data=combined.info.all[combined.info.all$Sex == 'Female',])

fh.group.col =c( 'No Reported Family History' = 'grey','Family History' ="#564D80")


plot2 = ggsurvplot(fit=km_trt_fit.female, conf.int = F,ylim = c(0.9,1),censor=F,legend = 'bottom',
                   palette = unname(fh.group.col),
                   ggtheme = theme_bw(),
                   legend.title = "Family History",
                   legend.labs = c("No Reported Family History\n(n=61,306)", "Family History\n(n=4,731)"),
                   xlab='Time to Cancer Event',
                   ylab='Cancer-Free Fraction') 


#plot.list = list(plot1,plot2,plot3)
setwd(figdir)

png(paste0('famhist.brca.only.cph.female.png'),height= 350, width= 600,type = 'cairo')
print(plot2)
dev.off()

#survival plots on cancer dx with grouped family history
km_trt_fit.female <- survfit(Surv(censorship_time/365, Event) ~ Family.history.all , data=combined.info.all[combined.info.all$Sex == 'Female',])

fh.group.col =c( 'No Reported Family History' = 'grey','Family History' ="#564D80")


plot2 = ggsurvplot(fit=km_trt_fit.female, conf.int = F,ylim = c(0.9,1),censor=F,legend = 'bottom',
                   palette = unname(fh.group.col),
                   ggtheme = theme_bw(),
                   legend.title = "Family History",
                   legend.labs = c("No Reported Family History\n(n=61,306)", "Family History\n(n=4,731)"),
                   xlab='Time to Cancer Event',
                   ylab='Cancer-Free Fraction') 


#plot.list = list(plot1,plot2,plot3)
setwd(figdir)

png(paste0('famhist.all.only.cph.male.png'),height= 350, width= 600,type = 'cairo')
print(plot2)
dev.off()


#calculating HR 
#smk
km_trt_fit.male.smk <- coxph(Surv(censorship_time/365, Event) ~ Smoking.Frequency , data=combined.info.all[combined.info.all$Sex == 'Female',])
tmp2 = summary(km_trt_fit.male.smk)
smk.hr.df = data.frame(Var = rep('Smoking Frequency',4),
                       Var.group = c('Never','Former','Occasional','Daily'),
                       HR = c(1,tmp2$coefficients[4:6]),
                       HRL =c(1,tmp2$conf.int[7:9]),
                       HRU =c(1,tmp2$conf.int[10:12]) , 
                       SE=c(0,tmp2$coefficients[7:9]),
                       pvalue = c(1,tmp2$coefficients[13:15]))
#fh
km_trt_fit.male.fh <- coxph(Surv(censorship_time/365, Event) ~ Family.history.breast , data=combined.info.all[combined.info.all$Sex == 'Female',])
tmp2 = summary(km_trt_fit.male.fh)
fh.hr.df = data.frame(Var = rep('Family History',2),
                      Var.group = c('No Reported Family History','Reported Family History'),
                      HR = c(1,tmp2$coefficients[2]),
                      HRL =c(1,tmp2$conf.int[3]),
                      HRU =c(1,tmp2$conf.int[4]) , 
                      SE=c(0,tmp2$coefficients[3]),
                      pvalue = c(1,tmp2$coefficients[5]))

#bmi

km_trt_fit.male.bmi.short <- coxph(Surv(censorship_time/365, Event) ~ bmi.group.short , data=combined.info.all[combined.info.all$Sex == 'Female',])
km_trt_fit.male.bmi.long <- coxph(Surv(censorship_time/365, Event) ~ bmi.group.long , data=combined.info.all[combined.info.all$Sex == 'Female',])

tmp2 = summary(km_trt_fit.male.bmi.short)
bmi.short.hr.df = data.frame(Var = rep('BMI Short',4),
                             Var.group =levels(combined.info.all$bmi.group.short),
                             HR = c(1,tmp2$coefficients[4:6]),
                             HRL =c(1,tmp2$conf.int[7:9]),
                             HRU =c(1,tmp2$conf.int[10:12]) , 
                             SE=c(0,tmp2$coefficients[7:9]),
                             pvalue = c(1,tmp2$coefficients[13:15]))
tmp2 = summary(km_trt_fit.male.bmi.long)
bmi.long.hr.df = data.frame(Var = rep('BMI Long',6),
                            Var.group =levels(combined.info.all$bmi.group.long),
                            HR = c(1,tmp2$coefficients[6:10]),
                            HRL =c(1,tmp2$conf.int[11:15]),
                            HRU =c(1,tmp2$conf.int[16:20]) , 
                            SE=c(0,tmp2$coefficients[11:15]),
                            pvalue = c(1,tmp2$coefficients[21:25]))

#age
km_trt_fit.male.age <- coxph(Surv(censorship_time/365, Event) ~ SDC_AGE_CALC , data=combined.info.all[combined.info.all$Sex == 'Female',])
tmp2 = summary(km_trt_fit.male.age)
age.hr.df = data.frame(Var = rep('Age',1),
                       Var.group = c('Age'),
                       HR = c(tmp2$coefficients[2]),
                       HRL =c(tmp2$conf.int[3]),
                       HRU =c(tmp2$conf.int[4]) , 
                       SE=c(0,tmp2$coefficients[3]),
                       pvalue = c(tmp2$coefficients[5]))

#alc
km_trt_fit.male.alc <- coxph(Surv(censorship_time/365, Event) ~ Alch.con.group , data=combined.info.all[combined.info.all$Sex == 'Female',])
tmp2 = summary(km_trt_fit.male.alc)
alc.hr.df = data.frame(Var = rep('Alcohol Frequency',5),
                       Var.group = c('Never','Former','Infrequent','Moderate','Frequent'),
                       HR = c(1,tmp2$coefficients[5:8]),
                       HRL =c(1,tmp2$conf.int[9:12]),
                       HRU =c(1,tmp2$conf.int[13:16]) , 
                       SE=c(0,tmp2$coefficients[9:12]),
                       pvalue = c(1,tmp2$coefficients[17:20]))

####methylation_score HR ####
discovery.set = readRDS('/.mounts/labs/awadallalab/private/ncheng/manuscripts/early.cancer.risk/discovery.set3.samples.RDS')
validation.set = readRDS('/.mounts/labs/awadallalab/private/ncheng/manuscripts/early.cancer.risk/validation.set3.samples.RDS')
sample.info =  rbind(discovery.set,validation.set)
pred.df.targ = readRDS('/.mounts/labs/awadallalab/private/ncheng/cfmedip_data/cptp_samples/fragmentation/aix13.updated1/all.sample.fragmentation/breast.cancer.cv/regulatory.regions.v10/validation.breast.hrsplit/enhancer.validation.final.90/validation.brca.scores.RDS')
#pred.df.targ = readRDS('/.mounts/labs/awadallalab/private/ncheng/cfmedip_data/cptp_samples/fragmentation/aix13.updated1/all.sample.fragmentation/breast.cancer.cv/regulatory.regions.v10.newsplit7//cv.figures.genhancer.final.90///discovery.brca.scores.RDS')

score.dmr = pred.df.targ
score.dmr$GRP_Id = gsub('_combined.*','',score.dmr$GRP_Id)
top.perf.df = pred.df.targ
top.perf.df$GRP_Id = gsub('_combined.*','',top.perf.df$GRP_Id)
pred.df.targ.collapse.all = ddply(top.perf.df[,c('GRP_Id','reported','methylation_score')],c('GRP_Id','reported'),numcolwise(mean))
pred.df.targ.collapse.all = pred.df.targ.collapse.all[pred.df.targ.collapse.all$GRP_Id %in% sample.info$GRP_Id,]
#pred.df.targ.collapse.all = readRDS('/.mounts/labs/awadallalab/private/ncheng/cfmedip_data/cptp_samples/fragmentation/aix13.updated1/methylation.insert.select1/predictions1/top.perf.discovery.cv.RDS')
cp.youden= cutpointr(pred.df.targ.collapse.all$methylation_score,pred.df.targ.collapse.all$reported, method = maximize_metric, metric = youden)$optimal_cutpoint

mean.perf.df=pred.df.targ
mean.perf.df.targ.tmp = pred.df.targ.collapse.all#[!pred.df.targ.collapse.all$GRP_Id %in% validation.samples$GRP_Id,]
cp.youden= cutpointr(mean.perf.df.targ.tmp$methylation_score,mean.perf.df.targ.tmp$reported, method = maximize_metric, metric = youden)$optimal_cutpoint

y.index = 0.40 #cp.youden
mean.perf.df.targ.tmp.merged = merge(mean.perf.df.targ.tmp, sample.info[,c('GRP_Id','censorship_time')],by='GRP_Id')
mean.perf.df.targ.tmp.merged$Event=ifelse(mean.perf.df.targ.tmp.merged$reported == 'Control',0,1)
mean.perf.df.targ.tmp.merged$group = mean.perf.df.targ.tmp.merged$reported
mean.perf.df.targ.tmp.merged = mean.perf.df.targ.tmp.merged[!is.na(mean.perf.df.targ.tmp.merged$methylation_score),]
mean.perf.df.targ.tmp.merged$Risk.group = ifelse(mean.perf.df.targ.tmp.merged$methylation_score >= y.index,'High Predicted Risk','Low Predicted Risk')
mean.perf.df.targ.tmp.merged$Risk.group = factor(as.character(mean.perf.df.targ.tmp.merged$Risk.group ),levels = c('Low Predicted Risk','High Predicted Risk'))
mean.perf.df.targ.tmp.merged$Event= ifelse(mean.perf.df.targ.tmp.merged$reported == 'Control',0,1 )

female.ohs.qx.weighting.bmiadjust = function(subject,combined.info.all) {
  #subject = merge(mean.perf.df.targ.tmp.merged,sample.info[,c('GRP_Id','age_group','Family.history.breast','Alch.con.group','Smoking.Frequency')],by='GRP_Id')
  combined.info.all.male = combined.info.all[combined.info.all$Sex == 'Female',]
  combined.info.all.male = combined.info.all.male[combined.info.all.male$Cancer %in% c('Control','Breast'),]
  combined.info.all.male$Smoking.Frequency = as.character(combined.info.all.male$Smoking.Frequency)
  combined.info.all.male[is.na(combined.info.all.male$Smoking.Frequency),'Smoking.Frequency'] = 'Never'
  bmi.short.groups = unique(combined.info.all.male$bmi.group.short)
  #bmi.long.groups = unique(combined.info.all.male$bmi.group.long)
  age.groups =  unique(combined.info.all.male$age_group)
  fh.groups = unique(combined.info.all.male$Family.history.breast)
  alc.groups = unique(combined.info.all.male$Alch.con.group)
  smk.groups = unique(combined.info.all.male$Smoking.Frequency)
  control.samples = subject[subject$reported == 'Control',]
  cancer.samples = subject[subject$reported != 'Control',]
  
  control.samples = combined.info.all[combined.info.all$GRP_Id %in% control.samples$GRP_Id,]
  cancer.samples = combined.info.all[combined.info.all$GRP_Id %in% cancer.samples$GRP_Id,]
  
  ohs.pop.samples = combined.info.all[combined.info.all$Cohort == 'EpiCan', ]
  control.weight.samples = NULL
  for (bmi.short in bmi.short.groups){
    for (age in age.groups ) {
      for (fh in fh.groups ){
        for (alc in alc.groups) {
          for (smk in smk.groups ){
            group.combination = paste(age,fh,alc,smk,sep='.')
            targ.control.cohort = control.samples[which(control.samples$age_group == age &
                                                          control.samples$Family.history.breast == fh &
                                                          control.samples$Alch.con.group == alc &
                                                          control.samples$Smoking.Frequency == smk &
                                                          control.samples$bmi.group.short == bmi.short) ,]
            if (nrow(targ.control.cohort)>0){
              cohort.freq =  ohs.pop.samples[which(ohs.pop.samples$age_group == age &
                                                     ohs.pop.samples$Family.history.breast == fh &
                                                     ohs.pop.samples$Alch.con.group == alc &
                                                     ohs.pop.samples$Smoking.Frequency == smk &
                                                     ohs.pop.samples$bmi.group.short == bmi.short ),]
              
              targ.control.cohort$weight= nrow(cohort.freq)/nrow(targ.control.cohort)
              
              control.weight.samples = rbind(control.weight.samples,targ.control.cohort)
            }
            
            
          }
        }
      }
    }
    
  }
  
  if (nrow(control.weight.samples) < nrow(control.samples)) {
    exc.samples = control.samples[!control.samples$GRP_Id %in% control.weight.samples$GRP_Id,]
    for (r in 1:nrow(exc.samples)) {
      print(r)
      targ.sample = exc.samples[r,]
      #fh, alc, smk, bmi
      non.na.vars = c('age_group','ALC_CUR_FREQ','SMK_CIG_STATUS','DIS_CANCER_F_EVER','bmi.group.short')
      non.na.vars = non.na.vars[which(!is.na(targ.sample[,non.na.vars]))]
      
      cohort.freq =  ohs.pop.samples
      for (v in non.na.vars) {
        cohort.freq =  cohort.freq[which(cohort.freq[,v] == targ.sample[,v]) ,]
        targ.freq = control.samples[which(control.samples[,v] == targ.sample[,v])  ,]
      }
      targ.sample$weight = nrow(cohort.freq)/nrow(targ.freq)
      control.weight.samples = rbind(control.weight.samples,targ.sample)
    }
  }
  cancer.samples$weight=1
  return.df = rbind(control.weight.samples,cancer.samples)
  return.df = return.df[order(match(return.df$GRP_Id,subject$GRP_Id)),]
  return(return.df$weight)
  
}
female.ohs.qx.weighting = function(subject,combined.info.all) {
  #subject = merge(mean.perf.df.targ.tmp.merged,sample.info[,c('GRP_Id','age_group','Family.history.breast','Alch.con.group','Smoking.Frequency')],by='GRP_Id')
  combined.info.all.male = combined.info.all[combined.info.all$Sex == 'Female',]
  combined.info.all.male = combined.info.all.male[combined.info.all.male$Cancer %in% c('Control','Breast'),]
  combined.info.all.male$Smoking.Frequency = as.character(combined.info.all.male$Smoking.Frequency)
  combined.info.all.male[is.na(combined.info.all.male$Smoking.Frequency),'Smoking.Frequency'] = 'Never'
  #bmi.short.groups = unique(combined.info.all.male$bmi.group.short)
  #bmi.long.groups = unique(combined.info.all.male$bmi.group.long)
  age.groups =  unique(combined.info.all.male$age_group)
  fh.groups = unique(combined.info.all.male$Family.history.breast)
  alc.groups = unique(combined.info.all.male$Alch.con.group)
  smk.groups = unique(combined.info.all.male$Smoking.Frequency)
  control.samples = subject[subject$reported == 'Control',]
  cancer.samples = subject[subject$reported != 'Control',]
  
  control.samples = combined.info.all[combined.info.all$GRP_Id %in% control.samples$GRP_Id,]
  cancer.samples = combined.info.all[combined.info.all$GRP_Id %in% cancer.samples$GRP_Id,]
  
  ohs.pop.samples = combined.info.all[combined.info.all$Cohort == 'EpiCan', ]
  control.weight.samples = NULL
  for (age in age.groups ) {
    for (fh in fh.groups ){
      for (alc in alc.groups) {
        for (smk in smk.groups ){
          group.combination = paste(age,fh,alc,smk,sep='.')
          targ.control.cohort = control.samples[which(control.samples$age_group == age &
                                                        control.samples$Family.history.breast == fh &
                                                        control.samples$Alch.con.group == alc  &
                                                        control.samples$Smoking.Frequency == smk ) ,]
          if (nrow(targ.control.cohort)>0){
            cohort.freq =  ohs.pop.samples[which(ohs.pop.samples$age_group == age &
                                                   ohs.pop.samples$Family.history.breast == fh &
                                                   ohs.pop.samples$Alch.con.group == alc &
                                                   ohs.pop.samples$Smoking.Frequency == smk ),]
            
            targ.control.cohort$weight= nrow(cohort.freq)/nrow(targ.control.cohort)
            
            control.weight.samples = rbind(control.weight.samples,targ.control.cohort)
          }
          
          
        }
      }
    }
  }
  
  if (nrow(control.weight.samples) < nrow(control.samples)) {
    exc.samples = control.samples[!control.samples$GRP_Id %in% control.weight.samples$GRP_Id,]
    for (r in 1:nrow(exc.samples)) {
      print(r)
      targ.sample = exc.samples[r,]
      #fh, alc, smk, bmi
      non.na.vars = c('age_group','ALC_CUR_FREQ','SMK_CIG_STATUS','DIS_CANCER_F_EVER')
      non.na.vars = non.na.vars[which(!is.na(targ.sample[,non.na.vars]))]
      
      cohort.freq =  ohs.pop.samples
      for (v in non.na.vars) {
        cohort.freq =  cohort.freq[which(cohort.freq[,v] == targ.sample[,v]) ,]
        targ.freq = control.samples[which(control.samples[,v] == targ.sample[,v])  ,]
      }
      targ.sample$weight = nrow(cohort.freq)/nrow(targ.freq)
      control.weight.samples = rbind(control.weight.samples,targ.sample)
    }
  }
  cancer.samples$weight=1
  return.df = rbind(control.weight.samples,cancer.samples)
  return.df = return.df[order(match(return.df$GRP_Id,subject$GRP_Id)),]
  return(return.df$weight)
  
}


#mean.perf.df.targ.tmp.merged = unique(mean.perf.df.targ.tmp.merged[!mean.perf.df.targ.tmp.merged$GRP_Id %in% validation.samples$GRP_Id,])
mean.perf.df.targ.tmp.merged = unique(mean.perf.df.targ.tmp.merged[mean.perf.df.targ.tmp.merged$GRP_Id %in% combined.info.all$GRP_Id,])
km_trt_fit.male.medip <- coxph(Surv(censorship_time/365, Event) ~ Risk.group , data=mean.perf.df.targ.tmp.merged,
                               weights = female.ohs.qx.weighting(unique(mean.perf.df.targ.tmp.merged), combined.info.all[combined.info.all$Cancer %in% c('Control','Breast'),]))
#km_trt_fit.male.medip <- coxph(Surv(censorship_time/365, Event) ~ Risk.group , data=mean.perf.df.targ.tmp.merged)
fit_spline <- cph(Surv(censorship_time/365, Event) ~ rcs(Risk.group,4) ,
                  data=mean.perf.df.targ.tmp.merged,
                  weights = female.ohs.qx.weighting(unique(mean.perf.df.targ.tmp.merged), combined.info.all[combined.info.all$Cancer %in% c('Control','Breast'),]))


tmp2 = summary(km_trt_fit.male.medip)
methscore.hr.df = data.frame(Var = rep('cfDNA Methylation Risk Score',2),
                             Var.group = c('Low Predicted Risk','High Predicted Risk'),
                             HR = c(1,tmp2$coefficients[1,2]),
                             HRL =c(1,tmp2$conf.int[1,3]),
                             HRU =c(1,tmp2$conf.int[1,4]) ,
                             SE=c(0,tmp2$coefficients[1,3]),
                             pvalue = c(1,tmp2$coefficients[1,5]))

combined.hr =rbind(smk.hr.df, alc.hr.df,fh.hr.df,age.hr.df, methscore.hr.df,bmi.long.hr.df,bmi.short.hr.df)
#saveRDS(combined.hr,'~/tmp/combined.covariate.hr.RDS')
saveRDS(combined.hr,'~/tmp/validation.covariate.hr.RDS')

combined.hr.qx = readRDS('/Users/ncheng/Desktop/ncheng/Desktop/Lab Work/Manuscripts/PhD Thesis /Chapter3/files/combined.covariate.hr.RDS')
combined.hr.qx = combined.hr.qx[combined.hr.qx$Var != "cfDNA Methylation Risk Score",]
setwd('/Users/ncheng/Desktop/ncheng/Desktop/Lab Work/Manuscripts/OHS_cancer_risk/figures/ml.predictions/hr.risk/')
#combined.hr = readRDS('combined.covariate.hr.RDS')
validation.hr = readRDS('validation.covariate.hr.RDS')
validation.hr = validation.hr[grepl('Methylation',validation.hr$Var),]
combined.hr = rbind(validation.hr,combined.hr.qx)
combined.hr$group = 'Discovery'
validation.hr$group='Validation'
combined.hr$Var = ifelse(combined.hr$Var == 'cfDNA Methylation Risk Score','cfDNA Methylation Risk Score (Discovery)',combined.hr$Var)
combined.hr$Var = ifelse(combined.hr$Var == 'cfDNA Methylation Risk Score','cfDNA Methylation Risk Score (Test Set)',combined.hr$Var)

validation.hr$Var = ifelse(validation.hr$Var == 'cfDNA Methylation Risk Score','cfDNA Methylation Risk Score (Test)',validation.hr$Var)

combined.hr= rbind(combined.hr,validation.hr)
combined.hr$SE1 = (combined.hr$HR - combined.hr$HRL)/1.96
combined.hr$est.log = log(combined.hr$HR)
require("survival")
library(ggforestplot)
library(ggplot2)
combined.hr = combined.hr[combined.hr$Var %in% c('Age','Family History','Alcohol Frequency','cfDNA Methylation Risk Score (Discovery)','cfDNA Methylation Risk Score (Test Set)'),]
combined.hr$Var = factor(as.character(combined.hr$Var), levels = c('Age','Family History','Smoking Frequency','Alcohol Frequency','BMI Short','cfDNA Methylation Risk Score (Discovery)','cfDNA Methylation Risk Score (Test Set)','BMI Long'))


plot1 =ggforestplot::forestplot(
  df = combined.hr,
  name = Var.group,
  estimate = est.log,
  pvalue =pvalue,
  se =SE,
  logodds=T,
  psignif = 0.05,
  xlab='Hazard Ratio (95% CI)',
  ylab='',
  #xlab = "1-SD increment in cardiometabolic trait\nper 1-SD increment in biomarker concentration",
  colour = Var
) +
  ggforce::facet_col(
    facets = ~Var,
    scales = "free_y",
    space = "free"
  ) + 
  #scale_x_continuous(breaks = seq(0,16,1),limits = c(0.00001,15)) +
  theme_bw()+
  theme(
    legend.position='none',
    text = element_text(size=8),
    axis.text=element_text(size=8, face = "bold"),
    axis.title=element_text(size=8,face="bold"),
    strip.background = element_rect(fill = "white"))+scale_color_manual(values =c('#143642','#0F8B8D','#EC9A29','#A8201A','#A270D6'))

print(plot1)
png('test.png',height=1500,width=1600,res=300)
print(plot1)
dev.off()
scipen(options = 1000)
model <- coxph(Surv(censorship_time/365, Event) ~ Risk.group +  Alch.con.group+SDC_AGE_CALC+Family.history.breast+Smoking.Frequency, 
               data=mean.perf.df.targ.tmp,
               weights = weightsf.males(mean.perf.df.targ.tmp))

model <- coxph( Surv(time, status) ~ sex + rx + adhere,
                data = colon )
ggforest(model)

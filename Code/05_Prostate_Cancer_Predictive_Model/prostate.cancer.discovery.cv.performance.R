#!/usr/bin/env Rscript
####loading packages#####
library(plyr)
library(ROCR)
library(rstatix)
library(DESeq2)
library(matrixTests)
library(gridExtra)
library(caret)
library(pROC)
library(stringr)
library(parallel)
library(survival)
library(ncvreg)
library(survminer)
library(RColorBrewer)
library(cutpointr)
library(dplyr)
library(plyr)
library(ggplot2)
library(ggh4x)
library(survcomp)
library(cenROC)

####functions to be used in analysis####
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
    #if male
    s <- all_cause_mortality_CAN$lx
    
  } else {
    #if female
    s <- all_cause_mortality_CAN$lx.1
    
  }
  f <- function(x) exp(splinefun(all_cause_mortality_CAN$x, log(s), method="mono")(x))
  return(f(age2) / f(age1))
  
  
}

aml_inc_cr <- Vectorize(function(gender, age1, age2) sum(diff(aml_inc(gender, seq(age1,age2,1) ))*all_surv(gender, age1, seq(age1,age2-1,1)) ), c("gender","age1","age2"))


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

####results analysis####
sample.info = readRDS('sample.info.RDS') 
sample.info = sample.info[sample.info$Sex == 'Male',]
wkdir='/wkdir/'
marker.list = c('silencer')

median.perf.list = list()
mean.perf.list = list()
for (m in marker.list){
  print(m)
  
  directions = list.files(pattern='predictions.*',path =paste0(wkdir,m,'/'))
  
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
            tmp$covariate = 'none'
            tmp$fc= 0.25
            tmp$direction = dir
            return(tmp)
            
            
          },mc.cores=5 )
          pred.df.all.tmp=do.call('rbind',pred.list)
          feature.size.count = data.frame(table(pred.df.all.tmp$n.features),stringsAsFactors = F)
          pred.df.all =rbind(pred.df.all, unique(pred.df.all.tmp))
          
        }
        
      }
      
    }
    pred.df.all$freq = 1
    #computing cross-validated classification score by averaging performance across CV repeats grouped by model, tuning paramters and feature numbers
    pred.df.all.freq =ddply(pred.df.all[,c('GRP_Id','reported','model','feature','test','n.features','direction','freq','covariate','dir','fc')],
                            c('GRP_Id','reported','model','feature','test','n.features','direction','covariate','dir','fc'),
                            numcolwise(sum))
    
    
    pred.df.collapse.mean = ddply(pred.df.all[,c('GRP_Id','reported','methylation_score','model','feature','test','n.features','direction','covariate','dir','fc')],
                                  c('GRP_Id','reported','model','feature','test','n.features','direction','covariate','dir','fc'),
                                  numcolwise(mean))
    
    pred.df.collapse.mean =merge(pred.df.collapse.mean,pred.df.all.freq,by=c('GRP_Id','reported','model','feature','test','n.features','mat','direction','covariate','dir','fc'))
    
  }
  
}
figdir = paste0(wkdir,'/discovery.cv.figures/')
dir.create(figdir,recursive = T)
perfdir= paste0(wkdir,'/')
perf.list = list('mean.cv'= pred.df.collapse.mean)


#computing concordance and auc from mean CV risk scores 
targ.perf = pred.df.collapse.mean
targ.perf.list = split(targ.perf, targ.perf[,c('model','feature','test','n.features','direction','covariate','fc')])
targ.perf.list = targ.perf.list[sapply(targ.perf.list,nrow)>0 ]
targ.perf.list =lapply(targ.perf.list, function(x) {
  tmp = x
  tmp$All = auc_calc(tmp,c('Control','Cancer'))
  concordance_calc = function(tmp,sample.info){
    x = merge(tmp,sample.info)
    x$event=ifelse(x$reported =='Cancer',1,0)
    x$reported.surv = ifelse(x$reported == 'Cancer',1,0)
    ci= concordance.index(x$methylation_score, x$'censorship_time', surv.event = x$event, comppairs=10, na.rm = FALSE)
    return(ci$c.index)
  }
  tmp$Concordance = concordance_calc(tmp,sample.info)
  return(tmp)
} )

overall.perf=do.call('rbind',targ.perf.list)
overall.perf = overall.perf[order(-overall.perf$All),]
rownames(overall.perf) = c(1:nrow(overall.perf))

####selecting top performing parameter according to best CV auroc

target.parameters = overall.perf#[overall.perf$n.features == 90 ,] 90 features used in manuscript
target.parameters = target.parameters[order(-target.parameters$All),][1,]

overall.perf.targ.prostate = overall.perf[overall.perf$n.features == target.parameters$n.features & 
                                            overall.perf$direction == target.parameters$direction & 
                                            overall.perf$feature == target.parameters$feature & 
                                            overall.perf$mat == target.parameters$mat &
                                            overall.perf$model == target.parameters$model, ]

score.cutoff = cutpointr(overall.perf.targ.prostate$methylation_score,overall.perf.targ.prostate$reported, method = maximize_metric, metric = youden)
score.cutoff.prostate = score.cutoff$optimal_cutpoint


####prostate cancer plotting#####
combined.info.all = readRDS('male.epican.weighting.info.RDS') #upload
pred.df.targ = overall.perf.targ.prostate
name = paste0(figdir,'prad.top');
merged.df.all= sample.inof;
score.cutoff=score.cutoff.prostate;
dx.all = F; 
cutpoint.use =T; #using Youden's cutoff
gs.plot= T #plotting gleason scores

#initally embedded as a function, but survminer has trouble reading in some data, so try to just run with the chunk of code in braces below
{

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
  
  #sample.info.filt.pretime = merged.df.all
  pred.df.targ.collapse.all = pred.df.targ
  pred.df.targ.collapse.all$model = gsub('.new','',pred.df.targ.collapse.all$model)

  pred.df.targ.collapse.list = split(pred.df.targ.collapse.all,pred.df.targ.collapse.all$model)
  pred.df.targ.collapse.list = lapply(pred.df.targ.collapse.list, function(x) {
    auc.plot.all = tpr.fpr.calc(x)
    auc.plot.all$auc = auc_calc(x,labels= c('Control','Cancer'))
    auc.plot.all$model = x$model[1]
    return(auc.plot.all)
  } )
  
  pred.df.targ.collapse.df = do.call('rbind',pred.df.targ.collapse.list)
  summary.calc = function(pred.df.targ.collapse,merged.df.all.tmp,weighting) {
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
        male.weights= weightsf.males(x)
        
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
      
      male.ohs.qx.weighting = function(subject1,combined.info.all) {
        subject = merge(subject1, all.sample.info[,c('GRP_Id','Smoking.Frequency','age_group','Family.history.prostate','Alch.con.group')],by='GRP_Id')
        combined.info.all.male = combined.info.all[combined.info.all$Sex == 'Male',]
        combined.info.all.male$Smoking.Frequency = as.character(combined.info.all.male$Smoking.Frequency)
        combined.info.all.male[is.na(combined.info.all.male$Smoking.Frequency),'Smoking.Frequency'] = 'Never'
        age.groups =  unique(combined.info.all.male$age_group)
        fh.groups = unique(combined.info.all.male$Family.history.prostate)
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
                                                        control.samples$Family.history.prostate == fh &
                                                        control.samples$Alch.con.group == alc &
                                                        control.samples$Smoking.Frequency == smk ,]
                if (nrow(targ.control.cohort)>0){
                  cohort.freq =  ohs.pop.samples[ohs.pop.samples$age_group == age &
                                                   ohs.pop.samples$Family.history.prostate == fh &
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
      
      mean.perf.df.targ.tmp.merged$weighting = male.ohs.qx.weighting(mean.perf.df.targ.tmp.merged, combined.info.all)
      
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
        male.weights= weightsf.males(x)
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
        male.weights= weightsf.males(x)
        
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
      subject1 = unique(mean.perf.df.targ.tmp.merged)
      male.ohs.qx.weighting = function(subject1,combined.info.all) {
        subject = merge(subject1, all.sample.info[,c('GRP_Id','Smoking.Frequency','age_group','Family.history.prostate','Alch.con.group')],by='GRP_Id')
        combined.info.all.male = combined.info.all[combined.info.all$Sex == 'Male',]
        combined.info.all.male$Smoking.Frequency = as.character(combined.info.all.male$Smoking.Frequency)
        combined.info.all.male[is.na(combined.info.all.male$Smoking.Frequency),'Smoking.Frequency'] = 'Never'
        age.groups =  unique(combined.info.all.male$age_group)
        fh.groups = unique(combined.info.all.male$Family.history.prostate)
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
                                                        control.samples$Family.history.prostate == fh &
                                                        control.samples$Alch.con.group == alc &
                                                        control.samples$Smoking.Frequency == smk ,]
                if (nrow(targ.control.cohort)>0){
                  cohort.freq =  ohs.pop.samples[ohs.pop.samples$age_group == age &
                                                   ohs.pop.samples$Family.history.prostate == fh &
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
      
      mean.perf.df.targ.tmp.merged$weighting = male.ohs.qx.weighting(mean.perf.df.targ.tmp.merged, combined.info.all)
      
      
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
  
  #grade
  #stage
  #last psa 
  #age
  
  #overall performance between models
  return.df.final.auroc = NULL
  return.df.final.auct = NULL
  for (c in unique(pred.df.targ.collapse.df$model)) {
    library(cutpointr)
    
    print('ck1')
    
    if (cutpoint.use == F) {
      cp.youden= cutpointr(pred.df.targ.collapse.all$methylation_score,pred.df.targ.collapse.all$reported, method = maximize_metric, metric = youden)$optimal_cutpoint
      cp.F1_score= cutpointr(pred.df.targ.collapse.all$methylation_score,pred.df.targ.collapse.all$reported, method = maximize_metric, metric = F1_score)$optimal_cutpoint
      print(cp.youden)
    } else {
      cp.youden= score.cutoff
      cp.F1_score= cutpointr(pred.df.targ.collapse.all$methylation_score,pred.df.targ.collapse.all$reported, method = maximize_metric, metric = F1_score)$optimal_cutpoint
      print(cp.youden)
    }
    combined.auroc = NULL
    auc.t.combined = NULL
    pred.df.targ.collapse =pred.df.targ.collapse.all#[pred.df.targ.collapse.all$model == as.character(c),]
    merged.df.all.tmp = merge(pred.df.targ.collapse, merged.df.all,by=c('GRP_Id','ResearchId'))
    
    if (dx.all == F) {
      #plotting dx time grouped by 2 years
      dx12 = merged.df.all.tmp[merged.df.all.tmp$Diagnosis_Time %in% c('Control','0-1','1-2','0-2'),]
      dx34 = merged.df.all.tmp[merged.df.all.tmp$Diagnosis_Time %in% c('Control','2-3','3-4','2-4'),]
      dx45 = merged.df.all.tmp[merged.df.all.tmp$Diagnosis_Time %in% c('Control','4-5','5+','4-6'),]
      
      auc.plot.all = tpr.fpr.calc(pred.df.targ.collapse)
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
      auc.plot.all.dx45$dx.time = '4+'
      auc.plot.all.dx45$auc = auc_calc(pred.df.targ.collapse[pred.df.targ.collapse$GRP_Id %in% dx45$GRP_Id,],c('Control','Cancer'))
      
      auc.res1 = rbind(auc.plot.all[1,c('auc','dx.time')],
                       auc.plot.all.dx12[1,c('auc','dx.time')],
                       auc.plot.all.dx34[1,c('auc','dx.time')],
                       auc.plot.all.dx45[1,c('auc','dx.time')])
      colnames(auc.res1) = c('auc','subgroup')
      diagnosis_time_colors1 = c('#7A797C',"#048BA8",'#AAF683','#FFD97D')
      names(diagnosis_time_colors1) = c('All','0-2','2-4','4+')
      combined.auc.plot.all = rbind(auc.plot.all,auc.plot.all.dx12,auc.plot.all.dx34,auc.plot.all.dx45)
      plot1 = ggplot(combined.auc.plot.all,aes(x = fpr, y = tpr, col = dx.time)) + geom_line(linewidth=1) +
        scale_color_manual(values = c(diagnosis_time_colors1))+ #+ ggtitle(title) +
        theme_bw()+ 
        theme(text = element_text(size=8),
              axis.text=element_text(size=8),
              axis.title=element_text(size=8,face="bold"),
              legend.position = "none")+ guides(colour = guide_legend(override.aes = list(size=4),nrow=4)) + xlab('False Positive Rate') + ylab('Sensitivity')
      
      png(paste0(name,'.',c,'.auroc.dxtime2.png'),height = 1000,width=1000,res = 300)
      print(plot1)
      dev.off()
      
      labels = combined.auc.plot.all[combined.auc.plot.all$l ==0,]
      labels$title = gsub('.new','',paste0(labels$dx.time, ' AUROC: ',round(labels$auc,digits=3)))
      colors.filt = diagnosis_time_colors1[names(diagnosis_time_colors1) %in% labels$dx.time]
      labels = labels[order(match(labels$dx.time,names(colors.filt) )),]  
      labels$dx.time = factor(labels$dx.time,levels = names(colors.filt))
      combined.auc.plot.all$dx.time = factor(as.character(combined.auc.plot.all$dx.time), levels= names(colors.filt))
      
      
      plot1 = ggplot(combined.auc.plot.all,aes(x = fpr, y = tpr, col = dx.time)) + geom_line(linewidth=1) +
        scale_color_manual(values = c(diagnosis_time_colors1),labels = labels$title)+ #+ ggtitle(title) +
        theme_bw()+ 
        theme(text = element_text(size=8),
              axis.text=element_text(size=8, face = "bold"),
              axis.title=element_text(size=8,face="bold"),
              legend.position = "right")+ guides(colour = guide_legend(override.aes = list(size=4),ncol=1)) + xlab('False Positive Rate') + ylab('Sensitivity')
      
      png(paste0(name,'.',c,'.auroc.dxtim2e.labs.png'),height = 500,width=1000,res = 300)
      print(plot1)
      dev.off()
      
    } else {
      #plotting dx time grouped by 2 years
      dx12 = merged.df.all.tmp[merged.df.all.tmp$Diagnosis_Time %in% c('Control','0-1','1-2','0-2'),]
      dx34 = merged.df.all.tmp[merged.df.all.tmp$Diagnosis_Time %in% c('Control','2-3','3-4','2-4'),]
      dx45 = merged.df.all.tmp[merged.df.all.tmp$Diagnosis_Time %in% c('Control','4-5','5+','5-6','4-6'),]
      dx67 = merged.df.all.tmp[merged.df.all.tmp$Diagnosis_Time %in% c('Control','6-7','7-8','6-8'),]
      dx89 = merged.df.all.tmp[merged.df.all.tmp$Diagnosis_Time %in% c('Control','8-9','9-10','8-10'),]
      print('dx tmie 2 years')
      
      auc.plot.all = summary.calc(pred.df.targ.collapse,merged.df.all.tmp)
      auc.plot.all[[2]]$dx.time = 'All'
      auc.plot.all[[2]]$var = 'All'
      auc.plot.all[[2]]$var.group = 'All'
      
      auc.plot.all.dx12 = summary.calc.dxtime(pred.df.targ.collapse[pred.df.targ.collapse$GRP_Id %in% dx12$GRP_Id,],merged.df.all.tmp)
      auc.plot.all.dx12[[2]]$dx.time = '0-2'
      auc.plot.all.dx12[[2]]$var = '0-2'
      auc.plot.all.dx12[[2]]$var.group = 'DxTime2'
      
      auc.plot.all.dx34 = summary.calc.dxtime(pred.df.targ.collapse[pred.df.targ.collapse$GRP_Id %in% dx34$GRP_Id,],merged.df.all.tmp)
      auc.plot.all.dx34[[2]]$dx.time = '2-4'
      auc.plot.all.dx34[[2]]$var = '2-4'
      auc.plot.all.dx34[[2]]$var.group = 'DxTime2'
      
      auc.plot.all.dx45 = summary.calc.dxtime(pred.df.targ.collapse[pred.df.targ.collapse$GRP_Id %in% dx45$GRP_Id,],merged.df.all.tmp)
      auc.plot.all.dx45[[2]]$dx.time = '4-6'
      auc.plot.all.dx45[[2]]$var = '4-6'
      auc.plot.all.dx45[[2]]$var.group = 'DxTime2'
      
      auc.plot.all.dx67 = summary.calc.dxtime(pred.df.targ.collapse[pred.df.targ.collapse$GRP_Id %in% dx67$GRP_Id,],merged.df.all.tmp)
      auc.plot.all.dx67[[2]]$dx.time =  '6-8'
      auc.plot.all.dx67[[2]]$var = '6-8'
      auc.plot.all.dx67[[2]]$var.group = 'DxTime2'
      
      auc.plot.all.dx89 = summary.calc.dxtime(pred.df.targ.collapse[pred.df.targ.collapse$GRP_Id %in% dx89$GRP_Id,],merged.df.all.tmp)
      auc.plot.all.dx89[[2]]$dx.time = '8-10'
      auc.plot.all.dx89[[2]]$var = '8-10'
      auc.plot.all.dx89[[2]]$var.group = 'DxTime2'
      
      
      diagnosis_time_colors1 = c('grey',"#048BA8",'#AAF683','#FFD97D','#FF9B85','#C8553D')
      names(diagnosis_time_colors1) = c('All','0-2','2-4','4-6','6-8','8-10')
      
      
      auc.res2 = rbind(auc.plot.all[[2]][1,],
                       auc.plot.all.dx12[[2]][1,],
                       auc.plot.all.dx34[[2]][1,],
                       auc.plot.all.dx45[[2]][1,],
                       auc.plot.all.dx67[[2]][1,])
      #  auc.plot.all.dx89[[2]][1,])
      auc.res2$title = 'DxTime2'
      combined.auroc = rbind(combined.auroc,auc.res2)
      
      combined.auc.plot.all = rbind(auc.plot.all[[2]],
                                    auc.plot.all.dx12[[2]],
                                    auc.plot.all.dx34[[2]],
                                    auc.plot.all.dx45[[2]],
                                    auc.plot.all.dx67[[2]])#,
      # auc.plot.all.dx89[[2]])
      
      plot1 = ggplot(combined.auc.plot.all,aes(x = fpr, y = tpr, col = var)) + geom_line(linewidth=1) +
        scale_color_manual(values = c(diagnosis_time_colors1))+ #+ ggtitle(title) +
        theme_bw()+ 
        theme(text = element_text(size=8),
              axis.text=element_text(size=8, face = "bold"),
              axis.title=element_text(size=8,face="bold"),
              legend.position = "none")+ guides(colour = guide_legend(override.aes = list(size=4),nrow=4)) + xlab('False Positive Rate') + ylab('Sensitivity')
      
      png(paste0(name,'.',c,'.auroc.dxtime2.png'),height = 1000,width=1000,res = 300)
      print(plot1)
      dev.off()
      
      labels = combined.auc.plot.all[combined.auc.plot.all$l ==0,]
      labels$title = gsub('.new','',paste0(labels$dx.time, ' AUROC: ',round(labels$auc,digits=3)))
      colors.filt = diagnosis_time_colors1[names(diagnosis_time_colors1) %in% labels$dx.time]
      labels = labels[order(match(labels$dx.time,names(colors.filt) )),]  
      labels$dx.time = factor(labels$dx.time,levels = names(colors.filt))
      combined.auc.plot.all$dx.time = factor(as.character(combined.auc.plot.all$dx.time), levels= names(colors.filt))
      
      
      plot1 = ggplot(combined.auc.plot.all,aes(x = fpr, y = tpr, col = var)) + geom_line(linewidth=1) +
        scale_color_manual(values = c(diagnosis_time_colors1),labels = labels$title)+ #+ ggtitle(title) +
        theme_bw()+ 
        theme(text = element_text(size=8),
              axis.text=element_text(size=8, face = "bold"),
              axis.title=element_text(size=8,face="bold"),
              legend.position = "right")+ guides(colour = guide_legend(override.aes = list(size=4),ncol=1)) + xlab('False Positive Rate') + ylab('Sensitivity')
      
      png(paste0(name,'.',c,'.auroc.dxtim2e.labs.png'),height = 500,width=1000,res = 300)
      print(plot1)
      dev.off()
      
      
      
      
      
      plot1 = ggplot(auc.res2,aes(x = dx.time, y = auc, col = dx.time)) +
        geom_errorbar(aes(ymin=auc.lower, ymax=auc.upper), width=.2,
                      position=position_dodge(0.05))+
        geom_point()+
        scale_y_continuous(limits=c(0,1))+
        scale_color_manual(values = c(diagnosis_time_colors1))+ #+ ggtitle(title) +
        theme_bw()+ 
        #facet_grid(.~title)+
        
        theme(text = element_text(size=8),
              axis.text=element_text(size=8, face = "bold"),
              axis.title=element_text(size=8,face="bold"),
              legend.position = "none")+ 
        guides(colour = guide_legend(override.aes = list(size=4),nrow=4)) +
        xlab('Sampling time prior to diagnosis') + ylab('CV AUROC (95% CI)')
      
      png(paste0(name,'.',c,'.auroc.ci.dxtime2.png'),height = 1100,width=1100,res = 300)
      print(plot1)
      dev.off()
      plot1 = ggplot(auc.res2,aes(x = dx.time, y = ci, col = dx.time)) +
        geom_errorbar(aes(ymin=ci.lower, ymax=ci.upper), width=.2,
                      position=position_dodge(0.05))+
        geom_point()+
        scale_y_continuous(limits=c(0,1))+
        scale_color_manual(values = c(diagnosis_time_colors1))+ #+ ggtitle(title) +
        theme_bw()+ 
        #facet_grid(.~title)+
        theme(text = element_text(size=8),
              axis.text=element_text(size=8, face = "bold"),
              axis.title=element_text(size=8,face="bold"),
              legend.position = "none")+ 
        guides(colour = guide_legend(override.aes = list(size=4),nrow=4)) +
        xlab('Time to Diagnosis (Years)') + ylab('CV Concordance Index (95% CI)')
      
      png(paste0(name,'.',c,'.ci.ci.dxtime2.png'),height = 1100,width=1100,res = 300)
      print(plot1)
      dev.off()
      plot1 = ggplot(auc.res2,aes(x = dx.time, y =se.spec.95 , col = dx.time)) +
        geom_errorbar(aes(ymin=se.spec.95.lower, ymax=se.spec.95.upper), width=.2,
                      position=position_dodge(0.05))+
        geom_point()+
        scale_y_continuous(limits=c(0,1))+
        scale_color_manual(values = c(diagnosis_time_colors1))+ #+ ggtitle(title) +
        theme_bw()+ 
        #facet_grid(.~title)+
        theme(text = element_text(size=8),
              axis.text=element_text(size=8, face = "bold"),
              axis.title=element_text(size=8,face="bold"),
              legend.position = "none")+ 
        guides(colour = guide_legend(override.aes = list(size=4),nrow=4)) +
        xlab('Time to Diagnosis (Years)') + ylab('Sensitivity at 95% Specificity (95% CI)')
      
      png(paste0(name,'.',c,'.sens.spec95.ci.dxtime2.png'),height = 1100,width=1100,res = 300)
      print(plot1)
      dev.off()
      
      plot1 = ggplot(auc.res2,aes(x = dx.time, y = se.spec.90, col = dx.time)) +
        geom_errorbar(aes(ymin=se.spec.90.lower, ymax=se.spec.90.upper), width=.2,
                      position=position_dodge(0.05))+
        geom_point()+
        scale_y_continuous(limits=c(0,1))+
        scale_color_manual(values = c(diagnosis_time_colors1))+ #+ ggtitle(title) +
        theme_bw()+ 
        #facet_grid(.~title)+
        theme(text = element_text(size=8),
              axis.text=element_text(size=8, face = "bold"),
              axis.title=element_text(size=8,face="bold"),
              legend.position = "none")+ 
        guides(colour = guide_legend(override.aes = list(size=4),nrow=4)) +
        xlab('Time to Diagnosis (Years)')+ ylab('Sensitivity at 90% Specificity (95% CI)')
      
      png(paste0(name,'.',c,'.sens.spec90.ci.dxtime2.png'),height = 1100,width=1100,res = 300)
      print(plot1)
      dev.off()
      
      
      #youdens index sens
      plot1 = ggplot(auc.res2,aes(x = dx.time, y = jcutpoint.sens, col = dx.time)) +
        geom_errorbar(aes(ymin=jcutpoint.sens.lower, ymax=jcutpoint.sens.upper), width=.2,
                      position=position_dodge(0.05))+
        geom_point()+
        scale_y_continuous(limits=c(0,1))+
        scale_color_manual(values = c(diagnosis_time_colors1))+ #+ ggtitle(title) +
        theme_bw()+ 
        #facet_grid(.~title)+
        theme(text = element_text(size=8),
              axis.text=element_text(size=8, face = "bold"),
              axis.title=element_text(size=8,face="bold"),
              legend.position = "none")+ 
        guides(colour = guide_legend(override.aes = list(size=4),nrow=4)) +
        xlab('Time to Diagnosis (Years)')+ ylab('Youden\'s Index Cutoff Sensitivity')
      
      png(paste0(name,'.',c,'.sens.youden.ci.dxtime2.png'),height = 1100,width=1100,res = 300)
      print(plot1)
      dev.off()
      
      #youdens index spec
      plot1 = ggplot(auc.res2,aes(x = dx.time, y = jcutpoint.spec, col = dx.time)) +
        geom_errorbar(aes(ymin=jcutpoint.spec.lower, ymax=jcutpoint.spec.upper), width=.2,
                      position=position_dodge(0.05))+
        geom_point()+
        scale_y_continuous(limits=c(0,1))+
        scale_color_manual(values = c(diagnosis_time_colors1))+ #+ ggtitle(title) +
        theme_bw()+ 
        #facet_grid(.~title)+
        theme(text = element_text(size=8),
              axis.text=element_text(size=8, face = "bold"),
              axis.title=element_text(size=8,face="bold"),
              legend.position = "none")+ 
        guides(colour = guide_legend(override.aes = list(size=4),nrow=4)) +
        xlab('Time to Diagnosis (Years)')+ ylab('Youden\'s Index Cutoff Specificity')
      
      png(paste0(name,'.',c,'.spec.youden.ci.dxtime2.png'),height = 1100,width=1100,res = 300)
      print(plot1)
      dev.off()
      
      
      #youdens index ppv
      plot1 = ggplot(auc.res2,aes(x = dx.time, y = jcutpoint.ppv, col = dx.time)) +
        geom_errorbar(aes(ymin=jcutpoint.ppv.lower, ymax=jcutpoint.ppv.upper), width=.2,
                      position=position_dodge(0.05))+
        geom_point()+
        scale_y_continuous(limits=c(0,1))+
        scale_color_manual(values = c(diagnosis_time_colors1))+ #+ ggtitle(title) +
        theme_bw()+ 
        #facet_grid(.~title)+
        theme(text = element_text(size=8),
              axis.text=element_text(size=8, face = "bold"),
              axis.title=element_text(size=8,face="bold"),
              legend.position = "none")+ 
        guides(colour = guide_legend(override.aes = list(size=4),nrow=4)) +
        xlab('Time to Diagnosis (Years)')+ ylab('Youden\'s Index Cutoff PPV')
      
      png(paste0(name,'.',c,'.ppv.youden.ci.dxtime2.png'),height = 1100,width=1100,res = 300)
      print(plot1)
      dev.off()
      
      #youdens index npv
      plot1 = ggplot(auc.res2,aes(x = dx.time, y = jcutpoint.npv, col = dx.time)) +
        geom_errorbar(aes(ymin=jcutpoint.npv.lower, ymax=jcutpoint.npv.upper), width=.2,
                      position=position_dodge(0.05))+
        geom_point()+
        scale_y_continuous(limits=c(0,1))+
        scale_color_manual(values = c(diagnosis_time_colors1))+ #+ ggtitle(title) +
        theme_bw()+ 
        #facet_grid(.~title)+
        theme(text = element_text(size=8),
              axis.text=element_text(size=8, face = "bold"),
              axis.title=element_text(size=8,face="bold"),
              legend.position = "none")+ 
        guides(colour = guide_legend(override.aes = list(size=4),nrow=4)) +
        xlab('Time to Diagnosis (Years)')+ ylab('Youden\'s Index Cutoff NPV')
      
      png(paste0(name,'.',c,'.npv.youden.ci.dxtime2.png'),height = 1100,width=1100,res = 300)
      print(plot1)
      dev.off()
      
      
      #f1 sens
      plot1 = ggplot(auc.res2,aes(x = dx.time, y = f1cutpoint.sens, col = dx.time)) +
        geom_errorbar(aes(ymin=f1cutpoint.sens.lower, ymax=f1cutpoint.sens.upper), width=.2,
                      position=position_dodge(0.05))+
        geom_point()+
        scale_y_continuous(limits=c(0,1))+
        scale_color_manual(values = c(diagnosis_time_colors1))+ #+ ggtitle(title) +
        theme_bw()+ 
        #facet_grid(.~title)+
        theme(text = element_text(size=8),
              axis.text=element_text(size=8, face = "bold"),
              axis.title=element_text(size=8,face="bold"),
              legend.position = "none")+ 
        guides(colour = guide_legend(override.aes = list(size=4),nrow=4)) +
        xlab('Time to Diagnosis (Years)')+ ylab('F1-Score Maximized Cutoff Sensitivity')
      
      png(paste0(name,'.',c,'.sens.f1.ci.dxtime2.png'),height = 1100,width=1100,res = 300)
      print(plot1)
      dev.off()
      
      #f1 index spec
      plot1 = ggplot(auc.res2,aes(x = dx.time, y = f1cutpoint.spec, col = dx.time)) +
        geom_errorbar(aes(ymin=f1cutpoint.spec.lower, ymax=f1cutpoint.spec.upper), width=.2,
                      position=position_dodge(0.05))+
        geom_point()+
        scale_y_continuous(limits=c(0,1))+
        scale_color_manual(values = c(diagnosis_time_colors1))+ #+ ggtitle(title) +
        theme_bw()+ 
        #facet_grid(.~title)+
        theme(text = element_text(size=8),
              axis.text=element_text(size=8, face = "bold"),
              axis.title=element_text(size=8,face="bold"),
              legend.position = "none")+ 
        guides(colour = guide_legend(override.aes = list(size=4),nrow=4)) +
        xlab('Time to Diagnosi (Years)s')+ ylab('F1-Score Maximized Cutoff Specificity')
      
      png(paste0(name,'.',c,'.spec.f1.ci.dxtime2.png'),height = 1100,width=1100,res = 300)
      print(plot1)
      dev.off()
      
      
      #youdens index ppv
      plot1 = ggplot(auc.res2,aes(x = dx.time, y = f1cutpoint.ppv, col = dx.time)) +
        geom_errorbar(aes(ymin=f1cutpoint.ppv.lower, ymax=f1cutpoint.ppv.upper), width=.2,
                      position=position_dodge(0.05))+
        geom_point()+
        scale_y_continuous(limits=c(0,1))+
        scale_color_manual(values = c(diagnosis_time_colors1))+ #+ ggtitle(title) +
        theme_bw()+ 
        #facet_grid(.~title)+
        theme(text = element_text(size=8),
              axis.text=element_text(size=8, face = "bold"),
              axis.title=element_text(size=8,face="bold"),
              legend.position = "none")+ 
        guides(colour = guide_legend(override.aes = list(size=4),nrow=4)) +
        xlab('Time to Diagnosis (Years)')+ ylab('F1-Score Maximized Cutoff PPV')
      
      png(paste0(name,'.',c,'.ppv.f1.ci.dxtime2.png'),height = 1100,width=1100,res = 300)
      print(plot1)
      dev.off()
      
      #f1 index npv
      plot1 = ggplot(auc.res2,aes(x = dx.time, y = f1cutpoint.npv, col = dx.time)) +
        geom_errorbar(aes(ymin=f1cutpoint.npv.lower, ymax=f1cutpoint.npv.upper), width=.2,
                      position=position_dodge(0.05))+
        geom_point()+
        scale_y_continuous(limits=c(0,1))+
        scale_color_manual(values = c(diagnosis_time_colors1))+ #+ ggtitle(title) +
        theme_bw()+ 
        #facet_grid(.~title)+
        theme(text = element_text(size=8),
              axis.text=element_text(size=8, face = "bold"),
              axis.title=element_text(size=8,face="bold"),
              legend.position = "none")+ 
        guides(colour = guide_legend(override.aes = list(size=4),nrow=4)) +
        xlab('Time to Diagnosis (Years)')+ ylab('F1-Score Maximized Cutoff NPV')
      
      png(paste0(name,'.',c,'.npv.f1.ci.dxtime2.png'),height = 1100,width=1100,res = 300)
      print(plot1)
      dev.off()
      
      my_comparisons = list(c('Control','0-2'),
                            c('Control','2-4'),
                            c('Control','4-6'),
                            c('Control','6-8'),
                            c('Control','8-10'))
      diagnosis_time_colors1 = c('grey',"#048BA8",'#AAF683','#FFD97D','#FF9B85','#C8553D')
      names(diagnosis_time_colors1) = c('All','0-2','2-4','4-6','6-8','8-10')
      
      options(scipen=2)
      
      diagnosis_time_grouping = function(diagnosis_time) {
        tmp = ifelse(diagnosis_time > 2920, '8-10', diagnosis_time)
        tmp = ifelse(diagnosis_time <= 2920, '6-8', tmp)
        tmp = ifelse(diagnosis_time <= 2190, '4-6', tmp)
        tmp = ifelse(diagnosis_time <= 1460, '2-4', tmp)
        tmp = ifelse(diagnosis_time <= 730, '0-2', tmp)
        tmp[is.na(diagnosis_time)] = 'Control'
        diagnosis_time_groups = c('0-2','2-4','4-6','6-8','8-10','Control')
        tmp = factor(tmp, levels = rev(diagnosis_time_groups))
        return(tmp)
      }
      pred.df.targ.collapse.annotated = merge(pred.df.targ.collapse, merged.df.all.tmp[,c('GRP_Id','diff_in_days')],by='GRP_Id')
      
      pred.df.targ.collapse.annotated$Diagnosis.time.ks = diagnosis_time_grouping(pred.df.targ.collapse.annotated$diff_in_days)
      print(kruskal.test(methylation_score ~ Diagnosis.time.ks,  data = pred.df.targ.collapse.annotated)  )
      
      plot1 = ggplot(pred.df.targ.collapse.annotated,aes(x = Diagnosis.time.ks, y = methylation_score, col = Diagnosis.time.ks)) + geom_boxplot(outlier.shape=NA)+geom_jitter(width=0.1)+
        scale_color_manual(values = c(diagnosis_time_colors1))+ #+ ggtitle(title) +
        theme_bw()+ 
        scale_y_continuous(limits=c(0,1.8),breaks = seq(0,1.85,0.25))+
        theme(text = element_text(size=8),
              axis.text=element_text(size=8, face = "bold"),
              axis.title=element_text(size=8,face="bold"),
              legend.position = "none")+ guides(colour = guide_legend(override.aes = list(size=4),nrow=4)) + xlab('Time to Diagnosis (Years)') + ylab('Methylation Score') +
        stat_compare_means(comparisons = my_comparisons,label.y=seq(1,1.85,0.1),size=3)
      
      
      png(paste0(name,'.dx2.methscore.',c,'.png'),height = 1100,width=1100,res = 300)
      print(plot1)
      dev.off()
      
      
      
    }
    
    
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
    
    #plotting by dx time
    merged.df.all.tmp$Diagnosis_Time = diagnosis_time_grouping(merged.df.all.tmp$diff_in_days)
    dx1 = merged.df.all.tmp[merged.df.all.tmp$Diagnosis_Time %in% c('Control','0-1'),]
    dx2 = merged.df.all.tmp[merged.df.all.tmp$Diagnosis_Time %in% c('Control','1-2'),]
    dx3 = merged.df.all.tmp[merged.df.all.tmp$Diagnosis_Time %in% c('Control','2-3'),]
    dx4 = merged.df.all.tmp[merged.df.all.tmp$Diagnosis_Time %in% c('Control','3-4'),]
    dx5 = merged.df.all.tmp[merged.df.all.tmp$Diagnosis_Time %in% c('Control','4-5'),]
    
    

    #dx time 1 year
    
    pred.df.targ.collapse$reported = factor(pred.df.targ.collapse$reported,levels = c('Control','Cancer'))
    auc.plot.all = summary.calc.dxtime(pred.df.targ.collapse,merged.df.all.tmp)
    auc.plot.all[[2]]$dx.time = 'All'
    auc.plot.all[[2]]$var = 'All'
    
    auc.plot.all[[2]]$var.group = 'Time to Diagnosis'
    
    auc.plot.all.dx1 = summary.calc.dxtime(pred.df.targ.collapse[pred.df.targ.collapse$GRP_Id %in% dx1$GRP_Id,],merged.df.all.tmp)
    auc.plot.all.dx1[[2]]$dx.time = '0-1'
    auc.plot.all.dx1[[2]]$var = '0-1'
    
    auc.plot.all.dx1[[2]]$var.group = 'Time to Diagnosis'
    
    
    auc.plot.all.dx2 = summary.calc.dxtime(pred.df.targ.collapse[pred.df.targ.collapse$GRP_Id %in% dx2$GRP_Id,],merged.df.all.tmp)
    auc.plot.all.dx2[[2]]$dx.time = '1-2'
    auc.plot.all.dx2[[2]]$var = '1-2'
    
    auc.plot.all.dx2[[2]]$var.group = 'Time to Diagnosis'
    
    auc.plot.all.dx3 = summary.calc.dxtime(pred.df.targ.collapse[pred.df.targ.collapse$GRP_Id %in% dx3$GRP_Id,],merged.df.all.tmp)
    auc.plot.all.dx3[[2]]$dx.time = '2-3'
    auc.plot.all.dx3[[2]]$var = '2-3'
    
    auc.plot.all.dx3[[2]]$var.group = 'Time to Diagnosis'
    
    auc.plot.all.dx4 = summary.calc.dxtime(pred.df.targ.collapse[pred.df.targ.collapse$GRP_Id %in% dx4$GRP_Id,],merged.df.all.tmp)
    auc.plot.all.dx4[[2]]$dx.time = '3-4'
    auc.plot.all.dx4[[2]]$var = '3-4'
    auc.plot.all.dx4[[2]]$var.group = 'Time to Diagnosis'
    
    auc.plot.all.dx5 = summary.calc.dxtime(pred.df.targ.collapse[pred.df.targ.collapse$GRP_Id %in% dx5$GRP_Id,],merged.df.all.tmp)
    auc.plot.all.dx5[[2]]$dx.time = '4-5'
    auc.plot.all.dx5[[2]]$var = '4-5'
    
    auc.plot.all.dx5[[2]]$var.group = 'Time to Diagnosis'
    
    if (dx.all == T){
      dx6 = merged.df.all.tmp[merged.df.all.tmp$Diagnosis_Time %in% c('Control','5-6'),]
      dx7 = merged.df.all.tmp[merged.df.all.tmp$Diagnosis_Time %in% c('Control','6-7'),]
      dx8 = merged.df.all.tmp[merged.df.all.tmp$Diagnosis_Time %in% c('Control','7-8'),]
      dx9 = merged.df.all.tmp[merged.df.all.tmp$Diagnosis_Time %in% c('Control','8-9','9+'),]
      
      auc.plot.all.dx6 = summary.calc.dxtime(pred.df.targ.collapse[pred.df.targ.collapse$GRP_Id %in% dx6$GRP_Id,],merged.df.all.tmp)
      auc.plot.all.dx6[[2]]$dx.time = '5-6'
      auc.plot.all.dx6[[2]]$var = '5-6'
      
      auc.plot.all.dx6[[2]]$var.group = 'Time to Diagnosis'
      
      auc.plot.all.dx7 = summary.calc.dxtime(pred.df.targ.collapse[pred.df.targ.collapse$GRP_Id %in% dx7$GRP_Id,],merged.df.all.tmp)
      auc.plot.all.dx7[[2]]$dx.time = '6-7'
      auc.plot.all.dx7[[2]]$var = '6-7'
      
      auc.plot.all.dx7[[2]]$var.group = 'Time to Diagnosis'
      
      
      auc.plot.all.dx8 = summary.calc.dxtime(pred.df.targ.collapse[pred.df.targ.collapse$GRP_Id %in% dx8$GRP_Id,],merged.df.all.tmp)
      auc.plot.all.dx8[[2]]$dx.time = '7-8'
      auc.plot.all.dx8[[2]]$var = '7-8'
      
      auc.plot.all.dx8[[2]]$var.group = 'Time to Diagnosis'
      
      
      auc.plot.all.dx9 = summary.calc.dxtime(pred.df.targ.collapse[pred.df.targ.collapse$GRP_Id %in% dx9$GRP_Id,],merged.df.all.tmp)
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
      
      diagnosis_time_colors1 = c('grey',"#048BA8",'#60D394','#AAF683','#FFD97D','#FF9B85','#C8553D','#F46197','#C3C3E6','#442B48')
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
      
      png(paste0(name,'.',c,'.auroc.dxtime1.png'),height = 1100,width=1100,res = 300)
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
      
      
      #annotating others
      diagnosis_time_colors1 = c('grey',"#048BA8",'#60D394','#AAF683','#FFD97D','#FF9B85','#C8553D')
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
      
      png(paste0(name,'.',c,'.auroc.dxtime1.png'),height = 1100,width=1100,res = 300)
      print(plot1)
      dev.off()
    }
    
    if (dx.all == T){
      merged.df.all.tmp
      dx6 = merged.df.all.tmp[merged.df.all.tmp$Diagnosis_Time %in% c('Control','5-6','6-7','7-8','8-9','9+'),]
      
      
      auc.plot.all.dx6 = summary.calc.dxtime(pred.df.targ.collapse[pred.df.targ.collapse$GRP_Id %in% dx6$GRP_Id,],merged.df.all.tmp)
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
      
      png(paste0(name,'.',c,'.auroc.dxtime1.png'),height = 1100,width=1100,res = 300)
      print(plot1)
      dev.off()
      
      
    } 
    
    
    
    #plotting auc with error bars
    auc.res2$title = 'Diagnosis Time'
    
    
    plot1 = ggplot(auc.res2,aes(x = dx.time, y = auc, col = dx.time)) +
      geom_errorbar(aes(ymin=auc.lower, ymax=auc.upper), width=.2,
                    position=position_dodge(0.05))+
      geom_point()+
      scale_y_continuous(limits=c(0,1))+
      scale_color_manual(values = c(diagnosis_time_colors1))+ #+ ggtitle(title) +
      theme_bw()+ 
      #facet_grid(.~title)+
      
      theme(text = element_text(size=8),
            axis.text=element_text(size=8, face = "bold"),
            axis.title=element_text(size=8,face="bold"),
            legend.position = "none")+ 
      guides(colour = guide_legend(override.aes = list(size=4),nrow=4)) +
      xlab('Sampling time prior to diagnosis') + ylab('CV AUROC (95% CI)')
    
    png(paste0(name,'.',c,'.auroc.ci.dxtime1.png'),height = 1100,width=1100,res = 300)
    print(plot1)
    dev.off()
    plot1 = ggplot(auc.res2,aes(x = dx.time, y = ci, col = dx.time)) +
      geom_errorbar(aes(ymin=ci.lower, ymax=ci.upper), width=.2,
                    position=position_dodge(0.05))+
      geom_point()+
      scale_y_continuous(limits=c(0,1))+
      scale_color_manual(values = c(diagnosis_time_colors1))+ #+ ggtitle(title) +
      theme_bw()+ 
      #facet_grid(.~title)+
      theme(text = element_text(size=8),
            axis.text=element_text(size=8, face = "bold"),
            axis.title=element_text(size=8,face="bold"),
            legend.position = "none")+ 
      guides(colour = guide_legend(override.aes = list(size=4),nrow=4)) +
      xlab('Time to Diagnosis (Years)') + ylab('CV Concordance Index (95% CI)')
    
    png(paste0(name,'.',c,'.ci.ci.dxtime1.png'),height = 1100,width=1100,res = 300)
    print(plot1)
    dev.off()
    plot1 = ggplot(auc.res2,aes(x = dx.time, y =se.spec.95 , col = dx.time)) +
      geom_errorbar(aes(ymin=se.spec.95.lower, ymax=se.spec.95.upper), width=.2,
                    position=position_dodge(0.05))+
      geom_point()+
      scale_y_continuous(limits=c(0,1))+
      scale_color_manual(values = c(diagnosis_time_colors1))+ #+ ggtitle(title) +
      theme_bw()+ 
      #facet_grid(.~title)+
      theme(text = element_text(size=8),
            axis.text=element_text(size=8, face = "bold"),
            axis.title=element_text(size=8,face="bold"),
            legend.position = "none")+ 
      guides(colour = guide_legend(override.aes = list(size=4),nrow=4)) +
      xlab('Time to Diagnosis (Years)') + ylab('Sensitivity at 95% Specificity (95% CI)')
    
    png(paste0(name,'.',c,'.sens.spec95.ci.dxtime1.png'),height = 1100,width=1100,res = 300)
    print(plot1)
    dev.off()
    
    plot1 = ggplot(auc.res2,aes(x = dx.time, y = se.spec.90, col = dx.time)) +
      geom_errorbar(aes(ymin=se.spec.90.lower, ymax=se.spec.90.upper), width=.2,
                    position=position_dodge(0.05))+
      geom_point()+
      scale_y_continuous(limits=c(0,1))+
      scale_color_manual(values = c(diagnosis_time_colors1))+ #+ ggtitle(title) +
      theme_bw()+ 
      #facet_grid(.~title)+
      theme(text = element_text(size=8),
            axis.text=element_text(size=8, face = "bold"),
            axis.title=element_text(size=8,face="bold"),
            legend.position = "none")+ 
      guides(colour = guide_legend(override.aes = list(size=4),nrow=4)) +
      xlab('Time to Diagnosis (Years)')+ ylab('Sensitivity at 90% Specificity (95% CI)')
    
    png(paste0(name,'.',c,'.sens.spec90.ci.dxtime1.png'),height = 1100,width=1100,res = 300)
    print(plot1)
    dev.off()
    
    
    #youdens index sens
    plot1 = ggplot(auc.res2,aes(x = dx.time, y = jcutpoint.sens, col = dx.time)) +
      geom_errorbar(aes(ymin=jcutpoint.sens.lower, ymax=jcutpoint.sens.upper), width=.2,
                    position=position_dodge(0.05))+
      geom_point()+
      scale_y_continuous(limits=c(0,1))+
      scale_color_manual(values = c(diagnosis_time_colors1))+ #+ ggtitle(title) +
      theme_bw()+ 
      #facet_grid(.~title)+
      theme(text = element_text(size=8),
            axis.text=element_text(size=8, face = "bold"),
            axis.title=element_text(size=8,face="bold"),
            legend.position = "none")+ 
      guides(colour = guide_legend(override.aes = list(size=4),nrow=4)) +
      xlab('Time to Diagnosis (Years)')+ ylab('Youden\'s Index Cutoff Sensitivity')
    
    png(paste0(name,'.',c,'.sens.youden.ci.dxtime1.png'),height = 1100,width=1100,res = 300)
    print(plot1)
    dev.off()
    
    #youdens index spec
    plot1 = ggplot(auc.res2,aes(x = dx.time, y = jcutpoint.spec, col = dx.time)) +
      geom_errorbar(aes(ymin=jcutpoint.spec.lower, ymax=jcutpoint.spec.upper), width=.2,
                    position=position_dodge(0.05))+
      geom_point()+
      scale_y_continuous(limits=c(0,1))+
      scale_color_manual(values = c(diagnosis_time_colors1))+ #+ ggtitle(title) +
      theme_bw()+ 
      #facet_grid(.~title)+
      theme(text = element_text(size=8),
            axis.text=element_text(size=8, face = "bold"),
            axis.title=element_text(size=8,face="bold"),
            legend.position = "none")+ 
      guides(colour = guide_legend(override.aes = list(size=4),nrow=4)) +
      xlab('Time to Diagnosis (Years)')+ ylab('Youden\'s Index Cutoff Specificity')
    
    png(paste0(name,'.',c,'.spec.youden.ci.dxtime1.png'),height = 1100,width=1100,res = 300)
    print(plot1)
    dev.off()
    
    
    #youdens index ppv
    plot1 = ggplot(auc.res2,aes(x = dx.time, y = jcutpoint.ppv, col = dx.time)) +
      geom_errorbar(aes(ymin=jcutpoint.ppv.lower, ymax=jcutpoint.ppv.upper), width=.2,
                    position=position_dodge(0.05))+
      geom_point()+
      scale_y_continuous(limits=c(0,1))+
      scale_color_manual(values = c(diagnosis_time_colors1))+ #+ ggtitle(title) +
      theme_bw()+ 
      #facet_grid(.~title)+
      theme(text = element_text(size=8),
            axis.text=element_text(size=8, face = "bold"),
            axis.title=element_text(size=8,face="bold"),
            legend.position = "none")+ 
      guides(colour = guide_legend(override.aes = list(size=4),nrow=4)) +
      xlab('Time to Diagnosis (Years)')+ ylab('Youden\'s Index Cutoff PPV')
    
    png(paste0(name,'.',c,'.ppv.youden.ci.dxtime1.png'),height = 1100,width=1100,res = 300)
    print(plot1)
    dev.off()
    
    #youdens index npv
    plot1 = ggplot(auc.res2,aes(x = dx.time, y = jcutpoint.npv, col = dx.time)) +
      geom_errorbar(aes(ymin=jcutpoint.npv.lower, ymax=jcutpoint.npv.upper), width=.2,
                    position=position_dodge(0.05))+
      geom_point()+
      scale_y_continuous(limits=c(0,1))+
      scale_color_manual(values = c(diagnosis_time_colors1))+ #+ ggtitle(title) +
      theme_bw()+ 
      #facet_grid(.~title)+
      theme(text = element_text(size=8),
            axis.text=element_text(size=8, face = "bold"),
            axis.title=element_text(size=8,face="bold"),
            legend.position = "none")+ 
      guides(colour = guide_legend(override.aes = list(size=4),nrow=4)) +
      xlab('Time to Diagnosis (Years)')+ ylab('Youden\'s Index Cutoff NPV')
    
    png(paste0(name,'.',c,'.npv.youden.ci.dxtime1.png'),height = 1100,width=1100,res = 300)
    print(plot1)
    dev.off()
    
    
    #f1 sens
    plot1 = ggplot(auc.res2,aes(x = dx.time, y = f1cutpoint.sens, col = dx.time)) +
      geom_errorbar(aes(ymin=f1cutpoint.sens.lower, ymax=f1cutpoint.sens.upper), width=.2,
                    position=position_dodge(0.05))+
      geom_point()+
      scale_y_continuous(limits=c(0,1))+
      scale_color_manual(values = c(diagnosis_time_colors1))+ #+ ggtitle(title) +
      theme_bw()+ 
      #facet_grid(.~title)+
      theme(text = element_text(size=8),
            axis.text=element_text(size=8, face = "bold"),
            axis.title=element_text(size=8,face="bold"),
            legend.position = "none")+ 
      guides(colour = guide_legend(override.aes = list(size=4),nrow=4)) +
      xlab('Time to Diagnosis (Years)')+ ylab('F1-Score Maximized Cutoff Sensitivity')
    
    png(paste0(name,'.',c,'.sens.f1.ci.dxtime1.png'),height = 1100,width=1100,res = 300)
    print(plot1)
    dev.off()
    
    #f1 index spec
    plot1 = ggplot(auc.res2,aes(x = dx.time, y = f1cutpoint.spec, col = dx.time)) +
      geom_errorbar(aes(ymin=f1cutpoint.spec.lower, ymax=f1cutpoint.spec.upper), width=.2,
                    position=position_dodge(0.05))+
      geom_point()+
      scale_y_continuous(limits=c(0,1))+
      scale_color_manual(values = c(diagnosis_time_colors1))+ #+ ggtitle(title) +
      theme_bw()+ 
      #facet_grid(.~title)+
      theme(text = element_text(size=8),
            axis.text=element_text(size=8, face = "bold"),
            axis.title=element_text(size=8,face="bold"),
            legend.position = "none")+ 
      guides(colour = guide_legend(override.aes = list(size=4),nrow=4)) +
      xlab('Time to Diagnosis (Years)')+ ylab('F1-Score Maximized Cutoff Specificity')
    
    png(paste0(name,'.',c,'.spec.f1.ci.dxtime1.png'),height = 1100,width=1100,res = 300)
    print(plot1)
    dev.off()
    
    
    #youdens index ppv
    plot1 = ggplot(auc.res2,aes(x = dx.time, y = f1cutpoint.ppv, col = dx.time)) +
      geom_errorbar(aes(ymin=f1cutpoint.ppv.lower, ymax=f1cutpoint.ppv.upper), width=.2,
                    position=position_dodge(0.05))+
      geom_point()+
      scale_y_continuous(limits=c(0,1))+
      scale_color_manual(values = c(diagnosis_time_colors1))+ #+ ggtitle(title) +
      theme_bw()+ 
      #facet_grid(.~title)+
      theme(text = element_text(size=8),
            axis.text=element_text(size=8, face = "bold"),
            axis.title=element_text(size=8,face="bold"),
            legend.position = "none")+ 
      guides(colour = guide_legend(override.aes = list(size=4),nrow=4)) +
      xlab('Time to Diagnosis (Years)')+ ylab('F1-Score Maximized Cutoff PPV')
    
    png(paste0(name,'.',c,'.ppv.f1.ci.dxtime1.png'),height = 1100,width=1100,res = 300)
    print(plot1)
    dev.off()
    
    #f1 index npv
    plot1 = ggplot(auc.res2,aes(x = dx.time, y = f1cutpoint.npv, col = dx.time)) +
      geom_errorbar(aes(ymin=f1cutpoint.npv.lower, ymax=f1cutpoint.npv.upper), width=.2,
                    position=position_dodge(0.05))+
      geom_point()+
      scale_y_continuous(limits=c(0,1))+
      scale_color_manual(values = c(diagnosis_time_colors1))+ #+ ggtitle(title) +
      theme_bw()+ 
      #facet_grid(.~title)+
      theme(text = element_text(size=8),
            axis.text=element_text(size=8, face = "bold"),
            axis.title=element_text(size=8,face="bold"),
            legend.position = "none")+ 
      guides(colour = guide_legend(override.aes = list(size=4),nrow=4)) +
      xlab('Time to Diagnosis (Years)')+ ylab('F1-Score Maximized Cutoff NPV')
    
    png(paste0(name,'.',c,'.npv.f1.ci.dxtime1.png'),height = 1100,width=1100,res = 300)
    print(plot1)
    dev.off()
    
    
    
    
    #dx1
    labels = combined.auc.plot.all[combined.auc.plot.all$l ==0,]
    labels$title = gsub('.new','',paste0(labels$dx.time, ' AUROC: ',round(labels$auc,digits=3)))
    colors.filt = diagnosis_time_colors1[names(diagnosis_time_colors1) %in% labels$dx.time]
    labels = labels[order(match(labels$dx.time,names(colors.filt) )),]  
    labels$dx.time = factor(labels$dx.time,levels = names(colors.filt))
    combined.auc.plot.all$dx.time = factor(as.character(combined.auc.plot.all$dx.time), levels= names(colors.filt))
    
    
    plot1 = ggplot(combined.auc.plot.all,aes(x = fpr, y = tpr, col = dx.time)) + geom_line(linewidth=1) +
      scale_color_manual(values = c(diagnosis_time_colors1),labels = labels$title)+ #+ ggtitle(title) +
      theme_bw()+ 
      theme(text = element_text(size=8),
            axis.text=element_text(size=8, face = "bold"),
            axis.title=element_text(size=8,face="bold"),
            legend.position = "right")+ guides(colour = guide_legend(override.aes = list(size=4),ncol=1)) + xlab('False Positive Rate') + ylab('Sensitivity')
    
    png(paste0(name,'.',c,'.auroc.dxtime1.labs.png'),height = 1100,width=1100,res = 300)
    print(plot1)
    dev.off()
    
    pred.df.targ.collapse.annotated = merge(pred.df.targ.collapse, merged.df.all.tmp[,c('GRP_Id','Diagnosis_Time')],by='GRP_Id')
    
    #boxplot of score
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
      diagnosis_time_colors1 = c('grey',"#048BA8",'#60D394','#AAF683','#FFD97D','#FF9B85','#C8553D','#F46197','#C3C3E6','#442B48')
      names(diagnosis_time_colors1) = c('Control','0-1','1-2','2-3','3-4','4-5','5-6','6-7','7-8','8-10')
      
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
    
    auc.t =rbind( #auc.plot.all[[3]],
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
    
    png(paste0(name,'.',c,'.timeauroc.stage.png'),height = 1100,width=1100,res = 300)
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
    
    png(paste0(name,'.',c,'.auroc.stage.png'),height = 1100,width=1100,res = 300)
    print(plot1)
    dev.off()
    
    
    #stage labels
    labels = combined.auc.plot.all[combined.auc.plot.all$l ==0,]
    labels$title = gsub('.new','',paste0(labels$stage, ' AUROC: ',round(labels$auc,digits=3)))
    colors.filt = stage_colors[names(stage_colors) %in% labels$stage]
    labels = labels[order(match(labels$stage,names(colors.filt) )),]  
    labels$stage = factor(labels$stage,levels = names(colors.filt))
    combined.auc.plot.all$stage = factor(as.character(combined.auc.plot.all$stage), levels= names(colors.filt))
    
    
    plot1 = ggplot(combined.auc.plot.all,aes(x = fpr, y = tpr, col = stage)) + geom_line(linewidth=1) +
      scale_color_manual(values = c(stage_colors),labels = labels$title)+ #+ ggtitle(title) +
      theme_bw()+ 
      theme(text = element_text(size=8),
            axis.text=element_text(size=8, face = "bold"),
            axis.title=element_text(size=8,face="bold"),
            legend.position = "right")+ guides(colour = guide_legend(override.aes = list(size=4),nrow=4)) + xlab('False Positive Rate') + ylab('Sensitivity')
    
    png(paste0(name,'.',c,'.auroc.stage.labs.png'),height = 1000,width=1000,res = 300)
    print(plot1)
    dev.off()
    
    stage_colors.sample = c('Not Reported' = '#4ECDC4','Control'= '#7A797C','0' = "#000004FF", 'I' = "#3B0F70FF",'II' = "#8C2981FF",'III' = "#DE4968FF",'IV' = "#FE9F6DFF")
    
    
    
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
    
    png(paste0(name,'.',c,'.auroc.ci.stage.png'),height = 1100,width=600,res = 300)
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
    
    png(paste0(name,'.',c,'.ci.ci.stage.png'),height = 1100,width=600,res = 300)
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
    
    png(paste0(name,'.',c,'.sens.spec95.ci.stage.png'),height = 1100,width=600,res = 300)
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
    
    png(paste0(name,'.',c,'.sens.spec90.ci.stage.png'),height = 1100,width=600,res = 300)
    print(plot1)
    dev.off()
    
    
    
    pred.df.targ.collapse.annotated = merge(pred.df.targ.collapse, merged.df.all.tmp[,c('GRP_Id','Diagnosis_Time','filler','Stage','GRADE_CD')],by='GRP_Id')
    stage_colors.score = c('Control' = "grey", 'I' = "#3B0F70FF",'II' = "#8C2981FF",'III' = "#DE4968FF",'IV' = '#FE9F6DFF','NR' = '#291F1E')
    
    #boxplot of score
    pred.df.targ.collapse.annotated$Stage = ifelse(pred.df.targ.collapse.annotated$Stage == 'Not Reported/Unknown','NR',pred.df.targ.collapse.annotated$Stage)
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
    
    
    
    
    #youdens index sens
    plot1 = ggplot(auc.res1,aes(x = var, y = jcutpoint.sens, col = var)) +
      geom_errorbar(aes(ymin=jcutpoint.sens.lower, ymax=jcutpoint.sens.upper), width=.2,
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
      xlab('Stage')+ ylab('Youden\'s Index Cutoff Sensitivity')
    
    png(paste0(name,'.',c,'.sens.youden.ci.stage.png'),height = 1100,width=600,res = 300)
    print(plot1)
    dev.off()
    
    #youdens index spec
    plot1 = ggplot(auc.res1,aes(x = var, y = jcutpoint.spec, col = var)) +
      geom_errorbar(aes(ymin=jcutpoint.spec.lower, ymax=jcutpoint.spec.upper), width=.2,
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
      xlab('Stage')+ ylab('Youden\'s Index Cutoff Specificity')
    
    png(paste0(name,'.',c,'.spec.youden.ci.stage.png'),height = 1100,width=600,res = 300)
    print(plot1)
    dev.off()
    
    
    #youdens index ppv
    plot1 = ggplot(auc.res1,aes(x = var, y = jcutpoint.ppv, col = var)) +
      geom_errorbar(aes(ymin=jcutpoint.ppv.lower, ymax=jcutpoint.ppv.upper), width=.2,
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
      xlab('Stage')+ ylab('Youden\'s Index Cutoff PPV')
    
    png(paste0(name,'.',c,'.ppv.youden.ci.stage.png'),height = 1100,width=600,res = 300)
    print(plot1)
    dev.off()
    
    #youdens index npv
    plot1 = ggplot(auc.res1,aes(x = var, y = jcutpoint.npv, col = var)) +
      geom_errorbar(aes(ymin=jcutpoint.npv.lower, ymax=jcutpoint.npv.upper), width=.2,
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
      xlab('Stage')+ ylab('Youden\'s Index Cutoff NPV')
    
    png(paste0(name,'.',c,'.npv.youden.ci.stage.png'),height = 1100,width=600,res = 300)
    print(plot1)
    dev.off()
    
    
    #f1 sens
    plot1 = ggplot(auc.res1,aes(x = var, y = f1cutpoint.sens, col = var)) +
      geom_errorbar(aes(ymin=f1cutpoint.sens.lower, ymax=f1cutpoint.sens.upper), width=.2,
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
      xlab('Stage')+ ylab('F1-Score Maximized Cutoff Sensitivity')
    
    png(paste0(name,'.',c,'.sens.f1.ci.stage.png'),height = 1100,width=600,res = 300)
    print(plot1)
    dev.off()
    
    #f1 index spec
    plot1 = ggplot(auc.res1,aes(x = var, y = f1cutpoint.spec, col = var)) +
      geom_errorbar(aes(ymin=f1cutpoint.spec.lower, ymax=f1cutpoint.spec.upper), width=.2,
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
      xlab('Stage')+ ylab('F1-Score Maximized Cutoff Specificity')
    
    png(paste0(name,'.',c,'.spec.f1.ci.stage.png'),height = 1100,width=600,res = 300)
    print(plot1)
    dev.off()
    
    
    #youdens index ppv
    plot1 = ggplot(auc.res1,aes(x = var, y = f1cutpoint.ppv, col = var)) +
      geom_errorbar(aes(ymin=f1cutpoint.ppv.lower, ymax=f1cutpoint.ppv.upper), width=.2,
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
      xlab('Stage')+ ylab('F1-Score Maximized Cutoff PPV')
    
    png(paste0(name,'.',c,'.ppv.f1.ci.stage.png'),height = 1100,width=600,res = 300)
    print(plot1)
    dev.off()
    
    #f1 index npv
    plot1 = ggplot(auc.res1,aes(x = var, y = f1cutpoint.npv, col = var)) +
      geom_errorbar(aes(ymin=f1cutpoint.npv.lower, ymax=f1cutpoint.npv.upper), width=.2,
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
      xlab('Stage')+ ylab('F1-Score Maximized Cutoff NPV')
    
    png(paste0(name,'.',c,'.npv.f1.ci.stage.png'),height = 1100,width=600,res = 300)
    print(plot1)
    dev.off()
    

    ####plotting grade#####
    print('grade')
    if (gs.plot == T){
      grade1 = merged.df.all.tmp[merged.df.all.tmp$Gleason.score %in% ('6') | merged.df.all.tmp$group =='Control',]
      grade2 = merged.df.all.tmp[merged.df.all.tmp$Gleason.score %in% c('7')| merged.df.all.tmp$group =='Control',]
      grade3 = merged.df.all.tmp[merged.df.all.tmp$Gleason.score %in% c('8','9')| merged.df.all.tmp$group =='Control',]
      
      
      
      auc.plot.all.g0 = summary.calc(pred.df.targ.collapse,merged.df.all.tmp )
      auc.plot.all.g0[[2]]$grade='All'
      auc.plot.all.g0[[2]]$var.group='Grade'
      auc.plot.all.g0[[2]]$var='All'
      auc.plot.all.g0[[3]]$var='All'
      
      
      auc.plot.all.g1 = summary.calc(pred.df.targ.collapse[pred.df.targ.collapse$GRP_Id %in% grade1$GRP_Id,],merged.df.all.tmp )
      auc.plot.all.g1[[2]]$grade='6'
      auc.plot.all.g1[[2]]$var.group='Grade'
      auc.plot.all.g1[[2]]$var='6'
      auc.plot.all.g1[[3]]$var='6'
      
      auc.plot.all.g2 = summary.calc(pred.df.targ.collapse[pred.df.targ.collapse$GRP_Id %in% grade2$GRP_Id,],merged.df.all.tmp )
      auc.plot.all.g2[[2]]$grade='7'
      auc.plot.all.g2[[2]]$var.group='Grade'
      auc.plot.all.g2[[2]]$var='7'
      auc.plot.all.g2[[3]]$var='7'
      
      auc.plot.all.g3 = summary.calc(pred.df.targ.collapse[pred.df.targ.collapse$GRP_Id %in% grade3$GRP_Id,],merged.df.all.tmp )
      auc.plot.all.g3[[2]]$grade='8/9'
      auc.plot.all.g3[[2]]$var.group='Grade'
      auc.plot.all.g3[[2]]$var='8/9'
      auc.plot.all.g3[[3]]$var='8/9'
      
      #
      
      auc.t =rbind( auc.plot.all.g0[[3]],
                    auc.plot.all.g1[[3]],
                    auc.plot.all.g2[[3]],
                    auc.plot.all.g3[[3]]
      )
      auc.t$title = 'Grade'
      auc.t.combined = rbind(auc.t.combined,auc.t)
      
      auc.res1 = rbind(
        auc.plot.all.g1[[2]][1,],
        auc.plot.all.g2[[2]][1,],
        auc.plot.all.g3[[2]][1,])
      auc.res1$title = 'Grade'
      
      combined.auroc = combined.auroc[,colnames(combined.auroc) %in% colnames(auc.res1)]
      
      combined.auroc = rbind(combined.auroc,auc.res1[,colnames(auc.res1) %in% colnames(combined.auroc)])
      
      
      
      #combined.auroc = combined.auroc[,colnames(combined.auroc) %in% colnames(auc.res1)]
      # auc.res1 = auc.res1[,colnames(auc.res1) %in% colnames(combined.auroc)]
      stage_colors = c('All' = "#000004FF", 'I' = "#3B0F70FF",'II' = "#8C2981FF",'III/IV' = "#DE4968FF")
      
      grade_colors = c('All' = "#000004FF",'6' = '#82A6B1','7' = '#76877D', '8/9' = '#5B3000')
      
      
      
      plot1 = ggplot(auc.t[auc.t$time >= 365*1,],aes(x = time/365, y = AUC, col = var)) + geom_line(linewidth=1) +
        scale_color_manual(values = c(grade_colors))+ #+ ggtitle(title) +
        scale_y_continuous(limits=c(0,1))+
        theme_bw()+ 
        #facet_grid2(.~title)+
        theme(text = element_text(size=8),
              axis.text=element_text(size=8, face = "bold"),
              axis.title=element_text(size=8,face="bold"),
              legend.position = "none")+ guides(colour = guide_legend(override.aes = list(size=4),nrow=4)) +
        xlab('Time (Years)') + ylab('Time-Dependent CV AUROC')
      
      png(paste0(name,'.',c,'.timeauroc.grade.png'),height = 1100,width=1100,res = 300)
      print(plot1)
      dev.off()
      
      
      combined.auc.plot.all = rbind(auc.plot.all.g1[[2]],auc.plot.all.g2[[2]],auc.plot.all.g3[[2]])
      plot1 = ggplot(combined.auc.plot.all,aes(x = fpr, y = tpr, col = grade)) + geom_line(linewidth=1) +
        scale_color_manual(values = c(grade_colors))+ #+ ggtitle(title) +
        theme_bw()+ 
        theme(text = element_text(size=8),
              axis.text=element_text(size=8, face = "bold"),
              axis.title=element_text(size=8,face="bold"),
              legend.position = "none")+ guides(colour = guide_legend(override.aes = list(size=4),nrow=4)) + xlab('False Positive Rate') + ylab('Sensitivity')
      
      png(paste0(name,'.',c,'.auroc.grade.png'),height = 1000,width=1000,res = 300)
      print(plot1)
      dev.off()
      
      #grade labels
      labels = combined.auc.plot.all[combined.auc.plot.all$l ==0,]
      labels$title = gsub('.new','',paste0(labels$grade, ' AUROC: ',round(labels$auc,digits=3)))
      colors.filt = grade_colors[names(grade_colors) %in% labels$grade]
      labels = labels[order(match(labels$grade,names(colors.filt) )),]  
      labels$grade = factor(labels$grade,levels = names(colors.filt))
      combined.auc.plot.all$grade = factor(as.character(combined.auc.plot.all$grade), levels= names(colors.filt))
      
      plot1 = ggplot(combined.auc.plot.all,aes(x = fpr, y = tpr, col = grade)) + geom_line(linewidth=1) +
        scale_color_manual(values = c(grade_colors),labels = labels$title)+ #+ ggtitle(title) +
        theme_bw()+ 
        theme(text = element_text(size=8),
              axis.text=element_text(size=8, face = "bold"),
              axis.title=element_text(size=8,face="bold"),
              legend.position = "right")+ guides(colour = guide_legend(override.aes = list(size=4),nrow=4)) + xlab('False Positive Rate') + ylab('Sensitivity')
      
      png(paste0(name,'.',c,'.auroc.grade.labs.png'),height = 1000,width=1000,res = 300)
      print(plot1)
      dev.off()
      grade_colors.sample = c('Control'= 'grey','6' = '#82A6B1','7' = '#35605A', '8' = '#2A324B','9'='#66462C')
      
      pred.df.targ.collapse$Gleason.score =factor(as.character(pred.df.targ.collapse$Gleason.score), levels = names(grade_colors.sample))
      
      
      my_comparisons= list()
      if (wilcox.test(methylation_score ~ Gleason.score, pred.df.targ.collapse.annotated[pred.df.targ.collapse.annotated$Gleason.score %in% c('Control','6'),])$p.value < 0.05 ) {
        my_comparisons[[length(my_comparisons)+1]] = c('Control','6')
      }
      
      if (wilcox.test(methylation_score ~ Gleason.score, pred.df.targ.collapse.annotated[pred.df.targ.collapse.annotated$Gleason.score %in% c('Control','7'),])$p.value < 0.05 ) {
        my_comparisons[[length(my_comparisons)+1]] = c('Control','7')
      }
      
      if (wilcox.test(methylation_score ~ Gleason.score, pred.df.targ.collapse.annotated[pred.df.targ.collapse.annotated$Gleason.score %in% c('Control','8'),])$p.value < 0.05 ) {
        my_comparisons[[length(my_comparisons)+1]] = c('Control','8')
      }
      
      if (wilcox.test(methylation_score ~ Gleason.score, pred.df.targ.collapse.annotated[pred.df.targ.collapse.annotated$Gleason.score %in% c('Control','9'),])$p.value < 0.05 ) {
        my_comparisons[[length(my_comparisons)+1]] = c('Control','9')
      }
      
      
      
      print(kruskal.test(methylation_score ~ Gleason.score,  data = pred.df.targ.collapse.annotated[as.character(pred.df.targ.collapse.annotated$Gleason.score) %in%  c('Control','6','7','8','9'),])  )
      
      
      plot1 = ggplot(pred.df.targ.collapse,aes(x = Gleason.score, y = methylation_score, col = Gleason.score)) + geom_boxplot(outlier.shape=NA)+geom_jitter(width=0.1)+
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
      
      
      #sens spec
      
      
      plot1 = ggplot(auc.res1,aes(x = var, y = auc, col = var)) +
        geom_errorbar(aes(ymin=auc.lower, ymax=auc.upper), width=.2,
                      position=position_dodge(0.05))+
        geom_point()+
        scale_y_continuous(limits=c(0,1))+
        scale_color_manual(values = c(grade_colors))+ #+ ggtitle(title) +
        theme_bw()+ 
        #facet_grid(.~title)+
        
        theme(text = element_text(size=8),
              axis.text=element_text(size=8, face = "bold"),
              axis.title=element_text(size=8,face="bold"),
              legend.position = "none")+ 
        guides(colour = guide_legend(override.aes = list(size=4),nrow=4)) +
        xlab('Grade at Diagnosis') + ylab('CV AUROC (95% CI)')
      
      png(paste0(name,'.',c,'.auroc.ci.grade.png'),height = 1100,width=600,res = 300)
      print(plot1)
      dev.off()
      plot1 = ggplot(auc.res1,aes(x = var, y = ci, col = var)) +
        geom_errorbar(aes(ymin=ci.lower, ymax=ci.upper), width=.2,
                      position=position_dodge(0.05))+
        geom_point()+
        scale_y_continuous(limits=c(0,1))+
        scale_color_manual(values = c(grade_colors))+ #+ ggtitle(title) +
        theme_bw()+ 
        #facet_grid(.~title)+
        theme(text = element_text(size=8),
              axis.text=element_text(size=8, face = "bold"),
              axis.title=element_text(size=8,face="bold"),
              legend.position = "none")+ 
        guides(colour = guide_legend(override.aes = list(size=4),nrow=4)) +
        xlab('Grade at Diagnosis') + ylab('CV Concordance Index (95% CI)')
      
      png(paste0(name,'.',c,'.ci.ci.grade.png'),height = 1100,width=600,res = 300)
      print(plot1)
      dev.off()
      plot1 = ggplot(auc.res1,aes(x = var, y =se.spec.95 , col = var)) +
        geom_errorbar(aes(ymin=se.spec.95.lower, ymax=se.spec.95.upper), width=.2,
                      position=position_dodge(0.05))+
        geom_point()+
        scale_y_continuous(limits=c(0,1))+
        scale_color_manual(values = c(grade_colors))+ #+ ggtitle(title) +
        theme_bw()+ 
        #facet_grid(.~title)+
        theme(text = element_text(size=8),
              axis.text=element_text(size=8, face = "bold"),
              axis.title=element_text(size=8,face="bold"),
              legend.position = "none")+ 
        guides(colour = guide_legend(override.aes = list(size=4),nrow=4)) +
        xlab('Grade at Diagnosis') + ylab('Sensitivity at 95% Specificity (95% CI)')
      
      png(paste0(name,'.',c,'.sens.spec95.ci.grade.png'),height = 1100,width=600,res = 300)
      print(plot1)
      dev.off()
      
      plot1 = ggplot(auc.res1,aes(x = grade, y = se.spec.90, col = grade)) +
        geom_errorbar(aes(ymin=se.spec.90.lower, ymax=se.spec.90.upper), width=.2,
                      position=position_dodge(0.05))+
        geom_point()+
        scale_y_continuous(limits=c(0,1))+
        scale_color_manual(values = c(grade_colors))+ #+ ggtitle(title) +
        theme_bw()+ 
        #facet_grid(.~title)+
        theme(text = element_text(size=8),
              axis.text=element_text(size=8, face = "bold"),
              axis.title=element_text(size=8,face="bold"),
              legend.position = "none")+ 
        guides(colour = guide_legend(override.aes = list(size=4),nrow=4)) +
        xlab('Grade at Diagnosis')+ ylab('Sensitivity at 90% Specificity (95% CI)')
      
      png(paste0(name,'.',c,'.sens.spec90.ci.grade.png'),height = 1100,width=600,res = 300)
      print(plot1)
      dev.off()
      
      
      #youdens index sens
      plot1 = ggplot(auc.res1,aes(x = var, y = jcutpoint.sens, col = var)) +
        geom_errorbar(aes(ymin=jcutpoint.sens.lower, ymax=jcutpoint.sens.upper), width=.2,
                      position=position_dodge(0.05))+
        geom_point()+
        scale_y_continuous(limits=c(0,1))+
        scale_color_manual(values = c(grade_colors))+ #+ ggtitle(title) +
        theme_bw()+ 
        #facet_grid(.~title)+
        theme(text = element_text(size=8),
              axis.text=element_text(size=8, face = "bold"),
              axis.title=element_text(size=8,face="bold"),
              legend.position = "none")+ 
        guides(colour = guide_legend(override.aes = list(size=4),nrow=4)) +
        xlab('Grade')+ ylab('Youden\'s Index Cutoff Sensitivity')
      
      png(paste0(name,'.',c,'.sens.youden.ci.grade.png'),height = 1100,width=600,res = 300)
      print(plot1)
      dev.off()
      
      #youdens index spec
      plot1 = ggplot(auc.res1,aes(x = var, y = jcutpoint.spec, col = var)) +
        geom_errorbar(aes(ymin=jcutpoint.spec.lower, ymax=jcutpoint.spec.upper), width=.2,
                      position=position_dodge(0.05))+
        geom_point()+
        scale_y_continuous(limits=c(0,1))+
        scale_color_manual(values = c(grade_colors))+ #+ ggtitle(title) +
        theme_bw()+ 
        #facet_grid(.~title)+
        theme(text = element_text(size=8),
              axis.text=element_text(size=8, face = "bold"),
              axis.title=element_text(size=8,face="bold"),
              legend.position = "none")+ 
        guides(colour = guide_legend(override.aes = list(size=4),nrow=4)) +
        xlab('Grade')+ ylab('Youden\'s Index Cutoff Specificity')
      
      png(paste0(name,'.',c,'.spec.youden.ci.grade.png'),height = 1100,width=600,res = 300)
      print(plot1)
      dev.off()
      
      
      #youdens index ppv
      plot1 = ggplot(auc.res1,aes(x = var, y = jcutpoint.ppv, col = var)) +
        geom_errorbar(aes(ymin=jcutpoint.ppv.lower, ymax=jcutpoint.ppv.upper), width=.2,
                      position=position_dodge(0.05))+
        geom_point()+
        scale_y_continuous(limits=c(0,1))+
        scale_color_manual(values = c(grade_colors))+ #+ ggtitle(title) +
        theme_bw()+ 
        #facet_grid(.~title)+
        theme(text = element_text(size=8),
              axis.text=element_text(size=8, face = "bold"),
              axis.title=element_text(size=8,face="bold"),
              legend.position = "none")+ 
        guides(colour = guide_legend(override.aes = list(size=4),nrow=4)) +
        xlab('Grade')+ ylab('Youden\'s Index Cutoff PPV')
      
      png(paste0(name,'.',c,'.ppv.youden.ci.grade.png'),height = 1100,width=600,res = 300)
      print(plot1)
      dev.off()
      
      #youdens index npv
      plot1 = ggplot(auc.res1,aes(x = var, y = jcutpoint.npv, col = var)) +
        geom_errorbar(aes(ymin=jcutpoint.npv.lower, ymax=jcutpoint.npv.upper), width=.2,
                      position=position_dodge(0.05))+
        geom_point()+
        scale_y_continuous(limits=c(0,1))+
        scale_color_manual(values = c(grade_colors))+ #+ ggtitle(title) +
        theme_bw()+ 
        #facet_grid(.~title)+
        theme(text = element_text(size=8),
              axis.text=element_text(size=8, face = "bold"),
              axis.title=element_text(size=8,face="bold"),
              legend.position = "none")+ 
        guides(colour = guide_legend(override.aes = list(size=4),nrow=4)) +
        xlab('Grade')+ ylab('Youden\'s Index Cutoff NPV')
      
      png(paste0(name,'.',c,'.npv.youden.ci.grade.png'),height = 1100,width=600,res = 300)
      print(plot1)
      dev.off()
      
      
      #f1 sens
      plot1 = ggplot(auc.res1,aes(x = var, y = f1cutpoint.sens, col = var)) +
        geom_errorbar(aes(ymin=f1cutpoint.sens.lower, ymax=f1cutpoint.sens.upper), width=.2,
                      position=position_dodge(0.05))+
        geom_point()+
        scale_y_continuous(limits=c(0,1))+
        scale_color_manual(values = c(grade_colors))+ #+ ggtitle(title) +
        theme_bw()+ 
        #facet_grid(.~title)+
        theme(text = element_text(size=8),
              axis.text=element_text(size=8, face = "bold"),
              axis.title=element_text(size=8,face="bold"),
              legend.position = "none")+ 
        guides(colour = guide_legend(override.aes = list(size=4),nrow=4)) +
        xlab('Grade')+ ylab('F1-Score Maximized Cutoff Sensitivity')
      
      png(paste0(name,'.',c,'.sens.f1.ci.grade.png'),height = 1100,width=600,res = 300)
      print(plot1)
      dev.off()
      
      #f1 index spec
      plot1 = ggplot(auc.res1,aes(x = var, y = f1cutpoint.spec, col = var)) +
        geom_errorbar(aes(ymin=f1cutpoint.spec.lower, ymax=f1cutpoint.spec.upper), width=.2,
                      position=position_dodge(0.05))+
        geom_point()+
        scale_y_continuous(limits=c(0,1))+
        scale_color_manual(values = c(grade_colors))+ #+ ggtitle(title) +
        theme_bw()+ 
        #facet_grid(.~title)+
        theme(text = element_text(size=8),
              axis.text=element_text(size=8, face = "bold"),
              axis.title=element_text(size=8,face="bold"),
              legend.position = "none")+ 
        guides(colour = guide_legend(override.aes = list(size=4),nrow=4)) +
        xlab('Grade')+ ylab('F1-Score Maximized Cutoff Specificity')
      
      png(paste0(name,'.',c,'.spec.f1.ci.grade.png'),height = 1100,width=600,res = 300)
      print(plot1)
      dev.off()
      
      
      #youdens index ppv
      plot1 = ggplot(auc.res1,aes(x = var, y = f1cutpoint.ppv, col = var)) +
        geom_errorbar(aes(ymin=f1cutpoint.ppv.lower, ymax=f1cutpoint.ppv.upper), width=.2,
                      position=position_dodge(0.05))+
        geom_point()+
        scale_y_continuous(limits=c(0,1))+
        scale_color_manual(values = c(grade_colors))+ #+ ggtitle(title) +
        theme_bw()+ 
        #facet_grid(.~title)+
        theme(text = element_text(size=8),
              axis.text=element_text(size=8, face = "bold"),
              axis.title=element_text(size=8,face="bold"),
              legend.position = "none")+ 
        guides(colour = guide_legend(override.aes = list(size=4),nrow=4)) +
        xlab('Grade')+ ylab('F1-Score Maximized Cutoff PPV')
      
      png(paste0(name,'.',c,'.ppv.f1.ci.Grade.png'),height = 1100,width=600,res = 300)
      print(plot1)
      dev.off()
      
      #f1 index npv
      plot1 = ggplot(auc.res1,aes(x = var, y = f1cutpoint.npv, col = var)) +
        geom_errorbar(aes(ymin=f1cutpoint.npv.lower, ymax=f1cutpoint.npv.upper), width=.2,
                      position=position_dodge(0.05))+
        geom_point()+
        scale_y_continuous(limits=c(0,1))+
        scale_color_manual(values = c(grade_colors))+ #+ ggtitle(title) +
        theme_bw()+ 
        #facet_grid(.~title)+
        theme(text = element_text(size=8),
              axis.text=element_text(size=8, face = "bold"),
              axis.title=element_text(size=8,face="bold"),
              legend.position = "none")+ 
        guides(colour = guide_legend(override.aes = list(size=4),nrow=4)) +
        xlab('Grade')+ ylab('F1-Score Maximized Cutoff NPV')
      
      png(paste0(name,'.',c,'.npv.f1.ci.grade.png'),height = 1100,width=600,res = 300)
      print(plot1)
      dev.off()
      
    }
    
  
    #plotting age##
    ####BASELINE AGE#####
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
    
    age2 = merged.df.all.tmp[merged.df.all.tmp$age_group2 %in% c('45-55') | merged.df.all.tmp$Cancer =='Control',]
    age3 = merged.df.all.tmp[merged.df.all.tmp$age_group2 %in% c('55-65') | merged.df.all.tmp$Cancer =='Control',]
    age4 = merged.df.all.tmp[merged.df.all.tmp$age_group2 %in% c('65-75') | merged.df.all.tmp$Cancer =='Control',]
    
    auc.plot.all = summary.calc(pred.df.targ.collapse,merged.df.all.tmp)
    auc.plot.all[[2]]$var = 'All'
    auc.plot.all[[2]]$var.group = 'Age'
    
    
    
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
    
    
    auc.t =rbind(  auc.plot.all.s1[[3]],
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
      #facet_grid2(.~title)+
      theme(text = element_text(size=8),
            axis.text=element_text(size=8, face = "bold"),
            axis.title=element_text(size=8,face="bold"),
            legend.position = "none")+ guides(colour = guide_legend(override.aes = list(size=4),nrow=4)) +
      xlab('Time (Years)') + ylab('Time-Dependent CV AUROC')
    
    png(paste0(name,'.',c,'.timeauroc.age.png'),height = 1100,width=1100,res = 300)
    print(plot1)
    dev.off()
    
    
    combined.auc.plot.all = rbind(
      
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
    
    png(paste0(name,'.',c,'.auroc.age.png'),height = 1100,width=1100,res = 300)
    print(plot1)
    dev.off()
    
    
    #age labels
    labels = combined.auc.plot.all[combined.auc.plot.all$l ==0,]
    labels$title = gsub('.new','',paste0(labels$var, ' AUROC: ',round(labels$auc,digits=3)))
    colors.filt = age_colors[names(age_colors) %in% labels$var]
    labels = labels[order(match(labels$var,names(colors.filt) )),]  
    labels$var = factor(labels$var,levels = names(colors.filt))
    combined.auc.plot.all$var = factor(as.character(combined.auc.plot.all$var), levels= names(colors.filt))
    
    
    
    plot1 = ggplot(combined.auc.plot.all,aes(x = fpr, y = tpr, col = var)) + geom_line(linewidth=1) +
      scale_color_manual(values = c(age_colors),labels =labels$title)+ #+ ggtitle(title) +
      theme_bw()+ 
      theme(text = element_text(size=8),
            axis.text=element_text(size=8, face = "bold"),
            axis.title=element_text(size=8,face="bold"),
            legend.position = "right")+ guides(colour = guide_legend(override.aes = list(size=4),nrow=4)) + xlab('False Positive Rate') + ylab('Sensitivity')
    
    png(paste0(name,'.',c,'.auroc.age.labs.png'),height = 1100,width=1100,res = 300)
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
      xlab('Baseline Age') + ylab('CV AUROC (95% CI)')
    
    png(paste0(name,'.',c,'.auroc.ci.age.png'),height = 1100,width=600,res = 300)
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
      xlab('Baseline Age') + ylab('CV Concordance Index (95% CI)')
    
    png(paste0(name,'.',c,'.ci.ci.age.png'),height = 1100,width=600,res = 300)
    print(plot1)
    dev.off()
    plot1 = ggplot(auc.res1,aes(x = var, y =se.spec.95 , col = var)) +
      geom_errorbar(aes(ymin=se.spec.95.lower, ymax=se.spec.95.upper), width=.2,
                    position=position_dodge(0.05))+
      geom_point()+
      scale_y_continuous(limits=c(0,1))+
      scale_color_manual(values = c(age_colors))+ #+ ggtitle(title) +
      theme_bw()+ 
      theme(text = element_text(size=8),
            axis.text=element_text(size=8, face = "bold"),
            axis.title=element_text(size=8,face="bold"),
            legend.position = "none")+ 
      guides(colour = guide_legend(override.aes = list(size=4),nrow=4)) +
      xlab('Baseline Age') + ylab('Sensitivity at 95% Specificity (95% CI)')
    
    png(paste0(name,'.',c,'.sens.spec95.ci.age.png'),height = 1100,width=600,res = 300)
    print(plot1)
    dev.off()
    
    plot1 = ggplot(auc.res1,aes(x = var, y = se.spec.90, col = var)) +
      geom_errorbar(aes(ymin=se.spec.90.lower, ymax=se.spec.90.upper), width=.2,
                    position=position_dodge(0.05))+
      geom_point()+
      scale_y_continuous(limits=c(0,1))+
      scale_color_manual(values = c(age_colors))+ #+ ggtitle(title) +
      theme_bw()+ 
      theme(text = element_text(size=8),
            axis.text=element_text(size=8, face = "bold"),
            axis.title=element_text(size=8,face="bold"),
            legend.position = "none")+ 
      guides(colour = guide_legend(override.aes = list(size=4),nrow=4)) +
      xlab('Baseline Age')+ ylab('Sensitivity at 90% Specificity (95% CI)')
    
    png(paste0(name,'.',c,'.sens.spec90.ci.age.png'),height = 1100,width=600,res = 300)
    print(plot1)
    dev.off()
    
    
    pred.df.targ.collapse.annotated = merge(pred.df.targ.collapse, merged.df.all.tmp[,c('GRP_Id','Cancer','Diagnosis_Time','filler','age_group2','GRADE_CD')],by='GRP_Id')
    
    #boxplot of score
    plot1 = ggplot(pred.df.targ.collapse.annotated[pred.df.targ.collapse.annotated$Cancer %in% c('Control','Breast'),],aes(x = age_group2, y = methylation_score, col = age_group2)) +
      geom_boxplot(outlier.shape=NA)+geom_jitter(width=0.1)+
      scale_color_manual(values = c(age_colors))+
      theme_bw()+ 
      scale_y_continuous(limits=c(0,1))+
      theme(text = element_text(size=8),
            axis.text=element_text(size=8, face = "bold"),
            axis.title=element_text(size=8,face="bold"),
            legend.position = "none")+ 
      guides(colour = guide_legend(override.aes = list(size=4),nrow=4)) +
      xlab('Baseline Age') + ylab('Methylation Score')
    
    png(paste0(name,'.age.methscore.',c,'.png'),height = 1100,width=600,res = 300)
    print(plot1)
    dev.off()
    
    
    
    
    #youdens index sens
    plot1 = ggplot(auc.res1,aes(x = var, y = jcutpoint.sens, col = var)) +
      geom_errorbar(aes(ymin=jcutpoint.sens.lower, ymax=jcutpoint.sens.upper), width=.2,
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
      xlab('Baseline Age')+ ylab('Youden\'s Index Cutoff Sensitivity')
    
    png(paste0(name,'.',c,'.sens.youden.ci.age.png'),height = 1100,width=600,res = 300)
    print(plot1)
    dev.off()
    
    #youdens index spec
    plot1 = ggplot(auc.res1,aes(x = var, y = jcutpoint.spec, col = var)) +
      geom_errorbar(aes(ymin=jcutpoint.spec.lower, ymax=jcutpoint.spec.upper), width=.2,
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
      xlab('Baseline Age')+ ylab('Youden\'s Index Cutoff Specificity')
    
    png(paste0(name,'.',c,'.spec.youden.ci.age.png'),height = 1100,width=600,res = 300)
    print(plot1)
    dev.off()
    
    
    #youdens index ppv
    plot1 = ggplot(auc.res1,aes(x = var, y = jcutpoint.ppv, col = var)) +
      geom_errorbar(aes(ymin=jcutpoint.ppv.lower, ymax=jcutpoint.ppv.upper), width=.2,
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
      xlab('Baseline Age')+ ylab('Youden\'s Index Cutoff PPV')
    
    png(paste0(name,'.',c,'.ppv.youden.ci.age.png'),height = 1100,width=600,res = 300)
    print(plot1)
    dev.off()
    
    #youdens index npv
    plot1 = ggplot(auc.res1,aes(x = var, y = jcutpoint.npv, col = var)) +
      geom_errorbar(aes(ymin=jcutpoint.npv.lower, ymax=jcutpoint.npv.upper), width=.2,
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
      xlab('Baseline Age')+ ylab('Youden\'s Index Cutoff NPV')
    
    png(paste0(name,'.',c,'.npv.youden.ci.age.png'),height = 1100,width=600,res = 300)
    print(plot1)
    dev.off()
    
    
    #f1 sens
    plot1 = ggplot(auc.res1,aes(x = var, y = f1cutpoint.sens, col = var)) +
      geom_errorbar(aes(ymin=f1cutpoint.sens.lower, ymax=f1cutpoint.sens.upper), width=.2,
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
      xlab('Baseline Age')+ ylab('F1-Score Maximized Cutoff Sensitivity')
    
    png(paste0(name,'.',c,'.sens.f1.ci.age.png'),height = 1100,width=600,res = 300)
    print(plot1)
    dev.off()
    
    #f1 index spec
    plot1 = ggplot(auc.res1,aes(x = var, y = f1cutpoint.spec, col = var)) +
      geom_errorbar(aes(ymin=f1cutpoint.spec.lower, ymax=f1cutpoint.spec.upper), width=.2,
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
      xlab('Baseline Age')+ ylab('F1-Score Maximized Cutoff Specificity')
    
    png(paste0(name,'.',c,'.spec.f1.ci.age.png'),height = 1100,width=600,res = 300)
    print(plot1)
    dev.off()
    
    
    #youdens index ppv
    plot1 = ggplot(auc.res1,aes(x = var, y = f1cutpoint.ppv, col = var)) +
      geom_errorbar(aes(ymin=f1cutpoint.ppv.lower, ymax=f1cutpoint.ppv.upper), width=.2,
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
      xlab('Baseline Age')+ ylab('F1-Score Maximized Cutoff PPV')
    
    png(paste0(name,'.',c,'.ppv.f1.ci.age.png'),height = 1100,width=600,res = 300)
    print(plot1)
    dev.off()
    
    #f1 index npv
    plot1 = ggplot(auc.res1,aes(x = var, y = f1cutpoint.npv, col = var)) +
      geom_errorbar(aes(ymin=f1cutpoint.npv.lower, ymax=f1cutpoint.npv.upper), width=.2,
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
      xlab('Baseline Age')+ ylab('F1-Score Maximized Cutoff NPV')
    
    png(paste0(name,'.',c,'.npv.f1.ci.age.png'),height = 1100,width=600,res = 300)
    print(plot1)
    dev.off()
    
    
    
    
    
    
    ####plotting dxage#####
    
    library(ggh4x)
    age_grouping3 = function(age) {
      tmp = ifelse(age >= 30 & age < 50, '30-50', age)
      tmp = ifelse(age >= 50 & age < 60, '50-60', tmp)
      tmp = ifelse(age >= 60 & age < 70, '60-70', tmp)
      tmp = ifelse(age >= 70, '70-80', tmp)
      
      age_groups = c('30-50','50-60','60-70','70-80')
      tmp = factor(tmp, levels = age_groups) 
    }
    #targline 
    merged.df.all.tmp$dx.age = ifelse(merged.df.all.tmp$Cancer != 'Control',merged.df.all.tmp$SDC_AGE_CALC + merged.df.all.tmp$diff_in_days/365,merged.df.all.tmp$SDC_AGE_CALC)
    merged.df.all.tmp$age_group3 = age_grouping3(merged.df.all.tmp$dx.age)
    
    age2 = merged.df.all.tmp[merged.df.all.tmp$age_group3 %in% c('50-60') | merged.df.all.tmp$Cancer =='Control',]
    age3 = merged.df.all.tmp[merged.df.all.tmp$age_group3 %in% c('60-70') | merged.df.all.tmp$Cancer =='Control',]
    age4 = merged.df.all.tmp[merged.df.all.tmp$age_group3 %in% c('70-80') | merged.df.all.tmp$Cancer =='Control',]
    
    
    
    auc.plot.all.s0 = summary.calc(pred.df.targ.collapse,merged.df.all.tmp)
    auc.plot.all.s0[[2]]$var = 'All'
    auc.plot.all.s0[[2]]$var.group = 'DxAge'
    
    
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
      auc.plot.all.s2[[2]][1,],
      auc.plot.all.s3[[2]][1,],
      auc.plot.all.s4[[2]][1,])
    auc.res1$title = 'DxAge'
    combined.auroc = rbind(combined.auroc,auc.res1[,colnames(auc.res1) %in% colnames(combined.auroc)])
    
    dxage_colors = c('All' = "#000004FF",'30-50' = '#C5C392','50-60'= '#FFB20F','60-70'='#FF4B3E', '70-80'='#972D07')
    
    
    #annotating auc t
    auc.plot.all.s0[[3]]$var = 'All'
    auc.plot.all.s2[[3]]$var = '50-60'
    auc.plot.all.s3[[3]]$var = '60-70'
    auc.plot.all.s4[[3]]$var = '70-80'
    
    
    auc.t =rbind( 
      #   auc.plot.all.s0[[3]],
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
    
    png(paste0(name,'.',c,'.timeauroc.dxage.png'),height = 1100,width=1100,res = 300)
    print(plot1)
    dev.off()
    
    
    combined.auc.plot.all = rbind(
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
    
    png(paste0(name,'.',c,'.auroc.dxage.png'),height = 1100,width=1100,res = 300)
    print(plot1)
    dev.off()
    
    #age labels
    labels = combined.auc.plot.all[combined.auc.plot.all$l ==0,]
    labels$title = gsub('.new','',paste0(labels$var, ' AUROC: ',round(labels$auc,digits=3)))
    colors.filt = dxage_colors[names(dxage_colors) %in% labels$var]
    labels = labels[order(match(labels$var,names(colors.filt) )),]  
    labels$var = factor(labels$var,levels = names(colors.filt))
    combined.auc.plot.all$var = factor(as.character(combined.auc.plot.all$var), levels= names(colors.filt))
    
    
    
    plot1 = ggplot(combined.auc.plot.all,aes(x = fpr, y = tpr, col = var)) + geom_line(linewidth=1) +
      scale_color_manual(values = c(dxage_colors),labels =labels$title)+ #+ ggtitle(title) +
      theme_bw()+ 
      theme(text = element_text(size=8),
            axis.text=element_text(size=8, face = "bold"),
            axis.title=element_text(size=8,face="bold"),
            legend.position = "right")+ guides(colour = guide_legend(override.aes = list(size=4),nrow=4)) + xlab('False Positive Rate') + ylab('Sensitivity')
    
    png(paste0(name,'.',c,'.auroc.dxage.labs.png'),height = 1100,width=1100,res = 300)
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
      xlab('Diagnosis Age') + ylab('CV AUROC (95% CI)')
    
    png(paste0(name,'.',c,'.auroc.ci.dxage.png'),height = 1100,width=600,res = 300)
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
    
    png(paste0(name,'.',c,'.ci.ci.dxage.png'),height = 1100,width=600,res = 300)
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
    
    png(paste0(name,'.',c,'.sens.spec95.ci.dxage.png'),height = 1100,width=600,res = 300)
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
    
    png(paste0(name,'.',c,'.sens.spec90.ci.dxage.png'),height = 1100,width=600,res = 300)
    print(plot1)
    dev.off()
    
    #age labels
    labels = combined.auc.plot.all[combined.auc.plot.all$l ==0,]
    labels$title = gsub('.new','',paste0(labels$var, ' AUROC: ',round(labels$auc,digits=3)))
    colors.filt = dxage_colors[names(dxage_colors) %in% labels$var]
    labels = labels[order(match(labels$var,names(colors.filt) )),]  
    labels$var = factor(labels$var,levels = names(colors.filt))
    combined.auc.plot.all$var = factor(as.character(combined.auc.plot.all$var), levels= names(colors.filt))
    
    
    plot1 = ggplot(combined.auc.plot.all,aes(x = fpr, y = tpr, col = var)) + geom_line(linewidth=1) +
      scale_color_manual(values = c(dxage_colors),labels = labels$title)+ #+ ggtitle(title) +
      theme_bw()+ 
      theme(text = element_text(size=8),
            axis.text=element_text(size=8, face = "bold"),
            axis.title=element_text(size=8,face="bold"),
            legend.position = "right")+ guides(colour = guide_legend(override.aes = list(size=4),ncol=1)) + 
      xlab('False Positive Rate') + ylab('Sensitivity')
    
    png(paste0(name,'.',c,'.auroc.dxage.labs.png'),height = 1100,width=600,res = 300)
    print(plot1)
    dev.off()
    
    pred.df.targ.collapse.annotated = merge(pred.df.targ.collapse, merged.df.all.tmp[,c('GRP_Id','Cancer','Diagnosis_Time','filler','age_group3','GRADE_CD')],by='GRP_Id')
    
    #boxplot of score
    plot1 = ggplot(pred.df.targ.collapse.annotated[pred.df.targ.collapse.annotated$Cancer %in% c('Control','Breast'),],aes(x = age_group3, y = methylation_score, col = age_group3)) +
      geom_boxplot(outlier.shape=NA)+geom_jitter(width=0.1)+
      scale_color_manual(values = c(dxage_colors))+ #+ ggtitle(title) +
      theme_bw()+ 
      scale_y_continuous(limits=c(0,1))+
      
      #facet_grid2(. ~ Cancer)+
      
      theme(text = element_text(size=8),
            axis.text=element_text(size=8, face = "bold"),
            axis.title=element_text(size=8,face="bold"),
            legend.position = "none")+ 
      guides(colour = guide_legend(override.aes = list(size=4),nrow=4)) +
      xlab('Diagnosis Age') + ylab('Methylation Score')
    
    png(paste0(name,'.dxage.methscore.',c,'.png'),height = 1100,width=600,res = 300)
    print(plot1)
    dev.off()
    
    
    
    
    #youdens index sens
    plot1 = ggplot(auc.res1,aes(x = var, y = jcutpoint.sens, col = var)) +
      geom_errorbar(aes(ymin=jcutpoint.sens.lower, ymax=jcutpoint.sens.upper), width=.2,
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
      xlab('Diagnosis Age')+ ylab('Youden\'s Index Cutoff Sensitivity')
    
    png(paste0(name,'.',c,'.sens.youden.ci.dxage.png'),height = 1100,width=600,res = 300)
    print(plot1)
    dev.off()
    
    #youdens index spec
    plot1 = ggplot(auc.res1,aes(x = var, y = jcutpoint.spec, col = var)) +
      geom_errorbar(aes(ymin=jcutpoint.spec.lower, ymax=jcutpoint.spec.upper), width=.2,
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
      xlab('Diagnosis Age')+ ylab('Youden\'s Index Cutoff Specificity')
    
    png(paste0(name,'.',c,'.spec.youden.ci.dxage.png'),height = 1100,width=600,res = 300)
    print(plot1)
    dev.off()
    
    
    #youdens index ppv
    plot1 = ggplot(auc.res1,aes(x = var, y = jcutpoint.ppv, col = var)) +
      geom_errorbar(aes(ymin=jcutpoint.ppv.lower, ymax=jcutpoint.ppv.upper), width=.2,
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
      xlab('Diagnosis Age')+ ylab('Youden\'s Index Cutoff PPV')
    
    png(paste0(name,'.',c,'.ppv.youden.ci.dxage.png'),height = 1100,width=600,res = 300)
    print(plot1)
    dev.off()
    
    #youdens index npv
    plot1 = ggplot(auc.res1,aes(x = var, y = jcutpoint.npv, col = var)) +
      geom_errorbar(aes(ymin=jcutpoint.npv.lower, ymax=jcutpoint.npv.upper), width=.2,
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
      xlab('Diagnosis Age')+ ylab('Youden\'s Index Cutoff NPV')
    
    png(paste0(name,'.',c,'.npv.youden.ci.dxage.png'),height = 1100,width=600,res = 300)
    print(plot1)
    dev.off()
    
    
    #f1 sens
    plot1 = ggplot(auc.res1,aes(x = var, y = f1cutpoint.sens, col = var)) +
      geom_errorbar(aes(ymin=f1cutpoint.sens.lower, ymax=f1cutpoint.sens.upper), width=.2,
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
      xlab('Diagnosis Age')+ ylab('F1-Score Maximized Cutoff Sensitivity')
    
    png(paste0(name,'.',c,'.sens.f1.ci.dxage.png'),height = 1100,width=600,res = 300)
    print(plot1)
    dev.off()
    
    #f1 index spec
    plot1 = ggplot(auc.res1,aes(x = var, y = f1cutpoint.spec, col = var)) +
      geom_errorbar(aes(ymin=f1cutpoint.spec.lower, ymax=f1cutpoint.spec.upper), width=.2,
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
      xlab('Diagnosis Age')+ ylab('F1-Score Maximized Cutoff Specificity')
    
    png(paste0(name,'.',c,'.spec.f1.ci.dxage.png'),height = 1100,width=600,res = 300)
    print(plot1)
    dev.off()
    
    
    #youdens index ppv
    plot1 = ggplot(auc.res1,aes(x = var, y = f1cutpoint.ppv, col = var)) +
      geom_errorbar(aes(ymin=f1cutpoint.ppv.lower, ymax=f1cutpoint.ppv.upper), width=.2,
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
      xlab('Diagnosis Age')+ ylab('F1-Score Maximized Cutoff PPV')
    
    png(paste0(name,'.',c,'.ppv.f1.ci.dxage.png'),height = 1100,width=600,res = 300)
    print(plot1)
    dev.off()
    
    #f1 index npv
    plot1 = ggplot(auc.res1,aes(x = var, y = f1cutpoint.npv, col = var)) +
      geom_errorbar(aes(ymin=f1cutpoint.npv.lower, ymax=f1cutpoint.npv.upper), width=.2,
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
      xlab('Diagnosis Age')+ ylab('F1-Score Maximized Cutoff NPV')
    
    png(paste0(name,'.',c,'.npv.f1.ci.dxage.png'),height = 1100,width=600,res = 300)
    print(plot1)
    dev.off()
    
    #specificty/npv plot youdens index cutoff
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
    
    png(paste0(name,'.',c,'.npv.spec.youdenoverall.png'),height = 1100,width=800,res = 300)
    print(plot1)
    dev.off()
    
    png(paste0(name,'.',c,'.npv.spec.youdenoverall.short.png'),height = 600,width=800,res = 300)
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
      theme_bw()+ 
      theme(text = element_text(size=8),
            axis.text=element_text(size=8, face = "bold"),
            axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
            axis.title.y=element_text(size=8,face="bold"),
            axis.title.x=element_blank(),
            legend.position = "none")+ 
      guides(colour = guide_legend(override.aes = list(size=4),nrow=4)) +
      xlab('')+ ylab('Value')
    
    png(paste0(name,'.',c,'.npv.spec.f1overall.png'),height = 1100,width=600,res = 300)
    print(plot1)
    dev.off()
    
    png(paste0(name,'.',c,'.npv.spec.f1overall.short.png'),height = 600,width=800,res = 300)
    print(plot1)
    dev.off()
    
    combined.auroc$model = c
    auc.t.combined$model = c
    return.df.final.auroc = rbind(return.df.final.auroc,combined.auroc)
    return.df.final.auct = rbind(return.df.final.auct,auc.t.combined)
    
    
    ####surv plots#####
    
    #km curves.all#
    pca.coords =merge(pred.df.targ.collapse.annotated,merged.df.all.tmp) # merge(pred.df.targ.collapse.annotated, merged.df.all.tmp
    cp.F1_score= cutpointr(pred.df.targ.collapse.all$methylation_score,pred.df.targ.collapse.all$reported, method = maximize_metric, metric = F1_score)$optimal_cutpoint
    #cp.youden
    #cp.F1_score
    if (dx.all == T) {
      max.time = 9
    } else {
      max.time=5
    }
    cp.median = median(pred.df.targ.collapse.all$methylation_score)
    
    cutoff.types = c(
      'median'=cp.median,
      'cutoff'=score.cutoff)
    for (s in names(cutoff.types)[2]) {
      score.cutoff.val = cutoff.types[s]
      pca.coords$median.g = factor(ifelse(pca.coords$methylation_score >= score.cutoff.val,'Higher Test Median','Lower Test Median'),levels=c('Lower Test Median','Higher Test Median'))
      
      pca.coords$Cancer.g= ifelse(pca.coords$reported == 'Control',0,1)
      male.weights= weightsf.males(pca.coords)
      cph = summary(coxph(Surv(censorship_time, Cancer.g) ~ median.g, data =  pca.coords,weights=male.weights))$logtest[3]
      hr = round(summary(coxph(Surv(censorship_time, Cancer.g) ~ median.g, data =  pca.coords,weights=male.weights))$coefficients[2],digits=3)
      median.g = pca.coords$median.g
      a=coxph(formula = Surv(censorship_time, Cancer.g) ~ median.g, data =  pca.coords,weights=male.weights)
      b = survfit(a,newdata=data.frame(365*seq(1:max.time)))
      plot4 = ggsurvplot(
        fit = survfit(Surv(censorship_time/365, Cancer.g) ~ median.g, data =  pca.coords,weights = male.weights), 
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
      
      png(paste0(name,".scores.surv_median.all.",c,'.',s,".png"),height = 1100, width = 1100,res=300)
      print(plot4)
      dev.off()
      
      plot4 = ggsurvplot(
        fit = survfit(Surv(censorship_time/365, Cancer.g) ~ median.g, data =  pca.coords,weights = male.weights), 
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
      
      png(paste0(name,".scores.surv_median.all.",c,'.',s,".full.png"),height = 1100, width = 1100,res=300)
      print(plot4)
      dev.off()
      
      pca.coords$quartile <- dplyr::ntile(pca.coords$methylation_score, 4) 
      quartile = pca.coords$quartile 
      cph = summary(coxph(Surv(censorship_time, Cancer.g) ~ quartile, data =  pca.coords,weights = male.weights))$logtest[3]
      hr = round(summary(coxph(Surv(censorship_time, Cancer.g) ~ quartile, data =  pca.coords,weights=male.weights))$coefficients[2],digits=3)
      survfit.obj = survfit(Surv(censorship_time/365, Cancer.g) ~ quartile, data =  pca.coords,weights = male.weights)
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
        ylim=c(0.7,1),
        ggtheme = theme_bw(),
        
        font.tickslab = c(8,'bold'), 
        font.main = c(8, "bold",'black'),
        font.x = c(8, "bold",'black'),
        font.y = c(8, "bold",'black'),
        alpha = 0.6,
        font.legend = c(8), risk.table.fontsize = 5, risk.table.col = "black", size = 1)  + 
        ggtitle(paste0(hr,'.',formatC(cph, format = "e", digits = 2)) )  #+ the #+ theme(legend.position = 'none')
      
      
      png(paste0(name,".scores.surv_quart.all.",c,'.',s,".cut.png"),height = 1100, width = 1100,res=300)
      print(plot4)
      dev.off()
      
      pca.coords$quartile <- ifelse(pca.coords$methylation_score >=0 & pca.coords$methylation_score < 0.25, c(1),
                                    ifelse(pca.coords$methylation_score >=0.25 & pca.coords$methylation_score < 0.5, c(2),
                                           ifelse(pca.coords$methylation_score >=0.5 & pca.coords$methylation_score < 0.75, c(3),
                                                  ifelse(pca.coords$methylation_score >=0.75 & pca.coords$methylation_score < 1, c(4),'Other'))))
      quartile = pca.coords$quartile 
      cph = summary(coxph(Surv(censorship_time, Cancer.g) ~ quartile, data =  pca.coords,weights = male.weights))$logtest[3]
      hr = round(summary(coxph(Surv(censorship_time, Cancer.g) ~ quartile, data =  pca.coords,weights=male.weights))$coefficients[2],digits=3)
      survfit.obj = survfit(Surv(censorship_time/365, Cancer.g) ~ quartile, data =  pca.coords,weights = male.weights)
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
        ylim=c(0.70,1),
        ggtheme = theme_bw(),
        
        font.tickslab = c(8,'bold'), 
        font.main = c(8, "bold",'black'),
        font.x = c(8, "bold",'black'),
        font.y = c(8, "bold",'black'),
        alpha = 0.6,
        font.legend = c(8), risk.table.fontsize = 5, risk.table.col = "black", size = 1)  + 
        ggtitle(paste0(hr,'.',formatC(cph, format = "e", digits = 2)) )  
      
      
      png(paste0(name,".scores.surv_quart.breaks.",c,'.',s,".cut.png"),height = 1100, width = 1100,res=300)
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
      
      png(paste0(name,".scores.surv_quart.all.",c,'.',s,".all.png"),height = 1100, width = 1100,res=300)
      print(plot4)
      dev.off()
      
      
      km=survfit(formula = Surv(censorship_time, Cancer.g) ~ median.g, data =  pca.coords,weights=male.weights)
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
        ylab('Prostate Cancer Probability (%)') + xlab('Time (Years)') 
      png(paste0(name,".scores.weighted.risk.all.",c,'.',s,".all.png"),height = 900, width = 1100,res=300)
      print(plot1)
      dev.off()
      
      asdf
      
      #split by gleason score
      if (gs.plot == T){
        for (gs in c('6','7','8/9')) {
          if (gs != '8/9') {
            pca.coords.filt =  pca.coords[pca.coords$Gleason.score %in% c('Control',gs),]
            male.weights.filt = male.weights[pca.coords$Gleason.score %in% c('Control',gs)]
          } else {
            pca.coords.filt =  pca.coords[pca.coords$Gleason.score %in% c('Control','8','9'),]
            male.weights.filt = male.weights[pca.coords$Gleason.score %in% c('Control','8','9')]
          }
          
          gs.name =gsub('\\/','',gs)
          cph = summary(coxph(Surv(censorship_time, Cancer.g) ~ median.g, data = pca.coords.filt,weights=male.weights.filt))$logtest[3]
          hr = round(summary(coxph(Surv(censorship_time, Cancer.g) ~ median.g, data =  pca.coords.filt,weights=male.weights.filt))$coefficients[2],digits=3)
          median.g = pca.coords.filt$median.g
          a=coxph(formula = Surv(censorship_time, Cancer.g) ~ median.g, data =  pca.coords.filt,weights=male.weights.filt)
          b = survfit(a,newdata=data.frame(365*seq(1:5)))
          plot4 = ggsurvplot(
            fit = survfit(Surv(censorship_time/365, Cancer.g) ~ median.g, data =  pca.coords.filt,weights = male.weights.filt), 
            risk.table = F,
            #cumevents = TRUE,
            palette = c('#2a9d8f','#e9c46a'),
            legend.labs=c("Higher Test Median","Lower Test Median"),
            xlim = c(0,10), xlab = "Time (Years)", ylab = c("Survival probability"), break.x.by = 1,
            ylim = c(0.85,1),
            ggtheme = theme_bw(),
            
            font.tickslab = c(8,'bold'), 
            font.main = c(8, "bold",'black'),
            font.x = c(8, "bold",'black'),
            font.y = c(8, "bold",'black'),
            alpha = 0.6,
            font.legend = c(8), risk.table.fontsize = 5, risk.table.col = "black", size = 1) + 
            ggtitle(paste(hr,'.',formatC(cph, format = "e", digits = 2))) #+ theme(legend.position = 'none')
          
          png(paste0(name,".scores.surv_median.gs.",gs.name,'.',c,'.',s,".cut.png"),height = 1100, width = 1100,res=300)
          print(plot4)
          dev.off()
          
          plot4 = ggsurvplot(
            fit = survfit(Surv(censorship_time/365, Cancer.g) ~ median.g, data =  pca.coords.filt,weights = male.weights.filt), 
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
          
          png(paste0(name,".scores.surv_median.gs.",gs.name,'.',c,'.',s,".all.png"),height = 1100, width = 1100,res=300)
          print(plot4)
          dev.off()
          
          
          pca.coords.filt$quartile <- dplyr::ntile(pca.coords.filt$methylation_score, 4) 
          quartile = pca.coords.filt$quartile 
          cph = summary(coxph(Surv(censorship_time, Cancer.g) ~ quartile, data =  pca.coords.filt,weights = male.weights.filt))$logtest[3]
          hr = round(summary(coxph(Surv(censorship_time, Cancer.g) ~ quartile, data =  pca.coords.filt,weights=male.weights.filt))$coefficients[2],digits=3)
          survfit.obj = survfit(Surv(censorship_time/365, Cancer.g) ~ quartile, data =  pca.coords.filt,weights = male.weights.filt)
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
            ylim=c(0.70,1),
            ggtheme = theme_bw(),
            
            font.tickslab = c(8,'bold'), 
            font.main = c(8, "bold",'black'),
            font.x = c(8, "bold",'black'),
            font.y = c(8, "bold",'black'),
            
            alpha = 0.6,
            font.legend = c(8), risk.table.fontsize = 5, risk.table.col = "black", size = 1)  + 
            ggtitle(paste0(hr,'.',formatC(cph, format = "e", digits = 2)) )  #+ the #+ theme(legend.position = 'none')
          
          
          png(paste0(name,".scores.surv_quart.gs.",gs.name,'.',c,'.',s,".png"),height = 1100, width = 1100,res=300)
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
          
          png(paste0(name,".scores.surv_quart.gs.",gs.name,'.',c,'.',s,".full.png"),height = 1100, width = 1100,res=300)
          print(plot4)
          dev.off()
          
          
          pca.coords$quartile <- ifelse(pca.coords$methylation_score >=0 & pca.coords$methylation_score < 0.25, c(1),
                                        ifelse(pca.coords$methylation_score >=0.25 & pca.coords$methylation_score < 0.5, c(2),
                                               ifelse(pca.coords$methylation_score >=0.5 & pca.coords$methylation_score < 0.75, c(3),
                                                      ifelse(pca.coords$methylation_score >=0.75 & pca.coords$methylation_score < 1, c(4),'Other'))))
          
          quartile = pca.coords.filt$quartile 
          cph = summary(coxph(Surv(censorship_time, Cancer.g) ~ quartile, data =  pca.coords.filt,weights = male.weights.filt))$logtest[3]
          hr = round(summary(coxph(Surv(censorship_time, Cancer.g) ~ quartile, data =  pca.coords.filt,weights=male.weights.filt))$coefficients[2],digits=3)
          survfit.obj = survfit(Surv(censorship_time/365, Cancer.g) ~ quartile, data =  pca.coords.filt,weights = male.weights.filt)
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
            ylim=c(0.7,1),
            ggtheme = theme_bw(),
            
            font.tickslab = c(8,'bold'), 
            font.main = c(8, "bold",'black'),
            font.x = c(8, "bold",'black'),
            font.y = c(8, "bold",'black'),
            
            alpha = 0.6,
            font.legend = c(8), risk.table.fontsize = 5, risk.table.col = "black", size = 1)  + 
            ggtitle(paste0(hr,'.',formatC(cph, format = "e", digits = 2)) )  #+ the #+ theme(legend.position = 'none')
          
          
          png(paste0(name,".scores.surv_quart.breaks.gs.",gs.name,'.',c,'.',s,".png"),height = 1100, width = 1100,res=300)
          print(plot4)
          dev.off()
          
          km=survfit(formula = Surv(censorship_time, Cancer.g) ~ median.g, data =  pca.coords.filt,weights=male.weights.filt)
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
            ylab('Prostate Cancer Probability (%)') + xlab('Time (Years)') 
          png(paste0(name,".scores.weighted.risk.gs.",gs.name,'.',c,'.',s,".all.png"),height = 900, width = 1100,res=300)
          print(plot1)
          dev.off()
          
          
          
        }
        
        
      }
      
      
      #split by diagnosis age
      for (dxage in c('50-60','60-70','70-80')) {
        pca.coords.filt =  pca.coords[as.character(pca.coords$age_group3) %in% c(dxage) | pca.coords$Cancer == 'Control',]
        male.weights.filt = male.weights[as.character(pca.coords$age_group3) %in% c(dxage)| pca.coords$Cancer == 'Control']
        
        
        cph = summary(coxph(Surv(censorship_time, Cancer.g) ~ median.g, data = pca.coords.filt,weights=male.weights.filt))$logtest[3]
        hr = round(summary(coxph(Surv(censorship_time, Cancer.g) ~ median.g, data =  pca.coords.filt,weights=male.weights.filt))$coefficients[2],digits=3)
        median.g = pca.coords.filt$median.g
        a=coxph(formula = Surv(censorship_time, Cancer.g) ~ median.g, data =  pca.coords.filt,weights=male.weights.filt)
        b = survfit(a,newdata=data.frame(365*seq(1:5)))
        
        if (length(unique(pca.coords.filt$median.g)) > 1) {
          plot4 = ggsurvplot(
            fit = survfit(Surv(censorship_time/365, Cancer.g) ~ median.g, data =  pca.coords.filt,weights = male.weights.filt), 
            risk.table = F,
            #cumevents = TRUE,
            palette = c('#2a9d8f','#e9c46a'),
            legend.labs=c("Higher Test Median","Lower Test Median"),
            xlim = c(0,10), xlab = "Time (Years)", ylab = c("Survival probability"), break.x.by = 1,
            ylim = c(0.85,1),
            ggtheme = theme_bw(),
            
            font.tickslab = c(8,'bold'), 
            font.main = c(8, "bold",'black'),
            font.x = c(8, "bold",'black'),
            font.y = c(8, "bold",'black'),
            alpha = 0.6,
            font.legend = c(8), risk.table.fontsize = 5, risk.table.col = "black", size = 1) + 
            ggtitle(paste(hr,'.',formatC(cph, format = "e", digits = 2))) #+ theme(legend.position = 'none')
          
          png(paste0(name,".scores.surv_median.dxage.",dxage,'.',c,'.',s,".cut.png"),height = 1100, width = 1100,res=300)
          print(plot4)
          dev.off()
          
          plot4 = ggsurvplot(
            fit = survfit(Surv(censorship_time/365, Cancer.g) ~ median.g, data =  pca.coords.filt,weights = male.weights.filt), 
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
          
          png(paste0(name,".scores.surv_median.dxage.",dxage,'.',c,'.',s,".all.png"),height = 1100, width = 1100,res=300)
          print(plot4)
          dev.off()
          
          
          pca.coords.filt$quartile <- dplyr::ntile(pca.coords.filt$methylation_score, 4) 
          quartile = pca.coords.filt$quartile 
          cph = summary(coxph(Surv(censorship_time, Cancer.g) ~ quartile, data =  pca.coords.filt,weights = male.weights.filt))$logtest[3]
          hr = round(summary(coxph(Surv(censorship_time, Cancer.g) ~ quartile, data =  pca.coords.filt,weights=male.weights.filt))$coefficients[2],digits=3)
          survfit.obj = survfit(Surv(censorship_time/365, Cancer.g) ~ quartile, data =  pca.coords.filt,weights = male.weights.filt)
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
            ylim=c(0.85,1),
            ggtheme = theme_bw(),
            
            font.tickslab = c(8,'bold'), 
            font.main = c(8, "bold",'black'),
            font.x = c(8, "bold",'black'),
            font.y = c(8, "bold",'black'),  alpha = 0.6,
            font.legend = c(8), risk.table.fontsize = 5, risk.table.col = "black", size = 1)  + 
            ggtitle(paste0(hr,'.',formatC(cph, format = "e", digits = 2)) )  #+ the #+ theme(legend.position = 'none')
          
          
          png(paste0(name,".scores.surv_quart.dxage.",dxage,'.',c,'.',s,".png"),height = 1100, width = 1100,res=300)
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
          
          png(paste0(name,".scores.surv_quart.dxage.",dxage,'.',c,'.',s,".full.png"),height = 1100, width = 1100,res=300)
          print(plot4)
          dev.off()
          
          pca.coords$quartile <- ifelse(pca.coords$methylation_score >=0 & pca.coords$methylation_score < 0.25, c(1),
                                        ifelse(pca.coords$methylation_score >=0.25 & pca.coords$methylation_score < 0.5, c(2),
                                               ifelse(pca.coords$methylation_score >=0.5 & pca.coords$methylation_score < 0.75, c(3),
                                                      ifelse(pca.coords$methylation_score >=0.75 & pca.coords$methylation_score < 1, c(4),'Other'))))
          quartile = pca.coords.filt$quartile 
          cph = summary(coxph(Surv(censorship_time, Cancer.g) ~ quartile, data =  pca.coords.filt,weights = male.weights.filt))$logtest[3]
          hr = round(summary(coxph(Surv(censorship_time, Cancer.g) ~ quartile, data =  pca.coords.filt,weights=male.weights.filt))$coefficients[2],digits=3)
          survfit.obj = survfit(Surv(censorship_time/365, Cancer.g) ~ quartile, data =  pca.coords.filt,weights = male.weights.filt)
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
            ylim=c(0.7,1),
            ggtheme = theme_bw(),
            
            font.tickslab = c(8,'bold'), 
            font.main = c(8, "bold",'black'),
            font.x = c(8, "bold",'black'),
            font.y = c(8, "bold",'black'),
            alpha = 0.6,
            font.legend = c(8), risk.table.fontsize = 5, risk.table.col = "black", size = 1)  + 
            ggtitle(paste0(hr,'.',formatC(cph, format = "e", digits = 2)) )  #+ the #+ theme(legend.position = 'none')
          
          
          png(paste0(name,".scores.surv_quart.break.dxage.",dxage,'.',c,'.',s,".png"),height = 1100, width = 1100,res=300)
          print(plot4)
          dev.off()
          
          
          km=survfit(formula = Surv(censorship_time, Cancer.g) ~ median.g, data =  pca.coords.filt,weights=male.weights.filt)
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
            ylab('Prostate Cancer Probability (%)') + xlab('Time (Years)') 
          png(paste0(name,".scores.weighted.risk.dxage.",dxage,'.',c,'.',s,".all.png"),height = 900, width = 1100,res=300)
          print(plot1)
          dev.off()
          
          
        }
        
        
        
        
      }
      
      
      
    }
    
    
    
  }
  
  
  
  
  
  
  return(list(return.df.final.auroc,return.df.final.auct))
  
}

combined.performance = list(return.df.final.auroc,return.df.final.auct) #performance stratified by subgroups



#!/usr/bin/env Rscript
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

auc_calc = function(prediction_table,labels = c('Control','Cancer')) {
  tmp = prediction_table
  tmp = tmp[order(-tmp$methylation_score),]
  #labels = c(tmp$reported)
  pred = prediction(predictions = c(tmp$methylation_score) ,labels =  tmp$reported, labels)
  perf_AUC=performance(pred,"auc") #Calculate the AUC value
  AUC=perf_AUC@y.values[[1]]
  return(AUC)
}
autosome_filt = function(medips.count_df) {
  tmp =gsub(':.*','',medips.count_df)
  return_df = medips.count_df[!tmp %in% c('chrY','chrX','chrM')]
  return(return_df)
}

options(expressions = 5e5)
library(DESeq2)
library("BiocParallel")
ncores = 20
register(MulticoreParam(ncores))


predictive.models.glm = function(targ.matrix1,
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
    targ.features1 = targ.features#[1:f]
    targ.matrix = targ.matrix1[,targ.features1]
    
    targ.matrix$GRP_Id = rownames(targ.matrix)
    
    targ.matrix=  merge(targ.matrix,merged.df,by='GRP_Id')
    
    
    targ.matrix.train = targ.matrix[targ.matrix$GRP_Id %in% train.set$GRP_Id,]
    targ.matrix.test = targ.matrix[targ.matrix$GRP_Id %in% test.set$GRP_Id,]
    
    rownames(targ.matrix.train)= targ.matrix.train$GRP_Id
    rownames(targ.matrix.test) = targ.matrix.test$GRP_Id
    
    targ.matrix.train$group = factor(ifelse(targ.matrix.train$group == 'Control','Control','Cancer'),levels = c('Control','Cancer'))
    targ.matrix.test$group = factor(ifelse(targ.matrix.test$group == 'Control','Control','Cancer'),levels=  c('Control','Cancer'))
    
    
    targ.matrix.train.model = model.matrix(~. -1,targ.matrix.train[,c(targ.features1)] )
    targ.matrix.test.model = model.matrix(~. -1,targ.matrix.test[,c(targ.features1)] )
    
    #logreg new
    ##
    
    
    #logreg old
    
    prediction.setup.lasso.prev= function(lambda.min, alpha,lasso.fit,test_set_matrix,features,model = 'logreg') {
      test_set_matrix.model= model.matrix(~. -1,targ.matrix.test[,features] )
      lasso.predict.prob = predict(lasso.fit,lambda = lambda.min,alpha=alpha ,newx=test_set_matrix.model, type="response") # check help(predict.glmnet)
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
    
    foldids = lapply(1:10, function(x) {
      tmp.list = split(targ.matrix.train,targ.matrix.train$group)
      tmp.list = lapply(tmp.list, function(x) {
        return.tmp = x
        return.tmp$foldid <- sample(rep(1:10, length.out = nrow(x)))
        return(return.tmp)
        
      } )
      return.df = do.call('rbind',tmp.list)
      rownames(return.df) = return.df$GRP_Id
      return.df = return.df[targ.matrix.train$GRP_Id,]
      return(return.df)
    } )
    start=Sys.time()
    tuning.df = NULL
    
    for (a in seq(0,1,0.25)) {
      print(a)
      MSEs = NULL
      for (i in 1:10){
        print(i)
        lasso.cv = cv.glmnet(targ.matrix.train.model, targ.matrix.train$group, 
                             family="binomial", type.measure="class",alpha=a,  nfolds = 10,
                             foldids = foldids[[i]]$foldid)
        MSEs <- cbind(MSEs, lasso.cv$cvm)
        
      }
      rownames(MSEs) <- lasso.cv$lambda
      lambda.min <- as.numeric(names(which.min(rowMeans(MSEs))))
      
      tuning.df.tmp = data.frame(MSE =min(rowMeans(MSEs)), lambda = lambda.min[1], alpha =a)
      tuning.df = rbind(tuning.df,tuning.df.tmp)
      
    }
    end=Sys.time()
    print(lambda.min)
    #nfolds=10,
    # lambda= lambda.min,
    #   type.measure = "class")
    alpha = tuning.df[which(tuning.df$MSE ==min(tuning.df$MSE)),'alpha'][1]
    lambda = tuning.df[which(tuning.df$MSE ==min(tuning.df$MSE)),'lambda'][1]
    
    lasso.fit = glmnet(as.matrix(targ.matrix.train.model), 
                       as.numeric(targ.matrix.train$group), 
                       family="binomial",
                       nfolds=10,
                       lambda= lambda,
                       alpha=alpha,
                       type.measure = "class")
    
    
    print(lasso.fit$lambda)
    
    prediction_table.logreg.old =prediction.setup.lasso.prev(lasso.fit$lambda,alpha,lasso.fit,targ.matrix.test,c(targ.features1),model = 'logreg')
    test_set_matrix.model= model.matrix(~. -1,targ.matrix.test[,targ.features1] )
    lasso.predict.prob = predict(lasso.fit, test_set_matrix.model, type="response") # check help(predict.glmnet)
    
    prediction_table.logreg.old$model = c('logreg.old.alphac')
    results.df = rbind(results.df,prediction_table.logreg.old)
    feature.coef =coef(lasso.fit,lambda.min)
    lasso.old.feature.weights =data.frame(PC =rownames(feature.coef), coef= feature.coef[1:length(feature.coef)])
    lasso.old.feature.weights$model = 'logreg.old.alphac'
    #combining feature weights
    feature.weights = lasso.old.feature.weights
    
    #alpha1 model
    lasso.fit = glmnet(as.matrix(targ.matrix.train.model), 
                       as.numeric(targ.matrix.train$group), 
                       family="binomial",
                       #nfolds=10,
                       lambda= tuning.df[which(tuning.df$alpha ==1),'lambda'],#tuning.df[which(tuning.df$alpha ==0),]$lambda,
                       alpha=1,
                       type.measure = "auc")
    
    

    prediction_table.logreg.old =prediction.setup.lasso.prev(lasso.fit$lambda,1,lasso.fit,targ.matrix.test,c(targ.features1),model = 'logreg')
    test_set_matrix.model= model.matrix(~. -1,targ.matrix.test[,targ.features1] )
    lasso.predict.prob = predict(lasso.fit, test_set_matrix.model, type="response") # check help(predict.glmnet)
    
    #
    prediction_table.logreg.old$model = c('logreg.old.alpha1')
    results.df = rbind(results.df,prediction_table.logreg.old)
    feature.coef =coef(lasso.fit,lambda.min)
    lasso.old.feature.weights =data.frame(PC =rownames(feature.coef), coef= feature.coef[1:length(feature.coef)])
    lasso.old.feature.weights$model = 'logreg.old.alpha1'
    #combining feature weights
    feature.weights = rbind(feature.weights,lasso.old.feature.weights)
    
    #alpha 0
    #alpha1 model
    lasso.fit = glmnet(as.matrix(targ.matrix.train.model), 
                       as.numeric(targ.matrix.train$group), 
                       family="binomial",
                       #nfolds=10,
                       lambda= tuning.df[which(tuning.df$alpha ==0),'lambda'],#tuning.df[which(tuning.df$alpha ==0),]$lambda,
                       alpha=0,
                       type.measure = "auc")
    
    

    prediction_table.logreg.old =prediction.setup.lasso.prev(lasso.fit$lambda,0,lasso.fit,targ.matrix.test,c(targ.features1),model = 'logreg')

    test_set_matrix.model= model.matrix(~. -1,targ.matrix.test[,targ.features1] )
    lasso.predict.prob = predict(lasso.fit, test_set_matrix.model, type="response") # check help(predict.glmnet)
    
    prediction_table.logreg.old$model = c('logreg.old.alpha0')
    results.df = rbind(results.df,prediction_table.logreg.old)
    feature.coef =coef(lasso.fit,lambda.min)
    lasso.old.feature.weights =data.frame(PC =rownames(feature.coef), coef= feature.coef[1:length(feature.coef)])
    lasso.old.feature.weights$model = 'logreg.old.alpha0'
    #combining feature weights
    feature.weights = rbind(feature.weights,lasso.old.feature.weights)
    
    
    
    
    #annotating scores
    results.df$feature = feature.type
    results.df$test = stat.test
    results.df$n.features = fl
    results.df.all = rbind(results.df.all,results.df)
    
    
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

####selecting top performing parameter#####
all.sample.info = readRDS('sample.info.RDS') 
discovery.set = all.sample.info[all.sample.info$data.partition == 'Discovery',]
ohs.test.set = all.sample.info[all.sample.info$data.partition != 'Discovery',]
savedir='/local/save/directory'
setwd(savedir)

dmrcalling = T
#female
if (dmrcalling == T) {
  dds = readRDS('ohs.1000.gehancer.dds.RDS')
  
  colData(dds)$total_reads= colData(dds)$total/1000000
  merged.df.filt.targ = discovery.set
  
  windows = rownames(dds)
  windows = windows[!grepl('chrX|chrY',windows)]
  
  
  targ.samples = discovery.set[discovery.set$Sex == 'Female',]
  targ.samples = unique(targ.samples[targ.samples$GRP_Id %in% colnames(dds),])
  dds = estimateSizeFactors(dds[windows,targ.samples$GRP_Id])
  dds$condition = dds$group
  
  colData(dds)$filler = factor(colData(dds)$filler,levels= c('MFiller','UFiller'))
  colData(dds)$group = factor(ifelse(colData(dds)$Cancer =='Control','Control','Cancer'),levels = c('Control','Cancer'))
  windows = rownames(dds)
  windows = windows[!grepl('chrX|chrY',windows)]
  
  ##breast cancer vs controls##
  merged.df.filt.targ = discovery.set[discovery.set$GRP_Id %in% colnames(dds),]
  merged.df.filt.targ = merged.df.filt.targ
  dds.filt = dds[windows,merged.df.filt.targ$GRP_Id]
  mm = model.matrix(~ filler + group, colData(dds.filt)) 
  
  dds.filt = dds.filt[rowSums(counts(dds.filt)) > 0,]
  dds.filt$condition = dds.filt$group
  
  
  ddssva <- DESeq(dds.filt,full = mm,parallel=T) #differential methylation analysis
  res.df = results(ddssva,contrast = list('groupCancer')) #generating results table
  saveRDS(res.df,'discovery.breast.genhancer.dmrs.RDS')
  
  
}


#####genhancer only ml####
results.df.all = NULL
feature = 'dmr'
res.df = readRDS('discovery.breast.genhancer.dmrs.RDS')

res.df$window =rownames(res.df)
directions = 'hyper'
dds = readRDS(paste0('ohs.1000.gehancer.dds.RDS'))
windows = rownames(dds)
windows = windows[!grepl('chrX|chrY',windows)]


targ.samples = all.sample.info
targ.samples = unique(targ.samples[targ.samples$GRP_Id %in% colnames(dds),])
dds = estimateSizeFactors(dds[windows,targ.samples$GRP_Id])
dds$condition = dds$group
sample.matrix = counts(dds, normalized = T)

fc.cutoff=0.25
results.df.all = NULL

for (fc in fc.cutoff) {
  for (d in directions) {
    
    res.df.targ = res.df[which(res.df$log2FoldChange> fc   & res.df$baseMean > 1),]
    res.df.targ = res.df.targ[rownames(res.df.targ) %in% rownames(sample.matrix),]
    res.df.targ = res.df.targ[order(res.df.targ$pvalue),]
    
    
    targ.samples = all.sample.info[all.sample.info$Sex == 'Female',]
    targ.samples = targ.samples[targ.samples$GRP_Id %in% colnames(sample.matrix),]
    
    std.matrix= data.frame(t(log2(sample.matrix[rownames(res.df.targ),unique(targ.samples$GRP_Id)]+1)),check.names=F) 
    rownames(std.matrix)= targ.samples$GRP_Id
    colnames(std.matrix) =rownames(res.df.targ)
    std.matrix = std.matrix[targ.samples$GRP_Id,]
    
    feature.weights.all = NULL
    targ.features = rownames(res.df.targ)
    
    train.set= discovery.set[discovery.set$Sex == 'Female',]
    test.set = ohs.test.set[ohs.test.set$GRP_Id %in% colnames(sample.matrix),]
    test.set = test.set[test.set$Sex == 'Female',]
    results.df = NULL
    n.features = 90
    for (fl in n.features){
      
      start= Sys.time()
      f = min(fl, length(targ.features))
      
      targ.features.df = res.df.targ
      targ.features.df$window = rownames(targ.features.df) 
      targ.features.base = targ.features.df$window[1:f]
      targ.matrix.base = std.matrix[,targ.features.base]
      targ.matrix.base$GRP_Id = rownames(targ.matrix.base)
      covariates = 'none'
      targ.matrix=  merge(targ.matrix.base,all.sample.info[,c('GRP_Id','group')],by='GRP_Id')
      targ.features1 =c(targ.features.base)  
      rownames(targ.matrix) = targ.matrix$GRP_Id
      
      
      
      
      tmp.res = mclapply(1:100, function(seednum) {
        set.seed(seednum)
        res.df.all = predictive.models.glm(targ.matrix,
                                           unique(all.sample.info[,c('GRP_Id','group','Cancer')]),
                                           train.set,
                                           test.set,
                                           targ.features =targ.features1,
                                           feats = fl,
                                           feature.type = 'gs',
                                           stat.test ='gs')
        perf.df = res.df.all[[1]] 
      },mc.cores=10)
      
      tmp.res = tmp.res[sapply(tmp.res,function(x) grepl('Error',x)[1] == F)]
      results.df.tmp = do.call('rbind',tmp.res)
      results.df.tmp = results.df.tmp[grepl('AIX.*',results.df.tmp$GRP_Id),]
      
      
      results.df.tmp$methylation_score=as.numeric(results.df.tmp$methylation_score)
      
      
      combined.collapse =ddply(results.df.tmp[,c('GRP_Id','reported','methylation_score','model')], c('GRP_Id','reported','model'),numcolwise(mean))
      combined.collapse = combined.collapse[combined.collapse$model == 'logreg.old.alpha1',]
      combined.collapse$auroc = auc_calc(combined.collapse)
      combined.collapse$predictions = ifelse(combined.collapse$methylation_score >= score.cutoff.breast,'Cancer','Control')
      combined.collapse$covariate=covar
      
      results.df.all = rbind(results.df.all,combined.collapse)
      
      results.df.all[results.df.all$model =='logreg.old.alpha1',]
      
      
      
    }

    saveRDS('validation.performance.RDS')

    
  }
  
  
}








#####ploting performance####
figdir ='/figures/directory/'

dir.create(figdir,recursive = T)
pred.df.targ = readRDS('validation.performance.RDS')
name = paste0(figdir,'brca.ohs.test');
merged.df.all = ohs.test.set
sample.info = ohs.test.set
dx.all = T
score.cutoff=0.425
cutpoint.use=T


#initially saved as a function, but R has trouble reading some survminer functions as a function
{
  c = 'All'
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




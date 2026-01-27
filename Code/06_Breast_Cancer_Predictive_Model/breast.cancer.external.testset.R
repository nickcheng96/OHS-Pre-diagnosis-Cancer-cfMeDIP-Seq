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
library(dplyr)


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

###reading in external test set files####

discovery.breast.genhancer.dmrs.RDS = readRDS('/.mounts/labs/awadallalab/private/ncheng/cfmedip_data/cptp_samples/fragmentation/aix13.updated1/all.sample.fragmentation/regulatory.counts/aix13.dv.octane.burgener.combined.genhancer.raw.counts.RDS') #upload

#if combining with octane etc
targ.samples.all = readRDS('/.mounts/labs/awadallalab/private/ncheng/cfmedip_data/cptp_samples/fragmentation/aix13.updated1/all.sample.fragmentation/regulatory.counts/aix.octane.burg.sample.info.RDS')
#creating dds object
library(DESeq2)
combined.sample.info = combined.sample.info[combined.sample.info$GRP_Id %in% colnames(combined.counts),]
dds <- DESeqDataSetFromMatrix(countData = combined.counts[,combined.sample.info$GRP_Id],
                              colData = combined.sample.info,
                              design= ~  Cancer ) #can add gender here but we're only looking at female samples currently



#####ml analysis######
#after loading in first half in validation.genhancer.R
matrices = c('Female.1000','All.1000','Female.400','All.400')
sf = c('AllChr.before','AutoChr.before')
sf = c('AllChr.before','AutoChr.before','SF.before')



#####genhancer only ml####

#setting wkdir
wkdir='/.mounts/labs/awadallalab/private/ncheng/cfmedip_data/cptp_samples/fragmentation/aix13.updated1/all.sample.fragmentation/breast.cancer.cv/regulatory.regions.v10/'
savedir=paste0(wkdir,'validation.breast.hrsplit/')

results.df.all = NULL
feature = 'dmr'

res.df = readRDS('discovery.breast.genhancer.dmrs.RDS')
targ.samples.all =readRDS('aix.octane.burg.sample.info.RDS')

res.df$window =rownames(res.df)
directions = c('hyper')

windows = rownames(dds)
windows = windows[!grepl('chrX|chrY',windows)]
targ.samples = unique(targ.samples.all[targ.samples.all$GRP_Id %in% colnames(dds),])
dds = estimateSizeFactors(dds[windows,targ.samples$GRP_Id])
dds$condition = dds$group
sample.matrix = counts(dds, normalized = T)


fc.cutoff=0.25
results.df.all = NULL
for (fc in fc.cutoff) {
  for (d in directions) {
    
    
    res.df.targ = res.df
    res.df.targ = res.df.targ[which(res.df.targ$log2FoldChange> fc   & res.df.targ$baseMean > 1),]
    res.df.targ = res.df.targ[rownames(res.df.targ) %in% rownames(sample.matrix),]
    res.df.targ = res.df.targ[order(res.df.targ$pvalue),]
    
    
    targ.samples = targ.samples.all[targ.samples.all$Sex == 'Female',]
    targ.samples = targ.samples[targ.samples$GRP_Id %in% colnames(sample.matrix),]
    
    mat = 'log2.std'
    std.matrix= data.frame(t(log2(sample.matrix[rownames(res.df.targ),unique(targ.samples$GRP_Id)]+1)),check.names=F) 
    
    
    rownames(std.matrix)= targ.samples$GRP_Id
    colnames(std.matrix) =rownames(res.df.targ)
    std.matrix = std.matrix[targ.samples$GRP_Id,]
    
    
    
    feature.weights.all = NULL
    targ.features = rownames(res.df.targ)
    
    train.set= discovery.set[discovery.set$Sex == 'Female',]
    test.set = targ.samples.all[!targ.samples.all$GRP_Id %in% train.set$GRP_Id,]
    test.set = test.set[test.set$Sex == 'Female',]
    results.df = NULL

    n.features=90
    score.cutoff.breast = 0.5
    for (filt in c('base')) {
      for (fl in n.features){

        start= Sys.time()
        f = min(fl, length(targ.features))
        
        #all features as features
        
        
        targ.features.df = res.df.targ
        targ.features.df$window = rownames(targ.features.df) 
        targ.features.base = targ.features.df$window[1:f]
        targ.matrix.base = std.matrix[,targ.features.base]
        targ.matrix.base$GRP_Id = rownames(targ.matrix.base)
        covariates = 'none'
        for (covar in covariates) {
          targ.samples$group = ifelse(targ.samples$Cancer == 'Control','Control','Cancer')
          targ.matrix=  merge(targ.matrix.base,targ.samples[,c('GRP_Id','group')],by='GRP_Id')
          targ.features1 =c(targ.features.base)  
          
          rownames(targ.matrix) = targ.matrix$GRP_Id
          
          
          targ.samples.all$group = ifelse(targ.samples.all$Cancer == 'Control','Control','Cancer')
          
          tmp.res = mclapply(1:100, function(seednum) {
            set.seed(seednum)
            
            res.df.all = predictive.models.glm(targ.matrix,
                                               unique(targ.samples.all[,c('GRP_Id','group','Cancer')]),
                                               train.set,
                                               test.set,
                                               targ.features =targ.features1,
                                               feats = fl,
                                               feature.type = 'genhancer',
                                               stat.test ='genhancer')
            perf.df = res.df.all[[1]] 
            return(perf.df)
          },mc.cores=10)
          
          tmp.res = tmp.res[sapply(tmp.res,function(x) grepl('Error',x)[1] == F)]
          results.df.tmp = do.call('rbind',tmp.res)
          results.df.tmp$methylation_score=as.numeric(results.df.tmp$methylation_score)
          combined.collapse =ddply(results.df.tmp[,c('GRP_Id','reported','methylation_score','model')], c('GRP_Id','reported','model'),numcolwise(mean))
          combined.collapse = combined.collapse[combined.collapse$model == 'logreg.old.alpha1',]
          combined.collapse$auroc = auc_calc(combined.collapse)
          combined.collapse$predictions = ifelse(combined.collapse$methylation_score >= score.cutoff.breast,'Cancer','Control')
          
          results.df.all = rbind(results.df.all,combined.collapse)
          

        }
        
        
        
        
      }
      
    }
  
    
    #
    
    
  }
  
  
}




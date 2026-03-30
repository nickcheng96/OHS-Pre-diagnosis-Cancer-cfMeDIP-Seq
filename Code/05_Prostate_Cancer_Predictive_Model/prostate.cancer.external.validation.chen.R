#!/usr/bin/env Rscript
library(DESeq2)
library(parallel)
library(caret)
library(glmnet)
library(ROCR)
library(plyr)
#####functions####
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
  cl <- makeCluster(5)
  registerDoParallel(cl)
  for (fl in feats){
    start= Sys.time()
    results.df = NULL
    f = min(fl, length(targ.features))
    
    #all features as features
    targ.features1 = targ.features#[1:f]
    targ.matrix = targ.matrix1[,targ.features1]
    
    targ.matrix$GRP_Id = rownames(targ.matrix)
    print('test')
    targ.matrix=  merge(targ.matrix,merged.df,by='GRP_Id')
    
    print('test1')
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
                       lambda= lambda,#tuning.df[which(tuning.df$alpha ==0),]$lambda,
                       alpha=alpha,
                       type.measure = "class")
    
    
    #lasso.fit = cv.glmnet(targ.matrix.train.model, targ.matrix.train$group, family="binomial")
    #lambda.min = lasso.fit$lambda.min
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
    
    
    #lasso.fit = cv.glmnet(targ.matrix.train.model, targ.matrix.train$group, family="binomial")
    #lambda.min = lasso.fit$lambda.min
    print(lasso.fit$lambda)
    
    prediction_table.logreg.old =prediction.setup.lasso.prev(lasso.fit$lambda,1,lasso.fit,targ.matrix.test,c(targ.features1),model = 'logreg')
    #
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
    
    
    #lasso.fit = cv.glmnet(targ.matrix.train.model, targ.matrix.train$group, family="binomial")
    #lambda.min = lasso.fit$lambda.min
    print(lasso.fit$lambda)
    
    prediction_table.logreg.old =prediction.setup.lasso.prev(lasso.fit$lambda,0,lasso.fit,targ.matrix.test,c(targ.features1),model = 'logreg')
    #
    test_set_matrix.model= model.matrix(~. -1,targ.matrix.test[,targ.features1] )
    lasso.predict.prob = predict(lasso.fit, test_set_matrix.model, type="response") # check help(predict.glmnet)
    
    #
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



#####merging with silencer combined set####
matrixdir='/matrix/directory/'
setwd(matrixdir)
#ohs sample info and count matrix
all.sample.info=readRDS('sample.info.RDS')
all.sample.info$Cancer.timing = ifelse(all.sample.info$Cancer != 'Control','Incident','Control')
ohs.dds = readRDS(paste0('ohs.1000.silencer.dds.RDS'))

#chen external counts and sample information
chen.raw.counts = readRDS('chen.hg38.silencer.raw.counts.RDS') 
chen.sample.info = readRDS('chen.medip.sample.info.RDS')

windows = rownames(ohs.dds)
windows = windows[!grepl('chrX|chrY',windows)]

#combining ohs counts with chen counts
targ.samples = all.sample.info
targ.samples = unique(targ.samples[targ.samples$GRP_Id %in% colnames(ohs.dds),])
ohs.raw.counts = counts(ohs.dds,normalize=F)

ohs.raw.counts.filt = data.frame(ohs.raw.counts[rownames(ohs.raw.counts) %in% rownames(chen.raw.counts),],stringsAsFactors = F,check.names = F)

chen.raw.counts.ordered = chen.raw.counts[,-1]

ohs.raw.counts.filt.merged = cbind(chen.raw.counts.ordered,ohs.raw.counts.filt)
combined.sample.info = rbind(chen.sample.info,all.sample.info[,colnames(chen.sample.info)])

#setting up counts
combined.sample.info.targ = combined.sample.info[combined.sample.info$GRP_Id %in% colnames(ohs.raw.counts.filt.merged),]
windows = rownames(ohs.raw.counts.filt.merged)
windows = windows[!grepl('chrX|chrY',windows)]
ohs.raw.counts.filt.merged = ohs.raw.counts.filt.merged[windows,]
combined.sample.info.targ = combined.sample.info.targ[combined.sample.info.targ$GRP_Id %in% colnames(ohs.raw.counts.filt.merged),]
dds <- DESeqDataSetFromMatrix(countData = ohs.raw.counts.filt.merged[,combined.sample.info.targ$GRP_Id],
                              colData = combined.sample.info.targ,
                              design= ~ Cancer ) #can add gender here but we're only looking at female samples currently


#normalized counts
dds = estimateSizeFactors(dds)
dds.matrix =counts(dds,normalize =T)
colnames(dds.matrix) =combined.sample.info.targ$GRP_Id
saveRDS(dds.matrix,'combined.chen.ohs.norm.counts.RDS')
saveRDS( ohs.raw.counts.filt.merged[,combined.sample.info$GRP_Id],'combined.chen.ohs.raw.counts.RDS') #combined count matrix for OHS + chen samples

#####applying ohs pre-diagnosis predictive model to chen et al. samples samples#####
res.df=readRDS('ohs.prostate.discovery.DMRs.RDS')

res.df$window =rownames(res.df)
directions = c('abs')
all.sample.info = combined.sample.info
sample.matrix = dds.matrix
fc.cutoff=0.25
results.df.all = NULL

for (fc in fc.cutoff) {
  for (d in directions) {
    sample.matrix = dds.matrix
    res.df.targ = res.df
    res.df.targ = res.df.targ[which(abs(res.df.targ$log2FoldChange) > 0.25 & res.df.targ$baseMean > 1),]
    res.df.targ = res.df.targ[rownames(res.df.targ) %in% rownames(sample.matrix),]
    res.df.targ = res.df.targ[order(res.df.targ$pvalue),]
    
    targ.samples = all.sample.info[all.sample.info$Sex == 'Male',]
    targ.samples = targ.samples[targ.samples$GRP_Id %in% colnames(sample.matrix),]
    
    
    std.matrix= data.frame(t((sample.matrix[rownames(res.df.targ),unique(targ.samples$GRP_Id)]+1)),check.names=F) 
    
    rownames(std.matrix)= targ.samples$GRP_Id
    colnames(std.matrix) =rownames(res.df.targ)
    std.matrix = std.matrix[targ.samples$GRP_Id,]
    
    
    feature.weights.all = NULL
    targ.features = rownames(res.df.targ)
    
    train.set= all.sample.info[all.sample.info$data.partition == 'Discovery' & all.sample.info$Sex == 'Male',]
    test.set = all.sample.info[!all.sample.info$GRP_Id %in% train.set$GRP_Id & all.sample.info$Sex == 'Male',]
    test.set = test.set[test.set$GRP_Id %in% colnames(sample.matrix),]
    test.set = test.set[test.set$Sex == 'Male',]
    results.df = NULL
    
    n.features=100
    score.cutoff = 0.112
    for (filt in c('base')) {
      print(filt)
      for (fl in n.features){
        print(fl)
        print(fl)
        start= Sys.time()
        f = min(fl, length(targ.features))
        
        #all features as features
        
        targ.features.df = res.df.targ
        targ.features.df$window = rownames(targ.features.df) 
        targ.features.base = targ.features.df$window[1:f]
        targ.matrix.base = std.matrix[,targ.features.base]
        targ.matrix.base$GRP_Id = rownames(targ.matrix.base)
        all.sample.info$group = factor(as.character(ifelse(all.sample.info$Cancer == 'Control','Control','Cancer')),levels = c('Control','Cancer'))
        targ.matrix=  merge(targ.matrix.base,all.sample.info[,c('GRP_Id','group')],by='GRP_Id')
        targ.features1 =c(targ.features.base)  
        
        rownames(targ.matrix) = targ.matrix$GRP_Id
        

        tmp.res = mclapply(1:100, function(seednum) {
          set.seed(seednum)
          res.df.all = predictive.models.glm(targ.matrix,
                                             unique(all.sample.info[all.sample.info$GRP_Id %in% rownames(targ.matrix),c('GRP_Id','group','Cancer')]),
                                             train.set,
                                             test.set,
                                             targ.features =targ.features1,
                                             feats = fl,
                                             feature.type = 'silencer',
                                             stat.test ='silencer')
          perf.df = res.df.all[[1]] 
          return(perf.df)
        },mc.cores=5)
        
        tmp.res = tmp.res[sapply(tmp.res,function(x) grepl('Error',x)[1] == F)]
        results.df.tmp = do.call('rbind',tmp.res)


        results.df.tmp$methylation_score=as.numeric(results.df.tmp$methylation_score)

        combined.collapse =ddply(results.df.tmp[,c('GRP_Id','reported','methylation_score','model')], c('GRP_Id','reported','model'),numcolwise(mean,na.rm=T))
        combined.collapse = combined.collapse[combined.collapse$model == 'logreg.old.alpha1',]
        combined.collapse$auroc = auc_calc(combined.collapse)
        combined.collapse$predictions = ifelse(combined.collapse$methylation_score >= score.cutoff,'Cancer','Control')
        
        
        results.df.all = rbind(results.df.all,combined.collapse)
        
        
        
        
        
        
      }
      
    }
    
    
    
  }
  
  
}

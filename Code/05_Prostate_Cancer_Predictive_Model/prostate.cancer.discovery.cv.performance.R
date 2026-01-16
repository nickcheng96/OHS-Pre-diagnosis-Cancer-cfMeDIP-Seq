#!/usr/bin/env Rscript
options(expressions = 5e5)
library(DESeq2)
library(parallel)
library(caret)
library(glmnet)
library(ROCR)
ncores = 10
register(MulticoreParam(ncores))

#####functions#####
train_test_partition_cv.singlesex = function(year2_diagnosis_samples, splits = 10, seed = 10) { 
  set.seed(seed)
  year2_diagnosis_samples$diagnosis_time_group =year2_diagnosis_samples$Cancer
  tmp_list = split(year2_diagnosis_samples, year2_diagnosis_samples[,c('Cancer','filler')])
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


auc_calc = function(prediction_table,labels = c('Control','Cancer')) {
  tmp = prediction_table
  tmp = tmp[order(-tmp$methylation_score),]
  pred = prediction(predictions = c(tmp$methylation_score) ,labels =  tmp$reported, labels)
  perf_AUC=performance(pred,"auc") #Calculate the AUC value
  AUC=perf_AUC@y.values[[1]]
  return(AUC)
}

predictive.models.glm = function(targ.matrix1,
                                 merged.df,
                                 train.set,
                                 test.set,
                                 targ.features,
                                 feats,
                                 feature.type,
                                 stat.test,no_cores=5, standardize = F) {

  auc_calc = function(prediction_table,labels = c('Control','Cancer')) {
    tmp = prediction_table
    tmp = tmp[order(-tmp$methylation_score),]
    #labels = c(tmp$reported)
    pred = prediction(predictions = c(tmp$methylation_score) ,labels =  tmp$reported, labels)
    perf_AUC=performance(pred,"auc") #Calculate the AUC value
    AUC=perf_AUC@y.values[[1]]
    return(AUC)
  }
  
  prediction.setup.lasso= function(lambda.min, alpha,lasso.fit,test_set_matrix,features,model = 'logreg') {
    test_set_matrix.model= model.matrix(~. -1,targ.matrix.test[,features] )
    lasso.predict.prob = predict(lasso.fit,lambda = lambda.min,alpha=alpha ,newx=test_set_matrix.model, type="response") # check help(predict.glmnet)
    predictions = ifelse(lasso.predict.prob > 0.5,'Cancer','Control')
    prediction_table = data.frame(GRP_Id=test_set_matrix$GRP_Id, 'predictions' = predictions[,1], reported = test_set_matrix$group, methylation_score = lasso.predict.prob[,1])
    prediction_table = prediction_table[order(-prediction_table$methylation_score),]
    prediction_table$auroc = auc_calc(prediction_table,c('Control','Cancer'))
    prediction_table$model = 'logreg'

    return(prediction_table)
  }
  
  results.df.all = NULL
  feature.weights.all = NULL

  for (fl in feats){
    start= Sys.time()
    results.df = NULL
    f = min(fl, length(targ.features))
    
    #all features as features
    targ.features1 = targ.features
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
    
    

    #glmnet iterating for optimal lambda and alpha
    foldids = lapply(1:10, function(x) {
      tmp.list = split(targ.matrix.train,targ.matrix.train$group)
      tmp.list = lapply(tmp.list, function(y) {
        return.tmp = y
        return.tmp$foldid <- sample(rep(1:10, length.out = nrow(y)))
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

    alpha = tuning.df[which(tuning.df$MSE ==min(tuning.df$MSE)),'alpha'][1]
    lambda = tuning.df[which(tuning.df$MSE ==min(tuning.df$MSE)),'lambda'][1]
    
    lasso.fit = glmnet(as.matrix(targ.matrix.train.model), 
                       as.numeric(targ.matrix.train$group), 
                       family="binomial",
                       nfolds=10,
                       lambda= lambda,#tuning.df[which(tuning.df$alpha ==0),]$lambda,
                       alpha=alpha,
                       type.measure = "class")

    
    prediction_table.logreg.old =prediction.setup.lasso(lasso.fit$lambda,alpha,lasso.fit,targ.matrix.test,c(targ.features1),model = 'logreg')
    test_set_matrix.model= model.matrix(~. -1,targ.matrix.test[,targ.features1] )
    lasso.predict.prob = predict(lasso.fit, test_set_matrix.model, type="response") # check help(predict.glmnet)
    
    prediction_table.logreg.old$model = c('logreg.old.alphac')
    results.df = rbind(results.df,prediction_table.logreg.old)
    feature.coef =coef(lasso.fit,lambda.min)
    lasso.old.feature.weights =data.frame(PC =rownames(feature.coef), coef= feature.coef[1:length(feature.coef)])
    lasso.old.feature.weights$model = 'logreg.old.alphac'
    feature.weights = lasso.old.feature.weights
    
    #glmnet fixing alpha at 1
    lasso.fit = glmnet(as.matrix(targ.matrix.train.model), 
                       as.numeric(targ.matrix.train$group), 
                       family="binomial",
                       lambda= tuning.df[which(tuning.df$alpha ==1),'lambda'],
                       alpha=1,
                       type.measure = "auc")
    
    
    prediction_table.logreg.old =prediction.setup.lasso(lasso.fit$lambda,1,lasso.fit,targ.matrix.test,c(targ.features1),model = 'logreg')
    test_set_matrix.model= model.matrix(~. -1,targ.matrix.test[,targ.features1] )
    lasso.predict.prob = predict(lasso.fit, test_set_matrix.model, type="response")
    prediction_table.logreg.old$model = c('logreg.old.alpha1')
    
    #combining prediction results
    results.df = rbind(results.df,prediction_table.logreg.old)
    feature.coef =coef(lasso.fit,lambda.min)
    lasso.old.feature.weights =data.frame(PC =rownames(feature.coef), coef= feature.coef[1:length(feature.coef)])
    lasso.old.feature.weights$model = 'logreg.old.alpha1'
    #combining feature weights
    feature.weights = rbind(feature.weights,lasso.old.feature.weights)
    
    #glmnet with alpha set at 0
    lasso.fit = glmnet(as.matrix(targ.matrix.train.model), 
                       as.numeric(targ.matrix.train$group), 
                       family="binomial",
                       lambda= tuning.df[which(tuning.df$alpha ==0),'lambda'],
                       alpha=0,
                       type.measure = "auc")
    
    prediction_table.logreg.old = prediction.setup.lasso(lasso.fit$lambda,0,lasso.fit,targ.matrix.test,c(targ.features1),model = 'logreg')
    test_set_matrix.model= model.matrix(~. -1,targ.matrix.test[,targ.features1] )
    lasso.predict.prob = predict(lasso.fit, test_set_matrix.model, type="response") 
    
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



####reading in annotations####
cpg_count = readRDS('hg38_cpg_window_300_count.RDS') #number of cpg sites across 300bp regions
cpg_count = cpg_count[cpg_count$count >= 5,] 

#loading sample information file 
sample.info.filt = readRDS('sample.info.RDS')
sample.info.filt = sample.info.filt[sample.info.filt$data.partition %in% c('Discovery') & sample.info.filt$Sex == 'Male',]
targets = list('Prostate' = sample.info.filt)


#setting directories
wkdir='/place/where/counts/are/stored/'
savedir='/place/to/save/files'
setwd(wkdir)
marker = c('silencer') 

#differential methylation calling
for (fold in foldno){
  for (i in names(targets)) {
    start = Sys.time()
    dds = readRDS(paste0(wkdir,'/ohs.1000.silencer.dds.RDS'))
    
    #selecting all male discovery set samples
    targ.samples = targets[[i]]
    targ.samples = targ.samples[targ.samples$GRP_Id %in% colnames(dds),]
    
    #partitioning discovery set into folds
    sample.split = train_test_partition_cv.singlesex(targ.samples, splits = 10, seed = seedno)
    
    #partitioning into train vs test folds based on partition
    train.set = sample.split[sample.split$fold != fold,]
    test.set = sample.split[sample.split$fold == fold,]

    #retaining only autosomes 
    windows = rownames(dds)
    windows = windows[!grepl('chrX|chrY',windows)]
    

    #making sure that filler and sample grouping are factorized
    colData(dds)$filler = factor(colData(dds)$filler,levels= c('MFiller','UFiller'))
    colData(dds)$group = factor(ifelse(colData(dds)$Cancer =='Control','Control','Cancer'),levels = c('Control','Cancer'))
    
    #ensuring correct samples and windows are used
    dds.filt = dds[windows,train.set$GRP_Id]

    #setting DMR calling model
    mm = model.matrix(~ filler + group, colData(dds.filt)) 
    
    #performing DMR calling on train set folds only
    ddssva <- DESeq(dds.filt,full = mm,parallel=T) 
    
    #extracting DMR results and storing as dataframe
    res.df = results(ddssva,contrast = list('groupCancer')) #generating results table
    res.df$feature =marker
    res.df$cancer= i
    res.df$window = rownames(res.df)
    res.df$seed = seedno
    res.df$fold = fold
    end = Sys.time()
    saveRDS(res.df,paste0(savedir,i,'.',seedno,'.',fold,'.dmr.RDS'))
    saveRDS(train.set,paste0(savedir,i,'.',seedno,'.',fold,'.samplesplit.RDS'))
    print(start- end)
    
    
    
  }
  
}


#machine learning training/testing
dds.matrix = readRDS(paste0(wkdir,'ohs.1000.silencer.norm.counts.RDS'))

for (fold in foldno)   {
  if (file.exists(paste0(savedir,i,'.',seedno,'.',fold,'.dmr.RDS')) == T) {
    #loading previously computed DMRs
    res.df = data.frame(readRDS(paste0(savedir,i,'.',seedno,'.',fold,'.dmr.RDS')),check.names=F)
    
    #selecting samples of interest (ie male prostate cancer/controls)
    targ.samples = targets[[i]]
    targ.samples = targ.samples[targ.samples$GRP_Id %in% colnames(dds.matrix),]
    targ.samples$group = ifelse(targ.samples$Cancer == 'Control','Control','Cancer')
    
    #loading train/test fold partition
    train.set = readRDS(paste0(savedir,i,'.',seedno,'.',fold,'.samplesplit.RDS'))
    test.set = targ.samples[!targ.samples$GRP_Id %in% train.set$GRP_Id,]
    
    
    #selecting features for ML models
    for (dir in c('abs')){ #abs regions used in study
      predir=paste0(savedir,'predictions.',dir,'/')
      dir.create(predir, recursive = T)
      
      if(file.exists(paste0(predir,i,'.',seedno,'.',fold,'.predictions.RDS')) == F) {
        
        #ordering and selecting top DMRs
        if (dir == 'hyper'){
          res.df.sig= res.df[which(res.df$log2FoldChange > 0.25  & res.df$baseMean  > 1),]
          
        } else {
          res.df.sig= res.df[which(abs(res.df$log2FoldChange) > 0.25  & res.df$baseMean  > 1),]
          
        }
        #ordering features across train set fold p-values
        res.df.sig = res.df.sig[order(res.df.sig$pvalue),]

        #ensuring groups are correct
        train.set$group = ifelse(train.set$Cancer=='Control','Control','Cancer')
        test.set$group = ifelse(test.set$Cancer=='Control','Control','Cancer')
        
        log2.filt= data.frame(t(log2(dds.matrix[res.df.sig$window,unique(targ.samples$GRP_Id)]+1)),check.names=F) 
        matrix.list= list(log2.std = log2.filt)
        
        complete.res.base = NULL
        complete.weight.base= NULL
        for (mat in names(matrix.list)) {
          targ.matrix1 = matrix.list[[mat]]
          targ.features = res.df.sig$window
          feature.sizes= c(seq(10,500,10)) #benchmarking across differing number of top ranking regions
          feature.sizes = feature.sizes[feature.sizes <= length(targ.features) ] #if number of DMRs with logFC > 0.25 is lower than the number of features interrogated
          
          for (fl in feature.sizes) {
            res.df.all = predictive.models.glm(targ.matrix1,
                                           unique(targ.samples[,c('GRP_Id','group','Cancer')]),
                                           train.set,
                                           test.set,
                                           targ.features = c(targ.features[1:fl]),
                                           feats = feature.sizes,
                                           feature.type = paste0(marker),
                                           stat.test = paste0('deseq.',mat))
            #test fold performance
            perf.df = res.df.all[[1]] 
            perf.df$direction = dir

            #feature weightings
            perf.df = res.df.all[[2]] 
            perf.df$direction = dir

            complete.res.base= rbind(complete.res.base,perf.df)
            complete.weight.base= rbind(complete.weight.base,perf.df)
            
          }
          
          
          
          
        }
        
        saveRDS(complete.weight.base,paste0(predir,i,'.',seedno,'.',fold,'.feature.weights.RDS'))
        saveRDS(complete.res.base,paste0(predir,i,'.',seedno,'.',fold,'.predictions.RDS'))

      }
    }
  }
}



########
library(ggplot2)
library(cutpointr)
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
autosome_filt = function(medips.count_df) {
  tmp =gsub(':.*','',medips.count_df)
  return_df = medips.count_df[!tmp %in% c('chrY','chrX','chrM')]
  return(return_df)
}

options(expressions = 5e5)
library(DESeq2)
library("BiocParallel")
ncores = 10
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

####selecting top performing parameter#####
gleason.score = readRDS('/.mounts/labs/awadallalab/private/ncheng/cfmedip_data/cptp_samples/pathology_records/prostate/gleason.score.RDS')
discovery.set = readRDS('/.mounts/labs/awadallalab/private/ncheng/manuscripts/early.cancer.risk/discovery.set2.samples.RDS')
validation.set = readRDS('/.mounts/labs/awadallalab/private/ncheng/manuscripts/early.cancer.risk/validation.set2.samples.RDS')
all.sample.info=readRDS('/.mounts/labs/awadallalab/private/ncheng/manuscripts/early.cancer.risk/combined.set2.samples.RDS')

#setting wkdir
wkdir='/.mounts/labs/awadallalab/private/ncheng/cfmedip_data/cptp_samples/fragmentation/aix13.updated1/all.sample.fragmentation/breast.cancer.cv/regulatory.regions.v10/'
savedir=paste0(wkdir,'validation.breast.test/')
savedir1 = paste0(savedir,'/genhancer.validation/')
dir.create(savedir1)
setwd(savedir)
matrixdir='/.mounts/labs/awadallalab/private/ncheng/cfmedip_data/cptp_samples/fragmentation/aix13.updated1/all.sample.fragmentation/regulatory.counts/'
####deseq dmrcalling#####
#dmr calling with 1000
dmrcalling = T
matrices = c('Female.1000','All.1000','Female.400','All.400')

#male
merged.df.filt.targ = discovery.set[discovery.set$Sex == 'Male',]
if (dmrcalling == T) {
  #dds = readRDS(paste0(matrixdir,'aix13.combined.',1000,'.','genhancer','.dds.RDS'))
  #dds = readRDS('/.mounts/labs/awadallalab/private/ncheng/cfmedip_data/cptp_samples/fragmentation/aix13.updated1/all.sample.fragmentation/regulatory.counts/aix13.discovery.1000.genhancer.dds.RDS')
  #dds = readRDS('/.mounts/labs/awadallalab/private/ncheng/cfmedip_data/cptp_samples/fragmentation/aix13.updated1/all.sample.fragmentation//AIX13.dv.Allall.samples.dds.inserts.1000.all.300.q20.RDS')
  #dds = readRDS('/.mounts/labs/awadallalab/private/ncheng/cfmedip_data/cptp_samples/fragmentation/aix13.updated1/all.sample.fragmentation/AIX13.discovery.Female.dds.inserts.400.all.300.q20.RDS')
  #dds = readRDS('/.mounts/labs/awadallalab/private/ncheng/cfmedip_data/cptp_samples/fragmentation/aix13.updated1/all.sample.fragmentation/regulatory.counts/aix13.all.dv.1000.silencer.dds.RDS')
  #dds = readRDS('/.mounts/labs/awadallalab/private/ncheng/cfmedip_data/cptp_samples/fragmentation/aix13.updated1/all.sample.fragmentation/regulatory.counts/aix13.maleall.1000.silencer.dds3.RDS')
  #matrix/proportions comparison#
  dds = readRDS('/.mounts/labs/awadallalab/private/ncheng/cfmedip_data/cptp_samples/fragmentation/aix13.updated1/all.sample.fragmentation/regulatory.counts/aix13.combined.1000.silencer.dds.RDS')
  
  
  colData(dds)$total_reads= colData(dds)$total/1000000
  
  colData(dds)$total_reads= colData(dds)$total/1000000
  windows = rownames(dds)
  windows = windows[!grepl('chrX|chrY',windows)]
  
  
  targ.samples = all.sample.info#[all.sample.info$Sex == 'Male',]
  targ.samples = unique(targ.samples[targ.samples$GRP_Id %in% colnames(dds),])
  #dds = estimateSizeFactors(dds[windows,targ.samples$GRP_Id])
  dds$condition = dds$group
  
  
  
  colData(dds)$total.fragment.count.m = colData(dds)$total.fragment.count/1000000
  colData(dds)$filler = factor(colData(dds)$filler,levels= c('MFiller','UFiller'))
  colData(dds)$group = factor(ifelse(colData(dds)$Cancer =='Control','Control','Cancer'),levels = c('Control','Cancer'))
  windows = rownames(dds)
  windows = windows[!grepl('chrX|chrY',windows)]
  ##breast cancer vs controls##
  merged.df.filt.targ = discovery.set[discovery.set$GRP_Id %in% colnames(dds),]
  merged.df.filt.targ = merged.df.filt.targ[merged.df.filt.targ$Sex == 'Male',]
  dds.filt = dds[windows,merged.df.filt.targ$GRP_Id]
  mm = model.matrix(~ THALIANA_BETA +  total.fragment.count.m + filler + group, colData(dds.filt)) # + filler:nested.batch +
  
  dds.filt = dds.filt[rowSums(counts(dds.filt)) > 0,]
  dds.filt$condition = dds.filt$group
  
  
  ddssva <- DESeq(dds.filt,full = mm,parallel=T) #differential methylation analysis
  res.df = results(ddssva,contrast = list('groupCancer')) #generating results table
  
  #saveRDS(res.df,paste0(savedir,'discovery.breast.nonage.adjusted.dmrs.genhancers.RDS'))
  #saveRDS(res.df,paste0(savedir1,'discovery.breast.femalenorm.ntotalTB.silencers.RDS'))
  #saveRDS(res.df,paste0(savedir1,'discovery.prostate.allnorm.ntotalTB.silencers.dvnorm.RDS'))
  savedir1 = paste0(savedir,'/genhancer.validation.prevmatrix/')
  #saveRDS(res.df,paste0(savedir1,'discovery.prostate.allnorm.ntotalTB.silencers.malenorm.prevmatrix.RDS'))
  saveRDS(res.df,paste0(savedir1,'discovery.prostate.allnorm.ntotalTB.silencers.malenorm.prevmatrix.ageadjusted.RDS'))
  
  
}

####setting tuning choices ####
matrices = c('Female.1000','All.1000','Female.400','All.400')
sf = c('AllChr.before','AutoChr.before')
sf = c('AllChr.before','AutoChr.before','SF.before')



#####silencer only ml####
results.df.all = NULL
feature = 'dmr'
#dvnorm wk2
#dnorm wk3
for (s in sf[2]) {
  for (m in matrices[2]) {
    savedir1 = paste0(savedir,'/genhancer.validation/')
    
    print(m)
    #new.split
    #res.df = readRDS(paste0(savedir1,'discovery.allbg.Breast.',m,'.',s,'.dmrs.RDS'))
    #res.df = readRDS(paste0(savedir1,'discovery.breast.allnorm.ntotalTB.genhancers.RDS'))
    res.df = readRDS(paste0(savedir1,'discovery.prostate.allnorm.ntotalTB.silencers.RDS'))
    
    res.df$window =rownames(res.df)
    directions = c('abs','hyper')[1]#[2]#[1]#[foldno]#[1]#[2]
    print('c1')
    if (m == 'All.1000') {
      #sample.matrix = readRDS('/.mounts/labs/awadallalab/private/ncheng/cfmedip_data/cptp_samples/fragmentation/aix13.updated1/all.sample.fragmentation//AIX13.dv.Allall.samples.deseq.normcounts.inserts.1000.all.300.q20.RDS')
      #sample.matrix = readRDS('/.mounts/labs/awadallalab/private/ncheng/cfmedip_data/cptp_samples/fragmentation/aix13.updated1/all.sample.fragmentation/regulatory.counts/aix13.all.dv.1000.genhancer.norm.counts.RDS')
      
      dds = readRDS(paste0('/.mounts/labs/awadallalab/private/ncheng/cfmedip_data/cptp_samples/fragmentation/aix13.updated1/all.sample.fragmentation/regulatory.counts/aix13.all.dv.1000.silencer.dds.RDS'))
      #dds.matrix = readRDS(paste0(wkdir,'aix13.combined.',inserts,'.',marker,'.norm.counts.RDS'))
      
      colData(dds)$total_reads= colData(dds)$total/1000000
      windows = rownames(dds)
      windows = windows[!grepl('chrX|chrY',windows)]
      
      
      targ.samples = all.sample.info#[discovery.set$Sex == 'Female',]
      targ.samples = unique(targ.samples[targ.samples$GRP_Id %in% colnames(dds),])
      dds = estimateSizeFactors(dds[windows,targ.samples$GRP_Id])
      dds$condition = dds$group
      sample.matrix = counts(dds, normalized = T)
      
      
      #sf.calc = F
      
      
      
      
      
    }
    
    fc.cutoff=0.25
    results.df.all = NULL
    print('c2')
    for (fc in fc.cutoff) {
      for (d in directions) {
        
        print(d)
        if(feature == 'dmr'){
          res.df.targ = res.df#[rownames(res.df) %in% background,]
          if (grepl('abs',d) == T) {
            res.df.targ = res.df.targ[which(abs(res.df.targ$log2FoldChange) > fc),]# & res.df.targ$baseMean > 1
          } else if (grepl('hyper',d) == T){
            res.df.targ = res.df.targ[which(res.df.targ$log2FoldChange> fc   ),]#& res.df.targ$baseMean > 1
          }
          print('c3')
          res.df.targ = res.df.targ[order(res.df.targ$pvalue),]
          res.df.targ = res.df.targ[rownames(res.df.targ) %in% rownames(sample.matrix),]
          res.df.targ = res.df.targ[order(res.df.targ$pvalue),]
          
        } else if (feature == 'coef') {
          res.df.targ = res.df.targ#[rownames(res.df) %in% background,]
          if (grepl('abs',d) == T) {
            res.df.targ = res.df.targ[which(abs(res.df.targ$coef) > fc ),]
          } else if (grepl('hyper',d) == T){
            res.df.targ = res.df.targ[which(res.df.targ$coef> fc),]
          }          
          print('c3')
          res.df.targ = res.df.targ[rownames(res.df.targ) %in% rownames(sample.matrix),]
          res.df.targ = res.df.targ[order(-abs(res.df.targ$coef)),]
          
        } else if (feature == 'freq') {
          res.df.targ = res.df.targ#[rownames(res.df) %in% background,]
          res.df.targ1 = res.df.targ[which(abs(res.df.targ$freq) >= 30 ),]
          
          print('c3')
          res.df.targ = res.df.targ[rownames(res.df.targ) %in% rownames(sample.matrix),]
          res.df.targ = res.df.targ[order(-res.df.targ$freq),]
          
        }
        
        #sample.matrix = readRDS('')
        targ.samples = all.sample.info[all.sample.info$Sex == 'Male',]
        targ.samples = targ.samples[targ.samples$GRP_Id %in% colnames(sample.matrix),]
        
        mat = 'log2.std'
        if (mat == 'log2.std') {
          #std.matrix= data.frame(t(log2(sample.matrix[rownames(res.df.targ),targ.samples]+1)),check.names=F) 
          std.matrix= data.frame(t((sample.matrix[rownames(res.df.targ),unique(targ.samples$GRP_Id)]+1)),check.names=F) 
          #std.matrix = data.frame(t(sample.matrix[rownames(res.df.targ),unique(targ.samples$GRP_Id)]),check.names=F)
          #plot.samples = colnames()
          #std.matrix= data.frame(apply(sample.matrix[rownames(res.df.targ),unique(targ.samples$GRP_Id)],1,scale),check.names=F)
          
          rownames(std.matrix)= targ.samples$GRP_Id#colnames(sample.matrix)
          colnames(std.matrix) =rownames(res.df.targ)
          std.matrix = std.matrix[targ.samples$GRP_Id,]
          #std.matrix
        } else {
          std.matrix = data.frame(t(sample.matrix[rownames(res.df.targ),unique(targ.samples$GRP_Id)]),check.names=F)
          rownames(std.matrix)= unique(targ.samples$GRP_Id)
          colnames(std.matrix) =rownames(res.df.targ)
          
        }
        
        print('c4')
        
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
        
        
        feature.weights.all = NULL
        targ.features = rownames(res.df.targ)
        
        train.set= discovery.set[discovery.set$Sex == 'Male',]
        test.set = validation.set[validation.set$GRP_Id %in% colnames(sample.matrix),]
        test.set = test.set[test.set$Sex == 'Male',]
        results.df = NULL
        #base
        #n.features = c(seq(10,100,10),seq(125,300,25))
        #n.features = c(seq(25,200,25))#,seq(250,400,50))
        n.features = seq(10,150,10)#c(80,90,100)
        n.features = n.features[n.features <= nrow(res.df.targ)]
        #n.features = c(25,50)
        print('c5')
        #n.features = 75#80#c(seq(50,90,10),nrow(res.df.targ))#,seq(250,400,50))
        n.features=seq(90,140,10)
        n.features=100
        #seq(100,130,10))
        score.cutoff.breast = 0.283
        for (filt in c('base')) {
          print(filt)
          for (fl in n.features){
            print(fl)
            print(fl)
            start= Sys.time()
            f = min(fl, length(targ.features))
            
            #all features as features
            
            
            targ.features.df = res.df.targ#
            #targ.features.df = targ.features.df[order(targ.features.df$pvalue),]
            targ.features.df$window = rownames(targ.features.df) 
            targ.features.base = targ.features.df$window[1:f]
            targ.matrix.base = std.matrix[,targ.features.base]
            targ.matrix.base$GRP_Id = rownames(targ.matrix.base)
            covariates = c('filler','none')[2]#,'filler.total','none')#[foldno]
            for (covar in covariates) {
              if (covar == 'filler') {
                targ.matrix=  merge(targ.matrix.base,all.sample.info[,c('GRP_Id','group','filler')],by='GRP_Id')
                targ.features1 =c('filler',targ.features.base)
              } else if (covar == 'filler.total') {
                targ.matrix=  merge(targ.matrix.base,all.sample.info[,c('GRP_Id','group','filler','total')],by='GRP_Id')
                targ.matrix$total = targ.matrix$total/1000000
                targ.features1 =c('filler','total',targ.features.base)            
              } else if (covar == 'TB') {
                targ.matrix=  merge(targ.matrix.base,all.sample.info[,c('GRP_Id','group','THALIANA_BETA')],by='GRP_Id')
                targ.features1 =c('THALIANA_BETA',targ.features.base)    
              } else if (covar == 'none'){
                targ.matrix=  merge(targ.matrix.base,all.sample.info[,c('GRP_Id','group')],by='GRP_Id')
                targ.features1 =c(targ.features.base)  
              } else if (covar == 'TB.filler') {
                targ.matrix=  merge(targ.matrix.base,all.sample.info[,c('GRP_Id','group','filler','THALIANA_BETA')],by='GRP_Id')
                targ.features1 =c('THALIANA_BETA','filler',targ.features.base)    
              } 
              rownames(targ.matrix) = targ.matrix$GRP_Id
              
              
              # results.df.tmp = NULL
              tmp.res = mclapply(1:50, function(seednum) {
                set.seed(seednum)
                res.df.all = predictive.models.glm(targ.matrix,
                                                   unique(all.sample.info[,c('GRP_Id','group','Cancer')]),
                                                   train.set,
                                                   test.set,
                                                   targ.features =targ.features1,
                                                   feats = fl,
                                                   feature.type = 'silencer',
                                                   stat.test ='silencer')
                perf.df = res.df.all[[1]] 
                # perf.df = perf.df[!]
                #results.df.tmp = rbind(results.df.tmp,perf.df)
                return(perf.df)
              },mc.cores=5)
              
              tmp.res = tmp.res[sapply(tmp.res,function(x) grepl('Error',x)[1] == F)]
              results.df.tmp = do.call('rbind',tmp.res)
              ##
              results.df.tmp = results.df.tmp[grepl('AIX.*',results.df.tmp$GRP_Id),]
              
              library(plyr)
              results.df.tmp$methylation_score=as.numeric(results.df.tmp$methylation_score)
              #results.df.tmp$auroc=as.numeric(results.df.tmp$auroc)
              
              combined.collapse =ddply(results.df.tmp[,c('GRP_Id','reported','methylation_score','model')], c('GRP_Id','reported','model'),numcolwise(mean,na.rm=T))
              combined.collapse = combined.collapse[combined.collapse$model == 'logreg.old.alpha1',]
              combined.collapse$auroc = auc_calc(combined.collapse)
              combined.collapse$predictions = ifelse(combined.collapse$methylation_score >= score.cutoff.breast,'Cancer','Control')
              combined.collapse$sf = s
              combined.collapse$matrix=m
              combined.collapse$direction= d
              combined.collapse$features = fl
              combined.collapse$filter = filt
              combined.collapse$fc = fc
              combined.collapse$covariate=covar
              
              results.df.all = rbind(results.df.all,combined.collapse)
              
              #
              
            }
            
            
            
            
          }
          
        }
        
        #results.df.all$fc = fc
        savedir2=paste0(savedir1,'/fc.assess.ntotal.genhancer.v3/')
        dir.create(savedir2,recursive = T)
        #saveRDS(results.df.all,paste0(savedir2,'discovery.Breast.',filt,'.',fc,'.',d,'.',s,'.',m,'.performance.covarbase.RDS'))
        #saveRDS(results.df.all[results.df.all$model =='logreg.old.alpha1',],paste0(savedir2,'discovery.Breast.validation.genhancer.final.performance.covarbase.RDS'))
        
        
        
        #
        
        
      }
      
      
    }
    
    
    
    
    
    
    
  }
  
}



wkdir='/.mounts/labs/awadallalab/private/ncheng/cfmedip_data/cptp_samples/fragmentation/aix13.updated1/all.sample.fragmentation/breast.cancer.cv/regulatory.regions.v10.newsplit5/'
savedir=paste0(wkdir,'validation.breast.test/')
figdir=paste0(savedir,'figures.prad.validation.feat100.bm0.final/')
dir.create(predir,recursive = T)
dir.create(figdir,recursive = T)
gleason.score = readRDS('/.mounts/labs/awadallalab/private/ncheng/cfmedip_data/cptp_samples/pathology_records/prostate/gleason.score.RDS')
discovery.set = readRDS('/.mounts/labs/awadallalab/private/ncheng/manuscripts/early.cancer.risk/discovery.set2.samples.RDS')
validation.set = readRDS('/.mounts/labs/awadallalab/private/ncheng/manuscripts/early.cancer.risk/validation.set2.samples.RDS')
all.sample.info=readRDS('/.mounts/labs/awadallalab/private/ncheng/manuscripts/early.cancer.risk/combined.set2.samples.RDS')
sample.info = all.sample.info
results.df.all = readRDS('/.mounts/labs/awadallalab/private/ncheng/cfmedip_data/cptp_samples/fragmentation/aix13.updated1/all.sample.fragmentation/breast.cancer.cv/regulatory.regions.v10.newsplit5/validation.breast.test/validation.prad.scores.RDS')

pred.df.targ = results.df.all[results.df.all$features == 100,];
name = paste0(figdir,'prad.validation');
merged.df.all= validation.set;
gleason.score = gleason.score;
score.cutoff=0.1511;
dx.all = T;
cutpoint.use =T;
gs.plot= T
#saveRDS(pred.df.targ,paste0(name,'.scores.RDS'))


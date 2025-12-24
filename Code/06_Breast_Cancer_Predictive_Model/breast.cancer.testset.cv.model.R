
colData(dds)$total_reads= colData(dds)$total/1000000
windows = rownames(dds)
windows = windows[!grepl('chrX|chrY',windows)]


targ.samples = all.sample.info#[all.sample.info$Sex == 'Female',]
targ.samples = unique(targ.samples[targ.samples$GRP_Id %in% colnames(dds),])
dds = estimateSizeFactors(dds[windows,targ.samples$GRP_Id])
dds$condition = dds$group
sample.matrix = counts(dds, normalized = T)
##

savedir='/.mounts/labs/awadallalab/private/ncheng/cfmedip_data/cptp_samples/fragmentation/aix13.updated1/all.sample.fragmentation/regulatory.counts/'
setwd(savedir)
#genhancer counts.ohs
gleason.score = readRDS('/.mounts/labs/awadallalab/private/ncheng/cfmedip_data/cptp_samples/pathology_records/prostate/gleason.score.RDS')
discovery.set = readRDS('/.mounts/labs/awadallalab/private/ncheng/manuscripts/early.cancer.risk/discovery.set3.samples.RDS')
validation.set = readRDS('/.mounts/labs/awadallalab/private/ncheng/manuscripts/early.cancer.risk/validation.set3.samples.RDS')
all.sample.info=readRDS('/.mounts/labs/awadallalab/private/ncheng/manuscripts/early.cancer.risk/combined.set3.samples.RDS')
#genhancer.matrix = readRDS(paste0('/.mounts/labs/awadallalab/private/ncheng/cfmedip_data/cptp_samples/fragmentation/aix13.updated1/all.sample.fragmentation/regulatory.counts/aix13.all.dv.1000.genhancer.dds.RDS'))
genhancer.dds = readRDS(paste0('/.mounts/labs/awadallalab/private/ncheng/cfmedip_data/cptp_samples/fragmentation/aix13.updated1/all.sample.fragmentation/regulatory.counts/aix13.all.dv.1000.genhancer.dds.RDS'))
genhancer.matrix =counts(genhancer.dds,normalize=F)
windows = rownames(genhancer.matrix)
windows.bed = data.frame(chr=gsub(':.*','',windows),
                         start = gsub('.*:','',gsub('-.*','',windows)),
                         end = gsub('.*-','',windows),
                         window = windows)
write.table(windows.bed,'genhancer.positions.hg38.bed',row.names=F,col.names=F,quote=F,sep='\t')

#dds.matrix = readRDS(paste0(wkdir,'aix13.combined.',inserts,'.',marker,'.norm.counts.RDS'))

#adding octane + burgener#
octane.counts = readRDS('/.mounts/labs/awadallalab/private/ncheng/cfmedip_data/other_data/BRCA_control_cescon/OCTANE_HG38_NonSex_MEDIPS_Counts.rds')
burgener.counts = readRDS('/.mounts/labs/awadallalab/private/ncheng/cfmedip_data/other_data/BRCA_control_cescon/Justin_HG38_NonSex_MEDIPS_Counts.rds')
octane.sample.info = readRDS('/.mounts/labs/awadallalab/private/ncheng/cfmedip_data/other_data/BRCA_control_cescon/OCTANE_Clinical_Info.rds')
burgener.sample.info = readRDS('/.mounts/labs/awadallalab/private/ncheng/cfmedip_data/other_data/BRCA_control_cescon/justin.sample.sex.RDS')

#octane.counts$window = gsub('1-','0-',sub('\\.','-',sub('\\.',':',octane.counts$window)))
#burgener.counts$window = gsub('1-','0-',sub('\\.','-',sub('\\.',':',burgener.counts$window)))

combined.cescon.counts = cbind(octane.counts[,-1],burgener.counts[,-1])
library(stringr)
windows = str_split_fixed(octane.counts$window,pattern = '\\.',n=3)
library(stringr)
windows.bed = data.frame(windows,check.names=F)
windows.bed$window = paste0(windows.bed[,1],':',windows.bed[,2],'-',windows.bed[,3])
#write.table(windows.bed,'medip.300.cescon.positions.hg38.bed',row.names=F,col.names=F,quote=F,sep='\t')


####
window.metadata = read.table('/.mounts/labs/awadallalab/private/ncheng/cfmedip_data/cptp_samples/fragmentation/aix13.updated1/all.sample.fragmentation/regulatory.counts/cescon.300.genhancer.overlap.bed',
                             sep='\t',stringsAsFactors = F,header=F)
window.metadata = window.metadata[,c(4,8)]
colnames(window.metadata) = c('window.genhancer','window.300')
combined.cescon.counts$window.300 = windows.bed$window
combined.cescon.counts.annotated = merge(window.metadata,combined.cescon.counts,by ='window.300')
library(plyr)
combined.cescon.counts.annotated.collapse = ddply(combined.cescon.counts.annotated[,-1],'window.genhancer',numcolwise(sum))
saveRDS(combined.cescon.counts.annotated.collapse,'/.mounts/labs/awadallalab/private/ncheng/cfmedip_data/cptp_samples/fragmentation/aix13.updated1/all.sample.fragmentation/regulatory.counts/octane.burgener.genhancer.inf.rawcounts.RDS')


#combining matrices

colnames(combined.cescon.counts.annotated.collapse)[1] = 'window'
genhancer.matrix = data.frame(genhancer.matrix,check.names=F)
genhancer.matrix$window = rownames(genhancer.matrix)
window.df= data.frame(window=genhancer.matrix$window)
octane.counts.combined=  merge(window.df, combined.cescon.counts.annotated.collapse,by='window',all.x=T)
octane.counts.combined[is.na(octane.counts.combined)] = 0
rownames(octane.counts.combined) = octane.counts.combined$window
octane.counts.combined = octane.counts.combined[genhancer.matrix$window,]
combined.counts = cbind(genhancer.matrix,octane.counts.combined[,-1])
saveRDS(combined.counts,'/.mounts/labs/awadallalab/private/ncheng/cfmedip_data/cptp_samples/fragmentation/aix13.updated1/all.sample.fragmentation/regulatory.counts/aix13.dv.octane.burgener.combined.genhancer.raw.counts.RDS')



#if combining with octane etc
octane.sample.info =data.frame(octane.sample.info)
octane.sample.info$GRP_Id = octane.sample.info$SampleID
octane.sample.info$Cancer = 'Breast'
octane.sample.info$Diagnosis_Time = 'Post-Diagnosis'
octane.sample.info$seq.date = 'OCTANE'
burgener.sample.info = data.frame(burgener.sample.info)
burgener.sample.info$GRP_Id = burgener.sample.info[,4]
burgener.sample.info$Cancer = ifelse(grepl('Norm',burgener.sample.info$GRP_Id ) == T, 'Control','HNSC')
burgener.sample.info$Diagnosis_Time = ifelse(burgener.sample.info$Cancer == 'Control','Control','Post-Diagnosis')
burgener.sample.info$seq.date='burgener'
burgener.sample.info$Sex = ifelse(burgener.sample.info[,2] == 'M','Male','Female')
octane.sample.info$Sex = 'Female'
all.sample.info$seq.date= as.character(all.sample.info$seq.date)
combined.sample.info =rbind(all.sample.info[,c('GRP_Id','Cancer','Diagnosis_Time','seq.date','Sex')],
                            burgener.sample.info[,c('GRP_Id','Cancer','Diagnosis_Time','seq.date','Sex')],
                            octane.sample.info[,c('GRP_Id','Cancer','Diagnosis_Time','seq.date','Sex')])
rownames(combined.sample.info) = combined.sample.info$GRP_Id
saveRDS(combined.sample.info[combined.sample.info$GRP_Id %in% colnames(combined.counts),],paste0('/.mounts/labs/awadallalab/private/ncheng/cfmedip_data/cptp_samples/fragmentation/aix13.updated1/all.sample.fragmentation/regulatory.counts/aix.octane.burg.sample.info.RDS'))

#creating dds object
library(DESeq2)
combined.sample.info = combined.sample.info[combined.sample.info$GRP_Id %in% colnames(combined.counts),]
dds <- DESeqDataSetFromMatrix(countData = combined.counts[,combined.sample.info$GRP_Id],
                              colData = combined.sample.info,
                              design= ~  Cancer ) #can add gender here but we're only looking at female samples currently
saveRDS(dds,'/.mounts/labs/awadallalab/private/ncheng/cfmedip_data/cptp_samples/fragmentation/aix13.updated1/all.sample.fragmentation/regulatory.counts/aix13.dv.octane.burgener.combined.genhancer.dds.RDS')



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
for (s in sf[2]) {
  for (m in matrices[2]) {
    savedir1 = paste0(savedir,'/genhancer.hrsplit.final/')
    
    print(m)
    #new.split
    #res.df = readRDS(paste0(savedir1,'discovery.allbg.Breast.',m,'.',s,'.dmrs.RDS'))
    res.df = readRDS(paste0(savedir1,'discovery.breast.allnorm.ntotalTB.genhancers.RDS'))
    #res.df = readRDS(paste0(savedir1,'discovery.breast.allnorm.base.genhancers.RDS'))
    
    res.df$window =rownames(res.df)
    directions = c('abs','hyper')[2]#[2]#[1]#[foldno]#[1]#[2]
    print('c1')
    {
      #sample.matrix = readRDS('/.mounts/labs/awadallalab/private/ncheng/cfmedip_data/cptp_samples/fragmentation/aix13.updated1/all.sample.fragmentation//AIX13.dv.Allall.samples.deseq.normcounts.inserts.1000.all.300.q20.RDS')
      #sample.matrix = readRDS('/.mounts/labs/awadallalab/private/ncheng/cfmedip_data/cptp_samples/fragmentation/aix13.updated1/all.sample.fragmentation/regulatory.counts/aix13.all.dv.1000.genhancer.norm.counts.RDS')
      
      dds = readRDS(paste0('/.mounts/labs/awadallalab/private/ncheng/cfmedip_data/cptp_samples/fragmentation/aix13.updated1/all.sample.fragmentation/regulatory.counts/aix13.dv.octane.burgener.combined.genhancer.dds.RDS'))
      #dds.matrix = readRDS(paste0(wkdir,'aix13.combined.',inserts,'.',marker,'.norm.counts.RDS'))
      targ.samples.all =readRDS('/.mounts/labs/awadallalab/private/ncheng/cfmedip_data/cptp_samples/fragmentation/aix13.updated1/all.sample.fragmentation/regulatory.counts//aix.octane.burg.sample.info.RDS')
      windows = rownames(dds)
      windows = windows[!grepl('chrX|chrY',windows)]
      
      
      #targ.samples = all.sample.info#[discovery.set$Sex == 'Female',]
      targ.samples = unique(targ.samples.all[targ.samples.all$GRP_Id %in% colnames(dds),])
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
            res.df.targ = res.df.targ[which(abs(res.df.targ$log2FoldChange) > fc & res.df.targ$baseMean > 1),]
          } else if (grepl('hyper',d) == T){
            res.df.targ = res.df.targ[which(res.df.targ$log2FoldChange> fc   & res.df.targ$baseMean > 1),]
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
        targ.samples = targ.samples.all[targ.samples.all$Sex == 'Female',]
        targ.samples = targ.samples[targ.samples$GRP_Id %in% colnames(sample.matrix),]
        
        mat = 'log2.std'
        if (mat == 'log2.std') {
          #std.matrix= data.frame(t(log2(sample.matrix[rownames(res.df.targ),]+1)),check.names=F) 
          std.matrix= data.frame(t(log2(sample.matrix[rownames(res.df.targ),unique(targ.samples$GRP_Id)]+1)),check.names=F) 
          #std.matrix = data.frame(t(sample.matrix[rownames(res.df.targ),unique(targ.samples$GRP_Id)]),check.names=F)
          #plot.samples = colnames()
          #std.matrix= data.frame(t(apply(std.matrix,1,scale)),check.names=F)
          
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
        
        train.set= discovery.set[discovery.set$Sex == 'Female',]
        test.set = targ.samples.all[!targ.samples.all$GRP_Id %in% train.set$GRP_Id,]
        test.set = test.set[test.set$Sex == 'Female',]
        results.df = NULL
        #base
        #n.features = c(seq(10,100,10),seq(125,300,25))
        #n.features = c(seq(25,200,25))#,seq(250,400,50))
        #n.features = c(seq(50,150,10))#c(80,90,100)
        #n.features = c(25,50)
        print('c5')
        #n.features = 75#80#c(seq(50,90,10),nrow(res.df.targ))#,seq(250,400,50))
        n.features=90
        score.cutoff.breast = 0.5
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
              targ.samples$group = ifelse(targ.samples$Cancer == 'Control','Control','Cancer')
              targ.matrix=  merge(targ.matrix.base,targ.samples[,c('GRP_Id','group')],by='GRP_Id')
              targ.features1 =c(targ.features.base)  
              
              rownames(targ.matrix) = targ.matrix$GRP_Id
              
              
              targ.samples.all$group = ifelse(targ.samples.all$Cancer == 'Control','Control','Cancer')
              
              tmp.res = mclapply(1:50, function(seednum) {
                set.seed(seednum)
                
                res.df.all = predictive.models.glm(targ.matrix,
                                                   unique(targ.samples.all[,c('GRP_Id','group','Cancer')]),
                                                   train.set,
                                                   test.set,
                                                   targ.features =targ.features1,
                                                   feats = fl,
                                                   feature.type = 'gs',
                                                   stat.test ='gs')
                perf.df = res.df.all[[1]] 
                # perf.df = perf.df[!]
                #results.df.tmp = rbind(results.df.tmp,perf.df)
                return(perf.df)
              },mc.cores=10)
              
              tmp.res = tmp.res[sapply(tmp.res,function(x) grepl('Error',x)[1] == F)]
              results.df.tmp = do.call('rbind',tmp.res)
              ##
              #results.df.tmp = results.df.tmp[grepl('AIX.*',results.df.tmp$GRP_Id),]
              
              library(plyr)
              results.df.tmp$methylation_score=as.numeric(results.df.tmp$methylation_score)
              
              
              ##
              
              library(plyr)
              combined.collapse =ddply(results.df.tmp[,c('GRP_Id','reported','methylation_score','model')], c('GRP_Id','reported','model'),numcolwise(mean))
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

####computing AUC plot + barplot######
saveRDS(combined.collapse,paste0(figdir,'octane.scores.RDS'))

#octane vs ohs controls
combined.collapse.oct.ohs = combined.collapse[grepl('OCT_',combined.collapse$GRP_Id) == T | grepl('AIX_',combined.collapse$GRP_Id) == T  & combined.collapse$reported =='Control',]

#plot auroc
library(pROC)
figdir='/.mounts/labs/awadallalab/private/ncheng/cfmedip_data/cptp_samples/fragmentation/aix13.updated1/all.sample.fragmentation/breast.cancer.cv/regulatory.regions.v10/validation.breast.test/octane.v/'
dir.create(figdir)
auc.plot.all.s2 = summary.calc(pred.df.targ.collapse[pred.df.targ.collapse$GRP_Id %in% age2$GRP_Id,],merged.df.all.tmp)
auc.plot.all.s2[[2]]$var = 'HR-\nHER2-'
auc.plot.all.s2[[2]]$var.group =  'Subtype'


octane.sample.info = readRDS('/.mounts/labs/awadallalab/private/ncheng/cfmedip_data/other_data/BRCA_control_cescon/OCTANE_Clinical_Info.rds')
octane.sample.info = data.frame(octane.sample.info,check.names=F,stringsAsFactors = F)
burgener.sample.info = readRDS('/.mounts/labs/awadallalab/private/ncheng/cfmedip_data/other_data/BRCA_control_cescon/justin.sample.sex.RDS')
octane.sample.info$GRP_Id = sub('-','_',octane.sample.info[,1])
#octane.sample.info =
octane.auc.plot = function(pred.df.targ,name,merged.df.all) {
  library(ggpubr)
  library(ggplot2)
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
  
  library(cenROC)
  library(plyr)
  library(ggplot2)
  library(ggh4x)
  pred.df.targ.collapse = ddply(pred.df.targ[,c('GRP_Id','reported','methylation_score')],c('GRP_Id','reported'),numcolwise(median))
  pred.df.targ.collapse.annotated = merge(pred.df.targ.collapse,merged.df.all,by='GRP_Id' )
  
  combined.auroc = NULL
  auc.t.combined = NULL
  for (c in c('All')) {
    
    #subtype
    
    dx1 = octane.sample.info[octane.sample.info$'ER STATUS' == 'Positive' & octane.sample.info$'HER-2 STATUS' == 'Positive' |
                               octane.sample.info$'PR STATUS' == 'Positive' & octane.sample.info$'HER-2 STATUS' == 'Positive',]
    dx2 = octane.sample.info[octane.sample.info$'ER STATUS' == 'Positive' & octane.sample.info$'HER-2 STATUS' == 'Negative' |
                               octane.sample.info$'PR STATUS' == 'Positive' & octane.sample.info$'HER-2 STATUS' == 'Negative',]
    dx3 = octane.sample.info[octane.sample.info$'ER STATUS' == 'Negative' & 
                               octane.sample.info$'PR STATUS' == 'Negative' & octane.sample.info$'HER-2 STATUS' == 'Positive',]
    dx4 = octane.sample.info[octane.sample.info$'ER STATUS' == 'Negative' & 
                               octane.sample.info$'PR STATUS' == 'Negative' & octane.sample.info$'HER-2 STATUS' == 'Negative',]
    
    
    library(pROC)
    #dx time 1 year
    
    summary.calc = function(pred.df.targ.collapse,merged.df.all.tmp) {
      if (length(unique(pred.df.targ.collapse$reported)) > 1){
        auc.plot.all = tpr.fpr.calc(pred.df.targ.collapse)
        auc.plot.all$auc = auc_calc(pred.df.targ.collapse,c('Control','Cancer'))
        roc.a= roc(pred.df.targ.collapse,'reported','methylation_score')
        a.confint = ci.auc(roc.a, conf.level=0.95, method=c( "bootstrap"), boot.n = 2000)
        auc.plot.all$auc.lower = a.confint[1]
        auc.plot.all$auc.upper = a.confint[3]
        
        spec.threshold=data.frame(ci.se(roc.a,c(0.95,0.9)))
        
        auc.plot.all$se.spec.95 = spec.threshold[1,2]
        auc.plot.all$se.spec.95.lower = spec.threshold[1,1]
        auc.plot.all$se.spec.95.upper = spec.threshold[1,3]
        
        
        auc.plot.all$se.spec.90 = spec.threshold[2,2]
        auc.plot.all$se.spec.90.lower = spec.threshold[2,1]
        auc.plot.all$se.spec.90.upper = spec.threshold[2,3]
        auc.plot.all$count = nrow(unique(pred.df.targ.collapse[which(pred.df.targ.collapse$reported =='Cancer'),]))
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
                                  se.spec.95=0,
                                  se.spec.95.lower=0,
                                  se.spec.95.upper=0,
                                  se.spec.90=0,
                                  se.spec.90.lower=0,
                                  se.spec.90.upper=0
        )
        
        return(list(NULL,auc.plot.all))
      }
      
    }
    merged.df.all.tmp = merged.df.all
    auc.plot.all = summary.calc(pred.df.targ.collapse,merged.df.all.tmp)
    auc.plot.all[[2]]$var = 'All'
    
    auc.plot.all[[2]]$var.group = 'Subtype'
    
    auc.plot.all.dx1 = summary.calc(pred.df.targ.collapse[pred.df.targ.collapse$GRP_Id %in% dx1$SampleID| pred.df.targ.collapse$reported == 'Control',],merged.df.all.tmp)
    auc.plot.all.dx1[[2]]$var = 'HR+\nHER2+'
    
    auc.plot.all.dx1[[2]]$var.group = 'Subtype'
    
    
    auc.plot.all.dx2 = summary.calc(pred.df.targ.collapse[pred.df.targ.collapse$GRP_Id %in% dx2$SampleID| pred.df.targ.collapse$reported == 'Control',],merged.df.all.tmp)
    auc.plot.all.dx2[[2]]$var = 'HR+\nHER2-'
    
    auc.plot.all.dx2[[2]]$var.group = 'Subtype'
    
    auc.plot.all.dx3 = summary.calc(pred.df.targ.collapse[pred.df.targ.collapse$GRP_Id %in% dx3$SampleID | pred.df.targ.collapse$reported == 'Control',],merged.df.all.tmp)
    auc.plot.all.dx3[[2]]$var = 'HR-\nHER2+'
    
    auc.plot.all.dx3[[2]]$var.group = 'Subtype'
    
    auc.plot.all.dx4 = summary.calc(pred.df.targ.collapse[pred.df.targ.collapse$GRP_Id %in% dx4$SampleID| pred.df.targ.collapse$reported == 'Control',],merged.df.all.tmp)
    
    auc.plot.all.dx4[[2]]$var = 'HR-\nHER2-'
    auc.plot.all.dx4[[2]]$var.group = 'Subtype'
    
    auc.res2 = rbind(auc.plot.all[[2]][1,],
                     auc.plot.all.dx1[[2]][1,],
                     auc.plot.all.dx2[[2]][1,],
                     auc.plot.all.dx3[[2]][1,],
                     auc.plot.all.dx4[[2]][1,])
    auc.res2$title = 'Subtype'
    combined.auroc = rbind(combined.auroc,auc.res2)
    
    
    #annotating others
    subtype_colors =  c('HR+\nHER2+' = '#8FD5A6','HR+\nHER2-' = '#BF1363',"HR-\nHER2+" = '#E09F3E','HR-\nHER2-'='#545E75','NR' ='grey','Control' = 'black')
    
    age_colors = c('All' = "#7f7f7f",'35-45' = '#C5C392','45-55'= '#FFB20F','55-65'='#FF4B3E', '65-75'='#972D07')
    stage_colors = c('All' = "#7f7f7f", 'I' = "#3B0F70FF",'II' = "#8C2981FF",'III/IV' = "#DE4968FF")
    combined.auc.plot.all = rbind(auc.plot.all[[2]],
                                  auc.plot.all.dx1[[2]],
                                  auc.plot.all.dx2[[2]],
                                  auc.plot.all.dx3[[2]],
                                  auc.plot.all.dx4[[2]])
    combined.auc.plot.all$title = 'Subtype'
    
    plot1 = ggplot(combined.auc.plot.all,aes(x = fpr, y = tpr, col = var)) + geom_line(linewidth=1) +
      scale_color_manual(values = c(subtype_colors))+ #+ ggtitle(title) +
      theme_bw()+ 
      #facet_grid2(.~title)+
      theme(text = element_text(size=8),
            axis.text=element_text(size=8),
            axis.title=element_text(size=8,face="bold"),
            legend.position = "none")+ guides(colour = guide_legend(override.aes = list(size=4),nrow=4)) + xlab('False Positive Rate') + ylab('Sensitivity')
    
    png(paste0(name,'.',c,'.auroc.subtype.png'),height = 1000,width=1000,res = 300)
    print(plot1)
    dev.off()
    
    #boxplot of score
    octane.hnsc = burgener.sample.info[burgener.sample.info[,3] == 'HN','id']
    octane.healthy = burgener.sample.info[burgener.sample.info[,3] == 'Norm','id']
    
    pred.df.targ.collapse.annotated = merge(pred.df.targ.collapse,merged.df.all,by='GRP_Id' )
    pred.df.targ.collapse.annotated = merge(pred.df.targ.collapse.annotated, octane.sample.info,by.x='GRP_Id',by.y = 'SampleID',all.x=T)
    pred.df.targ.collapse.annotated$subtype = ifelse(pred.df.targ.collapse.annotated$GRP_Id %in% dx1$SampleID,'HR+\nHER2+',
                                                     ifelse(pred.df.targ.collapse.annotated$GRP_Id %in% dx2$SampleID,'HR+\nHER2-',
                                                            ifelse(pred.df.targ.collapse.annotated$GRP_Id %in% dx3$SampleID,'HR-\nHER2+',
                                                                   ifelse(pred.df.targ.collapse.annotated$GRP_Id %in% dx4$SampleID,'HR-\nHER2-',
                                                                          ifelse(pred.df.targ.collapse.annotated$GRP_Id %in% octane.hnsc,'HNSC',
                                                                                 ifelse(pred.df.targ.collapse.annotated$GRP_Id %in% octane.healthy,'Healthy\nControl','Other'))))))
    
    
    
    plot1 = ggplot(pred.df.targ.collapse.annotated,aes(x = subtype, y = methylation_score, col = subtype)) + 
      geom_boxplot(outlier.shape=NA)+geom_jitter(width=0.1)+
      scale_color_manual(values = c(subtype_colors))+ #+ ggtitle(title) +
      theme_bw()+ 
      scale_y_continuous(limits=c(0,1))+
      theme(text = element_text(size=8),
            axis.text=element_text(size=8, face = "bold"),
            axis.title=element_text(size=8,face="bold"),
            legend.position = "none")+ guides(colour = guide_legend(override.aes = list(size=4),nrow=4)) + xlab('Time to Diagnosis') + ylab('Methylation Score')
    
    png(paste0(name,'.dx.methscore.',c,'.png'),height = 900,width=1100,res = 300)
    print(plot1)
    dev.off()
    
    pred.df.targ.collapse.annotated$Cancer.group =ifelse(pred.df.targ.collapse.annotated$subtype == 'HNSC','HNSC',
                                                         ifelse(pred.df.targ.collapse.annotated$subtype %in% c('Healthy\nControl','Control','Other'),'Healthy\nControl',
                                                                'Breast\nCancer'))
    octane.group.labs= c('Healthy\nControl' = 'grey','HNSC' = '#06D6A0','Breast\nCancer' = '#DB2955')
    pred.df.targ.collapse.annotated$Cancer.group = factor(pred.df.targ.collapse.annotated$Cancer.group,levels = names(octane.group.labs))
    
    non.cancer.vars = c('HNSC','Healthy\nControl')
    non.cancer.vars = non.cancer.vars[non.cancer.vars %in% as.character(pred.df.targ.collapse.annotated$Cancer.group)]
    my_comparisons=list()
    for (i in non.cancer.vars) {
      my_comparisons[[length(my_comparisons)+1]] = c('Breast\nCancer',i)
    }
    
    plot1 = ggplot(pred.df.targ.collapse.annotated,aes(x = Cancer.group, y = methylation_score, col = Cancer.group)) + geom_boxplot(outlier.shape=NA)+geom_jitter(width=0.1)+
      scale_color_manual(values = octane.group.labs)+ #+ ggtitle(title) +
      theme_bw()+ 
      scale_y_continuous(limits=c(0,1.2),breaks = seq(0,1,0.25))+
      theme(text = element_text(size=8),
            axis.text=element_text(size=8, face = "bold"),
            axis.title.y=element_text(size=8,face="bold"),
            axis.title.x=element_blank(),
            legend.position = "none")+ 
      guides(colour = guide_legend(override.aes = list(size=4),nrow=4)) + 
      xlab('Time to Diagnosis') + ylab('Methylation Score')+
      stat_compare_means(comparisons = my_comparisons,label.y=c(1,1.1,1.2,1.3,1.4,1.5)[1:length(my_comparisons)],size=3,method='wilcox.test')
    
    
    png(paste0(name,'.group.methscore.',c,'.png'),height = 1000,width=800,res = 300)
    print(plot1)
    dev.off()
    
    
    
    #age
    
    dx1 = octane.sample.info[octane.sample.info$'AGE' < 50 ,]
    dx2 = octane.sample.info[octane.sample.info$'AGE' >= 50 ,]
    
    library(pROC)
    #dx time 1 year
    summary.calc = function(pred.df.targ.collapse,merged.df.all.tmp) {
      if (length(unique(pred.df.targ.collapse$reported)) > 1){
        auc.plot.all = tpr.fpr.calc(pred.df.targ.collapse)
        auc.plot.all$auc = auc_calc(pred.df.targ.collapse,c('Control','Cancer'))
        roc.a= pROC::roc(pred.df.targ.collapse,'reported','methylation_score')
        a.confint = ci.auc(roc.a, conf.level=0.95, method=c( "bootstrap"), boot.n = 2000)
        auc.plot.all$auc.lower = a.confint[1]
        auc.plot.all$auc.upper = a.confint[3]
        
        spec.threshold=data.frame(ci.se(roc.a,c(0.95,0.9)))
        
        auc.plot.all$se.spec.95 = spec.threshold[1,2]
        auc.plot.all$se.spec.95.lower = spec.threshold[1,1]
        auc.plot.all$se.spec.95.upper = spec.threshold[1,3]
        
        
        auc.plot.all$se.spec.90 = spec.threshold[2,2]
        auc.plot.all$se.spec.90.lower = spec.threshold[2,1]
        auc.plot.all$se.spec.90.upper = spec.threshold[2,3]
        auc.plot.all$count = nrow(unique(pred.df.targ.collapse[pred.df.targ.collapse$reported =='Cancer',]))
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
                                  se.spec.95=0,
                                  se.spec.95.lower=0,
                                  se.spec.95.upper=0,
                                  se.spec.90=0,
                                  se.spec.90.lower=0,
                                  se.spec.90.upper=0
        )
        
        return(list(NULL,auc.plot.all))
      }
      
    }
    merged.df.all.tmp = merged.df.all
    auc.plot.all = summary.calc(pred.df.targ.collapse,merged.df.all.tmp)
    auc.plot.all[[2]]$var = 'All'
    
    auc.plot.all[[2]]$var.group ='Age at Diagnosis'
    
    auc.plot.all.dx1 = summary.calc(pred.df.targ.collapse[pred.df.targ.collapse$GRP_Id %in% dx1$SampleID| pred.df.targ.collapse$reported == 'Control',],merged.df.all.tmp)
    auc.plot.all.dx1[[2]]$var = '< 50'
    
    auc.plot.all.dx1[[2]]$var.group = 'Age at Diagnosis'
    
    
    auc.plot.all.dx2 = summary.calc(pred.df.targ.collapse[pred.df.targ.collapse$GRP_Id %in% dx2$SampleID| pred.df.targ.collapse$reported == 'Control',],merged.df.all.tmp)
    auc.plot.all.dx2[[2]]$var = '> 50'
    
    auc.plot.all.dx2[[2]]$var.group = 'Age at Diagnosis'
    
    auc.res2 = rbind(#auc.plot.all[[2]][1,],
      auc.plot.all.dx1[[2]][1,],
      auc.plot.all.dx2[[2]][1,])
    auc.res2$title = 'Age'
    combined.auroc = rbind(combined.auroc,auc.res2)
    
    
    #annotating others
    subtype_colors =  c('HR+\nHER2+' = '#8FD5A6','HR+\nHER2-' = '#BF1363',"HR-\nHER2+" = '#E09F3E','HR-\nHER2-'='#545E75','NR' ='grey','Control' = 'black')
    age50_colors = c('All'="#7f7f7f",'< 50' = '#C5C392','> 50' = '#972D07')
    age_colors = c('All' = "#7f7f7f",'35-45' = '#C5C392','45-55'= '#FFB20F','55-65'='#FF4B3E', '65-75'='#972D07')
    stage_colors = c('All' = "#7f7f7f", 'I' = "#3B0F70FF",'II' = "#8C2981FF",'III/IV' = "#DE4968FF")
    combined.auc.plot.all = rbind(#auc.plot.all[[2]],
      auc.plot.all.dx1[[2]],
      auc.plot.all.dx2[[2]])
    combined.auc.plot.all$title = 'Age'
    
    plot1 = ggplot(combined.auc.plot.all,aes(x = fpr, y = tpr, col = var)) + geom_line(linewidth=1) +
      scale_color_manual(values = c(age50_colors))+ #+ ggtitle(title) +
      theme_bw()+ 
      #facet_grid2(.~title)+
      theme(text = element_text(size=8),
            axis.text=element_text(size=8),
            axis.title=element_text(size=8,face="bold"),
            legend.position = "none")+ guides(colour = guide_legend(override.aes = list(size=4),nrow=4)) + xlab('False Positive Rate') + ylab('Sensitivity')
    
    png(paste0(name,'.',c,'.auroc.age.png'),height = 1000,width=1000,res = 300)
    print(plot1)
    dev.off()
    
    #boxplot of score
    
    pred.df.targ.collapse.annotated$age = ifelse(pred.df.targ.collapse.annotated$GRP_Id %in% dx1$SampleID,'< 50',
                                                 ifelse(pred.df.targ.collapse.annotated$GRP_Id %in% dx2$SampleID,'> 50','Other'))
    
    
    
    
    
    ###summary plot
    #sensitivity
    combined.auroc = combined.auroc[!duplicated(combined.auroc$var),]
    combined.auroc$var.group = ifelse(combined.auroc$var == 'All','Age at Diagnosis',combined.auroc$var.group)
    combined.auroc$var.group = factor(combined.auroc$var.group,levels = c('Subtype','Age at Diagnosis'))
    combined.auroc$var.plot = paste0(combined.auroc$var,'\n(n=',combined.auroc$count,')')
    
    plot1 = ggplot(combined.auroc[combined.auroc$var.group =='Subtype',],aes(x = var.plot, y = auc, fill = var,ymin=auc.lower, ymax=auc.upper)) + 
      geom_bar(stat='identity')+
      scale_fill_manual(values = c(age50_colors,subtype_colors))+ #+ ggtitle(title) +
      theme_bw()+ 
      geom_errorbar(width=.2,
                    position=position_dodge(.9)) +
      
      #facet_grid2(. ~ var.group,scales='free', independent = "x")+
      #force_panelsizes(cols = c(0.5, 0.5)) +
      scale_y_continuous(limits=c(0,1))+
      theme(text = element_text(size=8),
            axis.text=element_text(size=8),
            axis.title=element_text(size=8,face="bold"),
            legend.position = "none")+ 
      guides(colour = guide_legend(override.aes = list(size=4),nrow=4)) + 
      xlab('Subtype') + ylab('AUROC (95% CI)')
    
    png(paste0(name,'.combined.auc.subtype.',c,'.png'),height = 600,width=1000,res = 300)
    print(plot1)
    dev.off()
    
    
    plot1 = ggplot(combined.auroc[combined.auroc$var.group =='Age at Diagnosis',],aes(x = var.plot, y = auc, fill = var,ymin=auc.lower, ymax=auc.upper)) +
      geom_bar(stat='identity')+
      scale_fill_manual(values = c(age50_colors,subtype_colors))+ #+ ggtitle(title) +
      theme_bw()+ 
      geom_errorbar(width=.2,
                    position=position_dodge(.9)) +
      
      #facet_grid2(. ~ var.group,scales='free', independent = "x")+
      #force_panelsizes(cols = c(0.5, 0.5)) +
      scale_y_continuous(limits=c(0,1))+
      theme(text = element_text(size=8),
            axis.text=element_text(size=8),
            axis.title=element_text(size=8,face="bold"),
            legend.position = "none")+ 
      guides(colour = guide_legend(override.aes = list(size=4),nrow=4)) + 
      xlab('Age at Diagnosis') + ylab('AUROC (95% CI)')
    
    png(paste0(name,'.combined.auc.dxage.',c,'.png'),height = 600,width=1000,res = 300)
    print(plot1)
    dev.off()
    
    
    plot1 = ggplot(combined.auroc,aes(x = var.plot, y = se.spec.95,ymin=se.spec.95.lower, ymax=se.spec.95.upper, fill = var)) + geom_bar(stat='identity')+
      scale_fill_manual(values = c(age50_colors,subtype_colors))+ #+ ggtitle(title) +
      theme_bw()+ 
      geom_errorbar(width= 0.2, position=position_dodge(.9)) + 
      facet_grid2(. ~ var.group,scales='free', independent = "x")+
      force_panelsizes(cols = c(1, 0.5)) +
      
      scale_y_continuous(limits=c(0,1))+
      theme(text = element_text(size=8),
            axis.text=element_text(size=8),
            axis.title=element_text(size=8,face="bold"),
            legend.position = "none")+ 
      guides(colour = guide_legend(override.aes = list(size=4),nrow=4)) + 
      xlab('Time to Diagnosis') + ylab('Sensitivity at \n95% Specificity (95% CI)')
    
    png(paste0(name,'.combined.sens.overall.',c,'.png'),height = 600,width=1000,res = 300)
    print(plot1)
    dev.off()
    
    
    
    
    plot1 = ggplot(combined.auroc,aes(x = var.plot, y = auc, fill = var,ymin=auc.lower, ymax=auc.upper)) +
      geom_bar(stat='identity')+
      scale_fill_manual(values = c(age50_colors,subtype_colors))+ #+ ggtitle(title) +
      theme_bw()+ 
      geom_errorbar(width= 0.2, position=position_dodge(.9)) + 
      facet_grid(. ~ var.group,scales='free_x',space='fixed')+
      force_panelsizes(cols = c(1, 0.5)) +
      
      scale_y_continuous(limits=c(0,1))+
      theme(text = element_text(size=8),
            axis.text=element_text(size=8),
            axis.title=element_text(size=8,face="bold"),
            legend.position = "none",
            strip.background = element_rect(fill = "white"))+
      guides(colour = guide_legend(override.aes = list(size=4),nrow=4)) + 
      xlab('Time to Diagnosis') + ylab('AUROC (95% CI)')
    
    png(paste0(name,'.combined.auroc.overall.',c,'.png'),height = 600,width=1500,res = 300)
    print(plot1)
    dev.off()
    
    
    
    
  }
  
  
  
  return(combined.auroc)
  
}
targ.samples.all =readRDS('/.mounts/labs/awadallalab/private/ncheng/cfmedip_data/cptp_samples/fragmentation/aix13.updated1/all.sample.fragmentation/regulatory.counts//aix.octane.burg.sample.info.RDS')


#octane vs non-breast controls
combined.collapse.oct.burg = combined.collapse[!grepl('AIX',combined.collapse$GRP_Id) == T ,]
combined.collapse.oct.burg$reported =ifelse(grepl('OCT_',combined.collapse.oct.burg$GRP_Id) == T, 'Cancer','Control')


a = octane.auc.plot(combined.collapse.oct.ohs,paste0(figdir,'octane.ohs'),targ.samples.all)
b = octane.auc.plot(combined.collapse.oct.burg,paste0(figdir,'octane.burg'),targ.samples.all)
saveRDS(a,paste0(figdir,'octane.ohs.auroc.RDS'))
saveRDS(b,paste0(figdir,'octane.burg.auroc.RDS'))



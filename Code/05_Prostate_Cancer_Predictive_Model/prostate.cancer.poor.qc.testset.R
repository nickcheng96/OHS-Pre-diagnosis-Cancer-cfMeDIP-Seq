#!/usr/bin/env Rscript
library(DESeq2)
library(parallel)
library(caret)
library(glmnet)
library(ROCR)
library(plyr)



combined.sample.info = readRDS('combined.burgener.chen.ohs.silencer.sample.info.RDS') #upload


dv.sample.info=readRDS('sample.info.RDS')
sample.info.filt.ecx = readRDS('qcfail.sample.info.RDS') #upload with formatted sample info
sample.info.filt.ecx = sample.info.filt.ecx[sample.info.filt.ecx$Sex == 'Male',]
sample.info.filt.ecx$data.partition = 'Test'
combined.sample.info = rbind(dv.sample.info[,colnames(dv.sample.info) %in% colnames(sample.info.filt.ecx)], sample.info.filt.ecx[,colnames(sample.info.filt.ecx) %in% colnames(dv.sample.info)])
combined.sample.info.targ = combined.sample.info

combined.sample.info.targ$group = factor(as.character(ifelse(combined.sample.info.targ$Cancer == 'Control','Control','Cancer')),levels= c('Control','Cancer'))
all.sample.info = combined.sample.info.targ

#setting up counts
raw.counts = readRDS('Silencer.raw.counts.qcfail.RDS')
combined.sample.info.targ = combined.sample.info[combined.sample.info$GRP_Id %in% colnames(ohs.raw.counts.filt.merged),]
windows = rownames(raw.counts)
windows = windows[!grepl('chrX|chrY',windows)]
ohs.raw.counts.filt.merged = raw.counts[windows,]
dds <- DESeqDataSetFromMatrix(countData = ohs.raw.counts.filt.merged[,all.sample.info$GRP_Id],
                              colData = all.sample.info,
                              design= ~ Cancer ) #can add gender here but we're only looking at female samples currently


dds = estimateSizeFactors(dds[windows,targ.samples$GRP_Id])
dds$condition = dds$group
sample.matrix = counts(dds, normalized = T)


res.df=readRDS('ohs.prostate.discovery.DMRs.RDS')
res.df$window =rownames(res.df)
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

n.features=100
score.cutoff = 0.112
results.df.all = NULL

for (fl in n.features){
  
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


#
saveRDS(results.df.all,paste0('poor.qc.combined.results.RDS'))
results.df.validation = results.df.all[!results.df.all$GRP_Id %in% sample.info.filt.ecx$GRP_Id,]
results.df.excluded = results.df.all[results.df.all$GRP_Id %in% sample.info.filt.ecx$GRP_Id,]


#correlating scores with qc metrics

results.df.all.annotated = merge(results.df.all, all.sample.info,by='GRP_Id')
results.df.all.annotated$qc.group = ifelse(results.df.all.annotated$GRP_Id %in% sample.info.filt.ecx$GRP_Id,'Fail','Pass')

#total
plot1 = ggplot(results.df.all.annotated,aes(x = total.fragment.count, y = methylation_score, col = qc.group,shape = reported)) +
  geom_point()+
  scale_color_manual(values =c('Fail'='#DD7373','Pass'='#3B3561'))+ 
  scale_shape_manual(values = c('Cancer'=19,'Control'=1))+#+ ggtitle(title) +
  theme_bw()+ 
  scale_y_continuous(limits=c(0,1))+
  expand_limits(y = 0)+
  theme(text = element_text(size=8),
        axis.text=element_text(size=8),
        axis.title=element_text(size=8,face="bold"),
        legend.position = "none")+ 
  geom_vline(xintercept = 5000000,linetype = 'dotted')+
  
  guides(colour = guide_legend(override.aes = list(size=4),nrow=4)) +
  xlab('Total Deduplicated Reads') + ylab('Methylation Score') 


png(paste0(figdir,'total.methscore.correlation.png'),height = 1000,width=1000,res = 300)
print(plot1)
dev.off()

#total
plot1 = ggplot(results.df.all.annotated,aes(x = total.fragment.count, y = methylation_score, col = qc.group,shape = reported)) +
  geom_point()+
  scale_color_manual(values =c('Fail'='#DD7373','Pass'='#3B3561'))+ 
  scale_shape_manual(values = c('Cancer'=19,'Control'=1))+#+ ggtitle(title) +
  theme_bw()+ 
  scale_y_continuous(limits=c(0,1))+
  expand_limits(y = 0)+
  theme(text = element_text(size=8),
        axis.text=element_text(size=10),
        axis.title=element_text(size=12,face="bold"),
        legend.position = "right")+ 
  guides(colour = guide_legend(ncol=1)) + 
  geom_vline(xintercept = 5000000)+
  xlab('Total Deduplicated Reads') + ylab('Methylation Score') 


png(paste0(figdir,'total.methscore.correlation.labs.png'),height = 1000,width=1000,res = 300)
print(plot1)
dev.off()



#thalianabeta
plot1 = ggplot(results.df.all.annotated,aes(x = THALIANA_BETA, y = methylation_score, col = qc.group,shape = reported)) +
  geom_point()+
  scale_color_manual(values =c('Fail'='#DD7373','Pass'='#3B3561'))+ 
  scale_shape_manual(values = c('Cancer'=19,'Control'=1))+#+ ggtitle(title) +
  theme_bw()+ 
  scale_y_continuous(limits=c(0,1))+
  expand_limits(y = 0)+
  theme(text = element_text(size=8),
        axis.text=element_text(size=8),
        axis.title=element_text(size=8,face="bold"),
        legend.position = "none")+ 
  geom_vline(xintercept = 0.95,linetype = 'dotted')+
  
  guides(colour = guide_legend(override.aes = list(size=4),nrow=4)) +
  xlab('Methylated Spike-In Proportion') + ylab('Methylation Score') 


png(paste0(figdir,'tb.methscore.correlation.png'),height = 1000,width=1000,res = 300)
print(plot1)
dev.off()

#relH

plot1 = ggplot(results.df.all.annotated,aes(x = enrichment.score.relH, y = methylation_score, col = qc.group,shape = reported)) +
  geom_point()+
  scale_color_manual(values =c('Fail'='#DD7373','Pass'='#3B3561'))+ 
  scale_shape_manual(values = c('Cancer'=19,'Control'=1))+#+ ggtitle(title) +
  theme_bw()+ 
  scale_y_continuous(limits=c(0,1))+
  expand_limits(y = 0)+
  theme(text = element_text(size=8),
        axis.text=element_text(size=8),
        axis.title=element_text(size=8,face="bold"),
        legend.position = "none")+ 
  geom_vline(xintercept = 2.5,linetype = 'dotted')+
  
  guides(colour = guide_legend(override.aes = list(size=4),nrow=4)) +
  xlab('relH Enrichment Score') + ylab('Methylation Score') 


png(paste0(figdir,'relh.methscore.correlation.png'),height = 1000,width=1000,res = 300)
print(plot1)
dev.off()

#saturation

plot1 = ggplot(results.df.all.annotated,aes(x = maxTruCor, y = methylation_score, col = qc.group,shape = reported)) +
  geom_point()+
  scale_color_manual(values =c('Fail'='#DD7373','Pass'='#3B3561'))+ 
  scale_shape_manual(values = c('Cancer'=19,'Control'=1))+#+ ggtitle(title) +
  theme_bw()+ 
  scale_y_continuous(limits=c(0,1))+
  expand_limits(y = 0)+
  theme(text = element_text(size=8),
        axis.text=element_text(size=8),
        axis.title=element_text(size=8,face="bold"),
        legend.position = "none")+ 
  geom_vline(xintercept = 0.6,linetype = 'dotted')+
  
  guides(colour = guide_legend(override.aes = list(size=4),nrow=4)) +
  xlab('Saturation Score') + ylab('Methylation Score') 


png(paste0(figdir,'saturation.methscore.correlation.png'),height = 1000,width=1000,res = 300)
print(plot1)
dev.off()

##plotting performance##
savedir='/savedir/'

figdir=paste0(savedir,'/qc.poor.figures/')
dir.create(figdir,recursive = T)
pred.df.targ = results.df.all
name = paste0(figdir,'prad.qc.poor');
sample.info = all.sample.info;
merged.df.all= all.sample.info;
score.cutoff=0.15;
cutpoint.use =T; #using Youden's cutoff
gs.plot= T #plotting gleason scores
dx.all = T

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
      
      
      
      x = merge(mean.perf.df.targ.tmp.merged,sample.info)
      x$event=ifelse(x[,'reported'] =='Cancer',1,0)
      x$reported.surv = ifelse(x[,'reported'] == 'Cancer',1,0)
      x = x[order(match(x$GRP_Id,mean.perf.df.targ.tmp.merged$GRP_Id )),]
      male.weights= weightsf.males(x)
      mean.perf.df.targ.tmp.merged = x[,c(colnames(mean.perf.df.targ.tmp.merged))]
      mean.perf.df.targ.tmp.merged$weighting = male.weights
      
      
      
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
      
      
      x = merge(mean.perf.df.targ.tmp.merged,sample.info)
      x$event=ifelse(x[,'reported'] =='Cancer',1,0)
      x$reported.surv = ifelse(x[,'reported'] == 'Cancer',1,0)
      x = x[order(match(x$GRP_Id,mean.perf.df.targ.tmp.merged$GRP_Id )),]
      male.weights= weightsf.males(x)
      mean.perf.df.targ.tmp.merged = x[,c(colnames(mean.perf.df.targ.tmp.merged))]
      mean.perf.df.targ.tmp.merged$weighting = male.weights
      
      
      
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
    merged.df.all.tmp = merge(pred.df.targ.collapse, merged.df.all,by=c('GRP_Id'))
    
    
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
      
      png(paste0(name,".scores.surv_median.all.",c,'.',s,".png"),height = 1100, width = 950,res=300)
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
      png(paste0(name,".scores.weighted.risk.all.",c,'.',s,".all.png"),height = 1100, width = 1100,res=300)
      print(plot1)
      dev.off()
      
      
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
  
  
  
  
  
  return.df.final.auroc=combined.auroc
  return(list(return.df.final.auroc,return.df.final.auct))
  
}



combined.performance = list(return.df.final.auroc,return.df.final.auct) #performance stratified by subgroups



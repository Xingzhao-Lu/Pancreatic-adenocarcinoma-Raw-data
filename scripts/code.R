
if (T) {
  dir.create("scripts")
  dir.create("results")
  dir.create("files")
  dir.create("figures")
  dir.create("origin_datas/GEO",recursive = T)
  dir.create("origin_datas/TCGA")
}
library(RColorBrewer)
library(stringr)
library(openxlsx)
library(data.table)
library(reshape2)
library(dplyr)
library(tidydr)
library(tidyr)
library(tidyverse)
library(clusterProfiler)
library(pheatmap)
library(ComplexHeatmap)
library(GSVA)
library(GSEABase)
library(fgsea)
library(corrplot)
library(colorspace)
library(survival)
library(survminer)
library(maftools)
library(vegan)
library(forcats)
library(ggpubr)
library(ggsci)
library(ggplot2)
library(rstatix)
library(ggstatsplot)
library(survcomp)
library(customLayout)
#library(ggcor)
library(ggstance)
options(stringsAsFactors = F)
source('z:/projects/codes/mg_base.R')
MG_Grobal_baseFolder='z:/projects/codes'
merge_data_by_group=function(datExpr,anno=NULL,method=c('mean','median','max','min')[1],rm_muti=T){
  library(dplyr)
  
  anno=unique(anno[,1:2])
  #gp=table(anno[,1])
  anno=anno[which(!is.na(anno[,2])&!anno[,2]==''),]
  
  an.cmp=intersect(anno[,1],row.names(datExpr))
  if(rm_muti){
    if(ncol(datExpr)==1){
      test.data=cbind(datExpr[match(an.cmp,row.names(datExpr)),])
      row.names(test.data)=an.cmp
      colnames(test.data)=colnames(datExpr)
    }else{
      test.data=datExpr[match(an.cmp,row.names(datExpr)),]
    }
    test.data=as.data.frame(test.data,full = T)
    #class(test.data[,1])
    #head(test.data)
    #sum(is.na(test.data[,1]))
    #head(datExpr)
    test.data$MG_SXR_Group=anno[match(an.cmp,anno[,1]),2]
  }else{
    anno=anno[anno[,1]%in%an.cmp,]
    if(ncol(datExpr)==1){
      test.data=crbind2DataFrame(cbind(datExpr[match(anno[,1],row.names(datExpr)),]),full = T)
      colnames(test.data)=colnames(datExpr)
      row.names(test.data)=anno[,1]
    }else{
      test.data=crbind2DataFrame(datExpr[match(anno[,1],row.names(datExpr)),],full = T)
    }
    test.data$MG_SXR_Group=anno[,2]
  }
  
  
  #which(test.data[,1]=='Inf'&test.data[,1]=='NaN'&test.data[,1]=='NA')
  
  vd.test.data=NULL
  if(method=='mean'){
    vd.test.data=test.data %>% dplyr::group_by(MG_SXR_Group) %>% dplyr::summarise_each(dplyr::funs(mean_rm_na),vals=c(colnames(test.data)[1:(ncol(test.data)-1)]))
  }else if(method=='median'){
    vd.test.data=test.data %>% dplyr::group_by(MG_SXR_Group) %>% dplyr::summarise_each(dplyr::funs(median_rm_na),vals=c(colnames(test.data)[1:(ncol(test.data)-1)]))
  }else if(method=='max'){
    vd.test.data=test.data %>% dplyr::group_by(MG_SXR_Group) %>% dplyr::summarise_each(dplyr::funs(max_rm_na),vals=c(colnames(test.data)[1:(ncol(test.data)-1)]))
  }else if(method=='min'){
    vd.test.data=test.data %>% dplyr::group_by(MG_SXR_Group) %>% dplyr::summarise_each(dplyr::funs(min_rm_na),vals=c(colnames(test.data)[1:(ncol(test.data)-1)]))
  }
  if(!is.null(vd.test.data)){
    vd.test.data=as.data.frame(vd.test.data)
    #print(head(vd.test.data))
    #print(row.names(vd.test.data))
    #print(sum(is.na(vd.test.data[,1])))
    #print(sum(vd.test.data[,1]==''))
    if(ncol(vd.test.data)==2){
      vd.test.data.1=cbind(vd.test.data[,2])
      row.names(vd.test.data.1)=as.character(vd.test.data[,1])
      vd.test.data=vd.test.data.1
      #head(vd.test.data)
    }else{
      row.names(vd.test.data)=as.character(vd.test.data[,1])
      vd.test.data=vd.test.data[,-1]
    }
    colnames(vd.test.data)=colnames(test.data)[1:(ncol(test.data)-1)]
  }
  return(vd.test.data)
  detach('package:dplyr')
}

#TCGA################
tcga.pheno.all=read.delim('origin_datas/TCGA/TCGA.PAAD.sampleMap_PAAD_clinicalMatrix')
colnames(tcga.pheno.all)
table(tcga.pheno.all$X_primary_disease)
table(tcga.pheno.all$X_primary_site)

tcga.pheno=data.frame(Samples=tcga.pheno.all$sampleID,
                      Age=tcga.pheno.all$age_at_initial_pathologic_diagnosis,
                      Gender=tcga.pheno.all$gender,
                      tobacco_smoking_history=tcga.pheno.all$tobacco_smoking_history,
                      tcga.pheno.all[,c('pathologic_T','pathologic_N','pathologic_M','pathologic_stage')])
table(tcga.pheno$pathologic_T)
tcga.pheno$pathologic_T[tcga.pheno$pathologic_T %in% c('TX','[Discrepancy]')]=NA

table(tcga.pheno$pathologic_N)
tcga.pheno$pathologic_N=gsub('[abc]','',tcga.pheno$pathologic_N)
tcga.pheno$pathologic_N[tcga.pheno$pathologic_N %in% c('NX','')]=NA

table(tcga.pheno$pathologic_M)
tcga.pheno$pathologic_M[tcga.pheno$pathologic_M %in% c('MX','')]=NA

table(tcga.pheno$pathologic_stage)
tcga.pheno$pathologic_stage=gsub('[AB]','',tcga.pheno$pathologic_stage)
tcga.pheno$pathologic_stage[tcga.pheno$pathologic_stage %in% c('[Discrepancy]','')]=NA


tcga.survival=read.delim('origin_datas/TCGA/survival_PAAD_survival.txt')
head(tcga.survival)
tcga.survival=tcga.survival[,c(1,3:10)]
colnames(tcga.survival)[1]='Samples'
tcga.survival$OS.time
tcga.survival=tcga.survival %>% drop_na(OS.time) %>%subset(OS.time>0)
dim(tcga.survival)
#195

tcga.cli=merge(tcga.pheno,tcga.survival,by='Samples')
rownames(tcga.cli)=tcga.cli$Samples
dim(tcga.cli)
#195
head(tcga.cli)


tcga.exp=read.delim('origin_datas/TCGA/TCGA.PAAD.sampleMap_HiSeqV2.gz',row.names = 1,check.names = F)
tcga.exp[1:5,1:5]
dim(tcga.exp)
range(tcga.exp)

TCGA.samples=intersect(tcga.cli$Samples,colnames(tcga.exp))
length(TCGA.samples)
#182
tcga.cli=tcga.cli[TCGA.samples,]
tcga.exp=tcga.exp[,TCGA.samples]
dim(tcga.cli);dim(tcga.exp)
# [1] 182  16
# [1] 20530   182

save(tcga.exp,tcga.cli,file = 'results/00.pre_data/TCGA_data.RData')





#GSE57495############
load('origin_datas/GEO/GSE57495.RData')

GSE57495.cli=GSE57495$Sample
GSE57495.cli=data.frame(Samples=GSE57495.cli$Acc,Status=GSE57495.cli$vital.status,
                        OS.time=GSE57495.cli$`overall survival (month)`)
GSE57495.cli$OS.time
table(GSE57495.cli$Status)
GSE57495.cli$OS=ifelse(GSE57495.cli$Status=='ALIVE',0,1)


GSE57495.exp=GSE57495$Exp$GPL15048_60607_Data_col1
GSE57495.exp=exp_probe2symbol_v2(datExpr = GSE57495.exp,GPL = 'GPL15048')
GSE57495.exp[1:5,1:5]

save(GSE57495.exp,GSE57495.cli,file = 'results/00.pre_data/GSE57495_data.RData')





dir.create('results/01.TGFB_genes')

gmt_files <- list.files(path = "origin_datas/TGFB_geneset/",
                        pattern = "\\.gmt$", full.names = TRUE)


read_gmt <- function(file) {
  lines <- readLines(file)
  genes <- unlist(lapply(lines, function(x) {
    parts <- strsplit(x, "\t")[[1]]
    if (length(parts) > 2) {
      return(parts[-c(1,2)])   (GeneSetName, Description)
    } else {
      return(NULL)
    }
  }))
  return(genes)
}


TGFB.genes <- unique(unlist(lapply(gmt_files, read_gmt)))
TGFB.genes=intersect(rownames(tcga.exp),TGFB.genes)

length(TGFB.genes)   
#149
head(TGFB.genes)     


TGFB.genes.cox=cox_batch(dat = tcga.exp[TGFB.genes,tcga.cli$Samples],
                         time = tcga.cli$OS.time,event = tcga.cli$OS)
TGFB.genes.cox=na.omit(TGFB.genes.cox)
table(TGFB.genes.cox$p.value<0.01)


pre.genes=rownames(TGFB.genes.cox[TGFB.genes.cox$p.value<0.01,])
length(pre.genes)
#24
bioForest=function(rt=null,col){
  
  gene <- rownames(rt)
  hr <- sprintf("%.3f",rt[,2])
  hrLow  <- sprintf("%.3f",rt[,3])
  hrHigh <- sprintf("%.3f",rt[,4])
  Hazard.ratio <- paste0(hr,"(",hrLow,"-",hrHigh,")")
  pVal <- ifelse(rt[,1]<0.001, "<0.001", sprintf("%.3f", rt[,1]))
  

  n <- nrow(rt)
  nRow <- n+1
  ylim <- c(1,nRow)
  layout(matrix(c(1,2),nc=2),width=c(3,2.5))
  

  xlim = c(0,3)
  par(mar=c(4,2.5,2,1))
  plot(1,xlim=xlim,ylim=ylim,type="n",axes=F,xlab="",ylab="")
  text.cex=0.8
  text(0,n:1,gene,adj=0,cex=text.cex)
  text(1.5-0.5*0.2,n:1,pVal,adj=1,cex=text.cex);text(1.5-0.5*0.2,n+1,'pvalue',cex=text.cex,font=2,adj=1)
  text(3.1,n:1,Hazard.ratio,adj=1,cex=text.cex);text(3.1,n+1,'Hazard ratio',cex=text.cex,font=2,adj=1)
  

  par(mar=c(4,1,2,1),mgp=c(2,0.5,0))
  xlim = c(0,max(as.numeric(hrLow),as.numeric(hrHigh)))
  plot(1,xlim=xlim,ylim=ylim,type="n",axes=F,ylab="",xaxs="i",xlab="Hazard ratio")
  arrows(as.numeric(hrLow),n:1,as.numeric(hrHigh),n:1,angle=90,code=3,length=0.05,col="darkblue",lwd=2.5)
  abline(v=1,col="black",lty=2,lwd=2)
  boxcolor = ifelse(as.numeric(hr) > 1, col[2], col[1])
  points(as.numeric(hr), n:1, pch = 15, col = boxcolor, cex=1.5)
  axis(1)
}

pdf('results/01.TGFB_genes/Figure1.pdf',height = 6,width = 8)
bioForest(rt = TGFB.genes.cox[pre.genes,],col=c('skyblue','orange'))
dev.off()

write.csv(TGFB.genes.cox[pre.genes,],'results/01.TGFB_genes/TGFB.genes.cox.csv')

tcga_model_data <- cbind(tcga.cli[, c("OS.time", "OS")],
                         t(tcga.exp[pre.genes, tcga.cli$Samples]))
tcga_model_data=as.data.frame(tcga_model_data)


#02.LASSO###########
dir.create('results/02.LASSO')
library(glmnet)
set.seed(2025)
fit1=glmnet(as.matrix(tcga_model_data[,-c(1,2)])
            ,cbind(time=tcga_model_data$OS.time,
                   status=tcga_model_data$OS)
            ,family="cox"
            ,nlambda=100
            , alpha=1) 

cv.fit<-cv.glmnet(as.matrix( tcga_model_data[,-c(1,2)])
                  ,cbind(time=tcga_model_data$OS.time,
                         status=tcga_model_data$OS)
                  ,family="cox"
                  ,nfolds = 10
                  ,nlambda=100
                  , alpha=1)

sig.coef <- coefficients(cv.fit,s=cv.fit$lambda.min)[which(coefficients(cv.fit,s=cv.fit$lambda.min)[,1]!=0),1]
print(cv.fit$lambda.min)#0.04020895
length(names(sig.coef))#10
names(sig.coef)

pdf('results/02.LASSO/LASSO.pdf',height = 4,width = 8,onefile = F)
par(mfrow=c(1,2))
plot(fit1)
plot(cv.fit)
dev.off()


fmla <- as.formula(paste0("Surv(OS.time,OS) ~"
                          ,paste0(names(sig.coef),collapse = '+')))
cox <- coxph(fmla, data =as.data.frame(tcga_model_data))
cox=step(cox)
lan <- coef(cox)
lan
paste0(round(lan, 3), '*', names(lan),collapse = ')+(')


gene.coef=data.frame(gene=names(lan),coef=as.numeric(lan))
gene.coef
gene.coef$Type=ifelse(gene.coef$coef>0,'Risk','Protective')
gene.coef.fig=ggplot(gene.coef, aes(x = coef, y = reorder(gene,coef), fill =Type)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = c("#EF9708","#009B9F")) +
  labs(x = 'coefficient', y = "") +
  geom_text(aes(label = round(coef,3),hjust =2), data = subset(gene.coef, coef > 0))+ 
  geom_text(aes(label = round(coef,3), hjust = -1), data = subset(gene.coef, coef < 0))+  
  theme_bw()+ theme(text = element_text(size = 14,color='black'),
                    legend.position = 'top')
gene.coef.fig
ggsave('results/02.LASSO/gene_coef.pdf',gene.coef.fig,height = 3.5,width = 4.5)

module.coxforest=ggforest(cox, data = tcga_model_data,
                          main = "Hazardratio", fontsize =1,
                          noDigits = 2)
module.coxforest
ggsave('results/02.LASSO/module_coxforest.pdf',module.coxforest,height =5,width = 10)



dir.create('results/03.prognostic_model')

risktype.col=c( "#D691C1","#43B7B9")
tcga_model_data <- cbind(tcga.cli[, c("OS.time", "OS")],
                         t(tcga.exp[names(lan), tcga.cli$Samples]))
tcga_model_data=as.data.frame(tcga_model_data)

risk.tcga=as.numeric(lan%*%as.matrix(t(tcga_model_data[tcga.cli$Samples,names(lan)])))
tcga.risktype.cli=data.frame(tcga.cli,Riskscore=risk.tcga)
tcga.cutoff<-surv_cutpoint(tcga.risktype.cli,time="OS.time",event="OS",
                      variables='Riskscore')
tcga.risktype.cli$Risktype=ifelse(tcga.risktype.cli$Riskscore>tcga.cutoff$cutpoint$cutpoint,'High','Low')
table(tcga.risktype.cli$Risktype)

tcga.roc.OS=ggplotTimeROC(tcga.risktype.cli$OS.time,
                           tcga.risktype.cli$OS,
                           tcga.risktype.cli$Riskscore,mks = c(1:5))
tcga.roc.OS

tcga.km.OS=ggsurvplot(fit=survfit(Surv(OS.time/365, OS) ~ Risktype,
                                   data = tcga.risktype.cli),
                       data=tcga.risktype.cli,
                       conf.int = T,pval = T,risk.table = T, 
                       fun = "pct",size = 1,surv.median.line = 'hv',
                       title='TCGA',legend.title='Risktype',
                       legend.labs = c('High','Low'),
                       linetype = c("solid", "dashed","strata")[1],
                       palette = risktype.col,
                       legend.position='top')
tcga.km.OS=mg_merge_plot(tcga.km.OS$plot,tcga.km.OS$table,nrow=2,
                          heights = c(2.5,1),align = 'v')
tcga.km.OS


GSE57495_model_data <- cbind(GSE57495.cli[, c("OS.time", "OS")],
                             t(GSE57495.exp[names(lan), GSE57495.cli$Samples]))
GSE57495_model_data

risk.GSE57495=as.numeric(lan%*%as.matrix(t(GSE57495_model_data[GSE57495.cli$Samples,names(lan)])))
GSE57495.risktype.cli=data.frame(GSE57495.cli,Riskscore=risk.GSE57495)
GSE57495.cutoff<-surv_cutpoint(GSE57495.risktype.cli,time="OS.time",event="OS",
                      variables='Riskscore')
GSE57495.risktype.cli$Risktype=ifelse(GSE57495.risktype.cli$Riskscore>GSE57495.cutoff$cutpoint$cutpoint,'High','Low')

table(GSE57495.risktype.cli$Risktype)

GSE57495.roc.OS=ggplotTimeROC(time = GSE57495.risktype.cli$OS.time,
                              status = GSE57495.risktype.cli$OS,
                              score = GSE57495.risktype.cli$Riskscore,mks = c(1:5))
GSE57495.roc.OS

GSE57495.km.OS=ggsurvplot(fit=survfit( Surv(OS.time/12, OS) ~ Risktype,
                                       data = GSE57495.risktype.cli),
                          data=GSE57495.risktype.cli,
                          conf.int = T,pval = T,risk.table = T, 
                          fun = "pct",size = 1,surv.median.line = 'hv',
                          title='GSE57495',legend.title='Risktype',
                          legend.labs = c('High','Low'),
                          linetype = c("solid", "dashed","strata")[1],
                          palette = risktype.col,
                          legend.position='top')
GSE57495.km.OS=mg_merge_plot(GSE57495.km.OS$plot,GSE57495.km.OS$table,nrow=2,heights = c(2.5,1),align = 'v')
GSE57495.km.OS


df=cbind.data.frame(t(tcga.exp[names(lan),tcga.risktype.cli$Samples]),Risktype=tcga.risktype.cli$Risktype)
head(df)
df=melt(df)
head(df)

tcga.model.expr=ggviolin(df, x = "variable", y = "value", fill = "Risktype",add = "boxplot")+
  ggpubr::stat_compare_means(aes(group=Risktype), label = "p.signif", method = 't.test')+
  scale_fill_manual(values = risktype.col)+ylab('Expression')+xlab('')+
  theme(text = element_text(size=14))+ggtitle('TCGA')
tcga.model.expr




df=cbind.data.frame(t(GSE57495.exp[names(lan),GSE57495.risktype.cli$Samples]),Risktype=GSE57495.risktype.cli$Risktype)
head(df)
df=melt(df)
head(df)

GSE57495.model.expr=ggviolin(df, x = "variable", y = "value", fill = "Risktype",add = "boxplot")+
  ggpubr::stat_compare_means(aes(group=Risktype), label = "p.signif", method = 't.test')+
  scale_fill_manual(values = risktype.col)+ylab('Expression')+xlab('')+
  theme(text = element_text(size=14))+ggtitle('GSE57495')
GSE57495.model.expr

pdf('results/03.prognostic_model/Figure3.pdf',height = 10,width = 15)
mg_merge_plot(tcga.roc.OS,tcga.km.OS,tcga.model.expr,
              GSE57495.roc.OS,GSE57495.km.OS,GSE57495.model.expr,
              nrow=2,ncol=3,labels = LETTERS[1:6])
dev.off()
# 


# GSE71729_model_data <- cbind(GSE71729.cli[, c("OS.time", "OS")],
#                              t(GSE71729.exp[names(lan), GSE71729.cli$Samples]))
# # colnames(GSE71729_model_data) <- gsub('-', '_', colnames(GSE71729_model_data))
# GSE71729_model_data
# 
# risk.GSE71729=as.numeric(lan%*%as.matrix(t(GSE71729_model_data[GSE71729.cli$Samples,names(lan)])))
# GSE71729.risktype.cli=data.frame(GSE71729.cli,Riskscore=risk.GSE71729)
# GSE71729.risktype.cli=crbind2DataFrame(GSE71729.risktype.cli)
# GSE71729.cutoff<-surv_cutpoint(GSE71729.risktype.cli,time="OS.time",event="OS",
#                                variables='Riskscore')
# GSE71729.risktype.cli$Risktype=ifelse(GSE71729.risktype.cli$Riskscore>GSE71729.cutoff$cutpoint$cutpoint,'High','Low')
# 
# # GSE71729.risktype.cli$Risktype=ifelse(GSE71729.risktype.cli$Riskscore>median(risk.GSE71729),'High','Low')
# table(GSE71729.risktype.cli$Risktype)
# 
# GSE71729.roc.OS=ggplotTimeROC(time = GSE71729.risktype.cli$OS.time,
#                               status = GSE71729.risktype.cli$OS,
#                               score = GSE71729.risktype.cli$Riskscore,mks = c(1:5))
# GSE71729.roc.OS
# 
# GSE71729.km.OS=ggsurvplot(fit=survfit( Surv(OS.time/12, OS) ~ Risktype,
#                                        data = GSE71729.risktype.cli),
#                           data=GSE71729.risktype.cli,
#                           conf.int = T,pval = T,risk.table = T, 
#                           fun = "pct",size = 1,surv.median.line = 'hv',
#                           title='GSE71729',legend.title='Risktype',
#                           legend.labs = c('High','Low'),
#                           linetype = c("solid", "dashed","strata")[1],
#                           palette = risktype.col,
#                           legend.position='top')
# GSE71729.km.OS=mg_merge_plot(GSE71729.km.OS$plot,GSE71729.km.OS$table,nrow=2,heights = c(2.5,1),align = 'v')
# GSE71729.km.OS
# 
# 
# cptac_model_data <- cbind(cptac.cli[, c("OS.time", "OS")],
#                           t(cptac.exp[names(lan), cptac.cli$Samples]))
# # colnames(cptac_model_data) <- gsub('-', '_', colnames(cptac_model_data))
# cptac_model_data
# 
# risk.cptac=as.numeric(lan%*%as.matrix(t(cptac_model_data[cptac.cli$Samples,names(lan)])))
# cptac.risktype.cli=data.frame(cptac.cli,Riskscore=risk.cptac)
# 
# cptac.risktype.cli=crbind2DataFrame(cptac.risktype.cli)
# cptac.cutoff<-surv_cutpoint(cptac.risktype.cli,time="OS.time",event="OS",
#                             variables='Riskscore')
# cptac.risktype.cli$Risktype=ifelse(cptac.risktype.cli$Riskscore>cptac.cutoff$cutpoint$cutpoint,'High','Low')
# table(cptac.risktype.cli$Risktype)
# 
# cptac.roc.OS=ggplotTimeROC(time = cptac.risktype.cli$OS.time,
#                            status = cptac.risktype.cli$OS,
#                            score = cptac.risktype.cli$Riskscore,mks = c(1:5))
# cptac.roc.OS
# 
# cptac.km.OS=ggsurvplot(fit=survfit( Surv(OS.time/365, OS) ~ Risktype,
#                                     data = cptac.risktype.cli),
#                        data=cptac.risktype.cli,
#                        conf.int = T,pval = T,risk.table = T, 
#                        fun = "pct",size = 1,surv.median.line = 'hv',
#                        title='cptac',legend.title='Risktype',
#                        legend.labs = c('High','Low'),
#                        linetype = c("solid", "dashed","strata")[1],
#                        palette = risktype.col,
#                        legend.position='top')
# cptac.km.OS=mg_merge_plot(cptac.km.OS$plot,cptac.km.OS$table,nrow=2,heights = c(2.5,1),align = 'v')
# cptac.km.OS
# 
# 
# qcmg_model_data <- cbind(qcmg.cli[, c("OS.time", "OS")],
#                          t(qcmg.exp[names(lan), qcmg.cli$Samples]))
# # colnames(qcmg_model_data) <- gsub('-', '_', colnames(qcmg_model_data))
# qcmg_model_data
# 
# risk.qcmg=as.numeric(lan%*%as.matrix(t(qcmg_model_data[qcmg.cli$Samples,names(lan)])))
# qcmg.risktype.cli=data.frame(qcmg.cli,Riskscore=risk.qcmg)
# 
# qcmg.risktype.cli=crbind2DataFrame(qcmg.risktype.cli)
# qcmg.cutoff<-surv_cutpoint(qcmg.risktype.cli,time="OS.time",event="OS",
#                            variables='Riskscore')
# qcmg.risktype.cli$Risktype=ifelse(qcmg.risktype.cli$Riskscore>qcmg.cutoff$cutpoint$cutpoint,'High','Low')
# table(qcmg.risktype.cli$Risktype)
# 
# qcmg.roc.OS=ggplotTimeROC(time = qcmg.risktype.cli$OS.time,
#                           status = qcmg.risktype.cli$OS,
#                           score = qcmg.risktype.cli$Riskscore,mks = c(1:5))
# qcmg.roc.OS
# 
# qcmg.km.OS=ggsurvplot(fit=survfit( Surv(OS.time/365, OS) ~ Risktype,
#                                    data = qcmg.risktype.cli),
#                       data=qcmg.risktype.cli,
#                       conf.int = T,pval = T,risk.table = T, 
#                       fun = "pct",size = 1,surv.median.line = 'hv',
#                       title='qcmg',legend.title='Risktype',
#                       legend.labs = c('High','Low'),
#                       linetype = c("solid", "dashed","strata")[1],
#                       palette = risktype.col,
#                       legend.position='top')
# qcmg.km.OS=mg_merge_plot(qcmg.km.OS$plot,qcmg.km.OS$table,nrow=2,heights = c(2.5,1),align = 'v')
# qcmg.km.OS
# 
dir.create('results/04.TME')
tcga.estimate=read.delim('results/04.TME/TCGA_ESTIMATE_score.txt',row.names = 1,check.names = F)
head(tcga.estimate)
tcga.estimate.df=tcga.estimate[tcga.risktype.cli$Samples,-4]
tcga.estimate.df$group=tcga.risktype.cli$Risktype

tcga.estimate.df=melt(tcga.estimate.df,id.vars = 'group')
head(tcga.estimate.df)
fig4a=ggplot(data = tcga.estimate.df,aes(x = variable,y = value,fill = group))+
  geom_boxplot(position = position_dodge(width = 0.5),width = 0.4)+
  scale_fill_manual(values = risktype.col)+xlab('')+ylab('Score')+
  stat_compare_means(aes(group=group), label = "p.signif", method = 't.test')+
  theme_classic()+
  theme(text = element_text(size=14),legend.position = 'top',
        axis.text.x = element_text(angle = 45,hjust = 1,color='black',size=12),
        axis.text.y = element_text(color='black',size=12))
fig4a


tcga.cibersort=read.delim('results/04.TME/TCGA_CIBERSORT_Results.txt',row.names = 1,check.names = F)
head(tcga.cibersort)
tcga.cibersort=tcga.cibersort[tcga.risktype.cli$Samples,] %>% 
  mutate(group=tcga.risktype.cli$Risktype)  %>% 
  subset(`P-value`<0.05) %>% 
  select(c(26,1:22))
head(tcga.cibersort)
tcga.cibersort.df=melt(tcga.cibersort,id.vars = 'group')
head(tcga.cibersort.df)
fig4b=ggplot(data = tcga.cibersort.df,aes(x = variable,y = value,fill = group))+
  geom_boxplot(position = position_dodge(width = 0.5),width = 0.4)+
  scale_fill_manual(values = risktype.col)+xlab('')+ylab('Fraction')+
  stat_compare_means(aes(group=group), label = "p.signif", method = 't.test')+
  theme_classic()+
  theme(text = element_text(size=14),legend.position = 'top',
        axis.text.x = element_text(angle = 45,hjust = 1,color='black',size=12),
        axis.text.y = element_text(color='black',size=12))
fig4b

pdf('results/04.TME/fig4ab.pdf',height = 4.5,width = 12)
mg_merge_plot(fig4a,fig4b,widths = c(1,2.5),labels = c('A','B'),
              align = 'h',common.legend = T)
dev.off()

tcga.estimate.cor=data.frame(t(tcga.exp[names(lan),tcga.risktype.cli$Samples]),
                             tcga.estimate[tcga.risktype.cli$Samples,1:3])
cor_res <- Hmisc::rcorr(as.matrix(tcga.estimate.cor),type = 'pearson')
cor_res$P[is.na(cor_res$P)] <- 0
library(corrplot)
pdf('results/04.TME/fig4c.pdf',height = 5,width = 5)
corrplot(as.matrix(cor_res$r[names(lan),colnames(tcga.estimate)[1:3]]),
         p.mat = as.matrix(cor_res$P[names(lan),colnames(tcga.estimate)[1:3]]),
         mar = c(0,0,1,1),
         col=colorRampPalette(c('blue', 'white','red'))(100),
         tl.srt = 90,tl.cex = 1,tl.col = 'black',tl.offset = 0.5,
         cl.pos = c("b","r","n")[1],cl.align.text = 'l',cl.length = 5,
         cl.ratio = 0.1,cl.cex = 0.8,
         addgrid.col = 'white',
         method = c("circle", "square", "ellipse", "number", "shade", "color", "pie")[6],
         insig = 'label_sig',
         sig.level=c(0.001,0.01,0.05),
         pch.cex=1,is.corr=T,xpd=T)
dev.off()


tcga.cibersort.cor=data.frame(t(tcga.exp[names(lan),rownames(tcga.cibersort)]),
                              tcga.cibersort[,-1],check.names = F)
cor_res <- Hmisc::rcorr(as.matrix(tcga.cibersort.cor),type = 'pearson')
cor_res$P[is.na(cor_res$P)] <- 0
pdf('results/04.TME/fig4d.pdf',height = 5,width = 10)
corrplot(as.matrix(cor_res$r[names(lan),colnames(tcga.cibersort)[-1]]),
         p.mat = as.matrix(cor_res$P[names(lan),colnames(tcga.cibersort)[-1]]),
         mar = c(0,0,1,1),
         col=colorRampPalette(c('blue', 'white','red'))(100),
         tl.srt = 90,tl.cex = 1,tl.col = 'black',tl.offset = 0.5,
         cl.pos = c("b","r","n")[2],cl.align.text = 'l',cl.length = 5,
         cl.ratio = 0.1,cl.cex = 0.8,
         addgrid.col = 'white',
         method = c("circle", "square", "ellipse", "number", "shade", "color", "pie")[6],
         insig = 'label_sig',
         sig.level=c(0.001,0.01,0.05),
         pch.cex=1,is.corr=T,xpd=T)
dev.off()



dir.create('results/05.Treatment')
# tcga_tide_dat <- t(scale(t(tcga.exp),scale = F))
# dim(tcga_tide_dat)
# write.table(tcga_tide_dat,file = 'results/05.Treatment/tcga_tide_dat.txt',quote = F, sep = '\t')


tcga_tide_res<-read.csv('results/05.Treatment/tcga.tide.res.csv',row.names = 1,stringsAsFactors = F)
head(tcga_tide_res)
dim(tcga_tide_res)
TIDE.df=data.frame(tcga_tide_res[tcga.risktype.cli$Samples,],
                   tcga.risktype.cli)
head(TIDE.df)
table(TIDE.df$Responder,TIDE.df$Risktype)
chisq.test(table(TIDE.df$Responder,TIDE.df$Risktype))
# data:  table(TIDE.df$Responder, TIDE.df$Risktype)
# X-squared = 6.6611, df = 1, p-value = 0.009854
tide.responder=prop.table(table(TIDE.df$Responder,TIDE.df $Risktype),margin=2)
tide.responder=melt(tide.responder)
tide.responder
fig5a=ggplot(tide.responder, aes(x= Var2, y=value, fill=Var1))+
  geom_bar(stat = "identity")+xlab('Riskscore')+ylab('Percentage')+
  scale_fill_manual(values = pal_simpsons()(9)[5:6],name='Responder')+
  theme_bw()+geom_text(data=tide.responder,aes(label=paste(round(100*value,2),'%',sep='')))+
  theme(text = element_text(size = 12),legend.position = 'top',
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.text.x = element_text(size=14,color='black'),
        axis.text.y = element_text(size=14,color='black'))
fig5a
ggsave('results/05.Treatment/Fig5a.pdf',fig5a,height = 6,width = 4)

fig5b=TIDE.df %>%
  ggplot(aes(x=Risktype, y=TIDE,fill=Risktype)) +
  geom_boxplot()+
  scale_fill_manual(values =risktype.col)+   
  ggpubr::stat_compare_means(aes(group=Risktype), label = "p.format", method = 't.test')+
  labs(x="", y = "TIDE", fill = "Risktype") +theme_bw()+
  theme(legend.position = "none",text = element_text(size = 14),
        axis.text.x = element_text(size=14,color='black'),
        axis.text.y = element_text(size=14,color='black')) 
fig5b
ggsave('results/05.Treatment/Fig5b.pdf',fig5b,height = 6,width = 4)

##GSE91061【submap########
# path1='public/imv210_geo/'
# GSE91061_samples <- read.delim(paste0(path1,'GSE91061/GSE91061.cli2.txt'),
#                                row.names = 1, check.names = F, stringsAsFactors = F)
# GSE91061_samples=data.frame(Samples=rownames(GSE91061_samples),
#                             Visit=GSE91061_samples$`visit (pre or on treatment)`,
#                             Title=stringr::str_split_fixed(GSE91061_samples$Title,'_',2)[,1],
#                             Res=GSE91061_samples$response,
#                             Treatment_res=GSE91061_samples$response1)
# head(GSE91061_samples)
# table(GSE91061_samples$Treatment_res)
# table(GSE91061_samples$Res)
# GSE91061_samples$Res[which(GSE91061_samples$Res=='PD'|GSE91061_samples$Res=='SD')]='PD/SD'
# 
# GSE91061_exp <- read.delim(paste0(path1,'GSE91061/GSE91061.exp.txt'),
#                            row.names = 1, check.names = F, stringsAsFactors = F)
# dim(GSE91061_exp)
# GSE91061_exp=GSE91061_exp[,GSE91061_samples$Samples]
# range(GSE91061_exp)
# GSE91061_exp=log2(GSE91061_exp+1)
# 
# submap_test=mg_submap(original.data.A = tcga.exp[,tcga.risktype.cli$Samples],
#                       original.data.B = GSE91061_exp[,GSE91061_samples$Samples],
#                       original.data.A.group = as.character(tcga.risktype.cli$Risktype),
#                       original.data.B.group = as.character(GSE91061_samples$Treatment_res))
# 
# save(submap_test,file =  '~/project/2025/20251028_PAAD_TGFB/results/05.Treatment/submap_res.RData')

load('results/05.Treatment/submap_res.RData')
submap_test$Bonferroni.SA.matrix
submap_test$nominal.p.matrix.Fisher

treat_heatmap=rbind(submap_test$Bonferroni.SA.matrix,
                    submap_test$nominal.p.matrix.Fisher)
treat_heatmap_pval=data.frame(pval=rep(c('nominal p value','Bonferroni corrected'),
                                       c(2,2)))
rownames(treat_heatmap_pval)=rownames(treat_heatmap)

pdf('results/05.Treatment/Fig5c.pdf',height = 6,width = 8)
pheatmap(rbind(submap_test$Bonferroni.SA.matrix,
               submap_test$nominal.p.matrix.Fisher),
         scale = 'none',name = 'p.value',
         annotation_row = treat_heatmap_pval,
         cluster_cols = F,cluster_rows = F,
         gaps_row = 2,
         color = sequential_hcl(4, palette = "Viridis",rev = T),
         legend_breaks=c(0,0.3,0.6,0.9),
         display_numbers=T,number_color = 'black')
dev.off()




load('results/05.Treatment/tcga_durg_ic50_res.RData')

tcga_durg_ic50_res[1:5,1:5]
colnames(tcga_durg_ic50_res)[1]='Cisplatin'

library(ggcorrplot)
library(psych)
IC50_RS_cor <- corr.test(x =tcga.risktype.cli$Riskscore,
                         y = tcga_durg_ic50_res[tcga.risktype.cli$Samples,],
                         method = "pearson",adjust = "BH",ci = F)

IC50_RS_cor_res=data.frame(drugs=colnames(tcga_durg_ic50_res))
IC50_RS_cor_res$cor<-as.numeric(IC50_RS_cor$r)
IC50_RS_cor_res$p.adj<-as.numeric(IC50_RS_cor$p.adj)
head(IC50_RS_cor_res)
table(IC50_RS_cor_res$p.adj<0.05,abs(IC50_RS_cor_res$cor)>0.6)
IC50_RS_cor_res=IC50_RS_cor_res%>%subset(p.adj<0.05 & abs(cor)>0.6)
IC50_RS_cor_res=IC50_RS_cor_res[order(IC50_RS_cor_res$cor),]
head(IC50_RS_cor_res)


fig5d=ggplot(data=IC50_RS_cor_res,aes(x=cor,y=reorder(drugs,cor))) +
  geom_col(aes(fill=-log10(p.adj))) +
  scale_fill_gradient(low = "skyblue", high = "orange")+
  labs(x='Correlation',y='Drugs')+theme_bw()+
  theme(text = element_text(size=14,color='black'))
fig5d
ggsave('results/05.Treatment/Fig5d.pdf',fig5d,height = 4.5,width = 5.5)

IC50.df=data.frame(tcga_durg_ic50_res[tcga.risktype.cli$Samples,IC50_RS_cor_res$drugs],
                   Risktype=tcga.risktype.cli$Risktype,check.names = F)
IC50.df=melt(IC50.df)
head(IC50.df)
fig5e=ggplot(data = IC50.df,aes(x = variable,y = value,fill = Risktype ))+
  geom_boxplot(position = position_dodge(width = 0.5),width = 0.4)+
  scale_fill_manual(values = risktype.col)+xlab('')+ylab('IC50')+
  stat_compare_means(aes(group=Risktype ), label = "p.signif", method = 'wilcox.test')+
  theme_light()+
  theme(text = element_text(size=14),legend.position = 'top',
        axis.text.x = element_text(angle = 45,hjust = 1,color='black'))
fig5e
ggsave('results/05.Treatment/Fig5e.pdf',fig5e,height =4.5,width = 10)



#06.GSEA#######################
dir.create('results/06.GSEA')

getGeneFC=function(gene.exp,group,ulab=ulab,dlab=dlab){
  degs_C1_C3=mg_limma_DEG(gene.exp, 
                          group,
                          ulab=ulab,
                          dlab=dlab)
  
  degs_C1_C3_sig<-degs_C1_C3$DEG[which(degs_C1_C3$DEG$adj.P.Val <= 1),]
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(stringr)
  degs_C1_C3_sig_gene<-rownames(degs_C1_C3_sig)
  
  degs_gene_entz=bitr(degs_C1_C3_sig_gene,fromType="SYMBOL",toType="ENTREZID",OrgDb="org.Hs.eg.db")
  degs_gene_entz <- dplyr::distinct(degs_gene_entz,SYMBOL,.keep_all=TRUE)
  
  gene_df <- data.frame(logFC=degs_C1_C3_sig$logFC,
                        SYMBOL = rownames(degs_C1_C3_sig))
  gene_df <- merge(gene_df,degs_gene_entz,by="SYMBOL")
  head(gene_df)
  
  geneList<-gene_df$logFC
  names(geneList)=gene_df$ENTREZID 
  head(geneList)
  
  geneList=sort(geneList,decreasing = T)
  head(geneList) 
  return(geneList)
}
tcga.geneList=getGeneFC(gene.exp=tcga.exp[,tcga.risktype.cli$Samples],
                        group=tcga.risktype.cli$Risktype
                        ,ulab='High',dlab = 'Low')
set.seed(123)
pathway.gmt<-read.gmt("gmt/c2.cp.kegg.v7.5.1.entrez.gmt")
tcga.hallmark.gsea<-GSEA(tcga.geneList,TERM2GENE = pathway.gmt,seed=T)

write.csv(tcga.hallmark.gsea@result,'results/06.GSEA/tcga.hallmark.gsea.res.csv')

library(enrichplot)
library(ggplot2)
library(GseaVis)
gsea.dotplot=dotplotGsea(data = tcga.hallmark.gsea,topn = 10,
                         order.by = 'NES')

pdf('results/06.GSEA/Figure6.pdf',height = 7,width = 12)
gsea.dotplot$plot+#theme(text = element_text(family = 'Times'))+
  scale_y_discrete(labels=function(x) str_remove(x,"HALLMARK_"))
dev.off()


dir.create('results/07.Mutation')

load('results/07.Mutation/tcga.maf.RData')


tcga.TMB <- tmb(maf = tcga.maf,captureSize = 50,logScale = T)
head(tcga.TMB)
dim(tcga.TMB)
tcga.TMB$Tumor_Sample_Barcode=paste0(tcga.TMB$Tumor_Sample_Barcode,'-01')

tcga.TMB.merge=merge(tcga.TMB,tcga.risktype.cli,by.x='Tumor_Sample_Barcode',by.y='Samples')
head(tcga.TMB.merge)
dim(tcga.TMB.merge)
fig7a=tcga.TMB.merge %>%drop_na(total_perMB)%>% subset(total_perMB<100) %>%
  ggplot(aes(x=Risktype,y=total_perMB  ,fill=Risktype))+
  geom_boxplot()+ stat_compare_means(aes(group=Risktype), label = "p.format", method = 't.test')+
  scale_fill_manual(values =risktype.col)+theme_bw()+
  theme(text = element_text(size = 12),legend.position = 'none',
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.text.x = element_text(size=14,color='black'),
        axis.text.y = element_text(size=14,color='black'))
fig7a
ggsave('results/07.Mutation/Fig7a.pdf',fig7a,height = 4.5,width = 3.5)



tcga.risktype.use=tcga.risktype.cli[,c('Samples','Risktype')]
table(tcga.risktype.use$Risktype)
head(tcga.risktype.use)
colnames(tcga.risktype.use)[1]='Tumor_Sample_Barcode'
tcga.risktype.use$Tumor_Sample_Barcode=substr(tcga.risktype.use$Tumor_Sample_Barcode,1,12)
tcga.risktype.high=tcga.risktype.use[which(tcga.risktype.use$Risktype=='High'),]
tcga.risktype.low=tcga.risktype.use[which(tcga.risktype.use$Risktype=='Low'),]

write.table(tcga.risktype.high,file='results/07.Mutation/tcga.risktype.high.txt')
write.table(tcga.risktype.low ,file='results/07.Mutation/tcga.risktype.low.txt')


tcga.maf1=subsetMaf(tcga.maf,tsb=tcga.maf@data$Tumor_Sample_Barcode[substr(tcga.maf@data$Tumor_Sample_Barcode,1,12) %in% tcga.risktype.high$Tumor_Sample_Barcode])
tcga.maf1<-read.maf(tcga.maf1@data,isTCGA=T,
                    clinicalData = 'results/07.Mutation/tcga.risktype.high.txt')
tcga.maf1@clinical.data

tcga.maf2=subsetMaf(tcga.maf,tsb=tcga.maf@data$Tumor_Sample_Barcode[substr(tcga.maf@data$Tumor_Sample_Barcode,1,12) %in% tcga.risktype.low$Tumor_Sample_Barcode])
tcga.maf2<-read.maf(tcga.maf2@data,isTCGA=T,
                    clinicalData = 'results/07.Mutation/tcga.risktype.low.txt')
tcga.maf2@clinical.data



mdf_cp=mafCompare(tcga.maf1,tcga.maf2,m1Name = 'High',m2Name = 'Low')
maf.compare.res=mdf_cp$results

head(maf.compare.res)
table(maf.compare.res$pval<0.05)
maf.compare.res.filter=maf.compare.res[maf.compare.res$pval<0.05,]
head(maf.compare.res.filter)
write.csv(maf.compare.res.filter,'results/07.Mutation/maf_compare.csv')

pdf('results/07.Mutation/diff_mutation.pdf',height = 2.5,width = 5,onefile = F)
forestPlot(mafCompareRes = mdf_cp,color = c('royalblue', 'maroon'),pVal = 0.05)
dev.off()

pdf('results/07.Mutation/diff_mutation_oncoplot.pdf',height = 3.5,width = 7,onefile = F)
coOncoplot(m1=tcga.maf1, m2=tcga.maf2, m1Name="High", m2Name="Low",
           genes =  maf.compare.res.filter$Hugo_Symbol)
dev.off()


tcga.character=read.table('PMC5982584_supplement.txt')
table(tcga.character$`TCGA Study`)
tcga.character$Samples=paste0(rownames(tcga.character),'-01')
rownames(tcga.character)=tcga.character$Samples
tcga.character=merge(tcga.character,tcga.risktype.cli,by='Samples')
colnames(tcga.character)

tcga.character %>% drop_na(`Indel Neoantigens`)%>%
  ggplot(aes(x=Risktype, y=`Indel Neoantigens`,fill=Risktype)) +
  geom_boxplot()+
  # geom_jitter(aes(color = Risktype), position = position_jitter(width = 0.2, height = 0))+
  scale_fill_manual(values =risktype.col)+   
  ggpubr::stat_compare_means(aes(group=Risktype), label = "p.format", method = 'wilcox.test')+
  labs(x="", y = "Indel Neoantigens", fill = "Risktype") +theme_bw()+
  theme(legend.position = "none",text = element_text(family = 'Times',size = 14)) #  


tcga.character %>% drop_na(`Aneuploidy Score`)%>%
  ggplot(aes(x=Risktype, y=`Aneuploidy Score`,fill=Risktype)) +
  geom_boxplot()+
  # geom_jitter(aes(color = Risktype), position = position_jitter(width = 0.2, height = 0))+
  scale_fill_manual(values =risktype.col)+   
  ggpubr::stat_compare_means(aes(group=Risktype), label = "p.format", method = 'wilcox.test')+
  labs(x="", y = "Aneuploidy Score", fill = "Risktype") +theme_bw()+
  theme(legend.position = "none",text = element_text(family = 'Times',size = 14)) #     
ggsave('results/07.Mutation/Fig7b.pdf',height = 4.5,width = 3.5)

tcga.character %>% drop_na(`Number of Segments`)%>%
  ggplot(aes(x=Risktype, y=`Number of Segments`,fill=Risktype)) +
  geom_boxplot()+
  # geom_jitter(aes(color = Risktype), position = position_jitter(width = 0.2, height = 0))+
  scale_fill_manual(values =risktype.col)+   
  ggpubr::stat_compare_means(aes(group=Risktype), label = "p.format", method = 'wilcox.test')+
  labs(x="", y = "Number of Segments", fill = "Risktype") +theme_bw()+
  theme(legend.position = "none",text = element_text(family = 'Times',size = 14)) #  
ggsave('results/07.Mutation/Fig7c.pdf',height = 4.5,width = 3.5)

tcga.character %>% drop_na(`Homologous Recombination Defects`)%>%
  ggplot(aes(x=Risktype, y=`Homologous Recombination Defects`,fill=Risktype)) +
  geom_boxplot()+
  # geom_jitter(aes(color = Risktype), position = position_jitter(width = 0.2, height = 0))+
  scale_fill_manual(values =risktype.col)+   
  ggpubr::stat_compare_means(aes(group=Risktype), label = "p.format", method = 'wilcox.test')+
  labs(x="", y = "Homologous Recombination Defects", fill = "Risktype") +theme_bw()+
  theme(legend.position = "none",text = element_text(family = 'Times',size = 14)) #
ggsave('results/07.Mutation/Fig7d.pdf',height = 4.5,width = 3.5)




dir.create('results/08.Nomogram')
head(tcga.risktype.cli)
tcga_cox_datas=tcga.risktype.cli
tcga_cox_datas=as.data.frame(tcga_cox_datas)
colnames(tcga_cox_datas)
fivenum(tcga_cox_datas$Age)
tcga_cox_datas$Age1=ifelse(tcga_cox_datas$Age>65,'>65','<=65')
table(tcga_cox_datas$Age1)

table(tcga_cox_datas$pathologic_T)
tcga_cox_datas$pathologic_T[tcga_cox_datas$pathologic_T %in% c('T1','T2')]<-'T1+T2'
tcga_cox_datas$pathologic_T[tcga_cox_datas$pathologic_T %in% c('T3','T4')]<-'T3+T4'

table(tcga_cox_datas$pathologic_N)
table(tcga_cox_datas$pathologic_M)

table(tcga_cox_datas$pathologic_stage)
tcga_cox_datas$pathologic_stage[tcga_cox_datas$pathologic_stage %in% c('Stage I','Stage II')]<-'Stage I+II'
tcga_cox_datas$pathologic_stage[tcga_cox_datas$pathologic_stage %in% c('Stage III','Stage IV')]<-'Stage III+IV'




#Age
tcga_cox_datas=as.data.frame(tcga_cox_datas)
Age_sig_cox <- summary(coxph(formula=Surv(OS.time, OS)~Age,
                             data=tcga_cox_datas))
Age_sig_cox_dat <- data.frame(Names=rownames(Age_sig_cox[[8]]),
                              HR = round(Age_sig_cox[[7]][,2],3),
                              lower.95 = round(Age_sig_cox[[8]][,3],3),
                              upper.95 = round(Age_sig_cox[[8]][,4],3),
                              p.value=round(Age_sig_cox[[7]][,5],3))
Age_sig_cox_dat

##Sex
Gender_sig_cox <- summary(coxph(formula=Surv(OS.time, OS)~Gender,
                                data=tcga_cox_datas))
Gender_sig_cox_dat <- data.frame(Names=rownames(Gender_sig_cox[[8]]),
                                 HR = round(Gender_sig_cox[[7]][,2],3),
                                 lower.95 = round(Gender_sig_cox[[8]][,3],3),
                                 upper.95 = round(Gender_sig_cox[[8]][,4],3),
                                 p.value=round(Gender_sig_cox[[7]][,5],3))
Gender_sig_cox_dat

##smoke
smoke_sig_cox <- summary(coxph(formula=Surv(OS.time, OS)~tobacco_smoking_history,
                               data=tcga_cox_datas))
smoke_sig_cox_dat <- data.frame(Names=rownames(smoke_sig_cox[[8]]),
                                HR = round(smoke_sig_cox[[7]][,2],3),
                                lower.95 = round(smoke_sig_cox[[8]][,3],3),
                                upper.95 = round(smoke_sig_cox[[8]][,4],3),
                                p.value=round(smoke_sig_cox[[7]][,5],3))
smoke_sig_cox_dat

#stage
Stage_sig_cox <- summary(coxph(formula=Surv(OS.time, OS)~pathologic_stage,
                               data=tcga_cox_datas))
Stage_sig_cox_dat <- data.frame(Names=rownames(Stage_sig_cox[[8]]),
                                HR = round(Stage_sig_cox[[7]][,2],3),
                                lower.95 = round(Stage_sig_cox[[8]][,3],3),
                                upper.95 = round(Stage_sig_cox[[8]][,4],3),
                                p.value=round(Stage_sig_cox[[7]][,5],3))
Stage_sig_cox_dat


#T.stage
T.stage_sig_cox <- summary(coxph(formula=Surv(OS.time, OS)~pathologic_T,
                                 data=tcga_cox_datas))
T.stage_sig_cox_dat <- data.frame(Names=rownames(T.stage_sig_cox[[8]]),
                                  HR = round(T.stage_sig_cox[[7]][,2],3),
                                  lower.95 = round(T.stage_sig_cox[[8]][,3],3),
                                  upper.95 = round(T.stage_sig_cox[[8]][,4],3),
                                  p.value=round(T.stage_sig_cox[[7]][,5],3))
T.stage_sig_cox_dat

#N.stage
N.stage_sig_cox <- summary(coxph(formula=Surv(OS.time, OS)~pathologic_N,
                                 data=tcga_cox_datas))
N.stage_sig_cox_dat <- data.frame(Names=rownames(N.stage_sig_cox[[8]]),
                                  HR = round(N.stage_sig_cox[[7]][,2],3),
                                  lower.95 = round(N.stage_sig_cox[[8]][,3],3),
                                  upper.95 = round(N.stage_sig_cox[[8]][,4],3),
                                  p.value=round(N.stage_sig_cox[[7]][,5],3))
N.stage_sig_cox_dat

# M.stage
M.stage_sig_cox <- summary(coxph(formula=Surv(OS.time, OS)~pathologic_M,
                                 data=tcga_cox_datas))
M.stage_sig_cox_dat <- data.frame(Names=rownames(M.stage_sig_cox[[8]]),
                                  HR = round(M.stage_sig_cox[[7]][,2],3),
                                  lower.95 = round(M.stage_sig_cox[[8]][,3],3),
                                  upper.95 = round(M.stage_sig_cox[[8]][,4],3),
                                  p.value=round(M.stage_sig_cox[[7]][,5],3))
M.stage_sig_cox_dat


#riskscore
riskscore_sig_cox<-summary(coxph(formula=Surv(OS.time, OS)~Riskscore,
                                 data=tcga_cox_datas))
riskscore_sig_cox_dat <- data.frame(Names=rownames(riskscore_sig_cox[[8]]),
                                    HR = round(riskscore_sig_cox[[7]][,2],3),
                                    lower.95 = round(riskscore_sig_cox[[8]][,3],3),
                                    upper.95 = round(riskscore_sig_cox[[8]][,4],3),
                                    p.value=round(riskscore_sig_cox[[7]][,5],3))
riskscore_sig_cox_dat

sig_cox_dat <- rbind(Age_sig_cox_dat,
                     Gender_sig_cox_dat,
                     T.stage_sig_cox_dat,
                     N.stage_sig_cox_dat,
                     M.stage_sig_cox_dat,
                     Stage_sig_cox_dat,
                     smoke_sig_cox_dat,
                     riskscore_sig_cox_dat)
data.sig <- data.frame(Names=sig_cox_dat$Names,
                       p.value=sig_cox_dat$p.value,
                       sig_cox_dat$HR,
                       sig_cox_dat$lower.95,
                       sig_cox_dat$upper.95)
data.sig <- crbind2DataFrame(data.sig)
data.sig$Names
rownames(data.sig) <- c("Age",
                        "Gender",
                        "pathologic_T",
                        "pathologic_N",
                        "pathologic_M",
                        "pathologic_Stage",
                        "tobacco_smoking_history",
                        "RiskScore")
data.sig$Names <- rownames(data.sig)
data.sig
data.sig$p.value=ifelse(data.sig$p.value<0.001,'<0.001',data.sig$p.value)
pdf('results/08.Nomogram/Univariate.pdf',height = 5,width = 7.5,onefile = F)
mg_forestplot_v2(data.sig,xlog = T,colgap =5,lineheight = 9,zero = 1,
                 boxsize = .2,lwd.zero=1,lwd.ci=1.5,lwd.xaxis=1,
                 box_col='blue',summary_col="black",lines_col='black',zero_col='grey',
                 xlab='Hazard Ratio',lty.ci = 6,graph.pos =4)
dev.off()



#pathologic_T+pathologic_N+pathologic_M
muti_sig_cox <- summary(coxph(formula=Surv(OS.time, OS)~Age+pathologic_T+pathologic_N+Riskscore,
                              data=tcga_cox_datas))
muti_cox_dat <- data.frame(Names=rownames(muti_sig_cox[[8]]),
                           HR = round(muti_sig_cox[[7]][,2],3),
                           lower.95 = round(muti_sig_cox[[8]][,3],3),
                           upper.95 = round(muti_sig_cox[[8]][,4],3),
                           p.value=round(muti_sig_cox[[7]][,5],3))
data.muti <- data.frame(Features=muti_cox_dat$Names,
                        p.value=muti_cox_dat$p.value,
                        muti_cox_dat$HR,
                        muti_cox_dat$lower.95,
                        muti_cox_dat$upper.95)
data.muti <- as.data.frame(data.muti)
data.muti
rownames(data.muti) <- c('Age','pathologic_T','pathologic_N','Riskscore')
data.muti$Features=rownames(data.muti)
data.muti
data.muti$p.value=ifelse(data.muti$p.value<0.001,'<0.001',data.muti$p.value)
pdf('results/08.Nomogram/Multivariate.pdf',height = 3,width =7.5,onefile = F)
mg_forestplot_v2(data.muti,xlog = T,colgap =5,lineheight = 9,zero = 1,
                 boxsize = .2,lwd.zero=1,lwd.ci=1.5,lwd.xaxis=1,
                 box_col='blue',summary_col="black",lines_col='black',zero_col='grey',
                 xlab='Hazard Ratio',lty.ci = 6,graph.pos =4)

dev.off()



mg_nomogram=function(clinical_riskscore,os,status,title='Nomogram',
                     quick=T,mks = c(1,3,5)){
  #clinical_riskscore=dat1[,3:5]
  #os=dat1[,1]
  #status=dat1[,2]
  #sum(is.na(norm.stat.al[,3]))
  norm.stat.al=data.frame(clinical_riskscore,time=os,status=status)
  norm.stat.al=as.data.frame(norm.stat.al)
  library(rms)
  env <- globalenv()
  env$MG_Grobal_DDSet <- rms::datadist(norm.stat.al) 
  options(datadist='MG_Grobal_DDSet')
  fmla <- as.formula(paste0("Surv(time, status) ~",paste0(colnames(clinical_riskscore),collapse = '+')))
  cox2 <- cph(fmla, data = norm.stat.al,surv = T,x = T,y = T)
  #summary(cox2)
  #surv=Survival(cox2)
  
  fp <- predict(cox2)
  cindex=getC_index(fp,norm.stat.al$time,norm.stat.al$status)
  #cindex.orig=1-rcorr.cens(fp,Surv(norm.stat.al$time,norm.stat.al$status))
  cut.time=c()
  if(quantile(os[!is.na(os)])['75%']<12){
    cut.time=mks
  }else if(quantile(os[!is.na(os)])['75%']<365){
    cut.time=c(12*mks[1],12*mks[2],12*mks[3])
  }else{
    cut.time=c(365*mks[1],365*mks[2],365*mks[3])
  }
  cut.time=cut.time[which(cut.time<quantile(os,seq(0,1,0.01))['100%'])]
  print(cut.time)
  #regplot(cox2)
  #  print(regplot(cox3
  #,observation=pbc[2,] 

  #              ,title=title
  #              ,failtime = cut.time
  #              ,prfail = TRUE 
  #              ,showP = T 
  #              ,droplines = F#
  #,colors = mg_colors[1:3] 
  #,rank="decreasing") 
  #,interval="confidence"
  #,rank="decreasing"
  #,clickable=T
  #              ,points=TRUE)) 
  
  #  plot(nom)
  surv=Survival(cox2)
  survs=list()
  cal_all=list()
  for(i in 1:length(cut.time)){
    f1<-cph(formula = fmla,data=norm.stat.al,x=T,y=T,surv = T,na.action=na.delete,time.inc = cut.time[i]) 
    cal1<-calibrate(f1, cmethod="KM", method="boot",u=cut.time[i],m=floor(sum(f1$n)/3)) 
    cal_all=c(cal_all,list(cal1))
    #    surv0 <- function(x)surv(cut.time[i],lp=x) 
    #    survs=c(survs,list(surv0))
  }
  #f2<-cph(formula = fmla,data=norm.stat.al,x=T,y=T,surv = T,na.action=na.delete,time.inc = cut.time[2]) 
  #cal3<-calibrate(f2, cmethod="KM", method="boot",u=cut.time[2],m=100,B=200) 
  #f3<-cph(formula = fmla,data=norm.stat.al,x=T,y=T,surv = T,na.action=na.delete,time.inc = cut.time[3]) 
  #cal5<-calibrate(f3, cmethod="KM", method="boot",u=cut.time[3],m=100,B=200) 
  if(length(cut.time)==1){
    surv1 <- function(x)surv(cut.time[1],lp=x) 
    survs=list(surv1)
  }else if(length(cut.time)==2){
    surv1 <- function(x)surv(cut.time[1],lp=x) 
    surv2 <- function(x)surv(cut.time[2],lp=x) 
    survs=list(surv1,surv2)
  }else if(length(cut.time)==3){
    surv1 <- function(x)surv(cut.time[1],lp=x) 
    surv2 <- function(x)surv(cut.time[2],lp=x) 
    surv3 <- function(x)surv(cut.time[3],lp=x) 
    survs=list(surv1,surv2,surv3)
  }
  nom=nomogram(cox2,fun=survs,lp= F
               ,funlabel=c(paste0(mks[1], '-Year Survival'),
                           paste0(mks[2], '-Year Survival'),
                           paste0(mks[3], '-Year Survival'))[1:length(cut.time)]
               ,maxscale=100
               ,fun.at=seq(0,1,0.2)
  )
  
  if(!quick){
    cal_all=list()
    for(i in 1:length(cut.time)){
      cal1=get_best_calibrate(cox2,cut.time[i])
      cal_all=c(cal_all,list(cal1))
    }
    #cal3=get_best_calibrate(cox2,cut.time[2])
    #cal5=get_best_calibrate(cox2,cut.time[3])
  }
  lay2 <- customLayout::lay_new(matrix(1:2))
  lay1 <- customLayout::lay_new(matrix(1:1))
  cl <- customLayout::lay_bind_col(lay2, lay1, widths = c(1, 1.5)) 
  #customLayout::lay_show(cl)
  customLayout::lay_set(cl) 
  
  plot(cal_all[[1]],lwd = 2,lty = 1,errbar.col = mg_colors[1], bty = "l",xlim = c(0,1),ylim= c(0,1),xlab = "Nomogram-prediced OS (%)"
       ,ylab = "Observed OS (%)",col = mg_colors[1],cex.lab=1.2,cex.axis=1, cex.main=1.2, cex.sub=0.6,subtitles = F)
  #lines(cal_all[[1]][,c('mean.predicted',"KM")],type = 'b', lwd = 1, col = mg_colors[1], pch = 16)
  mtext("")
  if(length(cal_all)>1){
    for(i in 2:length(cal_all)){
      plot(cal_all[[i]],lwd = 2,lty = 1,errbar.col = mg_colors[i],xlim = c(0,1),ylim= c(0,1),col = mg_colors[i],add = T)
      #lines(cal_all[[i]][,c('mean.predicted',"KM")],type = 'b', lwd = 1, col = mg_colors[i], pch = 16)
    }
    #plot(cal3,lwd = 2,lty = 0,errbar.col = mg_colors[3],xlim = c(0,1),ylim= c(0,1),col = mg_colors[3],add = T)
    #lines(cal3[,c('mean.predicted',"KM")],type = 'b', lwd = 1, col = mg_colors[3], pch = 16)
  }
  abline(0,1, lwd = 2, lty = 3, col = 'black')
  legend("topleft", legend = c(paste0(mks[1], '-Year'),
                               paste0(mks[2], '-Year'),
                               paste0(mks[3], '-Year'))[1:length(cut.time)],col =mg_colors[1:length(cut.time)],lwd = 2,cex = 1.2,bty = "n")
  
  fp <- predict(cox2)
  cindex=getC_index(fp,norm.stat.al$time,norm.stat.al$status)
  dca_dat=data.frame(Nomogram=fp,status=norm.stat.al$status)
  fp.al=cbind(fp)
  for(i in 1:ncol(clinical_riskscore)){
    fmla1 <- as.formula(paste0("Surv(time, status) ~",colnames(clinical_riskscore)[i]))
    cox1 <- cph(fmla1, data = norm.stat.al,surv = T,x = T,y = T)
    fp1 <- predict(cox1)
    fp.al=cbind(fp.al,fp1)
  }
  colnames(fp.al)=c('Nomogram',colnames(clinical_riskscore))
  fp.al=as.data.frame(fp.al)
  fp.al$status=norm.stat.al$status
  mg_plotDCA(fp.al$status
             ,c('Nomogram',colnames(clinical_riskscore))
             ,c('Nomogram',colnames(clinical_riskscore)),fp.al)
  #plot(cal1,xlim=c(0,1),ylim=c(0,1))
  #plot(cal3,xlim=c(0,1),ylim=c(0,1))
  #plot(cal5,xlim=c(0,1),ylim=c(0,1))
  plot(nom)
  options(datadist=NULL)
  return(list(Mod=cox2,Cindex=cindex,CutTime=cut.time))
}

getC_index=function(riskscore,os,status){
  inds=which(!is.na(riskscore)&!is.na(os)&!is.na(status))
  riskscore=riskscore[inds]
  os=os[inds]
  status=status[inds]
  c1=survcomp::concordance.index(x=riskscore, surv.time=os, surv.event=status,
                                 method="noether")
  #c2=concordance.index(x=riskscore[order(rnorm(length(riskscore)))], surv.time=os, surv.event=status,
  #                     method="noether")
  #p=min(cindex.comp(c1, c2)$p.value,cindex.comp(c2, c1)$p.value)  
  return(c1)
}
pdf('results/08.Nomogram/nomogram.pdf', width = 12, height = 10)
nom.plot=mg_nomogram(data.frame(RiskScore=tcga_cox_datas$Riskscore,
                                Age=tcga_cox_datas$Age,
                                pathologic_N=tcga_cox_datas$pathologic_N),
                     os = tcga_cox_datas$OS.time,
                     status = tcga_cox_datas$OS,quick=T,
                     mks = c(1,2,3),title = 'Nomogram')
dev.off()

mg_nomogram_buti(nom.plot$Mod,cut.time = c(1*365,2*365,3*365))


sum1=summary(coxph(formula=Surv(OS.time, OS)~Age+pathologic_N+Riskscore, data=tcga_cox_datas))
sum1$concordance[1]
c_index=as.numeric(sum1$concordance[1])

clinical.feature=c('Age','Gender','tobacco_smoking_history','pathologic_T','pathologic_N','pathologic_M','pathologic_stage','Riskscore')
for (i in 1:length(clinical.feature)) {
  as.formula(paste0("Surv(OS.time, OS) ~",clinical.feature[i]))
  sum_res=summary(coxph(as.formula(paste0("Surv(OS.time, OS) ~",clinical.feature[i])), data=tcga_cox_datas))
  index=as.numeric(sum_res$concordance[1])
  c_index=append(c_index,index)
}
names(c_index)=c('Nomogram',clinical.feature)
c_index[order(c_index)]


c_index.df=data.frame(Feature=names(c_index),c_index=unname(c_index))
c_index.df=c_index.df[order(c_index.df$c_index,decreasing = T),]
c_index.df

pdf('results/08.Nomogram/c_index.pdf',height = 4.5,width = 7)

ggplot(c_index.df, aes(x = reorder(Feature,c_index), y = c_index)) +
  geom_bar(stat = "identity", fill = "steelblue", alpha = 0.8) +
  geom_text(aes(label = sprintf("%.3f", c_index)), 
            vjust = -0.5, size = 3.5, color = "black") +
  labs(title = "C-index of Different Features",
       x = "Feature",
       y = "C-index") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1,color='black',size=12),
        axis.text.y = element_text(color='black',size=12),
        plot.title = element_text(hjust = 0.5),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.y = element_blank()) +
  scale_y_continuous(limits = c(0, 0.75), expand = c(0, 0))

dev.off()

save.image(file = 'project.RData')




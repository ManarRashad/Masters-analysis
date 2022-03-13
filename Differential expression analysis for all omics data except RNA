###########################################
#  - IHMP   Project            #
#  - Get Differential features            #
#  - Main function to be called           #
#  - 2019-16-08                            #
#  - Copyright: Mohamed Hamed/ Manar Rashad  #
###########################################

#R programming environments:
# R version 3.4.1 (2017-06-30)
# Platform: x86_64-apple-darwin15.6.0 (64-bit)
# Running under: OS X El Capitan 10.11.6
# locale: en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8

library(readr)
library(genefilter)
library(multtest)

correctPvalueandReturnAll<-function(tt.pval,method)
{
  mt=mt.rawp2adjp(tt.pval,proc=method)
  adjp=mt$adjp[order(mt$index),]
  return(adjp[,2])
}

load("pheno.RDATA")


############## Specify the pheno you want to work on:##############
############## choose one of the following two lines:##############

#To compare between the cntrl and prediabetic along all healthy visits  
pheno=pheno.healthy.all
# To compare between the IR and IS in the prediabetic only along all healthy visits  
pheno=pheno.healthy.pre

############### Specify the contrast /column in the pheno table which has the two groups of comparisons
############### Choose one of the follwoing two lines ############### 
contrast=c("Class")
contrast=c("IRIS")

load("clinical_Class.RDATA")
load("clinical_IRIS.RDATA")

groups=as.character(unlist(unique(pheno[contrast])))
groups

file_directory="~/Dropbox/Rostock/projects/IHMP/Processed_Files"
datasets=c("RNAseq_abundance","cytokine_abundance","gut_16s_abundance","nares_16s_abundance",
           "proteome_abundance","metabolome_abundance","clinical_tests")

# dataset="RNAseq_abundance"
for (dataset in datasets)
{
  filepath=file.path(file_directory,dataset,paste(dataset,".txt",sep=""))
  df<- read_delim(filepath,"\t", escape_double = FALSE, trim_ws = TRUE)
  
  # get only samples which are common in the  pheno (healthy and iris) amd this specific omics dataset
  common.samples=intersect(pheno$SampleID, df$SampleID)
  pheno.df=pheno[pheno$SampleID %in% common.samples,]
  df=df[df$SampleID %in% common.samples, ]
  
  # remove duplicate samples
  df=df[!duplicated(df$SampleID),]
  sum(duplicated(df$SampleID))
  
  # remove the last 6 columns in the"proteome_abundance" and "metabolome_abundance"
  if(dataset %in% c("proteome_abundance","metabolome_abundance","clinical_tests"))
  {
    df=df[- c((dim(df)[2] -5) : dim(df)[2])]
  }
  
  
  # transpose and remove the sample ids
  df.samples=df$SampleID
  df=df[-1]
  df=t(df)
  colnames(df)=df.samples
  
  # check for complete cases in  clinical _tests
  if(dataset=="clinical_tests"){
    df=df[complete.cases(df),]
    #df[is.na(df)]=0  
    df=gsub(">","",df)
    df=gsub("<","",df)
    features=rownames(df)
    df=apply(df,2,as.numeric)
    rownames(df)=features
    }
  
  print(paste(dataset,">>>>",range(df)),sep="")
  
  
  # normalize the data only of the in case of cytokines and metabolome 
  if(dataset %in% c("metabolome_abundance","cytokine_abundance"))
  {
    df=log2(df+1)
  }

  samples.case=unlist(pheno.df[pheno.df[contrast]==groups[2],]$SampleID)
  samples.ctrl=unlist(pheno.df[pheno.df[contrast]==groups[1],]$SampleID)
  df.case=df[,colnames(df) %in%samples.case ]
  df.ctrl=df[,colnames(df) %in%samples.ctrl ]
  
  
  exp=cbind(df.case,df.ctrl)
  
  case.indecies=c(1: dim(df.case)[2])
  ctrl.indecies=c( (dim(df.case)[2] +1) : dim(exp)[2]) #### or  seq(from=9,to=dim(exp)[2])
  
  #make the rows of pheno table in the same order as the columns of the exp
  rownames(pheno.df)=pheno.df$SampleID
  pheno.df=pheno.df[colnames(exp),]
  
  
  ## calculating LFC
  lfc.diff=apply(exp,1, function(x)  mean(x[case.indecies]) -mean(x[ctrl.indecies]))
  x=as.data.frame(lfc.diff)
  
  ## calcualting p values
  f=factor( c( rep(1, length(case.indecies)) , rep(2, length(ctrl.indecies)) ))
  t.pval=rowttests(as.matrix(exp),f)$p.value
  t.pval.adj=correctPvalueandReturnAll(t.pval,"BH")
  
  res=cbind(lfc.diff,t.pval,t.pval.adj)
  
  #####  selection criteria for identifying DEGS : final crietria : pval<0.05 only 
  hits.res=res[t.pval<0.05,]  ##### identify DEGs based on teh significance level only
  # hits.res=res[t.pval.adj<0.05,]  ##### identify DEGs based on teh significance level only
  # hits.res=res[abs(lfc.diff) > log2(1.2),]  ##### identify DEGs based on the LFC only
  # hits.res=res[abs(lfc.diff) > log2(1.2)   & t.pval <0.05,]  ##### identify DEGs based on both  LFC and the significane level 
  hits.data=exp[rownames(hits.res),]
  hits.anno=c(  rep(groups[2],length(case.indecies)) , rep(groups[1],length(ctrl.indecies)) )
  
  # save the whole data in a list 
  l=list(res=as.data.frame(res),data=exp, hits.res=hits.res,hits.data=hits.data,hits.anno=hits.anno,this.pheno=pheno.df,contrast=contrast,groups=groups)
  # l=list(hits.res=hits.res,hits.data=hits.data,hits.anno=hits.anno,this.pheno=pheno.df,contrast=contrast,groups=groups)
  dataset.name= paste(as.character(strsplit(dataset,"_")[[1]][1]) ,contrast,sep="_")
  assign(dataset.name,l)
  save( list=dataset.name , file=paste(dataset.name,".RDATA",sep=""))
}

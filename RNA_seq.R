###########################################
#  - IHMP   Project                       #
#  - Get Differential features of RNA-Seq #
#  - Main function to be called           #
#  - 2019-16-08                           #
#  - Copyright: Mohamed Hamed             #
###########################################

#R programming environments:
# R version 3.4.1 (2017-06-30)
# Platform: x86_64-apple-darwin15.6.0 (64-bit)
# Running under: OS X El Capitan 10.11.6
# locale: en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8

library(readr)
library("vsn")
library(DESeq2)
library(pheatmap)


load("pheno.RDATA")

# To compare between the cntrl and prediabetic along all healthy visits  
pheno=pheno.healthy.all
# To compare between the IR and IS in the prediabetic only along all healthy visits  
pheno=pheno.healthy.pre

############### Specify the contrast /column in the pheno table which has the two groups of comparisons. choose only one line
contrast=c("Class")
contrast=c("IRIS")

groups=as.character(unlist(unique(pheno[contrast])))
groups


file_directory="~/Dropbox/Rostock/projects/IHMP/Processed_Files"
dataset= "RNAseq_abundance"
filepath=file.path(file_directory,dataset,paste(dataset,".txt",sep=""))
df<- read_delim(filepath,"\t", escape_double = FALSE, trim_ws = TRUE)
  
# get only samples which are common in the  pheno (healthy and iris) amd this specific omics dataset
common.samples=intersect(pheno$SampleID, df$SampleID)
pheno.df=pheno[pheno$SampleID %in% common.samples,]
pheno.df["Contrast.column"]=pheno.df[contrast] ## to be used later in DEseq package
df=df[df$SampleID %in% common.samples, ]
  
# remove duplicate samples, but here there are no duplicated
df=df[!duplicated(df$SampleID),]
sum(duplicated(df$SampleID))

  
# transpose and remove the sample ids
df.samples=df$SampleID
df=df[-1]
df=as.matrix(df)
df=t(df)
colnames(df)=df.samples
genes=rownames(df)  

# return the data of the RNAseq to its raw formats and coerce as integer ..a requirement for deseq2 package: log2(n+1) as stated in the paper
df=(2^df)-1
df=apply(df,2,as.integer)
rownames(df)=genes

# sort the df columns to be first the case group , then the cntrl group
# this is a useless step here with deseq, but it is just for my personal comfort. @MH
samples.case=unlist(pheno.df[pheno.df[contrast]==groups[2],]$SampleID)
samples.ctrl=unlist(pheno.df[pheno.df[contrast]==groups[1],]$SampleID)
df.case=df[,colnames(df) %in%samples.case ]
df.ctrl=df[,colnames(df) %in%samples.ctrl ]
df=cbind(df.case,df.ctrl)


#It is absolutely critical that the columns of the "df" and the rows of the "pheno.df" (information about samples) are in the same order. DESeq2 will not make guesses as to which column of the count matrix belongs to which row of the column data, these must be provided to DESeq2 already in consistent order. 
#put the "pheno.df" rows in the same order as columns of the "df"
rownames(pheno.df)=pheno.df$SampleID
#pheno.df=pheno.df[,"Contrast.column",drop=F]
pheno.df=pheno.df[colnames(df),]
match(colnames(df),rownames(pheno.df))
all(colnames(df)== rownames(pheno.df))

#### this is hardcoded..we wanna all genric ..don't do this..
# samples.pre=unlist(pheno.df[pheno.df$Class=="Prediabetic",]$SampleID)
# samples.ctrl=unlist(pheno.df[pheno.df$Class=="Control",]$SampleID)
# #for IRIS
# samples.pre.IR=unlist(pheno.df[pheno.df$IRIS=="IR",]$SampleID)
# samples.pre.IS=unlist(pheno.df[pheno.df$IRIS=="IS",]$SampleID)

# samples.case=unlist(pheno.df[pheno.df[contrast]==groups[2],]$SampleID)
# samples.ctrl=unlist(pheno.df[pheno.df[contrast]==groups[1],]$SampleID)
# df.case=df[,colnames(df) %in%samples.case ]
# df.ctrl=df[,colnames(df) %in%samples.ctrl ]
# exp=cbind(df.case,df.ctrl)

### same ..don't write any hard codes... we specify the comparison using the pheno type dataframe and the contrast column in the start of the script
# df.ctrl=df[,colnames(df) %in%samples.ctrl ]
# df.pred=df[,colnames(df) %in%samples.pre ]
# #for IRIS
# df.pred.IR =df[,colnames(df) %in%samples.pre.IR ]
# df.pred.IS =df[,colnames(df) %in%samples.pre.IS ]
#   exp=cbind(df.pred,df.ctrl)
# # # for IRIS
# # exp=cbind(df.pred.IR,df.pred.IS)
 
# pheno.class = pheno.df[,c("SampleID","Class")]#in case of class
# pheno.IRIS = pheno.df[,c("SampleID","IRIS")]#incase of IRIS
# case.indecies=c(1: dim(df.pred.IS)[2])
# ctrl.indecies=c( (dim(df.pred.IS)[2] +1) : dim(exp)[2])
  
###### DO the differential EXP analysis using DeSeq2
exp=df  
dds = DESeqDataSetFromMatrix( countData = exp , colData = pheno.df , design = ~ Contrast.column)
# dds = DESeqDataSetFromMatrix( countData = exp , colData = pheno.IRIS , design = ~ IRIS)
dds.run = DESeq(dds)
### direct results or specifying teh contrast (to make a res object based on two specific conditions/treatment)
res=results(dds.run)
res=results(dds.run, contrast = c("Contrast.column",groups[2] ,groups[1] ) )
  
#remove nulls
res.complete=res[complete.cases(res), ]
summary(res.complete)

all.res=as.data.frame(res.complete)
#plotMA(res, ylim=c(-1,1)) 
#summary (res)
#####  selection criteria for identifying DEGS : final crietria : pval<0.05 and logFC >2 
hits.res=all.res[all.res$padj< 0.05,] ### identify DEGs based on the adjusted pvalue
# hits.res=all.res[all.res$pvalue< 0.05 & abs(all.res$log2FoldChange)>log2(2),] ##### identify DEGs based on both  LFC and the significane level
# hits.res=all.res[all.res$pvalue< 0.05 & abs(all.res$log2FoldChange)>log2(1.5),] ##### identify DEGs based on both  LFC and the significane level
# hits.res=all.res[all.res$pvalue< 0.05,]##### identify DEGs based on the significance level only
# hits.res=all.res[all.res$padj< 0.05 & abs(all.res$log2FoldChange)>log2(2),]
# hits.res=all.res[all.res$padj< 0.05 & abs(all.res$log2FoldChange)>log2(1.2),]
#write.table(res.degs, file = "res.degs_deseq_class.txt", sep = "\t", quote = FALSE)

hits.res=hits.res[order(hits.res$padj),]  
# get the normalized/processed data  
ntd=normTransform(dds)
all.data= assay(ntd)  # or
#all.data2=log2(exp+1)
hits.data= all.data[ rownames(hits.res), ]
hits.anno=c(  rep(groups[2],dim(df.case)[2]) , rep(groups[1],dim(df.ctrl)[2]) )

# save the whole data in a list 
l=list(res=as.data.frame(res),data=all.data, hits.res=hits.res,hits.data=hits.data,hits.anno=hits.anno,this.pheno=pheno.df,contrast=contrast,groups=groups)
dataset.name= paste("DESeq" ,contrast,sep="_")
assign(dataset.name,l)
save( list=dataset.name , file=paste(dataset.name,".RDATA",sep=""))

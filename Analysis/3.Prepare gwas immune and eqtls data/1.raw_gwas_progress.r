library(readr)
library(stringr)
fl<-list.files('/share/pub/dengcy/Cancer_Gwas/CollectedData/Gwas_data/finngen_R7/rawdata')
setwd("/share/pub/dengcy/Cancer_Gwas/CollectedData/Gwas_data/finngen_R7/rawdata")
flx<-lapply(fl,function(x){
    GWAS_raw <-read_table(x)
    GWAS_raw<-GWAS_raw[,c(1,2,3,4,5,6,7,9,10,11,12,13)]
    colnames(GWAS_raw)<-c("chrom","pos","REF","ALT","rsid","nearest_genes","pval","beta","se","maf")
    x<-tolower(x)
    x<-str_replace(x,".gz","")
    write.table(GWAS_raw,file=paste("/share/pub/dengcy/Cancer_Gwas/CollectedData/Gwas_data/finngen_R7/",x,".txt",sep=""),row.names=F,quote=F)
    retrun(x)
    })

ab<-unlist(flx)
ab<-tolower(ab)
ab<-data.frame(ab)
write.table(ab,file="/share/pub/dengcy/Cancer_Gwas/CollectedData/Gwas_data/finngen_R7/finnngen_names_files.txt",quote=F,row.names=F,col.names=F)

library(readr)
setwd("/share/pub/dengcy/Cancer_Gwas/CollectedData/Gwas_data/IEU_opengwas_project/")

##bbj-a-140
Pancreas<-"bbj-a-140.vcf.gz"
#`#CHROM`    POS ID     REF   ALT   QUAL  FILTER INFO   FORMAT  `ieu-a-1008`  
GWAS_raw <-read_table(Pancreas,skip=109)
colnames(GWAS_raw )<-c("CHROM","POS","ID","REF","ALT","QUAL","FILTER","INFO","FORMAT","IEU")
#ES:SE:LP:AF:SS:ID
#ieu-a格式数据
if( GWAS_raw$FORMAT[1]=="ES:SE:LP:AF:ID"){
 value_list<-lapply(1:nrow(GWAS_raw),function(i){
  index <- strsplit(GWAS_raw$FORMAT[i], split = ":",fixed=T)
  index<-unlist(index)
  value <- strsplit(GWAS_raw$IEU[i], split = ":",fixed=T)
  value<-unlist(value)
  names(value)<-index
    
  value<-value[c("ES","SE","LP","AF")]
  return(value)
  }) 

 value_df<- matrix(as.numeric(unlist(value_list)),byrow=T,ncol=4)
 value_df<-as.data.frame(value_df)
 colnames(value_df)<-c("ES","SE","LP","AF")
 gwas_data<-GWAS_raw[,c(1:5)]
 gwas_data<- cbind(gwas_data,value_df)
 rm(value_df,GWAS_raw)
 gwas_data$LP <-  10^(-gwas_data$LP)
 gwas_data$LP[which(gwas_data$LP>1)]<-1
 #colnames(gwas_data)
 colnames(gwas_data)<-c("chrom","pos","rsid","REF","ALT","beta","se","p","maf")
 }
# if(is.na(gwas_data$maf[1])) 
    # gwas_data$maf<-0.1
write.table(gwas_data,file="/share/pub/dengcy/Cancer_Gwas/CollectedData/Gwas_data/IEU_opengwas_project/Pancreas_EastAsian_gwas_data.txt",row.names=F,quote=F)
gwas_data <-read_table("/share/pub/dengcy/Cancer_Gwas/CollectedData/Gwas_data/IEU_opengwas_project/Pancreas_EastAsian_gwas_data.txt")

gwas_data<-gwas_data[,c("rsid","chrom","pos","REF","ALT","beta","se","p","maf")]
colnames(gwas_data)<-c("SNP","chrom","pos","A1","A2","se","beta","pval","maf")
write.table(gwas_data,file="/share/pub/dengcy/Cancer_Gwas/Runtime1.0/3.LDSC/Ieu_gwasdata/Pancreas_EastAsian_gwas_LDSC.txt",row.names=F,quote=F)





##"ieu-a-822"
Pancreas<-"ieu-a-822.vcf.gz"
#`#CHROM`    POS ID     REF   ALT   QUAL  FILTER INFO   FORMAT  `ieu-a-1008`  
GWAS_raw <-read_table(Pancreas,skip=109)
colnames(GWAS_raw )<-c("CHROM","POS","ID","REF","ALT","QUAL","FILTER","INFO","FORMAT","IEU")
#ES:SE:LP:AF:SS:ID
#ieu-a格式数据
if( GWAS_raw$FORMAT[1]=="ES:SE:LP:SS:ID"){
 value_list<-lapply(1:nrow(GWAS_raw),function(i){
  index <- strsplit(GWAS_raw$FORMAT[i], split = ":",fixed=T)
  index<-unlist(index)
  value <- strsplit(GWAS_raw$IEU[i], split = ":",fixed=T)
  value<-unlist(value)
  names(value)<-index
    
  value<-value[c("ES","SE","LP","SS")]
  return(value)
  }) 

 value_df<- matrix(as.numeric(unlist(value_list)),byrow=T,ncol=4)
 value_df<-as.data.frame(value_df)
 colnames(value_df)<-c("ES","SE","LP","SS")
 gwas_data<-GWAS_raw[,c(1:5)]
 gwas_data<- cbind(gwas_data,value_df)
 rm(value_df,GWAS_raw)
 gwas_data$LP <-  10^(-gwas_data$LP)
 gwas_data$LP[which(gwas_data$LP>1)]<-1
 #colnames(gwas_data)

 colnames(gwas_data)<-c("chrom","pos","rsid","REF","ALT","beta","se","p","N")
  gwas_data$maf<-0.1
 }
# if(is.na(gwas_data$maf[1])) 
    # 
write.table(gwas_data,file="/share/pub/dengcy/Cancer_Gwas/CollectedData/Gwas_data/Pancreas_European_gwas_data.txt",row.names=F,quote=F)

gwas_data<-gwas_data[,c("rsid","chrom","pos","REF","ALT","beta","se","pval","maf")]
colnames(gwas_data)<-c("SNP","chrom","pos","A1","A2","beta","se","pval","maf")
write.table(gwas_data,file="/share/pub/dengcy/Cancer_Gwas/Runtime1.0/3.LDSC/Ieu_gwasdata/Pancreas_European_gwas_LDSC.txt",row.names=F,quote=F)
################################################################
#Lung
########################
library(readr)
setwd("/share/pub/dengcy/Cancer_Gwas/CollectedData/Gwas_data/raw_files/")

##"ieu-b-4954" European
Lung<-"ieu-b-4954.vcf.gz"
#`#CHROM`    POS ID     REF   ALT   QUAL  FILTER INFO   FORMAT  `ieu-a-1008`  
GWAS_raw <-read_table(Lung,skip=109)
colnames(GWAS_raw )<-c("CHROM","POS","ID","REF","ALT","QUAL","FILTER","INFO","FORMAT","IEU")
#ES:SE:LP:AF:SS:ID
#ieu-a格式数据
#if( GWAS_raw$FORMAT[1]=="ES:SE:LP:AF:ID"){
 value_list<-lapply(1:nrow(GWAS_raw),function(i){
  index <- strsplit(GWAS_raw$FORMAT[i], split = ":",fixed=T)
  index<-unlist(index)
  value <- strsplit(GWAS_raw$IEU[i], split = ":",fixed=T)
  value<-unlist(value)
  #names(value)<-index
    
  value<-value[1:4]
  return(value)
  }) 

 value_df<- matrix(as.numeric(unlist(value_list)),byrow=T,ncol=4)
 value_df<-as.data.frame(value_df)
 colnames(value_df)<-c("ES","SE","LP","AF")
 gwas_data<-GWAS_raw[,c(1:5)]
 gwas_data<- cbind(gwas_data,value_df)
 rm(value_df,GWAS_raw)
 gwas_data$LP <-  10^(-gwas_data$LP)
 gwas_data$LP[which(gwas_data$LP>1)]<-1
 #colnames(gwas_data)
 colnames(gwas_data)<-c("chrom","pos","rsid","REF","ALT","beta","se","p","maf")
 
# if(is.na(gwas_data$maf[1])) 
    # gwas_data$maf<-0.1
write.table(gwas_data,file="/share/pub/dengcy/Cancer_Gwas/CollectedData/Gwas_data/Lung_European_gwas_data.txt",row.names=F,quote=F)
gwas_data<-gwas_data[,c("rsid","chrom","pos","REF","ALT","beta","se","pval","maf")]
colnames(gwas_data)<-c("SNP","chrom","pos","A1","A2","beta","se","pval","maf")
write.table(gwas_data,file="/share/pub/dengcy/Cancer_Gwas/Runtime1.0/3.LDSC/Ieu_gwasdata/Lung_European_gwas_gwas_LDSC.txt",row.names=F,quote=F)



##"#"bbj-a-133" East Asian"
Lung<-"bbj-a-133.vcf.gz"
#`#CHROM`    POS ID     REF   ALT   QUAL  FILTER INFO   FORMAT  `ieu-a-1008`  
GWAS_raw <-read_table(Lung,skip=109)
colnames(GWAS_raw )<-c("CHROM","POS","ID","REF","ALT","QUAL","FILTER","INFO","FORMAT","IEU")

 value_list<-lapply(1:nrow(GWAS_raw),function(i){
  index <- strsplit(GWAS_raw$FORMAT[i], split = ":",fixed=T)
  index<-unlist(index)
  value <- strsplit(GWAS_raw$IEU[i], split = ":",fixed=T)
  value<-unlist(value)
  names(value)<-index
    
  value<-value[c("ES","SE","LP","SS")]
  return(value)
  }) 

 value_df<- matrix(as.numeric(unlist(value_list)),byrow=T,ncol=4)
 value_df<-as.data.frame(value_df)
 colnames(value_df)<-c("ES","SE","LP","SS")
 gwas_data<-GWAS_raw[,c(1:5)]
 gwas_data<- cbind(gwas_data,value_df)
 rm(value_df,GWAS_raw)
 gwas_data$LP <-  10^(-gwas_data$LP)
 gwas_data$LP[which(gwas_data$LP>1)]<-1
 #colnames(gwas_data)

 colnames(gwas_data)<-c("chrom","pos","rsid","REF","ALT","beta","se","p","N")
  gwas_data$maf<-0.1
 
write.table(gwas_data,file="/share/pub/dengcy/Cancer_Gwas/CollectedData/Gwas_data/Lung_EastAsian_gwas_data.txt",row.names=F,quote=F)

################################################################
#CRC
########################
library(readr)
setwd("/share/pub/dengcy/Cancer_Gwas/CollectedData/Gwas_data/raw_files/")

##"ieu-b-4965" European
CRC<-"ieu-b-4965.vcf.gz"
#`#CHROM`    POS ID     REF   ALT   QUAL  FILTER INFO   FORMAT  `ieu-a-1008`  
GWAS_raw <-read_table(CRC,skip=109)
colnames(GWAS_raw )<-c("CHROM","POS","ID","REF","ALT","QUAL","FILTER","INFO","FORMAT","IEU")
#ES:SE:LP:AF:SS:ID
#ieu-a格式数据
#if( GWAS_raw$FORMAT[1]=="ES:SE:LP:AF:ID"){
 value_list<-lapply(1:nrow(GWAS_raw),function(i){
  index <- strsplit(GWAS_raw$FORMAT[i], split = ":",fixed=T)
  index<-unlist(index)
  value <- strsplit(GWAS_raw$IEU[i], split = ":",fixed=T)
  value<-unlist(value)
  #names(value)<-index
    
  value<-value[1:4]
  return(value)
  }) 

 value_df<- matrix(as.numeric(unlist(value_list)),byrow=T,ncol=4)
 value_df<-as.data.frame(value_df)
 colnames(value_df)<-c("ES","SE","LP","AF")
 gwas_data<-GWAS_raw[,c(1:5)]
 gwas_data<- cbind(gwas_data,value_df)
 rm(value_df,GWAS_raw)
 gwas_data$LP <-  10^(-gwas_data$LP)
 gwas_data$LP[which(gwas_data$LP>1)]<-1
 #colnames(gwas_data)
 colnames(gwas_data)<-c("chrom","pos","rsid","REF","ALT","beta","se","p","maf")
 
# if(is.na(gwas_data$maf[1])) 
    # gwas_data$maf<-0.1
write.table(gwas_data,file="/share/pub/dengcy/Cancer_Gwas/CollectedData/Gwas_data/CRC_European_gwas_data.txt",row.names=F,quote=F)

##"bbj-a-107" 
CRC<-"bbj-a-107.vcf.gz"
#`#CHROM`    POS ID     REF   ALT   QUAL  FILTER INFO   FORMAT  `ieu-a-1008`  
GWAS_raw <-read_table(CRC,skip=109)
colnames(GWAS_raw )<-c("CHROM","POS","ID","REF","ALT","QUAL","FILTER","INFO","FORMAT","IEU")
#ES:SE:LP:AF:SS:ID
#ieu-a格式数据
#if( GWAS_raw$FORMAT[1]=="ES:SE:LP:AF:ID"){
 value_list<-lapply(1:nrow(GWAS_raw),function(i){
  index <- strsplit(GWAS_raw$FORMAT[i], split = ":",fixed=T)
  index<-unlist(index)
  value <- strsplit(GWAS_raw$IEU[i], split = ":",fixed=T)
  value<-unlist(value)
  #names(value)<-index
    
  value<-value[1:4]
  return(value)
  }) 

 value_df<- matrix(as.numeric(unlist(value_list)),byrow=T,ncol=4)
 value_df<-as.data.frame(value_df)
 colnames(value_df)<-c("ES","SE","LP","AF")
 gwas_data<-GWAS_raw[,c(1:5)]
 gwas_data<- cbind(gwas_data,value_df)
 rm(value_df,GWAS_raw)
 gwas_data$LP <-  10^(-gwas_data$LP)
 gwas_data$LP[which(gwas_data$LP>1)]<-1
 #colnames(gwas_data)
 colnames(gwas_data)<-c("chrom","pos","rsid","REF","ALT","beta","se","p","maf")
 gwas_data<-na.omit(gwas_data)
# if(is.na(gwas_data$maf[1])) 
    # gwas_data$maf<-0.1
write.table(gwas_data,file="/share/pub/dengcy/Cancer_Gwas/CollectedData/Gwas_data/CRC_EastAsian_gwas_data.txt",row.names=F,quote=F)

############################################
#breast cancer
########################

#ieu-b-4810 European
CRC<-"ieu-b-4810.vcf.gz"
#`#CHROM`    POS ID     REF   ALT   QUAL  FILTER INFO   FORMAT  `ieu-a-1008`  
GWAS_raw <-read_table(CRC,skip=109)
colnames(GWAS_raw )<-c("CHROM","POS","ID","REF","ALT","QUAL","FILTER","INFO","FORMAT","IEU")
#ES:SE:LP:AF:SS:ID
#ieu-a格式数据
#if( GWAS_raw$FORMAT[1]=="ES:SE:LP:AF:ID"){
 value_list<-lapply(1:nrow(GWAS_raw),function(i){
  index <- strsplit(GWAS_raw$FORMAT[i], split = ":",fixed=T)
  index<-unlist(index)
  value <- strsplit(GWAS_raw$IEU[i], split = ":",fixed=T)
  value<-unlist(value)
  #names(value)<-index
    
  value<-value[1:4]
  return(value)
  }) 

 value_df<- matrix(as.numeric(unlist(value_list)),byrow=T,ncol=4)
 value_df<-as.data.frame(value_df)
 colnames(value_df)<-c("ES","SE","LP","AF")
 gwas_data<-GWAS_raw[,c(1:5)]
 gwas_data<- cbind(gwas_data,value_df)
 rm(value_df,GWAS_raw)
 gwas_data$LP <-  10^(-gwas_data$LP)
 gwas_data$LP[which(gwas_data$LP>1)]<-1
 #colnames(gwas_data)
 colnames(gwas_data)<-c("chrom","pos","rsid","REF","ALT","beta","se","p","maf")
 
# if(is.na(gwas_data$maf[1])) 
    # gwas_data$maf<-0.1
write.table(gwas_data,file="/share/pub/dengcy/Cancer_Gwas/CollectedData/Gwas_data/Breast_European_gwas_data.txt",row.names=F,quote=F)

########################
library(readr)
setwd("/share/pub/dengcy/Cancer_Gwas/CollectedData/Gwas_data/IEU_opengwas_project")

#
a<-"ebi-a-GCST006464.vcf.gz"
#`#CHROM`    POS ID     REF   ALT   QUAL  FILTER INFO   FORMAT  `ieu-a-1008`  
GWAS_raw <-read_table(a,skip=109)
colnames(GWAS_raw )<-c("CHROM","POS","ID","REF","ALT","QUAL","FILTER","INFO","FORMAT","IEU")


###########################
#GCST007728 GCST90041886 
library(readr)
setwd("/share/pub/dengcy/Cancer_Gwas/CollectedData/Gwas_data/gwas_catalog")
a<-"GCST90041886_buildGRCh37.tsv.gz"
GWAS_raw <-read_table(a)
GWAS_raw<-GWAS_raw[,c(1:6,7,9,11)]
colnames(GWAS_raw)<-c("chrom","rsid","pos","REF","ALT","N","maf","se","beta")
write.table(GWAS_raw,file="/share/pub/dengcy/Cancer_Gwas/CollectedData/Gwas_data/GCST007728.txt",row.names=F,quote=F)

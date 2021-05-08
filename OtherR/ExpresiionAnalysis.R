##' convert gene ID among KH,KY,gene name, symbol,
##'
##'
##' @title Gff2GeneTable
##' @param gffFile GFF file
##' @param compress compress file or not
##' @return file save.
##' @export
##' @author Chen Zaohuang


# Expression Analysis

#allGeneList$operon<-ifelse(allGeneList$gene %in% operonID,"operon","nonoperon")
#allGeneList$HH<-ifelse(allGeneList$gene %in% ci.HH.genes,"H2H","")

#allGeneList$HHOperon<-ifelse(allGeneList$gene %in% ci.HHOperon.genes,"H2HOperon","")

#allGeneList %>% filter(operon=="operon" ) %>% select(1:6,8)->aa1
#allGeneList %>% filter(operon=="nonoperon") %>% select(1:6,8)->aa2
#allGeneList %>% filter(HH=="H2H") %>% select(1:5,CV,HH)->aa3
#allGeneList %>% filter(operon=="nonoperon" & HH!="H2H" ) %>% select(1:5,CV)->aa4
#aa4$type<-"Non-H2H non-operon"


#colnames(aa2)[6]<-"type"
#colnames(aa1)[6]<-"type"
#colnames(aa4)[7]<-"type"
#colnames(aa3)[7]<-"type"


#rbind(aa1,aa3,aa4,aa2)->aa
#aa$condition<-factor(aa$condition,levels = c("0h","24h","48h","96h","120h","144h"))
#aa$type<-factor(aa$type,levels = c("operon","H2H","Non-H2H non-operon","nonoperon"))
#ggplot(aa,aes(x=condition,y=log2(mean+1),fill=type))+geom_boxplot(notch = T)+
#  theme_bw()+guides(colour=guide_legend(title = ""))+theme(panel.grid.major = element_blank(),panel.background = element_blank(),panel.grid.minor = element_blank(),legend.position = "none",text = element_text(size = 20)) +
#  xlab(label = "")+ylab("Gene Expression (log2 TPM)")

#ggplot(aa,aes(x=condition,y=CV,fill=type))+geom_boxplot(notch = T)+
#  theme_bw()+guides(colour=guide_legend(title = ""))+theme(panel.grid.major = element_blank(),panel.background = element_blank(),panel.grid.minor = element_blank(),legend.position = "none",text = element_text(size = 20)) +
#  xlab(label = "")+ylab("Expression Variability (CV of TPM)")




#ggplot(aa,aes(x=condition,y=log2(mean+1),group=type,color=type))+stat_summary(fun.y = mean,geom = "line",position = position_dodge(width = 0.2)) +
#  stat_summary(fun.data = mean_cl_boot, geom = "errorbar",position = position_dodge(width = 0.2),width=0.7)+
##  theme_bw()+guides(colour=guide_legend(title = ""))+stat_summary(fun.y = mean,geom = "point",position = position_dodge(width = 0.2),size=3)+theme(panel.grid.major = element_blank(),panel.background = element_blank(),panel.grid.minor = element_blank(),text = element_text(size = 20),legend.position = "none") + xlab(label = "")+ylab("Mean Expression (TPM)")

#ggplot(aa,aes(x=condition,y=mean,group=type,color=type))+stat_summary(fun.y = mean,geom = "line",position = position_dodge(width = 0.2)) +
 # stat_summary(fun.data = mean_cl_boot, geom = "errorbar",position = position_dodge(width = 0.2),width=0.7)+
 # theme_bw()+guides(colour=guide_legend(title = ""))+stat_summary(fun.y = mean,geom = "point",position = position_dodge(width = 0.2),size=3)+theme(panel.grid.major = element_blank(),panel.background = element_blank(),panel.grid.minor = element_blank(),text = element_text(size = 20),) + xlab(label = "")+ylab("Mean Expression (TPM)")



#aa$mean<-log2(aa$mean+1)
#mean_normalizaed<-function(x){
#y<- ( x-mean(x))/(max(x)-min(x))

#  return(y)
#  }

#bb<-aggregate(aa$mean,list(aa$gene,aa$type),mean_normalizaed)
#as.numeric(bb$x)->aa$mean


#res$HH<-ifelse(res$gene %in% ci.HH.genes,"H2H","")
#res %>% filter(HH=="H2H") %>% select(1:6,8,HH)->aa1
#res %>% filter(operon=="operon") %>% select(1:8)->aa2
#res %>% filter(operon!="operon" & HH !="H2H") %>% select(1:6,8)->aa3
#aa3$type<-"Non-H2H non-operon"

#colnames(aa1)[8]<-"type"
#colnames(aa2)[7]<-"type"

#rbind(aa1,aa2,aa3)->aa
#aa$condition<-factor(aa$condition,levels = c("0h","24h","48h","96h","120h","144h"))
#aa$type<-factor(aa$type,levels = c("operon","H2H","Non-H2H non-operon"))
#ggplot(aa,aes(x=condition,y=log2FoldChange,group=type,color=type))+stat_summary(fun.y = mean,geom = "line",position = position_dodge(width = 0.2)) +
 # stat_summary(fun.data = mean_cl_boot, geom = "errorbar",position = position_dodge(width = 0.2),width=0.7)+
#  theme_bw()+guides(colour=guide_legend(title = ""))+stat_summary(fun.y = mean,geom = "point",position = position_dodge(width = 0.2),size=3)+theme(panel.grid.major = element_blank(),panel.background = element_blank(),panel.grid.minor = element_blank(),text = element_text(size = 20),legend.position = "none") + xlab(label = "")+ylab("Mean LFC")

#ggplot(aa,aes(x=condition,y=log2FoldChange,fill=type))+geom_boxplot(notch = T)+
#  theme_bw()+guides(colour=guide_legend(title = ""))+theme(panel.grid.major = element_blank(),panel.background = element_blank(),panel.grid.minor = element_blank(),legend.position = "none",text = element_text(size = 20)) +
#  xlab(label = "")+ylab("Gene Expression (LFC)")






gene_expression<-function(data,plot="boxplot",type="log2TPM"){
# Library required
  require(reshape2)
  require(ggplot2)
  require(tidyverse)



  tpm.melt<- melt(data)
  colnames(tpm.melt)<-c("gene","sample","tpm")
  tpm.melt$condition<-tpm.melt$sample
  tpm.melt$condition<-sub("CRHS96.*","96h",tpm.melt$condition)
  tpm.melt$condition<-sub("CRHS24.*","24h",tpm.melt$condition)
  tpm.melt$condition<-sub("CRHS120.*","120h",tpm.melt$condition)
  tpm.melt$condition<-sub("CRHS48.*","48h",tpm.melt$condition)
  tpm.melt$condition<-sub("CRHS144.*","144h",tpm.melt$condition)
  tpm.melt$condition<-sub("CRCK.*","0h",tpm.melt$condition)

  allGeneList  <-tpm.melt %>% group_by(condition,gene) %>% summarise(mean=mean(tpm),sd=sd(tpm),CV=sd(tpm)/mean(tpm))

 allGeneList$operon<-ifelse(allGeneList$gene %in% operonID,"operon","nonoperon")
 allGeneList$HH<-ifelse(allGeneList$gene %in% ci.HH.genes,"H2H","")


 allGeneList %>% filter(operon=="operon" ) %>% select(1:6)->aa1
 #allGeneList %>% filter(operon=="nonoperon") %>% select(1:6)->aa2
 allGeneList %>% filter(HH=="H2H") %>% select(1:5,HH)->aa3
 allGeneList %>% filter(operon!="operon" & HH!="H2H" ) %>% select(1:5,CV)->aa4
 aa4$type<-"Non-H2H non-operon"


 colnames(aa2)[6]<-"type"
 colnames(aa1)[6]<-"type"
 colnames(aa4)[6]<-"type"
 colnames(aa3)[6]<-"type"


 rbind(aa1,aa3,aa4)->aa
 aa$condition<-factor(aa$condition,levels = c("0h","24h","48h","96h","120h","144h"))
 aa$type<-factor(aa$type,levels = c("operon","H2H","Non-H2H non-operon"))

 if(plot=="boxplot" & type=="log2TPM"){
  p<- ggplot(aa,aes(x=condition,y=log2(mean+1),fill=type))+geom_boxplot(notch = T)+
     theme_bw()+guides(colour=guide_legend(title = ""))+theme(panel.grid.major = element_blank(),panel.background = element_blank(),panel.grid.minor = element_blank(),legend.position = "none",text = element_text(size = 20)) +
     xlab(label = "")+ylab("Gene Expression (log2 TPM)")
  return(p)

 } else if (plot=="boxplot" & type=='cvTPM'){
   p<-ggplot(aa,aes(x=condition,y=CV,fill=type))+geom_boxplot(notch = T)+
     theme_bw()+guides(colour=guide_legend(title = ""))+theme(panel.grid.major = element_blank(),panel.background = element_blank(),panel.grid.minor = element_blank(),legend.position = "none",text = element_text(size = 20)) +
     xlab(label = "")+ylab("Expression Variability (CV of TPM)")
   return(p)
 } else if(plot=="line" & type=="log2meanTPM"){
  p<- ggplot(aa,aes(x=condition,y=log2(mean+1),group=type,color=type))+stat_summary(fun.y = mean,geom = "line",position = position_dodge(width = 0.2)) +
     stat_summary(fun.data = mean_cl_boot, geom = "errorbar",position = position_dodge(width = 0.2),width=0.7)+
     theme_bw()+guides(colour=guide_legend(title = ""))+stat_summary(fun.y = mean,geom = "point",position = position_dodge(width = 0.2),size=3)+theme(panel.grid.major = element_blank(),panel.background = element_blank(),panel.grid.minor = element_blank(),text = element_text(size = 20),legend.position = "none") + xlab(label = "")+ylab("Mean Expression (log2 TPM)")
   return(p)
 } else if(plot=="line" & type=="meanTPM"){
  p<- ggplot(aa,aes(x=condition,y=mean,group=type,color=type))+stat_summary(fun.y = mean,geom = "line",position = position_dodge(width = 0.2)) +
     stat_summary(fun.data = mean_cl_boot, geom = "errorbar",position = position_dodge(width = 0.2),width=0.7)+
     theme_bw()+guides(colour=guide_legend(title = ""))+stat_summary(fun.y = mean,geom = "point",position = position_dodge(width = 0.2),size=3)+theme(panel.grid.major = element_blank(),panel.background = element_blank(),panel.grid.minor = element_blank(),text = element_text(size = 20),legend.position = "none") + xlab(label = "")+ylab("Mean Expression (TPM)")
   return(p)

 } else if (plot=="line" & type=="nor_meanTPM"){

   mean_normalizaed<-function(x){
     y<- ( x-mean(x))/(max(x)-min(x))

     return(y)
   }


   bb<-aa %>% group_by(type,gene) %>% summarise(log2nor_mean=mean_normalizaed(log2(mean+1)),nor_mean=mean_normalizaed(mean),condition=condition)



  p<- ggplot(bb,aes(x=condition,y=nor_mean,group=type,color=type))+stat_summary(fun.y = mean,geom = "line",position = position_dodge(width = 0.2)) +
    stat_summary(fun.data = mean_cl_boot, geom = "errorbar",position = position_dodge(width = 0.2),width=0.7)+
    theme_bw()+guides(colour=guide_legend(title = ""))+stat_summary(fun.y = mean,geom = "point",position = position_dodge(width = 0.2),size=3)+theme(panel.grid.major = element_blank(),panel.background = element_blank(),panel.grid.minor = element_blank(),text = element_text(size = 20),) +
    xlab(label = "")+ylab("Mean Normalization Gene Expression")
  return(p)

 } else if (plot=="line" & type=="log2nor_meanTPM"){
   bb<-aa %>% group_by(type,gene) %>% summarise(log2nor_mean=mean_normalizaed(log2(mean+1)),nor_mean=mean_normalizaed(mean),condition=condition)
   p<- ggplot(bb,aes(x=condition,y=log2nor_mean,group=type,color=type))+stat_summary(fun.y = mean,geom = "line",position = position_dodge(width = 0.2)) +
     stat_summary(fun.data = mean_cl_boot, geom = "errorbar",position = position_dodge(width = 0.2),width=0.7)+
     theme_bw()+guides(colour=guide_legend(title = ""))+stat_summary(fun.y = mean,geom = "point",position = position_dodge(width = 0.2),size=3)+theme(panel.grid.major = element_blank(),panel.background = element_blank(),panel.grid.minor = element_blank(),text = element_text(size = 20),) +
     xlab(label = "")+ylab("log2 Mean Normalization Gene Expression")
   return(p)
 }

}
















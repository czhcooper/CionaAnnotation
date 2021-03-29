
library(ggplot2)
library(tidyverse)

res1 %>% filter(HH =="HH" ) %>% select(log2FoldChange,condition)->aa1
aa1$Orientation<-"HH"
res1 %>% filter(HT =="HT") %>% select(log2FoldChange,condition)->aa2
aa2$Orientation<-"HT"
res1 %>% filter(TT =="TT") %>% select(log2FoldChange,condition)->aa3
aa3$Orientation<-"TT"
res1 %>% filter(Others =="Others") %>% select(log2FoldChange,condition)->aa4
aa4$Orientation<-"Others"
res1 %>% filter(operon =="operon") %>% select(log2FoldChange,condition)->aa5
aa5$Orientation<-"Operon"

rbind(aa1,aa2,aa3,aa4,aa5)->bb1

res2 %>% filter(HH =="HH") %>% select(log2FoldChange,condition)->aa1
aa1$Orientation<-"HH"
res2 %>% filter(HT =="HT") %>% select(log2FoldChange,condition)->aa2
aa2$Orientation<-"HT"
res2 %>% filter(TT =="TT") %>% select(log2FoldChange,condition)->aa3
aa3$Orientation<-"TT"
res2 %>% filter(Others =="Others") %>% select(log2FoldChange,condition)->aa4
aa4$Orientation<-"Others"
res2 %>% filter(operon =="operon") %>% select(log2FoldChange,condition)->aa5
aa5$Orientation<-"Operon"
rbind(aa1,aa2,aa3,aa4,aa5)->bb2

res3 %>% filter(HH =="HH") %>% select(log2FoldChange,condition)->aa1
aa1$Orientation<-"HH"
res3 %>% filter(HT =="HT") %>% select(log2FoldChange,condition)->aa2
aa2$Orientation<-"HT"
res3 %>% filter(TT =="TT") %>% select(log2FoldChange,condition)->aa3
aa3$Orientation<-"TT"
res3 %>% filter(Others =="Others") %>% select(log2FoldChange,condition)->aa4
aa4$Orientation<-"Others"
res3 %>% filter(operon =="operon") %>% select(log2FoldChange,condition)->aa5
aa5$Orientation<-"Operon"
rbind(aa1,aa2,aa3,aa4,aa5)->bb3

res4 %>% filter(HH =="HH") %>% select(log2FoldChange,condition)->aa1
aa1$Orientation<-"HH"
res4 %>% filter(HT =="HT") %>% select(log2FoldChange,condition)->aa2
aa2$Orientation<-"HT"
res4 %>% filter(TT =="TT") %>% select(log2FoldChange,condition)->aa3
aa3$Orientation<-"TT"
res4 %>% filter(Others =="Others") %>% select(log2FoldChange,condition)->aa4
aa4$Orientation<-"Others"
res4 %>% filter(operon =="operon") %>% select(log2FoldChange,condition)->aa5
aa5$Orientation<-"Operon"
rbind(aa1,aa2,aa3,aa4,aa5)->bb4

res5 %>% filter(HH =="HH") %>% select(log2FoldChange,condition)->aa1
aa1$Orientation<-"HH"
res5 %>% filter(HT =="HT") %>% select(log2FoldChange,condition)->aa2
aa2$Orientation<-"HT"
res5 %>% filter(TT =="TT") %>% select(log2FoldChange,condition)->aa3
aa3$Orientation<-"TT"
res5 %>% filter(Others =="Others") %>% select(log2FoldChange,condition)->aa4
aa4$Orientation<-"Others"
res5 %>% filter(operon =="operon") %>% select(log2FoldChange,condition)->aa5
aa5$Orientation<-"Operon"
rbind(aa1,aa2,aa3,aa4,aa5)->bb5


rbind(bb1,bb2,bb3,bb4,bb5)->bb

bb$condition<-factor(bb$condition,levels = c("24h","48h","96h","120h","144h"))
ggplot(bb,aes(x=abs(log2FoldChange),color=Orientation))+geom_density()+theme_bw()+guides(colour=guide_legend(title = ""))+theme(panel.grid.major = element_blank(),panel.background = element_blank(),panel.grid.minor = element_blank())

ggplot(bb,aes(x=condition,y=abs(log2FoldChange))) + geom_boxplot()

ggplot(bb,aes(x=condition,y=abs(log2FoldChange))) +stat_summary(fun.y = mean,geom = "line",position = position_dodge(width = 0.2)) +stat_summary(fun.data = mean_cl_boot, geom = "errorbar",position = position_dodge(width = 0.2),width=0.5)+xlab("condition")+ylab("Variability")+theme_bw()+guides(colour=guide_legend(title = ""))+theme(panel.grid.major = element_blank(),panel.background = element_blank(),panel.grid.minor = element_blank())

ggplot(bb[bb$Orientation!="Others",],aes(x=condition,y=abs(log2FoldChange),group=Orientation,color=Orientation))+stat_summary(fun.y = mean,geom = "line",position = position_dodge(width = 0.2)) +stat_summary(fun.data = mean_se, geom = "errorbar",position = position_dodge(width = 0.2),width=0.5)+xlab("condition")+ylab("Variability")+theme_bw()+guides(colour=guide_legend(title = ""))+theme(panel.grid.major = element_blank(),panel.background = element_blank(),panel.grid.minor = element_blank())+scale_color_manual(values = c("#d1495b" ,"#edae49","#00BFC4" ,"#66a182"))

ggplot(bb,aes(x=condition,y=log2FoldChange,group=Orientation,color=Orientation))+stat_summary(fun.y = mean,geom = "line",position = position_dodge(width = 0.2)) +stat_summary(fun.data = mean_cl_boot, geom = "errorbar",position = position_dodge(width = 0.2),width=0.5)+stat_summary(fun.y = mean,geom = "point",position = position_dodge(width = 0.2))+xlab("condition")+ylab("Mean LFC")+theme_bw()+guides(colour=guide_legend(title = ""))+theme(panel.grid.major = element_blank(),panel.background = element_blank(),panel.grid.minor = element_blank(),legend.position = "none")+scale_color_manual(values = c("#d1495b","#edae49","#00BFC4"  ,"#00A5FF" ,"#66a182"))


rbind(res1,res2,res3,res4,res5)->res
res$condition<-factor(res$condition,levels = c("24h","48h","96h","120h","144h"))
ggplot(res,aes(x=condition,y=log2FoldChange,group=operon,color=operon))+stat_summary(fun.y = mean,geom = "line",position = position_dodge(width = 0.2)) +stat_summary(fun.data = mean_se, geom = "errorbar",position = position_dodge(width = 0.2),width=0.5)+xlab("condition")+ylab("Mean LFC")+theme_bw()+guides(colour=guide_legend(title = ""))+theme(panel.grid.major = element_blank(),panel.background = element_blank(),panel.grid.minor = element_blank(),legend.position = "none")






rbind(res1,res2,res3,res4,res5)->res.all
sep1<-300
sep2<-600
sep3<-900

res.all %>% filter(HT =="HT" & operon!="operon" ) %>% dplyr::select(log2FoldChange,HTdistance,condition)->aa
aa$distance<-""
aa$distance[aa$HTdistance <=sep1]<-"a"
aa$distance[aa$HTdistance <=sep2 & aa$HTdistance >sep1]<-"b"
aa$distance[aa$HTdistance <=sep3 & aa$HTdistance >sep2 ]<-"c"
aa$distance[aa$HTdistance > sep3]<-"d"
aa->bb1



res.all %>% filter(HH =="HH" & operon!="operon" ) %>% dplyr::select(log2FoldChange,HTdistance,condition)->aa
aa$distance<-""
aa$distance[aa$HTdistance <=sep1]<-"a"
aa$distance[aa$HTdistance <=sep2 & aa$HTdistance >sep1]<-"b"
aa$distance[aa$HTdistance <=sep3 & aa$HTdistance >sep2 ]<-"c"
aa$distance[aa$HTdistance > sep3]<-"d"
aa->bb2

# combine all
colnames(bb1)[2]<-"Distance"
colnames(bb2)[2]<-"Distance"
bb1$Orientation<-"HT"
bb2$Orientation<-"HH"
rbind(bb1,bb2)->bb

na.omit(bb)->bb
bb$condition<-factor(bb$condition,levels = c("24h","48h","96h","120h","144h"))
ggplot(bb,aes(x=distance,y=abs(log2FoldChange),group=Orientation,color=Orientation))+stat_summary(fun.y = mean,geom = "line",position = position_dodge(width = 0.2)) +stat_summary(fun.data = mean_se, geom = "errorbar",position = position_dodge(width = 0.2),width=0.2, aes(shape=Orientation))+theme_bw()+guides(colour=guide_legend(title = ""))+theme(panel.grid.major = element_blank(),panel.background = element_blank(),panel.grid.minor = element_blank()) + facet_wrap(~condition,nrow = 1)

ggplot(bb,aes(x=distance,y=abs(log2FoldChange),group=Orientation,color=Orientation))+stat_summary(fun.y = mean,geom = "line",position = position_dodge(width = 0.2)) +stat_summary(fun.data = mean_se, geom = "errorbar",position = position_dodge(width = 0.2),width=0.2, aes(shape=Orientation))+theme_bw()+guides(colour=guide_legend(title = ""))+theme(panel.grid.major = element_blank(),panel.background = element_blank(),panel.grid.minor = element_blank()) +ylab("Variability")+scale_color_manual(values = c("#d1495b" ,"#edae49"))


pairwise.wilcox.test(abs(bb$log2FoldChange),g=bb$Orientation)


pairwise.wilcox.test(abs(bb$log2FoldChange[bb$Orientation=="HH"]),g=bb$distance[bb$Orientation=="HH"])

pairwise.wilcox.test(abs(bb$log2FoldChange[bb$Orientation=="HT"]),g=bb$distance[bb$Orientation=="HT"])

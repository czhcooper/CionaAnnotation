

ggplot(aaa,aes(x=condition,y=Esim,fill=type))+geom_boxplot(notch = T)+
  theme_bw()+guides(colour=guide_legend(title = ""))+theme(panel.grid.major = element_blank(),panel.background = element_blank(),panel.grid.minor = element_blank(),legend.position = "none",text = element_text(size = 20)) +
  xlab(label = "")+ylab("Correlation Coefficient")+scale_fill_manual(values = c("#F8766D","#39B600","#00B4EF","#BB9D00"))

thresfold<-0.5
all.Esim %>% filter(`0h`>=thresfold & type=="Operon" & `96h`>=thresfold & `24h`<=(-thresfold))->aa


thresfold<-0.5
all.Esim %>% filter(`0h`<=(-thresfold) & type=="Operon" & `96h`<=(-thresfold) )->aa


thresfold<-0.3
all.Esim %>% filter(`0h`<=(-thresfold) & type=="H2H" & `96h`<=(-thresfold) )->aa


all.Esim %>%filter(type=="Operon" )->aa.op
aa.op %>% filter(`0h`<(0.5) & `24h`>0.5 & `48h`>0.5)->up.op
 up.op.genes<-unlist(lapply(ci.operon_2gene,function(x){
  x[["GeneID"]][x[["Praent"]] %in% rownames(up.op)]
}),use.names = F)
up.op.GO<-GOenrich(up.op.genes,category = "BP")
GOplot(up.op.GO)


aa.op$CC<-ifelse(aa.op$`0h`>0.5,"HC",ifelse(aa.op$`0h`<(-0.5),"LC","NC"))
aa.op %>% gather(key="condition",value = "value",-type,-CC,-typeID)->aa
na.omit(aa)->aa
aa$condition<-factor(aa$condition,levels = c("0h","24h","48h","96h","120h","144h"))
ggplot(aa,aes(x=condition,y=value,fill=condition))+geom_boxplot()+
  theme_bw()+guides(colour=guide_legend(title = ""))+theme(panel.grid.major = element_blank(),panel.background = element_blank(),panel.grid.minor = element_blank(),legend.position = "none",text = element_text(size = 20)) +
  xlab(label = "")+ylab("Correlation Coefficient")+facet_grid(cols = vars(CC))


all.Esim %>% filter(type=="H2H")->aa.HH
aa.HH %>% filter(`0h`<0.5 & `24h`>0.5 &`48h`>0.5)->up.HH
aa.HH$CC<-ifelse(aa.HH$`0h`>0.5,"HC",ifelse(aa.HH$`0h`<(-0.5),"LC","NC"))
aa.HH %>% gather(key="condition",value = "value",-type,-CC,-typeID)->aa
na.omit(aa)->aa
aa$condition<-factor(aa$condition,levels = c("0h","24h","48h","96h","120h","144h"))
ggplot(aa,aes(x=condition,y=value,fill=condition))+geom_boxplot()+
  theme_bw()+guides(colour=guide_legend(title = ""))+theme(panel.grid.major = element_blank(),panel.background = element_blank(),panel.grid.minor = element_blank(),legend.position = "none",text = element_text(size = 20)) +
  xlab(label = "")+ylab("Correlation Coefficient")+facet_grid(cols = vars(CC))

    aa<- ci.HH.Ex[up.HH$typeID %in% names(ci.HH.Ex)]
 up.HH.genes<- unlist(lapply(aa,function(x){
    colnames(x)
  }),use.names = F)
up.HH.GO.BP<-GOenrich(up.HH.genes,category = "BP")
GOplot(up.HH.GO.BP)
up.HH.GO.MF<-GOenrich(up.HH.genes,category = "MF")
GOplot(up.HH.GO.MF)
up.HH.GO.CC<-GOenrich(up.HH.genes,category = "CC")
GOplot(up.HH.GO.CC)
up.HH.GO$GeneRatio<-up.HH.GO$Significant/up.HH.GO$Annotated
up.HH.GO %>%  arrange(classic)->bb
ggplot(bb[1:20,],aes(x=GeneRatio,y=Term))+geom_point(aes(size=Significant,color=classic))+facet_grid(rows = vars(type))

up.HH.GO.BP$GeneRatio<-up.HH.GO.BP$Significant/up.HH.GO.BP$Annotated
ggplot(up.HH.GO.BP[1:15,],aes(x=GeneRatio,y=Term))+geom_point(aes(size=Significant,color=classic))  +theme_bw(base_size = 14) +
  scale_colour_gradient(limits=c(0, 0.01), low="red") +labs(size="Count",color="p-value")+ylab("Biological Process")

up.op.GO.BP$GeneRatio<-up.op.GO.BP$Significant/up.op.GO.BP$Annotated
ggplot(up.op.GO.BP[1:15,],aes(x=GeneRatio,y=Term))+geom_point(aes(size=Significant,color=classic))  +theme_bw(base_size = 14) +
  scale_colour_gradient(limits=c(0, 0.05), low="red") +labs(size="Count",color="p-value")+ylab("Biological Process")



summary(ci.random.GO.sim)
 aa1<-data.frame(type="Random",GOsim=ci.random.GO.sim)
 aa2<-data.frame(type="Operon",GOsim=ci.operon.GO.sim)
 aa3<-data.frame(type="H2H",GOsim=ci.HH.GO.sim)
 aa4<-data.frame(type="H2T",GOsim=ci.HT.GO.sim)
 aa5<-data.frame(type="T2T",GOsim=ci.TT.GO.sim)

 rbind(aa1,aa2,aa3,aa4,aa5)->all.GOsim
 all.GOsim$type<-factor(all.GOsim$type,levels = c("Operon","H2H","H2T","T2T","Random"))
ggplot(all.GOsim,aes(x=type,y=GOsim,fill=type))+geom_boxplot(notch = T)  +ylab("GO semantic similarity")+
  theme(legend.position = "none",text = element_text(size = 20))+xlab("")




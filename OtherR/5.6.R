library(CionaAnnotation)
library(ggplot2)
library(topGO)
library(dplyr)
library(reshape2)
library(tidyverse)
library(GOSim)
library(ggdendro)
library("grid")

#functional analysis of operon and H2H
ci.operon.BP<-GOenrich(as.data.frame(operonID),category = "BP")
ci.HH.genes<-unique(unlist(lapply(ci.HH,function(x){
  x[,"Geneid"]
}),use.names = F))
ci.H2H.BP<-GOenrich(as.data.frame(ci.HH.genes),category = "BP")

aa1<-ci.operon.BP[c(3,5,10,11,12,13,14,15,16,17,18,19,20,21,23,24,25,27,28,29,30,32,33,34,57,58,49,53),]
aa1$type<-"Operon"

aa2<-ci.H2H.BP[c(5,6,8,9,10,11,12,13,14,15,17,18,19,20,21,22,23,24,27,28,29,30,33,34,36,37,38,40,41,47,49,56,57,58),]
aa2$type<-"H2H"

operon.H2H.BP<-rbind(aa1,aa2)

operon.H2H.BP.GOsim<-getTermSim(unique(operon.H2H.BP$GO.ID),method = "Resnik")

operon.H2H.BP.GOsim.scaled<-1- operon.H2H.BP.GOsim
operon.H2H.BP.GOsim.scaled<-as.dist(operon.H2H.BP.GOsim.scaled)

model<-hclust(operon.H2H.BP.GOsim.scaled,"ave")

p1<-ggdendrogram(hc,rotate = T,segments = T,labels = T)
dhc <- as.dendrogram(model)
# Rectangular lines
ddata <- dendro_data(dhc, type = "rectangle")
p1 <- ggplot(segment(ddata)) +
  geom_segment(aes(x = x, y = y, xend = xend, yend = yend)) +
  coord_flip() +
  scale_y_reverse(expand = c(0.2, 0))



#heatmap
operon.H2H.BP$type<-factor(operon.H2H.BP$type,levels = c("Operon","H2H"))
p2<-ggplot(operon.H2H.BP,aes(x=type,y=Term))+geom_point(aes(color= -log10(classic),size=Significant))+xlab("")+ylab("")+
  theme_classic()+theme()+
  scale_color_gradient(low = "blue",high = "red")

grid.newpage()
print(p2, vp = viewport(x = 0.4, y = 0.5, width = 0.8, height = 1.0))
print(p1, vp = viewport(x = 0.90, y = 0.445, width = 0.2, height = 1.0))

#re-order
otter.order<-order.dendrogram(dhc)
operon.H2H.BP$GO.ID<-factor(operon.H2H.BP$GO.ID,levels = c( unique(operon.H2H.BP$GO.ID)[otter.order]))
operon.H2H.BP$Term<-factor(operon.H2H.BP$Term,levels = unique(operon.H2H.BP$Term[operon.H2H.BP$GO.ID] ) )

p2<-ggplot(operon.H2H.BP,aes(x=type,y=Term))+geom_point(aes(color= -log10(classic),size=Significant))+xlab("")+ylab("")+
  theme_classic()+theme()+
  scale_color_gradient(low = "blue",high = "red")




#expression

test.tpm<-KH.tpm %>% select(CRCK.1:CRCK.30,CRHS24.2:CRHS48.6)
test.tpm$GendID<-rownames(test.tpm)
test.tpm->data

operon.H2H.tpm
ggplot(operon.H2H.tpm,aes(x=condition,y=log2(mean+1),fill=type))+geom_boxplot(notch = T)+
  theme_bw()+guides(colour=guide_legend(title = ""))+theme(panel.grid.major = element_blank(),panel.background = element_blank(),panel.grid.minor = element_blank(),legend.position = "none",text = element_text(size = 20)) +
  xlab(label = "")+ylab("Gene Expression (log2 TPM+1)")

ggplot(operon.H2H.tpm,aes(x=condition,y=CV,fill=type))+geom_boxplot(notch = T)+
  theme_bw()+guides(colour=guide_legend(title = ""))+theme(panel.grid.major = element_blank(),panel.background = element_blank(),panel.grid.minor = element_blank(),legend.position = "none",text = element_text(size = 20)) +
  xlab(label = "")+ylab("Gene Expression Variability (CV of TPM)")

operon2.gene<-unlist(lapply(ci.operon_2gene,function(x){
  x[,"GeneID"]
}),use.names = F)

operon2gene.tpm<-test.tpm %>% filter(GendID %in% operon2.gene)

split(operon2gene.tpm,operon2gene.tpm$GendID)


operon2.gene.cor<- lapply(ci.operon_2gene,function(x){
  y<-x %>% select(CRCK.1:CRCK.30,CRHS24.2:CRHS48.6)
  y<-apply(y,1,function(x){
    (x-mean(x))/sd(x)
  })
   h0<-cor(y[1:5,])
   a1<-data.frame(cor=h0[1,2],condition="0h")

   h24<-cor(y[6:9,])
   a2<-data.frame(cor=h24[1,2],condition="24h")


   h48<-cor(y[10:15,])
   a3<-data.frame(cor=h48[1,2],condition="48h")

   rbind(a1,a2,a3)->a
   return(a)
   })

operon2.gene.cor<-bind_rows(operon2.gene.cor,.id = "genePairsID")



H2H.gene.cor<-lapply(ci.HH,function(x){
  y<-x %>% select(CRCK01.bam:CRCK30.bam,CRHS24.2.bam:CRHS24.6.bam,CRHS48.1.bam:CRHS48.6.bam)
  y<-apply(y,1,function(x){
    (x-mean(x))/sd(x)
  })

  h0<-cor(y[1:5,])
  a1<-data.frame(cor=h0[1,2],condition="0h")

  h24<-cor(y[6:9,])
  a2<-data.frame(cor=h24[1,2],condition="24h")


  h48<-cor(y[10:15,])
  a3<-data.frame(cor=h48[1,2],condition="48h")

  rbind(a1,a2,a3)->a
  return(a)
})

H2H.gene.cor<-bind_rows(H2H.gene.cor,.id = "genePairsID")

operon.H2H.gene.cor<-rbind(operon2.gene.cor,H2H.gene.cor)


HC<-operon2.gene.cor$genePairsID[operon2.gene.cor$cor>0.5 & operon2.gene.cor$condition=="0h" ]
LC<-operon2.gene.cor$genePairsID[operon2.gene.cor$cor<(-0.5) & operon2.gene.cor$condition=="0h" ]
NC<-operon2.gene.cor$genePairsID[operon2.gene.cor$cor>=(-0.5) &  operon2.gene.cor$cor<=0.5 & operon2.gene.cor$condition=="0h" ]

operon2.gene.cor$type<-ifelse(operon2.gene.cor$genePairsID %in% HC,"HC",ifelse( operon2.gene.cor$genePairsID %in% LC,"LC", "NC"))


ggplot(na.omit(operon2.gene.cor),aes(x=condition,y=cor))+ geom_boxplot(notch = T)+facet_grid(cols = vars(type))+theme_classic()+
  theme(text = element_text(size = 20))+ylab("Expression Correlation")+xlab("")

HC2LC.operon<-operon2.gene.cor$genePairsID[operon2.gene.cor$type=="HC" &  operon2.gene.cor$cor[operon2.gene.cor$condition=="24h"] <(-0.5) & operon2.gene.cor$cor[operon2.gene.cor$condition=="48h"] <(-0.5) ]
HC2LC.operon.genes<- unlist(lapply(ci.operon[names(ci.operon) %in% HC2LC.operon],function(x){
   x[,"GeneID"]
 }),use.names = F)


a<-GOenrich(as.data.frame(HC2LC.operon.genes),category = "BP")

GOplot(a)



HC<-H2H.gene.cor$genePairsID[H2H.gene.cor$cor>0.5 & H2H.gene.cor$condition=="0h" ]
LC<-H2H.gene.cor$genePairsID[H2H.gene.cor$cor<(-0.5) & H2H.gene.cor$condition=="0h" ]
NC<-H2H.gene.cor$genePairsID[H2H.gene.cor$cor>=(-0.5) & H2H.gene.cor$cor<=0.5 & H2H.gene.cor$condition=="0h" ]

H2H.gene.cor$type<-ifelse(H2H.gene.cor$genePairsID %in% HC,"HC",ifelse(H2H.gene.cor$genePairsID %in% LC,"LC", "NC"))


ggplot(na.omit(H2H.gene.cor),aes(x=condition,y=cor))+ geom_boxplot(notch = T)+facet_grid(cols = vars(type))+theme_classic()+
  theme(text = element_text(size = 20))+ylab("Expression Correlation")+xlab("")







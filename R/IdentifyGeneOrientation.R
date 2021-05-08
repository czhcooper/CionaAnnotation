##' calculate the expression correlation between  gene pairs
##'
##'
##' @title Identify gene organization
##' @param
##' @param
##' @return A list that contains gene orientation and genes
##' @author Chen Zaohuang
##' @export


# Input should be a data.frame containing colnumn names "Geneid","Chr","Start", "End", "Strand" in order, others columns can exist but are not necessary.

geneOrientation<-function(data,Orientation,maxDistance=1000,minDistance=1,size=2){

  #Libaray required
  require(tidyverse)
  require(plyr)
  require(dplyr)


  # make a slide window function
  slidingwindow <- function(df, size){
    df <- data.frame(df)
    rownames(df) <- seq_len(nrow(df))
    windows <- alply(1:(nrow(df)-size+1), 1, function(x){c(x:(x+size-1))})

    # For each window, apply the function to the subset of indices
    outmat <- lapply(windows,
                     function(indices){
                       subset(df, as.numeric(rownames(df)) %in% indices)})

  # Sort the output into order of data frame rows
  }


  #identify head-to-hail orientation
  if(Orientation=="HH"){
   aa<- data %>% arrange(Chr,Start)

   slidingwindow(aa,size)->cc
   lapply(cc, function(x)  length(unique(x["Chr"]) ==1 ))->dd # make sure they are in a same chromosome
   unlist(dd)->ee
   ee[ee==FALSE]->ff    #remove these, because you may get a pair of genes from different chromosomes or scaffolds
   cc[names(ee[ee==TRUE])]->cc   #Now you have corrected data

   unlist(lapply(cc, function(x) x[5][1,]=="-" & x[5][2,]=="+" & x[3][2,]-x[4][1,]>=minDistance & x[3][2,]-x[4][1,] <=maxDistance))->HH.ee
   cc[names(HH.ee[HH.ee==TRUE])]->HH.cc



   return(HH.cc)


   # identify head-to-tail orientation
  } else if(Orientation=="HT"){

    aa<- data %>% arrange(Chr,Strand,Start)
    slidingwindow(aa,size)->cc
    lapply(cc, function(x)  length(unique(x[,"Chr"])  )==1)->HT.dd # make sure they are in a same chromosome
    unlist(HT.dd)->HT.ee
    cc[names(HT.ee[HT.ee==TRUE])]->HT.cc

    cc.dis<-lapply(cc, function(x){         # calculate intergenic distance within each group
      x$Dis<-NA
      for (i in 1:(dim(x)[1]-1)){
        x[i+1,"Dis"]<-x[i+1,"Start"]-x[i,"End"]
      }
      y<-na.omit(x)
      yx <-all(y[,"Dis"] <maxDistance & y[,"Dis"] > minDistance)
      return(yx)
    })


     HT.ee<-unlist(cc.dis)
     HT.cc[names(HT.ee[HT.ee==TRUE])]->HT.cc
      #names(HT.cc)<-paste0(Orientation,".",size,".",seq(1:length(names(HT.cc))))


    return(HT.cc)

    # identify tail-to-tail orientation
  } else if (Orientation=="TT"){
    aa<- data %>% arrange(Chr,Start)

    slidingwindow(aa,size)->cc
    lapply(cc, function(x)  length(unique(x["Chr"]) ==1 ))->dd # make sure they are in a same chromosome
    unlist(dd)->ee
    ee[ee==FALSE]->ff    #remove these
    cc[names(ee[ee==TRUE])]->cc   #Now you have true data
    ##head to head , without overlapping    1kb from where?
    unlist(lapply(cc, function(x) x[5][1,]=="+" & x[5][2,]=="-" & x[3][2,]-x[4][1,]>=minDistance & x[3][2,]-x[4][1,] <=maxDistance))->TT.ee
    cc[names(TT.ee[TT.ee==TRUE])]->TT.cc


    return(TT.cc)


  }


}



# identify Operon

#ee <-unlist(lapply(ci.HT.operon4,function(x){
#  all(x[[1]] %in% unique(unlist(lapply(ci.HT.operon5,function(x)x[[1]]),use.names = F)))

#}),use.names = F)

#ci.HT.operon4[!ee]->ci.HT.operon4

#ee <-unlist(lapply(ci.HT.operon3,function(x){
#  all(x[[1]] %in% c(unique(unlist(lapply(ci.HT.operon5,function(x)x[[1]]),use.names = F)),
#      unique(unlist(lapply(ci.HT.operon4,function(x)x[[1]]),use.names = F))))

#}),use.names = F)

#ci.HT.operon3[!ee]->ci.HT.operon3

#ee <-unlist(lapply(ci.HT.operon2,function(x){
#  all(x[[1]] %in% c(unique(unlist(lapply(ci.HT.operon5,function(x)x[[1]]),use.names = F)),
#                    unique(unlist(lapply(ci.HT.operon4,function(x)x[[1]]),use.names = F)),
#                    unique(unlist(lapply(ci.HT.operon3,function(x)x[[1]]),use.names = F)) ))

#}),use.names = F)

#ci.HT.operon2[!ee]->ci.HT.operon2







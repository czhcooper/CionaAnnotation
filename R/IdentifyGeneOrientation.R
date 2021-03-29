

# Input should be a data.frame containing colnumn names "Geneid","Chr","Start", "End", "Strand" in order, others columns can exist but are not necessary.

geneOrientation<-function(data,Orientation,maxDistance=1000,minDistance=1){
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

   slidingwindow(aa,2)->cc
   lapply(cc, function(x) identical(x[2][1,],x[2][2,]))->dd # make sure they are in a same chromosome
   unlist(dd)->ee
   ee[ee==FALSE]->ff    #remove these, because you may get a pair of genes from different chromosomes or scaffolds
   cc[names(ee[ee==TRUE])]->cc   #Now you have corrected data

   unlist(lapply(cc, function(x) x[5][1,]=="-" & x[5][2,]=="+" & x[3][2,]-x[4][1,]>=minDistance & x[3][2,]-x[4][1,] <=maxDistance))->HH.ee
   cc[names(HH.ee[HH.ee==TRUE])]->HH.cc



   return(HH.cc)


   # identify head-to-tail orientation
  } else if(Orientation=="HT"){

    aa<- data %>% arrange(Chr,Strand,Start)
    slidingwindow(aa,2)->cc
    lapply(cc, function(x) identical(x[2][1,],x[2][2,]))->HT.dd # make sure they are in a same chromosome
    unlist(HT.dd)->HT.ee
    cc[names(HT.ee[HT.ee==TRUE])]->HT.cc

    unlist(lapply(cc, function(x)  identical(x[2][1,],x[2][2,])& x[3][2,]-x[4][1,]>= minDistance & x[3][2,]-x[4][1,] <=maxDistance))->HT.ee
    HT.cc[names(HT.ee[HT.ee==TRUE])]->HT.cc


    return(HT.cc)

    # identify tail-to-tail orientation
  } else if (Orientation=="TT"){
    aa<- data %>% arrange(Chr,Start)

    slidingwindow(aa,2)->cc
    lapply(cc, function(x) identical(x[2][1,],x[2][2,]))->dd # make sure they are in a same chromosome
    unlist(dd)->ee
    ee[ee==FALSE]->ff    #remove these
    cc[names(ee[ee==TRUE])]->cc   #Now you have true data
    ##head to head , without overlapping    1kb from where?
    unlist(lapply(cc, function(x) x[5][1,]=="+" & x[5][2,]=="-" & x[3][2,]-x[4][1,]>=minDistance & x[3][2,]-x[4][1,] <=maxDistance))->TT.ee
    cc[names(TT.ee[TT.ee==TRUE])]->TT.cc


    return(TT.cc)


  }


}



rm(list = ls())

library(picante)
library(tidyverse)
library(magrittr)
library(ape)
library(foreach)
library(doParallel)
library(abind)
library(data.table)


### work in the tmp directory
df = read.csv('comm.tsv',sep='\t',header=T, row.names = 1)
saveRDS(df, "comm.rds")
mf = read.csv('meta_with_clinical.tsv',sep='\t', header=T, row.names =1)
saveRDS(mf,'traits.rds')

phy = read.tree("asv_iq_tree.treefile")
comm = readRDS("comm.rds")
traits = readRDS("traits.rds")

phylo<-phy
comun=t(comm)
dim(comun)


#combine tree and OTU table
phylocom = match.phylo.data(phylo, comun)
str(phylocom)
rm(phy, traits, comm, comun, phylo)

#calculate bMNTD
bMNTD = as.matrix(comdistnt(t(phylocom$data),cophenetic(phylocom$phy),abundance.weighted=T))
dim(bMNTD)

saveRDS(bMNTD,'wbMNTD.rds')

identical(colnames(phylocom$data),colnames(bMNTD))
identical(colnames(phylocom$data),rownames(bMNTD))

#definining functions

funcomdist = function(x) {
    as.matrix(comdistnt(t(x$data),taxaShuffle(cophenetic(x$phy)),abundance.weighted=T,exclude.conspecifics = F))
}

acomb <- function(...) abind(..., along=3)


#setup parallel

cores=detectCores()
cl <- makeCluster(cores[1]-1) #not to overload the computer
cl <- makeCluster(2)
registerDoParallel(cl)

null_bMNTD = foreach(i=1:999, .packages = c("picante","abind"), .combine = "acomb", .verbose = TRUE) %dopar% {
    funcomdist(phylocom)
}

stopCluster(cl)

dim(null_bMNTD)

#save the null distribution
saveRDS(null_bMNTD, "null_wbMTND.rds")

#put it together

wbNTI = matrix(c(NA),nrow=ncol(phylocom$data),ncol=ncol(phylocom$data))
dim(wbNTI)

for(columns in 1:(ncol(phylocom$data)-1)) {
    for(rows in (columns+1):ncol(phylocom$data)) {
        rand.vals = null_bMNTD[rows,columns,];
        wbNTI[rows,columns] = (bMNTD[rows,columns] - mean(rand.vals)) / sd(rand.vals)
        rm("rand.vals")
    }
}

rownames(wbNTI) = colnames(phylocom$data);
colnames(wbNTI) = colnames(phylocom$data);

saveRDS(wbNTI, "bNTI.rds")




##### raup crick

rm(list = ls())
library(picante)
library(magrittr)
library(data.table)
library(ecodist)
library(doParallel)
library(foreach)
library(NST)

setwd('/media/christos/ssd/work/Infants/publication/tmp')
df = read.csv('comm.tsv',sep='\t',header = T, row.names=1)
saveRDS(df, "comm.rds")
comm = readRDS("comm.rds") %>% as.matrix()




raup_crick_abu_par <- function(com, reps, ncore, classic_metric=FALSE, split_ties=TRUE){
  
  
  require("parallel")
  
  require("doSNOW")
  
  require("vegan")
  
  
  pb <- txtProgressBar(max =reps, style = 3)
  
  progress <- function(n) setTxtProgressBar(pb, n)
  
  opts <- list(progress = progress)
  
  cl <- makeCluster(ncore)
  
  registerDoSNOW(cl)
  
  bray.rand <- foreach(randomize = 1:reps,
                       
                       .options.snow = opts,
                       
                       .packages = c("vegan", "picante")) %dopar% {
                         
                         null.dist <- com*0
                         print(null.dist)
                         
                         for(i in 1:nrow(com)){

                           com.pa <- (com>0)*1
                           
                           gamma<-ncol(com)
                           
                           occur<-apply(com>0, MARGIN=2, FUN=sum)
                           
                           abundance<-apply(com, MARGIN=2, FUN=sum)
                           
                           com1 <- rep(0,gamma)

                           com1[sample(1:gamma, sum(com.pa[i,]), replace=FALSE, prob=occur)]<-1
                           
                           com1.samp.sp = sample(which(com1>0), (sum(com[i,])-sum(com1)),
                                                 
                                                 replace=TRUE,prob=abundance[which(com1>0)]);
                           
                           com1.samp.sp = cbind(com1.samp.sp,1)
                           
                           com1.sp.counts = as.data.frame(tapply(com1.samp.sp[,2],com1.samp.sp[,1],FUN=sum))
                           
                           colnames(com1.sp.counts) = 'counts'
                           
                           com1.sp.counts$sp = as.numeric(rownames(com1.sp.counts))
                           
                           com1[com1.sp.counts$sp] = com1[com1.sp.counts$sp] + com1.sp.counts$counts
                           
                           x <- com1
                           
                           null.dist[i,] <- x
                           
                           rm('com1.samp.sp','com1.sp.counts')
                           
                         }
                         
                         as.matrix(vegdist(null.dist, "bray"))
                         
                       }
  
  stopCluster(cl)
  
  
  
  ## Calculate beta-diversity for obs metacommunity
  
  bray.obs <- as.matrix(vegdist(com, "bray"))
  
  
  
  ##how many null observations is the observed value tied with?
  
  null_bray_curtis <- bray.rand
  
  num_exact_matching_in_null <- lapply(null_bray_curtis, function(x) x==bray.obs)
  
  num_exact_matching_in_null <- apply(simplify2array(num_exact_matching_in_null), 1:2, sum)
  
  
  
  ##how many null values are smaller than the observed dissimilarity?
  
  num_less_than_in_null <- lapply(null_bray_curtis, function(x) (x<bray.obs)*1)
  
  num_less_than_in_null <- apply(simplify2array(num_less_than_in_null), 1:2, sum)
  
  
  
  
  
  rc = (num_less_than_in_null)/reps; # rc;
  
  
  
  if(split_ties){
    
    
    
    rc = ((num_less_than_in_null +(num_exact_matching_in_null)/2)/reps)
    
  };
  
  
  
  
  
  if(!classic_metric){
    
    
    
    ##our modification of raup crick standardizes the metric to range from -1 to 1 instead of 0 to 1
    
    
    
    rc = (rc-.5)*2
    
  };
  
  
  
  return(rc)
  
  
  
}


RC.bray=raup_crick_abu_par(comm,reps=999,ncore = 4, classic_metric = FALSE)
x = rownames(comm)
x = cbind(x,"A") %>% as.matrix()
rownames(x) = x[,1]
x = x[,2]
x = x %>% as.matrix()

Sys.setenv(OMP_NUM_THREADS = 1)
Sys.setenv(MKL_NUM_THREADS = 1)

NST = tNST(comm, group = x, dist.method = "bray", abundance.weighted = TRUE, rand = 1000, nworker = 4, SES = TRUE, RC = TRUE, between.group = FALSE)
NST@index.pair
NST$index.pair$RC.bray %>% as.matrix
NST
#RC.mat=as.matrix(RC.bray)
saveRDS(NST, "RC_bray.rds")


#RC.bray=raup_crick_abu_par(comm,reps=999,ncore = 4, classic_metric = FALSE)
#RC.mat=as.matrix(RC.bray)
#saveRDS(RC.mat, "RC_bray.rds")

  

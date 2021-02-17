library(data.table)
library(vegan)
library(MASS)
data(varespec, varechem)

### work in the tmp directory

### RUN clinical_data.py before!!! 

M = read.table('meta_envfit.tsv', header=T, row.names=1, sep='\t')

C = read.table('CountTable.tsv',header=T, row.names = 1, sep='\t')
dim(M)

D = vegdist(C)
mds.stuff <- cmdscale(D, eig=TRUE,x.ret = TRUE)
mds.var.per <- round(mds.stuff$eig/sum(mds.stuff$eig)*100, 1)
mds.values <- mds.stuff$points
mds.data <- data.table(SampleID=rownames(mds.values),
                       X=mds.values[,1],
                       Y=mds.values[,2])
mds.values

expVar = c('Age','SEX', 'DM','GA','AntiCumIndex','oxygenCumIndex','il6WinIndex', 'Enteral_feeding')
preVar = c('ROP','IVH','BPD','Late.onset.sepsis')


## All gut samples 
Gmexp = subset(M, M$SampleType == 'Gut')[expVar]
Gmpre = subset(M, M$SampleType == 'Gut')[preVar]
Gc = C[rownames(Gmexp),]
Gc = Gc[,which(colSums(as.matrix(Gc)) > 0)]
envfit(ord=cmdscale(vegdist(as.matrix(Gc, rownames='SampleID'), 'bray'), eig=TRUE, x.ret = TRUE)$points, 
                 env=Gmexp, perm = 10000)
envfit(ord=cmdscale(vegdist(as.matrix(Gc, rownames='SampleID'), 'bray'), eig=TRUE, x.ret = TRUE)$points, 
                 env=Gmpre, perm = 10000)



## All oral samples 
Omexp = subset(M, M$SampleType == 'Oral_cavity')[expVar]
Ompre = subset(M, M$SampleType == 'Oral_cavity')[preVar]
Oc = C[rownames(Omexp),]
Oc = Oc[,which(colSums(as.matrix(Oc)) > 0)]
envfit(ord=cmdscale(vegdist(as.matrix(Oc, rownames='SampleID'), 'bray'), eig=TRUE, x.ret = TRUE)$points, 
       env=Omexp, perm = 10000)
envfit(ord=cmdscale(vegdist(as.matrix(Oc, rownames='SampleID'), 'bray'), eig=TRUE, x.ret = TRUE)$points, 
       env=Ompre, perm = 10000)


## All skin samples 
Smexp = subset(M, M$SampleType == 'Skin')[expVar]
Smpre = subset(M, M$SampleType == 'Skin')[preVar]
Sc = C[rownames(Smexp),]
Sc = Sc[,which(colSums(as.matrix(Sc)) > 0)]
envfit(ord=cmdscale(vegdist(as.matrix(Sc, rownames='SampleID'), 'bray'), eig=TRUE, x.ret = TRUE)$points, 
       env=Smexp, perm = 10000)
envfit(ord=cmdscale(vegdist(as.matrix(Sc, rownames='SampleID'), 'bray'), eig=TRUE, x.ret = TRUE)$points, 
       env=Smpre, perm = 10000)











  

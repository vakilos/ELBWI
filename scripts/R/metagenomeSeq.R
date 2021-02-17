library(metagenomeSeq)


#### work in metagenomeseq directory
## Load data, count data is already filtered for depth and samples are removed!
count = loadMeta(file.path("count_data.tsv"))
tax = read.delim(file.path('tax_data.tsv'), stringsAsFactors = FALSE)
meta = loadPhenoData(file.path('meta_data.tsv'), tran=TRUE, sep='\t')

### order data
ord = match(colnames(count$counts), rownames(meta))
meta = meta[ord, ]
meta
dim(count$counts)
### create MRexperiment object
phenotypeData = AnnotatedDataFrame(meta)
OTUdata = AnnotatedDataFrame(tax)

## create obj and normalize
obj = newMRexperiment(count$counts, phenoData = phenotypeData, featureData = OTUdata)
p = cumNormStatFast(obj)
data = cumNorm(obj, p = p)


taxa = sapply(fData(data)$ASVg,
              function(i){
                i[length(i)]
              })


### Subset data for twins
twindata = data[, which(pData(obj)$TWIN %in% c('1.0','4.0'))] ## include only 2 groups
#levels(pData(twindata)$)
pd <- pData(twindata)
mod <- model.matrix(~1 + TWIN, data=pd)
dataRes1 = fitFeatureModel(twindata, mod)
result <- MRscoefs(dataRes1, taxa=taxa, number = 50)
result[which(result$adjPvalues < 0.05),]


## Subset data for dbClusters
dbClusterData = data[, which(pData(obj)$dbCluster %in% c(2,3))]
levels(pData(dbClusterData)$dbCluster)
pd <- pData(dbClusterData)
mod <- model.matrix(~1 + dbCluster, data=pd)
dataRes1 = fitFeatureModel(dbClusterData, mod)
result <- MRcoefs(dataRes1, taxa=taxa, number=1000)
result[which(result$adjPvalues < 0.05),]

## subset data for Clusters
ClusterData = data[, which(pData(obj)$Cluster %in% c(3,2))]
levels(pData(ClusterData)$Cluster)
pd <- pData(ClusterData)
mod <- model.matrix(~1 + Cluster, data=pd)
dataRes1 = fitFeatureModel(ClusterData, mod)
result <- MRcoefs(dataRes1, taxa=taxa, number=1000)
result[which(result$adjPvalues < 0.05),]


d## Subset data for Sampletypes
sampleTypedata = data[, which(pData(obj)$SampleType %in% c('Skin','Gut'))]
levels(pData(sampleTypedata)$SampleType)
pd <- pData(sampleTypedata)
mod <- model.matrix(~1 + SampleType, data=pd)
dataRes1 = fitFeatureModel(sampleTypedata, mod)
result <- MRcoefs(dataRes1, taxa=taxa, number = 1000)
result[which(result$adjPvalues < 0.05),]

## Subset data for DM

pd <- pData(data)
mod <- model.matrix(~1 + DM, data=pd)
levels(pData(data)$DM)
dataRes1 = fitFeatureModel(data, mod)
result <- MRcoefs(dataRes1, taxa=taxa, number = 1000)
result[which(result$adjPvalues < 0.05),]

## subset data for BPD
BPDdata = data[, which(pData(obj)$BPD %in% c('2','0'))]
levels(pData(BPDdata)$BPD)
pd <- pData(BPDdata)
mod <- model.matrix(~1 + BPD, data=pd)
dataRes1 = fitFeatureModel(BPDdata, mod)
result <- MRcoefs(dataRes1, taxa=taxa, number = 1000)
result[which(result$adjPvalues < 0.05),]


## diff abundance for SEX
levels(pData(data)$SEX)
pd <- pData(data)
mod <- model.matrix(~1 + SEX, data=pd)
dataRes1 = fitFeatureModel(data, mod)
result <- MRcoefs(dataRes1, taxa=taxa, number = 1000)
result[which(result$adjPvalues < 0.05),]


## diff abundance for EOS (early onset sepsis)
levels(pData(data)$EOS)
pd <- pData(data)
mod <- model.matrix(~1 + EOS, data=pd)
dataRes1 = fitFeatureModel(data, mod)
result <- MRcoefs(dataRes1, taxa=taxa, number = 1000)
result[which(result$adjPvalues < 0.05),]

## subset data for ROP
ROPdata = data[, which(pData(obj)$ROP %in% c('0','1'))]
levels(pData(ROPdata)$ROP)
pd <- pData(ROPdata)
mod <- model.matrix(~1 + ROP, data=pd)
dataRes1 = fitFeatureModel(ROPdata, mod)
result <- MRcoefs(dataRes1, taxa=taxa, number = 50)
result[which(result$adjPvalues < 0.05),]


## subset data for IVH
IVHdata = data[, which(pData(obj)$IVH %in% c('1','2'))]
levels(pData(IVHdata)$IVH)
pd <- pData(IVHdata)
mod <- model.matrix(~1 + IVH, data=pd)
dataRes1 = fitFeatureModel(IVHdata, mod)
result <- MRcoefs(dataRes1, taxa=taxa, number = 1000)
result[which(result$adjPvalues < 0.05),]

## subset data for Patient
PAdata = data[, which(pData(obj)$Patient %in% c('P14','P15'))]
levels(pData(PAdata)$Patient)
pd <- pData(PAdata)
mod <- model.matrix(~1 + Patient, data=pd)
dataRes1 = fitFeatureModel(PAdata, mod)
result <- MRcoefs(dataRes1, taxa=taxa, number = 50)
result[which(result$adjPvalues < 0.05),]


### LOAT 

lodata = data[, which(pData(obj)$LOAT %in% c(7,15))]
levels(pData(lodata)$LOAT)
pd <-pData(lodata)
mod <- model.matrix(~1 + LOAT, data=pd)
dataRes1 = fitFeatureModel(lodata, mod)
result <- MRcoefs(dataRes1, taxa=taxa, number = 1000)
result[which(result$adjPvalues < 0.05),]

  

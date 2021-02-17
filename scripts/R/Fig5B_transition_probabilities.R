library("markovchain")
library(data.table)
library(ggplot2)
library(ggdag)
library(reshape2)
## work in the tmp directory
## save figures to exportDir

exportDir = ''
df = fread('Transitions.csv')
statesNames=c("EC","IC","LC","SC") ### order alphabetically

gut_matrix  = subset(df, df$site =='Gut')
gut_matrix = setDT(dcast(gut_matrix[,2:4], From~To))
gut_matrix
gut_matrix = as.matrix(setnafill(gut_matrix[,2:5], fill=0))
mcCG<-new("markovchain", states=statesNames, transitionMatrix=gut_matrix) 

oral_matrix  = subset(df, df$site =='Oral_cavity')
oral_matrix = setDT(dcast(oral_matrix[,2:4], From~To))
oral_matrix
#oral_matrix$LC = c(0,0,0,0) ### if no transition prob to LC
oral_matrix = as.matrix(setnafill(oral_matrix[,2:5], fill=0))
mcCO<-new("markovchain", states=statesNames, transitionMatrix=oral_matrix) 

skin_matrix  = subset(df, df$site =='Skin')
skin_matrix = setDT(dcast(skin_matrix[,2:4], From~To))
skin_matrix
skin_matrix = as.matrix(setnafill(skin_matrix[,2:5], fill=0))
mcCS<-new("markovchain", states=statesNames, transitionMatrix=skin_matrix) 

SSG<-steadyStates(mcCG)
SSO<-steadyStates(mcCO)
SSS<-steadyStates(mcCS)

SS<-rbind(SSG, SSO, SSS)
rownames(SS)<-c("Gut", "Oral cavity", "Skin")
df = setDT(as.data.frame(SS))
df$Bodysite = row.names(SS)
df = melt(df, variable.name = 'State', id.vars = 'Bodysite')
ggplot(df, aes(x=Bodysite, y=value, fill=State)) + geom_bar(sta='identity', width=.4, color='black') + 
  labs(y='Steady state frequency', x= '') + 
  scale_fill_manual(name='', 
                    values=c("EC"='#fca481',"SC"='#919ebb',
                             "LC"= '#7dc7b0', 'IC'='#eba1cf')) + coord_flip()+
  theme(axis.text.y = element_text(size=16), legend.text = element_text(size=12), axis.title.x = element_text(size=14), 
        aspect.ratio = .4)
ggsave(paste0(exportDir,'SteadyStateFrequencies.svg'))
##############################################################################
cf = fread("ClusterFreq_per_bodysite.csv")
ggplot(cf, aes(x=SampleType, y=ClusterFreqPerSite, fill=ClusterASV)) + geom_col( width=.4, color='black') + 
  labs(y='Empricial Cluster Frequency', x= '') + 
  scale_fill_manual(name='', 
                    values=c("EC"='#fca481',"SC"='#919ebb',
                             "LC"= '#7dc7b0', 'IC'='#eba1cf')) + coord_flip()+
  theme(axis.text.y = element_text(size=16), legend.text = element_text(size=12), axis.title.x = element_text(size=14), 
        aspect.ratio = .4)


#'-------------------------------------------------------------------------------------------------#
#' Title..........: Multi-environment Genomic Prediction using EnvRtype
#' Created........: 2022-18-08
#' Last update... : 2022-28-08
#' Goal...........: Demos on how to process ECs and run genomic prediction
#' Contents.......: (1) data processing demo ; (2) genomic prediction demo
#'                  
#' Author         : Germano Costa-Neto <gmc222@@cornell.edu>
#'-------------------------------------------------------------------------------------------------#

rm(list=ls())


require(EnvRtype)
require(tidyverse)
require(reshape2)
require(ggplot2)

load("./data/example_1_G2F.RData")


ECs = readRDS('./data/ECs_MET1.rds')

ite = 10E2
bur = 2E2
thi = 5

MET_data_set = MET_data1  # phenotypic data
M            = Marker_1   # molecular data


(my_environments = levels(MET_data_set$Environment ))

# standardization to easy the adaptation of the code
names(MET_data_set) = c('env','gid','value')


MxE = matrix(NA,nrow=length(my_environments),ncol=ncol(M))
colnames(MxE) = colnames(M)
rownames(MxE) = my_environments
SE_MxE = MxE

dim(MxE)

require(BGLR)

for(i in 1:length(my_environments))
{
  
  cat(paste0('doing environment....',i,' ',my_environments[i],'\n'))
  
  Y = 
    MET_data_set %>% 
    filter(env %in% my_environments[i]) %>% 
    droplevels()
  
  
  
  gid = levels(Y$gid)
  M_0 = M[rownames(M) %in% gid,]
  
  
  ETA = list(Marker = list(X = M_0 ,model = 'BRR')) # you can change BRR for BL, BayesA, BayesB etc
  fit = BGLR::BGLR(y = Y$value,ETA = ETA,nIter = ite,burnIn = bur,thin = thi,verbose = F)

  
  MxE[i,]    = fit$ETA$Marker$b # marker effects for a given environment
  SE_MxE[i,] = fit$ETA$Marker$SD.b # SE for the marker effects
  
}

MxE[1:5,1:5]



#### Reaction-norm analysis of the marker effects ####
# this is an idea: interpret marker effects using env data
# for this, and as example, I am using PLS (partial least squares)
# however this can be done by RF, GAM, XGBoost entre outros

require(plsdepot)
pls = plsreg2(predictors = ECs,responses = MxE[1:6,1:6],crosval = F)
plot(pls) 
pls$expvar # X = ECs and Y = MxE using 2 latent vectors (LV, similar to principal components)

pls = plsreg2(predictors = ECs,responses = MxE[1:6,1:6],crosval = F,comps = 3)

pls$expvar # X = ECs and Y = MxE using 2 latent vectors (LV, similar to principal components)
pls$std.coefs # betas/intercept ()
pls$reg.coefs # intercept + betas

# example, using only chromsossome 1
markers_chr = map_G2F %>% filter(chr %in% '1')
markers_chr=markers_chr$marker

pls = plsreg2(predictors = ECs,responses = MxE[,which(colnames(MxE) %in% markers_chr)],crosval = F)

#pls$y.pred
#pls$VIP
#pls$std.coefs

heatmap(pls$std.coefs) # how the markers respond to the environment

pls$expvar # X = ECs (69% explained) and Y = MxE (52% explained) using 2 components

sort(pls$VIP[,1],decreasing = T)[1:10] # most important variables in LV 1

sort(pls$VIP[,2],decreasing = T)[1:10] # most important variables in LV 2


STD_Coef = 
  pls$std.coefs %>% melt(varname=c('EC','marker'),value.name = 'reaction_norm') %>% 
  merge(map_G2F,by='marker')
head(STD_Coef)


STD_Coef %>% 
  filter(marker %in% markers_chr) %>% 
  filter(EC %in% names(sort(pls$VIP[,1],decreasing = T)[1:5])) %>% 
  ggplot(aes(x=pos,y=reaction_norm,fill=EC,colour=EC))+stat_smooth(method = 'loess')+
  facet_grid(~EC)


STD_Coef %>% 
  filter(marker %in% markers_chr) %>% 
  filter(EC %in% names(sort(pls$VIP[,1],decreasing = T)[1:5])) %>% 
  ggplot(aes(x=pos,y=reaction_norm,fill=EC,colour=EC))+stat_smooth(method = 'loess')


STD_Coef %>% 
  filter(marker %in% markers_chr) %>% 
  filter(EC %in% names(sort(pls$VIP[,1],decreasing = F)[1:5])) %>% 
  ggplot(aes(x=pos,y=reaction_norm,colour=EC))+stat_smooth(method = 'loess')



## one more example, chr 5

markers_chr = map_G2F %>% filter(chr %in% '5')
markers_chr=markers_chr$marker

pls=plsreg2(predictors = ECs,responses = MxE[,which(colnames(MxE) %in% markers_chr)],crosval = F)

#pls$y.pred
#pls$VIP
#pls$std.coefs

heatmap(pls$std.coefs) # how the markers respond to the environment

pls$expvar # X = ECs (69% explained) and Y = MxE (52% explained) using 2 components

sort(pls$VIP[,1],decreasing = T)[1:10] # most important variables in LV 1

sort(pls$VIP[,2],decreasing = T)[1:10] # most important variables in LV 2


STD_Coef = 
  pls$std.coefs %>% melt(varname=c('EC','marker'),value.name = 'reaction_norm') %>% 
  merge(map_G2F,by='marker')
head(STD_Coef)


# you see here a different pattern than those observed in chr 1
# it seems the same genes are responding to the same variables


STD_Coef %>% 
  filter(marker %in% markers_chr) %>% 
  filter(EC %in% names(sort(pls$VIP[,1],decreasing = T)[1:5])) %>% 
  ggplot(aes(x=pos,y=reaction_norm,fill=EC,colour=EC))+stat_smooth(method = 'loess')+
  facet_grid(~EC)


STD_Coef %>% 
  filter(marker %in% markers_chr) %>% 
  filter(EC %in% names(sort(pls$VIP[,1],decreasing = T)[1:5])) %>% 
  ggplot(aes(x=pos,y=reaction_norm,fill=EC,colour=EC))+stat_smooth(method = 'loess')


STD_Coef %>% 
  filter(marker %in% markers_chr) %>% 
  filter(EC %in% names(sort(pls$VIP[,1],decreasing = F)[1:5])) %>% 
  ggplot(aes(x=pos,y=reaction_norm,colour=EC))+stat_smooth(method = 'loess')


#### lets group markers and envs ###

# again, still using only chr at the time as example
MxE_4_PCA= 
  STD_Coef %>% 
  filter(marker %in% markers_chr) %>% 
  acast(marker~EC,value.var = 'reaction_norm') %>% 
  FactoMineR::PCA(graph = F) %>% 
  FactoMineR::HCPC(nb.clust = -1)

# 3 groups of markers for this chr
# probably is: arms of the chr and centromere
# I am just showing some possibilities of working with this data
# not sure if this will let us for any "new" pathway
MxE_4_PCA$data.clust$clust

EC_4_PCA= 
  STD_Coef %>% 
  filter(marker %in% markers_chr) %>% 
  acast(EC~marker,value.var = 'reaction_norm') %>% 
  FactoMineR::PCA(graph = F) %>% 
  FactoMineR::HCPC(nb.clust = -1)


# 3 groups of covariables 
EC_4_PCA$data.clust$clust



STD_Coef %>% 
  filter(marker %in% markers_chr) %>% 
  filter(EC %in% rownames(EC_4_PCA$data.clust)[EC_4_PCA$data.clust$clust %in% '1']) %>% 
  ggplot(aes(x=pos,y=reaction_norm,fill=EC,colour=EC))+stat_smooth(method = 'loess')+
  theme(legend.position = 'none')

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




(my_environments = levels(MET_data1$Environment))

data("G2F_2014_17")

head(G2F_2014_17)
G2F_2014_17$env = gsub(G2F_2014_17$env,pattern = '@',replacement = '_')

data = G2F_2014_17 %>% filter(env %in% my_environments) %>% droplevels()
data

env   = data$env
lon   = data$lon
lat   = data$lat
env   = data$env
start = data$start
end = data$end
country = rep('USA1',length(lon))


env.data =
  EnvRtype::get_weather(env.id = env,
                        lat = lat,
                        lon = lon,
                        start.day = start,
                        end.day = end,
                        country = country,
                        parallel = TRUE)

head(env.data)

env.data <-   
  env.data %>% 
  param_temperature(merge = T) %>% 
  param_atmospheric(merge=T) %>% 
  param_radiation(merge=T)

head(env.data)

(var_names = names(env.data)[c(2,11:13,15:18,21:23,26,28)])


## I suggest to build a EC matrix considering all possible environments
# Why? Because you are centering your env.data based on the core of environmental conditions
#      the plant breeders will face (and the which are the targets environments)
#      Then, if you are doing one-a-one predictions you will not face that issues for centering
#      mens based on two envs

ECs = EnvRtype::W_matrix(env.data = env.data,var.id = var_names)
heatmap(ECs)

# I suggest not use the "mean" but the quantiles
ECs = EnvRtype::W_matrix(env.data = env.data,var.id = var_names,statistic = 'quantile')
heatmap(ECs)

# and also create development stages. Here I will assume fixed days after emergence
# I suggest you to use your expertise about the growing conditions faced in china and which
# development stages you believe are more important.
# here is an example:

time.intervals = c(0,20,35,65,90)
name_stages    = c('VE-V3','V3-V6','V6-R1','R1-R5','R5-Harvest')

ECs = EnvRtype::W_matrix(env.data = env.data,
                         var.id = var_names,
                         statistic = 'quantile',by.interval = T,
                         time.window = time.intervals,names.window = name_stages)

heatmap(ECs)

# after R5, to the best of my knowledge, the env variables will not help too much
# at least with the maize in Brazil
# so I will remove it

ECs = ECs[,!grepl(colnames(ECs),pattern = 'R5-Harvest')]
heatmap(ECs)

dim(ECs) # 156 environmental covariables!

## once you have the environmental covariables done, you don't need now to create the kenrels
# an interesting fact is: now you can use the ECs for any other appraoch (XGboost, RF, PLS, for instance)
# and, of course, for creating the relative environmental similarity matrices

saveRDS(object = ECs,file = './data/ECs_MET1.rds')

### Organizing the G matrix
require(AGHmatrix) # I like to use this package becuase of the flexibility for creating different matrices

# Example: Additive Relationship Matrix
G_a = Gmatrix(SNPmatrix = Marker_1,method = 'VanRaden',ploidy = 2)

# Example: Dominance Relationship Matrix
G_d = Gmatrix(SNPmatrix = Marker_1,method = 'Vitezica',ploidy = 2)

## as the EC, here I will not create a GRM using all environments because the goal of your research
# is to predict one-to-one environments. I will let it to be built further in this code


# end of the data processing demo
#######################################
#'-------------------------------------------------------------------------------------------------#
# (2) GP demo

# parameters for bayesian models
ite = 10E2
bur = 2E2
thi = 5

MET_data_set = MET_data1  # phenotypic data
M            = Marker_1   # molecular data

# standardization to easy the adaptation of the code
names(MET_data_set) = c('env','gid','value')

(my_environments = levels(MET_data_set$env))


## creating the sets (training/testing environments)
sets = 
  expand.grid(tst=my_environments,trt=my_environments) %>% 
  filter(!tst == trt)


#### Running GP ####

## EXAMPLE Witout the loop

ind = 1
MODEL = 1


tst = sets$tst[ind]
trt = sets$trt[ind]


# selecting phenotypic data
Y = 
  MET_data_set %>% 
  filter(env %in% c(tst,trt)) %>% 
  droplevels()

head(Y)

# selecting environmental data
ECs_0 = ECs[rownames(ECs) %in% c(tst,trt), ]


# selecting molecular data
gid = levels(Y$gid)
M_0 = M[rownames(M) %in% gid,]


G_0 = AGHmatrix::Gmatrix(SNPmatrix = M_0, method = 'VanRaden',ploidy = 2)

ERM_W = env_kernel(env.data = ECs_0,gaussian = TRUE )[[2]] # I suggest to use gaussian = TRUE. But you can test also gaussian = F
# gaussian = FALSE (Jarquin et al., 2014)
# gaussian = TRUE (He et al., 2019; Costa-Neto et al, 2021, 2022 proved that nonlinear is better than linear)
dim(ERM_W)


ERM_I = ERM_W *0+diag(1,nrow = nrow(ERM_W ),ncol = ncol(ERM_W )) # identity matrix


ERM_W_by_stage = env_kernel(env.data = ECs_0,gaussian = TRUE,stages = name_stages[1:4])[[2]]

ERM_W_by_stage$`VE-V3` # as we have only 2 environments, this will not work so well
ERM_W_by_stage$`V3-V6` # however, for more than 2 envs this works nicely because you are able to see the differences in the environmentam similarity due to dev. stages
ERM_W_by_stage$`V6-R1`


# thus, I will use here: Identity matrix (KE_I), W-matrix (ERM_I)
KE_W = list(W = ERM_W)
KE_I = list(E = ERM_I)
KG   = list(G = G_0)



# y = E+G, e~N(0,I)
M01 = EnvRtype::get_kernel(K_E = KE_I,K_G = KG,
                           data = Y,
                           y = 'value',env = 'env',gid = 'gid',
                           model = 'EMM')

# y = E+G+GE, e~N(0,I),  GE = E kronecker G
M02 = EnvRtype::get_kernel(K_E = KE_I,K_G = KG,
                           data = Y,
                           y = 'value',env = 'env',gid = 'gid',
                           model = 'EMDs')

# y = E+G, e~N(0,W)
M03 = EnvRtype::get_kernel(K_E = KE_W,K_G = KG,
                           data = Y,
                           y = 'value',env = 'env',gid = 'gid',
                           model = 'EMM')

# y = E+G+GE, e~N(0,W),  GE = W kronecker G
M04 = EnvRtype::get_kernel(K_E = KE_W,K_G = KG,
                           data = Y,
                           y = 'value',env = 'env',gid = 'gid',
                           model = 'RNMM')

MODEL_LIST = list(M01,M02,M03,M04)
Models = paste0('Model_',1:length(MODEL_LIST))

yNA      <- Y

trainingset = which(yNA$env %in% as.character(trt))
yNA$value[-trainingset] = NA

MODEL = 1 # try model 2, 3, and 4 to see how it works

fit =  EnvRtype::kernel_model(y = 'value',data = yNA,random = MODEL_LIST[[MODEL]],
                              env = 'env',gid = 'gid',tol = 1e-20,
                              iterations = ite,thining = thi,burnin = bur)

fit$VarComp # you can save the VarComps, however I suggest not to run a CV for it
            # more examples in the code "code_GS_envrtype_2.R"

cor(fit$yHat[trainingset],Y$value[trainingset]) # you can save the cor now
# however I prefer to save every pred x obs value for further analysis


output <- data.frame(obs=Y$value,pred=fit$yHat,
                     gid=Y$gid, env=Y$env,
                     Model = Models [MODEL],training_env = trt,testing_env = tst,set=NA)

output$set[ trainingset] <- 'training'
output$set[-trainingset] <- 'testing'

# then I save everything

#return(output)

#### now putting everything with the loop ###


MET_data_set = MET_data1  # phenotypic data
M            = Marker_1   # molecular data

# standardization to easy the adaptation of the code
names(MET_data_set) = c('env','gid','value')

(my_environments = levels(MET_data_set$Environment))


## creating the sets (training/testing environments)
sets = 
  expand.grid(tst=my_environments,trt=my_environments) %>% 
  filter(!tst == trt)

#'-------------------------------------------------------------------------------------------------#
#### Running GP ####

## EXAMPLE Witout the loop
#'-------------------------------------------------------------------------------------------------#
load("./data/example_1_G2F.RData")
ECs = readRDS('./data/ECs_MET1.rds')

# parameters for bayesian models
ite = 10E2
bur = 2E2
thi = 5
seed = 9821

MET_data_set = MET_data1  # phenotypic data
M            = Marker_1   # molecular data

# standardization to easy the adaptation of the code
names(MET_data_set) = c('env','gid','value')

(my_environments = levels(MET_data_set$Environment))


## creating the sets (training/testing environments)
sets = 
  expand.grid(tst=my_environments,trt=my_environments) %>% 
  filter(!tst == trt)

# will run 30 scenarios x 4 models (120 models)
path = './output'
dir.create(path)

output_1 = c()
for(sets_to_run in 1:nrow(sets))
{
  tst = sets$tst[sets_to_run]
  trt = sets$trt[sets_to_run]
  
  
  # selecting phenotypic data
  Y = 
    MET_data_set %>% 
    filter(env %in% c(tst,trt)) %>% 
    droplevels()
  
  head(Y)
  
  # selecting environmental data
  ECs_0 = ECs[rownames(ECs) %in% c(tst,trt), ]
  
  
  # selecting molecular data
  gid = levels(Y$gid)
  M_0 = M[rownames(M) %in% gid,]
  
  
  G_0 = Gmatrix(SNPmatrix = M_0, method = 'VanRaden',ploidy = 2)
  
  ERM_W = env_kernel(env.data = ECs_0,gaussian = TRUE )[[2]] 
  ERM_I = ERM_W *0+diag(1,nrow = nrow(ERM_W ),ncol = ncol(ERM_W )) # identity matrix
  
  
  KE_W = list(W = ERM_W)
  KE_I = list(E = ERM_I)
  KG   = list(G = G_0)
  
  
  # y = E+G, e~N(0,I)
  M01 = EnvRtype::get_kernel(K_E = KE_I,K_G = KG,
                             data = Y,
                             y = 'value',env = 'env',gid = 'gid',
                             model = 'EMM')
  
  # y = E+G+GE, e~N(0,I),  GE = E kronecker G
  M02 = EnvRtype::get_kernel(K_E = KE_I,K_G = KG,
                             data = Y,
                             y = 'value',env = 'env',gid = 'gid',
                             model = 'EMDs')
  
  # y = E+G, e~N(0,W)
  M03 = EnvRtype::get_kernel(K_E = KE_W,K_G = KG,
                             data = Y,
                             y = 'value',env = 'env',gid = 'gid',
                             model = 'EMM')
  
  # y = E+G+GE, e~N(0,W),  GE = W kronecker G
  M04 = EnvRtype::get_kernel(K_E = KE_W,K_G = KG,
                             data = Y,
                             y = 'value',env = 'env',gid = 'gid',
                             model = 'RNMM')
  
  MODEL_LIST = list(M01,M02,M03,M04)
  Models = paste0('Model_',1:length(MODEL_LIST))
  
  yNA      <- Y
  
  trainingset = which(yNA$env %in% as.character(trt))
  yNA$value[-trainingset] = NA
  
  output_0 = c()
  
  for(MODEL in 1:length(Models))
  {
    set.seed(seed)
    fit =  EnvRtype::kernel_model(y = 'value',data = yNA,random = MODEL_LIST[[MODEL]],
                                  env = 'env',gid = 'gid',#tol = 1e-20,
                                  iterations = ite,thining = thi,burnin = bur)
    
    output <- data.frame(obs=Y$value,pred=fit$yHat,
                         gid=Y$gid, env=Y$env,
                         Model = Models [MODEL],
                         training_env = trt,testing_env = tst,set=NA)
    
    output$set[ trainingset] <- 'training'
    output$set[-trainingset] <- 'testing'
    
    write.table(x = data.frame(Model=Models [MODEL],training_env = trt,testing_env = tst,
                               r = cor(output$obs[-trainingset],output$pred[-trainingset])),
                file ='predictive_ability.txt',append = T )
    
    saveRDS(output,file = paste0(path,'/',paste(Models [MODEL],training_env = trt,testing_env = tst,sep = '@')))
    output_0 = rbind(output_0,output)
  }
  output_1 = rbind(output_0,output_1)

}

require(plyr)

mse = function(x,y) sum((x-y)^2)/length(y)

output_1 = c()
for(i in 1:length(list.files(path = path)))
{
  output_1=rbind(output_1,readRDS(paste0('./output/',list.files(path = path)[i])))
}

output_1 %>% filter(set %in% 'testing') %>% ggplot(aes(x=obs,y=pred,colour=training_env))+geom_point()+
  facet_grid(training_env~Model)


output_1 %>% 
  ddply(.(Model,training_env,testing_env,set),summarise,
        pa= round(cor(obs,pred),3),
        mse = round(mse(obs,pred),3)) %>%
  filter(set %in% 'testing') %>% 
  ggplot(aes(x=training_env, y=testing_env,fill=pa))+
  geom_tile()+facet_grid(~Model)+
  xlab('Training Env')+ylab("Testing Env")+
  scale_fill_gradientn(colours = rainbow(5))+
  theme(axis.text.x = element_text(angle=90,hjust=1))




## OBS: THis funciton works well for gathering data for multiple environments
# I guess for two-by-two this is not the best situation
# Also, this data set has a more simple GxE pattern -- so it is expect that the use of ECs
# will not increase too much accuracy in GxE prediction

# However we can replace the kernel_model() by a BGLR funcion


head(output_1)
    
### END ##################


## you can also run it in parallel. However I have not organized it yet for this code
# there is better ways for doing it. I guess the code below is not working, but can be modified for someone with more time 
# I aways use it, just like i did here: https://raw.githubusercontent.com/allogamous/EnvRtype/master/EnvRtype_full_tutorial.R


require(doParallel)
require(foreach)

cl <- makeCluster(4)
registerDoParallel(cl)


results <-
  foreach(sets_to_run = 1:length(sets$tst),.packages = c('tidyverse','AGHmatrix','reshape2'), .combine = "rbind")%:%
  foreach(MODEL = 1:4, .packages = c('tidyverse','AGHmatrix','reshape2'),.combine = "rbind") %dopar% # I will run 4 models
  {
    tst = sets$tst[sets_to_run]
    trt = sets$trt[sets_to_run]
    
    
    # selecting phenotypic data
    Y = 
      MET_data_set %>% 
      filter(env %in% c(tst,trt)) %>% 
      droplevels()
    
    head(Y)
    
    # selecting environmental data
    ECs_0 = ECs[rownames(ECs) %in% c(tst,trt), ]
    
    
    # selecting molecular data
    gid = levels(Y$gid)
    M_0 = M[rownames(M) %in% gid,]
    
    
    G_0 = Gmatrix(SNPmatrix = M_0, method = 'VanRaden',ploidy = 2)
    
    ERM_W = env_kernel(env.data = ECs_0,gaussian = TRUE )[[2]] 
    ERM_I = ERM_W *0+diag(1,nrow = nrow(ERM_W ),ncol = ncol(ERM_W )) # identity matrix
    
    
    KE_W = list(W = ERM_W)
    KE_I = list(E = ERM_I)
    KG   = list(G = G_0)
    
    
    # y = E+G, e~N(0,I)
    M01 = EnvRtype::get_kernel(K_E = KE_I,K_G = KG,
                               data = Y,
                               y = 'value',env = 'env',gid = 'gid',
                               model = 'EMM')
    
    # y = E+G+GE, e~N(0,I),  GE = E kronecker G
    M02 = EnvRtype::get_kernel(K_E = KE_I,K_G = KG,
                               data = Y,
                               y = 'value',env = 'env',gid = 'gid',
                               model = 'EMDs')
    
    # y = E+G, e~N(0,W)
    M03 = EnvRtype::get_kernel(K_E = KE_W,K_G = KG,
                               data = Y,
                               y = 'value',env = 'env',gid = 'gid',
                               model = 'EMM')
    
    # y = E+G+GE, e~N(0,W),  GE = W kronecker G
    M04 = EnvRtype::get_kernel(K_E = KE_W,K_G = KG,
                               data = Y,
                               y = 'value',env = 'env',gid = 'gid',
                               model = 'RNMM')
    
    MODEL_LIST = list(M01,M02,M03,M04)
    Models = paste0('Model_',1:length(MODEL_LIST))
    
    yNA      <- Y
    
    trainingset = which(yNA$env %in% as.character(trt))
    yNA$value[-trainingset] = NA
    
    set.seed(seed)
    fit =  EnvRtype::kernel_model(y = 'value',data = yNA,random = MODEL_LIST[[MODEL]],
                                  env = 'env',gid = 'gid',tol = 1e-20,
                                  iterations = ite,thining = thi,burnin = bur)
    
    output <- data.frame(obs=Y$value,pred=fit$yHat,
                         gid=Y$gid, env=Y$env,
                         Model = Models [MODEL],
                         training_env = trt,testing_env = tst,set=NA)
    
    output$set[ trainingset] <- 'training'
    output$set[-trainingset] <- 'testing'
    
    write.table(x = data.frame(Model=Models [MODEL],training_env = trt,testing_env = tst,
                               r = cor(output$obs[trainingset],output$pred[trainingset])),
                file ='predictive_ability.txt',append = T )
    
    saveRDS(output,file = paste0(path,'/',paste(Models [MODEL],training_env = trt,testing_env = tst,sep = '@')))
    
    
    return(output)
  }
    
stopCluster(cl)

# this is a more complex data set (diverse years).You can also test it, 
#   just replace Marker_1 by Marker_2 and MET_data1 by MET_data2
load("./data/example_2_G2F.RData")


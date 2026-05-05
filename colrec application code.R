library(survivalNET)
library(survivalPLANN)
library(relsurv)
library(parallel)

library(ggplot2)
library(dplyr)
library(tidyr)

data(colrec)
data(slopop)
colrec$agey <- colrec$age/365.241

colrec$sexchara <- "male"
colrec$sexchara[colrec$sex==2] <- "female"
colrec$sex01 <- (colrec$sex)-1

colrec$colon <- 1*(as.character(colrec$site)=="colon")

colrec <- colrec[colrec$stage!=99,] # suppresion des patients avec valeur manquante

colrec$stage1 <- 1*(colrec$stage==1)
colrec$stage2 <- 1*(colrec$stage==2)
colrec$stage3 <- 1*(colrec$stage==3)

table(colrec$stat)

colrec$agey10 <- colrec$agey/10

colrec$sex.organ <- paste0(colrec$sexchara, colrec$colon)
colrec$sex.organ <- as.numeric(factor(colrec$sex.organ, levels = c("male0", "male1", "female0", "female1")))

### Objectif : prédire la survie nette en fonction des covariables (âge, sexe, localisation)

##  On ne comparera pas à PP car ce n'est qu'un estimateur et pas un modèle

######## ATTENTION, colrec à 5971 data et non juste 5578 on a supprimé les individus sans données sur une colonne 

set.seed(203512)

ind <- rbinom(dim(colrec)[1],1,0.5)

data_train <- colrec[ind ==0,] 
data_valid <- colrec[ind ==1,] 

pro.time <- max(data_train$time)



cv_para <- function(i){
  set.seed(110226 +i)
  if(i ==1){
    
    cv <- cvFLEXNET(formula = Surv(time, stat) ~ stage2 + stage3 + agey10 + sex01 + colon+
                      ratetable(age, diag, sexchara), 
                    ratetable=slopop , pro.time = pro.time,
                    data = data_train, cv = 30, m = c(0,1,2,3))
    
    
  }else if(i ==2){
    
    cv <- cvFLEXNET(formula = Surv(time, stat) ~ stage2 + stage3 + agey10 + strata(sex.organ)+
                      ratetable(age, diag, sexchara), ratetable=slopop , pro.time = pro.time,
                    data = data_train, cv = 30, m = c(0,1,2,3), m_s = c(0,1))
    
  }else if(i ==3){
    
    
    cv <- cvPLANN(formula = Surv(time, stat) ~ stage2 + stage3 + agey10 + sex01 + colon, pro.time = pro.time,
                  data = data_train, cv = 30, inter= 365.241/12, size=c(6,12,24,32),
                  decay=c(0.01, 0.05, 0.1), maxit=1000, MaxNWts=10000, metric = "ibs")
    
  }
  return(cv)
}

res_cv <- mclapply(1:3, cv_para, mc.cores = 3)
names(res_cv) <- c("flexnet","flexnet2","plann")
######### paramètres des modèles 

tune.flex <- res_cv$flexnet
tune.flex2 <- res_cv$flexnet2
tune.plann <- res_cv$plann

##flex1 

m1 <- tune.flex$optimal$m

##flex2 

m2 <- tune.flex2$optimal$m
m_s <- tune.flex2$optimal$m_s

## PLANN

inter <- tune.plann$optimal$inter

size <- tune.plann$optimal$size

decay <- tune.plann$optimal$decay

maxit <- tune.plann$optimal$maxit

MaxNWts <- tune.plann$optimal$MaxNWts 

save(list = ls(), file = paste0("/home/thomas/Documents/Rstudio/Application/Application colrec/CV_4strates_save.Rdata"))

################## Estimation des modèles 

## flex1

flex.model1 <- survivalFLEXNET(
  formula = Surv(time, stat) ~ stage2 + stage3 + agey10 + sex01 + colon +
    ratetable(age, diag, sexchara),
  data = data_train,
  ratetable = slopop,
  m = m1
)
##flex2
flex.model2 <-survivalFLEXNET(formula = Surv(time, stat) ~ stage2 + stage3 + agey10 + strata(sex.organ) +
                                ratetable(age, diag, sexchara), data = data_train,
                              ratetable=slopop, m = m2, m_s = m_s)
##PLANN

plann.model <- sPLANN(formula = Surv(time, stat) ~ stage2 + stage3 + agey10 + sex01 + colon,
                      pro.time = pro.time, data = data_train, inter= inter, size=size, 
                      decay=decay, maxit=maxit, MaxNWts=MaxNWts)

#################

newtimes = c(3,5,10)*365.241
newtimespred <- sort(unique(c(seq(0, 365.241*10, 365.241/12), newtimes)))[-1]
#### prédictions

plannpred <- predictRS(plann.model, data = data_train, newtimes = newtimespred, ratetable=slopop, 
                       age = "age", year = "diag", sex = "sexchara")

mean.plann <- plannpred$mpredictions$net_survival 

###flexnet

### modèle PH
##2 noeuds

flexpred1 <- predict(flex.model1, newtimes = newtimespred)$predictions

mean.flex1 <- apply(flexpred1, FUN="mean", MARGIN=2)


#### modèle NPH

##2 noeuds

flexpred2 <- predict(flex.model2 , newtimes = newtimespred)$predictions

mean.flex2  <- apply(flexpred2 , FUN="mean", MARGIN=2)

## pohar perme

PP <- summary( rs.surv(Surv(time, stat) ~ 1, data = data_train,
                       ratetable = slopop, method = "pohar-perme", add.times = newtimespred,
                       rmap = list(age = age, sex = sexchara, year = diag)), times = newtimespred, extend = TRUE) 

poharpred <- PP$surv 

##### Indicateurs

## ROC.net

library(RISCA)
##choix de l'espacement des cutoffs 
# c.off <- seq(0,1, by = 0.2)
# sort(unique(1-ind_estimP[, l + 1]))
# unique( quantile( 1-ind_estimP[, l + 1],  probs= seq(0,1, by = 0.01)) ) #pas en dessous de 20 points 0.025

##### calcul de ROC net


datatrain <- data_train
ind_estimP <- plannpred$ipredictions$relative_survival[,-1]
ind_estimF1.2 <- flexpred1
ind_estimF2.2 <- flexpred2


# Initialiser les vecteurs de résultats
hold_P <- hold_F1.2 <- hold_F2.2 <- c()
holdmod_P <- holdmod_F1.2 <- holdmod_F2.2 <- list()


library(doParallel)
library(foreach)

cl <- makeCluster(parallel::detectCores() - 5)
registerDoParallel(cl)

res <- foreach(l = seq_along(newtimes)) %dopar% {
  library(RISCA)
  
  rocP <- roc.net(datatrain$time, datatrain$stat, 1-ind_estimP[, l + 1],
                  datatrain$age, datatrain$sexchara, datatrain$diag, slopop,
                  pro.time = newtimes[l],
                  cut.off = unique(quantile(1-ind_estimP[, l + 1],
                                            probs = seq(0, 1, .01))))
  
  rocF1.2 <- roc.net(datatrain$time, datatrain$stat, 1-ind_estimF1.2[, l],
                     datatrain$age, datatrain$sexchara, datatrain$diag,
                     slopop, pro.time = newtimes[l],
                     cut.off = unique(quantile(1-ind_estimF1.2[, l],
                                               probs = seq(0, 1, .01))))
  
  rocF2.2 <- roc.net(datatrain$time, datatrain$stat, 1-ind_estimF2.2[, l],
                     datatrain$age, datatrain$sexchara, datatrain$diag,
                     slopop, pro.time = newtimes[l],
                     cut.off = unique(quantile(1-ind_estimF2.2[, l],
                                               probs = seq(0, 1, .01))))
  
  list(
    aucP = rocP$auc,
    aucF1.2 = rocF1.2$auc,
    aucF2.2 = rocF2.2$auc,
    rocP = rocP,
    rocF1.2 = rocF1.2,
    rocF2.2 = rocF2.2
  )
}

stopCluster(cl)

hold_P <- c(res[[1]]$aucP, res[[2]]$aucP, res[[3]]$aucP) 
hold_F1.2 <- c(res[[1]]$aucF1.2, res[[2]]$aucF1.2, res[[3]]$aucF1.2) 
hold_F2.2 <- c(res[[1]]$aucF2.2, res[[2]]$aucF2.2, res[[3]]$aucF2.2) 

for(i in 1:length(newtimes)){
  holdmod_P[[i]] <- res[[i]]$rocP
  holdmod_F1.2[[i]] <- res[[i]]$rocF1.2
  holdmod_F2.2[[i]] <- res[[i]]$rocF2.2
  
}
results_list = list(
  P     = hold_P,
  F1.2  = hold_F1.2,
  F2.2  = hold_F2.2
)

results_mod_list = list(
  P     = holdmod_P,
  F1.2  = holdmod_F1.2,
  F2.2  = holdmod_F2.2
)
# Transformer en matrices pour chaque méthode :
AUC_WHOLE_P    <- results_list$P
AUC_WHOLE_F1.2 <- results_list$F1.2
AUC_WHOLE_F2.2 <- results_list$F2.2

save(list = ls(), file = paste0("/home/thomas/Documents/Rstudio/Application/Application colrec/AUC_train_save.Rdata"))

############# PAR STRATES

strata_names = c("HC", "HR", "FC", "FR")

for(j in strata_names){
  if(j == "HC"){
    data_s <- data_train[data_train$sex ==1 & data_train$site == "colon",  ]
  }else if(j == "HR"){
    data_s <- data_train[data_train$sex ==1 & data_train$site == "rectum",  ]
    
  }else if(j == "FC"){
    data_s <- data_train[data_train$sex ==2 & data_train$site == "colon",  ]
  }else if(j == "FR"){
    data_s <- data_train[data_train$sex ==2 & data_train$site == "rectum",  ]
  }
  
  plannpred_s <- predictRS(plann.model, data = data_s, newtimes = newtimespred, ratetable=slopop, 
                           age = "age", year = "diag", sex = "sexchara")
  
  mean.plann_s <- plannpred_s$mpredictions$net_survival 
  
  ###flexnet
  
  ### modèle PH
  ##2 noeuds
  
  flexpred1_s <- predict(flex.model1, newtimes = newtimespred, newdata = data_s)$predictions
  
  mean.flex1_s <- apply(flexpred1_s, FUN="mean", MARGIN=2)
  
  
  #### modèle NPH
  
  ##2 noeuds
  
  flexpred2_s <- predict(flex.model2 , newtimes = newtimespred, newdata = data_s)$predictions
  
  mean.flex2_s  <- apply(flexpred2_s , FUN="mean", MARGIN=2)
  
  ## pohar perme
  
  PP_s <- summary( rs.surv(Surv(time, stat) ~ 1, data = data_s,
                           ratetable = slopop, method = "pohar-perme", add.times = newtimespred,
                           rmap = list(age = age, sex = sexchara, year = diag)), times = newtimespred, extend = TRUE) 
  
  poharpred_s <- PP_s$surv 
  
  if(j == "HC"){
    plannpred_HC <<- plannpred_s
    mean.plann_HC <<- mean.plann_s
    
    flexpred1_HC <<- flexpred1_s
    mean.flex1_HC <<- mean.flex1_s
    
    flexpred2_HC <<- flexpred2_s
    mean.flex2_HC <<- mean.flex2_s
    
    PP_HC <<- PP_s
    poharpred_HC <<- poharpred_s
  }else if(j == "HR"){
    plannpred_HR <<- plannpred_s
    mean.plann_HR <<- mean.plann_s
    
    flexpred1_HR <<- flexpred1_s
    mean.flex1_HR <<- mean.flex1_s
    
    flexpred2_HR <<- flexpred2_s
    mean.flex2_HR <<- mean.flex2_s
    
    PP_HR <<- PP_s
    poharpred_HR <<- poharpred_s
  }else if(j == "FC"){
    plannpred_FC <<- plannpred_s
    mean.plann_FC <<- mean.plann_s
    
    flexpred1_FC <<- flexpred1_s
    mean.flex1_FC <<- mean.flex1_s
    
    flexpred2_FC <<- flexpred2_s
    mean.flex2_FC <<- mean.flex2_s
    
    PP_FC <<- PP_s
    poharpred_FC <<- poharpred_s
  }else if(j == "FR"){
    plannpred_FR <<- plannpred_s
    mean.plann_FR <<- mean.plann_s
    
    flexpred1_FR <<- flexpred1_s
    mean.flex1_FR <<- mean.flex1_s
    
    flexpred2_FR <<- flexpred2_s
    mean.flex2_FR <<- mean.flex2_s
    
    PP_FR <<- PP_s
    poharpred_FR <<- poharpred_s
  }
  
  ind_estimP <- plannpred_s$ipredictions$relative_survival[,-1]
  ind_estimF1.2 <- flexpred1_s
  ind_estimF2.2 <- flexpred2_s
  
  
  # Initialiser les vecteurs de résultats
  hold_P <- hold_F1.2 <- hold_F2.2 <- c()
  holdmod_P <- holdmod_F1.2 <- holdmod_F2.2 <- list()
  
  
  library(doParallel)
  library(foreach)
  
  cl <- makeCluster(parallel::detectCores() - 5)
  registerDoParallel(cl)
  
  res_s <- foreach(l = seq_along(newtimes)) %dopar% {
    library(RISCA)
    
    rocP <- roc.net(data_s$time, data_s$stat, 1-ind_estimP[, l + 1],
                    data_s$age, data_s$sexchara, data_s$diag, slopop,
                    pro.time = newtimes[l],
                    cut.off = unique(quantile(1-ind_estimP[, l + 1],
                                              probs = seq(0, 1, .01))))
    
    rocF1.2 <- roc.net(data_s$time, data_s$stat, 1-ind_estimF1.2[, l],
                       data_s$age, data_s$sexchara, data_s$diag,
                       slopop, pro.time = newtimes[l],
                       cut.off = unique(quantile(1-ind_estimF1.2[, l],
                                                 probs = seq(0, 1, .01))))
    
    rocF2.2 <- roc.net(data_s$time, data_s$stat, 1-ind_estimF2.2[, l],
                       data_s$age, data_s$sexchara, data_s$diag,
                       slopop, pro.time = newtimes[l],
                       cut.off = unique(quantile(1-ind_estimF2.2[, l],
                                                 probs = seq(0, 1, .01))))
    
    list(
      aucP = rocP$auc,
      aucF1.2 = rocF1.2$auc,
      aucF2.2 = rocF2.2$auc,
      rocP = rocP,
      rocF1.2 = rocF1.2,
      rocF2.2 = rocF2.2
    )
  }
  
  stopCluster(cl)
  
  hold_P <- c(res_s[[1]]$aucP, res_s[[2]]$aucP, res_s[[3]]$aucP) 
  hold_F1.2 <- c(res_s[[1]]$aucF1.2, res_s[[2]]$aucF1.2, res_s[[3]]$aucF1.2) 
  hold_F2.2 <- c(res_s[[1]]$aucF2.2, res_s[[2]]$aucF2.2, res_s[[3]]$aucF2.2) 
  
  for(i in 1:length(newtimes)){
    holdmod_P[[i]] <- res_s[[i]]$rocP
    holdmod_F1.2[[i]] <- res_s[[i]]$rocF1.2
    holdmod_F2.2[[i]] <- res_s[[i]]$rocF2.2
    
  }
  results_list_s = list(
    P     = hold_P,
    F1.2  = hold_F1.2,
    F2.2  = hold_F2.2
  )
  
  results_mod_list_s = list(
    P     = holdmod_P,
    F1.2  = holdmod_F1.2,
    F2.2  = holdmod_F2.2
  )
  
  if( j == "HC"){
    results_list_HC <- results_list_s
    results_mod_list_HC <- results_mod_list_s
  }else if( j == "HR"){
    results_list_HR <- results_list_s
    results_mod_list_HR <- results_mod_list_s
  }else if( j == "FC"){
    results_list_FC <- results_list_s
    results_mod_list_FC <- results_mod_list_s
  }else if( j == "FR"){
    results_list_FR <- results_list_s
    results_mod_list_FR <- results_mod_list_s
  }
}
# Transformer en matrices pour chaque méthode :
AUC_WHOLE_P_HC    <- results_list_HC$P
AUC_WHOLE_F1.2_HC <- results_list_HC$F1.2
AUC_WHOLE_F2.2_HC <- results_list_HC$F2.2

AUC_WHOLE_P_HR    <- results_list_HR$P
AUC_WHOLE_F1.2_HR <- results_list_HR$F1.2
AUC_WHOLE_F2.2_HR <- results_list_HR$F2.2

AUC_WHOLE_P_FC    <- results_list_FC$P
AUC_WHOLE_F1.2_FC <- results_list_FC$F1.2
AUC_WHOLE_F2.2_FC <- results_list_FC$F2.2

AUC_WHOLE_P_FR    <- results_list_FR$P
AUC_WHOLE_F1.2_FR <- results_list_FR$F1.2
AUC_WHOLE_F2.2_FR <- results_list_FR$F2.2


save(list = ls(), file = paste0("/home/thomas/Documents/Rstudio/Application/Application colrec/AUC_train_strata_save.Rdata"))

#################################################################################################################################
##                                           VALIDATION 
#################################################################################################################################

plannpredval <- predictRS(plann.model, data = data_valid, newtimes = newtimespred, ratetable=slopop, 
                          age = "age", year = "diag", sex = "sexchara")

mean.plannval <- plannpredval$mpredictions$net_survival 

###flexnet

### modèle PH
##2 noeuds

flexpredval1 <- predict(flex.model1, newdata = data_valid ,newtimes = newtimespred)$predictions

mean.flexval1 <- apply(flexpredval1, FUN="mean", MARGIN=2)


#### modèle NPH

##2 noeuds

flexpredval2 <- predict(flex.model2 , newdata = data_valid, newtimes = newtimespred)$predictions

mean.flexval2  <- apply(flexpredval2 , FUN="mean", MARGIN=2)

## pohar perme

PPval <- summary( rs.surv(Surv(time, stat) ~ 1, data = data_valid,
                          ratetable = slopop, method = "pohar-perme", add.times = newtimespred,
                          rmap = list(age = age, sex = sexchara, year = diag)), times = newtimespred, extend = TRUE) 

poharpredval <- PPval$surv 

##### Indicateurs

## ROC.net

library(RISCA)
##choix de l'espacement des cutoffs 
# c.off <- seq(0,1, by = 0.2)
# sort(unique(1-ind_estimP[, l + 1]))
# unique( quantile( 1-ind_estimP[, l + 1],  probs= seq(0,1, by = 0.01)) ) #pas en dessous de 20 points 0.025

##### calcul de ROC net


datatrain <- data_valid
ind_estimP <- plannpredval$ipredictions$relative_survival
ind_estimF1.2 <- flexpredval1
ind_estimF2.2 <- flexpredval2


# Initialiser les vecteurs de résultats
hold_P <- hold_F1.2 <- hold_F2.2 <- c()
holdmod_P <- holdmod_F1.2 <- holdmod_F2.2 <- list()

cl <- makeCluster(parallel::detectCores() - 5)
registerDoParallel(cl)

res_val <- foreach(l = seq_along(newtimes)) %dopar% {
  library(RISCA)
  
  rocP <- roc.net(datatrain$time, datatrain$stat, 1-ind_estimP[, l + 1],
                  datatrain$age, datatrain$sexchara, datatrain$diag, slopop,
                  pro.time = newtimes[l],
                  cut.off = unique(quantile(1-ind_estimP[, l + 1],
                                            probs = seq(0, 1, .01))))
  
  rocF1.2 <- roc.net(datatrain$time, datatrain$stat, 1-ind_estimF1.2[, l],
                     datatrain$age, datatrain$sexchara, datatrain$diag,
                     slopop, pro.time = newtimes[l],
                     cut.off = unique(quantile(1-ind_estimF1.2[, l],
                                               probs = seq(0, 1, .01))))
  
  rocF2.2 <- roc.net(datatrain$time, datatrain$stat, 1-ind_estimF2.2[, l],
                     datatrain$age, datatrain$sexchara, datatrain$diag,
                     slopop, pro.time = newtimes[l],
                     cut.off = unique(quantile(1-ind_estimF2.2[, l],
                                               probs = seq(0, 1, .01))))
  
  list(
    aucP = rocP$auc,
    aucF1.2 = rocF1.2$auc,
    aucF2.2 = rocF2.2$auc,
    rocP = rocP,
    rocF1.2 = rocF1.2,
    rocF2.2 = rocF2.2
  )
}

stopCluster(cl)

hold_P <- c(res_val[[1]]$aucP, res_val[[2]]$aucP, res_val[[3]]$aucP) 
hold_F1.2 <- c(res_val[[1]]$aucF1.2, res_val[[2]]$aucF1.2, res_val[[3]]$aucF1.2) 
hold_F2.2 <- c(res_val[[1]]$aucF2.2, res_val[[2]]$aucF2.2, res_val[[3]]$aucF2.2) 

for(i in 1:length(newtimes)){
  holdmod_P[[i]] <- res_val[[i]]$rocP
  holdmod_F1.2[[i]] <- res_val[[i]]$rocF1.2
  holdmod_F2.2[[i]] <- res_val[[i]]$rocF2.2
  
}

resultsval_list = list(
  P     = hold_P,
  F1.2  = hold_F1.2,
  F2.2  = hold_F2.2
)
resultsval_mod_list = list(
  P     = holdmod_P,
  F1.2  = holdmod_F1.2,
  F2.2  = holdmod_F2.2
)
# Transformer en matrices pour chaque méthode :
AUC_WHOLEval_P    <- resultsval_list$P
AUC_WHOLEval_F1.2 <- resultsval_list$F1.2
AUC_WHOLEval_F2.2 <- resultsval_list$F2.2

save(list = ls(), file = paste0("/home/thomas/Documents/Rstudio/Application/Application colrec/AUC_valid_save.Rdata"))

############# PAR STRATES



for(j in strata_names){
  if(j == "HC"){
    data_s <- data_valid[data_valid$sex ==1 & data_valid$site == "colon",  ]
  }else if(j == "HR"){
    data_s <- data_valid[data_valid$sex ==1 & data_valid$site == "rectum",  ]
    
  }else if(j == "FC"){
    data_s <- data_valid[data_valid$sex ==2 & data_valid$site == "colon",  ]
  }else if(j == "FR"){
    data_s <- data_valid[data_valid$sex ==2 & data_valid$site == "rectum",  ]
  }
  
  plannpred_s <- predictRS(plann.model, data = data_s, newtimes = newtimespred, ratetable=slopop, 
                           age = "age", year = "diag", sex = "sexchara")
  
  mean.plann_s <- plannpred_s$mpredictions$net_survival 
  
  ###flexnet
  
  ### modèle PH
  ##2 noeuds
  
  flexpred1_s <- predict(flex.model1, newtimes = newtimespred, newdata = data_s)$predictions
  
  mean.flex1_s <- apply(flexpred1_s, FUN="mean", MARGIN=2)
  
  
  #### modèle NPH
  
  ##2 noeuds
  
  flexpred2_s <- predict(flex.model2 , newtimes = newtimespred, newdata = data_s)$predictions
  
  mean.flex2_s  <- apply(flexpred2_s , FUN="mean", MARGIN=2)
  
  ## pohar perme
  
  PP_s <- summary( rs.surv(Surv(time, stat) ~ 1, data = data_s,
                           ratetable = slopop, method = "pohar-perme", add.times = newtimespred,
                           rmap = list(age = age, sex = sexchara, year = diag)), times = newtimespred, extend = TRUE) 
  
  poharpred_s <- PP_s$surv 
  
  if(j == "HC"){
    plannpredval_HC <<- plannpred_s
    mean.plannval_HC <<- mean.plann_s
    
    flexpred1val_HC <<- flexpred1_s
    mean.flex1val_HC <<- mean.flex1_s
    
    flexpred2val_HC <<- flexpred2_s
    mean.flex2val_HC <<- mean.flex2_s
    
    PPval_HC <<- PP_s
    poharpredval_HC <<- poharpred_s
  }else if(j == "HR"){
    plannpredval_HR <<- plannpred_s
    mean.plannval_HR <<- mean.plann_s
    
    flexpred1val_HR <<- flexpred1_s
    mean.flex1val_HR <<- mean.flex1_s
    
    flexpred2val_HR <<- flexpred2_s
    mean.flex2val_HR <<- mean.flex2_s
    
    PPval_HR <<- PP_s
    poharpredval_HR <<- poharpred_s
  }else if(j == "FC"){
    plannpredval_FC <<- plannpred_s
    mean.plannval_FC <<- mean.plann_s
    
    flexpred1val_FC <<- flexpred1_s
    mean.flex1val_FC <<- mean.flex1_s
    
    flexpred2val_FC <<- flexpred2_s
    mean.flex2val_FC <<- mean.flex2_s
    
    PPval_FC <<- PP_s
    poharpredval_FC <<- poharpred_s
  }else if(j == "FR"){
    plannpredval_FR <<- plannpred_s
    mean.plannval_FR <<- mean.plann_s
    
    flexpred1val_FR <<- flexpred1_s
    mean.flex1val_FR <<- mean.flex1_s
    
    flexpred2val_FR <<- flexpred2_s
    mean.flex2val_FR <<- mean.flex2_s
    
    PPval_FR <<- PP_s
    poharpredval_FR <<- poharpred_s
  }
  
  ind_estimP <- plannpred_s$ipredictions$relative_survival[,-1]
  ind_estimF1.2 <- flexpred1_s
  ind_estimF2.2 <- flexpred2_s
  
  
  # Initialiser les vecteurs de résultats
  hold_P <- hold_F1.2 <- hold_F2.2 <- c()
  holdmod_P <- holdmod_F1.2 <- holdmod_F2.2 <- list()
  
  
  library(doParallel)
  library(foreach)
  
  cl <- makeCluster(parallel::detectCores() - 5)
  registerDoParallel(cl)
  
  res_s <- foreach(l = seq_along(newtimes)) %dopar% {
    library(RISCA)
    
    rocP <- roc.net(data_s$time, data_s$stat, 1-ind_estimP[, l + 1],
                    data_s$age, data_s$sexchara, data_s$diag, slopop,
                    pro.time = newtimes[l],
                    cut.off = unique(quantile(1-ind_estimP[, l + 1],
                                              probs = seq(0, 1, .01))))
    
    rocF1.2 <- roc.net(data_s$time, data_s$stat, 1-ind_estimF1.2[, l],
                       data_s$age, data_s$sexchara, data_s$diag,
                       slopop, pro.time = newtimes[l],
                       cut.off = unique(quantile(1-ind_estimF1.2[, l],
                                                 probs = seq(0, 1, .01))))
    
    rocF2.2 <- roc.net(data_s$time, data_s$stat, 1-ind_estimF2.2[, l],
                       data_s$age, data_s$sexchara, data_s$diag,
                       slopop, pro.time = newtimes[l],
                       cut.off = unique(quantile(1-ind_estimF2.2[, l],
                                                 probs = seq(0, 1, .01))))
    
    list(
      aucP = rocP$auc,
      aucF1.2 = rocF1.2$auc,
      aucF2.2 = rocF2.2$auc,
      rocP = rocP,
      rocF1.2 = rocF1.2,
      rocF2.2 = rocF2.2
    )
  }
  
  stopCluster(cl)
  
  hold_P <- c(res_s[[1]]$aucP, res_s[[2]]$aucP, res_s[[3]]$aucP) 
  hold_F1.2 <- c(res_s[[1]]$aucF1.2, res_s[[2]]$aucF1.2, res_s[[3]]$aucF1.2) 
  hold_F2.2 <- c(res_s[[1]]$aucF2.2, res_s[[2]]$aucF2.2, res_s[[3]]$aucF2.2) 
  
  for(i in 1:length(newtimes)){
    holdmod_P[[i]] <- res_s[[i]]$rocP
    holdmod_F1.2[[i]] <- res_s[[i]]$rocF1.2
    holdmod_F2.2[[i]] <- res_s[[i]]$rocF2.2
    
  }
  resultsval_list_s = list(
    P     = hold_P,
    F1.2  = hold_F1.2,
    F2.2  = hold_F2.2
  )
  
  resultsval_mod_list_s = list(
    P     = holdmod_P,
    F1.2  = holdmod_F1.2,
    F2.2  = holdmod_F2.2
  )
  
  if( j == "HC"){
    resultsval_list_HC <- resultsval_list_s
    resultsval_mod_list_HC <- resultsval_mod_list_s
  }else if( j == "HR"){
    resultsval_list_HR <- resultsval_list_s
    resultsval_mod_list_HR <- resultsval_mod_list_s
  }else if( j == "FC"){
    resultsval_list_FC <- resultsval_list_s
    resultsval_mod_list_FC <- resultsval_mod_list_s
  }else if( j == "FR"){
    resultsval_list_FR <- resultsval_list_s
    resultsval_mod_list_FR <- resultsval_mod_list_s
  }
}
# Transformer en matrices pour chaque méthode :
AUC_WHOLEval_P_HC    <- resultsval_list_HC$P
AUC_WHOLEval_F1.2_HC <- resultsval_list_HC$F1.2
AUC_WHOLEval_F2.2_HC <- resultsval_list_HC$F2.2

AUC_WHOLEval_P_HR    <- resultsval_list_HR$P
AUC_WHOLEval_F1.2_HR <- resultsval_list_HR$F1.2
AUC_WHOLEval_F2.2_HR <- resultsval_list_HR$F2.2

AUC_WHOLEval_P_FC    <- resultsval_list_FC$P
AUC_WHOLEval_F1.2_FC <- resultsval_list_FC$F1.2
AUC_WHOLEval_F2.2_FC <- resultsval_list_FC$F2.2

AUC_WHOLEval_P_FR    <- resultsval_list_FR$P
AUC_WHOLEval_F1.2_FR <- resultsval_list_FR$F1.2
AUC_WHOLEval_F2.2_FR <- resultsval_list_FR$F2.2


save(list = ls(), file = paste0("/home/thomas/Documents/Rstudio/Application/Application colrec/AUC_valid_strata.Rdata"))


#############################################################################################

##################################           PLOT

################################################################################################################

## TRAIN


#########################################################################################################################
############ Plot des valeurs de survie nette moyennes pour WHOLE
plot_df <- data.frame(
  time = c(0, newtimespred) / 365.241,
  PoharPerme = c(1, PP$surv),
  PLANN = mean.plann,
  'Spline PH' = c(1, mean.flex1),
  'Spline NPH' = c(1, mean.flex2),
  check.names = FALSE
  
)

plot_df_long <- plot_df %>%
  pivot_longer(cols = -time, names_to = "Method", values_to = "Survival")

plot_df_long$Method <- factor(
  plot_df_long$Method,
  levels = c("PoharPerme", "PLANN", "Spline PH", "Spline NPH")  # desired order
)

labels = function(x) {
  sapply(x, function(v) {
    if (v == round(v)) {
      # integer → remove .00
      format(v, trim = TRUE, scientific = FALSE, nsmall = 0)
    } else {
      # decimal → keep two digits (e.g. 0.50, 0.25, 0.75)
      format(v, trim = TRUE, scientific = FALSE, nsmall = 1)
    }
  })
}

library(scales)
p_mean_whole <- ggplot(plot_df_long, aes(x = time, y = Survival, color = Method, linetype = Method)) +
  geom_line(size = 1.2) +
  geom_line(data = subset(plot_df_long, Method == "Spline PH"),
            size = 1.2) +
  scale_linetype_manual(values=c( "PoharPerme" = "solid", "Spline PH" ="dotted", "Spline NPH" ="longdash", "PLANN" ="twodash"))+ 
  scale_color_manual(
    values = c("PoharPerme" = "black", "Spline PH" = "brown1", "Spline NPH" = "darkolivegreen", "PLANN" = "aquamarine"),
  )  +
  labs(
    x = "Time (years)",
    y = "Net survival"
  ) +
  theme_minimal(base_family = "CMU Serif",base_size = 14) +
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
        plot.title = element_text(hjust = 0.5),
        axis.ticks = element_line(color = "black", size = 0.5),  # tick line
        axis.ticks.length = unit(0.25, "cm"),
        legend.position =  "none", #c(0.98, 0.10)
        legend.justification = c(1, 0),
        legend.background = element_rect(fill = "white", color = "black")
  ) +
  scale_x_continuous(limits = c(0, 10.5), expand = expansion(mult = 0), labels = labels) +
  scale_y_continuous(limits = c(0, 1), expand = expansion(mult = 0)
                     ,
                     breaks = function(x) {
                       b <- pretty(x)
                       b[b != 0]           
                     }, labels = labels) 

print(p_mean_whole)
ggsave(
  filename = paste0("/home/thomas/Documents/Articles/Article survivalPLANN/Figures/Figures 4 strates NEW/MEAN_WHOLE_TRAIN.png"),
  plot = p_mean_whole,
  width = 6,
  height = 6,
  dpi = 800
)

# print(p_mean_whole)
# ggsave(
#   filename = paste0("/home/thomas/Documents/Articles/Article survivalPLANN/Figures/Figures 4 strates NEW/MEAN_WHOLE_TRAIN.eps"),
#   plot = p_mean_whole,
#   width = 6,
#   height = 6,
#   dpi = 800
# )
######### plots par strates
l_mean_train <- list()
l_mean_train[[1]] <- p_mean_whole
l = 2
strata_names = c("HC", "HR", "FC", "FR")
for(j in strata_names){
  if(j == "HC"){
    data_s <- data_train[data_train$sex ==1 & data_train$site == "colon",  ]
  }else if(j == "HR"){
    data_s <- data_train[data_train$sex ==1 & data_train$site == "rectum",  ]
    
  }else if(j == "FC"){
    data_s <- data_train[data_train$sex ==2 & data_train$site == "colon",  ]
  }else if(j == "FR"){
    data_s <- data_train[data_train$sex ==2 & data_train$site == "rectum",  ]
  }
  
  plannpred_s <- predictRS(plann.model, data = data_s, newtimes = newtimespred, ratetable=slopop, 
                           age = "age", year = "diag", sex = "sexchara")
  
  mean.plann_s <- plannpred_s$mpredictions$net_survival 
  
  ###flexnet
  
  ### modèle PH
  ##2 noeuds
  
  flexpred1_s <- predict(flex.model1, newtimes = newtimespred, newdata = data_s)$predictions
  
  mean.flex1_s <- apply(flexpred1_s, FUN="mean", MARGIN=2)
  
  
  #### modèle NPH
  
  ##2 noeuds
  
  flexpred2_s <- predict(flex.model2 , newtimes = newtimespred, newdata = data_s)$predictions
  
  mean.flex2_s  <- apply(flexpred2_s , FUN="mean", MARGIN=2)
  
  ## pohar perme
  
  PP_s <- summary( rs.surv(Surv(time, stat) ~ 1, data = data_s,
                           ratetable = slopop, method = "pohar-perme", add.times = newtimespred,
                           rmap = list(age = age, sex = sexchara, year = diag)), times = newtimespred, extend = TRUE) 
  
  poharpred_s <- PP_s$surv 
  
  if(j == "HC"){
    plannpred_HC <<- plannpred_s
    mean.plann_HC <<- mean.plann_s
    
    flexpred1_HC <<- flexpred1_s
    mean.flex1_HC <<- mean.flex1_s
    
    flexpred2_HC <<- flexpred2_s
    mean.flex2_HC <<- mean.flex2_s
    
    PP_HC <<- PP_s
    poharpred_HC <<- poharpred_s
  }else 
    if(j == "HR"){
      plannpred_HR <<- plannpred_s
      mean.plann_HR <<- mean.plann_s
      
      flexpred1_HR <<- flexpred1_s
      mean.flex1_HR <<- mean.flex1_s
      
      flexpred2_HR <<- flexpred2_s
      mean.flex2_HR <<- mean.flex2_s
      
      PP_HR <<- PP_s
      poharpred_HR <<- poharpred_s
    }else 
      if(j == "FC"){
        plannpred_FC <<- plannpred_s
        mean.plann_FC <<- mean.plann_s
        
        flexpred1_FC <<- flexpred1_s
        mean.flex1_FC <<- mean.flex1_s
        
        flexpred2_FC <<- flexpred2_s
        mean.flex2_FC <<- mean.flex2_s
        
        PP_FC <<- PP_s
        poharpred_FC <<- poharpred_s
      }else 
        if(j == "FR"){
          plannpred_FR <<- plannpred_s
          mean.plann_FR <<- mean.plann_s
          
          flexpred1_FR <<- flexpred1_s
          mean.flex1_FR <<- mean.flex1_s
          
          flexpred2_FR <<- flexpred2_s
          mean.flex2_FR <<- mean.flex2_s
          
          PP_FR <<- PP_s
          poharpred_FR <<- poharpred_s
        }
  
  plot_df <- data.frame(
    time = c(0, newtimespred) / 365.241,
    PoharPerme = c(1, get(paste0("PP_",j))$surv),
    PLANN = get(paste0("mean.plann_",j)),
    RP = c(1, get(paste0("mean.flex1_",j))),
    RP2 = c(1, get(paste0("mean.flex2_",j)))
  )
  
  plot_df_long <- plot_df %>%
    pivot_longer(cols = -time, names_to = "Method", values_to = "Survival")
  
  plot_df_long$Method <- factor(
    plot_df_long$Method,
    levels = c("PoharPerme", "PLANN", "RP", "RP2")  # desired order
  )
  
  p <- ggplot(plot_df_long, aes(x = time, y = Survival, color = Method, linetype = Method)) +
    geom_line(size = 1.2) +
    geom_line(data = subset(plot_df_long, Method == "RP"),
              size = 1.2) +
    scale_linetype_manual(values=c("PoharPerme" = "solid", "RP" ="dotted", "RP2" ="longdash", "PLANN" ="twodash")) + 
    scale_color_manual(
      values = c("PoharPerme" = "black", "RP" = "brown1", "RP2" = "darkolivegreen", "PLANN" = "aquamarine")
    ) +
    labs(
      x = "Time (years)",
      y = "Net survival"
    ) +
    theme_minimal(base_family = "CMU Serif", base_size = 14) +
    theme(
      legend.position = "none",        # <<---------------- remove legend
      panel.border = element_blank(),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.line = element_line(colour = "black"),
      axis.ticks = element_line(color = "black", size = 0.5),  # tick line
      axis.ticks.length = unit(0.25, "cm")
    ) +
    scale_x_continuous(limits = c(0, 10.5), expand = expansion(mult = 0), labels = labels) +
    scale_y_continuous(limits = c(0, 1), expand = expansion(mult = 0), labels = labels,
                       breaks = function(x) {
                         b <- pretty(x)
                         b[b != 0]
                       })
  
  print(p)
  l_mean_train[[l]] <- p
  ggsave(
    filename = paste0("/home/thomas/Documents/Articles/Article survivalPLANN/Figures/Figures 4 strates NEW/MEAN_",j,"_TRAIN.png"),
    plot = p,
    width = 6,
    height = 6,
    dpi = 800
  )
  l = l+1
}

################################################################
##############        PLOT AUC     #############################
################################################################

labels = function(x) {
  sapply(x, function(v) {
    if (v == round(v)) {
      # integer → remove .00
      format(v, trim = TRUE, scientific = FALSE, nsmall = 0)
    } else {
      # decimal → keep two digits (e.g. 0.50, 0.25, 0.75)
      format(v, trim = TRUE, scientific = FALSE, nsmall = 1)
    }
  })
}

main_vec <- c("3 years", "5 years", "10 years")

if(FALSE){
  ########### BOUCLE PLOTS AUC
  for(p in c(1,2,3)){
    main_name <- main_vec[p]
    time_map <- c(3,5,10)
    df <- rbind(
      data.frame(FPR = 1 - results_mod_list$P[[p]]$table$sp,TPR = results_mod_list$P[[p]]$table$se,
                 Model = "PLANN"),
      data.frame(FPR   = 1 - results_mod_list$F1.2[[p]]$table$sp,TPR = results_mod_list$F1.2[[p]]$table$se,
                 Model = "Spline PH"),
      data.frame(FPR = 1 - results_mod_list$F2.2[[p]]$table$sp, TPR = results_mod_list$F2.2[[p]]$table$se,
                 Model = "Spline NPH"))
    
    plot_to_print <- ggplot(df, aes(x = FPR, y = TPR, color = Model, linetype = Model)) +
      geom_line(size = 1.2) +
      geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
      scale_color_manual(limits = c("PLANN", "Spline PH", "Spline NPH"),
                         values = c("Spline PH" = "brown1", "Spline NPH" = "darkolivegreen", "PLANN" = "aquamarine"),labels = c(
                           PLANN = paste0("PLANN            (AUC : ", sprintf("%.4f", AUC_WHOLE_P[p]), ")"),
                           'Spline PH'    = paste0("Spline PH          (AUC : ", sprintf("%.4f", AUC_WHOLE_F1.2[p]), ")"),
                           'Spline NPH'   = paste0("Spline NPH        (AUC : ", sprintf("%.4f", AUC_WHOLE_F2.2[p]), ")"))
      )+
      scale_linetype_manual(limits = c("PLANN", "Spline PH", "Spline NPH"), values=c("Spline PH" ="dotted", "Spline NPH" ="longdash", "PLANN" ="twodash"), labels = c(
        PLANN = paste0("PLANN            (AUC : ", sprintf("%.4f", AUC_WHOLE_P[p]), ")"),
        'Spline PH'    = paste0("Spline PH          (AUC : ", sprintf("%.4f", AUC_WHOLE_F1.2[p]), ")"),
        'Spline NPH'   = paste0("Spline NPH        (AUC : ", sprintf("%.4f", AUC_WHOLE_F2.2[p]), ")"))) + 
      labs(
        x = "1 - Specificity",
        y = "Sensitivity",
        title = "ROC Curve"
      )  +
      labs(title = main_name, x = "1 - Specificity", y = "Sensitivity") +
      theme_minimal(base_family = "CMU Serif", base_size = 14) +
      theme(panel.border = element_blank(), panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
            plot.title = element_text(hjust = 0.5),
            axis.ticks = element_line(color = "black", size = 0.5),  # tick line
            axis.ticks.length = unit(0.25, "cm"),
            legend.position = c(0.98, 0.10),
            legend.justification = c(1, 0),
            legend.background = element_rect(fill = "white", color = "black")
      ) + scale_x_continuous(limits = c(-0.0025,1.05), expand = c(0,0), breaks = c(0, 0.25, 0.50, 0.75, 1), labels = labels) +
      scale_y_continuous(limits = c(0,1), expand = expansion(mult = 0), breaks = c(0.25, 0.50, 0.75, 1), labels = labels
      )
    
    print(plot_to_print)
    ggsave(
      filename = paste0("/home/thomas/Documents/Articles/Article survivalPLANN/Figures/Figures 4 strates NEW/AUC_TRAIN_WHOLE_", time_map[p],".png"),
      plot = plot_to_print,
      width = 6,
      height = 6,
      dpi = 800
    )
  }
  
  
  ######### plots par strates
  strata_names = c("HC", "HR", "FC", "FR")
  for(j in strata_names){
    
    results_mod_list_s <- get(paste0("results_mod_list_",j))
    AUC_WHOLE_P_s <- get(paste0("AUC_WHOLE_P_",j))
    AUC_WHOLE_F1.2_s <- get(paste0("AUC_WHOLE_F1.2_",j)) 
    AUC_WHOLE_F2.2_s <- get(paste0("AUC_WHOLE_F2.2_",j)) 
    
    for(p in c(1,2,3)){
      
      main_name <- main_vec[p]
      time_map <- c(3,5,10)
      
      df <- rbind(
        data.frame(FPR = 1 - results_mod_list_s$P[[p]]$table$sp,TPR = results_mod_list_s$P[[p]]$table$se,
                   Model = "PLANN"),
        data.frame(FPR   = 1 - results_mod_list_s$F1.2[[p]]$table$sp,TPR = results_mod_list_s$F1.2[[p]]$table$se,
                   Model = "Spline PH"),
        data.frame(FPR = 1 - results_mod_list_s$F2.2[[p]]$table$sp, TPR = results_mod_list_s$F2.2[[p]]$table$se,
                   Model = "Spline NPH"))
      
      plot_to_print <- ggplot(df, aes(x = FPR, y = TPR, color = Model, linetype = Model)) +
        geom_line(size = 1.2) +
        geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
        scale_color_manual(limits = c("PLANN", "Spline PH", "Spline NPH"),
                           values = c("Spline PH" = "brown1", "Spline NPH" = "darkolivegreen", "PLANN" = "aquamarine"),labels = c(
                             PLANN = paste0("PLANN            (AUC : ", sprintf("%.4f", AUC_WHOLE_P_s[p]), ")"),
                             'Spline PH'    = paste0("Spline PH          (AUC : ", sprintf("%.4f", AUC_WHOLE_F1.2_s[p]), ")"),
                             'Spline NPH'   = paste0("Spline NPH        (AUC : ", sprintf("%.4f", AUC_WHOLE_F2.2_s[p]), ")"))
        )+
        scale_linetype_manual(limits = c("PLANN", "Spline PH", "Spline NPH"), values=c("Spline PH" ="dotted", "Spline NPH" ="longdash", "PLANN" ="twodash"), labels = c(
          PLANN = paste0("PLANN            (AUC : ", sprintf("%.4f", AUC_WHOLE_P_s[p]), ")"),
          'Spline PH'    = paste0("Spline PH          (AUC : ", sprintf("%.4f", AUC_WHOLE_F1.2_s[p]), ")"),
          'Spline NPH'   = paste0("Spline NPH        (AUC : ", sprintf("%.4f", AUC_WHOLE_F2.2_s[p]), ")"))) + 
        labs(
          x = "1 - Specificity",
          y = "Sensitivity",
          title = "ROC Curve"
        )  +
        labs(title = main_name, x = "1 - Specificity", y = "Sensitivity") +
        theme_minimal(base_family = "CMU Serif", base_size = 14) +
        theme(panel.border = element_blank(), panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
              plot.title = element_text(hjust = 0.5),
              axis.ticks = element_line(color = "black", size = 0.5),  # tick line
              axis.ticks.length = unit(0.25, "cm"),
              legend.position = c(0.98, 0.10),
              legend.justification = c(1, 0),
              legend.background = element_rect(fill = "white", color = "black")
        ) + scale_x_continuous(limits = c(-0.0025,1.05), expand = c(0,0), breaks = c(0, 0.25, 0.50, 0.75, 1), labels = labels) +
        scale_y_continuous(limits = c(0,1), expand = expansion(mult = 0), breaks = c(0.25, 0.50, 0.75, 1), labels = labels
        )
      
      print(plot_to_print)
      ggsave(
        filename = paste0("/home/thomas/Documents/Articles/Article survivalPLANN/Figures/Figures 4 strates NEW/AUC_TRAIN_",j,"_", time_map[p],".png"),
        plot = plot_to_print,
        width = 6,
        height = 6,
        dpi = 800
      )
    }#fin boucle temps 
  }# fin boucle strata
}
##############################################################################################
###### plots de calibration 
l_cali_train <- list()
l = 1
for(nt in newtimes){
  
  main_name <- main_vec[match(nt, newtimes)]
  mat <- match(newtimes, newtimespred)
  correstab <- setNames(mat, newtimes)
  
  .predP <- plannpred$ipredictions$relative_survival[,-1][, as.numeric(correstab[as.character(nt)])]
  .predF1 <- flexpred1[, as.numeric(correstab[as.character(nt)])]
  .predF2 <- flexpred2[, as.numeric(correstab[as.character(nt)])]
  
  n.groups <- 5
  
  .grpsP <- as.numeric(cut(.predP, breaks = c(-Inf, quantile(.predP, seq(1/n.groups, 1, 1/n.groups))),
                           labels = 1:n.groups))
  .grpsF1 <- as.numeric(cut(.predF1, breaks = c(-Inf, quantile(.predF1, seq(1/n.groups, 1, 1/n.groups))),
                            labels = 1:n.groups))
  .grpsF2 <- as.numeric(cut(.predF2, breaks = c(-Inf, quantile(.predF2, seq(1/n.groups, 1, 1/n.groups))),
                            labels = 1:n.groups))
  
  .estP <- sapply(1:n.groups, FUN = function(x) { mean(.predP[.grpsP==x]) } )
  .estF1 <- sapply(1:n.groups, FUN = function(x) { mean(.predF1[.grpsF1==x]) } )
  .estF2 <- sapply(1:n.groups, FUN = function(x) { mean(.predF2[.grpsF2==x]) } )
  
  time <- data_train$time 
  event <- data_train$stat 
  
  .dataP <- data.frame(time = time, event = event, age = data_train$age, sexchara = data_train$sexchara, diag = data_train$diag, grps = .grpsP)
  .dataF1 <- data.frame(time = time, event = event, age = data_train$age, sexchara = data_train$sexchara, diag = data_train$diag, grps = .grpsF1)
  .dataF2 <- data.frame(time = time, event = event, age = data_train$age, sexchara = data_train$sexchara, diag = data_train$diag, grps = .grpsF2)
  
  .survfitP <- summary( rs.surv(Surv(time, event) ~ .grpsP, data = .dataP,
                                ratetable = slopop, method = "pohar-perme", add.times = nt,
                                rmap = list(age = age, sex = sexchara, year = diag)), times = nt, extend = TRUE) 
  .survfitF1 <- summary( rs.surv(Surv(time, event) ~ .grpsF1, data = .dataF1,
                                 ratetable = slopop, method = "pohar-perme", add.times = nt,
                                 rmap = list(age = age, sex = sexchara, year = diag)), times = nt, extend = TRUE)
  .survfitF2 <- summary( rs.surv(Surv(time, event) ~ .grpsF2, data = .dataF2,
                                 ratetable = slopop, method = "pohar-perme", add.times = nt,
                                 rmap = list(age = age, sex = sexchara, year = diag)), times = nt, extend = TRUE)
  
  .obsP <- .survfitP$surv
  .lowerP <- .survfitP$lower
  .upperP <- .survfitP$upper
  
  .obsF1 <- .survfitF1$surv
  .lowerF1 <-.survfitF1$lower
  .upperF1 <- .survfitF1$upper 
  
  .obsF2 <- .survfitF2$surv
  .lowerF2 <- .survfitF2$lower
  .upperF2 <- .survfitF2$upper
  
  if(hasArg(cex)==FALSE) {cex <-1} else {cex <- list(...)$cex}
  if(hasArg(cex.lab)==FALSE) {cex.lab <- 1} else {cex.lab <- list(...)$cex.lab}
  if(hasArg(cex.axis)==FALSE) {cex.axis <- 1} else {cex.axis <- list(...)$cex.axis}
  if(hasArg(cex.main)==FALSE) {cex.main <- 1} else {cex.main <- list(...)$cex.main}
  if(hasArg(type)==FALSE) {type <- "b"} else {type <- list(...)$type}
  if(hasArg(col)==FALSE) {col <- 1} else {col <- list(...)$col}
  if(hasArg(lty)==FALSE) {lty <- 1} else {lty <- list(...)$lty}
  if(hasArg(lwd)==FALSE) {lwd <- 1} else {lwd <- list(...)$lwd}
  if(hasArg(pch)==FALSE) {pch <- 16} else {pch <- list(...)$pch}
  
  if(hasArg(ylim)==FALSE) {ylim <- c(0,1)} else {ylim <- list(...)$ylim}
  if(hasArg(xlim)==FALSE) {xlim  <- c(0,1)} else {xlim <- list(...)$xlim}
  
  if(hasArg(ylab)==FALSE) {ylab <- "Pohar-Perme estimations"} else {ylab <- list(...)$ylab}
  if(hasArg(xlab)==FALSE) {xlab <- "Models estimations"} else {xlab <- list(...)$xlab}
  if(hasArg(main)==FALSE) {main <- paste0(nt/365.241, " years")} else {main <- list(...)$main}
  
  df_plot <- bind_rows(
    data.frame(Method = "PLANN", est = .estP, obs = .obsP, lower = .lowerP, upper = .upperP
    ),
    data.frame(Method = "Spline PH", est = .estF1, obs = .obsF1, lower = .lowerF1, upper = .upperF1
    ),
    data.frame(Method = "Spline NPH", est = .estF2, obs = .obsF2, lower = .lowerF2, upper = .upperF2
    )
  )
  
  df_plot$Method <- factor(
    df_plot$Method,
    levels = c("PLANN", "Spline PH", "Spline NPH")
  )
  # ---- 2. ggplot equivalent of the base R graph ----
  plot_to_print <- ggplot(df_plot, aes(x = est, y = obs, color = Method, linetype = Method)) +
    geom_line(linewidth = 1) +
    geom_point(size = 2) +
    geom_segment(aes(x = est, xend = est, y = lower, yend = upper)) +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed") +
    scale_linetype_manual(limits = c("PLANN", "Spline PH", "Spline NPH"), values=c("PLANN" ="twodash", "Spline PH" ="dotted", "Spline NPH" ="longdash" ))+ 
    scale_color_manual(limits = c("PLANN", "Spline PH", "Spline NPH"),
                       values = c( "PLANN" = "aquamarine", "Spline PH" = "brown1", "Spline NPH" = "darkolivegreen" ),
    )+
    labs(
      x = "Model estimates",
      y = "Pohar-Perme estimates",
      title = main_name
    ) +
    theme_minimal(base_family = "CMU Serif",base_size = 14) +
    theme(panel.border = element_blank(), panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
          plot.title = element_text(hjust = 0.5),
          axis.ticks = element_line(color = "black", size = 0.5),  # tick line
          axis.ticks.length = unit(0.25, "cm"),
          legend.position = "none", #c(0.98, 0.10),
          legend.justification = c(1, 0),
          legend.background = element_rect(fill = "white", color = "black")
    ) + scale_x_continuous(limits = c(-0.0025,1.05), expand = c(0,0), breaks = c(0, 0.25, 0.50, 0.75, 1), labels = labels) +
    scale_y_continuous(limits = c(0,1), expand = expansion(mult = 0), breaks = c(0.25, 0.50, 0.75, 1), labels = labels
    )
  
  time_map <- c(3,5,10)
  print(plot_to_print)
  l_cali_train[[l]] <- plot_to_print
  ggsave(
    filename = paste0("/home/thomas/Documents/Articles/Article survivalPLANN/Figures/Figures 4 strates NEW/Calibration_TRAIN_WHOLE_",  time_map[which(time_map == nt/365.241)]
                      ,".png"),
    plot = plot_to_print,
    width = 6,
    height = 6,
    dpi = 800
  )
  l = l+1
}
############################

#####Par strates 
l_cali_train_strata <- list()
l = 1
for(j in strata_names){
  
  if(j == "HC"){
    data_s <- data_train[data_train$sex ==1 & data_train$site == "colon",  ]
    plannpred_s <- plannpred_HC
    flexpred1_s <- flexpred1_HC
    flexpred2_s <- flexpred2_HC
  }else
    if(j == "HR"){
      data_s <- data_train[data_train$sex ==1 & data_train$site == "rectum",  ]
      plannpred_s <- plannpred_HR
      flexpred1_s <- flexpred1_HR
      flexpred2_s <- flexpred2_HR
    }else 
      if(j == "FC"){
        data_s <- data_train[data_train$sex ==2 & data_train$site == "colon",  ]
        plannpred_s <- plannpred_FC
        flexpred1_s <- flexpred1_FC
        flexpred2_s <- flexpred2_FC
      }else 
        if(j == "FR"){
          data_s <- data_train[data_train$sex ==2 & data_train$site == "rectum",  ]
          plannpred_s <- plannpred_FR
          flexpred1_s <- flexpred1_FR
          flexpred2_s <- flexpred2_FR  
        }
  
  
  for(nt in newtimes){
    
    main_name <- main_vec[match(nt, newtimes)]
    
    mat <- match(newtimes, newtimespred)
    correstab <- setNames(mat, newtimes)
    
    .predP <- plannpred_s$ipredictions$relative_survival[,-1][, as.numeric(correstab[as.character(nt)])]
    .predF1 <- flexpred1_s[, as.numeric(correstab[as.character(nt)])]
    .predF2 <- flexpred2_s[, as.numeric(correstab[as.character(nt)])]
    
    n.groups <- 5
    
    .grpsP <- as.numeric(cut(.predP, breaks = c(-Inf, quantile(.predP, seq(1/n.groups, 1, 1/n.groups))),
                             labels = 1:n.groups))
    .grpsF1 <- as.numeric(cut(.predF1, breaks = c(-Inf, quantile(.predF1, seq(1/n.groups, 1, 1/n.groups))),
                              labels = 1:n.groups))
    .grpsF2 <- as.numeric(cut(.predF2, breaks = c(-Inf, quantile(.predF2, seq(1/n.groups, 1, 1/n.groups))),
                              labels = 1:n.groups))
    
    .estP <- sapply(1:n.groups, FUN = function(x) { mean(.predP[.grpsP==x]) } )
    .estF1 <- sapply(1:n.groups, FUN = function(x) { mean(.predF1[.grpsF1==x]) } )
    .estF2 <- sapply(1:n.groups, FUN = function(x) { mean(.predF2[.grpsF2==x]) } )
    
    time <- data_s$time 
    event <- data_s$stat 
    
    .dataP <- data.frame(time = time, event = event, age = data_s$age, sexchara = data_s$sexchara, diag = data_s$diag, grps = .grpsP)
    .dataF1 <- data.frame(time = time, event = event, age = data_s$age, sexchara = data_s$sexchara, diag = data_s$diag, grps = .grpsF1)
    .dataF2 <- data.frame(time = time, event = event, age = data_s$age, sexchara = data_s$sexchara, diag = data_s$diag, grps = .grpsF2)
    
    .survfitP <- summary( rs.surv(Surv(time, event) ~ .grpsP, data = .dataP,
                                  ratetable = slopop, method = "pohar-perme", add.times = nt,
                                  rmap = list(age = age, sex = sexchara, year = diag)), times = nt, extend = TRUE) 
    .survfitF1 <- summary( rs.surv(Surv(time, event) ~ .grpsF1, data = .dataF1,
                                   ratetable = slopop, method = "pohar-perme", add.times = nt,
                                   rmap = list(age = age, sex = sexchara, year = diag)), times = nt, extend = TRUE)
    .survfitF2 <- summary( rs.surv(Surv(time, event) ~ .grpsF2, data = .dataF2,
                                   ratetable = slopop, method = "pohar-perme", add.times = nt,
                                   rmap = list(age = age, sex = sexchara, year = diag)), times = nt, extend = TRUE)
    
    .obsP <- .survfitP$surv
    .lowerP <- .survfitP$lower
    .upperP <- .survfitP$upper
    
    .obsF1 <- .survfitF1$surv
    .lowerF1 <-.survfitF1$lower
    .upperF1 <- .survfitF1$upper 
    
    .obsF2 <- .survfitF2$surv
    .lowerF2 <- .survfitF2$lower
    .upperF2 <- .survfitF2$upper
    
    if(hasArg(cex)==FALSE) {cex <-1} else {cex <- list(...)$cex}
    if(hasArg(cex.lab)==FALSE) {cex.lab <- 1} else {cex.lab <- list(...)$cex.lab}
    if(hasArg(cex.axis)==FALSE) {cex.axis <- 1} else {cex.axis <- list(...)$cex.axis}
    if(hasArg(cex.main)==FALSE) {cex.main <- 1} else {cex.main <- list(...)$cex.main}
    if(hasArg(type)==FALSE) {type <- "b"} else {type <- list(...)$type}
    if(hasArg(col)==FALSE) {col <- 1} else {col <- list(...)$col}
    if(hasArg(lty)==FALSE) {lty <- 1} else {lty <- list(...)$lty}
    if(hasArg(lwd)==FALSE) {lwd <- 1} else {lwd <- list(...)$lwd}
    if(hasArg(pch)==FALSE) {pch <- 16} else {pch <- list(...)$pch}
    
    if(hasArg(ylim)==FALSE) {ylim <- c(0,1)} else {ylim <- list(...)$ylim}
    if(hasArg(xlim)==FALSE) {xlim  <- c(0,1)} else {xlim <- list(...)$xlim}
    
    if(hasArg(ylab)==FALSE) {ylab <- "Pohar-Perme estimations"} else {ylab <- list(...)$ylab}
    if(hasArg(xlab)==FALSE) {xlab <- "Models estimations"} else {xlab <- list(...)$xlab}
    if(hasArg(main)==FALSE) {main <- paste0(nt/365.241, " years")} else {main <- list(...)$main}
    
    df_plot <- bind_rows(
      data.frame(Method = "PLANN", est = .estP, obs = .obsP, lower = .lowerP, upper = .upperP
      ),
      data.frame(Method = "RP", est = .estF1, obs = .obsF1, lower = .lowerF1, upper = .upperF1
      ),
      data.frame(Method = "RP2", est = .estF2, obs = .obsF2, lower = .lowerF2, upper = .upperF2
      )
    )
    maxy <- ifelse(max(df_plot$upper) > 1,
                   max(df_plot$upper) + 0.01,
                   1) 
    # ---- 2. ggplot equivalent of the base R graph ----
    p <- ggplot(df_plot, aes(x = est, y = obs, color = Method, linetype = Method)) +
      geom_line(linewidth = 1) +
      geom_point(size = 2) +
      geom_segment(aes(x = est, xend = est, y = lower, yend = upper)) +
      geom_abline(intercept = 0, slope = 1, linetype = "dashed") +
      scale_color_manual(
        values = c(
          "PLANN" = "aquamarine",
          "RP" = "brown1",
          "RP2" = "darkolivegreen"
        )
      ) +
      labs(
        x = "Model estimates",
        y = "Pohar-Perme estimates",
        title = main_name
      ) +
      scale_linetype_manual(values=c("RP" ="dotted", "RP2" ="longdash", "PLANN" ="twodash")) + 
      theme_bw() +
      theme_minimal(base_family = "CMU Serif",base_size = 14) +
      theme(legend.position = "none",
            panel.border = element_blank(),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            plot.title = element_text(hjust = 0.5),
            axis.ticks = element_line(color = "black", size = 0.5),  # tick line
            axis.ticks.length = unit(0.25, "cm"),
            axis.line = element_line(colour = "black")
      ) + scale_x_continuous(limits = c(-0.0025,1.05), expand = c(0,0), breaks = c(0, 0.25, 0.50, 0.75, 1), labels = labels) +
      scale_y_continuous( expand = expansion(mult = 0), breaks = c(0.25, 0.50, 0.75, 1), labels = labels
      ) +
      coord_cartesian(ylim = c(0, maxy))
    
    print(p)
    l_cali_train_strata[[l]] <- p
    time_map <- c(3,5,10)
    ggsave(
      filename = paste0("/home/thomas/Documents/Articles/Article survivalPLANN/Figures/Figures 4 strates NEW/Calibration_TRAIN_", j,"_",  time_map[which(time_map == nt/365.241)]
                        ,".png"),
      plot = p,
      width = 6,
      height = 6,
      dpi = 800
    )
    l = l+1
  }
  
}


########################### VALIDATION 
############### Plot des valeurs moyennes base entière validation

plot_df <- data.frame(
  time = c(0, newtimespred) / 365.241,
  PoharPerme = c(1, PPval$surv),
  PLANN = mean.plannval,
  RP = c(1, mean.flexval1),
  RP2 = c(1, mean.flexval2)
)

plot_df_long <- plot_df %>%
  pivot_longer(cols = -time, names_to = "Method", values_to = "Survival")

plot_df_long$Method <- factor(
  plot_df_long$Method,
  levels = c("PoharPerme", "PLANN", "RP", "RP2")  # desired order
)

labels = function(x) {
  sapply(x, function(v) {
    if (v == round(v)) {
      # integer → remove .00
      format(v, trim = TRUE, scientific = FALSE, nsmall = 0)
    } else {
      # decimal → keep two digits (e.g. 0.50, 0.25, 0.75)
      format(v, trim = TRUE, scientific = FALSE, nsmall = 1)
    }
  })
}

library(scales)
p_mean_whole_val <-ggplot(plot_df_long, aes(x = time, y = Survival, color = Method, linetype = Method)) +
  geom_line(size = 1.2) +
  geom_line(data = subset(plot_df_long, Method == "RP"),
            size = 1.2) +
  scale_linetype_manual(values=c( "PoharPerme" = "solid", "RP" ="dotted", "RP2" ="longdash", "PLANN" ="twodash"))+ 
  scale_color_manual(
    values = c("PoharPerme" = "black", "RP" = "brown1", "RP2" = "darkolivegreen", "PLANN" = "aquamarine"),
  )  +
  labs(
    x = "Time (years)",
    y = "Net survival"
  ) +
  theme_minimal(base_family = "CMU Serif",base_size = 14) +
  theme(
    legend.position = "none",
    legend.background = element_rect(
      fill = "white",                # background color of the box
      color = "black",               # border color of the box
      size = 0.5,                    # border thickness
      linetype = "solid"
    ),
    panel.border = element_blank(), panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
    axis.ticks = element_line(color = "black", size = 0.5),  # tick line
    axis.ticks.length = unit(0.25, "cm")) +
  scale_x_continuous(limits = c(0, 10.5), expand = expansion(mult = 0), labels = labels) +
  scale_y_continuous(limits = c(0, 1), expand = expansion(mult = 0)
                     ,
                     breaks = function(x) {
                       b <- pretty(x)
                       b[b != 0]           
                     }, labels = labels) 

print(p_mean_whole_val)
ggsave(
  filename = paste0("/home/thomas/Documents/Articles/Article survivalPLANN/Figures/Figures 4 strates NEW/MEAN_WHOLE_VALID.png"),
  plot = p_mean_whole_val,
  width = 6,
  height = 6,
  dpi = 800
)
######### plots par strates
l_mean_valid <- list()
l_mean_valid[[1]] <- p_mean_whole_val
l = 2
strata_names = c("HC", "HR", "FC", "FR")
for(j in strata_names){
  if(j == "HC"){
    data_s <- data_valid[data_valid$sex ==1 & data_valid$site == "colon",  ]
  }else if(j == "HR"){
    data_s <- data_valid[data_valid$sex ==1 & data_valid$site == "rectum",  ]
    
  }else if(j == "FC"){
    data_s <- data_valid[data_valid$sex ==2 & data_valid$site == "colon",  ]
  }else if(j == "FR"){
    data_s <- data_valid[data_valid$sex ==2 & data_valid$site == "rectum",  ]
  }
  
  plannpred_s <- predictRS(plann.model, data = data_s, newtimes = newtimespred, ratetable=slopop, 
                           age = "age", year = "diag", sex = "sexchara")
  
  mean.plann_s <- plannpred_s$mpredictions$net_survival 
  
  ###flexnet
  
  ### modèle PH
  ##2 noeuds
  
  flexpred1_s <- predict(flex.model1, newtimes = newtimespred, newdata = data_s)$predictions
  
  mean.flex1_s <- apply(flexpred1_s, FUN="mean", MARGIN=2)
  
  
  #### modèle NPH
  
  ##2 noeuds
  
  flexpred2_s <- predict(flex.model2 , newtimes = newtimespred, newdata = data_s)$predictions
  
  mean.flex2_s  <- apply(flexpred2_s , FUN="mean", MARGIN=2)
  
  ## pohar perme
  
  PP_s <- summary( rs.surv(Surv(time, stat) ~ 1, data = data_s,
                           ratetable = slopop, method = "pohar-perme", add.times = newtimespred,
                           rmap = list(age = age, sex = sexchara, year = diag)), times = newtimespred, extend = TRUE) 
  
  poharpred_s <- PP_s$surv 
  
  if(j == "HC"){
    plannpredval_HC <<- plannpred_s
    mean.plannval_HC <<- mean.plann_s
    
    flexpred1val_HC <<- flexpred1_s
    mean.flex1val_HC <<- mean.flex1_s
    
    flexpred2val_HC <<- flexpred2_s
    mean.flex2val_HC <<- mean.flex2_s
    
    PPval_HC <<- PP_s
    poharpredval_HC <<- poharpred_s
  }else if(j == "HR"){
    plannpredval_HR <<- plannpred_s
    mean.plannval_HR <<- mean.plann_s
    
    flexpred1val_HR <<- flexpred1_s
    mean.flex1val_HR <<- mean.flex1_s
    
    flexpred2val_HR <<- flexpred2_s
    mean.flex2val_HR <<- mean.flex2_s
    
    PPval_HR <<- PP_s
    poharpredval_HR <<- poharpred_s
  }else if(j == "FC"){
    plannpredval_FC <<- plannpred_s
    mean.plannval_FC <<- mean.plann_s
    
    flexpred1val_FC <<- flexpred1_s
    mean.flex1val_FC <<- mean.flex1_s
    
    flexpred2val_FC <<- flexpred2_s
    mean.flex2val_FC <<- mean.flex2_s
    
    PPval_FC <<- PP_s
    poharpredval_FC <<- poharpred_s
  }else if(j == "FR"){
    plannpredval_FR <<- plannpred_s
    mean.plannval_FR <<- mean.plann_s
    
    flexpred1val_FR <<- flexpred1_s
    mean.flex1val_FR <<- mean.flex1_s
    
    flexpred2val_FR <<- flexpred2_s
    mean.flex2val_FR <<- mean.flex2_s
    
    PPval_FR <<- PP_s
    poharpredval_FR <<- poharpred_s
  }
  
  plot_df <- data.frame(
    time = c(0, newtimespred) / 365.241,
    PoharPerme = c(1, get(paste0("PP_",j))$surv),
    PLANN = get(paste0("mean.plann_",j)),
    RP = c(1, get(paste0("mean.flex1_",j))),
    RP2 = c(1, get(paste0("mean.flex2_",j)))
  )
  
  plot_df_long <- plot_df %>%
    pivot_longer(cols = -time, names_to = "Method", values_to = "Survival")
  
  plot_df_long$Method <- factor(
    plot_df_long$Method,
    levels = c("PoharPerme", "PLANN", "RP", "RP2")  # desired order
  )
  
  p <- ggplot(plot_df_long, aes(x = time, y = Survival, color = Method, linetype = Method)) +
    geom_line(size = 1.2) +
    geom_line(data = subset(plot_df_long, Method == "RP"),
              size = 1.2) +
    scale_linetype_manual(values=c("PoharPerme" = "solid", "RP" ="dotted", "RP2" ="longdash", "PLANN" ="twodash")) + 
    scale_color_manual(
      values = c("PoharPerme" = "black", "RP" = "brown1", "RP2" = "darkolivegreen", "PLANN" = "aquamarine")
    ) +
    labs(
      x = "Time (years)",
      y = "Net survival"
    ) +
    theme_minimal(base_family = "CMU Serif", base_size = 14) +
    theme(
      legend.position = "none",        # <<---------------- remove legend
      panel.border = element_blank(),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.line = element_line(colour = "black"),
      axis.ticks = element_line(color = "black", size = 0.5),  # tick line
      axis.ticks.length = unit(0.25, "cm")
    ) +
    scale_x_continuous(limits = c(0, 10.5), expand = expansion(mult = 0), labels = labels) +
    scale_y_continuous(limits = c(0, 1), expand = expansion(mult = 0), labels = labels,
                       breaks = function(x) {
                         b <- pretty(x)
                         b[b != 0]
                       })
  
  print(p)
  l_mean_valid[[l]] <- p
  ggsave(
    filename = paste0("/home/thomas/Documents/Articles/Article survivalPLANN/Figures/Figures 4 strates NEW/MEAN_",j,"_VALID.png"),
    plot = p,
    width = 6,
    height = 6,
    dpi = 800
  )
  l = l+1
}

##############"PLOT AUC

labels = function(x) {
  sapply(x, function(v) {
    if (v == round(v)) {
      # integer → remove .00
      format(v, trim = TRUE, scientific = FALSE, nsmall = 0)
    } else {
      # decimal → keep two digits (e.g. 0.50, 0.25, 0.75)
      format(v, trim = TRUE, scientific = FALSE, nsmall = 1)
    }
  })
}

main_vec <- c("3 years", "5 years", "10 years")

############## PLOTS AUC VAL
if(FALSE){
  for(p in c(1,2,3)){
    main_name <- main_vec[p]
    time_map <- c(3,5,10)
    
    df <- rbind(
      data.frame(FPR = 1 - resultsval_mod_list$P[[p]]$table$sp,TPR = resultsval_mod_list$P[[p]]$table$se,
                 Model = "PLANN"),
      data.frame(FPR   = 1 - resultsval_mod_list$F1.2[[p]]$table$sp,TPR = resultsval_mod_list$F1.2[[p]]$table$se,
                 Model = "Spline PH"),
      data.frame(FPR = 1 - resultsval_mod_list$F2.2[[p]]$table$sp, TPR = resultsval_mod_list$F2.2[[p]]$table$se,
                 Model = "Spline NPH"))
    
    plot_to_print <- ggplot(df, aes(x = FPR, y = TPR, color = Model, linetype = Model)) +
      geom_line(size = 1.2) +
      geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
      scale_color_manual(limits = c("PLANN", "Spline PH", "Spline NPH"),
                         values = c("Spline PH" = "brown1", "Spline NPH" = "darkolivegreen", "PLANN" = "aquamarine"),labels = c(
                           PLANN = paste0("PLANN            (AUC : ", sprintf("%.4f", AUC_WHOLEval_P[p]), ")"),
                           'Spline PH'    = paste0("Spline PH          (AUC : ", sprintf("%.4f", AUC_WHOLEval_F1.2[p]), ")"),
                           'Spline NPH'   = paste0("Spline NPH        (AUC : ", sprintf("%.4f", AUC_WHOLEval_F2.2[p]), ")"))
      )+
      scale_linetype_manual(limits = c("PLANN", "Spline PH", "Spline NPH"), values=c("Spline PH" ="dotted", "Spline NPH" ="longdash", "PLANN" ="twodash"), labels = c(
        PLANN = paste0("PLANN            (AUC : ", sprintf("%.4f", AUC_WHOLEval_P[p]), ")"),
        'Spline PH'    = paste0("Spline PH          (AUC : ", sprintf("%.4f", AUC_WHOLEval_F1.2[p]), ")"),
        'Spline NPH'   = paste0("Spline NPH        (AUC : ", sprintf("%.4f", AUC_WHOLEval_F2.2[p]), ")"))) + 
      labs(
        x = "1 - Specificity",
        y = "Sensitivity",
        title = "ROC Curve"
      )  +
      labs(title = main_name, x = "1 - Specificity", y = "Sensitivity") +
      theme_minimal(base_family = "CMU Serif", base_size = 14) +
      theme(panel.border = element_blank(), panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
            plot.title = element_text(hjust = 0.5),
            axis.ticks = element_line(color = "black", size = 0.5),  # tick line
            axis.ticks.length = unit(0.25, "cm"),
            legend.position = c(0.98, 0.10),
            legend.justification = c(1, 0),
            legend.background = element_rect(fill = "white", color = "black")
      ) + scale_x_continuous(limits = c(-0.0025,1.05), expand = c(0,0), breaks = c(0, 0.25, 0.50, 0.75, 1), labels = labels) +
      scale_y_continuous(limits = c(0,1), expand = expansion(mult = 0), breaks = c(0.25, 0.50, 0.75, 1), labels = labels
      )
    
    print(plot_to_print)
    ggsave(
      filename = paste0("/home/thomas/Documents/Articles/Article survivalPLANN/Figures/Figures 4 strates/AUC_VALID_WHOLE_", time_map[p],".png"),
      plot = plot_to_print,
      width = 6,
      height = 6,
      dpi = 800
    )
  }
  
  ######### plots par strates
  strata_names = c("HC", "HR", "FC", "FR")
  for(j in strata_names){
    
    results_mod_list_s <- get(paste0("resultsval_mod_list_",j))
    AUC_WHOLE_P_s <- get(paste0("AUC_WHOLEval_P_",j))
    AUC_WHOLE_F1.2_s <- get(paste0("AUC_WHOLEval_F1.2_",j)) 
    AUC_WHOLE_F2.2_s <- get(paste0("AUC_WHOLEval_F2.2_",j)) 
    
    for(p in c(1,2,3)){
      
      main_name <- main_vec[p]
      time_map <- c(3,5,10)
      
      df <- rbind(
        data.frame(FPR = 1 - results_mod_list_s$P[[p]]$table$sp,TPR = results_mod_list_s$P[[p]]$table$se,
                   Model = "PLANN"),
        data.frame(FPR   = 1 - results_mod_list_s$F1.2[[p]]$table$sp,TPR = results_mod_list_s$F1.2[[p]]$table$se,
                   Model = "Spline PH"),
        data.frame(FPR = 1 - results_mod_list_s$F2.2[[p]]$table$sp, TPR = results_mod_list_s$F2.2[[p]]$table$se,
                   Model = "Spline NPH"))
      
      plot_to_print <- ggplot(df, aes(x = FPR, y = TPR, color = Model, linetype = Model)) +
        geom_line(size = 1.2) +
        geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
        scale_color_manual(limits = c("PLANN", "Spline PH", "Spline NPH"),
                           values = c("Spline PH" = "brown1", "Spline NPH" = "darkolivegreen", "PLANN" = "aquamarine"),labels = c(
                             PLANN = paste0("PLANN            (AUC : ", sprintf("%.4f", AUC_WHOLE_P_s[p]), ")"),
                             'Spline PH'    = paste0("Spline PH          (AUC : ", sprintf("%.4f", AUC_WHOLE_F1.2_s[p]), ")"),
                             'Spline NPH'   = paste0("Spline NPH        (AUC : ", sprintf("%.4f", AUC_WHOLE_F2.2_s[p]), ")"))
        )+
        scale_linetype_manual(limits = c("PLANN", "Spline PH", "Spline NPH"), values=c("Spline PH" ="dotted", "Spline NPH" ="longdash", "PLANN" ="twodash"), labels = c(
          PLANN = paste0("PLANN            (AUC : ", sprintf("%.4f", AUC_WHOLE_P_s[p]), ")"),
          'Spline PH'    = paste0("Spline PH          (AUC : ", sprintf("%.4f", AUC_WHOLE_F1.2_s[p]), ")"),
          'Spline NPH'   = paste0("Spline NPH        (AUC : ", sprintf("%.4f", AUC_WHOLE_F2.2_s[p]), ")"))) + 
        labs(
          x = "1 - Specificity",
          y = "Sensitivity",
          title = "ROC Curve"
        )  +
        labs(title = main_name, x = "1 - Specificity", y = "Sensitivity") +
        theme_minimal(base_family = "CMU Serif", base_size = 14) +
        theme(panel.border = element_blank(), panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
              plot.title = element_text(hjust = 0.5),
              axis.ticks = element_line(color = "black", size = 0.5),  # tick line
              axis.ticks.length = unit(0.25, "cm"),
              legend.position = c(0.98, 0.10),
              legend.justification = c(1, 0),
              legend.background = element_rect(fill = "white", color = "black")
        ) + scale_x_continuous(limits = c(-0.0025,1.05), expand = c(0,0), breaks = c(0, 0.25, 0.50, 0.75, 1), labels = labels) +
        scale_y_continuous(limits = c(0,1), expand = expansion(mult = 0), breaks = c(0.25, 0.50, 0.75, 1), labels = labels
        )
      
      print(plot_to_print)
      ggsave(
        filename = paste0("/home/thomas/Documents/Articles/Article survivalPLANN/Figures/Figures 4 strates NEW/AUC_VALID_",j,"_", time_map[p],".png"),
        plot = plot_to_print,
        width = 6,
        height = 6,
        dpi = 800
      )
      
    }#fin boucle temps 
  }# fin boucle strata
}
#####################################################################################

###### plots de calibration 
l_cali_valid <- list()
l = 1
for(nt in newtimes){
  
  main_name <- main_vec[match(nt, newtimes)]
  mat <- match(newtimes, newtimespred)
  correstab <- setNames(mat, newtimes)
  
  .predP <- plannpredval$ipredictions$relative_survival[,-1][, as.numeric(correstab[as.character(nt)])]
  .predF1 <- flexpredval1[, as.numeric(correstab[as.character(nt)])]
  .predF2 <- flexpredval2[, as.numeric(correstab[as.character(nt)])]
  
  n.groups <- 5
  
  .grpsP <- as.numeric(cut(.predP,
                           breaks = c(-Inf, quantile(.predP, seq(1/n.groups, 1, 1/n.groups))),
                           labels = 1:n.groups))
  .grpsF1 <- as.numeric(cut(.predF1,
                            breaks = c(-Inf, quantile(.predF1, seq(1/n.groups, 1, 1/n.groups))),
                            labels = 1:n.groups))
  .grpsF2 <- as.numeric(cut(.predF2,
                            breaks = c(-Inf, quantile(.predF2, seq(1/n.groups, 1, 1/n.groups))),
                            labels = 1:n.groups))
  
  .estP <- sapply(1:n.groups, FUN = function(x) { mean(.predP[.grpsP==x]) } )
  .estF1 <- sapply(1:n.groups, FUN = function(x) { mean(.predF1[.grpsF1==x]) } )
  .estF2 <- sapply(1:n.groups, FUN = function(x) { mean(.predF2[.grpsF2==x]) } )
  
  time <- data_valid$time 
  event <- data_valid$stat 
  
  .dataP <- data.frame(time = time, event = event, age = data_valid$age, sexchara = data_valid$sexchara, diag = data_valid$diag, grps = .grpsP)
  .dataF1 <- data.frame(time = time, event = event, age = data_valid$age, sexchara = data_valid$sexchara, diag = data_valid$diag, grps = .grpsF1)
  .dataF2 <- data.frame(time = time, event = event, age = data_valid$age, sexchara = data_valid$sexchara, diag = data_valid$diag, grps = .grpsF2)
  
  .survfitP <- summary( rs.surv(Surv(time, event) ~ .grpsP, data = .dataP,
                                ratetable = slopop, method = "pohar-perme", add.times = nt,
                                rmap = list(age = age, sex = sexchara, year = diag)), times = nt, extend = TRUE) 
  .survfitF1 <- summary( rs.surv(Surv(time, event) ~ .grpsF1, data = .dataF1,
                                 ratetable = slopop, method = "pohar-perme", add.times = nt,
                                 rmap = list(age = age, sex = sexchara, year = diag)), times = nt, extend = TRUE)
  .survfitF2 <- summary( rs.surv(Surv(time, event) ~ .grpsF2, data = .dataF2,
                                 ratetable = slopop, method = "pohar-perme", add.times = nt,
                                 rmap = list(age = age, sex = sexchara, year = diag)), times = nt, extend = TRUE)
  
  .obsP <- .survfitP$surv
  
  .lowerP <- .survfitP$lower
  
  .upperP <- .survfitP$upper
  
  .obsF1 <- .survfitF1$surv
  
  .lowerF1 <-.survfitF1$lower
  
  .upperF1 <- .survfitF1$upper 
  
  .obsF2 <- .survfitF2$surv
  
  .lowerF2 <- .survfitF2$lower
  
  .upperF2 <- .survfitF2$upper
  
  if(hasArg(cex)==FALSE) {cex <-1} else {cex <- list(...)$cex}
  if(hasArg(cex.lab)==FALSE) {cex.lab <- 1} else {cex.lab <- list(...)$cex.lab}
  if(hasArg(cex.axis)==FALSE) {cex.axis <- 1} else {cex.axis <- list(...)$cex.axis}
  if(hasArg(cex.main)==FALSE) {cex.main <- 1} else {cex.main <- list(...)$cex.main}
  if(hasArg(type)==FALSE) {type <- "b"} else {type <- list(...)$type}
  if(hasArg(col)==FALSE) {col <- 1} else {col <- list(...)$col}
  if(hasArg(lty)==FALSE) {lty <- 1} else {lty <- list(...)$lty}
  if(hasArg(lwd)==FALSE) {lwd <- 1} else {lwd <- list(...)$lwd}
  if(hasArg(pch)==FALSE) {pch <- 16} else {pch <- list(...)$pch}
  
  if(hasArg(ylim)==FALSE) {ylim <- c(0,1)} else {ylim <- list(...)$ylim}
  if(hasArg(xlim)==FALSE) {xlim  <- c(0,1)} else {xlim <- list(...)$xlim}
  
  if(hasArg(ylab)==FALSE) {ylab <- "Pohar-Perme estimations"} else {ylab <- list(...)$ylab}
  if(hasArg(xlab)==FALSE) {xlab <- "Models estimations"} else {xlab <- list(...)$xlab}
  if(hasArg(main)==FALSE) {main <- paste0(nt/365.241, " years")} else {main <- list(...)$main}
  
  df_plot <- bind_rows(
    data.frame(
      Method = "PLANN",
      est = .estP,
      obs = .obsP,
      lower = .lowerP,
      upper = .upperP
    ),
    data.frame(
      Method = "RP",
      est = .estF1,
      obs = .obsF1,
      lower = .lowerF1,
      upper = .upperF1
    ),
    data.frame(
      Method = "RP2",
      est = .estF2,
      obs = .obsF2,
      lower = .lowerF2,
      upper = .upperF2
    )
  )
  
  # ---- 2. ggplot equivalent of the base R graph ----
  plot_to_print <- ggplot(df_plot, aes(x = est, y = obs, color = Method, linetype = Method)) +
    geom_line(linewidth = 1) +
    geom_point(size = 2) +
    geom_segment(aes(x = est, xend = est, y = lower, yend = upper)) +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed") +
    scale_color_manual(
      values = c(
        "PLANN" = "aquamarine",
        "RP" = "brown1",
        "RP2" = "darkolivegreen"
      )
    ) +
    labs(
      x = "Model estimates",
      y = "Pohar-Perme estimates",
      title = main_name
    ) +
    scale_linetype_manual(values=c("RP" ="dotted", "RP2" ="longdash", "PLANN" ="twodash")) + 
    theme_minimal(base_family = "CMU Serif",base_size = 14) +
    theme(legend.position = "none",
          panel.border = element_blank(), panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
          plot.title = element_text(hjust = 0.5),
          axis.ticks = element_line(color = "black", size = 0.5),  # tick line
          axis.ticks.length = unit(0.25, "cm")) + scale_x_continuous(limits = c(-0.0025,1.05), expand = c(0,0), breaks = c(0, 0.25, 0.50, 0.75, 1), labels = labels) +
    scale_y_continuous(limits = c(0,1), expand = expansion(mult = 0), breaks = c(0.25, 0.50, 0.75, 1), labels = labels
    )
  assign(paste0("cali_WHOLE_val_",nt), p)
  print(p)
  
  time_map <- c(3,5,10)
  print(plot_to_print)
  l_cali_valid[[l]] <- plot_to_print
  ggsave(
    filename = paste0("/home/thomas/Documents/Articles/Article survivalPLANN/Figures/Figures 4 strates NEW/Calibration_VALID_WHOLE_",  time_map[which(time_map == nt/365.241)]
                      ,".png"),
    plot = plot_to_print,
    width = 6,
    height = 6,
    dpi = 800
  )
  l = l+1
  
}
############################

#####Par strates 
l_cali_valid_strata <- list()
l = 1
for(j in strata_names){
  
  if(j == "HC"){
    data_s <- data_valid[data_valid$sex ==1 & data_valid$site == "colon",  ]
    plannpred_s <- plannpredval_HC
    flexpred1_s <- flexpred1val_HC
    flexpred2_s <- flexpred2val_HC
  }else if(j == "HR"){
    data_s <- data_valid[data_valid$sex ==1 & data_valid$site == "rectum",  ]
    plannpred_s <- plannpredval_HR
    flexpred1_s <- flexpred1val_HR
    flexpred2_s <- flexpred2val_HR
  }else if(j == "FC"){
    data_s <- data_valid[data_valid$sex ==2 & data_valid$site == "colon",  ]
    plannpred_s <- plannpredval_FC
    flexpred1_s <- flexpred1val_FC
    flexpred2_s <- flexpred2val_FC
  }else if(j == "FR"){
    data_s <- data_valid[data_valid$sex ==2 & data_valid$site == "rectum",  ]
    plannpred_s <- plannpredval_FR
    flexpred1_s <- flexpred1val_FR
    flexpred2_s <- flexpred2val_FR  
  }
  
  
  for(nt in newtimes){
    
    main_name <- main_vec[match(nt, newtimes)]
    
    mat <- match(newtimes, newtimespred)
    correstab <- setNames(mat, newtimes)
    
    .predP <- plannpred_s$ipredictions$relative_survival[,-1][, as.numeric(correstab[as.character(nt)])]
    .predF1 <- flexpred1_s[, as.numeric(correstab[as.character(nt)])]
    .predF2 <- flexpred2_s[, as.numeric(correstab[as.character(nt)])]
    
    n.groups <- 5
    
    .grpsP <- as.numeric(cut(.predP, breaks = c(-Inf, quantile(.predP, seq(1/n.groups, 1, 1/n.groups))),
                             labels = 1:n.groups))
    .grpsF1 <- as.numeric(cut(.predF1, breaks = c(-Inf, quantile(.predF1, seq(1/n.groups, 1, 1/n.groups))),
                              labels = 1:n.groups))
    .grpsF2 <- as.numeric(cut(.predF2, breaks = c(-Inf, quantile(.predF2, seq(1/n.groups, 1, 1/n.groups))),
                              labels = 1:n.groups))
    
    .estP <- sapply(1:n.groups, FUN = function(x) { mean(.predP[.grpsP==x]) } )
    .estF1 <- sapply(1:n.groups, FUN = function(x) { mean(.predF1[.grpsF1==x]) } )
    .estF2 <- sapply(1:n.groups, FUN = function(x) { mean(.predF2[.grpsF2==x]) } )
    
    time <- data_s$time 
    event <- data_s$stat 
    
    .dataP <- data.frame(time = time, event = event, age = data_s$age, sexchara = data_s$sexchara, diag = data_s$diag, grps = .grpsP)
    .dataF1 <- data.frame(time = time, event = event, age = data_s$age, sexchara = data_s$sexchara, diag = data_s$diag, grps = .grpsF1)
    .dataF2 <- data.frame(time = time, event = event, age = data_s$age, sexchara = data_s$sexchara, diag = data_s$diag, grps = .grpsF2)
    
    .survfitP <- summary( rs.surv(Surv(time, event) ~ .grpsP, data = .dataP,
                                  ratetable = slopop, method = "pohar-perme", add.times = nt,
                                  rmap = list(age = age, sex = sexchara, year = diag)), times = nt, extend = TRUE) 
    .survfitF1 <- summary( rs.surv(Surv(time, event) ~ .grpsF1, data = .dataF1,
                                   ratetable = slopop, method = "pohar-perme", add.times = nt,
                                   rmap = list(age = age, sex = sexchara, year = diag)), times = nt, extend = TRUE)
    .survfitF2 <- summary( rs.surv(Surv(time, event) ~ .grpsF2, data = .dataF2,
                                   ratetable = slopop, method = "pohar-perme", add.times = nt,
                                   rmap = list(age = age, sex = sexchara, year = diag)), times = nt, extend = TRUE)
    
    .obsP <- .survfitP$surv
    .lowerP <- .survfitP$lower
    .upperP <- .survfitP$upper
    
    .obsF1 <- .survfitF1$surv
    .lowerF1 <-.survfitF1$lower
    .upperF1 <- .survfitF1$upper 
    
    .obsF2 <- .survfitF2$surv
    .lowerF2 <- .survfitF2$lower
    .upperF2 <- .survfitF2$upper
    
    if(hasArg(cex)==FALSE) {cex <-1} else {cex <- list(...)$cex}
    if(hasArg(cex.lab)==FALSE) {cex.lab <- 1} else {cex.lab <- list(...)$cex.lab}
    if(hasArg(cex.axis)==FALSE) {cex.axis <- 1} else {cex.axis <- list(...)$cex.axis}
    if(hasArg(cex.main)==FALSE) {cex.main <- 1} else {cex.main <- list(...)$cex.main}
    if(hasArg(type)==FALSE) {type <- "b"} else {type <- list(...)$type}
    if(hasArg(col)==FALSE) {col <- 1} else {col <- list(...)$col}
    if(hasArg(lty)==FALSE) {lty <- 1} else {lty <- list(...)$lty}
    if(hasArg(lwd)==FALSE) {lwd <- 1} else {lwd <- list(...)$lwd}
    if(hasArg(pch)==FALSE) {pch <- 16} else {pch <- list(...)$pch}
    
    if(hasArg(ylim)==FALSE) {ylim <- c(0,1)} else {ylim <- list(...)$ylim}
    if(hasArg(xlim)==FALSE) {xlim  <- c(0,1)} else {xlim <- list(...)$xlim}
    
    if(hasArg(ylab)==FALSE) {ylab <- "Pohar-Perme estimations"} else {ylab <- list(...)$ylab}
    if(hasArg(xlab)==FALSE) {xlab <- "Models estimations"} else {xlab <- list(...)$xlab}
    if(hasArg(main)==FALSE) {main <- paste0(nt/365.241, " years")} else {main <- list(...)$main}
    
    df_plot <- bind_rows(
      data.frame(Method = "PLANN", est = .estP, obs = .obsP, lower = .lowerP, upper = .upperP
      ),
      data.frame(Method = "RP", est = .estF1, obs = .obsF1, lower = .lowerF1, upper = .upperF1
      ),
      data.frame(Method = "RP2", est = .estF2, obs = .obsF2, lower = .lowerF2, upper = .upperF2
      )
    )
    maxy <- ifelse(max(df_plot$upper) > 1,
                   max(df_plot$upper) + 0.01,
                   1)     
    # ---- 2. ggplot equivalent of the base R graph ----
    p <- ggplot(df_plot, aes(x = est, y = obs, color = Method, linetype = Method)) +
      geom_line(linewidth = 1) +
      geom_point(size = 2) +
      geom_segment(aes(x = est, xend = est, y = lower, yend = upper)) +
      geom_abline(intercept = 0, slope = 1, linetype = "dashed") +
      scale_color_manual(
        values = c(
          "PLANN" = "aquamarine",
          "RP" = "brown1",
          "RP2" = "darkolivegreen"
        )
      ) +
      labs(
        x = "Model estimates",
        y = "Pohar-Perme estimates",
        title = main_name
      ) +
      scale_linetype_manual(values=c("RP" ="dotted", "RP2" ="longdash", "PLANN" ="twodash")) + 
      theme_bw() +
      theme_minimal(base_family = "CMU Serif",base_size = 14) +
      theme(legend.position = "none",
            panel.border = element_blank(),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            plot.title = element_text(hjust = 0.5),
            axis.ticks = element_line(color = "black", size = 0.5),  # tick line
            axis.ticks.length = unit(0.25, "cm"),
            axis.line = element_line(colour = "black")
      ) + scale_x_continuous(limits = c(-0.0025,1.05), expand = c(0,0), breaks = c(0, 0.25, 0.50, 0.75, 1), labels = labels) +
      scale_y_continuous(expand = expansion(mult = 0), breaks = c(0.25, 0.50, 0.75, 1), labels = labels
      )+
      coord_cartesian(ylim = c(0, maxy))
    
    print(p)
    l_cali_valid_strata[[l]] <- p
    time_map <- c(3,5,10)
    ggsave(
      filename = paste0("/home/thomas/Documents/Articles/Article survivalPLANN/Figures/Figures 4 strates NEW/Calibration_VALID_", j,"_",  time_map[which(time_map == nt/365.241)]
                        ,".png"),
      plot = p,
      width = 6,
      height = 6,
      dpi = 800
    )
    l=l+1
  }
  
}

########################## multi plot
library(patchwork)

## train 


# mean cali
leg_mean <- theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                  panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
                  plot.title = element_text(hjust = 0.5),
                  legend.position =  c(0.25, 0.02),
                  legend.justification = c(1, 0),
                  legend.title = element_text(size=8), 
                  legend.text = element_text(size=8),
                  legend.background = element_rect(fill = "white", color = "black"))

leg_cali <- theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                  panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
                  plot.title = element_text(hjust = 0.5),
                  axis.ticks = element_line(color = "black", size = 0.5),  # tick line
                  axis.ticks.length = unit(0.25, "cm"),
                  legend.position = c(0.96, 0.02),
                  legend.justification = c(1, 0),
                  legend.background = element_rect(fill = "white", color = "black") )

combined <- ( (l_mean_train[[1]] + leg_mean + ggtitle("Whole")) + (l_cali_train[[1]]+leg_cali) ) /
  ((l_cali_train[[2]]+leg_cali) + (l_cali_train[[3]]+leg_cali))
combined
ggsave(
  filename = paste0("/home/thomas/Documents/Articles/Article survivalPLANN/Figures/Figures 4 strates NEW/mean_cali_train.png"),
  plot = combined, width = 12, height = 8, dpi = 800)
#mean strata
combined <- ( (l_mean_train[[2]] + leg_mean + ggtitle("Male & Colon") ) + (l_mean_train[[3]] + leg_mean + ggtitle("Male & Rectum")) ) /
  ((l_mean_train[[4]] + leg_mean + ggtitle("Female & Colon") ) + (l_mean_train[[5]] + leg_mean + ggtitle("Female & Rectum")))
combined
ggsave(
  filename = paste0("/home/thomas/Documents/Articles/Article survivalPLANN/Figures/Figures 4 strates NEW/multi_mean_cali_train_whole.png"),
  plot = combined, width = 12, height = 8, dpi = 800)
#means strata cali

#HC
combined <- ( (l_mean_train[[2]] + leg_mean+ ggtitle("Male & Colon")) + (l_cali_train_strata[[1]]+leg_cali) ) /
  ((l_cali_train_strata[[2]]+leg_cali) + (l_cali_train_strata[[3]]+leg_cali))
combined
ggsave(
  filename = paste0("/home/thomas/Documents/Articles/Article survivalPLANN/Figures/Figures 4 strates NEW/mean_cali_train_HC.png"),
  plot = combined, width = 12, height = 8, dpi = 800)

#HR
combined <- ( (l_mean_train[[3]] + leg_mean+ ggtitle("Male & Rectum")) + (l_cali_train_strata[[4]]+leg_cali) ) /
  ((l_cali_train_strata[[5]]+leg_cali) + (l_cali_train_strata[[6]]+leg_cali))
combined
ggsave(
  filename = paste0("/home/thomas/Documents/Articles/Article survivalPLANN/Figures/Figures 4 strates NEW/mean_cali_train_HR.png"),
  plot = combined, width = 12, height = 8, dpi = 800)
#FC
combined <- ( (l_mean_train[[4]] + leg_mean + ggtitle("Female & Colon")) + (l_cali_train_strata[[7]]+leg_cali) ) /
  ((l_cali_train_strata[[8]]+leg_cali) + (l_cali_train_strata[[9]]+leg_cali))
combined
ggsave(
  filename = paste0("/home/thomas/Documents/Articles/Article survivalPLANN/Figures/Figures 4 strates NEW/mean_cali_train_FC.png"),
  plot = combined, width = 12, height = 8, dpi = 800)
#FR
combined <- ( (l_mean_train[[5]] + leg_mean + ggtitle("Female & Rectum")) + (l_cali_train_strata[[10]]+leg_cali) ) /
  ((l_cali_train_strata[[11]]+leg_cali) + (l_cali_train_strata[[12]]+leg_cali))
combined
ggsave(
  filename = paste0("/home/thomas/Documents/Articles/Article survivalPLANN/Figures/Figures 4 strates NEW/mean_cali_train_FR.png"),
  plot = combined, width = 12, height = 8, dpi = 800)

## valid 


# mean cali

combined <- ( (l_mean_valid[[1]] + leg_mean + ggtitle("Whole")) + (l_cali_valid[[1]]+leg_cali) ) /
  ((l_cali_valid[[2]]+leg_cali) + (l_cali_valid[[3]]+leg_cali))
combined
ggsave(
  filename = paste0("/home/thomas/Documents/Articles/Article survivalPLANN/Figures/Figures 4 strates NEW/mean_cali_valid.png"),
  plot = combined, width = 12, height = 8, dpi = 800)
#mean strata
combined <- ( (l_mean_valid[[2]] + leg_mean + ggtitle("Male & Colon") ) + (l_mean_valid[[3]] + leg_mean + ggtitle("Male & Rectum")) ) /
  ((l_mean_valid[[4]] + leg_mean + ggtitle("Female & Colon") ) + (l_mean_valid[[5]] + leg_mean + ggtitle("Female & Rectum")))
combined
ggsave(
  filename = paste0("/home/thomas/Documents/Articles/Article survivalPLANN/Figures/Figures 4 strates NEW/multi_mean_cali_valid_whole.png"),
  plot = combined, width = 12, height = 8, dpi = 800)
#means strata cali

#HC
combined <- ( (l_mean_valid[[2]] + leg_mean + ggtitle("Male & Colon")) + (l_cali_valid_strata[[1]]+leg_cali) ) /
  ((l_cali_valid_strata[[2]]+leg_cali) + (l_cali_valid_strata[[3]]+leg_cali))
combined
ggsave(
  filename = paste0("/home/thomas/Documents/Articles/Article survivalPLANN/Figures/Figures 4 strates NEW/mean_cali_valid_HC.png"),
  plot = combined, width = 12, height = 8, dpi = 800)

#HR
combined <- ( (l_mean_valid[[3]] + leg_mean + ggtitle("Male & Rectum") ) + (l_cali_valid_strata[[4]]+leg_cali) ) /
  ((l_cali_valid_strata[[5]]+leg_cali) + (l_cali_valid_strata[[6]]+leg_cali))
combined
ggsave(
  filename = paste0("/home/thomas/Documents/Articles/Article survivalPLANN/Figures/Figures 4 strates NEW/mean_cali_valid_HR.png"),
  plot = combined, width = 12, height = 8, dpi = 800)
#FC
combined <- ( (l_mean_valid[[4]] + leg_mean + ggtitle("Female & Colon")) + (l_cali_valid_strata[[7]]+leg_cali) ) /
  ((l_cali_valid_strata[[8]]+leg_cali) + (l_cali_valid_strata[[9]]+leg_cali))
combined
ggsave(
  filename = paste0("/home/thomas/Documents/Articles/Article survivalPLANN/Figures/Figures 4 strates NEW/mean_cali_valid_FC.png"),
  plot = combined, width = 12, height = 8, dpi = 800)
#FR
combined <- ( (l_mean_valid[[5]] + leg_mean + ggtitle("Female & Rectum") ) + (l_cali_valid_strata[[10]]+leg_cali) ) /
  ((l_cali_valid_strata[[11]]+leg_cali) + (l_cali_valid_strata[[12]]+leg_cali))
combined
ggsave(
  filename = paste0("/home/thomas/Documents/Articles/Article survivalPLANN/Figures/Figures 4 strates NEW/mean_cali_valid_FR.png"),
  plot = combined, width = 12, height = 8, dpi = 800)

### courbes log(-log()) pour NPH
library(relsurv)

HC = summary( rs.surv(Surv(time, stat) ~ 1, data = colrec[colrec$site == "colon" & colrec$sex ==1,],
                      ratetable = slopop, method = "pohar-perme", add.times = newtimespred,
                      rmap = list(age = age, sex = sexchara, year = diag)), times = newtimespred, extend = TRUE) 

HR = summary( rs.surv(Surv(time, stat) ~ 1, data = colrec[colrec$site == "rectum" & colrec$sex ==1 ,],
                      ratetable = slopop, method = "pohar-perme", add.times = newtimespred,
                      rmap = list(age = age, sex = sexchara, year = diag)), times = newtimespred, extend = TRUE) 



FC = summary( rs.surv(Surv(time, stat) ~ 1, data = colrec[colrec$site == "colon"  & colrec$sex ==2  ,],
                      ratetable = slopop, method = "pohar-perme", add.times = newtimespred,
                      rmap = list(age = age, sex = sexchara, year = diag)), times = newtimespred, extend = TRUE) 



FR = summary( rs.surv(Surv(time, stat) ~ 1, data = colrec[colrec$site == "rectum" & colrec$sex ==2 ,],
                      ratetable = slopop, method = "pohar-perme", add.times = newtimespred,
                      rmap = list(age = age, sex = sexchara, year = diag)), times = newtimespred, extend = TRUE) 

plot(newtimespred/365.241, log(-log(HC$surv)),  type = 'l', ylab = "log(-log(Se(t)))", xlab = "Time (years)")
lines(newtimespred/365.241, log(-log(HR$surv)), col = 'red', lty = 2)
lines(newtimespred/365.241, log(-log(FC$surv)), col = 'green', lty = 3)
lines(newtimespred/365.241, log(-log(FR$surv)), col = 'blue', lty = 4)


legend("bottomright", legend = c("Male & Colon", "Male & Rectum", "Female & Colon", "Female & Rectum"),
       col = c("black", "red", "green", "blue"), lty = c(1,2,3,4))

df <- data.frame( time = newtimespred / 365.241, Male_Colon = log(-log(HC$surv)), 
                  Male_Rectum   = log(-log(HR$surv)), Female_Colon  = log(-log(FC$surv)),
                  Female_Rectum = log(-log(FR$surv)) )

# Convert to long format for ggplot
library(tidyr)
df_long <- pivot_longer(df,cols = -time,names_to = "Group",values_to = "Value")
df_long$Group <- gsub("_", " & ", df_long$Group)
ggplot(df_long, aes(x = time, y = Value, color = Group, linetype = Group)) +
  geom_line(size = 1) +
  labs(x = "Time (years)",y = "log(-log(Se(t)))") +
  scale_color_manual(values = c("Male & Colon" = "black","Male & Rectum" = "red",
                                "Female & Colon"  = "green", "Female & Rectum" = "blue")) +
  scale_linetype_manual(values = c("Male & Colon" = 1, "Male & Rectum" = 2,
                                   "Female & Colon" = 3, "Female & Rectum" = 4) ) +
  scale_x_continuous(
    breaks = scales::pretty_breaks(),  # or: breaks = 0:10
    labels = function(x) format(round(x), nsmall = 0)  # integers only
  ) + 
  theme_minimal(base_family = "Helvetica", base_size = 14) +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(color = "black", linewidth = 0.6),
    axis.ticks = element_line(color = "black", linewidth = 0.5),
    axis.ticks.length = unit(0.20, "cm"),
    plot.title = element_text(hjust = 0.5)  # center title
  ) + theme(legend.title = element_blank(), legend.position = c(0.80, 0.20),   # bottom-right *inside* the panel
            legend.background = element_rect(fill = "white", color = "black", linewidth = 0.5),
            legend.box.background = element_rect(fill = "white", color = "black")
  )
############### TRACER LES COURBES ROCS et donner les valeurs dans un tableau

newtimespred <- sort(unique(c(seq(0, 365.241*10, 365.241/12), newtimes)))[-1]

H = summary( rs.surv(Surv(time, stat) ~ 1, data = colrec[colrec$site =="colon",],
                     ratetable = slopop, method = "pohar-perme", add.times = newtimespred,
                     rmap = list(age = age, sex = sexchara, year = diag)), times = newtimespred, extend = TRUE) 

Fe = summary( rs.surv(Surv(time, stat) ~ 1, data = colrec[colrec$site =="rectum" ,],
                      ratetable = slopop, method = "pohar-perme", add.times = newtimespred,
                      rmap = list(age = age, sex = sexchara, year = diag)), times = newtimespred, extend = TRUE) 

plot(newtimespred/365.241, log(-log(H$surv)),  type = 'l', ylab = "log(-log(Se(t)))", xlab = "Time (years)")
lines(newtimespred/365.241, log(-log(Fe$surv)), col = 'green', lty = 3)

df <- data.frame( time = newtimespred / 365.241, Colon = log(-log(H$surv)), 
                  Rectum  = log(-log(Fe$surv)))

# Convert to long format for ggplot
library(tidyr)
df_long <- pivot_longer(df,cols = -time,names_to = "Group",values_to = "Value")
df_long$Group <- gsub("_", " & ", df_long$Group)
ggplot(df_long, aes(x = time, y = Value, color = Group, linetype = Group)) +
  geom_line(size = 1) +
  labs(x = "Time (years)",y = "log(-log(Se(t)))") +
  scale_color_manual(values = c("Colon" = "black","Rectum" = "red")) +
  scale_linetype_manual(values = c("Colon" = 1, "Rectum" = 2)) +
  scale_x_continuous(
    breaks = scales::pretty_breaks(),  # or: breaks = 0:10
    labels = function(x) format(round(x), nsmall = 0)  # integers only
  ) + 
  theme_minimal(base_family = "CMU Serif", base_size = 14) +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(color = "black", linewidth = 0.6),
    axis.ticks = element_line(color = "black", linewidth = 0.5),
    axis.ticks.length = unit(0.20, "cm"),
    plot.title = element_text(hjust = 0.5)  # center title
  ) + theme(legend.title = element_blank(), legend.position = c(0.80, 0.20),   # bottom-right *inside* the panel
            legend.background = element_rect(fill = "white", color = "black", linewidth = 0.5),
            legend.box.background = element_rect(fill = "white", color = "black")
  )

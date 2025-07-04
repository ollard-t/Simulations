# remotes::install_github("ollard-t/survivalNET")
# remotes::install_github("ollard-t/survivalPLANN")

# install.packages("~/Documents/GitHub/survivalNET", repos = NULL, type = "source")
# install.packages("~/Documents/GitHub/survivalPLANN", repos = NULL, type = "source")

##CTRL+F balises

## 1000SiSTA : N1000 début simulations
## 1000SiEND : N1000 fin simulations#
## 1000IndSTA : N1000 début calcul indicateurs
##     1000ROCSTA
## 1000IndEND : N1000 fin calcul indicateurs

## 3000SiSTA : N3000 début simulations
## 3000SiEND : N3000 fin simulations#
## 3000IndSTA : N3000 début calcul indicateurs
##     3000ROCSTA
## 3000IndEND : N3000 fin calcul indicateurs

## 5000SiSTA : N5000 début simulations
## 5000SiEND : N5000 fin simulations#
## 5000IndSTA : N5000 début calcul indicateurs
##     5000ROCSTA
## 5000IndEND : N5000 fin calcul indicateurs

#################################################################################################
#################################################################################################
#################################################################################################
                                            ##1000##
#################################################################################################
#################################################################################################
#################################################################################################

## 1000SiSTA

library(survivalNET)
library(survivalPLANN)
library(parallel)
library(doParallel)


path0 <- "~/Documents/Simulations/Simulations mai 2025/N1000/Output simulations_N1000_NPH_HC/"
path1 <- "~/Documents/Simulations/BASES/"
############

#importation des données
data("fr.ratetable")
data("slopop")
data(colrec)

# data management


colrec$agey <- colrec$age/365.241

colrec$sexchara <- "male"
colrec$sexchara[colrec$sex==2] <- "female"

colrec$colon <- 1*(as.character(colrec$site)=="colon")

colrec <- colrec[colrec$stage!=99,] # suppresion des patients avec valeur manquante
table(colrec$stage)

colrec$stage1 <- 1*(colrec$stage==1)
colrec$stage2 <- 1*(colrec$stage==2)
colrec$stage3 <- 1*(colrec$stage==3)

table(colrec$stat)

colrec$agey10 <- colrec$agey/10
names(colrec)


###################################################################################################
### Fonctions de survie nette, attendue, observée et de simulation de temps
###################################################################################################


Sn <- function(time, sigma, nu, theta, beta, covariates)
{
  exp((1-(1+(time/(exp(sigma)))^exp(nu))^(1/exp(theta))) * exp(sum(covariates*beta)) )
}


simulate_iteration <- function(i){
  
  set.seed(i)
  
  start <- Sys.time()
  
  data <- read.csv(file = paste0(path1,"ind1000NPH_HC/",i,"_df1000.csv"), sep = ";")
  Tsigma <- 12.6 
  Tnu <- -0.5
  Ttheta <- 0
  betaZ <- c(0.9,2.7,0.3,-0.1,-0.1)
  
  TsigmaHC <- 13.0
  TnuHC <- -0.55
  TthetaHC <- 0
  betaZHC <- c(1,2.9,0.2)
  
  #### Hommes Rectum
  TsigmaHR <- 12.75
  TnuHR <- -0.3
  TthetaHR <- 0                      
  betaZHR <- c(1.2,3,0.35)
  
  #### Femmes Colon
  TsigmaFC <- 14
  TnuFC <- -0.55
  TthetaFC <- 0
  betaZFC <- c(0.6,2.6,-0.1)
  
  #### Femmes Rectum
  TsigmaFR <- 11.5
  TnuFR <- -0.4
  TthetaFR <- 0
  betaZFR <- c(0.9,2.7,0.16)
  ############### validation sample
  N <- dim(data)[1]
  N_train <- N/2 #500
  data_train = data[1:N_train,]
  data_valid = data[(N_train+1):dim(data)[1],]
  
  ###############################
  ########## STRATES ############
  ###############################
  
  ### training
  #################### Hommes ####################
  
  ### hommes & colon
  
  data_train_HC <- data_train[data_train$sex == 1 & data_train$colon == 1,]
  
  ### hommes & rectum
  
  data_train_HR <- data_train[data_train$sex == 1 & data_train$colon == 0,]
  
  #################### Femmes ####################
  
  ### femmes & colon
  
  data_train_FC <- data_train[data_train$sex == 2 & data_train$colon == 1,]
  
  ### femmes & rectum
  
  data_train_FR <- data_train[data_train$sex == 2 & data_train$colon == 0,]
  
  
  ########################################
  ### validation
  
  #################### Hommes ####################
  
  ### hommes & colon
  
  data_valid_HC <- data_valid[data_valid$sex == 1 & data_valid$colon == 1,]
  
  ### hommes & rectum
  
  data_valid_HR <- data_valid[data_valid$sex == 1 & data_valid$colon == 0,]
  
  #################### Femmes ####################
  
  ### femmes & colon
  
  data_valid_FC <- data_valid[data_valid$sex == 2 & data_valid$colon == 1,]
  
  ### femmes & rectum
  
  data_valid_FR <- data_valid[data_valid$sex == 2 & data_valid$colon == 0,]
  
  ########################################
  
  ######################valeurs issues du modèle theorique
  
  pro.time <- max(data_train$times[data_train$status==1])
  
  ###temps d'évaluations (t =1, 3, 5, 10 ans)
  
  newtimes <- c(365.241, 365.241*3, 365.241*5, 365.241*10) 
  # @THOMAS : 80% 50% 20% sujets à risque 
  
  Zdata_train <- cbind(data_train$stage2, data_train$stage3, data_train$agey10, data_train$sex, data_train$colon)
  colnames(Zdata_train) <- c("stage2", "stage3", "agey10", "sex", "colon")
  
  
  funtheo <- function(i){
    if(Zdata_train[i,"sex"] == 1 & Zdata_train[i,"colon"] == 1){
      Tsigma = TsigmaHC
      Tnu = TnuHC
      Ttheta = TthetaHC
      betaZ = betaZHC
    }
    # Hommes rectum
    if(Zdata_train[i,"sex"] == 1 & Zdata_train[i,"colon"] == 0){
      Tsigma = TsigmaHR
      Tnu = TnuHR
      Ttheta = TthetaHR
      betaZ = betaZHR
    }
    # Femmes colon
    if(Zdata_train[i,"sex"] == 2 & Zdata_train[i,"colon"] == 1){
      Tsigma = TsigmaFC
      Tnu = TnuFC
      Ttheta = TthetaFC
      betaZ = betaZFC
    }
    # Femmes rectum
    if(Zdata_train[i,"sex"] == 2 & Zdata_train[i,"colon"] == 0){
      Tsigma = TsigmaFR
      Tnu = TnuFR
      Ttheta = TthetaFR
      betaZ = betaZFR
    }
    Sn(time = newtimes, sigma = Tsigma, nu = Tnu, theta = Ttheta, beta = betaZ, covariates = Zdata_train[i,-c(4,5)]) }
  
  theopred <- t(sapply(1:dim(Zdata_train)[1], FUN = funtheo ))
  colnames(theopred) <- newtimes
  
  mean.theo <- apply(theopred, FUN="mean", MARGIN=2)
  
  ###strates
  #####
  
  strata_names <- c("HC", "HR", "FC", "FR")
  
  for(j in strata_names){
    
    data_train_var <- get(paste0("data_train_", j))
    
    Zdata_train <- cbind(data_train_var$stage2, data_train_var$stage3, data_train_var$agey10, data_train_var$sex, data_train_var$colon)
    colnames(Zdata_train) <- c("stage2", "stage3", "agey10", "sex", "colon")
    
    funtheostrat <- function(i){
      if(Zdata_train[i,"sex"] == 1 & Zdata_train[i,"colon"] == 1){
        Tsigma = TsigmaHC
        Tnu = TnuHC
        Ttheta = TthetaHC
        betaZ = betaZHC
      }
      # Hommes rectum
      if(Zdata_train[i,"sex"] == 1 & Zdata_train[i,"colon"] == 0){
        Tsigma = TsigmaHR
        Tnu = TnuHR
        Ttheta = TthetaHR
        betaZ = betaZHR
      }
      # Femmes colon
      if(Zdata_train[i,"sex"] == 2 & Zdata_train[i,"colon"] == 1){
        Tsigma = TsigmaFC
        Tnu = TnuFC
        Ttheta = TthetaFC
        betaZ = betaZFC
      }
      # Femmes rectum
      if(Zdata_train[i,"sex"] == 2 & Zdata_train[i,"colon"] == 0){
        Tsigma = TsigmaFR
        Tnu = TnuFR
        Ttheta = TthetaFR
        betaZ = betaZFR
      }
      Sn(time = newtimes, sigma = Tsigma, nu = Tnu, theta = Ttheta, beta = betaZ, covariates = Zdata_train[i,-c(4,5)]) }
    
    theopredstrat <- t(sapply(1:dim(Zdata_train)[1], FUN = funtheostrat ))
    
    colnames(theopredstrat) <- newtimes
    
    mean.theo_strat <- apply(theopredstrat, FUN = "mean", MARGIN = 2)
    
    assign(paste0("Zdata_train_", j), Zdata_train)
    assign(paste0("theopred_", j), theopredstrat)
    assign(paste0("mean.theo_", j), mean.theo_strat)  
    
  }
  
  #####
  
  ####### sur la base de validation
  
  Zdata_val <- cbind(data_valid$stage2, data_valid$stage3, data_valid$agey10, data_valid$sex, data_valid$colon)
  
  theopredval <- t(sapply(1:dim(Zdata_val)[1], FUN = function(i){
    Sn(time = newtimes, sigma = Tsigma, nu = Tnu, theta = Ttheta, beta = betaZ,  covariates = Zdata_val[i,]) } ) )
  
  colnames(theopredval) <- newtimes
  
  mean.theoval <- apply(theopredval, FUN="mean", MARGIN=2) 
  
  ###strates
  #####
  
  for(j in strata_names){
    
    data_valid_var <- get(paste0("data_valid_", j))
    
    Zdata_valid <- cbind(data_valid_var$stage2, data_valid_var$stage3, data_valid_var$agey10, data_valid_var$sex, data_valid_var$colon)
    
    theopredvalstrat <- t(sapply(1:nrow(Zdata_valid), FUN = function(i) {
      Sn(time = newtimes, sigma = Tsigma, nu = Tnu, theta = Ttheta, beta = betaZ, covariates = Zdata_valid[i,])
    }))
    
    colnames(theopredvalstrat) <- newtimes
    
    mean.theoval_strat <- apply(theopredvalstrat, FUN = "mean", MARGIN = 2)
    
    assign(paste0("Zdata_val_", j), Zdata_valid)
    assign(paste0("theopredval_", j), theopredvalstrat)
    assign(paste0("mean.theoval_", j), mean.theoval_strat)  
    
  }
  
  
  #####
  #######################################################
  ###### choix des hyper-paramètres par cross-validation
  #######################################################
  
  
  # 
  # ## Modèle PLANN
  # # 
  
  startCVplann <- Sys.time()
  
  tune.plann <- cvPLANN(formula = Surv(times, status) ~ stage2 + stage3 + agey10 + sex01 + colon, pro.time = pro.time,
                        data = data_train, cv = 10, inter= 365.241/12, size=c(2, 4, 8, 10),
                        decay=c(0.01, 0.1), maxit=1000, MaxNWts=10000)
  
  TimeCVplann <- Sys.time() - startCVplann
  
  # 
  # ####################
  # #paramètres optimaux
  # ####################
  # 
  
  # inter <- 365.241/12
  # size <- 4
  # decay <- 0.01
  # maxit <- 1000
  # MaxNWts <- 10000
  
  ##### PLANN
  
  inter <- tune.plann$optimal$inter
  
  size <- tune.plann$optimal$size
  
  decay <- tune.plann$optimal$decay
  
  maxit <- tune.plann$optimal$maxit
  
  MaxNWts <- tune.plann$optimal$MaxNWts  
  
  # ##### FLEXNET
  
  #flex1 
  
  m1.2 <- 2
  m1.4 <- 4
  failedindic1.2 <- NULL
  failedindic1.4 <- NULL
  
  #flex2  
  
  m2.2 <- 2
  m2.4 <- 4
  failedindic2.2 <- NULL
  failedindic2.4 <- NULL
  
  
  ##################
  # Modèles fittés
  ##################
  # 
  
  ### plann
  
  plann.model <- sPLANN(formula = Surv(times, status) ~ stage2 + stage3 + agey10 + sex01 + colon,
                        pro.time = pro.time, data = data_train, inter= inter, size=size, 
                        decay=decay, maxit=maxit, MaxNWts=MaxNWts)
  
  ### Weibull Généralisé 
  
  WG.model <- survivalNET(formula = Surv(times, status) ~ stage2 + stage3 + agey10 + sex01 + colon +
                            ratetable(age, year, sexchara), data = data_train,
                          ratetable=slopop, dist = "genweibull")
  
  logcoeffWG <- tail(WG.model$coefficients, 3)
  
  ################## flex1
  ########### 2 noeuds
  m1.2_original <- m1.2
  m1.2_changed <- FALSE
  
  converged <- FALSE
  
  while (m1.2 >= 0 && !converged) {
    tryCatch({
      flex.model1.2 <- survivalFLEXNET(
        formula = Surv(times, status) ~ stage2 + stage3 + agey10 + sex01 + colon +
          ratetable(age, year, sexchara),
        data = data_train,
        ratetable = slopop,
        m = m1.2
      )
      # If the model succeeds, set converged to TRUE and break out of the loop
      converged <- TRUE
    }, error = function(e) {
      # If an error occurs, decrement m1
      m1.2_changed <<- TRUE
      m1.2 <<- m1.2 - 1
    })
  }
  
  if (m1.2_changed & converged) {
    failedindic1.2 <- paste0( " : this value was decreased from ", m1.2_original)
  }
  
  if (!converged) {
    flex.model1.2 <- WG.model
    failedindic1.2 <- "flex.model1.2 didn't converge, even after decreasing the number of knots incrementaly. WG.model was used instead."
  }
  
  ########### 4 noeuds
  m1.4_original <- m1.4
  m1.4_changed <- FALSE
  
  converged <- FALSE
  
  while (m1.4 >= 0 && !converged) {
    tryCatch({
      flex.model1.4 <- survivalFLEXNET(
        formula = Surv(times, status) ~ stage2 + stage3 + agey10 + sex01 + colon +
          ratetable(age, year, sexchara),
        data = data_train,
        ratetable = slopop,
        m = m1.4
      )
      # If the model succeeds, set converged to TRUE and break out of the loop
      converged <- TRUE
    }, error = function(e) {
      # If an error occurs, decrement m1
      m1.4_changed <<- TRUE
      m1.4 <<- m1.4 - 1
    })
  }
  
  if (m1.4_changed & converged) {
    failedindic1.4 <- paste0( " : this value was decreased from ", m1.4_original)
  }
  
  if (!converged) {
    flex.model1.4 <- WG.model
    failedindic1.4 <- "flex.model1.4 didn't converge, even after decreasing the number of knots incrementaly. WG.model was used instead."
  }
  ### flex2
  ## 2 noeuds 
  m2.2_original <- m2.2
  m2.2_changed <- FALSE
  
  converged <- FALSE
  
  while (m2.2 >= 0 && !converged) {
    tryCatch({flex.model2.2 <-survivalFLEXNET(formula = Surv(times, status) ~ stage2 + stage3 + agey10 + strata(sex.organ) +
                                                ratetable(age, year, sexchara), data = data_train,
                                              ratetable=slopop, m = m2.2)
    # If the model succeeds, set converged to TRUE and break out of the loop
    converged <- TRUE
    }, error = function(e) {
      # If an error occurs, decrement m2
      m2.2_changed <<- TRUE
      m2.2 <<- m2.2 - 1
    })
  }
  
  if (m2.2_changed & converged) {
    failedindic2.2 <- paste0(" this value was decreased from ", m2.2_original)
  }
  
  if (!converged) {
    flex.model2.2 <- flex.model1.2
    failedindic2.2 <- "flex.model2.2 didn't converge, even after decreasing the number of knots incrementaly. flex.model1.2 was used instead."
    
  }
  ## 4 noeuds 
  
  m2.4_original <- m2.4
  m2.4_changed <- FALSE
  
  converged <- FALSE
  
  while (m2.4 >= 0 && !converged) {
    tryCatch({flex.model2.4 <-survivalFLEXNET(formula = Surv(times, status) ~ stage2 + stage3 + agey10 + strata(sex.organ) +
                                                ratetable(age, year, sexchara), data = data_train,
                                              ratetable=slopop, m = m2.4)
    # If the model succeeds, set converged to TRUE and break out of the loop
    converged <- TRUE
    }, error = function(e) {
      # If an error occurs, decrement m2
      m2.4_changed <<- TRUE
      m2.4 <<- m2.4 - 1
    })
  }
  
  if (m2.4_changed & converged) {
    failedindic2.4 <- paste0(" this value was decreased from ", m2.4_original)
  }
  
  if (!converged) {
    flex.model2.4 <- flex.model1.4
    failedindic2.4 <- "flex.model2.4 didn't converge, even after decreasing the number of knots incrementaly. flex.model1.4 was used instead."
    
  }
  
  ##############
  #predictions #
  ##############
  
  #bases entières (entraînement et validation)
  ##################################################################################################
  ########### BASE DE TRAIN
  ###plann
  
  plannpred <- predictRS(plann.model, data = data_train, newtimes = newtimes, ratetable=slopop, 
                         age = "age", year = "year", sex = "sexchara")
  
  mean.plann <- plannpred$mpredictions$net_survival 
  
  ###flexnet
  
  ### modèle PH
  ##2 noeuds
  
  flexpred1.2 <- predict(flex.model1.2, newtimes = newtimes)$predictions
  
  mean.flex1.2 <- apply(flexpred1.2, FUN="mean", MARGIN=2)
  
  ##4 noeuds
  
  flexpred1.4 <- predict(flex.model1.4, newtimes = newtimes)$predictions
  
  mean.flex1.4 <- apply(flexpred1.4, FUN="mean", MARGIN=2)
  
  #### modèle NPH
  
  ##2 noeuds
  
  flexpred2.2 <- predict(flex.model2.2 , newtimes = newtimes)$predictions
  
  mean.flex2.2  <- apply(flexpred2.2 , FUN="mean", MARGIN=2)
  
  ##4 noeuds
  
  flexpred2.4 <- predict(flex.model2.4 , newtimes = newtimes)$predictions
  
  mean.flex2.4  <- apply(flexpred2.4 , FUN="mean", MARGIN=2)
  ### genweibull
  
  WGpred <- predict(WG.model, newtimes = newtimes)$predictions
  
  mean.WG <- apply(WGpred, FUN="mean", MARGIN=2)
  
  
  ### Estimateur de Pohar-Perme 
  
  PP <- summary( rs.surv(Surv(times, status) ~ 1, data = data_train,
                         ratetable = slopop, method = "pohar-perme", add.times = newtimes,
                         rmap = list(age = age, sex = sexchara, year = year)), times = newtimes, extend = TRUE) 
  
  poharpred <- PP$surv 
  
  ###### par strates ######
  ######
  #plann
  
  for(j in strata_names){
    
    data_train_var <- get(paste0("data_train_", j))
    
    plannpredS <- predictRS(plann.model, data = data_train_var, newtimes = newtimes, 
                            ratetable = slopop, age = "age", year = "year", sex = "sexchara")
    
    mean.plannS <- plannpredS$mpredictions$net_survival
    
    assign(paste0("plannpred_", j), plannpredS)
    assign(paste0("mean.plann_", j), mean.plannS)
  }
  
  ###flexnet
  
  ### modèle PH
  
  ## 2 noeuds
  
  for(j in strata_names){
    
    data_train_var <- get(paste0("data_train_", j))
    
    flexpred1.2S <- predict(flex.model1.2, newtimes = newtimes, newdata = data_train_var)$predictions
    
    mean.flex1.2S <- apply(flexpred1.2S, FUN="mean", MARGIN=2)
    
    assign(paste0("flexpred1.2_", j), flexpred1.2S)
    assign(paste0("mean.flex1.2_", j), mean.flex1.2S)
  }
  
  ## 4 noeuds
  
  for(j in strata_names){
    
    data_train_var <- get(paste0("data_train_", j))
    
    flexpred1.4S <- predict(flex.model1.4, newtimes = newtimes, newdata = data_train_var)$predictions
    
    mean.flex1.4S <- apply(flexpred1.4S, FUN="mean", MARGIN=2)
    
    assign(paste0("flexpred1.4_", j), flexpred1.4S)
    assign(paste0("mean.flex1.4_", j), mean.flex1.4S)
  }
  
  #### modèle NPH
  
  ##2 noeuds
  
  for(j in strata_names){
    
    data_train_var <- get(paste0("data_train_", j))
    
    flexpred2.2S <- predict(flex.model2.2, newtimes = newtimes, newdata = data_train_var)$predictions
    
    mean.flex2.2S <- apply(flexpred2.2S, FUN="mean", MARGIN=2)
    
    assign(paste0("flexpred2.2_", j), flexpred2.2S)
    assign(paste0("mean.flex2.2_", j), mean.flex2.2S)
  }
  
  ##4 noeuds
  
  for(j in strata_names){
    
    data_train_var <- get(paste0("data_train_", j))
    
    flexpred2.4S <- predict(flex.model2.4, newtimes = newtimes, newdata = data_train_var)$predictions
    
    mean.flex2.4S <- apply(flexpred2.4S, FUN="mean", MARGIN=2)
    
    assign(paste0("flexpred2.4_", j), flexpred2.4S)
    assign(paste0("mean.flex2.4_", j), mean.flex2.4S)
  }
  
  ### genweibull
  
  for(j in strata_names){
    
    data_train_var <- get(paste0("data_train_", j))
    
    WGpredS <- predict(WG.model, newtimes = newtimes, newdata = data_train_var)$predictions
    
    mean.WGS <- apply(WGpredS, FUN="mean", MARGIN=2)
    
    assign(paste0("WGpred_", j), WGpredS)
    assign(paste0("mean.WG_", j), mean.WGS)
  }
  
  
  ##Estimateur de Pohar-Perme 
  
  for(j in strata_names){
    
    data_train_var <- get(paste0("data_train_", j))
    
    PPS <- summary( rs.surv(Surv(times, status) ~ 1, data = data_train_var,
                            ratetable = slopop, method = "pohar-perme", add.times = newtimes,
                            rmap = list(age = age, sex = sexchara, year = year)), times = newtimes, extend = TRUE) 
    
    poharpredS <- PPS$surv
    
    assign(paste0("PP_", j), PPS)
    assign(paste0("poharpred_", j), poharpredS)
  }
  
  ###### 
  
  ############################################
  ########### BASE DE VALIDATION
  ###plann
  
  plannpredval <-  predictRS(plann.model, data = data_valid, newtimes = newtimes, ratetable=slopop, 
                             age = "age", year = "year", sex = "sexchara")
  
  mean.plannval <- plannpredval$mpredictions$net_survival 
  
  ###flexnet
  
  ### modèle PH
  
  ### 2 noeuds
  
  flexpredval1.2 <- predict(flex.model1.2, newtimes = newtimes, newdata = data_valid)$predictions
  
  mean.flexval1.2 <- apply(flexpredval1.2, FUN="mean", MARGIN=2)
  
  ### 4 noeuds
  
  flexpredval1.4 <- predict(flex.model1.4, newtimes = newtimes, newdata = data_valid)$predictions
  
  mean.flexval1.4 <- apply(flexpredval1.4, FUN="mean", MARGIN=2)
  
  #### modèle NPH
  
  ### 2 noeuds
  
  flexpredval2.2 <- predict(flex.model2.2, newtimes = newtimes, newdata = data_valid)$predictions
  
  mean.flexval2.2 <- apply(flexpredval2.2, FUN="mean", MARGIN=2)
  
  ### 4 noeuds
  
  flexpredval2.4 <- predict(flex.model2.4, newtimes = newtimes, newdata = data_valid)$predictions
  
  mean.flexval2.4 <- apply(flexpredval2.4, FUN="mean", MARGIN=2)
  
  ### genweibull
  
  WGpredval <- predict(WG.model, newtimes = newtimes, newdata = data_valid)$predictions
  
  mean.WGval <- apply(WGpredval, FUN="mean", MARGIN=2)
  
  ### Estimateur de Pohar-Perme 
  
  PPval <- summary( rs.surv(Surv(times, status) ~ 1, data = data_valid,
                            ratetable = slopop, method = "pohar-perme", add.times = newtimes,
                            rmap = list(age = age, sex = sexchara, year = year)), times = newtimes, extend = TRUE) 
  
  poharpredval <- PPval$surv 
  
  ###### par strates ######
  
  
  ###plann
  
  for(j in strata_names){
    
    data_valid_var <- get(paste0("data_valid_", j))
    
    plannpredvalS <- predictRS(plann.model, data = data_valid_var, newtimes = newtimes, 
                               ratetable = slopop, age = "age", year = "year", sex = "sexchara")
    
    mean.plannvalS <- plannpredvalS$mpredictions$net_survival
    
    assign(paste0("plannpredval_", j), plannpredvalS)
    assign(paste0("mean.plannval_", j), mean.plannvalS)
  }
  
  ###flexnet
  
  ### modèle PH
  
  ### 2 noeuds 
  
  for(j in strata_names){
    
    data_valid_var <- get(paste0("data_valid_", j))
    
    flexpredval1.2S <- predict(flex.model1.2, newtimes = newtimes, newdata = data_valid_var)$predictions
    
    mean.flexval1.2S <- apply(flexpredval1.2S, FUN="mean", MARGIN=2)
    
    assign(paste0("flexpredval1.2_", j), flexpredval1.2S)
    assign(paste0("mean.flexval1.2_", j), mean.flexval1.2S)
  }
  
  ### 4 noeuds 
  
  for(j in strata_names){
    
    data_valid_var <- get(paste0("data_valid_", j))
    
    flexpredval1.4S <- predict(flex.model1.4, newtimes = newtimes, newdata = data_valid_var)$predictions
    
    mean.flexval1.4S <- apply(flexpredval1.4S, FUN="mean", MARGIN=2)
    
    assign(paste0("flexpredval1.4_", j), flexpredval1.4S)
    assign(paste0("mean.flexval1.4_", j), mean.flexval1.4S)
  }
  #### modèle NPH
  
  ## 2 noeuds
  
  for(j in strata_names){
    
    data_valid_var <- get(paste0("data_valid_", j))
    
    flexpredval2.2S <- predict(flex.model2.2, newtimes = newtimes, newdata = data_valid_var)$predictions
    
    mean.flexval2.2S <- apply(flexpredval2.2S, FUN="mean", MARGIN=2)
    
    assign(paste0("flexpredval2.2_", j), flexpredval2.2S)
    assign(paste0("mean.flexval2.2_", j), mean.flexval2.2S)
  }
  
  ## 4 noeuds
  
  for(j in strata_names){
    
    data_valid_var <- get(paste0("data_valid_", j))
    
    flexpredval2.4S <- predict(flex.model2.4, newtimes = newtimes, newdata = data_valid_var)$predictions
    
    mean.flexval2.4S <- apply(flexpredval2.4S, FUN="mean", MARGIN=2)
    
    assign(paste0("flexpredval2.4_", j), flexpredval2.4S)
    assign(paste0("mean.flexval2.4_", j), mean.flexval2.4S)
  }
  
  ### genweibull
  
  for(j in strata_names){
    
    data_valid_var <- get(paste0("data_valid_", j))
    
    WGpredvalS <- predict(WG.model, newtimes = newtimes, newdata = data_valid_var)$predictions
    
    mean.WGvalS <- apply(WGpredvalS, FUN="mean", MARGIN=2)
    
    assign(paste0("WGpredval_", j), WGpredvalS)
    assign(paste0("mean.WGval_", j), mean.WGvalS)
  }
  
  ##Estimateur de Pohar-Perme 
  
  for(j in strata_names){
    
    data_valid_var <- get(paste0("data_valid_", j))
    
    PPvalS <- summary( rs.surv(Surv(times, status) ~ 1, data = data_valid_var,
                               ratetable = slopop, method = "pohar-perme", add.times = newtimes,
                               rmap = list(age = age, sex = sexchara, year = year)), times = newtimes, extend = TRUE) 
    
    poharpredvalS <- PPvalS$surv
    
    assign(paste0("PPval_", j), PPvalS)
    assign(paste0("poharpredval_", j), poharpredvalS)
  }
  ##################################################################################################
  
  ######
  #time
  ######
  
  time <- Sys.time()-start
  time <- as.numeric(time, units = 'mins')
  ########################
  #       outputs        #
  ########################
  
  #####bases####
  
  ### TRAIN
  write.table(data_train,  paste0(path0, "DATAFRAMES/TRAIN/WHOLE/",i,"_dftrain.csv"), sep = ";", row.names = F, col.names = T)
  ##strates
  for(j in strata_names){
    
    data_train_var <- get(paste0("data_train_", j))
    
    file_path <- paste0(path0, "DATAFRAMES/TRAIN/", j, "/", i, "_dftrain_", j, ".csv")
    
    write.table(data_train_var, file = file_path, sep = ";", row.names = FALSE, col.names = TRUE)
  }
  ### VALID
  write.table(data_valid,  paste0(path0, "DATAFRAMES/VALID/WHOLE/",i,"_dfvalid.csv"), sep = ";", row.names = F, col.names = T)
  ##strates
  
  for(j in strata_names){
    
    data_valid_var <- get(paste0("data_valid_", j))
    
    file_path <- paste0(path0, "DATAFRAMES/VALID/", j, "/", i, "_dfvalid_", j, ".csv")
    
    write.table(data_valid_var, file = file_path, sep = ";", row.names = FALSE, col.names = TRUE)
  }
  
  ############## paramètres estimés par CV ###########
  
  m1.2value <- paste0(m1.2, failedindic1.2)
  m1.4value <- paste0(m1.4, failedindic1.4)
  
  m2.2value <- paste0(m2.2, failedindic2.2)
  m2.4value <- paste0(m2.4, failedindic2.4)
  
  paramsCV <- as.data.frame(c(inter = inter,size = size,decay = decay, maxit = maxit,
                              MaxNWts = MaxNWts, m1.2 = m1.2value, m1.4 = m1.4value, 
                              m2.2 = m2.2value, m2.4 = m2.4value, logcoeffWG = logcoeffWG))
  
  write.table(paramsCV,  paste0(path0, "PARAM/",i,"_param.csv"))
  
  
  ############# temps écoulé total & par CV de chaque modèle et taille de la base #########################
  
  timeComputation <- as.data.frame(c(time = time, plann = TimeCVplann, N = dim(data)[1]))
  
  write.table(timeComputation,  paste0(path0, "TIME/",i,"_time.csv"))
  
  ##############Log(LIKELIHOOD)####################### 
  
  plann.loglik <- plannpred$loglik  
  
  flex1.2.loglik <- unname(flex.model1.2$loglik[1])
  flex1.4.loglik <- unname(flex.model1.4$loglik[1])
  
  flex2.2.loglik <- unname(flex.model2.2$loglik[1])
  flex2.4.loglik <- unname(flex.model2.4$loglik[1])
  
  WG.loglik <- unname(WG.model$loglik[1])
  
  logliks <- as.data.frame(c(plann = plann.loglik, flex1.2 = flex1.2.loglik,
                             flex1.4 = flex1.4.loglik, flex2.2 = flex2.2.loglik,
                             flex2.4 = flex2.4.loglik, WG = WG.loglik))
  
  write.table(logliks,  paste0(path0, "LOGLIK/",i,"_loglik.csv"))
  
  ##################DATA TRAIN
  #################################### MOYENNES
  #######THEO
  
  mean_theo <- t(data.frame(THEO = mean.theo))
  write.table(mean_theo, paste0(path0, "TRAIN/WHOLE/THEO/mean_survival/",i,"_mean_theo.csv"), sep = ";", row.names = F, col.names = T)
  
  
  ######PLANN
  
  mean_PLANN <- t(data.frame(PLANN = mean.plann))
  write.table(mean_PLANN, paste0(path0, "TRAIN/WHOLE/PLANN/mean_survival/",i,"_mean_PLANN.csv"), sep = ";", row.names = F, col.names = T)
  
  #######FLEX1
  
  ##2 noeuds
  mean_flex1.2 <- t(data.frame(flex1.2=mean.flex1.2))
  write.table(mean_flex1.2, paste0(path0, "TRAIN/WHOLE/FLEX1.2/mean_survival/",i,"_mean_flex12.csv"), sep = ";", row.names = F, col.names = T)
  
  ##4 noeuds
  mean_flex1.4 <- t(data.frame(flex1.4=mean.flex1.4))
  write.table(mean_flex1.4, paste0(path0, "TRAIN/WHOLE/FLEX1.4/mean_survival/",i,"_mean_flex14.csv"), sep = ";", row.names = F, col.names = T)
  
  ######FLEX2
  
  ##2 noeuds
  mean_flex2.2 <- t(data.frame(flex2.2=mean.flex2.2))
  write.table(mean_flex2.2, paste0(path0, "TRAIN/WHOLE/FLEX2.2/mean_survival/",i,"_mean_flex22.csv"), sep = ";", row.names = F, col.names = T)
  
  ##4 noeuds
  mean_flex2.4 <- t(data.frame(flex2.4=mean.flex2.4))
  write.table(mean_flex2.4, paste0(path0, "TRAIN/WHOLE/FLEX2.4/mean_survival/",i,"_mean_flex24.csv"), sep = ";", row.names = F, col.names = T)
  
  ######WG
  
  mean_flexWG <- t(data.frame(flexWG=mean.WG))
  write.table(mean_flexWG, paste0(path0, "TRAIN/WHOLE/WG/mean_survival/",i,"_mean_flexWG.csv"), sep = ";",row.names = F, col.names = T)
  
  ######Estimateur de Pohar-Perme 
  
  write.table(poharpred, paste0(path0, "TRAIN/WHOLE/PP/",i,"PP.csv"), sep = ";",row.names = T, col.names = T)
  ###strates 
  ######
  
  ##theo
  
  mean.theo_STRATES = c()
  
  for (j in strata_names) {
    
    mean.theo_STRATES <- rbind(mean.theo_STRATES, get(paste0("mean.theo_", j) ) )
    rownames(mean.theo_STRATES)[nrow(mean.theo_STRATES)] <- paste0("mean.theo_", j)
  }
  
  mean.theo_STRATES <- data.frame(strat = rownames(mean.theo_STRATES), mean.theo_STRATES)
  write.table(mean.theo_STRATES, paste0(path0, "TRAIN/STRATA/THEO/mean_survival/",i,"mean_theo_strates.csv"), sep = ";",row.names = F, col.names = T)
  
  #PLANN
  
  mean.plann_STRATES = c()
  
  for (j in strata_names) {
    
    mean.plann_STRATES <- rbind(mean.plann_STRATES, get(paste0("mean.plann_", j) ) )
    rownames(mean.plann_STRATES)[nrow(mean.plann_STRATES)] <- paste0("mean.plann_", j)
  }
  
  mean.plann_STRATES <- data.frame(strat = rownames(mean.plann_STRATES), mean.plann_STRATES)
  write.table(mean.plann_STRATES, paste0(path0, "TRAIN/STRATA/PLANN/mean_survival/",i,"mean_plann_strates.csv"), sep = ";",row.names = F, col.names = T)
  
  
  #flex1
  
  ##2 noeuds
  mean.flex1.2_STRATES = c()
  
  for (j in strata_names) {
    
    mean.flex1.2_STRATES <- rbind(mean.flex1.2_STRATES, get(paste0("mean.flex1.2_", j) ) )
    rownames(mean.flex1.2_STRATES)[nrow(mean.flex1.2_STRATES)] <- paste0("mean.flex1.2_", j)
  }
  
  mean.flex1.2_STRATES <- data.frame(strat = rownames(mean.flex1.2_STRATES), mean.flex1.2_STRATES)
  write.table(mean.flex1.2_STRATES, paste0(path0, "TRAIN/STRATA/FLEX1.2/mean_survival/",i,"mean_flex12_strates.csv"), sep = ";",row.names = F, col.names = T)
  
  ##4 noeuds
  
  mean.flex1.4_STRATES = c()
  
  for (j in strata_names) {
    
    mean.flex1.4_STRATES <- rbind(mean.flex1.4_STRATES, get(paste0("mean.flex1.4_", j) ) )
    rownames(mean.flex1.4_STRATES)[nrow(mean.flex1.4_STRATES)] <- paste0("mean.flex1.4_", j)
  }
  
  mean.flex1.4_STRATES <- data.frame(strat = rownames(mean.flex1.4_STRATES), mean.flex1.4_STRATES)
  write.table(mean.flex1.4_STRATES, paste0(path0, "TRAIN/STRATA/FLEX1.4/mean_survival/",i,"mean_flex14_strates.csv"), sep = ";",row.names = F, col.names = T)
  
  #flex2
  
  ## 2 noeuds
  
  mean.flex2.2_STRATES = c()
  
  for (j in strata_names) {
    
    mean.flex2.2_STRATES <- rbind(mean.flex2.2_STRATES, get(paste0("mean.flex2.2_", j) ) )
    rownames(mean.flex2.2_STRATES)[nrow(mean.flex2.2_STRATES)] <- paste0("mean.flex2.2_", j)
  }
  
  mean.flex2.2_STRATES <- data.frame(strat = rownames(mean.flex2.2_STRATES), mean.flex2.2_STRATES)
  write.table(mean.flex2.2_STRATES, paste0(path0, "TRAIN/STRATA/FLEX2.2/mean_survival/",i,"mean_flex22_strates.csv"), sep = ";",row.names = F, col.names = T)
  
  ## 4 noeuds
  
  mean.flex2.4_STRATES = c()
  
  for (j in strata_names) {
    
    mean.flex2.4_STRATES <- rbind(mean.flex2.4_STRATES, get(paste0("mean.flex2.4_", j) ) )
    rownames(mean.flex2.4_STRATES)[nrow(mean.flex2.4_STRATES)] <- paste0("mean.flex2.4_", j)
  }
  
  mean.flex2.4_STRATES <- data.frame(strat = rownames(mean.flex2.4_STRATES), mean.flex2.4_STRATES)
  write.table(mean.flex2.4_STRATES, paste0(path0, "TRAIN/STRATA/FLEX2.4/mean_survival/",i,"mean_flex24_strates.csv"), sep = ";",row.names = F, col.names = T)
  
  #WG
  mean.WG_STRATES = c()
  
  for (j in strata_names) {
    
    mean.WG_STRATES <- rbind(mean.WG_STRATES, get(paste0("mean.WG_", j) ) )
    rownames(mean.WG_STRATES)[nrow(mean.WG_STRATES)] <- paste0("mean.WG_", j)
  }
  
  mean.WG_STRATES <- data.frame(strat = rownames(mean.WG_STRATES), mean.WG_STRATES)
  write.table(mean.WG_STRATES, paste0(path0, "TRAIN/STRATA/WG/mean_survival/",i,"mean_WG_strates.csv"), sep = ";",row.names = F, col.names = T)
  
  
  ####Estimateur de Pohar-Perme
  
  poharpred_strates = c()
  
  for (j in strata_names) {
    
    poharpred_strates <- rbind(poharpred_strates, get(paste0("poharpred_", j) ) )
    rownames(poharpred_strates)[nrow(poharpred_strates)] <- paste0("poharpred_", j)
  }
  
  poharpred_strates <- data.frame(strat = rownames(poharpred_strates), poharpred_strates)
  write.table(poharpred_strates, paste0(path0, "TRAIN/STRATA/PP/",i,"PP_strates.csv"), sep = ";",row.names = F, col.names = T)
  
  ######
  
  #################################### SURVIES INDIVIDUELLES
  
  #######THEO
  
  ind_theo <- data.frame(theopred)
  write.table(ind_theo, paste0(path0, "TRAIN/WHOLE/THEO/individual_survival/",
                               i,"_ind_theo.csv"), sep = ";", row.names = F, col.names = T)
  
  ######PLANN
  
  ind_PLANN <- data.frame(plannpred$ipredictions$relative_survival)
  write.table(ind_PLANN, paste0(path0, "TRAIN/WHOLE/PLANN/individual_survival/",
                                i,"_ind_PLANN.csv"), sep = ";", row.names = F, col.names = T)
  
  #######FLEX1
  
  ##2 noeuds
  ind_flex1.2 <- data.frame(flexpred1.2)
  write.table(ind_flex1.2, paste0(path0, "TRAIN/WHOLE/FLEX1.2/individual_survival/"
                                  ,i,"_ind_flex12.csv"), sep = ";", row.names = F, col.names = T)
  
  ##4 noeuds
  ind_flex1.4 <- data.frame(flexpred1.4)
  write.table(ind_flex1.4, paste0(path0, "TRAIN/WHOLE/FLEX1.4/individual_survival/"
                                  ,i,"_ind_flex14.csv"), sep = ";", row.names = F, col.names = T)
  ######FLEX2
  
  ##2 noeuds
  ind_flex2.2 <- data.frame(flexpred2.2)
  write.table(ind_flex2.2, paste0(path0, "TRAIN/WHOLE/FLEX2.2/individual_survival/"
                                  ,i,"_ind_flex22.csv"), sep = ";", row.names = F, col.names = T)
  
  ##4 noeuds
  ind_flex2.4 <- data.frame(flexpred2.4)
  write.table(ind_flex2.4, paste0(path0, "TRAIN/WHOLE/FLEX2.4/individual_survival/"
                                  ,i,"_ind_flex24.csv"), sep = ";", row.names = F, col.names = T)
  
  ######WG
  
  ind_WG <- data.frame(WGpred)
  write.table(ind_WG, paste0(path0, "TRAIN/WHOLE/WG/individual_survival/"
                             ,i,"_ind_WG.csv"), sep = ";",row.names = F, col.names = T)
  
  
  ###strates 
  ######
  
  #theo
  
  for (j in strata_names) {
    
    ind_theo <- data.frame(get(paste0("theopred_", j)))
    
    write.table(ind_theo, paste0(path0, "TRAIN/STRATA/THEO/individual_survival/", i, "_ind_theo_", j, ".csv"),
                sep = ";",row.names = F, col.names = T)
  }
  
  ##plann
  
  for (j in strata_names) {
    
    ind_PLANN <- data.frame(get(paste0("plannpred_", j))$ipredictions$relative_survival)
    
    write.table(ind_PLANN, paste0(path0, "TRAIN/STRATA/PLANN/individual_survival/", i, "_ind_PLANN_", j, ".csv"),
                sep = ";",row.names = F, col.names = T)
  }
  
  #flex1
  
  ## 2 noeuds
  for (j in strata_names) {
    
    ind_flex1.2 <- data.frame(get(paste0("flexpred1.2_", j)))
    
    write.table(ind_flex1.2, paste0(path0, "TRAIN/STRATA/FLEX1.2/individual_survival/", i, "_ind_flex12_", j, ".csv"),
                sep = ";",row.names = F, col.names = T)
  }
  
  ## 4 noeuds
  for (j in strata_names) {
    
    ind_flex1.4 <- data.frame(get(paste0("flexpred1.4_", j)))
    
    write.table(ind_flex1.4, paste0(path0, "TRAIN/STRATA/FLEX1.4/individual_survival/", i, "_ind_flex14_", j, ".csv"),
                sep = ";",row.names = F, col.names = T)
  }
  
  #flex2
  
  ##2 noeuds
  for (j in strata_names) {
    
    ind_flex2.2 <- data.frame(get(paste0("flexpred2.2_", j)))
    
    write.table(ind_flex2.2, paste0(path0, "TRAIN/STRATA/FLEX2.2/individual_survival/", i, "_ind_flex22_", j, ".csv"),
                sep = ";",row.names = F, col.names = T)
  }
  
  ##4 noeuds
  for (j in strata_names) {
    
    ind_flex2.4 <- data.frame(get(paste0("flexpred2.4_", j)))
    
    write.table(ind_flex2.4, paste0(path0, "TRAIN/STRATA/FLEX2.4/individual_survival/", i, "_ind_flex24_", j, ".csv"),
                sep = ";",row.names = F, col.names = T)
  }
  #WG
  
  for (j in strata_names) {
    
    ind_WG <- data.frame(get(paste0("WGpred_", j)))
    
    write.table(ind_WG, paste0(path0, "TRAIN/STRATA/WG/individual_survival/", i, "_ind_WG_", j, ".csv"),
                sep = ";",row.names = F, col.names = T)
  }
  
  ######
  
  ##################DATA VALID
  #################################### MOYENNES
  
  #######THEO
  
  mean_theoval <- t(data.frame(THEO = mean.theoval))
  write.table(mean_theoval, paste0(path0, "VALID/WHOLE/THEO/mean_survival/",
                                   i,"_mean_theoval.csv"), sep = ";", row.names = F, col.names = T)
  
  ######PLANN
  
  mean_PLANNval <- t(data.frame(PLANN = mean.plannval))
  write.table(mean_PLANNval, paste0(path0, "VALID/WHOLE/PLANN/mean_survival/",
                                    i,"_mean_PLANNval.csv"), sep = ";", row.names = F, col.names = T)
  
  #######FLEX1
  
  ##2 noeuds
  mean_flexval1.2 <- t(data.frame(flex1.2=mean.flexval1.2))
  write.table(mean_flexval1.2, paste0(path0, "VALID/WHOLE/FLEX1.2/mean_survival/"
                                      ,i,"_mean_flexval12.csv"), sep = ";", row.names = F, col.names = T)
  ##4 noeuds
  mean_flexval1.4 <- t(data.frame(flex1.4=mean.flexval1.4))
  write.table(mean_flexval1.4, paste0(path0, "VALID/WHOLE/FLEX1.4/mean_survival/"
                                      ,i,"_mean_flexval14.csv"), sep = ";", row.names = F, col.names = T)
  
  ######FLEX2
  
  ## 2 noeuds
  mean_flexval2.2 <- t(data.frame(flex2.2=mean.flexval2.2))
  write.table(mean_flexval2.2, paste0(path0, "VALID/WHOLE/FLEX2.2/mean_survival/"
                                      ,i,"_mean_flexval22.csv"), sep = ";", row.names = F, col.names = T)
  
  ## 4 noeuds
  mean_flexval2.4 <- t(data.frame(flex2.4=mean.flexval2.4))
  write.table(mean_flexval2.4, paste0(path0, "VALID/WHOLE/FLEX2.4/mean_survival/"
                                      ,i,"_mean_flexval24.csv"), sep = ";", row.names = F, col.names = T)
  ######WG
  
  mean_flexWGval <- t(data.frame(flexWG=mean.WGval))
  write.table(mean_flexWGval, paste0(path0, "VALID/WHOLE/WG/mean_survival/"
                                     ,i,"_mean_flexWGval.csv"), sep = ";",row.names = F, col.names = T)
  
  ###### Estimateur de pohar perme
  
  write.table(poharpred, paste0(path0, "VALID/WHOLE/PP/",i,"PPval.csv"), sep = ";",row.names = T, col.names = T)
  
  ###strates 
  ######
  
  ##theo
  
  mean.theoval_STRATES <- c()
  
  for (j in strata_names) {
    
    mean.theoval_STRATES <- rbind(mean.theoval_STRATES, get(paste0("mean.theoval_", j) ) )
    rownames(mean.theoval_STRATES)[nrow(mean.theoval_STRATES)] <- paste0("mean.theoval_", j)
  }
  
  mean.theoval_STRATES <- data.frame(strat = rownames(mean.theo_STRATES), mean.theo_STRATES)
  write.table(mean.theoval_STRATES, paste0(path0, "VALID/STRATA/THEO/mean_survival/",i,"mean_theoval_strates.csv"), sep = ";",row.names = F, col.names = T)
  
  
  #PLANN
  
  mean.plannval_STRATES = c()
  
  for (j in strata_names) {
    
    mean.plannval_STRATES <- rbind(mean.plannval_STRATES, get(paste0("mean.plannval_", j) ) )
    rownames(mean.plannval_STRATES)[nrow(mean.plannval_STRATES)] <- paste0("mean.plannval_", j)
  }
  
  mean.plannval_STRATES <- data.frame(strat = rownames(mean.plannval_STRATES), mean.plannval_STRATES)
  write.table(mean.plannval_STRATES, paste0(path0, "VALID/STRATA/PLANN/mean_survival/",i,"mean_plannval_strates.csv"), sep = ";",row.names = F, col.names = T)
  
  #flex1
  ## 2 noeuds
  mean.flexval1.2_STRATES = c()
  
  for (j in strata_names) {
    
    mean.flexval1.2_STRATES <- rbind(mean.flexval1.2_STRATES, get(paste0("mean.flexval1.2_", j) ) )
    rownames(mean.flexval1.2_STRATES)[nrow(mean.flexval1.2_STRATES)] <- paste0("mean.flexval1.2_", j)
  }
  
  mean.flexval1.2_STRATES <- data.frame(strat = rownames(mean.flexval1.2_STRATES), mean.flexval1.2_STRATES)
  write.table(mean.flexval1.2_STRATES, paste0(path0, "VALID/STRATA/FLEX1.2/mean_survival/",i,"mean_flexval12_strates.csv"), sep = ";",row.names = F, col.names = T)
  
  ##4 noeuds
  mean.flexval1.4_STRATES = c()
  
  for (j in strata_names) {
    
    mean.flexval1.4_STRATES <- rbind(mean.flexval1.4_STRATES, get(paste0("mean.flexval1.4_", j) ) )
    rownames(mean.flexval1.4_STRATES)[nrow(mean.flexval1.4_STRATES)] <- paste0("mean.flexval1.4_", j)
  }
  
  mean.flexval1.4_STRATES <- data.frame(strat = rownames(mean.flexval1.4_STRATES), mean.flexval1.4_STRATES)
  write.table(mean.flexval1.4_STRATES, paste0(path0, "VALID/STRATA/FLEX1.4/mean_survival/",i,"mean_flexval14_strates.csv"), sep = ";",row.names = F, col.names = T)
  
  
  #flex2
  
  ##2 noeuds
  mean.flexval2.2_STRATES = c()
  
  for (j in strata_names) {
    
    mean.flexval2.2_STRATES <- rbind(mean.flexval2.2_STRATES, get(paste0("mean.flexval2.2_", j) ) )
    rownames(mean.flexval2.2_STRATES)[nrow(mean.flexval2.2_STRATES)] <- paste0("mean.flexval2.2_", j)
  }
  
  mean.flexval2.2_STRATES <- data.frame(strat = rownames(mean.flexval2.2_STRATES), mean.flexval2.2_STRATES)
  write.table(mean.flexval2.2_STRATES, paste0(path0, "VALID/STRATA/FLEX2.2/mean_survival/",i,"mean_flexval22_strates.csv"), sep = ";",row.names = F, col.names = T)
  
  ##4 noeuds
  mean.flexval2.4_STRATES = c()
  
  for (j in strata_names) {
    
    mean.flexval2.4_STRATES <- rbind(mean.flexval2.4_STRATES, get(paste0("mean.flexval2.4_", j) ) )
    rownames(mean.flexval2.4_STRATES)[nrow(mean.flexval2.4_STRATES)] <- paste0("mean.flexval2.4_", j)
  }
  
  mean.flexval2.4_STRATES <- data.frame(strat = rownames(mean.flexval2.4_STRATES), mean.flexval2.4_STRATES)
  write.table(mean.flexval2.4_STRATES, paste0(path0, "VALID/STRATA/FLEX2.4/mean_survival/",i,"mean_flexval24_strates.csv"), sep = ";",row.names = F, col.names = T)
  
  #WG
  mean.WGval_STRATES = c()
  
  for (j in strata_names) {
    
    mean.WGval_STRATES <- rbind(mean.WGval_STRATES, get(paste0("mean.WGval_", j) ) )
    rownames(mean.WGval_STRATES)[nrow(mean.WGval_STRATES)] <- paste0("mean.WGval_", j)
  }
  
  mean.WGval_STRATES <- data.frame(strat = rownames(mean.WGval_STRATES), mean.WGval_STRATES)
  write.table(mean.WGval_STRATES, paste0(path0, "VALID/STRATA/WG/mean_survival/",i,"mean_WGval_strates.csv"), sep = ";",row.names = F, col.names = T)
  
  ####Estimateur de Pohar-Perme
  
  poharpredval_strates = c()
  
  for (j in strata_names) {
    
    poharpredval_strates <- rbind(poharpredval_strates, get(paste0("poharpredval_", j) ) )
    rownames(poharpredval_strates)[nrow(poharpredval_strates)] <- paste0("poharpredval_", j)
  }
  
  poharpredval_strates <- data.frame(strat = rownames(poharpredval_strates), poharpredval_strates)
  write.table(poharpredval_strates, paste0(path0, "VALID/STRATA/PP/",i,"PPval_strates.csv"), sep = ";",row.names = F, col.names = T)
  
  ######
  #################################### SURVIES INDIVIDUELLES
  
  #######THEO
  ind_theoval <- data.frame(theopredval)
  write.table(ind_theoval, 
              paste0(path0, "VALID/WHOLE/THEO/individual_survival/",i
                     ,"_ind_theoval.csv"), sep = ";", row.names = F, col.names = T)
  
  ######PLANN
  
  ind_PLANNval <- data.frame(plannpredval$ipredictions$relative_survival)
  write.table(ind_PLANNval, paste0(path0, "VALID/WHOLE/PLANN/individual_survival/",
                                   i,"_ind_PLANNval.csv"), sep = ";", row.names = F, col.names = T)
  
  #######FLEX1
  
  ##2 noeuds
  
  ind_flexval1.2 <- data.frame(flexpredval1.2)
  write.table(ind_flexval1.2, paste0(path0, "VALID/WHOLE/FLEX1.2/individual_survival/"
                                     ,i,"_ind_flexval12.csv"), sep = ";", row.names = F, col.names = T)
  
  ##4 noeuds
  
  ind_flexval1.4 <- data.frame(flexpredval1.4)
  write.table(ind_flexval1.4, paste0(path0, "VALID/WHOLE/FLEX1.4/individual_survival/"
                                     ,i,"_ind_flexval14.csv"), sep = ";", row.names = F, col.names = T)
  ######FLEX2
  
  ##2 noeuds
  ind_flexval2.2 <- data.frame(flexpredval2.2)
  write.table(ind_flexval2.2, paste0(path0, "VALID/WHOLE/FLEX2.2/individual_survival/"
                                     ,i,"_ind_flexval22.csv"), sep = ";", row.names = F, col.names = T)
  
  ##4 noeuds
  ind_flexval2.4 <- data.frame(flexpredval2.4)
  write.table(ind_flexval2.4, paste0(path0, "VALID/WHOLE/FLEX2.4/individual_survival/"
                                     ,i,"_ind_flexval24.csv"), sep = ";", row.names = F, col.names = T)
  ######WG
  
  ind_flexWGval <- data.frame(WGpredval)
  write.table(ind_flexWGval, paste0(path0, "VALID/WHOLE/WG/individual_survival/"
                                    ,i,"_ind_flexWGval.csv"), sep = ";",row.names = F, col.names = T)
  
  
  ###strates ####
  ######
  
  #theo
  
  for (j in strata_names) {
    
    ind_theoval <- data.frame(get(paste0("theopredval_", j)))
    
    write.table(ind_theoval, paste0(path0, "VALID/STRATA/THEO/individual_survival/", i, "_ind_theoval_", j, ".csv"),
                sep = ";",row.names = F, col.names = T)
  }
  
  ##plann
  
  for (j in strata_names) {
    
    ind_PLANNval <- data.frame(get(paste0("plannpredval_", j))$ipredictions$relative_survival)
    
    write.table(ind_PLANNval, paste0(path0, "VALID/STRATA/PLANN/individual_survival/", i, "_ind_PLANNval_", j, ".csv"),
                sep = ";",row.names = F, col.names = T)
  }
  
  #flex1
  
  ##2 noeuds
  for (j in strata_names) {
    
    ind_flexval1.2 <- data.frame(get(paste0("flexpredval1.2_", j)))
    
    write.table(ind_flexval1.2, paste0(path0, "VALID/STRATA/FLEX1.2/individual_survival/", i, "_ind_flexval12_", j, ".csv"),
                sep = ";",row.names = F, col.names = T)
  }
  ##4 noeuds
  for (j in strata_names) {
    
    ind_flexval1.4 <- data.frame(get(paste0("flexpredval1.4_", j)))
    
    write.table(ind_flexval1.4, paste0(path0, "VALID/STRATA/FLEX1.4/individual_survival/", i, "_ind_flexval14_", j, ".csv"),
                sep = ";",row.names = F, col.names = T)
  }
  #flex2
  
  ##2 noeuds
  for (j in strata_names) {
    
    ind_flexval2.2 <- data.frame(get(paste0("flexpredval2.2_", j)))
    
    write.table(ind_flexval2.2, paste0(path0, "VALID/STRATA/FLEX2.2/individual_survival/", i, "_ind_flexval22_", j, ".csv"),
                sep = ";",row.names = F, col.names = T)
  }
  
  ##4 noeuds
  for (j in strata_names) {
    
    ind_flexval2.4 <- data.frame(get(paste0("flexpredval2.4_", j)))
    
    write.table(ind_flexval2.4, paste0(path0, "VALID/STRATA/FLEX2.4/individual_survival/", i, "_ind_flexval24_", j, ".csv"),
                sep = ";",row.names = F, col.names = T)
  }
  
  #WG
  
  for (j in strata_names) {
    
    ind_WGval <- data.frame(get(paste0("WGpredval_", j)))
    
    write.table(ind_WGval, paste0(path0, "VALID/STRATA/WG/individual_survival/", i, "_ind_WGval_", j, ".csv"),
                sep = ";",row.names = F, col.names = T)
  }
  ######
  
}

iterations <- 1:1000

mclapply(iterations, simulate_iteration, mc.cores = detectCores() - 14)

max_attempts <- 3
attempt <- 1
while (TRUE) {
  indic <- c()  
  
  for (i in 1:1000) {
    if (!file.exists(paste0(path0, "DATAFRAMES/TRAIN/WHOLE/", i, "_dftrain.csv"))) {
      indic <- c(indic, i)
    }
  }
  
  if (length(indic) == 0) {
    break
  }
  if (attempt > max_attempts) {
    break
  }
  mclapply(indic, simulate_iteration, mc.cores = detectCores() - 14)
  attempt <- attempt + 1
}

rm(list = ls(envir = .GlobalEnv), envir = .GlobalEnv)

## 1000SiEND

############################

## 1000IndSTA
library(parallel)
library(doParallel)

path0 <- "~/Documents/Simulations/Simulations mai 2025/N1000/Output simulations_N1000_NPH_HC/"

indic = c()
for(i in 1:1000){
  if(!file.exists(paste0(path0,"DATAFRAMES/TRAIN/WHOLE/",i,"_dftrain.csv"))){ indic <<- c(indic, i)} 
}

if(!is.null(indic)){
  iterations <- c(1:1000)[-indic]
}else{
  iterations <- c(1:1000)
}
strata_names <- c("HC", "HR", "FC", "FR")
newtimes <- c(365.241, 365.241*3, 365.241*5, 365.241*10) 

### FIN INDICATEURS ~ligne 3113

start <- Sys.time()

### Calculs des indicateurs

### RMSE & BIAIS
############################################ TRAIN #################################################
#WHOLE

###on récupère les valeurs théoriques moyennes Se(t) 
ALL_theo <- data.frame()

for(k in iterations){
  assign(paste0("meanTrainTheo_", k) , read.csv(paste0(path0, "TRAIN/WHOLE/THEO/mean_survival/",k,"_mean_theo.csv"),
                                                sep = ";") )
  
  ALL_theo <-rbind(ALL_theo, get(paste0("meanTrainTheo_", k)))
  rm(list =paste0("meanTrainTheo_", k))
}

#moyenne théorique de la valur de survie
S_theo <- colMeans(ALL_theo)
############################################
#importation des données des différents modèles 
#########
##plann
ALL_PLANN <- data.frame() 
for(k in iterations){
  assign(paste0("meanTrainPLANN_", k) , read.csv(paste0(path0, "TRAIN/WHOLE/PLANN/mean_survival/",k,"_mean_PLANN.csv"), sep = ";") )
  ALL_PLANN <-rbind(ALL_PLANN, get(paste0("meanTrainPLANN_", k)))
  rm(list =paste0("meanTrainPLANN_", k))
}

###flex1
##2 noeuds
ALL_flex1.2 <- data.frame()
for(k in iterations){
  assign(paste0("meanTrainflex1.2_", k) , read.csv(paste0(path0, "TRAIN/WHOLE/FLEX1.2/mean_survival/",k,"_mean_flex12.csv"), sep = ";") )
  ALL_flex1.2 <- rbind(ALL_flex1.2, get(paste0("meanTrainflex1.2_",k)))
  rm(list =paste0("meanTrainflex1.2_", k))
}

##4 noeuds
ALL_flex1.4 <- data.frame()
for(k in iterations){
  assign(paste0("meanTrainflex1.4_", k) , read.csv(paste0(path0, "TRAIN/WHOLE/FLEX1.4/mean_survival/",k,"_mean_flex14.csv"), sep = ";") )
  ALL_flex1.4 <- rbind(ALL_flex1.4, get(paste0("meanTrainflex1.4_",k)))
  rm(list =paste0("meanTrainflex1.4_", k))
}

###flex2

##2 noeuds
ALL_flex2.2 <- data.frame()
for(k in iterations){
  assign(paste0("meanTrainflex2.2_", k) , read.csv(paste0(path0, "TRAIN/WHOLE/FLEX2.2/mean_survival/",k,"_mean_flex22.csv"), sep = ";") )
  ALL_flex2.2 <- rbind(ALL_flex2.2, get(paste0("meanTrainflex2.2_",k)))
  rm(list =paste0("meanTrainflex2.2_", k))
}

##4 noeuds
ALL_flex2.4 <- data.frame()
for(k in iterations){
  assign(paste0("meanTrainflex2.4_", k) , read.csv(paste0(path0, "TRAIN/WHOLE/FLEX2.4/mean_survival/",k,"_mean_flex24.csv"), sep = ";") )
  ALL_flex2.4 <- rbind(ALL_flex2.4, get(paste0("meanTrainflex2.4_",k)))
  rm(list =paste0("meanTrainflex2.4_", k))
}
### WG
ALL_WG <- data.frame()
for(k in iterations){
  assign(paste0("meanTrainWG_", k) , read.csv(paste0(path0, "TRAIN/WHOLE/WG/mean_survival/",k,"_mean_flexWG.csv"), sep = ";") )
  ALL_WG <- rbind(ALL_WG,get(paste0("meanTrainWG_",k))) 
  rm(list =paste0("meanTrainWG_", k))
}

ALL_PP <- data.frame() 
for(k in iterations){
  assign(paste0("meanTrainPP_", k) , read.csv(paste0(path0, "TRAIN/WHOLE/PP/",k,"PP.csv"), sep = ";") )
  ALL_PP <-rbind(ALL_PP, get(paste0("meanTrainPP_", k)))
  rm(list =paste0("meanTrainPP_", k))
}
ALL_PP <- as.data.frame(matrix(ALL_PP$x, ncol = 4, byrow = TRUE))

biais_PLANN <- colMeans(ALL_PLANN[,-1] - S_theo[col( ALL_PLANN[,-1] )] )
RMSE_PLANN <- sqrt( colMeans((ALL_PLANN[,-1] - S_theo[col( ALL_PLANN[,-1] )])^2) )

biais_flex1.2 <- colMeans( ALL_flex1.2 - S_theo[col( ALL_flex1.2)] )
RMSE_flex1.2 <- sqrt( colMeans( (ALL_flex1.2 - S_theo[col( ALL_flex1.2)])^2 ) )

biais_flex1.4 <- colMeans( ALL_flex1.4 - S_theo[col( ALL_flex1.4)] )
RMSE_flex1.4 <- sqrt( colMeans( (ALL_flex1.4 - S_theo[col( ALL_flex1.4)])^2 ) )

biais_flex2.2 <- colMeans( ALL_flex2.2 - S_theo[col( ALL_flex2.2)] )
RMSE_flex2.2 <- sqrt( colMeans( (ALL_flex2.2 - S_theo[col( ALL_flex2.2)])^2 ) )

biais_flex2.4 <- colMeans( ALL_flex2.4 - S_theo[col( ALL_flex2.4)] )
RMSE_flex2.4 <- sqrt( colMeans( (ALL_flex2.4 - S_theo[col( ALL_flex2.4)])^2 ) )

biais_WG <- colMeans( ALL_WG - S_theo[col( ALL_WG)] )
RMSE_WG <- sqrt( colMeans( (ALL_WG - S_theo[col( ALL_WG)])^2 ) )

biais_PP <- colMeans( ALL_PP - S_theo[col( ALL_PP)] )
RMSE_PP <- sqrt( colMeans( (ALL_PP - S_theo[col( ALL_PP)])^2 ) )
#STRATES

###on récupère les valeurs théoriques moyennes Se(t) 
for(j in strata_names){
  assign(paste0("ALL_",j, "_theo"),  data.frame()) 
}

for(k in iterations){
  assign(paste0("meanTrainTheoS_", k) , read.csv(paste0(path0, "TRAIN/STRATA/THEO/mean_survival/",k,"mean_theo_strates.csv"),
                                                 sep = ";", row.names = 1) ) 
  for(j in strata_names){
    
    current_row <- get(paste0("meanTrainTheoS_", k))[paste0("mean.theo_", j), , drop = FALSE]
    rownames(current_row) <- paste0("mean.theo_", j, k)
    assign(paste0("ALL_", j, "_theo"), rbind(get(paste0("ALL_", j, "_theo")), current_row))
    
  }
  
  rm(list =paste0("meanTrainTheoS_", k))
}

#moyenne théorique de la valeur de survie
for(j in strata_names){
  assign(paste0("S_theo", j ),colMeans(get(paste0("ALL_",j,"_theo") ) ) )
}
############################################
#importation des données des différents modèles 
############################################

# PLANN

for(j in strata_names){
  assign(paste0("ALL_PLANN_", j),data.frame()) 
}

for(k in iterations){
  assign(paste0("meanTrainPLANNstr_", k) , read.csv(paste0(path0, "TRAIN/STRATA/PLANN/mean_survival/",k,"mean_plann_strates.csv"), sep = ";", row.names = 1) )
  
  for(j in strata_names){
    assign(paste0("ALL_PLANN_",j) ,rbind(get(paste0("ALL_PLANN_",j)),get(paste0("meanTrainPLANNstr_",k))[paste0("mean.plann_", j),-1]  ) )
  }
  rm(list =paste0("meanTrainPLANNstr_", k))
}

### flex1

##2 noeuds
for(j in strata_names){
  assign(paste0("ALL_FLEX1.2_", j),data.frame()) 
}

for(k in iterations){
  assign(paste0("meanTrainFLEX1.2str_", k) , read.csv(paste0(path0, "TRAIN/STRATA/FLEX1.2/mean_survival/",k,"mean_flex12_strates.csv"), sep = ";", row.names = 1) )
  
  for(j in strata_names){
    assign(paste0("ALL_FLEX1.2_",j) ,rbind(get(paste0("ALL_FLEX1.2_",j)), get(paste0("meanTrainFLEX1.2str_",k))[paste0("mean.flex1.2_", j),]  ) )
  }
  rm(list =paste0("meanTrainFLEX1.2str_", k))
}

##4 noeuds
for(j in strata_names){
  assign(paste0("ALL_FLEX1.4_", j),data.frame()) 
}

for(k in iterations){
  assign(paste0("meanTrainFLEX1.4str_", k) , read.csv(paste0(path0, "TRAIN/STRATA/FLEX1.4/mean_survival/",k,"mean_flex14_strates.csv"), sep = ";", row.names = 1) )
  
  for(j in strata_names){
    assign(paste0("ALL_FLEX1.4_",j) ,rbind(get(paste0("ALL_FLEX1.4_",j)), get(paste0("meanTrainFLEX1.4str_",k))[paste0("mean.flex1.4_", j),]  ) )
  }
  rm(list =paste0("meanTrainFLEX1.4str_", k))
}

### flex2
##2 noeuds
for(j in strata_names){
  assign(paste0("ALL_FLEX2.2_", j),data.frame()) 
}

for(k in iterations){
  assign(paste0("meanTrainFLEX2.2str_", k) , read.csv(paste0(path0, "TRAIN/STRATA/FLEX2.2/mean_survival/",k,"mean_flex22_strates.csv"), sep = ";", row.names = 1) )
  
  for(j in strata_names){
    assign(paste0("ALL_FLEX2.2_",j) ,rbind(get(paste0("ALL_FLEX2.2_",j)), get(paste0("meanTrainFLEX2.2str_",k))[paste0("mean.flex2.2_", j),]  ) )
  }
  rm(list =paste0("meanTrainFLEX2.2str_", k))
}

##4 noeuds
for(j in strata_names){
  assign(paste0("ALL_FLEX2.4_", j),data.frame()) 
}

for(k in iterations){
  assign(paste0("meanTrainFLEX2.4str_", k) , read.csv(paste0(path0, "TRAIN/STRATA/FLEX2.4/mean_survival/",k,"mean_flex24_strates.csv"), sep = ";", row.names = 1) )
  
  for(j in strata_names){
    assign(paste0("ALL_FLEX2.4_",j) ,rbind(get(paste0("ALL_FLEX2.4_",j)), get(paste0("meanTrainFLEX2.4str_",k))[paste0("mean.flex2.4_", j),]  ) )
  }
  rm(list =paste0("meanTrainFLEX2.4str_", k))
}
### WG

for(j in strata_names){
  assign(paste0("ALL_WG_", j),data.frame()) 
}

for(k in iterations){
  assign(paste0("meanTrainWGstr_", k) , read.csv(paste0(path0, "TRAIN/STRATA/WG/mean_survival/",k,"mean_WG_strates.csv"), sep = ";", row.names = 1) )
  
  for(j in strata_names){
    assign(paste0("ALL_WG_",j) ,rbind(get(paste0("ALL_WG_",j)), get(paste0("meanTrainWGstr_",k))[paste0("mean.WG_", j),]  ) )
  }
  rm(list =paste0("meanTrainWGstr_", k))
}

### PP

for(j in strata_names){
  assign(paste0("ALL_PP_", j),data.frame()) 
}

for(k in iterations){
  assign(paste0("meanTrainPPstr_", k) , read.csv(paste0(path0, "TRAIN/STRATA/PP/",k,"PP_strates.csv"), sep = ";", row.names = 1) )
  
  for(j in strata_names){
    assign(paste0("ALL_PP_",j) ,rbind(get(paste0("ALL_PP_",j)), get(paste0("meanTrainPPstr_",k))[paste0("poharpred_", j),]  ) )
  }
  rm(list =paste0("meanTrainPPstr_", k))
}

#### Calcul des indicateurs
##PLANN
biais_PLANN_strates <- data.frame()
RMSE_PLANN_strates <- data.frame()
for(j in strata_names){
  assign(paste0("biais_PLANN_", j), 
         colMeans( get( paste0("ALL_PLANN_", j) ) - get(paste0("S_theo", j) )[col(get( paste0("ALL_PLANN_", j)) )] ) )
  #on concatene tous les résultats dans une dataframe pour ne pas saturer l'evt 
  biais_PLANN_strates <- rbind(biais_PLANN_strates, get(paste0("biais_PLANN_", j)) )
  rownames(biais_PLANN_strates)[dim(biais_PLANN_strates)[1]] <- paste0(j)
  rm(list = paste0("biais_PLANN_", j))
  
  assign(paste0("RMSE_PLANN_",  j) ,
         sqrt( colMeans( ( get( paste0("ALL_PLANN_", j) ) - 
                             get(paste0("S_theo", j) )[col(get( paste0("ALL_PLANN_", j)) )] )^2 ) ) )
  RMSE_PLANN_strates <- rbind(RMSE_PLANN_strates, get(paste0("RMSE_PLANN_", j)) )
  rownames(RMSE_PLANN_strates)[dim(RMSE_PLANN_strates)[1]] <- paste0(j)
  rm(list = paste0("RMSE_PLANN_", j))
}

colnames(biais_PLANN_strates) <- colnames(ALL_flex1.4)
colnames(RMSE_PLANN_strates) <- colnames(ALL_flex1.4)

# FLEX 1
#2 noeuds
biais_FLEX1.2_strates <- data.frame()
RMSE_FLEX1.2_strates <- data.frame()
for(j in strata_names){
  assign(paste0("biais_FLEX1.2_", j), 
         colMeans( get( paste0("ALL_FLEX1.2_", j) ) - get(paste0("S_theo", j) )[col(get( paste0("ALL_FLEX1.2_", j)) )] ) )
  
  #on concatene tous les résultats dans une dataframe pour ne pas saturer l'evt 
  biais_FLEX1.2_strates <- rbind(biais_FLEX1.2_strates, get(paste0("biais_FLEX1.2_", j)) )
  rownames(biais_FLEX1.2_strates)[dim(biais_FLEX1.2_strates)[1]] <- paste0(j)
  rm(list = paste0("biais_FLEX1.2_", j))
  
  assign(paste0("RMSE_FLEX1.2_",  j) ,
         sqrt( colMeans( ( get( paste0("ALL_FLEX1.2_", j) ) - 
                             get(paste0("S_theo", j) )[col(get( paste0("ALL_FLEX1.2_", j)) )] )^2 ) ) )
  
  RMSE_FLEX1.2_strates <- rbind(RMSE_FLEX1.2_strates, get(paste0("RMSE_FLEX1.2_", j)) )
  rownames(RMSE_FLEX1.2_strates)[dim(RMSE_FLEX1.2_strates)[1]] <- paste0(j)
  rm(list = paste0("RMSE_FLEX1.2_", j))
  
}

colnames(biais_FLEX1.2_strates) <- colnames(ALL_flex1.2)
colnames(RMSE_FLEX1.2_strates) <- colnames(ALL_flex1.2)

#4 noeuds
biais_FLEX1.4_strates <- data.frame()
RMSE_FLEX1.4_strates <- data.frame()
for(j in strata_names){
  assign(paste0("biais_FLEX1.4_", j), 
         colMeans( get( paste0("ALL_FLEX1.4_", j) ) - get(paste0("S_theo", j) )[col(get( paste0("ALL_FLEX1.4_", j)) )] ) )
  
  #on concatene tous les résultats dans une dataframe pour ne pas saturer l'evt 
  biais_FLEX1.4_strates <- rbind(biais_FLEX1.4_strates, get(paste0("biais_FLEX1.4_", j)) )
  rownames(biais_FLEX1.4_strates)[dim(biais_FLEX1.4_strates)[1]] <- paste0(j)
  rm(list = paste0("biais_FLEX1.4_", j))
  
  assign(paste0("RMSE_FLEX1.4_",  j) ,
         sqrt( colMeans( ( get( paste0("ALL_FLEX1.4_", j) ) - 
                             get(paste0("S_theo", j) )[col(get( paste0("ALL_FLEX1.4_", j)) )] )^2 ) ) )
  
  RMSE_FLEX1.4_strates <- rbind(RMSE_FLEX1.4_strates, get(paste0("RMSE_FLEX1.4_", j)) )
  rownames(RMSE_FLEX1.4_strates)[dim(RMSE_FLEX1.4_strates)[1]] <- paste0(j)
  rm(list = paste0("RMSE_FLEX1.4_", j))
  
}

colnames(biais_FLEX1.4_strates) <- colnames(ALL_flex1.4)
colnames(RMSE_FLEX1.4_strates) <- colnames(ALL_flex1.4)


# FLEX 2
##2 noeuds
biais_FLEX2.2_strates <- data.frame()
RMSE_FLEX2.2_strates <- data.frame()
for(j in strata_names){
  assign(paste0("biais_FLEX2.2_", j), 
         colMeans( get( paste0("ALL_FLEX2.2_", j) ) - get(paste0("S_theo", j) )[col(get( paste0("ALL_FLEX2.2_", j)) )] ) )
  
  biais_FLEX2.2_strates <- rbind(biais_FLEX2.2_strates, get(paste0("biais_FLEX2.2_", j)) )
  rownames(biais_FLEX2.2_strates)[dim(biais_FLEX2.2_strates)[1]] <- paste0(j)
  rm(list = paste0("biais_FLEX2.2_", j))
  
  assign(paste0("RMSE_FLEX2.2_",  j) ,
         sqrt( colMeans( ( get( paste0("ALL_FLEX2.2_", j) ) - 
                             get(paste0("S_theo", j) )[col(get( paste0("ALL_FLEX2.2_", j)) )] )^2 ) ) )
  
  RMSE_FLEX2.2_strates <- rbind(RMSE_FLEX2.2_strates, get(paste0("RMSE_FLEX2.2_", j)) )
  rownames(RMSE_FLEX2.2_strates)[dim(RMSE_FLEX2.2_strates)[1]] <- paste0(j)
  rm(list = paste0("RMSE_FLEX2.2_", j))
}
colnames(biais_FLEX2.2_strates) <- colnames(ALL_flex2.2)
colnames(RMSE_FLEX2.2_strates) <- colnames(ALL_flex2.2)

##4 noeuds
biais_FLEX2.4_strates <- data.frame()
RMSE_FLEX2.4_strates <- data.frame()
for(j in strata_names){
  assign(paste0("biais_FLEX2.4_", j), 
         colMeans( get( paste0("ALL_FLEX2.4_", j) ) - get(paste0("S_theo", j) )[col(get( paste0("ALL_FLEX2.4_", j)) )] ) )
  
  biais_FLEX2.4_strates <- rbind(biais_FLEX2.4_strates, get(paste0("biais_FLEX2.4_", j)) )
  rownames(biais_FLEX2.4_strates)[dim(biais_FLEX2.4_strates)[1]] <- paste0(j)
  rm(list = paste0("biais_FLEX2.4_", j))
  
  assign(paste0("RMSE_FLEX2.4_",  j) ,
         sqrt( colMeans( ( get( paste0("ALL_FLEX2.4_", j) ) - 
                             get(paste0("S_theo", j) )[col(get( paste0("ALL_FLEX2.4_", j)) )] )^2 ) ) )
  
  RMSE_FLEX2.4_strates <- rbind(RMSE_FLEX2.4_strates, get(paste0("RMSE_FLEX2.4_", j)) )
  rownames(RMSE_FLEX2.4_strates)[dim(RMSE_FLEX2.4_strates)[1]] <- paste0(j)
  rm(list = paste0("RMSE_FLEX2.4_", j))
}
colnames(biais_FLEX2.4_strates) <- colnames(ALL_flex2.4)
colnames(RMSE_FLEX2.4_strates) <- colnames(ALL_flex2.4)

# WG 
biais_WG_strates <- data.frame()
RMSE_WG_strates <- data.frame()
for(j in strata_names){
  assign(paste0("biais_WG_", j), 
         colMeans( get( paste0("ALL_WG_", j) ) - get(paste0("S_theo", j) )[col(get( paste0("ALL_WG_", j)) )] ) )
  
  biais_WG_strates <- rbind(biais_WG_strates, get(paste0("biais_WG_", j)) )
  rownames(biais_WG_strates)[dim(biais_WG_strates)[1]] <- paste0(j)
  rm(list = paste0("biais_WG_", j))
  
  assign(paste0("RMSE_WG_",  j) ,
         sqrt( colMeans( ( get( paste0("ALL_WG_", j) ) - 
                             get(paste0("S_theo", j) )[col(get( paste0("ALL_WG_", j)) )] )^2 ) ) )
  
  RMSE_WG_strates <- rbind(RMSE_WG_strates, get(paste0("RMSE_WG_", j)) )
  rownames(RMSE_WG_strates)[dim(RMSE_WG_strates)[1]] <- paste0(j)
  rm(list = paste0("RMSE_WG_", j))
}

colnames(biais_WG_strates) <- colnames(ALL_WG)
colnames(RMSE_WG_strates) <- colnames(ALL_WG)

# PP
biais_PP_strates <- data.frame()
RMSE_PP_strates <- data.frame()
for(j in strata_names){
  assign(paste0("biais_PP_", j), 
         colMeans( get( paste0("ALL_PP_", j) ) - get(paste0("S_theo", j) )[col(get( paste0("ALL_PP_", j)) )] ) )
  
  biais_PP_strates <- rbind(biais_PP_strates, get(paste0("biais_PP_", j)) )
  rownames(biais_PP_strates)[dim(biais_PP_strates)[1]] <- paste0(j)
  rm(list = paste0("biais_PP_", j))
  
  assign(paste0("RMSE_PP_",  j) ,
         sqrt( colMeans( ( get( paste0("ALL_PP_", j) ) - 
                             get(paste0("S_theo", j) )[col(get( paste0("ALL_PP_", j)) )] )^2 ) ) )
  
  RMSE_PP_strates <- rbind(RMSE_PP_strates, get(paste0("RMSE_PP_", j)) )
  rownames(RMSE_PP_strates)[dim(RMSE_PP_strates)[1]] <- paste0(j)
  rm(list = paste0("RMSE_PP_", j))
}

colnames(biais_PP_strates) <- colnames(ALL_PP)
colnames(RMSE_PP_strates) <- colnames(ALL_PP)

###concaten des S_theo, là aussi pour libérer de la place dans l'evt global
S_theo_strates <- data.frame()
for(j in strata_names){
  S_theo_strates <- rbind(S_theo_strates, get(paste0("S_theo",j)))
  rownames(S_theo_strates)[dim(S_theo_strates)[1]] <- paste0(j)
  rm(list = paste0("S_theo",j))
}
colnames(S_theo_strates) <- colnames(ALL_WG)

############################################ VALID #################################################

# WHOLE
###on récupère les valeurs théoriques moyennes Se(t) 
ALL_theoval <- data.frame()

for(k in iterations){
  assign(paste0("meanValidTheo_", k) , read.csv(paste0(path0, "VALID/WHOLE/THEO/mean_survival/",k,"_mean_theoval.csv"),
                                                sep = ";") )
  
  ALL_theoval <-rbind(ALL_theoval, get(paste0("meanValidTheo_", k)))
  rm(list =paste0("meanValidTheo_", k))
}

#moyenne théorique de la valur de survie
S_theoval <- colMeans(ALL_theoval)

############################################
#importation des données des différents modèles 
#########
##plann
ALL_PLANNval <- data.frame() 
for(k in iterations){
  assign(paste0("meanValidPLANN_", k) , read.csv(paste0(path0, "VALID/WHOLE/PLANN/mean_survival/",k,"_mean_PLANNval.csv"), sep = ";") )
  ALL_PLANNval <-rbind(ALL_PLANNval, get(paste0("meanValidPLANN_", k)))
  rm(list =paste0("meanValidPLANN_", k))
}

###flex1
##2 noeuds
ALL_flex1.2val <- data.frame()
for(k in iterations){
  assign(paste0("meanValidflex1.2_", k) , read.csv(paste0(path0, "VALID/WHOLE/FLEX1.2/mean_survival/",k,"_mean_flexval12.csv"), sep = ";") )
  ALL_flex1.2val <- rbind(ALL_flex1.2val, get(paste0("meanValidflex1.2_",k)))
  rm(list =paste0("meanValidflex1.2_", k))
}
##4 noeuds
ALL_flex1.4val <- data.frame()
for(k in iterations){
  assign(paste0("meanValidflex1.4_", k) , read.csv(paste0(path0, "VALID/WHOLE/FLEX1.4/mean_survival/",k,"_mean_flexval14.csv"), sep = ";") )
  ALL_flex1.4val <- rbind(ALL_flex1.4val, get(paste0("meanValidflex1.4_",k)))
  rm(list =paste0("meanValidflex1.4_", k))
}
### flex2

##2 noeuds
ALL_flex2.2val <- data.frame()
for(k in iterations){
  assign(paste0("meanValidflex2.2_", k) , read.csv(paste0(path0, "VALID/WHOLE/FLEX2.2/mean_survival/",k,"_mean_flexval22.csv"), sep = ";") )
  ALL_flex2.2val <- rbind(ALL_flex2.2val, get(paste0("meanValidflex2.2_",k)))
  rm(list =paste0("meanValidflex2.2_", k))
}

##4 noeuds
ALL_flex2.4val <- data.frame()
for(k in iterations){
  assign(paste0("meanValidflex2.4_", k) , read.csv(paste0(path0, "VALID/WHOLE/FLEX2.4/mean_survival/",k,"_mean_flexval24.csv"), sep = ";") )
  ALL_flex2.4val <- rbind(ALL_flex2.4val, get(paste0("meanValidflex2.4_",k)))
  rm(list =paste0("meanValidflex2.4_", k))
}
### WG
ALL_WGval <- data.frame()
for(k in iterations){
  assign(paste0("meanValidWG_", k) , read.csv(paste0(path0, "VALID/WHOLE/WG/mean_survival/",k,"_mean_flexWGval.csv"), sep = ";") )
  ALL_WGval <- rbind(ALL_WGval, get(paste0("meanValidWG_",k))) 
  rm(list =paste0("meanValidWG_", k))
}

##PP
ALL_PPval <- data.frame()
for(k in iterations){
  assign(paste0("meanValidPP_", k) , read.csv(paste0(path0, "VALID/WHOLE/PP/",k,"PPval.csv"), sep = ";") )
  ALL_PPval <- rbind(ALL_PPval, get(paste0("meanValidPP_",k))) 
  rm(list =paste0("meanValidPP_", k))
}
ALL_PPval <- as.data.frame(matrix(ALL_PPval$x, ncol = 4, byrow = TRUE))

biais_PLANNval <- colMeans(ALL_PLANNval[,-1] - S_theoval[col( ALL_PLANNval[,-1] )] )
RMSE_PLANNval <- sqrt( colMeans((ALL_PLANNval[,-1] - S_theoval[col( ALL_PLANNval[,-1] )])^2) )

biais_flex1.2val <- colMeans( ALL_flex1.2val - S_theoval[col( ALL_flex1.2val)] )
RMSE_flex1.2val <- sqrt( colMeans( (ALL_flex1.2val - S_theoval[col( ALL_flex1.2val)])^2 ) )

biais_flex1.4val <- colMeans( ALL_flex1.4val - S_theoval[col( ALL_flex1.4val)] )
RMSE_flex1.4val <- sqrt( colMeans( (ALL_flex1.4val - S_theoval[col( ALL_flex1.4val)])^2 ) )

biais_flex2.2val <- colMeans( ALL_flex2.2val - S_theoval[col( ALL_flex2.2val)] )
RMSE_flex2.2val <- sqrt( colMeans( (ALL_flex2.2val - S_theoval[col( ALL_flex2.2val)])^2 ) )

biais_flex2.4val <- colMeans( ALL_flex2.4val - S_theoval[col( ALL_flex2.4val)] )
RMSE_flex2.4val <- sqrt( colMeans( (ALL_flex2.4val - S_theoval[col( ALL_flex2.4val)])^2 ) )

biais_WGval <- colMeans( ALL_WGval - S_theoval[col( ALL_WGval)] )
RMSE_WGval <- sqrt( colMeans( (ALL_WGval - S_theoval[col( ALL_WGval)])^2 ) )

biais_PPval <- colMeans( ALL_PP - S_theoval[col( ALL_PPval)] )
RMSE_PPval <- sqrt( colMeans( (ALL_PP - S_theoval[col( ALL_PPval)])^2 ) )
#STRATES

###on récupère les valeurs théoriques moyennes Se(t) 
for(j in strata_names){
  assign(paste0("ALL_",j, "_theoval"),  data.frame()) 
}

for(k in iterations){
  assign(paste0("meanValidTheoS_", k) , read.csv(paste0(path0, "VALID/STRATA/THEO/mean_survival/",k,"mean_theoval_strates.csv"),
                                                 sep = ";", row.names = 1)[,-1] ) ### [,-1] car quand j'ai enregistré les résultats, j'ai dupliué la premiere colonne donc il faut l'enlever
  for(j in strata_names){
    
    current_row <- get(paste0("meanValidTheoS_", k))[paste0("mean.theo_", j), , drop = FALSE]
    rownames(current_row) <- paste0("mean.theoval_", j, k)
    assign(paste0("ALL_", j, "_theoval"), rbind(get(paste0("ALL_", j, "_theoval")), current_row))
    
  }
  
  rm(list =paste0("meanValidTheoS_", k))
}
rm(current_row)
#moyenne théorique de la valeur de survie
for(j in strata_names){
  assign(paste0("S_theoval", j ),colMeans(get(paste0("ALL_",j,"_theoval") ) ) )
}
############################################
#importation des données des différents modèles 
############################################

# PLANN

for(j in strata_names){
  assign(paste0("ALL_PLANNval_", j),data.frame()) 
}

for(k in iterations){
  assign(paste0("meanValidPLANNstr_", k) , read.csv(paste0(path0, "VALID/STRATA/PLANN/mean_survival/",k,"mean_plannval_strates.csv"), sep = ";", row.names = 1) )
  
  for(j in strata_names){
    assign(paste0("ALL_PLANNval_",j) ,rbind(get(paste0("ALL_PLANNval_",j)),get(paste0("meanValidPLANNstr_",k))[paste0("mean.plannval_", j),-1]  ) )
  }
  rm(list =paste0("meanValidPLANNstr_", k))
}

### flex1
##2 noeuds
for(j in strata_names){
  assign(paste0("ALL_FLEX1.2val_", j),data.frame()) 
}

for(k in iterations){
  assign(paste0("meanValidFLEX1.2str_", k) , read.csv(paste0(path0, "VALID/STRATA/FLEX1.2/mean_survival/",k,"mean_flexval12_strates.csv"), sep = ";", row.names = 1) )
  
  for(j in strata_names){
    assign(paste0("ALL_FLEX1.2val_",j) ,rbind(get(paste0("ALL_FLEX1.2val_",j)), get(paste0("meanValidFLEX1.2str_",k))[paste0("mean.flexval1.2_", j),]  ) )
  }
  rm(list =paste0("meanValidFLEX1.2str_", k))
}

##4 noeuds
for(j in strata_names){
  assign(paste0("ALL_FLEX1.4val_", j),data.frame()) 
}

for(k in iterations){
  assign(paste0("meanValidFLEX1.4str_", k) , read.csv(paste0(path0, "VALID/STRATA/FLEX1.4/mean_survival/",k,"mean_flexval14_strates.csv"), sep = ";", row.names = 1) )
  
  for(j in strata_names){
    assign(paste0("ALL_FLEX1.4val_",j) ,rbind(get(paste0("ALL_FLEX1.4val_",j)), get(paste0("meanValidFLEX1.4str_",k))[paste0("mean.flexval1.4_", j),]  ) )
  }
  rm(list =paste0("meanValidFLEX1.4str_", k))
}
### flex2

##2 noeuds
for(j in strata_names){
  assign(paste0("ALL_FLEX2.2val_", j),data.frame()) 
}

for(k in iterations){
  assign(paste0("meanValidFLEX2.2str_", k) , read.csv(paste0(path0, "VALID/STRATA/FLEX2.2/mean_survival/",k,"mean_flexval22_strates.csv"), sep = ";", row.names = 1) )
  
  for(j in strata_names){
    assign(paste0("ALL_FLEX2.2val_",j) ,rbind(get(paste0("ALL_FLEX2.2val_",j)), get(paste0("meanValidFLEX2.2str_",k))[paste0("mean.flexval2.2_", j),]  ) )
  }
  rm(list =paste0("meanValidFLEX2.2str_", k))
}

##4 noeuds
for(j in strata_names){
  assign(paste0("ALL_FLEX2.4val_", j),data.frame()) 
}

for(k in iterations){
  assign(paste0("meanValidFLEX2.4str_", k) , read.csv(paste0(path0, "VALID/STRATA/FLEX2.4/mean_survival/",k,"mean_flexval24_strates.csv"), sep = ";", row.names = 1) )
  
  for(j in strata_names){
    assign(paste0("ALL_FLEX2.4val_",j) ,rbind(get(paste0("ALL_FLEX2.4val_",j)), get(paste0("meanValidFLEX2.4str_",k))[paste0("mean.flexval2.4_", j),]  ) )
  }
  rm(list =paste0("meanValidFLEX2.4str_", k))
}
### WG

for(j in strata_names){
  assign(paste0("ALL_WGval_", j),data.frame()) 
}

for(k in iterations){
  assign(paste0("meanValidWGstr_", k) , read.csv(paste0(path0, "VALID/STRATA/WG/mean_survival/",k,"mean_WGval_strates.csv"), sep = ";", row.names = 1) )
  
  for(j in strata_names){
    assign(paste0("ALL_WGval_",j) ,rbind(get(paste0("ALL_WGval_",j)), get(paste0("meanValidWGstr_",k))[paste0("mean.WGval_", j),]  ) )
  }
  rm(list =paste0("meanValidWGstr_", k))
}

### PP 

for(j in strata_names){
  assign(paste0("ALL_PPval_", j),data.frame()) 
}

for(k in iterations){
  assign(paste0("meanValidPPstr_", k) , read.csv(paste0(path0, "VALID/STRATA/PP/",k,"PPval_strates.csv"), sep = ";", row.names = 1) )
  
  for(j in strata_names){
    assign(paste0("ALL_PPval_",j) ,rbind(get(paste0("ALL_PPval_",j)), get(paste0("meanValidPPstr_",k))[paste0("poharpredval_", j),]  ) )
  }
  rm(list =paste0("meanValidPPstr_", k))
}
#### Calcul des indicateurs
##PLANN
biais_PLANNval_strates <- data.frame()
RMSE_PLANNval_strates <- data.frame()
for(j in strata_names){
  assign(paste0("biais_PLANNval_", j), 
         colMeans( get( paste0("ALL_PLANNval_", j) ) - get(paste0("S_theoval", j) )[col(get( paste0("ALL_PLANNval_", j)) )] ) )
  
  biais_PLANNval_strates <- rbind(biais_PLANNval_strates, get(paste0("biais_PLANNval_", j)) )
  rownames(biais_PLANNval_strates)[dim(biais_PLANNval_strates)[1]] <- paste0(j)
  rm(list = paste0("biais_PLANNval_", j))
  
  assign(paste0("RMSE_PLANNval_",  j) ,
         sqrt( colMeans( ( get( paste0("ALL_PLANNval_", j) ) - 
                             get(paste0("S_theoval", j) )[col(get( paste0("ALL_PLANNval_", j)) )] )^2 ) ) )
  
  RMSE_PLANNval_strates <- rbind(RMSE_PLANNval_strates, get(paste0("RMSE_PLANNval_", j)) )
  rownames(RMSE_PLANNval_strates)[dim(RMSE_PLANNval_strates)[1]] <- paste0(j)
  rm(list = paste0("RMSE_PLANNval_", j))
  
}

colnames(biais_PLANNval_strates) <- colnames(ALL_flex1.2)
colnames(RMSE_PLANNval_strates) <- colnames(ALL_flex1.2)

# FLEX 1

##2 noeuds
biais_FLEX1.2val_strates <- data.frame()
RMSE_FLEX1.2val_strates <- data.frame()
for(j in strata_names){
  assign(paste0("biais_FLEX1.2val_", j), 
         colMeans( get( paste0("ALL_FLEX1.2val_", j) ) - get(paste0("S_theoval", j) )[col(get( paste0("ALL_FLEX1.2val_", j)) )] ) )
  
  biais_FLEX1.2val_strates <- rbind(biais_FLEX1.2val_strates, get(paste0("biais_FLEX1.2val_", j)) )
  rownames(biais_FLEX1.2val_strates)[dim(biais_FLEX1.2val_strates)[1]] <- paste0(j)
  rm(list = paste0("biais_FLEX1.2val_", j))
  
  assign(paste0("RMSE_FLEX1.2val_",  j) ,
         sqrt( colMeans( ( get( paste0("ALL_FLEX1.2val_", j) ) - 
                             get(paste0("S_theoval", j) )[col(get( paste0("ALL_FLEX1.2val_", j)) )] )^2 ) ) )
  
  RMSE_FLEX1.2val_strates <- rbind(RMSE_FLEX1.2val_strates, get(paste0("RMSE_FLEX1.2val_", j)) )
  rownames(RMSE_FLEX1.2val_strates)[dim(RMSE_FLEX1.2val_strates)[1]] <- paste0(j)
  rm(list = paste0("RMSE_FLEX1.2val_", j))
}

colnames(biais_FLEX1.2val_strates) <- colnames(ALL_flex1.2)
colnames(RMSE_FLEX1.2val_strates) <- colnames(ALL_flex1.2)

##4 noeuds
biais_FLEX1.4val_strates <- data.frame()
RMSE_FLEX1.4val_strates <- data.frame()
for(j in strata_names){
  assign(paste0("biais_FLEX1.4val_", j), 
         colMeans( get( paste0("ALL_FLEX1.4val_", j) ) - get(paste0("S_theoval", j) )[col(get( paste0("ALL_FLEX1.4val_", j)) )] ) )
  
  biais_FLEX1.4val_strates <- rbind(biais_FLEX1.4val_strates, get(paste0("biais_FLEX1.4val_", j)) )
  rownames(biais_FLEX1.4val_strates)[dim(biais_FLEX1.4val_strates)[1]] <- paste0(j)
  rm(list = paste0("biais_FLEX1.4val_", j))
  
  assign(paste0("RMSE_FLEX1.4val_",  j) ,
         sqrt( colMeans( ( get( paste0("ALL_FLEX1.4val_", j) ) - 
                             get(paste0("S_theoval", j) )[col(get( paste0("ALL_FLEX1.4val_", j)) )] )^2 ) ) )
  
  RMSE_FLEX1.4val_strates <- rbind(RMSE_FLEX1.4val_strates, get(paste0("RMSE_FLEX1.4val_", j)) )
  rownames(RMSE_FLEX1.4val_strates)[dim(RMSE_FLEX1.4val_strates)[1]] <- paste0(j)
  rm(list = paste0("RMSE_FLEX1.4val_", j))
}

colnames(biais_FLEX1.4val_strates) <- colnames(ALL_flex1.4)
colnames(RMSE_FLEX1.4val_strates) <- colnames(ALL_flex1.4)

# FLEX 2

##2 noeuds
biais_FLEX2.2val_strates <- data.frame()
RMSE_FLEX2.2val_strates <- data.frame()

for(j in strata_names){
  assign(paste0("biais_FLEX2.2val_", j), 
         colMeans( get( paste0("ALL_FLEX2.2val_", j) ) - get(paste0("S_theoval", j) )[col(get( paste0("ALL_FLEX2.2val_", j)) )] ) )
  
  biais_FLEX2.2val_strates <- rbind(biais_FLEX2.2val_strates, get(paste0("biais_FLEX2.2val_", j)) )
  rownames(biais_FLEX2.2val_strates)[dim(biais_FLEX2.2val_strates)[1]] <- paste0(j)
  rm(list = paste0("biais_FLEX2.2val_", j))
  
  assign(paste0("RMSE_FLEX2.2val_",  j) ,
         sqrt( colMeans( ( get( paste0("ALL_FLEX2.2val_", j) ) - 
                             get(paste0("S_theoval", j) )[col(get( paste0("ALL_FLEX2.2val_", j)) )] )^2 ) ) )
  
  RMSE_FLEX2.2val_strates <- rbind(RMSE_FLEX2.2val_strates, get(paste0("RMSE_FLEX2.2val_", j)) )
  rownames(RMSE_FLEX2.2val_strates)[dim(RMSE_FLEX2.2val_strates)[1]] <- paste0(j)
  rm(list = paste0("RMSE_FLEX2.2val_", j))
  
}

colnames(biais_FLEX2.2val_strates) <- colnames(ALL_flex1.2)
colnames(RMSE_FLEX2.2val_strates) <- colnames(ALL_flex1.2)

##4 noeuds
biais_FLEX2.4val_strates <- data.frame()
RMSE_FLEX2.4val_strates <- data.frame()

for(j in strata_names){
  assign(paste0("biais_FLEX2.4val_", j), 
         colMeans( get( paste0("ALL_FLEX2.4val_", j) ) - get(paste0("S_theoval", j) )[col(get( paste0("ALL_FLEX2.4val_", j)) )] ) )
  
  biais_FLEX2.4val_strates <- rbind(biais_FLEX2.4val_strates, get(paste0("biais_FLEX2.4val_", j)) )
  rownames(biais_FLEX2.4val_strates)[dim(biais_FLEX2.4val_strates)[1]] <- paste0(j)
  rm(list = paste0("biais_FLEX2.4val_", j))
  
  assign(paste0("RMSE_FLEX2.4val_",  j) ,
         sqrt( colMeans( ( get( paste0("ALL_FLEX2.4val_", j) ) - 
                             get(paste0("S_theoval", j) )[col(get( paste0("ALL_FLEX2.4val_", j)) )] )^2 ) ) )
  
  RMSE_FLEX2.4val_strates <- rbind(RMSE_FLEX2.4val_strates, get(paste0("RMSE_FLEX2.4val_", j)) )
  rownames(RMSE_FLEX2.4val_strates)[dim(RMSE_FLEX2.4val_strates)[1]] <- paste0(j)
  rm(list = paste0("RMSE_FLEX2.4val_", j))
  
}

colnames(biais_FLEX2.4val_strates) <- colnames(ALL_flex1.2)
colnames(RMSE_FLEX2.4val_strates) <- colnames(ALL_flex1.2)
# WG 

biais_WGval_strates <- data.frame()
RMSE_WGval_strates <- data.frame()

for(j in strata_names){
  assign(paste0("biais_WGval_", j), 
         colMeans( get( paste0("ALL_WGval_", j) ) - get(paste0("S_theoval", j) )[col(get( paste0("ALL_WGval_", j)) )] ) )
  
  biais_WGval_strates <- rbind(biais_WGval_strates, get(paste0("biais_WGval_", j)) )
  rownames(biais_WGval_strates)[dim(biais_WGval_strates)[1]] <- paste0(j)
  rm(list = paste0("biais_WGval_", j))
  
  assign(paste0("RMSE_WGval_",  j) ,
         sqrt( colMeans( ( get( paste0("ALL_WGval_", j) ) - 
                             get(paste0("S_theoval", j) )[col(get( paste0("ALL_WGval_", j)) )] )^2 ) ) )
  
  RMSE_WGval_strates <- rbind(RMSE_WGval_strates, get(paste0("RMSE_WGval_", j)) )
  rownames(RMSE_WGval_strates)[dim(RMSE_WGval_strates)[1]] <- paste0(j)
  rm(list = paste0("RMSE_WGval_", j))
  
}

colnames(biais_WGval_strates) <- colnames(ALL_flex1.2)
colnames(RMSE_WGval_strates) <- colnames(ALL_flex1.2)

# PP ############################################## A CHANGER 
biais_PPval_strates <- data.frame()
RMSE_PPval_strates <- data.frame()
for(j in strata_names){
  assign(paste0("biais_PPval_", j), 
         colMeans( get( paste0("ALL_PP_", j) ) - get(paste0("S_theoval", j) )[col(get( paste0("ALL_PP_", j)) )] ) )
  
  biais_PPval_strates <- rbind(biais_PPval_strates, get(paste0("biais_PPval_", j)) )
  rownames(biais_PPval_strates)[dim(biais_PPval_strates)[1]] <- paste0(j)
  rm(list = paste0("biais_PPval_", j))
  
  assign(paste0("RMSE_PPval_",  j) ,
         sqrt( colMeans( ( get( paste0("ALL_PP_", j) ) - 
                             get(paste0("S_theoval", j) )[col(get( paste0("ALL_PP_", j)) )] )^2 ) ) )
  
  RMSE_PPval_strates <- rbind(RMSE_PPval_strates, get(paste0("RMSE_PPval_", j)) )
  rownames(RMSE_PPval_strates)[dim(RMSE_PPval_strates)[1]] <- paste0(j)
  rm(list = paste0("RMSE_PPval_", j))
}

colnames(biais_PPval_strates) <- colnames(ALL_PPval)
colnames(RMSE_PPval_strates) <- colnames(ALL_PPval)

###concaten des S_theo, là aussi pour libérer de la place dans l'evt global
S_theoval_strates <- data.frame()
for(j in strata_names){
  S_theoval_strates <- rbind(S_theoval_strates, get(paste0("S_theoval",j)))
  rownames(S_theoval_strates)[dim(S_theoval_strates)[1]] <- paste0(j)
  rm(list = paste0("S_theoval",j))
}
colnames(S_theoval_strates) <- colnames(ALL_WG)


## mise en forme de l'evt global 
##############################################
##concatenation STHEO/ STHEOSTRATES
##TRAIN
S_theo <- rbind(S_theo, S_theo_strates)
rownames(S_theo)[1] <- "ALL" 
rm(S_theo_strates)
##VALID

S_theoval <- rbind(S_theoval, S_theoval_strates)
rownames(S_theoval)[1] <- "ALL" 
rm(S_theoval_strates)
## concatenation des résultats sur la base totale avec ceux fait en strates 
##TRAIN
##PLANN
biais_PLANN <- rbind(biais_PLANN, biais_PLANN_strates)
RMSE_PLANN <- rbind(RMSE_PLANN, RMSE_PLANN_strates)
rownames(biais_PLANN)[1] <- "ALL"
rownames(RMSE_PLANN)[1] <- "ALL"
rm(biais_PLANN_strates)
rm(RMSE_PLANN_strates)
##FLEX1
##2noeuds
biais_FLEX1.2 <- rbind(biais_flex1.2, biais_FLEX1.2_strates)
RMSE_FLEX1.2 <- rbind(RMSE_flex1.2, RMSE_FLEX1.2_strates)
rownames(biais_FLEX1.2)[1] <- "ALL"
rownames(RMSE_FLEX1.2)[1] <- "ALL"
rm(biais_flex1.2, biais_FLEX1.2_strates)
rm(RMSE_flex1.2, RMSE_FLEX1.2_strates)

##4noeuds
biais_FLEX1.4 <- rbind(biais_flex1.4, biais_FLEX1.4_strates)
RMSE_FLEX1.4 <- rbind(RMSE_flex1.4, RMSE_FLEX1.4_strates)
rownames(biais_FLEX1.4)[1] <- "ALL"
rownames(RMSE_FLEX1.4)[1] <- "ALL"
rm(biais_flex1.4, biais_FLEX1.4_strates)
rm(RMSE_flex1.4, RMSE_FLEX1.4_strates)

##FLEX2
##2 noeuds
biais_FLEX2.2 <- rbind(biais_flex2.2, biais_FLEX2.2_strates)
RMSE_FLEX2.2 <- rbind(RMSE_flex2.2, RMSE_FLEX2.2_strates)
rownames(biais_FLEX2.2)[1] <- "ALL"
rownames(RMSE_FLEX2.2)[1] <- "ALL"
rm(biais_flex2.2, biais_FLEX2.2_strates)
rm(RMSE_flex2.2, RMSE_FLEX2.2_strates)
##4 noeuds
biais_FLEX2.4 <- rbind(biais_flex2.4, biais_FLEX2.4_strates)
RMSE_FLEX2.4 <- rbind(RMSE_flex2.4, RMSE_FLEX2.4_strates)
rownames(biais_FLEX2.4)[1] <- "ALL"
rownames(RMSE_FLEX2.4)[1] <- "ALL"
rm(biais_flex2.4, biais_FLEX2.4_strates)
rm(RMSE_flex2.4, RMSE_FLEX2.4_strates)
##WG

biais_WG <- rbind(biais_WG, biais_WG_strates)
RMSE_WG <- rbind(RMSE_WG, RMSE_WG_strates)
rownames(biais_WG)[1] <- "ALL"
rownames(RMSE_WG)[1] <- "ALL"
rm(biais_WG_strates)
rm(RMSE_WG_strates)

##PP

biais_PP <- rbind(biais_PP, biais_PP_strates)
RMSE_PP <- rbind(RMSE_PP, RMSE_PP_strates)
rownames(biais_PP)[1] <- "ALL"
rownames(RMSE_PP)[1] <- "ALL"
rm(biais_PP_strates)
rm(RMSE_PP_strates)

######VALID

biais_PLANNval <- rbind(biais_PLANNval, biais_PLANNval_strates)
RMSE_PLANNval <- rbind(RMSE_PLANNval, RMSE_PLANNval_strates)
rownames(biais_PLANNval)[1] <- "ALL"
rownames(RMSE_PLANNval)[1] <- "ALL"
rm(biais_PLANNval_strates)
rm(RMSE_PLANNval_strates)
##FLEX1
##2 noeuds
biais_FLEX1.2val <- rbind(biais_flex1.2val, biais_FLEX1.2val_strates)
RMSE_FLEX1.2val <- rbind(RMSE_flex1.2val, RMSE_FLEX1.2val_strates)
rownames(biais_FLEX1.2val)[1] <- "ALL"
rownames(RMSE_FLEX1.2val)[1] <- "ALL"
rm(biais_flex1.2val, biais_FLEX1.2val_strates)
rm(RMSE_flex1.2val, RMSE_FLEX1.2val_strates)
##4 noeuds
biais_FLEX1.4val <- rbind(biais_flex1.4val, biais_FLEX1.4val_strates)
RMSE_FLEX1.4val <- rbind(RMSE_flex1.4val, RMSE_FLEX1.4val_strates)
rownames(biais_FLEX1.4val)[1] <- "ALL"
rownames(RMSE_FLEX1.4val)[1] <- "ALL"
rm(biais_flex1.4val, biais_FLEX1.4val_strates)
rm(RMSE_flex1.4val, RMSE_FLEX1.4val_strates)
##FLEX2
##2 noeuds
biais_FLEX2.2val <- rbind(biais_flex2.2val, biais_FLEX2.2val_strates)
RMSE_FLEX2.2val <- rbind(RMSE_flex2.2val, RMSE_FLEX2.2val_strates)
rownames(biais_FLEX2.2val)[1] <- "ALL"
rownames(RMSE_FLEX2.2val)[1] <- "ALL"
rm(biais_flex2.2val, biais_FLEX2.2val_strates)
rm(RMSE_flex2.2val, RMSE_FLEX2.2val_strates)
##4 noeuds
biais_FLEX2.4val <- rbind(biais_flex2.4val, biais_FLEX2.4val_strates)
RMSE_FLEX2.4val <- rbind(RMSE_flex2.4val, RMSE_FLEX2.4val_strates)
rownames(biais_FLEX2.4val)[1] <- "ALL"
rownames(RMSE_FLEX2.4val)[1] <- "ALL"
rm(biais_flex2.4val, biais_FLEX2.4val_strates)
rm(RMSE_flex2.4val, RMSE_FLEX2.4val_strates)
##WG

biais_WGval <- rbind(biais_WGval, biais_WGval_strates)
RMSE_WGval <- rbind(RMSE_WGval, RMSE_WGval_strates)
rownames(biais_WGval)[1] <- "ALL"
rownames(RMSE_WGval)[1] <- "ALL"
rm(biais_WGval_strates)
rm(RMSE_WGval_strates)

##PP

biais_PPval <- rbind(biais_PPval, biais_PPval_strates)
RMSE_PPval <- rbind(RMSE_PPval, RMSE_PPval_strates)
rownames(biais_PPval)[1] <- "ALL"
rownames(RMSE_PPval)[1] <- "ALL"
rm(biais_PPval_strates)
rm(RMSE_PPval_strates)

### les dataframes des moyennes sont enregistrées sous forme de liste par rapport au modèle utilisé
########
##TRAIN
THEO_means_names <- c()
for(j in strata_names){
  THEO_means_names <- c(THEO_means_names, paste0("ALL_",j,"_theo"))
}
THEO_means <- c(mget("ALL_theo"), mget(THEO_means_names))
rm(ALL_theo)
rm(THEO_means_names)
for(j in strata_names){
  rm(list = paste0("ALL_",j,"_theo") )
}

##PLANN
PLANN_means_names <- c()
for(j in strata_names){
  PLANN_means_names <- c(PLANN_means_names, paste0("ALL_PLANN_",j))
}
PLANN_means <- c(mget("ALL_PLANN"), mget(PLANN_means_names))
rm(ALL_PLANN)
rm(PLANN_means_names)
for(j in strata_names){
  rm(list = paste0("ALL_PLANN_",j) )
}

##FLEX1
##2 noeuds
FLEX1.2_means_names <- c()
for(j in strata_names){
  FLEX1.2_means_names <- c(FLEX1.2_means_names, paste0("ALL_FLEX1.2_",j))
}
FLEX1.2_means <- c(mget("ALL_flex1.2"), mget(FLEX1.2_means_names))
rm(ALL_flex1.2)
rm(FLEX1.2_means_names)
for(j in strata_names){
  rm(list = paste0("ALL_FLEX1.2_",j) )
}

##4 noeuds
FLEX1.4_means_names <- c()
for(j in strata_names){
  FLEX1.4_means_names <- c(FLEX1.4_means_names, paste0("ALL_FLEX1.4_",j))
}
FLEX1.4_means <- c(mget("ALL_flex1.4"), mget(FLEX1.4_means_names))
rm(ALL_flex1.4)
rm(FLEX1.4_means_names)
for(j in strata_names){
  rm(list = paste0("ALL_FLEX1.4_",j) )
}

##FLEX2
##2 noeuds
FLEX2.2_means_names <- c()
for(j in strata_names){
  FLEX2.2_means_names <- c(FLEX2.2_means_names, paste0("ALL_FLEX2.2_",j))
}
FLEX2.2_means <- c(mget("ALL_flex2.2"), mget(FLEX2.2_means_names))
rm(ALL_flex2.2)
rm(FLEX2.2_means_names)
for(j in strata_names){
  rm(list = paste0("ALL_FLEX2.2_",j) )
}
##4 noeuds
FLEX2.4_means_names <- c()
for(j in strata_names){
  FLEX2.4_means_names <- c(FLEX2.4_means_names, paste0("ALL_FLEX2.4_",j))
}
FLEX2.4_means <- c(mget("ALL_flex2.4"), mget(FLEX2.4_means_names))
rm(ALL_flex2.4)
rm(FLEX2.4_means_names)
for(j in strata_names){
  rm(list = paste0("ALL_FLEX2.4_",j) )
}
##WG
WG_means_names <- c()
for(j in strata_names){
  WG_means_names <- c(WG_means_names, paste0("ALL_WG_",j))
}
WG_means <- c(mget("ALL_WG"), mget(WG_means_names))
rm(ALL_WG)
rm(WG_means_names)
for(j in strata_names){
  rm(list = paste0("ALL_WG_",j) )
}

PP_means_names <- c()
for(j in strata_names){
  PP_means_names <- c(PP_means_names, paste0("ALL_PP_",j))
}
PP_means <- c(mget("ALL_PP"), mget(PP_means_names))
rm(ALL_PP)
rm(PP_means_names)
for(j in strata_names){
  rm(list = paste0("ALL_PP_",j) )
}
############VALID
THEOval_means_names <- c()
for(j in strata_names){
  THEOval_means_names <- c(THEOval_means_names, paste0("ALL_",j,"_theoval"))
}
THEOval_means <- c(mget("ALL_theoval"), mget(THEOval_means_names))
rm(ALL_theoval)
rm(THEOval_means_names)
for(j in strata_names){
  rm(list = paste0("ALL_",j,"_theoval") )
}

##PLANN
PLANNval_means_names <- c()
for(j in strata_names){
  PLANNval_means_names <- c(PLANNval_means_names, paste0("ALL_PLANNval_",j))
}
PLANNval_means <- c(mget("ALL_PLANNval"), mget(PLANNval_means_names))
rm(ALL_PLANNval)
rm(PLANNval_means_names)
for(j in strata_names){
  rm(list = paste0("ALL_PLANNval_",j) )
}

##FLEX1
##2 noeuds
FLEX1.2val_means_names <- c()
for(j in strata_names){
  FLEX1.2val_means_names <- c(FLEX1.2val_means_names, paste0("ALL_FLEX1.2val_",j))
}
FLEX1.2val_means <- c(mget("ALL_flex1.2val"), mget(FLEX1.2val_means_names))
rm(ALL_flex1.2val)
rm(FLEX1.2val_means_names)
for(j in strata_names){
  rm(list = paste0("ALL_FLEX1.2val_",j) )
}

##4 noeuds
FLEX1.4val_means_names <- c()
for(j in strata_names){
  FLEX1.4val_means_names <- c(FLEX1.4val_means_names, paste0("ALL_FLEX1.4val_",j))
}
FLEX1.4val_means <- c(mget("ALL_flex1.4val"), mget(FLEX1.4val_means_names))
rm(ALL_flex1.4val)
rm(FLEX1.4val_means_names)
for(j in strata_names){
  rm(list = paste0("ALL_FLEX1.4val_",j) )
}
##FLEX2
##2 noeuds

FLEX2.2val_means_names <- c()
for(j in strata_names){
  FLEX2.2val_means_names <- c(FLEX2.2val_means_names, paste0("ALL_FLEX2.2val_",j))
}
FLEX2.2val_means <- c(mget("ALL_flex2.2val"), mget(FLEX2.2val_means_names))
rm(ALL_flex2.2val)
rm(FLEX2.2val_means_names)
for(j in strata_names){
  rm(list = paste0("ALL_FLEX2.2val_",j) )
}

##4 noeuds

FLEX2.4val_means_names <- c()
for(j in strata_names){
  FLEX2.4val_means_names <- c(FLEX2.4val_means_names, paste0("ALL_FLEX2.4val_",j))
}
FLEX2.4val_means <- c(mget("ALL_flex2.4val"), mget(FLEX2.4val_means_names))
rm(ALL_flex2.4val)
rm(FLEX2.4val_means_names)
for(j in strata_names){
  rm(list = paste0("ALL_FLEX2.4val_",j) )
}
##WG
WGval_means_names <- c()
for(j in strata_names){
  WGval_means_names <- c(WGval_means_names, paste0("ALL_WGval_",j))
}
WGval_means <- c(mget("ALL_WGval"), mget(WGval_means_names))
rm(ALL_WGval)
rm(WGval_means_names)
for(j in strata_names){
  rm(list = paste0("ALL_WGval_",j) )
}

##PP
PPval_means_names <- c()
for(j in strata_names){
  PPval_means_names <- c(PPval_means_names, paste0("ALL_PPval_",j))
}
PPval_means <- c(mget("ALL_PPval"), mget(PPval_means_names))
rm(ALL_PPval)
rm(PPval_means_names)
for(j in strata_names){
  rm(list = paste0("ALL_PPval_",j) )
}
##     1000ROCSTA
################################################### ROC.net ##############################################################

library(RISCA)
##choix de l'espacement des cutoffs 
# c.off <- seq(0,1, by = 0.2)
# sort(unique(1-ind_estimP[, l + 1]))
# unique( quantile( 1-ind_estimP[, l + 1],  probs= seq(0,1, by = 0.01)) ) #pas en dessous de 20 points 0.025
####TRAIN
### WHOLE

##importation des bases de données / survies individuelles prédites par les modèles 
for(k in iterations){
  ##bases de données 
  assign(paste0("DATATrain_", k) , read.csv(paste0(path0, "DATAFRAMES/TRAIN/WHOLE/",k,"_dftrain.csv"),
                                            sep = ";") )
  
  ##survies individuelles 
  #plann
  assign(paste0("indTrainPLANN_", k) , read.csv(paste0(path0, "TRAIN/WHOLE/PLANN/individual_survival/",k,"_ind_PLANN.csv"),
                                                sep = ";") )
  #flex1
  assign(paste0("indTrainFLEX1.2_", k) , read.csv(paste0(path0, "TRAIN/WHOLE/FLEX1.2/individual_survival/",k,"_ind_flex12.csv"),
                                                  sep = ";") )
  assign(paste0("indTrainFLEX1.4_", k) , read.csv(paste0(path0, "TRAIN/WHOLE/FLEX1.4/individual_survival/",k,"_ind_flex14.csv"),
                                                  sep = ";") )
  #flex2
  assign(paste0("indTrainFLEX2.2_", k) , read.csv(paste0(path0, "TRAIN/WHOLE/FLEX2.2/individual_survival/",k,"_ind_flex22.csv"),
                                                  sep = ";") )
  assign(paste0("indTrainFLEX2.4_", k) , read.csv(paste0(path0, "TRAIN/WHOLE/FLEX2.4/individual_survival/",k,"_ind_flex24.csv"),
                                                  sep = ";") )
  #WG
  assign(paste0("indTrainWG_", k) , read.csv(paste0(path0, "TRAIN/WHOLE/WG/individual_survival/",k,"_ind_WG.csv"),
                                             sep = ";") )
}

##### calcul de ROC net

AUC_calc_W <- function(k) {
  datatrain <- get(paste0("DATATrain_", k))
  ind_estimP <- get(paste0("indTrainPLANN_", k))
  ind_estimF1.2 <- get(paste0("indTrainFLEX1.2_", k))
  ind_estimF1.4 <- get(paste0("indTrainFLEX1.4_", k))
  ind_estimF2.2 <- get(paste0("indTrainFLEX2.2_", k))
  ind_estimF2.4 <- get(paste0("indTrainFLEX2.4_", k))
  ind_estimWG <- get(paste0("indTrainWG_", k))
  
  rm(list = c(paste0("DATATrain_",k), paste0("indTrainPLANN_",k),
              paste0("indTrainFLEX1.2_",k), paste0("indTrainFLEX1.4_",k),
              paste0("indTrainFLEX2.2_",k), paste0("indTrainFLEX2.4_",k),
              paste0("indTrainWG_",k)))
  
  # Initialiser les vecteurs de résultats
  hold_P <- hold_F1.2 <- hold_F1.4 <- hold_F2.2 <- hold_F2.4 <- hold_WG <- c()
  # sort(unique(1-ind_estimP[, l + 1]))
  # unique( quantile( 1-ind_estimP[, l + 1],  probs= seq(0,1, by = 0.01)) ) #pas en dessous de 20 points 0.025
  
  for (l in seq_along(newtimes)) {
    hold_P    <- c(hold_P,    roc.net(datatrain$times, datatrain$status, 1-ind_estimP[,l + 1], datatrain$age, datatrain$sexchara,
                                      datatrain$year, slopop, pro.time = newtimes[l], 
                                      cut.off = unique( quantile(1-ind_estimP[,l+1],probs= seq(0,1, by = 0.01)) ) )$auc)
    hold_F1.2 <- c(hold_F1.2, roc.net(datatrain$times, datatrain$status, 1-ind_estimF1.2[, l], datatrain$age, datatrain$sexchara,
                                      datatrain$year, slopop, pro.time = newtimes[l],
                                      cut.off = unique( quantile( 1-ind_estimF1.2[, l],  probs= seq(0,1, by = 0.01)) ))$auc)
    hold_F1.4 <- c(hold_F1.4, roc.net(datatrain$times, datatrain$status, 1-ind_estimF1.4[, l], datatrain$age, datatrain$sexchara,
                                      datatrain$year, slopop, pro.time = newtimes[l],
                                      cut.off = unique( quantile( 1-ind_estimF1.4[, l],  probs= seq(0,1, by = 0.01)) ))$auc)
    hold_F2.2 <- c(hold_F2.2, roc.net(datatrain$times, datatrain$status, 1-ind_estimF2.2[, l], datatrain$age, datatrain$sexchara,
                                      datatrain$year, slopop, pro.time = newtimes[l], 
                                      cut.off = unique( quantile( 1-ind_estimF2.2[, l],  probs= seq(0,1, by = 0.01)) ))$auc)
    hold_F2.4 <- c(hold_F2.4, roc.net(datatrain$times, datatrain$status, 1-ind_estimF2.4[, l], datatrain$age, datatrain$sexchara,
                                      datatrain$year, slopop, pro.time = newtimes[l],
                                      cut.off = unique( quantile( 1-ind_estimF2.4[, l],  probs= seq(0,1, by = 0.01)) ))$auc)
    hold_WG   <- c(hold_WG,   roc.net(datatrain$times, datatrain$status, 1-ind_estimWG[, l],   datatrain$age, datatrain$sexchara,
                                      datatrain$year, slopop, pro.time = newtimes[l],
                                      cut.off = unique( quantile( 1-ind_estimWG[, l],  probs= seq(0,1, by = 0.01)) ))$auc)
  }
  
  return(list(
    P     = hold_P,
    F1.2  = hold_F1.2,
    F1.4  = hold_F1.4,
    F2.2  = hold_F2.2,
    F2.4  = hold_F2.4,
    WG    = hold_WG
  ))
}

results_list <- mclapply(iterations, AUC_calc_W, mc.cores = detectCores() - 4)

#supression des bases inutiles 
for(k in iterations){
  rm(list = c(paste0("DATATrain_",k), paste0("indTrainPLANN_",k),
              paste0("indTrainFLEX1.2_",k), paste0("indTrainFLEX1.4_",k),
              paste0("indTrainFLEX2.2_",k), paste0("indTrainFLEX2.4_",k),
              paste0("indTrainWG_",k)))
}
# Transformer en matrices pour chaque méthode :
AUC_WHOLE_P    <- do.call(rbind, lapply(results_list, function(x) x$P))
AUC_WHOLE_F1.2 <- do.call(rbind, lapply(results_list, function(x) x$F1.2))
AUC_WHOLE_F1.4 <- do.call(rbind, lapply(results_list, function(x) x$F1.4))
AUC_WHOLE_F2.2 <- do.call(rbind, lapply(results_list, function(x) x$F2.2))
AUC_WHOLE_F2.4 <- do.call(rbind, lapply(results_list, function(x) x$F2.4))
AUC_WHOLE_WG   <- do.call(rbind, lapply(results_list, function(x) x$WG))


### pour renommer les colonnes avec les temps de prognostic
for(l in c("P","F1.2","F1.4","F2.2","F2.4","WG")){
  dathold <- get(paste0("AUC_WHOLE_", l))
  colnames(dathold) <- c("1 year", "3 years", "5 years", "10 years")
  rownames(dathold) <- c(iterations)
  assign(paste0("AUC_WHOLE_", l),dathold)
}
## moyennes par temps 

for(l in c("P","F1.2","F1.4","F2.2","F2.4","WG")){
  
  assign(paste0("AUC_means_",l) ,  colMeans( get(paste0("AUC_WHOLE_",l)) ) )
  
}
################## STRATES

##importation des bases de données / survies individuelles prédites par les modèles 

for(k in iterations){
  ##bases de données 
  for (j in strata_names){
    assign(paste0("DATATrain_",j,"_", k) , read.csv(paste0(path0, "DATAFRAMES/TRAIN/",j,"/",k,"_dftrain_",j,".csv"),
                                                    sep = ";") )
    ##survies individuelles 
    #plann
    assign(paste0("indTrainPLANN_",j,"_", k) , read.csv(paste0(path0, "TRAIN/STRATA/PLANN/individual_survival/",k,"_ind_PLANN_",j,".csv"),
                                                        sep = ";") )
    #flex1
    assign(paste0("indTrainFLEX1.2_",j,"_", k) , read.csv(paste0(path0, "TRAIN/STRATA/FLEX1.2/individual_survival/",k,"_ind_flex12_",j,".csv"),
                                                          sep = ";") )
    assign(paste0("indTrainFLEX1.4_",j,"_", k) , read.csv(paste0(path0, "TRAIN/STRATA/FLEX1.4/individual_survival/",k,"_ind_flex14_",j,".csv"),
                                                          sep = ";") )
    #flex2
    assign(paste0("indTrainFLEX2.2_",j,"_", k) , read.csv(paste0(path0, "TRAIN/STRATA/FLEX2.2/individual_survival/",k,"_ind_flex22_",j,".csv"),
                                                          sep = ";") )
    assign(paste0("indTrainFLEX2.4_",j,"_", k) , read.csv(paste0(path0, "TRAIN/STRATA/FLEX2.4/individual_survival/",k,"_ind_flex24_",j,".csv"),
                                                          sep = ";") )
    #WG
    assign(paste0("indTrainWG_",j,"_", k) , read.csv(paste0(path0, "TRAIN/STRATA/WG/individual_survival/",k,"_ind_WG_",j,".csv"),
                                                     sep = ";") )
  }
}

##### calcul de ROC net

AUC_calc_S <- function(k) {
  auc_results <- list()
  
  for (j in strata_names) {
    datatrain <- get(paste0("DATATrain_", j, "_", k))
    ind_estimP <- get(paste0("indTrainPLANN_", j, "_", k))
    ind_estimF1.2 <- get(paste0("indTrainFLEX1.2_", j, "_", k))
    ind_estimF1.4 <- get(paste0("indTrainFLEX1.4_", j, "_", k))
    ind_estimF2.2 <- get(paste0("indTrainFLEX2.2_", j, "_", k))
    ind_estimF2.4 <- get(paste0("indTrainFLEX2.4_", j, "_", k))
    ind_estimWG <- get(paste0("indTrainWG_", j, "_", k))
    
    hold_P <- hold_F1.2 <- hold_F1.4 <- hold_F2.2 <- hold_F2.4 <- hold_WG <- c()
    
    for (l in seq_along(newtimes)) {
      hold_P    <- c(hold_P,    roc.net(datatrain$times, datatrain$status, 1-ind_estimP[,l + 1], datatrain$age, datatrain$sexchara,
                                        datatrain$year, slopop, pro.time = newtimes[l], 
                                        cut.off = unique( quantile(1-ind_estimP[,l+1],probs= seq(0,1, by = 0.01)) ) )$auc)
      hold_F1.2 <- c(hold_F1.2, roc.net(datatrain$times, datatrain$status, 1-ind_estimF1.2[, l], datatrain$age, datatrain$sexchara,
                                        datatrain$year, slopop, pro.time = newtimes[l],
                                        cut.off = unique( quantile( 1-ind_estimF1.2[, l],  probs= seq(0,1, by = 0.01)) ))$auc)
      hold_F1.4 <- c(hold_F1.4, roc.net(datatrain$times, datatrain$status, 1-ind_estimF1.4[, l], datatrain$age, datatrain$sexchara,
                                        datatrain$year, slopop, pro.time = newtimes[l],
                                        cut.off = unique( quantile( 1-ind_estimF1.4[, l],  probs= seq(0,1, by = 0.01)) ))$auc)
      hold_F2.2 <- c(hold_F2.2, roc.net(datatrain$times, datatrain$status, 1-ind_estimF2.2[, l], datatrain$age, datatrain$sexchara,
                                        datatrain$year, slopop, pro.time = newtimes[l], 
                                        cut.off = unique( quantile( 1-ind_estimF2.2[, l],  probs= seq(0,1, by = 0.01)) ))$auc)
      hold_F2.4 <- c(hold_F2.4, roc.net(datatrain$times, datatrain$status, 1-ind_estimF2.4[, l], datatrain$age, datatrain$sexchara,
                                        datatrain$year, slopop, pro.time = newtimes[l],
                                        cut.off = unique( quantile( 1-ind_estimF2.4[, l],  probs= seq(0,1, by = 0.01)) ))$auc)
      hold_WG   <- c(hold_WG,   roc.net(datatrain$times, datatrain$status, 1-ind_estimWG[, l],   datatrain$age, datatrain$sexchara,
                                        datatrain$year, slopop, pro.time = newtimes[l],
                                        cut.off = unique( quantile( 1-ind_estimWG[, l],  probs= seq(0,1, by = 0.01)) ))$auc)
    }
    
    auc_results[[j]] <- list(
      k = k,
      P = hold_P,
      F1.2 = hold_F1.2,
      F1.4 = hold_F1.4,
      F2.2 = hold_F2.2,
      F2.4 = hold_F2.4,
      WG = hold_WG
    )
  }
  
  return(auc_results)
}

res_all <- mclapply(iterations, AUC_calc_S, mc.cores = detectCores() - 4)

##supression des bases inutiles pour faire de la place
for(k in iterations){
  for(j in strata_names){
    rm(list = c(
      paste0("DATATrain_", j, "_", k),
      paste0("indTrainPLANN_", j, "_", k),
      paste0("indTrainFLEX1.2_", j, "_", k),
      paste0("indTrainFLEX1.4_", j, "_", k),
      paste0("indTrainFLEX2.2_", j, "_", k),
      paste0("indTrainFLEX2.4_", j, "_", k),
      paste0("indTrainWG_", j, "_", k)
    ))
  }
}

# Transformer en matrices pour chaque méthode :

for(j in strata_names){
  assign(paste0("AUC_",j,"_P"), do.call(rbind, lapply(res_all, function(x) x[[j]]$P)) )
  assign(paste0("AUC_",j,"_F1.2"), do.call(rbind, lapply(res_all, function(x) x[[j]]$F1.2)) )
  assign(paste0("AUC_",j,"_F1.4"), do.call(rbind, lapply(res_all, function(x) x[[j]]$F1.4)) )
  assign(paste0("AUC_",j,"_F2.2"), do.call(rbind, lapply(res_all, function(x) x[[j]]$F2.2)) )
  assign(paste0("AUC_",j,"_F2.4"), do.call(rbind, lapply(res_all, function(x) x[[j]]$F2.4)) )
  assign(paste0("AUC_",j,"_WG"), do.call(rbind, lapply(res_all, function(x) x[[j]]$WG)) )
}
### pour renommer les colonnes avec les temps de prognostic
for(l in c("P","F1.2","F1.4","F2.2","F2.4","WG")){
  for(j in strata_names){
    dathold <- get(paste0("AUC_", j,"_",l))
    colnames(dathold) <- c("1 year", "3 years", "5 years", "10 years")
    rownames(dathold) <- c(iterations)
    assign(paste0("AUC_",j,"_", l),dathold)
  }
}
## moyennes par temps 

for(l in c("P","F1.2","F1.4","F2.2","F2.4","WG")){
  for(j in strata_names){
    assign(paste0("AUC_means_",j,"_",l) ,  colMeans( get(paste0("AUC_",j,"_",l)) ) )
  }
}


### sous forme de liste pour libérer de l'espace dans l'evt (pour les moyennes)

ROCmean_results <- list()

for(l in c("P","F1.2","F1.4","F2.2","F2.4","WG")){
  ROCmean_results[['WHOLE']][[l]] <- list()
  ROCmean_results[['WHOLE']][[l]] <- get( paste0("AUC_means_",l) )
  rm(list = paste0("AUC_means_",l))
  
  for(j in strata_names){
    ROCmean_results[[j]][[l]] <- get(paste0("AUC_means_", j, "_", l)) 
    rm(list = paste0("AUC_means_", j, "_", l))
  }
  
}
### meme chose pour les valeurs pour chaque itération 


ROC_results <- list()

for(l in c("P","F1.2","F1.4","F2.2","F2.4","WG")){
  ROC_results[['WHOLE']][[l]] <- list()
  ROC_results[['WHOLE']][[l]] <- get( paste0("AUC_WHOLE_",l) )
  rm(list = paste0("AUC_WHOLE_",l) )
  for(j in strata_names){
    ROC_results[[j]][[l]] <- get(paste0("AUC_", j, "_", l)) 
    rm(list = paste0("AUC_", j, "_", l))
  }
  
}
#################### VALIDATION ######################

### WHOLE

##importation des bases de données / survies individuelles prédites par les modèles 
for(k in iterations){
  ##bases de données 
  assign(paste0("DATAValid_", k) , read.csv(paste0(path0, "DATAFRAMES/VALID/WHOLE/",k,"_dfvalid.csv"),
                                            sep = ";") )
  
  ##survies individuelles 
  #plann
  assign(paste0("indValidPLANN_", k) , read.csv(paste0(path0, "VALID/WHOLE/PLANN/individual_survival/",k,"_ind_PLANNval.csv"),
                                                sep = ";") )
  
  #flex1
  assign(paste0("indValidFLEX1.2_",k) , read.csv(paste0(path0, "VALID/WHOLE/FLEX1.2/individual_survival/",k,"_ind_flexval12.csv"),
                                                 sep = ";") )
  assign(paste0("indValidFLEX1.4_",k) , read.csv(paste0(path0, "VALID/WHOLE/FLEX1.4/individual_survival/",k,"_ind_flexval14.csv"),
                                                 sep = ";") )
  #flex2
  assign(paste0("indValidFLEX2.2_",k) , read.csv(paste0(path0, "VALID/WHOLE/FLEX2.2/individual_survival/",k,"_ind_flexval22.csv"),
                                                 sep = ";") )
  assign(paste0("indValidFLEX2.4_",k) , read.csv(paste0(path0, "VALID/WHOLE/FLEX2.4/individual_survival/",k,"_ind_flexval24.csv"),
                                                 sep = ";") )
  #WG
  assign(paste0("indValidWG_", k) , read.csv(paste0(path0, "VALID/WHOLE/WG/individual_survival/",k,"_ind_flexWGval.csv"),
                                             sep = ";") )
}

##### calcul de ROC net

##creation d'objets pour enregistrer les AUC de chaque itération aux différents temps et pour les différents modèles 

AUC_calc_val_W <- function(k){ 
  
  dataval <- get(paste0("DATAValid_",k))
  ind_estimP <- get(paste0("indValidPLANN_",k))
  ind_estimF1.2 <- get(paste0("indValidFLEX1.2_",k))
  ind_estimF1.4 <- get(paste0("indValidFLEX1.4_", k))
  ind_estimF2.2 <- get(paste0("indValidFLEX2.2_", k))
  ind_estimF2.4 <- get(paste0("indValidFLEX2.4_", k))
  ind_estimWG <- get(paste0("indValidWG_",k))
  
  rm(list = c(paste0("DATAValid_",k), paste0("indValidPLANN_",k),paste0("indValidFLEX1.2_",k), paste0("indValidFLEX1.4_", k),
              paste0("indValidFLEX2.2_", k),paste0("indValidFLEX2.4_", k),paste0("indValidWG_",k)))
  
  # Initialiser les vecteurs de résultats
  hold_P <- hold_F1.2 <- hold_F1.4 <- hold_F2.2 <- hold_F2.4 <- hold_WG <- c()
  
  ###PLANN
  for (l in seq_along(newtimes)) {
    hold_P    <- c(hold_P,    roc.net(dataval$times, dataval$status, 1-ind_estimP[,l + 1], dataval$age, dataval$sexchara,
                                      dataval$year, slopop, pro.time = newtimes[l], 
                                      cut.off = unique( quantile(1-ind_estimP[,l+1],probs= seq(0,1, by = 0.01)) ) )$auc)
    hold_F1.2 <- c(hold_F1.2, roc.net(dataval$times, dataval$status, 1-ind_estimF1.2[, l], dataval$age, dataval$sexchara,
                                      dataval$year, slopop, pro.time = newtimes[l],
                                      cut.off = unique( quantile( 1-ind_estimF1.2[, l],  probs= seq(0,1, by = 0.01)) ))$auc)
    hold_F1.4 <- c(hold_F1.4, roc.net(dataval$times, dataval$status, 1-ind_estimF1.4[, l], dataval$age, dataval$sexchara,
                                      dataval$year, slopop, pro.time = newtimes[l],
                                      cut.off = unique( quantile( 1-ind_estimF1.4[, l],  probs= seq(0,1, by = 0.01)) ))$auc)
    hold_F2.2 <- c(hold_F2.2, roc.net(dataval$times, dataval$status, 1-ind_estimF2.2[, l], dataval$age, dataval$sexchara,
                                      dataval$year, slopop, pro.time = newtimes[l], 
                                      cut.off = unique( quantile( 1-ind_estimF2.2[, l],  probs= seq(0,1, by = 0.01)) ))$auc)
    hold_F2.4 <- c(hold_F2.4, roc.net(dataval$times, dataval$status, 1-ind_estimF2.4[, l], dataval$age, dataval$sexchara,
                                      dataval$year, slopop, pro.time = newtimes[l],
                                      cut.off = unique( quantile( 1-ind_estimF2.4[, l],  probs= seq(0,1, by = 0.01)) ))$auc)
    hold_WG   <- c(hold_WG,   roc.net(dataval$times, dataval$status, 1-ind_estimWG[, l],   dataval$age, dataval$sexchara,
                                      dataval$year, slopop, pro.time = newtimes[l],
                                      cut.off = unique( quantile( 1-ind_estimWG[, l],  probs= seq(0,1, by = 0.01)) ))$auc)
  }
  
  return(list(
    P     = hold_P,
    F1.2  = hold_F1.2,
    F1.4  = hold_F1.4,
    F2.2  = hold_F2.2,
    F2.4  = hold_F2.4,
    WG    = hold_WG
  ))
}

resultval_list <- mclapply(iterations, AUC_calc_val_W, mc.cores = detectCores() - 4)

#supression des bases inutiles 
for(k in iterations){
  rm(list = c(paste0("DATAValid_",k), paste0("indValidPLANN_",k),
              paste0("indValidFLEX1.2_",k), paste0("indValidFLEX1.4_",k),
              paste0("indValidFLEX2.2_",k), paste0("indValidFLEX2.4_",k),
              paste0("indValidWG_",k)))
}
# Transformer en matrices pour chaque méthode :
AUCval_WHOLE_P    <- do.call(rbind, lapply(resultval_list, function(x) x$P))
AUCval_WHOLE_F1.2 <- do.call(rbind, lapply(resultval_list, function(x) x$F1.2))
AUCval_WHOLE_F1.4 <- do.call(rbind, lapply(resultval_list, function(x) x$F1.4))
AUCval_WHOLE_F2.2 <- do.call(rbind, lapply(resultval_list, function(x) x$F2.2))
AUCval_WHOLE_F2.4 <- do.call(rbind, lapply(resultval_list, function(x) x$F2.4))
AUCval_WHOLE_WG   <- do.call(rbind, lapply(resultval_list, function(x) x$WG))

### pour renommer les colonnes avec les temps de prognostic
for(l in c("P","F1.2","F1.4","F2.2","F2.4","WG")){
  dathold <- get(paste0("AUCval_WHOLE_", l))
  colnames(dathold) <- c("1 year", "3 years", "5 years", "10 years")
  rownames(dathold) <- c(iterations)
  assign(paste0("AUCval_WHOLE_", l),dathold)
}
## moyennes par temps 

for(l in c("P","F1.2","F1.4","F2.2","F2.4","WG")){
  
  assign(paste0("AUCval_means_",l) ,  colMeans( get(paste0("AUCval_WHOLE_",l)) ) )
  
}
################## STRATES

##importation des bases de données / survies individuelles prédites par les modèles 

for(k in iterations){
  ##bases de données 
  for (j in strata_names){
    assign(paste0("DATAValid_",j,"_", k) , read.csv(paste0(path0, "DATAFRAMES/VALID/",j,"/",k,"_dfvalid_",j,".csv"),
                                                    sep = ";") )
    ##survies individuelles 
    #plann
    assign(paste0("indValidPLANN_",j,"_", k) , read.csv(paste0(path0, "VALID/STRATA/PLANN/individual_survival/",k,"_ind_PLANNval_",j,".csv"),
                                                        sep = ";") )
    #flex1
    assign(paste0("indValidFLEX1.2_",j,"_", k) , read.csv(paste0(path0, "VALID/STRATA/FLEX1.2/individual_survival/",k,"_ind_flexval12_",j,".csv"),
                                                          sep = ";") )
    assign(paste0("indValidFLEX1.4_",j,"_", k) , read.csv(paste0(path0, "VALID/STRATA/FLEX1.4/individual_survival/",k,"_ind_flexval14_",j,".csv"),
                                                          sep = ";") )
    #flex2
    assign(paste0("indValidFLEX2.2_",j,"_", k) , read.csv(paste0(path0, "VALID/STRATA/FLEX2.2/individual_survival/",k,"_ind_flexval22_",j,".csv"),
                                                          sep = ";") )
    assign(paste0("indValidFLEX2.4_",j,"_", k) , read.csv(paste0(path0, "VALID/STRATA/FLEX2.4/individual_survival/",k,"_ind_flexval24_",j,".csv"),
                                                          sep = ";") )
    #WG
    assign(paste0("indValidWG_",j,"_", k) , read.csv(paste0(path0, "VALID/STRATA/WG/individual_survival/",k,"_ind_WGval_",j,".csv"),
                                                     sep = ";") )
  }
}

##### calcul de ROC net

AUC_calc_val_S <- function(k){ 
  auc_results <- list()
  for(j in strata_names){
    dataval <- get(paste0("DATAValid_",j,"_", k))
    ind_estimP <- get(paste0("indValidPLANN_",j,"_", k))
    ind_estimF1.2 <- get(paste0("indValidFLEX1.2_",j,"_", k))
    ind_estimF1.4 <- get(paste0("indValidFLEX1.4_",j,"_", k))
    ind_estimF2.2 <- get(paste0("indValidFLEX2.2_",j,"_", k))
    ind_estimF2.4 <- get(paste0("indValidFLEX2.4_",j,"_", k))
    ind_estimWG <- get(paste0("indValidWG_",j,"_", k))
    
    rm(list = c(paste0("DATAValid_",j,"_", k), paste0("indValidPLANN_",j,"_", k), paste0("indValidFLEX1.2_",j,"_", k),paste0("indValidFLEX1.4_",j,"_", k),
                paste0("indValidFLEX2.2_",j,"_", k), paste0("indValidFLEX2.2_",j,"_", k), paste0("indValidWG_",j,"_", k)))
    
    hold_P <- hold_F1.2 <- hold_F1.4 <- hold_F2.2 <- hold_F2.4 <- hold_WG <- c()
    
    for (l in seq_along(newtimes)) {
      hold_P    <- c(hold_P,    roc.net(dataval$times, dataval$status, 1-ind_estimP[,l + 1], dataval$age, dataval$sexchara,
                                        dataval$year, slopop, pro.time = newtimes[l], 
                                        cut.off = unique( quantile(1-ind_estimP[,l+1],probs= seq(0,1, by = 0.01)) ) )$auc)
      hold_F1.2 <- c(hold_F1.2, roc.net(dataval$times, dataval$status, 1-ind_estimF1.2[, l], dataval$age, dataval$sexchara,
                                        dataval$year, slopop, pro.time = newtimes[l],
                                        cut.off = unique( quantile( 1-ind_estimF1.2[, l],  probs= seq(0,1, by = 0.01)) ))$auc)
      hold_F1.4 <- c(hold_F1.4, roc.net(dataval$times, dataval$status, 1-ind_estimF1.4[, l], dataval$age, dataval$sexchara,
                                        dataval$year, slopop, pro.time = newtimes[l],
                                        cut.off = unique( quantile( 1-ind_estimF1.4[, l],  probs= seq(0,1, by = 0.01)) ))$auc)
      hold_F2.2 <- c(hold_F2.2, roc.net(dataval$times, dataval$status, 1-ind_estimF2.2[, l], dataval$age, dataval$sexchara,
                                        dataval$year, slopop, pro.time = newtimes[l], 
                                        cut.off = unique( quantile( 1-ind_estimF2.2[, l],  probs= seq(0,1, by = 0.01)) ))$auc)
      hold_F2.4 <- c(hold_F2.4, roc.net(dataval$times, dataval$status, 1-ind_estimF2.4[, l], dataval$age, dataval$sexchara,
                                        dataval$year, slopop, pro.time = newtimes[l],
                                        cut.off = unique( quantile( 1-ind_estimF2.4[, l],  probs= seq(0,1, by = 0.01)) ))$auc)
      hold_WG   <- c(hold_WG,   roc.net(dataval$times, dataval$status, 1-ind_estimWG[, l],   dataval$age, dataval$sexchara,
                                        dataval$year, slopop, pro.time = newtimes[l],
                                        cut.off = unique( quantile( 1-ind_estimWG[, l],  probs= seq(0,1, by = 0.01)) ))$auc)
    }
    
    auc_results[[j]] <- list(
      k = k,
      P = hold_P,
      F1.2 = hold_F1.2,
      F1.4 = hold_F1.4,
      F2.2 = hold_F2.2,
      F2.4 = hold_F2.4,
      WG = hold_WG
    )
  }
  
  return(auc_results)
}

resval_all <- mclapply(iterations, AUC_calc_val_S, mc.cores = detectCores() - 4)

for(k in iterations){
  for(j in strata_names){
    rm(list = c(
      paste0("DATAValid_", j, "_", k),
      paste0("indValidPLANN_", j, "_", k),
      paste0("indValidFLEX1.2_", j, "_", k),
      paste0("indValidFLEX1.4_", j, "_", k),
      paste0("indValidFLEX2.2_", j, "_", k),
      paste0("indValidFLEX2.4_", j, "_", k),
      paste0("indValidWG_", j, "_", k)
    ))
  }
}

# Transformer en matrices pour chaque méthode :

for(j in strata_names){
  assign(paste0("AUCval_",j,"_P"), do.call(rbind, lapply(resval_all, function(x) x[[j]]$P)) )
  assign(paste0("AUCval_",j,"_F1.2"), do.call(rbind, lapply(resval_all, function(x) x[[j]]$F1.2)) )
  assign(paste0("AUCval_",j,"_F1.4"), do.call(rbind, lapply(resval_all, function(x) x[[j]]$F1.4)) )
  assign(paste0("AUCval_",j,"_F2.2"), do.call(rbind, lapply(resval_all, function(x) x[[j]]$F2.2)) )
  assign(paste0("AUCval_",j,"_F2.4"), do.call(rbind, lapply(resval_all, function(x) x[[j]]$F2.4)) )
  assign(paste0("AUCval_",j,"_WG"), do.call(rbind, lapply(resval_all, function(x) x[[j]]$WG)) )
}
### pour renommer les colonnes avec les temps de prognostic
for(l in c("P","F1.2","F1.4","F2.2","F2.4","WG")){
  for(j in strata_names){
    dathold <- get(paste0("AUCval_", j,"_",l))
    colnames(dathold) <- c("1 year", "3 years", "5 years", "10 years")
    rownames(dathold) <- c(iterations)
    assign(paste0("AUCval_",j,"_", l),dathold)
  }
}
## moyennes par temps 

for(l in c("P","F1.2","F1.4","F2.2","F2.4","WG")){
  for(j in strata_names){
    assign(paste0("AUCval_means_",j,"_",l) ,  colMeans( get(paste0("AUCval_",j,"_",l)) ) )
  }
}


### sous forme de liste pour libérer de l'espace dans l'evt (pour les moyennes)

ROCmeanval_results <- list()

for(l in c("P","F1.2","F1.4","F2.2","F2.4","WG")){
  ROCmeanval_results[['WHOLE']][[l]] <- list()
  ROCmeanval_results[['WHOLE']][[l]] <- get( paste0("AUCval_means_",l) )
  rm(list = paste0("AUCval_means_",l))
  
  for(j in strata_names){
    ROCmeanval_results[[j]][[l]] <- get(paste0("AUCval_means_", j, "_", l)) 
    rm(list = paste0("AUCval_means_", j, "_", l))
  }
  
}
### meme chose pour les valeurs pour chaque itération 


ROCval_results <- list()

for(l in c("P","F1.2","F1.4","F2.2","F2.4","WG")){
  ROCval_results[['WHOLE']][[l]] <- list()
  ROCval_results[['WHOLE']][[l]] <- get( paste0("AUCval_WHOLE_",l) )
  rm(list = paste0("AUCval_WHOLE_",l) )
  for(j in strata_names){
    ROCval_results[[j]][[l]] <- get(paste0("AUCval_", j, "_", l)) 
    rm(list = paste0("AUCval_", j, "_", l))
  }
  
}

print(Sys.time()-start)


save.image("~/Documents/Simulations/Simulations mai 2025/Résultats/N1000_results/NPH/1000ite_1000ind_NPH_HC_04_07.Rdata")

rm(list = ls(envir = .GlobalEnv), envir = .GlobalEnv)

### 1000IndEND

##################################################################################################
##################################################################################################
##################################################################################################
                                          ##3000##
##################################################################################################
##################################################################################################
##################################################################################################

### 3000SiSTA

library(survivalNET)
library(survivalPLANN)
library(parallel)
library(doParallel)

#################################################################################################
### path ###

path0 <- "~/Documents/Simulations/Simulations mai 2025/N3000/Output simulations_N3000_NPH_HC/"
path1 <- "~/Documents/Simulations/BASES/"
# path0 <- paste0(getwd(),"/")
############

#importation des données
data("fr.ratetable")
data("slopop")
data(colrec)

# data management


colrec$agey <- colrec$age/365.241

colrec$sexchara <- "male"
colrec$sexchara[colrec$sex==2] <- "female"

colrec$colon <- 1*(as.character(colrec$site)=="colon")

colrec <- colrec[colrec$stage!=99,] # suppresion des patients avec valeur manquante
table(colrec$stage)

colrec$stage1 <- 1*(colrec$stage==1)
colrec$stage2 <- 1*(colrec$stage==2)
colrec$stage3 <- 1*(colrec$stage==3)

table(colrec$stat)

colrec$agey10 <- colrec$agey/10
names(colrec)


###################################################################################################
### Fonctions de survie nette, attendue, observée et de simulation de temps
###################################################################################################


Sn <- function(time, sigma, nu, theta, beta, covariates)
{
  exp((1-(1+(time/(exp(sigma)))^exp(nu))^(1/exp(theta))) * exp(sum(covariates*beta)) )
}


###################################################################################################
##############
# simulation#
#############

simulate_iteration <- function(i){
  
  set.seed(i)
  
  start <- Sys.time()
  
  data <- read.csv(file = paste0(path1,"ind3000NPH_HC/",i,"_df3000.csv"), sep = ";")
  Tsigma <- 12.6 
  Tnu <- -0.5
  Ttheta <- 0
  betaZ <- c(0.9,2.7,0.3,-0.1,-0.1)
  
  TsigmaHC <- 13.0
  TnuHC <- -0.55
  TthetaHC <- 0
  betaZHC <- c(1,2.9,0.2)
  
  #### Hommes Rectum
  TsigmaHR <- 12.75
  TnuHR <- -0.3
  TthetaHR <- 0                      
  betaZHR <- c(1.2,3,0.35)
  
  #### Femmes Colon
  TsigmaFC <- 14
  TnuFC <- -0.55
  TthetaFC <- 0
  betaZFC <- c(0.6,2.6,-0.1)
  
  #### Femmes Rectum
  TsigmaFR <- 11.5
  TnuFR <- -0.4
  TthetaFR <- 0
  betaZFR <- c(0.9,2.7,0.16)
  ############### validation sample
  N <- dim(data)[1]
  N_train <- N/2 #500
  data_train = data[1:N_train,]
  data_valid = data[(N_train+1):dim(data)[1],]
  
  ###############################
  ########## STRATES ############
  ###############################
  
  ### training
  #################### Hommes ####################
  
  ### hommes & colon
  
  data_train_HC <- data_train[data_train$sex == 1 & data_train$colon == 1,]
  
  ### hommes & rectum
  
  data_train_HR <- data_train[data_train$sex == 1 & data_train$colon == 0,]
  
  #################### Femmes ####################
  
  ### femmes & colon
  
  data_train_FC <- data_train[data_train$sex == 2 & data_train$colon == 1,]
  
  ### femmes & rectum
  
  data_train_FR <- data_train[data_train$sex == 2 & data_train$colon == 0,]
  
  
  ########################################
  ### validation
  
  #################### Hommes ####################
  
  ### hommes & colon
  
  data_valid_HC <- data_valid[data_valid$sex == 1 & data_valid$colon == 1,]
  
  ### hommes & rectum
  
  data_valid_HR <- data_valid[data_valid$sex == 1 & data_valid$colon == 0,]
  
  #################### Femmes ####################
  
  ### femmes & colon
  
  data_valid_FC <- data_valid[data_valid$sex == 2 & data_valid$colon == 1,]
  
  ### femmes & rectum
  
  data_valid_FR <- data_valid[data_valid$sex == 2 & data_valid$colon == 0,]
  
  ########################################
  
  ######################valeurs issues du modèle theorique
  
  pro.time <- max(data_train$times[data_train$status==1])
  
  ###temps d'évaluations (t =1, 3, 5, 10 ans)
  
  newtimes <- c(365.241, 365.241*3, 365.241*5, 365.241*10) 
  # @THOMAS : 80% 50% 20% sujets à risque 
  
  Zdata_train <- cbind(data_train$stage2, data_train$stage3, data_train$agey10, data_train$sex, data_train$colon)
  colnames(Zdata_train) <- c("stage2", "stage3", "agey10", "sex", "colon")
  
  funtheo <- function(i){
    if(Zdata_train[i,"sex"] == 1 & Zdata_train[i,"colon"] == 1){
      Tsigma = TsigmaHC
      Tnu = TnuHC
      Ttheta = TthetaHC
      betaZ = betaZHC
    }
    # Hommes rectum
    if(Zdata_train[i,"sex"] == 1 & Zdata_train[i,"colon"] == 0){
      Tsigma = TsigmaHR
      Tnu = TnuHR
      Ttheta = TthetaHR
      betaZ = betaZHR
    }
    # Femmes colon
    if(Zdata_train[i,"sex"] == 2 & Zdata_train[i,"colon"] == 1){
      Tsigma = TsigmaFC
      Tnu = TnuFC
      Ttheta = TthetaFC
      betaZ = betaZFC
    }
    # Femmes rectum
    if(Zdata_train[i,"sex"] == 2 & Zdata_train[i,"colon"] == 0){
      Tsigma = TsigmaFR
      Tnu = TnuFR
      Ttheta = TthetaFR
      betaZ = betaZFR
    }
    Sn(time = newtimes, sigma = Tsigma, nu = Tnu, theta = Ttheta, beta = betaZ, covariates = Zdata_train[i,-c(4,5)]) }
  
  theopred <- t(sapply(1:dim(Zdata_train)[1], FUN = funtheo ))
  
  colnames(theopred) <- newtimes
  
  mean.theo <- apply(theopred, FUN="mean", MARGIN=2)
  
  ###strates
  #####
  
  strata_names <- c("HC", "HR", "FC", "FR")
  
  for(j in strata_names){
    
    data_train_var <- get(paste0("data_train_", j))
    
    Zdata_train <- cbind(data_train_var$stage2, data_train_var$stage3, data_train_var$agey10, data_train_var$sex, data_train_var$colon)
    colnames(Zdata_train) <- c("stage2", "stage3", "agey10", "sex", "colon")
    
    funtheostrat <- function(i){
      if(Zdata_train[i,"sex"] == 1 & Zdata_train[i,"colon"] == 1){
        Tsigma = TsigmaHC
        Tnu = TnuHC
        Ttheta = TthetaHC
        betaZ = betaZHC
      }
      # Hommes rectum
      if(Zdata_train[i,"sex"] == 1 & Zdata_train[i,"colon"] == 0){
        Tsigma = TsigmaHR
        Tnu = TnuHR
        Ttheta = TthetaHR
        betaZ = betaZHR
      }
      # Femmes colon
      if(Zdata_train[i,"sex"] == 2 & Zdata_train[i,"colon"] == 1){
        Tsigma = TsigmaFC
        Tnu = TnuFC
        Ttheta = TthetaFC
        betaZ = betaZFC
      }
      # Femmes rectum
      if(Zdata_train[i,"sex"] == 2 & Zdata_train[i,"colon"] == 0){
        Tsigma = TsigmaFR
        Tnu = TnuFR
        Ttheta = TthetaFR
        betaZ = betaZFR
      }
      Sn(time = newtimes, sigma = Tsigma, nu = Tnu, theta = Ttheta, beta = betaZ, covariates = Zdata_train[i,-c(4,5)]) }
    
    theopredstrat <- t(sapply(1:dim(Zdata_train)[1], FUN = funtheostrat ))
    
    colnames(theopredstrat) <- newtimes
    
    mean.theo_strat <- apply(theopredstrat, FUN = "mean", MARGIN = 2)
    
    assign(paste0("Zdata_train_", j), Zdata_train)
    assign(paste0("theopred_", j), theopredstrat)
    assign(paste0("mean.theo_", j), mean.theo_strat)  
    
  }
  
  #####
  
  ####### sur la base de validation
  
  Zdata_val <- cbind(data_valid$stage2, data_valid$stage3, data_valid$agey10, data_valid$sex, data_valid$colon)
  
  theopredval <- t(sapply(1:dim(Zdata_val)[1], FUN = function(i){
    Sn(time = newtimes, sigma = Tsigma, nu = Tnu, theta = Ttheta, beta = betaZ,  covariates = Zdata_val[i,]) } ) )
  
  colnames(theopredval) <- newtimes
  
  mean.theoval <- apply(theopredval, FUN="mean", MARGIN=2) 
  
  ###strates
  #####
  
  for(j in strata_names){
    
    data_valid_var <- get(paste0("data_valid_", j))
    
    Zdata_valid <- cbind(data_valid_var$stage2, data_valid_var$stage3, data_valid_var$agey10, data_valid_var$sex, data_valid_var$colon)
    
    theopredvalstrat <- t(sapply(1:nrow(Zdata_valid), FUN = function(i) {
      Sn(time = newtimes, sigma = Tsigma, nu = Tnu, theta = Ttheta, beta = betaZ, covariates = Zdata_valid[i,])
    }))
    
    colnames(theopredvalstrat) <- newtimes
    
    mean.theoval_strat <- apply(theopredvalstrat, FUN = "mean", MARGIN = 2)
    
    assign(paste0("Zdata_val_", j), Zdata_valid)
    assign(paste0("theopredval_", j), theopredvalstrat)
    assign(paste0("mean.theoval_", j), mean.theoval_strat)  
    
  }
  
  
  #####
  #######################################################
  ###### choix des hyper-paramètres par cross-validation
  #######################################################
  
  
  # 
  # ## Modèle PLANN
  # # 
  
  startCVplann <- Sys.time()
  
  tune.plann <- cvPLANN(formula = Surv(times, status) ~ stage2 + stage3 + agey10 + sex01 + colon, pro.time = pro.time,
                        data = data_train, cv = 10, inter= 365.241/12, size=c(2, 4, 8, 10),
                        decay=c(0.01, 0.1), maxit=1000, MaxNWts=10000)
  
  TimeCVplann <- Sys.time() - startCVplann
  
  # 
  # ####################
  # #paramètres optimaux
  # ####################
  # 
  
  # inter <- 365.241/12
  # size <- 4
  # decay <- 0.01
  # maxit <- 1000
  # MaxNWts <- 10000
  
  ##### PLANN
  
  inter <- tune.plann$optimal$inter
  
  size <- tune.plann$optimal$size
  
  decay <- tune.plann$optimal$decay
  
  maxit <- tune.plann$optimal$maxit
  
  MaxNWts <- tune.plann$optimal$MaxNWts  
  
  # ##### FLEXNET
  
  #flex1 
  
  m1.2 <- 2
  m1.4 <- 4
  failedindic1.2 <- NULL
  failedindic1.4 <- NULL
  
  #flex2  
  
  m2.2 <- 2
  m2.4 <- 4
  failedindic2.2 <- NULL
  failedindic2.4 <- NULL
  
  
  ##################
  # Modèles fittés
  ##################
  # 
  
  ### plann
  
  plann.model <- sPLANN(formula = Surv(times, status) ~ stage2 + stage3 + agey10 + sex01 + colon,
                        pro.time = pro.time, data = data_train, inter= inter, size=size, 
                        decay=decay, maxit=maxit, MaxNWts=MaxNWts)
  
  ### Weibull Généralisé 
  
  WG.model <- survivalNET(formula = Surv(times, status) ~ stage2 + stage3 + agey10 + sex01 + colon +
                            ratetable(age, year, sexchara), data = data_train,
                          ratetable=slopop, dist = "genweibull")
  
  logcoeffWG <- tail(WG.model$coefficients, 3)
  
  ################## flex1
  ########### 2 noeuds
  m1.2_original <- m1.2
  m1.2_changed <- FALSE
  
  converged <- FALSE
  
  while (m1.2 >= 0 && !converged) {
    tryCatch({
      flex.model1.2 <- survivalFLEXNET(
        formula = Surv(times, status) ~ stage2 + stage3 + agey10 + sex01 + colon +
          ratetable(age, year, sexchara),
        data = data_train,
        ratetable = slopop,
        m = m1.2
      )
      # If the model succeeds, set converged to TRUE and break out of the loop
      converged <- TRUE
    }, error = function(e) {
      # If an error occurs, decrement m1
      m1.2_changed <<- TRUE
      m1.2 <<- m1.2 - 1
    })
  }
  
  if (m1.2_changed & converged) {
    failedindic1.2 <- paste0( " : this value was decreased from ", m1.2_original)
  }
  
  if (!converged) {
    flex.model1.2 <- WG.model
    failedindic1.2 <- "flex.model1.2 didn't converge, even after decreasing the number of knots incrementaly. WG.model was used instead."
  }
  
  ########### 4 noeuds
  m1.4_original <- m1.4
  m1.4_changed <- FALSE
  
  converged <- FALSE
  
  while (m1.4 >= 0 && !converged) {
    tryCatch({
      flex.model1.4 <- survivalFLEXNET(
        formula = Surv(times, status) ~ stage2 + stage3 + agey10 + sex01 + colon +
          ratetable(age, year, sexchara),
        data = data_train,
        ratetable = slopop,
        m = m1.4
      )
      # If the model succeeds, set converged to TRUE and break out of the loop
      converged <- TRUE
    }, error = function(e) {
      # If an error occurs, decrement m1
      m1.4_changed <<- TRUE
      m1.4 <<- m1.4 - 1
    })
  }
  
  if (m1.4_changed & converged) {
    failedindic1.4 <- paste0( " : this value was decreased from ", m1.4_original)
  }
  
  if (!converged) {
    flex.model1.4 <- WG.model
    failedindic1.4 <- "flex.model1.4 didn't converge, even after decreasing the number of knots incrementaly. WG.model was used instead."
  }
  ### flex2
  ## 2 noeuds 
  m2.2_original <- m2.2
  m2.2_changed <- FALSE
  
  converged <- FALSE
  
  while (m2.2 >= 0 && !converged) {
    tryCatch({flex.model2.2 <-survivalFLEXNET(formula = Surv(times, status) ~ stage2 + stage3 + agey10 + strata(sex.organ) +
                                                ratetable(age, year, sexchara), data = data_train,
                                              ratetable=slopop, m = m2.2)
    # If the model succeeds, set converged to TRUE and break out of the loop
    converged <- TRUE
    }, error = function(e) {
      # If an error occurs, decrement m2
      m2.2_changed <<- TRUE
      m2.2 <<- m2.2 - 1
    })
  }
  
  if (m2.2_changed & converged) {
    failedindic2.2 <- paste0(" this value was decreased from ", m2.2_original)
  }
  
  if (!converged) {
    flex.model2.2 <- flex.model1.2
    failedindic2.2 <- "flex.model2.2 didn't converge, even after decreasing the number of knots incrementaly. flex.model1.2 was used instead."
    
  }
  ## 4 noeuds 
  
  m2.4_original <- m2.4
  m2.4_changed <- FALSE
  
  converged <- FALSE
  
  while (m2.4 >= 0 && !converged) {
    tryCatch({flex.model2.4 <-survivalFLEXNET(formula = Surv(times, status) ~ stage2 + stage3 + agey10 + strata(sex.organ) +
                                                ratetable(age, year, sexchara), data = data_train,
                                              ratetable=slopop, m = m2.4)
    # If the model succeeds, set converged to TRUE and break out of the loop
    converged <- TRUE
    }, error = function(e) {
      # If an error occurs, decrement m2
      m2.4_changed <<- TRUE
      m2.4 <<- m2.4 - 1
    })
  }
  
  if (m2.4_changed & converged) {
    failedindic2.4 <- paste0(" this value was decreased from ", m2.4_original)
  }
  
  if (!converged) {
    flex.model2.4 <- flex.model1.4
    failedindic2.4 <- "flex.model2.4 didn't converge, even after decreasing the number of knots incrementaly. flex.model1.4 was used instead."
    
  }
  
  ##############
  #predictions #
  ##############
  
  #bases entières (entraînement et validation)
  ##################################################################################################
  ########### BASE DE TRAIN
  ###plann
  
  plannpred <- predictRS(plann.model, data = data_train, newtimes = newtimes, ratetable=slopop, 
                         age = "age", year = "year", sex = "sexchara")
  
  mean.plann <- plannpred$mpredictions$net_survival 
  
  ###flexnet
  
  ### modèle PH
  ##2 noeuds
  
  flexpred1.2 <- predict(flex.model1.2, newtimes = newtimes)$predictions
  
  mean.flex1.2 <- apply(flexpred1.2, FUN="mean", MARGIN=2)
  
  ##4 noeuds
  
  flexpred1.4 <- predict(flex.model1.4, newtimes = newtimes)$predictions
  
  mean.flex1.4 <- apply(flexpred1.4, FUN="mean", MARGIN=2)
  
  #### modèle NPH
  
  ##2 noeuds
  
  flexpred2.2 <- predict(flex.model2.2 , newtimes = newtimes)$predictions
  
  mean.flex2.2  <- apply(flexpred2.2 , FUN="mean", MARGIN=2)
  
  ##4 noeuds
  
  flexpred2.4 <- predict(flex.model2.4 , newtimes = newtimes)$predictions
  
  mean.flex2.4  <- apply(flexpred2.4 , FUN="mean", MARGIN=2)
  ### genweibull
  
  WGpred <- predict(WG.model, newtimes = newtimes)$predictions
  
  mean.WG <- apply(WGpred, FUN="mean", MARGIN=2)
  
  
  ### Estimateur de Pohar-Perme 
  
  PP <- summary( rs.surv(Surv(times, status) ~ 1, data = data_train,
                         ratetable = slopop, method = "pohar-perme", add.times = newtimes,
                         rmap = list(age = age, sex = sexchara, year = year)), times = newtimes, extend = TRUE) 
  
  poharpred <- PP$surv 
  
  ###### par strates ######
  ######
  #plann
  
  for(j in strata_names){
    
    data_train_var <- get(paste0("data_train_", j))
    
    plannpredS <- predictRS(plann.model, data = data_train_var, newtimes = newtimes, 
                            ratetable = slopop, age = "age", year = "year", sex = "sexchara")
    
    mean.plannS <- plannpredS$mpredictions$net_survival
    
    assign(paste0("plannpred_", j), plannpredS)
    assign(paste0("mean.plann_", j), mean.plannS)
  }
  
  ###flexnet
  
  ### modèle PH
  
  ## 2 noeuds
  
  for(j in strata_names){
    
    data_train_var <- get(paste0("data_train_", j))
    
    flexpred1.2S <- predict(flex.model1.2, newtimes = newtimes, newdata = data_train_var)$predictions
    
    mean.flex1.2S <- apply(flexpred1.2S, FUN="mean", MARGIN=2)
    
    assign(paste0("flexpred1.2_", j), flexpred1.2S)
    assign(paste0("mean.flex1.2_", j), mean.flex1.2S)
  }
  
  ## 4 noeuds
  
  for(j in strata_names){
    
    data_train_var <- get(paste0("data_train_", j))
    
    flexpred1.4S <- predict(flex.model1.4, newtimes = newtimes, newdata = data_train_var)$predictions
    
    mean.flex1.4S <- apply(flexpred1.4S, FUN="mean", MARGIN=2)
    
    assign(paste0("flexpred1.4_", j), flexpred1.4S)
    assign(paste0("mean.flex1.4_", j), mean.flex1.4S)
  }
  
  #### modèle NPH
  
  ##2 noeuds
  
  for(j in strata_names){
    
    data_train_var <- get(paste0("data_train_", j))
    
    flexpred2.2S <- predict(flex.model2.2, newtimes = newtimes, newdata = data_train_var)$predictions
    
    mean.flex2.2S <- apply(flexpred2.2S, FUN="mean", MARGIN=2)
    
    assign(paste0("flexpred2.2_", j), flexpred2.2S)
    assign(paste0("mean.flex2.2_", j), mean.flex2.2S)
  }
  
  ##4 noeuds
  
  for(j in strata_names){
    
    data_train_var <- get(paste0("data_train_", j))
    
    flexpred2.4S <- predict(flex.model2.4, newtimes = newtimes, newdata = data_train_var)$predictions
    
    mean.flex2.4S <- apply(flexpred2.4S, FUN="mean", MARGIN=2)
    
    assign(paste0("flexpred2.4_", j), flexpred2.4S)
    assign(paste0("mean.flex2.4_", j), mean.flex2.4S)
  }
  
  ### genweibull
  
  for(j in strata_names){
    
    data_train_var <- get(paste0("data_train_", j))
    
    WGpredS <- predict(WG.model, newtimes = newtimes, newdata = data_train_var)$predictions
    
    mean.WGS <- apply(WGpredS, FUN="mean", MARGIN=2)
    
    assign(paste0("WGpred_", j), WGpredS)
    assign(paste0("mean.WG_", j), mean.WGS)
  }
  
  
  ##Estimateur de Pohar-Perme 
  
  for(j in strata_names){
    
    data_train_var <- get(paste0("data_train_", j))
    
    PPS <- summary( rs.surv(Surv(times, status) ~ 1, data = data_train_var,
                            ratetable = slopop, method = "pohar-perme", add.times = newtimes,
                            rmap = list(age = age, sex = sexchara, year = year)), times = newtimes, extend = TRUE) 
    
    poharpredS <- PPS$surv
    
    assign(paste0("PP_", j), PPS)
    assign(paste0("poharpred_", j), poharpredS)
  }
  
  ###### 
  
  ############################################
  ########### BASE DE VALIDATION
  ###plann
  
  plannpredval <-  predictRS(plann.model, data = data_valid, newtimes = newtimes, ratetable=slopop, 
                             age = "age", year = "year", sex = "sexchara")
  
  mean.plannval <- plannpredval$mpredictions$net_survival 
  
  ###flexnet
  
  ### modèle PH
  
  ### 2 noeuds
  
  flexpredval1.2 <- predict(flex.model1.2, newtimes = newtimes, newdata = data_valid)$predictions
  
  mean.flexval1.2 <- apply(flexpredval1.2, FUN="mean", MARGIN=2)
  
  ### 4 noeuds
  
  flexpredval1.4 <- predict(flex.model1.4, newtimes = newtimes, newdata = data_valid)$predictions
  
  mean.flexval1.4 <- apply(flexpredval1.4, FUN="mean", MARGIN=2)
  
  #### modèle NPH
  
  ### 2 noeuds
  
  flexpredval2.2 <- predict(flex.model2.2, newtimes = newtimes, newdata = data_valid)$predictions
  
  mean.flexval2.2 <- apply(flexpredval2.2, FUN="mean", MARGIN=2)
  
  ### 4 noeuds
  
  flexpredval2.4 <- predict(flex.model2.4, newtimes = newtimes, newdata = data_valid)$predictions
  
  mean.flexval2.4 <- apply(flexpredval2.4, FUN="mean", MARGIN=2)
  
  ### genweibull
  
  WGpredval <- predict(WG.model, newtimes = newtimes, newdata = data_valid)$predictions
  
  mean.WGval <- apply(WGpredval, FUN="mean", MARGIN=2)
  
  ### Estimateur de Pohar-Perme 
  
  PPval <- summary( rs.surv(Surv(times, status) ~ 1, data = data_valid,
                            ratetable = slopop, method = "pohar-perme", add.times = newtimes,
                            rmap = list(age = age, sex = sexchara, year = year)), times = newtimes, extend = TRUE) 
  
  poharpredval <- PPval$surv 
  
  ###### par strates ######
  
  
  ###plann
  
  for(j in strata_names){
    
    data_valid_var <- get(paste0("data_valid_", j))
    
    plannpredvalS <- predictRS(plann.model, data = data_valid_var, newtimes = newtimes, 
                               ratetable = slopop, age = "age", year = "year", sex = "sexchara")
    
    mean.plannvalS <- plannpredvalS$mpredictions$net_survival
    
    assign(paste0("plannpredval_", j), plannpredvalS)
    assign(paste0("mean.plannval_", j), mean.plannvalS)
  }
  
  ###flexnet
  
  ### modèle PH
  
  ### 2 noeuds 
  
  for(j in strata_names){
    
    data_valid_var <- get(paste0("data_valid_", j))
    
    flexpredval1.2S <- predict(flex.model1.2, newtimes = newtimes, newdata = data_valid_var)$predictions
    
    mean.flexval1.2S <- apply(flexpredval1.2S, FUN="mean", MARGIN=2)
    
    assign(paste0("flexpredval1.2_", j), flexpredval1.2S)
    assign(paste0("mean.flexval1.2_", j), mean.flexval1.2S)
  }
  
  ### 4 noeuds 
  
  for(j in strata_names){
    
    data_valid_var <- get(paste0("data_valid_", j))
    
    flexpredval1.4S <- predict(flex.model1.4, newtimes = newtimes, newdata = data_valid_var)$predictions
    
    mean.flexval1.4S <- apply(flexpredval1.4S, FUN="mean", MARGIN=2)
    
    assign(paste0("flexpredval1.4_", j), flexpredval1.4S)
    assign(paste0("mean.flexval1.4_", j), mean.flexval1.4S)
  }
  #### modèle NPH
  
  ## 2 noeuds
  
  for(j in strata_names){
    
    data_valid_var <- get(paste0("data_valid_", j))
    
    flexpredval2.2S <- predict(flex.model2.2, newtimes = newtimes, newdata = data_valid_var)$predictions
    
    mean.flexval2.2S <- apply(flexpredval2.2S, FUN="mean", MARGIN=2)
    
    assign(paste0("flexpredval2.2_", j), flexpredval2.2S)
    assign(paste0("mean.flexval2.2_", j), mean.flexval2.2S)
  }
  
  ## 4 noeuds
  
  for(j in strata_names){
    
    data_valid_var <- get(paste0("data_valid_", j))
    
    flexpredval2.4S <- predict(flex.model2.4, newtimes = newtimes, newdata = data_valid_var)$predictions
    
    mean.flexval2.4S <- apply(flexpredval2.4S, FUN="mean", MARGIN=2)
    
    assign(paste0("flexpredval2.4_", j), flexpredval2.4S)
    assign(paste0("mean.flexval2.4_", j), mean.flexval2.4S)
  }
  
  ### genweibull
  
  for(j in strata_names){
    
    data_valid_var <- get(paste0("data_valid_", j))
    
    WGpredvalS <- predict(WG.model, newtimes = newtimes, newdata = data_valid_var)$predictions
    
    mean.WGvalS <- apply(WGpredvalS, FUN="mean", MARGIN=2)
    
    assign(paste0("WGpredval_", j), WGpredvalS)
    assign(paste0("mean.WGval_", j), mean.WGvalS)
  }
  
  ##Estimateur de Pohar-Perme 
  
  for(j in strata_names){
    
    data_valid_var <- get(paste0("data_valid_", j))
    
    PPvalS <- summary( rs.surv(Surv(times, status) ~ 1, data = data_valid_var,
                               ratetable = slopop, method = "pohar-perme", add.times = newtimes,
                               rmap = list(age = age, sex = sexchara, year = year)), times = newtimes, extend = TRUE) 
    
    poharpredvalS <- PPvalS$surv
    
    assign(paste0("PPval_", j), PPvalS)
    assign(paste0("poharpredval_", j), poharpredvalS)
  }
  ##################################################################################################
  
  ######
  #time
  ######
  
  time <- Sys.time()-start
  time <- as.numeric(time, units = 'mins')
  ########################
  #       outputs        #
  ########################
  
  #####bases####
  
  ### TRAIN
  write.table(data_train,  paste0(path0, "DATAFRAMES/TRAIN/WHOLE/",i,"_dftrain.csv"), sep = ";", row.names = F, col.names = T)
  ##strates
  for(j in strata_names){
    
    data_train_var <- get(paste0("data_train_", j))
    
    file_path <- paste0(path0, "DATAFRAMES/TRAIN/", j, "/", i, "_dftrain_", j, ".csv")
    
    write.table(data_train_var, file = file_path, sep = ";", row.names = FALSE, col.names = TRUE)
  }
  ### VALID
  write.table(data_valid,  paste0(path0, "DATAFRAMES/VALID/WHOLE/",i,"_dfvalid.csv"), sep = ";", row.names = F, col.names = T)
  ##strates
  
  for(j in strata_names){
    
    data_valid_var <- get(paste0("data_valid_", j))
    
    file_path <- paste0(path0, "DATAFRAMES/VALID/", j, "/", i, "_dfvalid_", j, ".csv")
    
    write.table(data_valid_var, file = file_path, sep = ";", row.names = FALSE, col.names = TRUE)
  }
  
  ############## paramètres estimés par CV ###########
  
  m1.2value <- paste0(m1.2, failedindic1.2)
  m1.4value <- paste0(m1.4, failedindic1.4)
  
  m2.2value <- paste0(m2.2, failedindic2.2)
  m2.4value <- paste0(m2.4, failedindic2.4)
  
  paramsCV <- as.data.frame(c(inter = inter,size = size,decay = decay, maxit = maxit,
                              MaxNWts = MaxNWts, m1.2 = m1.2value, m1.4 = m1.4value, 
                              m2.2 = m2.2value, m2.4 = m2.4value, logcoeffWG = logcoeffWG))
  
  write.table(paramsCV,  paste0(path0, "PARAM/",i,"_param.csv"))
  
  
  ############# temps écoulé total & par CV de chaque modèle et taille de la base #########################
  
  timeComputation <- as.data.frame(c(time = time, plann = TimeCVplann, N = dim(data)[1]))
  
  write.table(timeComputation,  paste0(path0, "TIME/",i,"_time.csv"))
  
  ##############Log(LIKELIHOOD)####################### 
  
  plann.loglik <- plannpred$loglik  
  
  flex1.2.loglik <- unname(flex.model1.2$loglik[1])
  flex1.4.loglik <- unname(flex.model1.4$loglik[1])
  
  flex2.2.loglik <- unname(flex.model2.2$loglik[1])
  flex2.4.loglik <- unname(flex.model2.4$loglik[1])
  
  WG.loglik <- unname(WG.model$loglik[1])
  
  logliks <- as.data.frame(c(plann = plann.loglik, flex1.2 = flex1.2.loglik,
                             flex1.4 = flex1.4.loglik, flex2.2 = flex2.2.loglik,
                             flex2.4 = flex2.4.loglik, WG = WG.loglik))
  
  write.table(logliks,  paste0(path0, "LOGLIK/",i,"_loglik.csv"))
  
  ##################DATA TRAIN
  #################################### MOYENNES
  #######THEO
  
  mean_theo <- t(data.frame(THEO = mean.theo))
  write.table(mean_theo, paste0(path0, "TRAIN/WHOLE/THEO/mean_survival/",i,"_mean_theo.csv"), sep = ";", row.names = F, col.names = T)
  
  
  ######PLANN
  
  mean_PLANN <- t(data.frame(PLANN = mean.plann))
  write.table(mean_PLANN, paste0(path0, "TRAIN/WHOLE/PLANN/mean_survival/",i,"_mean_PLANN.csv"), sep = ";", row.names = F, col.names = T)
  
  #######FLEX1
  
  ##2 noeuds
  mean_flex1.2 <- t(data.frame(flex1.2=mean.flex1.2))
  write.table(mean_flex1.2, paste0(path0, "TRAIN/WHOLE/FLEX1.2/mean_survival/",i,"_mean_flex12.csv"), sep = ";", row.names = F, col.names = T)
  
  ##4 noeuds
  mean_flex1.4 <- t(data.frame(flex1.4=mean.flex1.4))
  write.table(mean_flex1.4, paste0(path0, "TRAIN/WHOLE/FLEX1.4/mean_survival/",i,"_mean_flex14.csv"), sep = ";", row.names = F, col.names = T)
  
  ######FLEX2
  
  ##2 noeuds
  mean_flex2.2 <- t(data.frame(flex2.2=mean.flex2.2))
  write.table(mean_flex2.2, paste0(path0, "TRAIN/WHOLE/FLEX2.2/mean_survival/",i,"_mean_flex22.csv"), sep = ";", row.names = F, col.names = T)
  
  ##4 noeuds
  mean_flex2.4 <- t(data.frame(flex2.4=mean.flex2.4))
  write.table(mean_flex2.4, paste0(path0, "TRAIN/WHOLE/FLEX2.4/mean_survival/",i,"_mean_flex24.csv"), sep = ";", row.names = F, col.names = T)
  
  ######WG
  
  mean_flexWG <- t(data.frame(flexWG=mean.WG))
  write.table(mean_flexWG, paste0(path0, "TRAIN/WHOLE/WG/mean_survival/",i,"_mean_flexWG.csv"), sep = ";",row.names = F, col.names = T)
  
  ######Estimateur de Pohar-Perme 
  
  write.table(poharpred, paste0(path0, "TRAIN/WHOLE/PP/",i,"PP.csv"), sep = ";",row.names = T, col.names = T)
  ###strates 
  ######
  
  ##theo
  
  mean.theo_STRATES = c()
  
  for (j in strata_names) {
    
    mean.theo_STRATES <- rbind(mean.theo_STRATES, get(paste0("mean.theo_", j) ) )
    rownames(mean.theo_STRATES)[nrow(mean.theo_STRATES)] <- paste0("mean.theo_", j)
  }
  
  mean.theo_STRATES <- data.frame(strat = rownames(mean.theo_STRATES), mean.theo_STRATES)
  write.table(mean.theo_STRATES, paste0(path0, "TRAIN/STRATA/THEO/mean_survival/",i,"mean_theo_strates.csv"), sep = ";",row.names = F, col.names = T)
  
  #PLANN
  
  mean.plann_STRATES = c()
  
  for (j in strata_names) {
    
    mean.plann_STRATES <- rbind(mean.plann_STRATES, get(paste0("mean.plann_", j) ) )
    rownames(mean.plann_STRATES)[nrow(mean.plann_STRATES)] <- paste0("mean.plann_", j)
  }
  
  mean.plann_STRATES <- data.frame(strat = rownames(mean.plann_STRATES), mean.plann_STRATES)
  write.table(mean.plann_STRATES, paste0(path0, "TRAIN/STRATA/PLANN/mean_survival/",i,"mean_plann_strates.csv"), sep = ";",row.names = F, col.names = T)
  
  
  #flex1
  
  ##2 noeuds
  mean.flex1.2_STRATES = c()
  
  for (j in strata_names) {
    
    mean.flex1.2_STRATES <- rbind(mean.flex1.2_STRATES, get(paste0("mean.flex1.2_", j) ) )
    rownames(mean.flex1.2_STRATES)[nrow(mean.flex1.2_STRATES)] <- paste0("mean.flex1.2_", j)
  }
  
  mean.flex1.2_STRATES <- data.frame(strat = rownames(mean.flex1.2_STRATES), mean.flex1.2_STRATES)
  write.table(mean.flex1.2_STRATES, paste0(path0, "TRAIN/STRATA/FLEX1.2/mean_survival/",i,"mean_flex12_strates.csv"), sep = ";",row.names = F, col.names = T)
  
  ##4 noeuds
  
  mean.flex1.4_STRATES = c()
  
  for (j in strata_names) {
    
    mean.flex1.4_STRATES <- rbind(mean.flex1.4_STRATES, get(paste0("mean.flex1.4_", j) ) )
    rownames(mean.flex1.4_STRATES)[nrow(mean.flex1.4_STRATES)] <- paste0("mean.flex1.4_", j)
  }
  
  mean.flex1.4_STRATES <- data.frame(strat = rownames(mean.flex1.4_STRATES), mean.flex1.4_STRATES)
  write.table(mean.flex1.4_STRATES, paste0(path0, "TRAIN/STRATA/FLEX1.4/mean_survival/",i,"mean_flex14_strates.csv"), sep = ";",row.names = F, col.names = T)
  
  #flex2
  
  ## 2 noeuds
  
  mean.flex2.2_STRATES = c()
  
  for (j in strata_names) {
    
    mean.flex2.2_STRATES <- rbind(mean.flex2.2_STRATES, get(paste0("mean.flex2.2_", j) ) )
    rownames(mean.flex2.2_STRATES)[nrow(mean.flex2.2_STRATES)] <- paste0("mean.flex2.2_", j)
  }
  
  mean.flex2.2_STRATES <- data.frame(strat = rownames(mean.flex2.2_STRATES), mean.flex2.2_STRATES)
  write.table(mean.flex2.2_STRATES, paste0(path0, "TRAIN/STRATA/FLEX2.2/mean_survival/",i,"mean_flex22_strates.csv"), sep = ";",row.names = F, col.names = T)
  
  ## 4 noeuds
  
  mean.flex2.4_STRATES = c()
  
  for (j in strata_names) {
    
    mean.flex2.4_STRATES <- rbind(mean.flex2.4_STRATES, get(paste0("mean.flex2.4_", j) ) )
    rownames(mean.flex2.4_STRATES)[nrow(mean.flex2.4_STRATES)] <- paste0("mean.flex2.4_", j)
  }
  
  mean.flex2.4_STRATES <- data.frame(strat = rownames(mean.flex2.4_STRATES), mean.flex2.4_STRATES)
  write.table(mean.flex2.4_STRATES, paste0(path0, "TRAIN/STRATA/FLEX2.4/mean_survival/",i,"mean_flex24_strates.csv"), sep = ";",row.names = F, col.names = T)
  
  #WG
  mean.WG_STRATES = c()
  
  for (j in strata_names) {
    
    mean.WG_STRATES <- rbind(mean.WG_STRATES, get(paste0("mean.WG_", j) ) )
    rownames(mean.WG_STRATES)[nrow(mean.WG_STRATES)] <- paste0("mean.WG_", j)
  }
  
  mean.WG_STRATES <- data.frame(strat = rownames(mean.WG_STRATES), mean.WG_STRATES)
  write.table(mean.WG_STRATES, paste0(path0, "TRAIN/STRATA/WG/mean_survival/",i,"mean_WG_strates.csv"), sep = ";",row.names = F, col.names = T)
  
  
  ####Estimateur de Pohar-Perme
  
  poharpred_strates = c()
  
  for (j in strata_names) {
    
    poharpred_strates <- rbind(poharpred_strates, get(paste0("poharpred_", j) ) )
    rownames(poharpred_strates)[nrow(poharpred_strates)] <- paste0("poharpred_", j)
  }
  
  poharpred_strates <- data.frame(strat = rownames(poharpred_strates), poharpred_strates)
  write.table(poharpred_strates, paste0(path0, "TRAIN/STRATA/PP/",i,"PP_strates.csv"), sep = ";",row.names = F, col.names = T)
  
  ######
  
  #################################### SURVIES INDIVIDUELLES
  
  #######THEO
  
  ind_theo <- data.frame(theopred)
  write.table(ind_theo, paste0(path0, "TRAIN/WHOLE/THEO/individual_survival/",
                               i,"_ind_theo.csv"), sep = ";", row.names = F, col.names = T)
  
  ######PLANN
  
  ind_PLANN <- data.frame(plannpred$ipredictions$relative_survival)
  write.table(ind_PLANN, paste0(path0, "TRAIN/WHOLE/PLANN/individual_survival/",
                                i,"_ind_PLANN.csv"), sep = ";", row.names = F, col.names = T)
  
  #######FLEX1
  
  ##2 noeuds
  ind_flex1.2 <- data.frame(flexpred1.2)
  write.table(ind_flex1.2, paste0(path0, "TRAIN/WHOLE/FLEX1.2/individual_survival/"
                                  ,i,"_ind_flex12.csv"), sep = ";", row.names = F, col.names = T)
  
  ##4 noeuds
  ind_flex1.4 <- data.frame(flexpred1.4)
  write.table(ind_flex1.4, paste0(path0, "TRAIN/WHOLE/FLEX1.4/individual_survival/"
                                  ,i,"_ind_flex14.csv"), sep = ";", row.names = F, col.names = T)
  ######FLEX2
  
  ##2 noeuds
  ind_flex2.2 <- data.frame(flexpred2.2)
  write.table(ind_flex2.2, paste0(path0, "TRAIN/WHOLE/FLEX2.2/individual_survival/"
                                  ,i,"_ind_flex22.csv"), sep = ";", row.names = F, col.names = T)
  
  ##4 noeuds
  ind_flex2.4 <- data.frame(flexpred2.4)
  write.table(ind_flex2.4, paste0(path0, "TRAIN/WHOLE/FLEX2.4/individual_survival/"
                                  ,i,"_ind_flex24.csv"), sep = ";", row.names = F, col.names = T)
  
  ######WG
  
  ind_WG <- data.frame(WGpred)
  write.table(ind_WG, paste0(path0, "TRAIN/WHOLE/WG/individual_survival/"
                             ,i,"_ind_WG.csv"), sep = ";",row.names = F, col.names = T)
  
  
  ###strates 
  ######
  
  #theo
  
  for (j in strata_names) {
    
    ind_theo <- data.frame(get(paste0("theopred_", j)))
    
    write.table(ind_theo, paste0(path0, "TRAIN/STRATA/THEO/individual_survival/", i, "_ind_theo_", j, ".csv"),
                sep = ";",row.names = F, col.names = T)
  }
  
  ##plann
  
  for (j in strata_names) {
    
    ind_PLANN <- data.frame(get(paste0("plannpred_", j))$ipredictions$relative_survival)
    
    write.table(ind_PLANN, paste0(path0, "TRAIN/STRATA/PLANN/individual_survival/", i, "_ind_PLANN_", j, ".csv"),
                sep = ";",row.names = F, col.names = T)
  }
  
  #flex1
  
  ## 2 noeuds
  for (j in strata_names) {
    
    ind_flex1.2 <- data.frame(get(paste0("flexpred1.2_", j)))
    
    write.table(ind_flex1.2, paste0(path0, "TRAIN/STRATA/FLEX1.2/individual_survival/", i, "_ind_flex12_", j, ".csv"),
                sep = ";",row.names = F, col.names = T)
  }
  
  ## 4 noeuds
  for (j in strata_names) {
    
    ind_flex1.4 <- data.frame(get(paste0("flexpred1.4_", j)))
    
    write.table(ind_flex1.4, paste0(path0, "TRAIN/STRATA/FLEX1.4/individual_survival/", i, "_ind_flex14_", j, ".csv"),
                sep = ";",row.names = F, col.names = T)
  }
  
  #flex2
  
  ##2 noeuds
  for (j in strata_names) {
    
    ind_flex2.2 <- data.frame(get(paste0("flexpred2.2_", j)))
    
    write.table(ind_flex2.2, paste0(path0, "TRAIN/STRATA/FLEX2.2/individual_survival/", i, "_ind_flex22_", j, ".csv"),
                sep = ";",row.names = F, col.names = T)
  }
  
  ##4 noeuds
  for (j in strata_names) {
    
    ind_flex2.4 <- data.frame(get(paste0("flexpred2.4_", j)))
    
    write.table(ind_flex2.4, paste0(path0, "TRAIN/STRATA/FLEX2.4/individual_survival/", i, "_ind_flex24_", j, ".csv"),
                sep = ";",row.names = F, col.names = T)
  }
  #WG
  
  for (j in strata_names) {
    
    ind_WG <- data.frame(get(paste0("WGpred_", j)))
    
    write.table(ind_WG, paste0(path0, "TRAIN/STRATA/WG/individual_survival/", i, "_ind_WG_", j, ".csv"),
                sep = ";",row.names = F, col.names = T)
  }
  
  ######
  
  ##################DATA VALID
  #################################### MOYENNES
  
  #######THEO
  
  mean_theoval <- t(data.frame(THEO = mean.theoval))
  write.table(mean_theoval, paste0(path0, "VALID/WHOLE/THEO/mean_survival/",
                                   i,"_mean_theoval.csv"), sep = ";", row.names = F, col.names = T)
  
  ######PLANN
  
  mean_PLANNval <- t(data.frame(PLANN = mean.plannval))
  write.table(mean_PLANNval, paste0(path0, "VALID/WHOLE/PLANN/mean_survival/",
                                    i,"_mean_PLANNval.csv"), sep = ";", row.names = F, col.names = T)
  
  #######FLEX1
  
  ##2 noeuds
  mean_flexval1.2 <- t(data.frame(flex1.2=mean.flexval1.2))
  write.table(mean_flexval1.2, paste0(path0, "VALID/WHOLE/FLEX1.2/mean_survival/"
                                      ,i,"_mean_flexval12.csv"), sep = ";", row.names = F, col.names = T)
  ##4 noeuds
  mean_flexval1.4 <- t(data.frame(flex1.4=mean.flexval1.4))
  write.table(mean_flexval1.4, paste0(path0, "VALID/WHOLE/FLEX1.4/mean_survival/"
                                      ,i,"_mean_flexval14.csv"), sep = ";", row.names = F, col.names = T)
  
  ######FLEX2
  
  ## 2 noeuds
  mean_flexval2.2 <- t(data.frame(flex2.2=mean.flexval2.2))
  write.table(mean_flexval2.2, paste0(path0, "VALID/WHOLE/FLEX2.2/mean_survival/"
                                      ,i,"_mean_flexval22.csv"), sep = ";", row.names = F, col.names = T)
  
  ## 4 noeuds
  mean_flexval2.4 <- t(data.frame(flex2.4=mean.flexval2.4))
  write.table(mean_flexval2.4, paste0(path0, "VALID/WHOLE/FLEX2.4/mean_survival/"
                                      ,i,"_mean_flexval24.csv"), sep = ";", row.names = F, col.names = T)
  ######WG
  
  mean_flexWGval <- t(data.frame(flexWG=mean.WGval))
  write.table(mean_flexWGval, paste0(path0, "VALID/WHOLE/WG/mean_survival/"
                                     ,i,"_mean_flexWGval.csv"), sep = ";",row.names = F, col.names = T)
  
  ###### Estimateur de pohar perme
  
  write.table(poharpred, paste0(path0, "VALID/WHOLE/PP/",i,"PPval.csv"), sep = ";",row.names = T, col.names = T)
  
  ###strates 
  ######
  
  ##theo
  
  mean.theoval_STRATES <- c()
  
  for (j in strata_names) {
    
    mean.theoval_STRATES <- rbind(mean.theoval_STRATES, get(paste0("mean.theoval_", j) ) )
    rownames(mean.theoval_STRATES)[nrow(mean.theoval_STRATES)] <- paste0("mean.theoval_", j)
  }
  
  mean.theoval_STRATES <- data.frame(strat = rownames(mean.theo_STRATES), mean.theo_STRATES)
  write.table(mean.theoval_STRATES, paste0(path0, "VALID/STRATA/THEO/mean_survival/",i,"mean_theoval_strates.csv"), sep = ";",row.names = F, col.names = T)
  
  
  #PLANN
  
  mean.plannval_STRATES = c()
  
  for (j in strata_names) {
    
    mean.plannval_STRATES <- rbind(mean.plannval_STRATES, get(paste0("mean.plannval_", j) ) )
    rownames(mean.plannval_STRATES)[nrow(mean.plannval_STRATES)] <- paste0("mean.plannval_", j)
  }
  
  mean.plannval_STRATES <- data.frame(strat = rownames(mean.plannval_STRATES), mean.plannval_STRATES)
  write.table(mean.plannval_STRATES, paste0(path0, "VALID/STRATA/PLANN/mean_survival/",i,"mean_plannval_strates.csv"), sep = ";",row.names = F, col.names = T)
  
  #flex1
  ## 2 noeuds
  mean.flexval1.2_STRATES = c()
  
  for (j in strata_names) {
    
    mean.flexval1.2_STRATES <- rbind(mean.flexval1.2_STRATES, get(paste0("mean.flexval1.2_", j) ) )
    rownames(mean.flexval1.2_STRATES)[nrow(mean.flexval1.2_STRATES)] <- paste0("mean.flexval1.2_", j)
  }
  
  mean.flexval1.2_STRATES <- data.frame(strat = rownames(mean.flexval1.2_STRATES), mean.flexval1.2_STRATES)
  write.table(mean.flexval1.2_STRATES, paste0(path0, "VALID/STRATA/FLEX1.2/mean_survival/",i,"mean_flexval12_strates.csv"), sep = ";",row.names = F, col.names = T)
  
  ##4 noeuds
  mean.flexval1.4_STRATES = c()
  
  for (j in strata_names) {
    
    mean.flexval1.4_STRATES <- rbind(mean.flexval1.4_STRATES, get(paste0("mean.flexval1.4_", j) ) )
    rownames(mean.flexval1.4_STRATES)[nrow(mean.flexval1.4_STRATES)] <- paste0("mean.flexval1.4_", j)
  }
  
  mean.flexval1.4_STRATES <- data.frame(strat = rownames(mean.flexval1.4_STRATES), mean.flexval1.4_STRATES)
  write.table(mean.flexval1.4_STRATES, paste0(path0, "VALID/STRATA/FLEX1.4/mean_survival/",i,"mean_flexval14_strates.csv"), sep = ";",row.names = F, col.names = T)
  
  
  #flex2
  
  ##2 noeuds
  mean.flexval2.2_STRATES = c()
  
  for (j in strata_names) {
    
    mean.flexval2.2_STRATES <- rbind(mean.flexval2.2_STRATES, get(paste0("mean.flexval2.2_", j) ) )
    rownames(mean.flexval2.2_STRATES)[nrow(mean.flexval2.2_STRATES)] <- paste0("mean.flexval2.2_", j)
  }
  
  mean.flexval2.2_STRATES <- data.frame(strat = rownames(mean.flexval2.2_STRATES), mean.flexval2.2_STRATES)
  write.table(mean.flexval2.2_STRATES, paste0(path0, "VALID/STRATA/FLEX2.2/mean_survival/",i,"mean_flexval22_strates.csv"), sep = ";",row.names = F, col.names = T)
  
  ##4 noeuds
  mean.flexval2.4_STRATES = c()
  
  for (j in strata_names) {
    
    mean.flexval2.4_STRATES <- rbind(mean.flexval2.4_STRATES, get(paste0("mean.flexval2.4_", j) ) )
    rownames(mean.flexval2.4_STRATES)[nrow(mean.flexval2.4_STRATES)] <- paste0("mean.flexval2.4_", j)
  }
  
  mean.flexval2.4_STRATES <- data.frame(strat = rownames(mean.flexval2.4_STRATES), mean.flexval2.4_STRATES)
  write.table(mean.flexval2.4_STRATES, paste0(path0, "VALID/STRATA/FLEX2.4/mean_survival/",i,"mean_flexval24_strates.csv"), sep = ";",row.names = F, col.names = T)
  
  #WG
  mean.WGval_STRATES = c()
  
  for (j in strata_names) {
    
    mean.WGval_STRATES <- rbind(mean.WGval_STRATES, get(paste0("mean.WGval_", j) ) )
    rownames(mean.WGval_STRATES)[nrow(mean.WGval_STRATES)] <- paste0("mean.WGval_", j)
  }
  
  mean.WGval_STRATES <- data.frame(strat = rownames(mean.WGval_STRATES), mean.WGval_STRATES)
  write.table(mean.WGval_STRATES, paste0(path0, "VALID/STRATA/WG/mean_survival/",i,"mean_WGval_strates.csv"), sep = ";",row.names = F, col.names = T)
  
  ####Estimateur de Pohar-Perme
  
  poharpredval_strates = c()
  
  for (j in strata_names) {
    
    poharpredval_strates <- rbind(poharpredval_strates, get(paste0("poharpredval_", j) ) )
    rownames(poharpredval_strates)[nrow(poharpredval_strates)] <- paste0("poharpredval_", j)
  }
  
  poharpredval_strates <- data.frame(strat = rownames(poharpredval_strates), poharpredval_strates)
  write.table(poharpredval_strates, paste0(path0, "VALID/STRATA/PP/",i,"PPval_strates.csv"), sep = ";",row.names = F, col.names = T)
  
  ######
  #################################### SURVIES INDIVIDUELLES
  
  #######THEO
  ind_theoval <- data.frame(theopredval)
  write.table(ind_theoval, 
              paste0(path0, "VALID/WHOLE/THEO/individual_survival/",i
                     ,"_ind_theoval.csv"), sep = ";", row.names = F, col.names = T)
  
  ######PLANN
  
  ind_PLANNval <- data.frame(plannpredval$ipredictions$relative_survival)
  write.table(ind_PLANNval, paste0(path0, "VALID/WHOLE/PLANN/individual_survival/",
                                   i,"_ind_PLANNval.csv"), sep = ";", row.names = F, col.names = T)
  
  #######FLEX1
  
  ##2 noeuds
  
  ind_flexval1.2 <- data.frame(flexpredval1.2)
  write.table(ind_flexval1.2, paste0(path0, "VALID/WHOLE/FLEX1.2/individual_survival/"
                                     ,i,"_ind_flexval12.csv"), sep = ";", row.names = F, col.names = T)
  
  ##4 noeuds
  
  ind_flexval1.4 <- data.frame(flexpredval1.4)
  write.table(ind_flexval1.4, paste0(path0, "VALID/WHOLE/FLEX1.4/individual_survival/"
                                     ,i,"_ind_flexval14.csv"), sep = ";", row.names = F, col.names = T)
  ######FLEX2
  
  ##2 noeuds
  ind_flexval2.2 <- data.frame(flexpredval2.2)
  write.table(ind_flexval2.2, paste0(path0, "VALID/WHOLE/FLEX2.2/individual_survival/"
                                     ,i,"_ind_flexval22.csv"), sep = ";", row.names = F, col.names = T)
  
  ##4 noeuds
  ind_flexval2.4 <- data.frame(flexpredval2.4)
  write.table(ind_flexval2.4, paste0(path0, "VALID/WHOLE/FLEX2.4/individual_survival/"
                                     ,i,"_ind_flexval24.csv"), sep = ";", row.names = F, col.names = T)
  ######WG
  
  ind_flexWGval <- data.frame(WGpredval)
  write.table(ind_flexWGval, paste0(path0, "VALID/WHOLE/WG/individual_survival/"
                                    ,i,"_ind_flexWGval.csv"), sep = ";",row.names = F, col.names = T)
  
  
  ###strates ####
  ######
  
  #theo
  
  for (j in strata_names) {
    
    ind_theoval <- data.frame(get(paste0("theopredval_", j)))
    
    write.table(ind_theoval, paste0(path0, "VALID/STRATA/THEO/individual_survival/", i, "_ind_theoval_", j, ".csv"),
                sep = ";",row.names = F, col.names = T)
  }
  
  ##plann
  
  for (j in strata_names) {
    
    ind_PLANNval <- data.frame(get(paste0("plannpredval_", j))$ipredictions$relative_survival)
    
    write.table(ind_PLANNval, paste0(path0, "VALID/STRATA/PLANN/individual_survival/", i, "_ind_PLANNval_", j, ".csv"),
                sep = ";",row.names = F, col.names = T)
  }
  
  #flex1
  
  ##2 noeuds
  for (j in strata_names) {
    
    ind_flexval1.2 <- data.frame(get(paste0("flexpredval1.2_", j)))
    
    write.table(ind_flexval1.2, paste0(path0, "VALID/STRATA/FLEX1.2/individual_survival/", i, "_ind_flexval12_", j, ".csv"),
                sep = ";",row.names = F, col.names = T)
  }
  ##4 noeuds
  for (j in strata_names) {
    
    ind_flexval1.4 <- data.frame(get(paste0("flexpredval1.4_", j)))
    
    write.table(ind_flexval1.4, paste0(path0, "VALID/STRATA/FLEX1.4/individual_survival/", i, "_ind_flexval14_", j, ".csv"),
                sep = ";",row.names = F, col.names = T)
  }
  #flex2
  
  ##2 noeuds
  for (j in strata_names) {
    
    ind_flexval2.2 <- data.frame(get(paste0("flexpredval2.2_", j)))
    
    write.table(ind_flexval2.2, paste0(path0, "VALID/STRATA/FLEX2.2/individual_survival/", i, "_ind_flexval22_", j, ".csv"),
                sep = ";",row.names = F, col.names = T)
  }
  
  ##4 noeuds
  for (j in strata_names) {
    
    ind_flexval2.4 <- data.frame(get(paste0("flexpredval2.4_", j)))
    
    write.table(ind_flexval2.4, paste0(path0, "VALID/STRATA/FLEX2.4/individual_survival/", i, "_ind_flexval24_", j, ".csv"),
                sep = ";",row.names = F, col.names = T)
  }
  
  #WG
  
  for (j in strata_names) {
    
    ind_WGval <- data.frame(get(paste0("WGpredval_", j)))
    
    write.table(ind_WGval, paste0(path0, "VALID/STRATA/WG/individual_survival/", i, "_ind_WGval_", j, ".csv"),
                sep = ";",row.names = F, col.names = T)
  }
  ######
  
}
#fin itération


iterations <- 1:1000

mclapply(iterations, simulate_iteration, mc.cores = detectCores() - 14)

max_attempts <- 3
attempt <- 1
while (TRUE) {
  indic <- c()  
  
  for (i in 1:1000) {
    if (!file.exists(paste0(path0, "DATAFRAMES/TRAIN/WHOLE/", i, "_dftrain.csv"))) {
      indic <- c(indic, i)
    }
  }
  
  if (length(indic) == 0) {
    break
  }
  if (attempt > max_attempts) {
    break
  }
  mclapply(indic, simulate_iteration, mc.cores = detectCores() - 14)
  attempt <- attempt + 1
}
# nb.cluster = detectCores() - 4   
# cl <- makeCluster(nb.cluster, type="FORK") 
# registerDoParallel(cl)
# 
# foreach(i = 1:iterations,  .inorder = FALSE, .verbose = T) %dopar% 
#   { simulate_iteration(i) }
# stopCluster(cl)


rm(list = ls(envir = .GlobalEnv), envir = .GlobalEnv)
### 3000SiEND

########### CALC_INDICATEURS

### 3000IndSTA
library(parallel)
library(doParallel)

path0 <- "~/Documents/Simulations/Simulations mai 2025/N3000/Output simulations_N3000_NPH_HC/"

indic = c()
for(i in 1:1000){
  if(!file.exists(paste0(path0,"DATAFRAMES/TRAIN/WHOLE/",i,"_dftrain.csv"))){ indic <<- c(indic, i)} 
}


if(!is.null(indic)){
  iterations <- c(1:1000)[-indic]
}else{
  iterations <- c(1:1000)
}
strata_names <- c("HC", "HR", "FC", "FR")
newtimes <- c(365.241, 365.241*3, 365.241*5, 365.241*10) 

### FIN INDICATEURS ~ligne 3113

start <- Sys.time()

### Calculs des indicateurs

### RMSE & BIAIS
############################################ TRAIN #################################################
#WHOLE

###on récupère les valeurs théoriques moyennes Se(t) 
ALL_theo <- data.frame()

for(k in iterations){
  assign(paste0("meanTrainTheo_", k) , read.csv(paste0(path0, "TRAIN/WHOLE/THEO/mean_survival/",k,"_mean_theo.csv"),
                                                sep = ";") )
  
  ALL_theo <-rbind(ALL_theo, get(paste0("meanTrainTheo_", k)))
  rm(list =paste0("meanTrainTheo_", k))
}

#moyenne théorique de la valur de survie
S_theo <- colMeans(ALL_theo)
############################################
#importation des données des différents modèles 
#########
##plann
ALL_PLANN <- data.frame() 
for(k in iterations){
  assign(paste0("meanTrainPLANN_", k) , read.csv(paste0(path0, "TRAIN/WHOLE/PLANN/mean_survival/",k,"_mean_PLANN.csv"), sep = ";") )
  ALL_PLANN <-rbind(ALL_PLANN, get(paste0("meanTrainPLANN_", k)))
  rm(list =paste0("meanTrainPLANN_", k))
}

###flex1
##2 noeuds
ALL_flex1.2 <- data.frame()
for(k in iterations){
  assign(paste0("meanTrainflex1.2_", k) , read.csv(paste0(path0, "TRAIN/WHOLE/FLEX1.2/mean_survival/",k,"_mean_flex12.csv"), sep = ";") )
  ALL_flex1.2 <- rbind(ALL_flex1.2, get(paste0("meanTrainflex1.2_",k)))
  rm(list =paste0("meanTrainflex1.2_", k))
}

##4 noeuds
ALL_flex1.4 <- data.frame()
for(k in iterations){
  assign(paste0("meanTrainflex1.4_", k) , read.csv(paste0(path0, "TRAIN/WHOLE/FLEX1.4/mean_survival/",k,"_mean_flex14.csv"), sep = ";") )
  ALL_flex1.4 <- rbind(ALL_flex1.4, get(paste0("meanTrainflex1.4_",k)))
  rm(list =paste0("meanTrainflex1.4_", k))
}

###flex2

##2 noeuds
ALL_flex2.2 <- data.frame()
for(k in iterations){
  assign(paste0("meanTrainflex2.2_", k) , read.csv(paste0(path0, "TRAIN/WHOLE/FLEX2.2/mean_survival/",k,"_mean_flex22.csv"), sep = ";") )
  ALL_flex2.2 <- rbind(ALL_flex2.2, get(paste0("meanTrainflex2.2_",k)))
  rm(list =paste0("meanTrainflex2.2_", k))
}

##4 noeuds
ALL_flex2.4 <- data.frame()
for(k in iterations){
  assign(paste0("meanTrainflex2.4_", k) , read.csv(paste0(path0, "TRAIN/WHOLE/FLEX2.4/mean_survival/",k,"_mean_flex24.csv"), sep = ";") )
  ALL_flex2.4 <- rbind(ALL_flex2.4, get(paste0("meanTrainflex2.4_",k)))
  rm(list =paste0("meanTrainflex2.4_", k))
}
### WG
ALL_WG <- data.frame()
for(k in iterations){
  assign(paste0("meanTrainWG_", k) , read.csv(paste0(path0, "TRAIN/WHOLE/WG/mean_survival/",k,"_mean_flexWG.csv"), sep = ";") )
  ALL_WG <- rbind(ALL_WG,get(paste0("meanTrainWG_",k))) 
  rm(list =paste0("meanTrainWG_", k))
}

ALL_PP <- data.frame() 
for(k in iterations){
  assign(paste0("meanTrainPP_", k) , read.csv(paste0(path0, "TRAIN/WHOLE/PP/",k,"PP.csv"), sep = ";") )
  ALL_PP <-rbind(ALL_PP, get(paste0("meanTrainPP_", k)))
  rm(list =paste0("meanTrainPP_", k))
}
ALL_PP <- as.data.frame(matrix(ALL_PP$x, ncol = 4, byrow = TRUE))

biais_PLANN <- colMeans(ALL_PLANN[,-1] - S_theo[col( ALL_PLANN[,-1] )] )
RMSE_PLANN <- sqrt( colMeans((ALL_PLANN[,-1] - S_theo[col( ALL_PLANN[,-1] )])^2) )

biais_flex1.2 <- colMeans( ALL_flex1.2 - S_theo[col( ALL_flex1.2)] )
RMSE_flex1.2 <- sqrt( colMeans( (ALL_flex1.2 - S_theo[col( ALL_flex1.2)])^2 ) )

biais_flex1.4 <- colMeans( ALL_flex1.4 - S_theo[col( ALL_flex1.4)] )
RMSE_flex1.4 <- sqrt( colMeans( (ALL_flex1.4 - S_theo[col( ALL_flex1.4)])^2 ) )

biais_flex2.2 <- colMeans( ALL_flex2.2 - S_theo[col( ALL_flex2.2)] )
RMSE_flex2.2 <- sqrt( colMeans( (ALL_flex2.2 - S_theo[col( ALL_flex2.2)])^2 ) )

biais_flex2.4 <- colMeans( ALL_flex2.4 - S_theo[col( ALL_flex2.4)] )
RMSE_flex2.4 <- sqrt( colMeans( (ALL_flex2.4 - S_theo[col( ALL_flex2.4)])^2 ) )

biais_WG <- colMeans( ALL_WG - S_theo[col( ALL_WG)] )
RMSE_WG <- sqrt( colMeans( (ALL_WG - S_theo[col( ALL_WG)])^2 ) )

biais_PP <- colMeans( ALL_PP - S_theo[col( ALL_PP)] )
RMSE_PP <- sqrt( colMeans( (ALL_PP - S_theo[col( ALL_PP)])^2 ) )
#STRATES

###on récupère les valeurs théoriques moyennes Se(t) 
for(j in strata_names){
  assign(paste0("ALL_",j, "_theo"),  data.frame()) 
}

for(k in iterations){
  assign(paste0("meanTrainTheoS_", k) , read.csv(paste0(path0, "TRAIN/STRATA/THEO/mean_survival/",k,"mean_theo_strates.csv"),
                                                 sep = ";", row.names = 1) ) 
  for(j in strata_names){
    
    current_row <- get(paste0("meanTrainTheoS_", k))[paste0("mean.theo_", j), , drop = FALSE]
    rownames(current_row) <- paste0("mean.theo_", j, k)
    assign(paste0("ALL_", j, "_theo"), rbind(get(paste0("ALL_", j, "_theo")), current_row))
    
  }
  
  rm(list =paste0("meanTrainTheoS_", k))
}

#moyenne théorique de la valeur de survie
for(j in strata_names){
  assign(paste0("S_theo", j ),colMeans(get(paste0("ALL_",j,"_theo") ) ) )
}
############################################
#importation des données des différents modèles 
############################################

# PLANN

for(j in strata_names){
  assign(paste0("ALL_PLANN_", j),data.frame()) 
}

for(k in iterations){
  assign(paste0("meanTrainPLANNstr_", k) , read.csv(paste0(path0, "TRAIN/STRATA/PLANN/mean_survival/",k,"mean_plann_strates.csv"), sep = ";", row.names = 1) )
  
  for(j in strata_names){
    assign(paste0("ALL_PLANN_",j) ,rbind(get(paste0("ALL_PLANN_",j)),get(paste0("meanTrainPLANNstr_",k))[paste0("mean.plann_", j),-1]  ) )
  }
  rm(list =paste0("meanTrainPLANNstr_", k))
}

### flex1

##2 noeuds
for(j in strata_names){
  assign(paste0("ALL_FLEX1.2_", j),data.frame()) 
}

for(k in iterations){
  assign(paste0("meanTrainFLEX1.2str_", k) , read.csv(paste0(path0, "TRAIN/STRATA/FLEX1.2/mean_survival/",k,"mean_flex12_strates.csv"), sep = ";", row.names = 1) )
  
  for(j in strata_names){
    assign(paste0("ALL_FLEX1.2_",j) ,rbind(get(paste0("ALL_FLEX1.2_",j)), get(paste0("meanTrainFLEX1.2str_",k))[paste0("mean.flex1.2_", j),]  ) )
  }
  rm(list =paste0("meanTrainFLEX1.2str_", k))
}

##4 noeuds
for(j in strata_names){
  assign(paste0("ALL_FLEX1.4_", j),data.frame()) 
}

for(k in iterations){
  assign(paste0("meanTrainFLEX1.4str_", k) , read.csv(paste0(path0, "TRAIN/STRATA/FLEX1.4/mean_survival/",k,"mean_flex14_strates.csv"), sep = ";", row.names = 1) )
  
  for(j in strata_names){
    assign(paste0("ALL_FLEX1.4_",j) ,rbind(get(paste0("ALL_FLEX1.4_",j)), get(paste0("meanTrainFLEX1.4str_",k))[paste0("mean.flex1.4_", j),]  ) )
  }
  rm(list =paste0("meanTrainFLEX1.4str_", k))
}

### flex2
##2 noeuds
for(j in strata_names){
  assign(paste0("ALL_FLEX2.2_", j),data.frame()) 
}

for(k in iterations){
  assign(paste0("meanTrainFLEX2.2str_", k) , read.csv(paste0(path0, "TRAIN/STRATA/FLEX2.2/mean_survival/",k,"mean_flex22_strates.csv"), sep = ";", row.names = 1) )
  
  for(j in strata_names){
    assign(paste0("ALL_FLEX2.2_",j) ,rbind(get(paste0("ALL_FLEX2.2_",j)), get(paste0("meanTrainFLEX2.2str_",k))[paste0("mean.flex2.2_", j),]  ) )
  }
  rm(list =paste0("meanTrainFLEX2.2str_", k))
}

##4 noeuds
for(j in strata_names){
  assign(paste0("ALL_FLEX2.4_", j),data.frame()) 
}

for(k in iterations){
  assign(paste0("meanTrainFLEX2.4str_", k) , read.csv(paste0(path0, "TRAIN/STRATA/FLEX2.4/mean_survival/",k,"mean_flex24_strates.csv"), sep = ";", row.names = 1) )
  
  for(j in strata_names){
    assign(paste0("ALL_FLEX2.4_",j) ,rbind(get(paste0("ALL_FLEX2.4_",j)), get(paste0("meanTrainFLEX2.4str_",k))[paste0("mean.flex2.4_", j),]  ) )
  }
  rm(list =paste0("meanTrainFLEX2.4str_", k))
}
### WG

for(j in strata_names){
  assign(paste0("ALL_WG_", j),data.frame()) 
}

for(k in iterations){
  assign(paste0("meanTrainWGstr_", k) , read.csv(paste0(path0, "TRAIN/STRATA/WG/mean_survival/",k,"mean_WG_strates.csv"), sep = ";", row.names = 1) )
  
  for(j in strata_names){
    assign(paste0("ALL_WG_",j) ,rbind(get(paste0("ALL_WG_",j)), get(paste0("meanTrainWGstr_",k))[paste0("mean.WG_", j),]  ) )
  }
  rm(list =paste0("meanTrainWGstr_", k))
}

### PP

for(j in strata_names){
  assign(paste0("ALL_PP_", j),data.frame()) 
}

for(k in iterations){
  assign(paste0("meanTrainPPstr_", k) , read.csv(paste0(path0, "TRAIN/STRATA/PP/",k,"PP_strates.csv"), sep = ";", row.names = 1) )
  
  for(j in strata_names){
    assign(paste0("ALL_PP_",j) ,rbind(get(paste0("ALL_PP_",j)), get(paste0("meanTrainPPstr_",k))[paste0("poharpred_", j),]  ) )
  }
  rm(list =paste0("meanTrainPPstr_", k))
}

#### Calcul des indicateurs
##PLANN
biais_PLANN_strates <- data.frame()
RMSE_PLANN_strates <- data.frame()
for(j in strata_names){
  assign(paste0("biais_PLANN_", j), 
         colMeans( get( paste0("ALL_PLANN_", j) ) - get(paste0("S_theo", j) )[col(get( paste0("ALL_PLANN_", j)) )] ) )
  #on concatene tous les résultats dans une dataframe pour ne pas saturer l'evt 
  biais_PLANN_strates <- rbind(biais_PLANN_strates, get(paste0("biais_PLANN_", j)) )
  rownames(biais_PLANN_strates)[dim(biais_PLANN_strates)[1]] <- paste0(j)
  rm(list = paste0("biais_PLANN_", j))
  
  assign(paste0("RMSE_PLANN_",  j) ,
         sqrt( colMeans( ( get( paste0("ALL_PLANN_", j) ) - 
                             get(paste0("S_theo", j) )[col(get( paste0("ALL_PLANN_", j)) )] )^2 ) ) )
  RMSE_PLANN_strates <- rbind(RMSE_PLANN_strates, get(paste0("RMSE_PLANN_", j)) )
  rownames(RMSE_PLANN_strates)[dim(RMSE_PLANN_strates)[1]] <- paste0(j)
  rm(list = paste0("RMSE_PLANN_", j))
}

colnames(biais_PLANN_strates) <- colnames(ALL_flex1.4)
colnames(RMSE_PLANN_strates) <- colnames(ALL_flex1.4)

# FLEX 1
#2 noeuds
biais_FLEX1.2_strates <- data.frame()
RMSE_FLEX1.2_strates <- data.frame()
for(j in strata_names){
  assign(paste0("biais_FLEX1.2_", j), 
         colMeans( get( paste0("ALL_FLEX1.2_", j) ) - get(paste0("S_theo", j) )[col(get( paste0("ALL_FLEX1.2_", j)) )] ) )
  
  #on concatene tous les résultats dans une dataframe pour ne pas saturer l'evt 
  biais_FLEX1.2_strates <- rbind(biais_FLEX1.2_strates, get(paste0("biais_FLEX1.2_", j)) )
  rownames(biais_FLEX1.2_strates)[dim(biais_FLEX1.2_strates)[1]] <- paste0(j)
  rm(list = paste0("biais_FLEX1.2_", j))
  
  assign(paste0("RMSE_FLEX1.2_",  j) ,
         sqrt( colMeans( ( get( paste0("ALL_FLEX1.2_", j) ) - 
                             get(paste0("S_theo", j) )[col(get( paste0("ALL_FLEX1.2_", j)) )] )^2 ) ) )
  
  RMSE_FLEX1.2_strates <- rbind(RMSE_FLEX1.2_strates, get(paste0("RMSE_FLEX1.2_", j)) )
  rownames(RMSE_FLEX1.2_strates)[dim(RMSE_FLEX1.2_strates)[1]] <- paste0(j)
  rm(list = paste0("RMSE_FLEX1.2_", j))
  
}

colnames(biais_FLEX1.2_strates) <- colnames(ALL_flex1.2)
colnames(RMSE_FLEX1.2_strates) <- colnames(ALL_flex1.2)

#4 noeuds
biais_FLEX1.4_strates <- data.frame()
RMSE_FLEX1.4_strates <- data.frame()
for(j in strata_names){
  assign(paste0("biais_FLEX1.4_", j), 
         colMeans( get( paste0("ALL_FLEX1.4_", j) ) - get(paste0("S_theo", j) )[col(get( paste0("ALL_FLEX1.4_", j)) )] ) )
  
  #on concatene tous les résultats dans une dataframe pour ne pas saturer l'evt 
  biais_FLEX1.4_strates <- rbind(biais_FLEX1.4_strates, get(paste0("biais_FLEX1.4_", j)) )
  rownames(biais_FLEX1.4_strates)[dim(biais_FLEX1.4_strates)[1]] <- paste0(j)
  rm(list = paste0("biais_FLEX1.4_", j))
  
  assign(paste0("RMSE_FLEX1.4_",  j) ,
         sqrt( colMeans( ( get( paste0("ALL_FLEX1.4_", j) ) - 
                             get(paste0("S_theo", j) )[col(get( paste0("ALL_FLEX1.4_", j)) )] )^2 ) ) )
  
  RMSE_FLEX1.4_strates <- rbind(RMSE_FLEX1.4_strates, get(paste0("RMSE_FLEX1.4_", j)) )
  rownames(RMSE_FLEX1.4_strates)[dim(RMSE_FLEX1.4_strates)[1]] <- paste0(j)
  rm(list = paste0("RMSE_FLEX1.4_", j))
  
}

colnames(biais_FLEX1.4_strates) <- colnames(ALL_flex1.4)
colnames(RMSE_FLEX1.4_strates) <- colnames(ALL_flex1.4)


# FLEX 2
##2 noeuds
biais_FLEX2.2_strates <- data.frame()
RMSE_FLEX2.2_strates <- data.frame()
for(j in strata_names){
  assign(paste0("biais_FLEX2.2_", j), 
         colMeans( get( paste0("ALL_FLEX2.2_", j) ) - get(paste0("S_theo", j) )[col(get( paste0("ALL_FLEX2.2_", j)) )] ) )
  
  biais_FLEX2.2_strates <- rbind(biais_FLEX2.2_strates, get(paste0("biais_FLEX2.2_", j)) )
  rownames(biais_FLEX2.2_strates)[dim(biais_FLEX2.2_strates)[1]] <- paste0(j)
  rm(list = paste0("biais_FLEX2.2_", j))
  
  assign(paste0("RMSE_FLEX2.2_",  j) ,
         sqrt( colMeans( ( get( paste0("ALL_FLEX2.2_", j) ) - 
                             get(paste0("S_theo", j) )[col(get( paste0("ALL_FLEX2.2_", j)) )] )^2 ) ) )
  
  RMSE_FLEX2.2_strates <- rbind(RMSE_FLEX2.2_strates, get(paste0("RMSE_FLEX2.2_", j)) )
  rownames(RMSE_FLEX2.2_strates)[dim(RMSE_FLEX2.2_strates)[1]] <- paste0(j)
  rm(list = paste0("RMSE_FLEX2.2_", j))
}
colnames(biais_FLEX2.2_strates) <- colnames(ALL_flex2.2)
colnames(RMSE_FLEX2.2_strates) <- colnames(ALL_flex2.2)

##4 noeuds
biais_FLEX2.4_strates <- data.frame()
RMSE_FLEX2.4_strates <- data.frame()
for(j in strata_names){
  assign(paste0("biais_FLEX2.4_", j), 
         colMeans( get( paste0("ALL_FLEX2.4_", j) ) - get(paste0("S_theo", j) )[col(get( paste0("ALL_FLEX2.4_", j)) )] ) )
  
  biais_FLEX2.4_strates <- rbind(biais_FLEX2.4_strates, get(paste0("biais_FLEX2.4_", j)) )
  rownames(biais_FLEX2.4_strates)[dim(biais_FLEX2.4_strates)[1]] <- paste0(j)
  rm(list = paste0("biais_FLEX2.4_", j))
  
  assign(paste0("RMSE_FLEX2.4_",  j) ,
         sqrt( colMeans( ( get( paste0("ALL_FLEX2.4_", j) ) - 
                             get(paste0("S_theo", j) )[col(get( paste0("ALL_FLEX2.4_", j)) )] )^2 ) ) )
  
  RMSE_FLEX2.4_strates <- rbind(RMSE_FLEX2.4_strates, get(paste0("RMSE_FLEX2.4_", j)) )
  rownames(RMSE_FLEX2.4_strates)[dim(RMSE_FLEX2.4_strates)[1]] <- paste0(j)
  rm(list = paste0("RMSE_FLEX2.4_", j))
}
colnames(biais_FLEX2.4_strates) <- colnames(ALL_flex2.4)
colnames(RMSE_FLEX2.4_strates) <- colnames(ALL_flex2.4)

# WG 
biais_WG_strates <- data.frame()
RMSE_WG_strates <- data.frame()
for(j in strata_names){
  assign(paste0("biais_WG_", j), 
         colMeans( get( paste0("ALL_WG_", j) ) - get(paste0("S_theo", j) )[col(get( paste0("ALL_WG_", j)) )] ) )
  
  biais_WG_strates <- rbind(biais_WG_strates, get(paste0("biais_WG_", j)) )
  rownames(biais_WG_strates)[dim(biais_WG_strates)[1]] <- paste0(j)
  rm(list = paste0("biais_WG_", j))
  
  assign(paste0("RMSE_WG_",  j) ,
         sqrt( colMeans( ( get( paste0("ALL_WG_", j) ) - 
                             get(paste0("S_theo", j) )[col(get( paste0("ALL_WG_", j)) )] )^2 ) ) )
  
  RMSE_WG_strates <- rbind(RMSE_WG_strates, get(paste0("RMSE_WG_", j)) )
  rownames(RMSE_WG_strates)[dim(RMSE_WG_strates)[1]] <- paste0(j)
  rm(list = paste0("RMSE_WG_", j))
}

colnames(biais_WG_strates) <- colnames(ALL_WG)
colnames(RMSE_WG_strates) <- colnames(ALL_WG)

# PP
biais_PP_strates <- data.frame()
RMSE_PP_strates <- data.frame()
for(j in strata_names){
  assign(paste0("biais_PP_", j), 
         colMeans( get( paste0("ALL_PP_", j) ) - get(paste0("S_theo", j) )[col(get( paste0("ALL_PP_", j)) )] ) )
  
  biais_PP_strates <- rbind(biais_PP_strates, get(paste0("biais_PP_", j)) )
  rownames(biais_PP_strates)[dim(biais_PP_strates)[1]] <- paste0(j)
  rm(list = paste0("biais_PP_", j))
  
  assign(paste0("RMSE_PP_",  j) ,
         sqrt( colMeans( ( get( paste0("ALL_PP_", j) ) - 
                             get(paste0("S_theo", j) )[col(get( paste0("ALL_PP_", j)) )] )^2 ) ) )
  
  RMSE_PP_strates <- rbind(RMSE_PP_strates, get(paste0("RMSE_PP_", j)) )
  rownames(RMSE_PP_strates)[dim(RMSE_PP_strates)[1]] <- paste0(j)
  rm(list = paste0("RMSE_PP_", j))
}

colnames(biais_PP_strates) <- colnames(ALL_PP)
colnames(RMSE_PP_strates) <- colnames(ALL_PP)

###concaten des S_theo, là aussi pour libérer de la place dans l'evt global
S_theo_strates <- data.frame()
for(j in strata_names){
  S_theo_strates <- rbind(S_theo_strates, get(paste0("S_theo",j)))
  rownames(S_theo_strates)[dim(S_theo_strates)[1]] <- paste0(j)
  rm(list = paste0("S_theo",j))
}
colnames(S_theo_strates) <- colnames(ALL_WG)

############################################ VALID #################################################

# WHOLE
###on récupère les valeurs théoriques moyennes Se(t) 
ALL_theoval <- data.frame()

for(k in iterations){
  assign(paste0("meanValidTheo_", k) , read.csv(paste0(path0, "VALID/WHOLE/THEO/mean_survival/",k,"_mean_theoval.csv"),
                                                sep = ";") )
  
  ALL_theoval <-rbind(ALL_theoval, get(paste0("meanValidTheo_", k)))
  rm(list =paste0("meanValidTheo_", k))
}

#moyenne théorique de la valur de survie
S_theoval <- colMeans(ALL_theoval)

############################################
#importation des données des différents modèles 
#########
##plann
ALL_PLANNval <- data.frame() 
for(k in iterations){
  assign(paste0("meanValidPLANN_", k) , read.csv(paste0(path0, "VALID/WHOLE/PLANN/mean_survival/",k,"_mean_PLANNval.csv"), sep = ";") )
  ALL_PLANNval <-rbind(ALL_PLANNval, get(paste0("meanValidPLANN_", k)))
  rm(list =paste0("meanValidPLANN_", k))
}

###flex1
##2 noeuds
ALL_flex1.2val <- data.frame()
for(k in iterations){
  assign(paste0("meanValidflex1.2_", k) , read.csv(paste0(path0, "VALID/WHOLE/FLEX1.2/mean_survival/",k,"_mean_flexval12.csv"), sep = ";") )
  ALL_flex1.2val <- rbind(ALL_flex1.2val, get(paste0("meanValidflex1.2_",k)))
  rm(list =paste0("meanValidflex1.2_", k))
}
##4 noeuds
ALL_flex1.4val <- data.frame()
for(k in iterations){
  assign(paste0("meanValidflex1.4_", k) , read.csv(paste0(path0, "VALID/WHOLE/FLEX1.4/mean_survival/",k,"_mean_flexval14.csv"), sep = ";") )
  ALL_flex1.4val <- rbind(ALL_flex1.4val, get(paste0("meanValidflex1.4_",k)))
  rm(list =paste0("meanValidflex1.4_", k))
}
### flex2

##2 noeuds
ALL_flex2.2val <- data.frame()
for(k in iterations){
  assign(paste0("meanValidflex2.2_", k) , read.csv(paste0(path0, "VALID/WHOLE/FLEX2.2/mean_survival/",k,"_mean_flexval22.csv"), sep = ";") )
  ALL_flex2.2val <- rbind(ALL_flex2.2val, get(paste0("meanValidflex2.2_",k)))
  rm(list =paste0("meanValidflex2.2_", k))
}

##4 noeuds
ALL_flex2.4val <- data.frame()
for(k in iterations){
  assign(paste0("meanValidflex2.4_", k) , read.csv(paste0(path0, "VALID/WHOLE/FLEX2.4/mean_survival/",k,"_mean_flexval24.csv"), sep = ";") )
  ALL_flex2.4val <- rbind(ALL_flex2.4val, get(paste0("meanValidflex2.4_",k)))
  rm(list =paste0("meanValidflex2.4_", k))
}
### WG
ALL_WGval <- data.frame()
for(k in iterations){
  assign(paste0("meanValidWG_", k) , read.csv(paste0(path0, "VALID/WHOLE/WG/mean_survival/",k,"_mean_flexWGval.csv"), sep = ";") )
  ALL_WGval <- rbind(ALL_WGval, get(paste0("meanValidWG_",k))) 
  rm(list =paste0("meanValidWG_", k))
}

##PP
ALL_PPval <- data.frame()
for(k in iterations){
  assign(paste0("meanValidPP_", k) , read.csv(paste0(path0, "VALID/WHOLE/PP/",k,"PPval.csv"), sep = ";") )
  ALL_PPval <- rbind(ALL_PPval, get(paste0("meanValidPP_",k))) 
  rm(list =paste0("meanValidPP_", k))
}
ALL_PPval <- as.data.frame(matrix(ALL_PPval$x, ncol = 4, byrow = TRUE))

biais_PLANNval <- colMeans(ALL_PLANNval[,-1] - S_theoval[col( ALL_PLANNval[,-1] )] )
RMSE_PLANNval <- sqrt( colMeans((ALL_PLANNval[,-1] - S_theoval[col( ALL_PLANNval[,-1] )])^2) )

biais_flex1.2val <- colMeans( ALL_flex1.2val - S_theoval[col( ALL_flex1.2val)] )
RMSE_flex1.2val <- sqrt( colMeans( (ALL_flex1.2val - S_theoval[col( ALL_flex1.2val)])^2 ) )

biais_flex1.4val <- colMeans( ALL_flex1.4val - S_theoval[col( ALL_flex1.4val)] )
RMSE_flex1.4val <- sqrt( colMeans( (ALL_flex1.4val - S_theoval[col( ALL_flex1.4val)])^2 ) )

biais_flex2.2val <- colMeans( ALL_flex2.2val - S_theoval[col( ALL_flex2.2val)] )
RMSE_flex2.2val <- sqrt( colMeans( (ALL_flex2.2val - S_theoval[col( ALL_flex2.2val)])^2 ) )

biais_flex2.4val <- colMeans( ALL_flex2.4val - S_theoval[col( ALL_flex2.4val)] )
RMSE_flex2.4val <- sqrt( colMeans( (ALL_flex2.4val - S_theoval[col( ALL_flex2.4val)])^2 ) )

biais_WGval <- colMeans( ALL_WGval - S_theoval[col( ALL_WGval)] )
RMSE_WGval <- sqrt( colMeans( (ALL_WGval - S_theoval[col( ALL_WGval)])^2 ) )

biais_PPval <- colMeans( ALL_PP - S_theoval[col( ALL_PPval)] )
RMSE_PPval <- sqrt( colMeans( (ALL_PP - S_theoval[col( ALL_PPval)])^2 ) )
#STRATES

###on récupère les valeurs théoriques moyennes Se(t) 
for(j in strata_names){
  assign(paste0("ALL_",j, "_theoval"),  data.frame()) 
}

for(k in iterations){
  assign(paste0("meanValidTheoS_", k) , read.csv(paste0(path0, "VALID/STRATA/THEO/mean_survival/",k,"mean_theoval_strates.csv"),
                                                 sep = ";", row.names = 1)[,-1] ) ### [,-1] car quand j'ai enregistré les résultats, j'ai dupliué la premiere colonne donc il faut l'enlever
  for(j in strata_names){
    
    current_row <- get(paste0("meanValidTheoS_", k))[paste0("mean.theo_", j), , drop = FALSE]
    rownames(current_row) <- paste0("mean.theoval_", j, k)
    assign(paste0("ALL_", j, "_theoval"), rbind(get(paste0("ALL_", j, "_theoval")), current_row))
    
  }
  
  rm(list =paste0("meanValidTheoS_", k))
}
rm(current_row)
#moyenne théorique de la valeur de survie
for(j in strata_names){
  assign(paste0("S_theoval", j ),colMeans(get(paste0("ALL_",j,"_theoval") ) ) )
}
############################################
#importation des données des différents modèles 
############################################

# PLANN

for(j in strata_names){
  assign(paste0("ALL_PLANNval_", j),data.frame()) 
}

for(k in iterations){
  assign(paste0("meanValidPLANNstr_", k) , read.csv(paste0(path0, "VALID/STRATA/PLANN/mean_survival/",k,"mean_plannval_strates.csv"), sep = ";", row.names = 1) )
  
  for(j in strata_names){
    assign(paste0("ALL_PLANNval_",j) ,rbind(get(paste0("ALL_PLANNval_",j)),get(paste0("meanValidPLANNstr_",k))[paste0("mean.plannval_", j),-1]  ) )
  }
  rm(list =paste0("meanValidPLANNstr_", k))
}

### flex1
##2 noeuds
for(j in strata_names){
  assign(paste0("ALL_FLEX1.2val_", j),data.frame()) 
}

for(k in iterations){
  assign(paste0("meanValidFLEX1.2str_", k) , read.csv(paste0(path0, "VALID/STRATA/FLEX1.2/mean_survival/",k,"mean_flexval12_strates.csv"), sep = ";", row.names = 1) )
  
  for(j in strata_names){
    assign(paste0("ALL_FLEX1.2val_",j) ,rbind(get(paste0("ALL_FLEX1.2val_",j)), get(paste0("meanValidFLEX1.2str_",k))[paste0("mean.flexval1.2_", j),]  ) )
  }
  rm(list =paste0("meanValidFLEX1.2str_", k))
}

##4 noeuds
for(j in strata_names){
  assign(paste0("ALL_FLEX1.4val_", j),data.frame()) 
}

for(k in iterations){
  assign(paste0("meanValidFLEX1.4str_", k) , read.csv(paste0(path0, "VALID/STRATA/FLEX1.4/mean_survival/",k,"mean_flexval14_strates.csv"), sep = ";", row.names = 1) )
  
  for(j in strata_names){
    assign(paste0("ALL_FLEX1.4val_",j) ,rbind(get(paste0("ALL_FLEX1.4val_",j)), get(paste0("meanValidFLEX1.4str_",k))[paste0("mean.flexval1.4_", j),]  ) )
  }
  rm(list =paste0("meanValidFLEX1.4str_", k))
}
### flex2

##2 noeuds
for(j in strata_names){
  assign(paste0("ALL_FLEX2.2val_", j),data.frame()) 
}

for(k in iterations){
  assign(paste0("meanValidFLEX2.2str_", k) , read.csv(paste0(path0, "VALID/STRATA/FLEX2.2/mean_survival/",k,"mean_flexval22_strates.csv"), sep = ";", row.names = 1) )
  
  for(j in strata_names){
    assign(paste0("ALL_FLEX2.2val_",j) ,rbind(get(paste0("ALL_FLEX2.2val_",j)), get(paste0("meanValidFLEX2.2str_",k))[paste0("mean.flexval2.2_", j),]  ) )
  }
  rm(list =paste0("meanValidFLEX2.2str_", k))
}

##4 noeuds
for(j in strata_names){
  assign(paste0("ALL_FLEX2.4val_", j),data.frame()) 
}

for(k in iterations){
  assign(paste0("meanValidFLEX2.4str_", k) , read.csv(paste0(path0, "VALID/STRATA/FLEX2.4/mean_survival/",k,"mean_flexval24_strates.csv"), sep = ";", row.names = 1) )
  
  for(j in strata_names){
    assign(paste0("ALL_FLEX2.4val_",j) ,rbind(get(paste0("ALL_FLEX2.4val_",j)), get(paste0("meanValidFLEX2.4str_",k))[paste0("mean.flexval2.4_", j),]  ) )
  }
  rm(list =paste0("meanValidFLEX2.4str_", k))
}
### WG

for(j in strata_names){
  assign(paste0("ALL_WGval_", j),data.frame()) 
}

for(k in iterations){
  assign(paste0("meanValidWGstr_", k) , read.csv(paste0(path0, "VALID/STRATA/WG/mean_survival/",k,"mean_WGval_strates.csv"), sep = ";", row.names = 1) )
  
  for(j in strata_names){
    assign(paste0("ALL_WGval_",j) ,rbind(get(paste0("ALL_WGval_",j)), get(paste0("meanValidWGstr_",k))[paste0("mean.WGval_", j),]  ) )
  }
  rm(list =paste0("meanValidWGstr_", k))
}

### PP 

for(j in strata_names){
  assign(paste0("ALL_PPval_", j),data.frame()) 
}

for(k in iterations){
  assign(paste0("meanValidPPstr_", k) , read.csv(paste0(path0, "VALID/STRATA/PP/",k,"PPval_strates.csv"), sep = ";", row.names = 1) )
  
  for(j in strata_names){
    assign(paste0("ALL_PPval_",j) ,rbind(get(paste0("ALL_PPval_",j)), get(paste0("meanValidPPstr_",k))[paste0("poharpredval_", j),]  ) )
  }
  rm(list =paste0("meanValidPPstr_", k))
}
#### Calcul des indicateurs
##PLANN
biais_PLANNval_strates <- data.frame()
RMSE_PLANNval_strates <- data.frame()
for(j in strata_names){
  assign(paste0("biais_PLANNval_", j), 
         colMeans( get( paste0("ALL_PLANNval_", j) ) - get(paste0("S_theoval", j) )[col(get( paste0("ALL_PLANNval_", j)) )] ) )
  
  biais_PLANNval_strates <- rbind(biais_PLANNval_strates, get(paste0("biais_PLANNval_", j)) )
  rownames(biais_PLANNval_strates)[dim(biais_PLANNval_strates)[1]] <- paste0(j)
  rm(list = paste0("biais_PLANNval_", j))
  
  assign(paste0("RMSE_PLANNval_",  j) ,
         sqrt( colMeans( ( get( paste0("ALL_PLANNval_", j) ) - 
                             get(paste0("S_theoval", j) )[col(get( paste0("ALL_PLANNval_", j)) )] )^2 ) ) )
  
  RMSE_PLANNval_strates <- rbind(RMSE_PLANNval_strates, get(paste0("RMSE_PLANNval_", j)) )
  rownames(RMSE_PLANNval_strates)[dim(RMSE_PLANNval_strates)[1]] <- paste0(j)
  rm(list = paste0("RMSE_PLANNval_", j))
  
}

colnames(biais_PLANNval_strates) <- colnames(ALL_flex1.2)
colnames(RMSE_PLANNval_strates) <- colnames(ALL_flex1.2)

# FLEX 1

##2 noeuds
biais_FLEX1.2val_strates <- data.frame()
RMSE_FLEX1.2val_strates <- data.frame()
for(j in strata_names){
  assign(paste0("biais_FLEX1.2val_", j), 
         colMeans( get( paste0("ALL_FLEX1.2val_", j) ) - get(paste0("S_theoval", j) )[col(get( paste0("ALL_FLEX1.2val_", j)) )] ) )
  
  biais_FLEX1.2val_strates <- rbind(biais_FLEX1.2val_strates, get(paste0("biais_FLEX1.2val_", j)) )
  rownames(biais_FLEX1.2val_strates)[dim(biais_FLEX1.2val_strates)[1]] <- paste0(j)
  rm(list = paste0("biais_FLEX1.2val_", j))
  
  assign(paste0("RMSE_FLEX1.2val_",  j) ,
         sqrt( colMeans( ( get( paste0("ALL_FLEX1.2val_", j) ) - 
                             get(paste0("S_theoval", j) )[col(get( paste0("ALL_FLEX1.2val_", j)) )] )^2 ) ) )
  
  RMSE_FLEX1.2val_strates <- rbind(RMSE_FLEX1.2val_strates, get(paste0("RMSE_FLEX1.2val_", j)) )
  rownames(RMSE_FLEX1.2val_strates)[dim(RMSE_FLEX1.2val_strates)[1]] <- paste0(j)
  rm(list = paste0("RMSE_FLEX1.2val_", j))
}

colnames(biais_FLEX1.2val_strates) <- colnames(ALL_flex1.2)
colnames(RMSE_FLEX1.2val_strates) <- colnames(ALL_flex1.2)

##4 noeuds
biais_FLEX1.4val_strates <- data.frame()
RMSE_FLEX1.4val_strates <- data.frame()
for(j in strata_names){
  assign(paste0("biais_FLEX1.4val_", j), 
         colMeans( get( paste0("ALL_FLEX1.4val_", j) ) - get(paste0("S_theoval", j) )[col(get( paste0("ALL_FLEX1.4val_", j)) )] ) )
  
  biais_FLEX1.4val_strates <- rbind(biais_FLEX1.4val_strates, get(paste0("biais_FLEX1.4val_", j)) )
  rownames(biais_FLEX1.4val_strates)[dim(biais_FLEX1.4val_strates)[1]] <- paste0(j)
  rm(list = paste0("biais_FLEX1.4val_", j))
  
  assign(paste0("RMSE_FLEX1.4val_",  j) ,
         sqrt( colMeans( ( get( paste0("ALL_FLEX1.4val_", j) ) - 
                             get(paste0("S_theoval", j) )[col(get( paste0("ALL_FLEX1.4val_", j)) )] )^2 ) ) )
  
  RMSE_FLEX1.4val_strates <- rbind(RMSE_FLEX1.4val_strates, get(paste0("RMSE_FLEX1.4val_", j)) )
  rownames(RMSE_FLEX1.4val_strates)[dim(RMSE_FLEX1.4val_strates)[1]] <- paste0(j)
  rm(list = paste0("RMSE_FLEX1.4val_", j))
}

colnames(biais_FLEX1.4val_strates) <- colnames(ALL_flex1.4)
colnames(RMSE_FLEX1.4val_strates) <- colnames(ALL_flex1.4)

# FLEX 2

##2 noeuds
biais_FLEX2.2val_strates <- data.frame()
RMSE_FLEX2.2val_strates <- data.frame()

for(j in strata_names){
  assign(paste0("biais_FLEX2.2val_", j), 
         colMeans( get( paste0("ALL_FLEX2.2val_", j) ) - get(paste0("S_theoval", j) )[col(get( paste0("ALL_FLEX2.2val_", j)) )] ) )
  
  biais_FLEX2.2val_strates <- rbind(biais_FLEX2.2val_strates, get(paste0("biais_FLEX2.2val_", j)) )
  rownames(biais_FLEX2.2val_strates)[dim(biais_FLEX2.2val_strates)[1]] <- paste0(j)
  rm(list = paste0("biais_FLEX2.2val_", j))
  
  assign(paste0("RMSE_FLEX2.2val_",  j) ,
         sqrt( colMeans( ( get( paste0("ALL_FLEX2.2val_", j) ) - 
                             get(paste0("S_theoval", j) )[col(get( paste0("ALL_FLEX2.2val_", j)) )] )^2 ) ) )
  
  RMSE_FLEX2.2val_strates <- rbind(RMSE_FLEX2.2val_strates, get(paste0("RMSE_FLEX2.2val_", j)) )
  rownames(RMSE_FLEX2.2val_strates)[dim(RMSE_FLEX2.2val_strates)[1]] <- paste0(j)
  rm(list = paste0("RMSE_FLEX2.2val_", j))
  
}

colnames(biais_FLEX2.2val_strates) <- colnames(ALL_flex1.2)
colnames(RMSE_FLEX2.2val_strates) <- colnames(ALL_flex1.2)

##4 noeuds
biais_FLEX2.4val_strates <- data.frame()
RMSE_FLEX2.4val_strates <- data.frame()

for(j in strata_names){
  assign(paste0("biais_FLEX2.4val_", j), 
         colMeans( get( paste0("ALL_FLEX2.4val_", j) ) - get(paste0("S_theoval", j) )[col(get( paste0("ALL_FLEX2.4val_", j)) )] ) )
  
  biais_FLEX2.4val_strates <- rbind(biais_FLEX2.4val_strates, get(paste0("biais_FLEX2.4val_", j)) )
  rownames(biais_FLEX2.4val_strates)[dim(biais_FLEX2.4val_strates)[1]] <- paste0(j)
  rm(list = paste0("biais_FLEX2.4val_", j))
  
  assign(paste0("RMSE_FLEX2.4val_",  j) ,
         sqrt( colMeans( ( get( paste0("ALL_FLEX2.4val_", j) ) - 
                             get(paste0("S_theoval", j) )[col(get( paste0("ALL_FLEX2.4val_", j)) )] )^2 ) ) )
  
  RMSE_FLEX2.4val_strates <- rbind(RMSE_FLEX2.4val_strates, get(paste0("RMSE_FLEX2.4val_", j)) )
  rownames(RMSE_FLEX2.4val_strates)[dim(RMSE_FLEX2.4val_strates)[1]] <- paste0(j)
  rm(list = paste0("RMSE_FLEX2.4val_", j))
  
}

colnames(biais_FLEX2.4val_strates) <- colnames(ALL_flex1.2)
colnames(RMSE_FLEX2.4val_strates) <- colnames(ALL_flex1.2)
# WG 

biais_WGval_strates <- data.frame()
RMSE_WGval_strates <- data.frame()

for(j in strata_names){
  assign(paste0("biais_WGval_", j), 
         colMeans( get( paste0("ALL_WGval_", j) ) - get(paste0("S_theoval", j) )[col(get( paste0("ALL_WGval_", j)) )] ) )
  
  biais_WGval_strates <- rbind(biais_WGval_strates, get(paste0("biais_WGval_", j)) )
  rownames(biais_WGval_strates)[dim(biais_WGval_strates)[1]] <- paste0(j)
  rm(list = paste0("biais_WGval_", j))
  
  assign(paste0("RMSE_WGval_",  j) ,
         sqrt( colMeans( ( get( paste0("ALL_WGval_", j) ) - 
                             get(paste0("S_theoval", j) )[col(get( paste0("ALL_WGval_", j)) )] )^2 ) ) )
  
  RMSE_WGval_strates <- rbind(RMSE_WGval_strates, get(paste0("RMSE_WGval_", j)) )
  rownames(RMSE_WGval_strates)[dim(RMSE_WGval_strates)[1]] <- paste0(j)
  rm(list = paste0("RMSE_WGval_", j))
  
}

colnames(biais_WGval_strates) <- colnames(ALL_flex1.2)
colnames(RMSE_WGval_strates) <- colnames(ALL_flex1.2)

# PP ############################################## A CHANGER 
biais_PPval_strates <- data.frame()
RMSE_PPval_strates <- data.frame()
for(j in strata_names){
  assign(paste0("biais_PPval_", j), 
         colMeans( get( paste0("ALL_PP_", j) ) - get(paste0("S_theoval", j) )[col(get( paste0("ALL_PP_", j)) )] ) )
  
  biais_PPval_strates <- rbind(biais_PPval_strates, get(paste0("biais_PPval_", j)) )
  rownames(biais_PPval_strates)[dim(biais_PPval_strates)[1]] <- paste0(j)
  rm(list = paste0("biais_PPval_", j))
  
  assign(paste0("RMSE_PPval_",  j) ,
         sqrt( colMeans( ( get( paste0("ALL_PP_", j) ) - 
                             get(paste0("S_theoval", j) )[col(get( paste0("ALL_PP_", j)) )] )^2 ) ) )
  
  RMSE_PPval_strates <- rbind(RMSE_PPval_strates, get(paste0("RMSE_PPval_", j)) )
  rownames(RMSE_PPval_strates)[dim(RMSE_PPval_strates)[1]] <- paste0(j)
  rm(list = paste0("RMSE_PPval_", j))
}

colnames(biais_PPval_strates) <- colnames(ALL_PPval)
colnames(RMSE_PPval_strates) <- colnames(ALL_PPval)

###concaten des S_theo, là aussi pour libérer de la place dans l'evt global
S_theoval_strates <- data.frame()
for(j in strata_names){
  S_theoval_strates <- rbind(S_theoval_strates, get(paste0("S_theoval",j)))
  rownames(S_theoval_strates)[dim(S_theoval_strates)[1]] <- paste0(j)
  rm(list = paste0("S_theoval",j))
}
colnames(S_theoval_strates) <- colnames(ALL_WG)


## mise en forme de l'evt global 
##############################################
##concatenation STHEO/ STHEOSTRATES
##TRAIN
S_theo <- rbind(S_theo, S_theo_strates)
rownames(S_theo)[1] <- "ALL" 
rm(S_theo_strates)
##VALID

S_theoval <- rbind(S_theoval, S_theoval_strates)
rownames(S_theoval)[1] <- "ALL" 
rm(S_theoval_strates)
## concatenation des résultats sur la base totale avec ceux fait en strates 
##TRAIN
##PLANN
biais_PLANN <- rbind(biais_PLANN, biais_PLANN_strates)
RMSE_PLANN <- rbind(RMSE_PLANN, RMSE_PLANN_strates)
rownames(biais_PLANN)[1] <- "ALL"
rownames(RMSE_PLANN)[1] <- "ALL"
rm(biais_PLANN_strates)
rm(RMSE_PLANN_strates)
##FLEX1
##2noeuds
biais_FLEX1.2 <- rbind(biais_flex1.2, biais_FLEX1.2_strates)
RMSE_FLEX1.2 <- rbind(RMSE_flex1.2, RMSE_FLEX1.2_strates)
rownames(biais_FLEX1.2)[1] <- "ALL"
rownames(RMSE_FLEX1.2)[1] <- "ALL"
rm(biais_flex1.2, biais_FLEX1.2_strates)
rm(RMSE_flex1.2, RMSE_FLEX1.2_strates)

##4noeuds
biais_FLEX1.4 <- rbind(biais_flex1.4, biais_FLEX1.4_strates)
RMSE_FLEX1.4 <- rbind(RMSE_flex1.4, RMSE_FLEX1.4_strates)
rownames(biais_FLEX1.4)[1] <- "ALL"
rownames(RMSE_FLEX1.4)[1] <- "ALL"
rm(biais_flex1.4, biais_FLEX1.4_strates)
rm(RMSE_flex1.4, RMSE_FLEX1.4_strates)

##FLEX2
##2 noeuds
biais_FLEX2.2 <- rbind(biais_flex2.2, biais_FLEX2.2_strates)
RMSE_FLEX2.2 <- rbind(RMSE_flex2.2, RMSE_FLEX2.2_strates)
rownames(biais_FLEX2.2)[1] <- "ALL"
rownames(RMSE_FLEX2.2)[1] <- "ALL"
rm(biais_flex2.2, biais_FLEX2.2_strates)
rm(RMSE_flex2.2, RMSE_FLEX2.2_strates)
##4 noeuds
biais_FLEX2.4 <- rbind(biais_flex2.4, biais_FLEX2.4_strates)
RMSE_FLEX2.4 <- rbind(RMSE_flex2.4, RMSE_FLEX2.4_strates)
rownames(biais_FLEX2.4)[1] <- "ALL"
rownames(RMSE_FLEX2.4)[1] <- "ALL"
rm(biais_flex2.4, biais_FLEX2.4_strates)
rm(RMSE_flex2.4, RMSE_FLEX2.4_strates)
##WG

biais_WG <- rbind(biais_WG, biais_WG_strates)
RMSE_WG <- rbind(RMSE_WG, RMSE_WG_strates)
rownames(biais_WG)[1] <- "ALL"
rownames(RMSE_WG)[1] <- "ALL"
rm(biais_WG_strates)
rm(RMSE_WG_strates)

##PP

biais_PP <- rbind(biais_PP, biais_PP_strates)
RMSE_PP <- rbind(RMSE_PP, RMSE_PP_strates)
rownames(biais_PP)[1] <- "ALL"
rownames(RMSE_PP)[1] <- "ALL"
rm(biais_PP_strates)
rm(RMSE_PP_strates)

######VALID

biais_PLANNval <- rbind(biais_PLANNval, biais_PLANNval_strates)
RMSE_PLANNval <- rbind(RMSE_PLANNval, RMSE_PLANNval_strates)
rownames(biais_PLANNval)[1] <- "ALL"
rownames(RMSE_PLANNval)[1] <- "ALL"
rm(biais_PLANNval_strates)
rm(RMSE_PLANNval_strates)
##FLEX1
##2 noeuds
biais_FLEX1.2val <- rbind(biais_flex1.2val, biais_FLEX1.2val_strates)
RMSE_FLEX1.2val <- rbind(RMSE_flex1.2val, RMSE_FLEX1.2val_strates)
rownames(biais_FLEX1.2val)[1] <- "ALL"
rownames(RMSE_FLEX1.2val)[1] <- "ALL"
rm(biais_flex1.2val, biais_FLEX1.2val_strates)
rm(RMSE_flex1.2val, RMSE_FLEX1.2val_strates)
##4 noeuds
biais_FLEX1.4val <- rbind(biais_flex1.4val, biais_FLEX1.4val_strates)
RMSE_FLEX1.4val <- rbind(RMSE_flex1.4val, RMSE_FLEX1.4val_strates)
rownames(biais_FLEX1.4val)[1] <- "ALL"
rownames(RMSE_FLEX1.4val)[1] <- "ALL"
rm(biais_flex1.4val, biais_FLEX1.4val_strates)
rm(RMSE_flex1.4val, RMSE_FLEX1.4val_strates)
##FLEX2
##2 noeuds
biais_FLEX2.2val <- rbind(biais_flex2.2val, biais_FLEX2.2val_strates)
RMSE_FLEX2.2val <- rbind(RMSE_flex2.2val, RMSE_FLEX2.2val_strates)
rownames(biais_FLEX2.2val)[1] <- "ALL"
rownames(RMSE_FLEX2.2val)[1] <- "ALL"
rm(biais_flex2.2val, biais_FLEX2.2val_strates)
rm(RMSE_flex2.2val, RMSE_FLEX2.2val_strates)
##4 noeuds
biais_FLEX2.4val <- rbind(biais_flex2.4val, biais_FLEX2.4val_strates)
RMSE_FLEX2.4val <- rbind(RMSE_flex2.4val, RMSE_FLEX2.4val_strates)
rownames(biais_FLEX2.4val)[1] <- "ALL"
rownames(RMSE_FLEX2.4val)[1] <- "ALL"
rm(biais_flex2.4val, biais_FLEX2.4val_strates)
rm(RMSE_flex2.4val, RMSE_FLEX2.4val_strates)
##WG

biais_WGval <- rbind(biais_WGval, biais_WGval_strates)
RMSE_WGval <- rbind(RMSE_WGval, RMSE_WGval_strates)
rownames(biais_WGval)[1] <- "ALL"
rownames(RMSE_WGval)[1] <- "ALL"
rm(biais_WGval_strates)
rm(RMSE_WGval_strates)

##PP

biais_PPval <- rbind(biais_PPval, biais_PPval_strates)
RMSE_PPval <- rbind(RMSE_PPval, RMSE_PPval_strates)
rownames(biais_PPval)[1] <- "ALL"
rownames(RMSE_PPval)[1] <- "ALL"
rm(biais_PPval_strates)
rm(RMSE_PPval_strates)

### les dataframes des moyennes sont enregistrées sous forme de liste par rapport au modèle utilisé
########
##TRAIN
THEO_means_names <- c()
for(j in strata_names){
  THEO_means_names <- c(THEO_means_names, paste0("ALL_",j,"_theo"))
}
THEO_means <- c(mget("ALL_theo"), mget(THEO_means_names))
rm(ALL_theo)
rm(THEO_means_names)
for(j in strata_names){
  rm(list = paste0("ALL_",j,"_theo") )
}

##PLANN
PLANN_means_names <- c()
for(j in strata_names){
  PLANN_means_names <- c(PLANN_means_names, paste0("ALL_PLANN_",j))
}
PLANN_means <- c(mget("ALL_PLANN"), mget(PLANN_means_names))
rm(ALL_PLANN)
rm(PLANN_means_names)
for(j in strata_names){
  rm(list = paste0("ALL_PLANN_",j) )
}

##FLEX1
##2 noeuds
FLEX1.2_means_names <- c()
for(j in strata_names){
  FLEX1.2_means_names <- c(FLEX1.2_means_names, paste0("ALL_FLEX1.2_",j))
}
FLEX1.2_means <- c(mget("ALL_flex1.2"), mget(FLEX1.2_means_names))
rm(ALL_flex1.2)
rm(FLEX1.2_means_names)
for(j in strata_names){
  rm(list = paste0("ALL_FLEX1.2_",j) )
}

##4 noeuds
FLEX1.4_means_names <- c()
for(j in strata_names){
  FLEX1.4_means_names <- c(FLEX1.4_means_names, paste0("ALL_FLEX1.4_",j))
}
FLEX1.4_means <- c(mget("ALL_flex1.4"), mget(FLEX1.4_means_names))
rm(ALL_flex1.4)
rm(FLEX1.4_means_names)
for(j in strata_names){
  rm(list = paste0("ALL_FLEX1.4_",j) )
}

##FLEX2
##2 noeuds
FLEX2.2_means_names <- c()
for(j in strata_names){
  FLEX2.2_means_names <- c(FLEX2.2_means_names, paste0("ALL_FLEX2.2_",j))
}
FLEX2.2_means <- c(mget("ALL_flex2.2"), mget(FLEX2.2_means_names))
rm(ALL_flex2.2)
rm(FLEX2.2_means_names)
for(j in strata_names){
  rm(list = paste0("ALL_FLEX2.2_",j) )
}
##4 noeuds
FLEX2.4_means_names <- c()
for(j in strata_names){
  FLEX2.4_means_names <- c(FLEX2.4_means_names, paste0("ALL_FLEX2.4_",j))
}
FLEX2.4_means <- c(mget("ALL_flex2.4"), mget(FLEX2.4_means_names))
rm(ALL_flex2.4)
rm(FLEX2.4_means_names)
for(j in strata_names){
  rm(list = paste0("ALL_FLEX2.4_",j) )
}
##WG
WG_means_names <- c()
for(j in strata_names){
  WG_means_names <- c(WG_means_names, paste0("ALL_WG_",j))
}
WG_means <- c(mget("ALL_WG"), mget(WG_means_names))
rm(ALL_WG)
rm(WG_means_names)
for(j in strata_names){
  rm(list = paste0("ALL_WG_",j) )
}

PP_means_names <- c()
for(j in strata_names){
  PP_means_names <- c(PP_means_names, paste0("ALL_PP_",j))
}
PP_means <- c(mget("ALL_PP"), mget(PP_means_names))
rm(ALL_PP)
rm(PP_means_names)
for(j in strata_names){
  rm(list = paste0("ALL_PP_",j) )
}
############VALID
THEOval_means_names <- c()
for(j in strata_names){
  THEOval_means_names <- c(THEOval_means_names, paste0("ALL_",j,"_theoval"))
}
THEOval_means <- c(mget("ALL_theoval"), mget(THEOval_means_names))
rm(ALL_theoval)
rm(THEOval_means_names)
for(j in strata_names){
  rm(list = paste0("ALL_",j,"_theoval") )
}

##PLANN
PLANNval_means_names <- c()
for(j in strata_names){
  PLANNval_means_names <- c(PLANNval_means_names, paste0("ALL_PLANNval_",j))
}
PLANNval_means <- c(mget("ALL_PLANNval"), mget(PLANNval_means_names))
rm(ALL_PLANNval)
rm(PLANNval_means_names)
for(j in strata_names){
  rm(list = paste0("ALL_PLANNval_",j) )
}

##FLEX1
##2 noeuds
FLEX1.2val_means_names <- c()
for(j in strata_names){
  FLEX1.2val_means_names <- c(FLEX1.2val_means_names, paste0("ALL_FLEX1.2val_",j))
}
FLEX1.2val_means <- c(mget("ALL_flex1.2val"), mget(FLEX1.2val_means_names))
rm(ALL_flex1.2val)
rm(FLEX1.2val_means_names)
for(j in strata_names){
  rm(list = paste0("ALL_FLEX1.2val_",j) )
}

##4 noeuds
FLEX1.4val_means_names <- c()
for(j in strata_names){
  FLEX1.4val_means_names <- c(FLEX1.4val_means_names, paste0("ALL_FLEX1.4val_",j))
}
FLEX1.4val_means <- c(mget("ALL_flex1.4val"), mget(FLEX1.4val_means_names))
rm(ALL_flex1.4val)
rm(FLEX1.4val_means_names)
for(j in strata_names){
  rm(list = paste0("ALL_FLEX1.4val_",j) )
}
##FLEX2
##2 noeuds

FLEX2.2val_means_names <- c()
for(j in strata_names){
  FLEX2.2val_means_names <- c(FLEX2.2val_means_names, paste0("ALL_FLEX2.2val_",j))
}
FLEX2.2val_means <- c(mget("ALL_flex2.2val"), mget(FLEX2.2val_means_names))
rm(ALL_flex2.2val)
rm(FLEX2.2val_means_names)
for(j in strata_names){
  rm(list = paste0("ALL_FLEX2.2val_",j) )
}

##4 noeuds

FLEX2.4val_means_names <- c()
for(j in strata_names){
  FLEX2.4val_means_names <- c(FLEX2.4val_means_names, paste0("ALL_FLEX2.4val_",j))
}
FLEX2.4val_means <- c(mget("ALL_flex2.4val"), mget(FLEX2.4val_means_names))
rm(ALL_flex2.4val)
rm(FLEX2.4val_means_names)
for(j in strata_names){
  rm(list = paste0("ALL_FLEX2.4val_",j) )
}
##WG
WGval_means_names <- c()
for(j in strata_names){
  WGval_means_names <- c(WGval_means_names, paste0("ALL_WGval_",j))
}
WGval_means <- c(mget("ALL_WGval"), mget(WGval_means_names))
rm(ALL_WGval)
rm(WGval_means_names)
for(j in strata_names){
  rm(list = paste0("ALL_WGval_",j) )
}

##PP
PPval_means_names <- c()
for(j in strata_names){
  PPval_means_names <- c(PPval_means_names, paste0("ALL_PPval_",j))
}
PPval_means <- c(mget("ALL_PPval"), mget(PPval_means_names))
rm(ALL_PPval)
rm(PPval_means_names)
for(j in strata_names){
  rm(list = paste0("ALL_PPval_",j) )
}
########
##     3000ROCSTA
################################################### ROC.net ##############################################################

library(RISCA)
##choix de l'espacement des cutoffs 
# c.off <- seq(0,1, by = 0.2)
# sort(unique(1-ind_estimP[, l + 1]))
# unique( quantile( 1-ind_estimP[, l + 1],  probs= seq(0,1, by = 0.01)) ) #pas en dessous de 20 points 0.025
####TRAIN
### WHOLE

##importation des bases de données / survies individuelles prédites par les modèles 
for(k in iterations){
  ##bases de données 
  assign(paste0("DATATrain_", k) , read.csv(paste0(path0, "DATAFRAMES/TRAIN/WHOLE/",k,"_dftrain.csv"),
                                            sep = ";") )
  
  ##survies individuelles 
  #plann
  assign(paste0("indTrainPLANN_", k) , read.csv(paste0(path0, "TRAIN/WHOLE/PLANN/individual_survival/",k,"_ind_PLANN.csv"),
                                                sep = ";") )
  #flex1
  assign(paste0("indTrainFLEX1.2_", k) , read.csv(paste0(path0, "TRAIN/WHOLE/FLEX1.2/individual_survival/",k,"_ind_flex12.csv"),
                                                  sep = ";") )
  assign(paste0("indTrainFLEX1.4_", k) , read.csv(paste0(path0, "TRAIN/WHOLE/FLEX1.4/individual_survival/",k,"_ind_flex14.csv"),
                                                  sep = ";") )
  #flex2
  assign(paste0("indTrainFLEX2.2_", k) , read.csv(paste0(path0, "TRAIN/WHOLE/FLEX2.2/individual_survival/",k,"_ind_flex22.csv"),
                                                  sep = ";") )
  assign(paste0("indTrainFLEX2.4_", k) , read.csv(paste0(path0, "TRAIN/WHOLE/FLEX2.4/individual_survival/",k,"_ind_flex24.csv"),
                                                  sep = ";") )
  #WG
  assign(paste0("indTrainWG_", k) , read.csv(paste0(path0, "TRAIN/WHOLE/WG/individual_survival/",k,"_ind_WG.csv"),
                                             sep = ";") )
}

##### calcul de ROC net

AUC_calc_W <- function(k) {
  datatrain <- get(paste0("DATATrain_", k))
  ind_estimP <- get(paste0("indTrainPLANN_", k))
  ind_estimF1.2 <- get(paste0("indTrainFLEX1.2_", k))
  ind_estimF1.4 <- get(paste0("indTrainFLEX1.4_", k))
  ind_estimF2.2 <- get(paste0("indTrainFLEX2.2_", k))
  ind_estimF2.4 <- get(paste0("indTrainFLEX2.4_", k))
  ind_estimWG <- get(paste0("indTrainWG_", k))
  
  rm(list = c(paste0("DATATrain_",k), paste0("indTrainPLANN_",k),
              paste0("indTrainFLEX1.2_",k), paste0("indTrainFLEX1.4_",k),
              paste0("indTrainFLEX2.2_",k), paste0("indTrainFLEX2.4_",k),
              paste0("indTrainWG_",k)))
  
  # Initialiser les vecteurs de résultats
  hold_P <- hold_F1.2 <- hold_F1.4 <- hold_F2.2 <- hold_F2.4 <- hold_WG <- c()
  # sort(unique(1-ind_estimP[, l + 1]))
  # unique( quantile( 1-ind_estimP[, l + 1],  probs= seq(0,1, by = 0.01)) ) #pas en dessous de 20 points 0.025
  
  for (l in seq_along(newtimes)) {
    hold_P    <- c(hold_P,    roc.net(datatrain$times, datatrain$status, 1-ind_estimP[,l + 1], datatrain$age, datatrain$sexchara,
                                      datatrain$year, slopop, pro.time = newtimes[l], 
                                      cut.off = unique( quantile(1-ind_estimP[,l+1],probs= seq(0,1, by = 0.01)) ) )$auc)
    hold_F1.2 <- c(hold_F1.2, roc.net(datatrain$times, datatrain$status, 1-ind_estimF1.2[, l], datatrain$age, datatrain$sexchara,
                                      datatrain$year, slopop, pro.time = newtimes[l],
                                      cut.off = unique( quantile( 1-ind_estimF1.2[, l],  probs= seq(0,1, by = 0.01)) ))$auc)
    hold_F1.4 <- c(hold_F1.4, roc.net(datatrain$times, datatrain$status, 1-ind_estimF1.4[, l], datatrain$age, datatrain$sexchara,
                                      datatrain$year, slopop, pro.time = newtimes[l],
                                      cut.off = unique( quantile( 1-ind_estimF1.4[, l],  probs= seq(0,1, by = 0.01)) ))$auc)
    hold_F2.2 <- c(hold_F2.2, roc.net(datatrain$times, datatrain$status, 1-ind_estimF2.2[, l], datatrain$age, datatrain$sexchara,
                                      datatrain$year, slopop, pro.time = newtimes[l], 
                                      cut.off = unique( quantile( 1-ind_estimF2.2[, l],  probs= seq(0,1, by = 0.01)) ))$auc)
    hold_F2.4 <- c(hold_F2.4, roc.net(datatrain$times, datatrain$status, 1-ind_estimF2.4[, l], datatrain$age, datatrain$sexchara,
                                      datatrain$year, slopop, pro.time = newtimes[l],
                                      cut.off = unique( quantile( 1-ind_estimF2.4[, l],  probs= seq(0,1, by = 0.01)) ))$auc)
    hold_WG   <- c(hold_WG,   roc.net(datatrain$times, datatrain$status, 1-ind_estimWG[, l],   datatrain$age, datatrain$sexchara,
                                      datatrain$year, slopop, pro.time = newtimes[l],
                                      cut.off = unique( quantile( 1-ind_estimWG[, l],  probs= seq(0,1, by = 0.01)) ))$auc)
  }
  
  return(list(
    P     = hold_P,
    F1.2  = hold_F1.2,
    F1.4  = hold_F1.4,
    F2.2  = hold_F2.2,
    F2.4  = hold_F2.4,
    WG    = hold_WG
  ))
}

results_list <- mclapply(iterations, AUC_calc_W, mc.cores = detectCores() - 4)

#supression des bases inutiles 
for(k in iterations){
  rm(list = c(paste0("DATATrain_",k), paste0("indTrainPLANN_",k),
              paste0("indTrainFLEX1.2_",k), paste0("indTrainFLEX1.4_",k),
              paste0("indTrainFLEX2.2_",k), paste0("indTrainFLEX2.4_",k),
              paste0("indTrainWG_",k)))
}
# Transformer en matrices pour chaque méthode :
AUC_WHOLE_P    <- do.call(rbind, lapply(results_list, function(x) x$P))
AUC_WHOLE_F1.2 <- do.call(rbind, lapply(results_list, function(x) x$F1.2))
AUC_WHOLE_F1.4 <- do.call(rbind, lapply(results_list, function(x) x$F1.4))
AUC_WHOLE_F2.2 <- do.call(rbind, lapply(results_list, function(x) x$F2.2))
AUC_WHOLE_F2.4 <- do.call(rbind, lapply(results_list, function(x) x$F2.4))
AUC_WHOLE_WG   <- do.call(rbind, lapply(results_list, function(x) x$WG))


### pour renommer les colonnes avec les temps de prognostic
for(l in c("P","F1.2","F1.4","F2.2","F2.4","WG")){
  dathold <- get(paste0("AUC_WHOLE_", l))
  colnames(dathold) <- c("1 year", "3 years", "5 years", "10 years")
  rownames(dathold) <- c(iterations)
  assign(paste0("AUC_WHOLE_", l),dathold)
}
## moyennes par temps 

for(l in c("P","F1.2","F1.4","F2.2","F2.4","WG")){
  
  assign(paste0("AUC_means_",l) ,  colMeans( get(paste0("AUC_WHOLE_",l)) ) )
  
}
################## STRATES

##importation des bases de données / survies individuelles prédites par les modèles 

for(k in iterations){
  ##bases de données 
  for (j in strata_names){
    assign(paste0("DATATrain_",j,"_", k) , read.csv(paste0(path0, "DATAFRAMES/TRAIN/",j,"/",k,"_dftrain_",j,".csv"),
                                                    sep = ";") )
    ##survies individuelles 
    #plann
    assign(paste0("indTrainPLANN_",j,"_", k) , read.csv(paste0(path0, "TRAIN/STRATA/PLANN/individual_survival/",k,"_ind_PLANN_",j,".csv"),
                                                        sep = ";") )
    #flex1
    assign(paste0("indTrainFLEX1.2_",j,"_", k) , read.csv(paste0(path0, "TRAIN/STRATA/FLEX1.2/individual_survival/",k,"_ind_flex12_",j,".csv"),
                                                          sep = ";") )
    assign(paste0("indTrainFLEX1.4_",j,"_", k) , read.csv(paste0(path0, "TRAIN/STRATA/FLEX1.4/individual_survival/",k,"_ind_flex14_",j,".csv"),
                                                          sep = ";") )
    #flex2
    assign(paste0("indTrainFLEX2.2_",j,"_", k) , read.csv(paste0(path0, "TRAIN/STRATA/FLEX2.2/individual_survival/",k,"_ind_flex22_",j,".csv"),
                                                          sep = ";") )
    assign(paste0("indTrainFLEX2.4_",j,"_", k) , read.csv(paste0(path0, "TRAIN/STRATA/FLEX2.4/individual_survival/",k,"_ind_flex24_",j,".csv"),
                                                          sep = ";") )
    #WG
    assign(paste0("indTrainWG_",j,"_", k) , read.csv(paste0(path0, "TRAIN/STRATA/WG/individual_survival/",k,"_ind_WG_",j,".csv"),
                                                     sep = ";") )
  }
}

##### calcul de ROC net

AUC_calc_S <- function(k) {
  auc_results <- list()
  
  for (j in strata_names) {
    datatrain <- get(paste0("DATATrain_", j, "_", k))
    ind_estimP <- get(paste0("indTrainPLANN_", j, "_", k))
    ind_estimF1.2 <- get(paste0("indTrainFLEX1.2_", j, "_", k))
    ind_estimF1.4 <- get(paste0("indTrainFLEX1.4_", j, "_", k))
    ind_estimF2.2 <- get(paste0("indTrainFLEX2.2_", j, "_", k))
    ind_estimF2.4 <- get(paste0("indTrainFLEX2.4_", j, "_", k))
    ind_estimWG <- get(paste0("indTrainWG_", j, "_", k))
    
    hold_P <- hold_F1.2 <- hold_F1.4 <- hold_F2.2 <- hold_F2.4 <- hold_WG <- c()
    
    for (l in seq_along(newtimes)) {
      hold_P    <- c(hold_P,    roc.net(datatrain$times, datatrain$status, 1-ind_estimP[,l + 1], datatrain$age, datatrain$sexchara,
                                        datatrain$year, slopop, pro.time = newtimes[l], 
                                        cut.off = unique( quantile(1-ind_estimP[,l+1],probs= seq(0,1, by = 0.01)) ) )$auc)
      hold_F1.2 <- c(hold_F1.2, roc.net(datatrain$times, datatrain$status, 1-ind_estimF1.2[, l], datatrain$age, datatrain$sexchara,
                                        datatrain$year, slopop, pro.time = newtimes[l],
                                        cut.off = unique( quantile( 1-ind_estimF1.2[, l],  probs= seq(0,1, by = 0.01)) ))$auc)
      hold_F1.4 <- c(hold_F1.4, roc.net(datatrain$times, datatrain$status, 1-ind_estimF1.4[, l], datatrain$age, datatrain$sexchara,
                                        datatrain$year, slopop, pro.time = newtimes[l],
                                        cut.off = unique( quantile( 1-ind_estimF1.4[, l],  probs= seq(0,1, by = 0.01)) ))$auc)
      hold_F2.2 <- c(hold_F2.2, roc.net(datatrain$times, datatrain$status, 1-ind_estimF2.2[, l], datatrain$age, datatrain$sexchara,
                                        datatrain$year, slopop, pro.time = newtimes[l], 
                                        cut.off = unique( quantile( 1-ind_estimF2.2[, l],  probs= seq(0,1, by = 0.01)) ))$auc)
      hold_F2.4 <- c(hold_F2.4, roc.net(datatrain$times, datatrain$status, 1-ind_estimF2.4[, l], datatrain$age, datatrain$sexchara,
                                        datatrain$year, slopop, pro.time = newtimes[l],
                                        cut.off = unique( quantile( 1-ind_estimF2.4[, l],  probs= seq(0,1, by = 0.01)) ))$auc)
      hold_WG   <- c(hold_WG,   roc.net(datatrain$times, datatrain$status, 1-ind_estimWG[, l],   datatrain$age, datatrain$sexchara,
                                        datatrain$year, slopop, pro.time = newtimes[l],
                                        cut.off = unique( quantile( 1-ind_estimWG[, l],  probs= seq(0,1, by = 0.01)) ))$auc)
    }
    
    auc_results[[j]] <- list(
      k = k,
      P = hold_P,
      F1.2 = hold_F1.2,
      F1.4 = hold_F1.4,
      F2.2 = hold_F2.2,
      F2.4 = hold_F2.4,
      WG = hold_WG
    )
  }
  
  return(auc_results)
}

res_all <- mclapply(iterations, AUC_calc_S, mc.cores = detectCores() - 4)

##supression des bases inutiles pour faire de la place
for(k in iterations){
  for(j in strata_names){
    rm(list = c(
      paste0("DATATrain_", j, "_", k),
      paste0("indTrainPLANN_", j, "_", k),
      paste0("indTrainFLEX1.2_", j, "_", k),
      paste0("indTrainFLEX1.4_", j, "_", k),
      paste0("indTrainFLEX2.2_", j, "_", k),
      paste0("indTrainFLEX2.4_", j, "_", k),
      paste0("indTrainWG_", j, "_", k)
    ))
  }
}

# Transformer en matrices pour chaque méthode :

for(j in strata_names){
  assign(paste0("AUC_",j,"_P"), do.call(rbind, lapply(res_all, function(x) x[[j]]$P)) )
  assign(paste0("AUC_",j,"_F1.2"), do.call(rbind, lapply(res_all, function(x) x[[j]]$F1.2)) )
  assign(paste0("AUC_",j,"_F1.4"), do.call(rbind, lapply(res_all, function(x) x[[j]]$F1.4)) )
  assign(paste0("AUC_",j,"_F2.2"), do.call(rbind, lapply(res_all, function(x) x[[j]]$F2.2)) )
  assign(paste0("AUC_",j,"_F2.4"), do.call(rbind, lapply(res_all, function(x) x[[j]]$F2.4)) )
  assign(paste0("AUC_",j,"_WG"), do.call(rbind, lapply(res_all, function(x) x[[j]]$WG)) )
}
### pour renommer les colonnes avec les temps de prognostic
for(l in c("P","F1.2","F1.4","F2.2","F2.4","WG")){
  for(j in strata_names){
    dathold <- get(paste0("AUC_", j,"_",l))
    colnames(dathold) <- c("1 year", "3 years", "5 years", "10 years")
    rownames(dathold) <- c(iterations)
    assign(paste0("AUC_",j,"_", l),dathold)
  }
}
## moyennes par temps 

for(l in c("P","F1.2","F1.4","F2.2","F2.4","WG")){
  for(j in strata_names){
    assign(paste0("AUC_means_",j,"_",l) ,  colMeans( get(paste0("AUC_",j,"_",l)) ) )
  }
}


### sous forme de liste pour libérer de l'espace dans l'evt (pour les moyennes)

ROCmean_results <- list()

for(l in c("P","F1.2","F1.4","F2.2","F2.4","WG")){
  ROCmean_results[['WHOLE']][[l]] <- list()
  ROCmean_results[['WHOLE']][[l]] <- get( paste0("AUC_means_",l) )
  rm(list = paste0("AUC_means_",l))
  
  for(j in strata_names){
    ROCmean_results[[j]][[l]] <- get(paste0("AUC_means_", j, "_", l)) 
    rm(list = paste0("AUC_means_", j, "_", l))
  }
  
}
### meme chose pour les valeurs pour chaque itération 


ROC_results <- list()

for(l in c("P","F1.2","F1.4","F2.2","F2.4","WG")){
  ROC_results[['WHOLE']][[l]] <- list()
  ROC_results[['WHOLE']][[l]] <- get( paste0("AUC_WHOLE_",l) )
  rm(list = paste0("AUC_WHOLE_",l) )
  for(j in strata_names){
    ROC_results[[j]][[l]] <- get(paste0("AUC_", j, "_", l)) 
    rm(list = paste0("AUC_", j, "_", l))
  }
  
}
#################### VALIDATION ######################

### WHOLE

##importation des bases de données / survies individuelles prédites par les modèles 
for(k in iterations){
  ##bases de données 
  assign(paste0("DATAValid_", k) , read.csv(paste0(path0, "DATAFRAMES/VALID/WHOLE/",k,"_dfvalid.csv"),
                                            sep = ";") )
  
  ##survies individuelles 
  #plann
  assign(paste0("indValidPLANN_", k) , read.csv(paste0(path0, "VALID/WHOLE/PLANN/individual_survival/",k,"_ind_PLANNval.csv"),
                                                sep = ";") )
  
  #flex1
  assign(paste0("indValidFLEX1.2_",k) , read.csv(paste0(path0, "VALID/WHOLE/FLEX1.2/individual_survival/",k,"_ind_flexval12.csv"),
                                                 sep = ";") )
  assign(paste0("indValidFLEX1.4_",k) , read.csv(paste0(path0, "VALID/WHOLE/FLEX1.4/individual_survival/",k,"_ind_flexval14.csv"),
                                                 sep = ";") )
  #flex2
  assign(paste0("indValidFLEX2.2_",k) , read.csv(paste0(path0, "VALID/WHOLE/FLEX2.2/individual_survival/",k,"_ind_flexval22.csv"),
                                                 sep = ";") )
  assign(paste0("indValidFLEX2.4_",k) , read.csv(paste0(path0, "VALID/WHOLE/FLEX2.4/individual_survival/",k,"_ind_flexval24.csv"),
                                                 sep = ";") )
  #WG
  assign(paste0("indValidWG_", k) , read.csv(paste0(path0, "VALID/WHOLE/WG/individual_survival/",k,"_ind_flexWGval.csv"),
                                             sep = ";") )
}

##### calcul de ROC net

##creation d'objets pour enregistrer les AUC de chaque itération aux différents temps et pour les différents modèles 

AUC_calc_val_W <- function(k){ 
  
  dataval <- get(paste0("DATAValid_",k))
  ind_estimP <- get(paste0("indValidPLANN_",k))
  ind_estimF1.2 <- get(paste0("indValidFLEX1.2_",k))
  ind_estimF1.4 <- get(paste0("indValidFLEX1.4_", k))
  ind_estimF2.2 <- get(paste0("indValidFLEX2.2_", k))
  ind_estimF2.4 <- get(paste0("indValidFLEX2.4_", k))
  ind_estimWG <- get(paste0("indValidWG_",k))
  
  rm(list = c(paste0("DATAValid_",k), paste0("indValidPLANN_",k),paste0("indValidFLEX1.2_",k), paste0("indValidFLEX1.4_", k),
              paste0("indValidFLEX2.2_", k),paste0("indValidFLEX2.4_", k),paste0("indValidWG_",k)))
  
  # Initialiser les vecteurs de résultats
  hold_P <- hold_F1.2 <- hold_F1.4 <- hold_F2.2 <- hold_F2.4 <- hold_WG <- c()
  
  ###PLANN
  for (l in seq_along(newtimes)) {
    hold_P    <- c(hold_P,    roc.net(dataval$times, dataval$status, 1-ind_estimP[,l + 1], dataval$age, dataval$sexchara,
                                      dataval$year, slopop, pro.time = newtimes[l], 
                                      cut.off = unique( quantile(1-ind_estimP[,l+1],probs= seq(0,1, by = 0.01)) ) )$auc)
    hold_F1.2 <- c(hold_F1.2, roc.net(dataval$times, dataval$status, 1-ind_estimF1.2[, l], dataval$age, dataval$sexchara,
                                      dataval$year, slopop, pro.time = newtimes[l],
                                      cut.off = unique( quantile( 1-ind_estimF1.2[, l],  probs= seq(0,1, by = 0.01)) ))$auc)
    hold_F1.4 <- c(hold_F1.4, roc.net(dataval$times, dataval$status, 1-ind_estimF1.4[, l], dataval$age, dataval$sexchara,
                                      dataval$year, slopop, pro.time = newtimes[l],
                                      cut.off = unique( quantile( 1-ind_estimF1.4[, l],  probs= seq(0,1, by = 0.01)) ))$auc)
    hold_F2.2 <- c(hold_F2.2, roc.net(dataval$times, dataval$status, 1-ind_estimF2.2[, l], dataval$age, dataval$sexchara,
                                      dataval$year, slopop, pro.time = newtimes[l], 
                                      cut.off = unique( quantile( 1-ind_estimF2.2[, l],  probs= seq(0,1, by = 0.01)) ))$auc)
    hold_F2.4 <- c(hold_F2.4, roc.net(dataval$times, dataval$status, 1-ind_estimF2.4[, l], dataval$age, dataval$sexchara,
                                      dataval$year, slopop, pro.time = newtimes[l],
                                      cut.off = unique( quantile( 1-ind_estimF2.4[, l],  probs= seq(0,1, by = 0.01)) ))$auc)
    hold_WG   <- c(hold_WG,   roc.net(dataval$times, dataval$status, 1-ind_estimWG[, l],   dataval$age, dataval$sexchara,
                                      dataval$year, slopop, pro.time = newtimes[l],
                                      cut.off = unique( quantile( 1-ind_estimWG[, l],  probs= seq(0,1, by = 0.01)) ))$auc)
  }
  
  return(list(
    P     = hold_P,
    F1.2  = hold_F1.2,
    F1.4  = hold_F1.4,
    F2.2  = hold_F2.2,
    F2.4  = hold_F2.4,
    WG    = hold_WG
  ))
}

resultval_list <- mclapply(iterations, AUC_calc_val_W, mc.cores = detectCores() - 4)

#supression des bases inutiles 
for(k in iterations){
  rm(list = c(paste0("DATAValid_",k), paste0("indValidPLANN_",k),
              paste0("indValidFLEX1.2_",k), paste0("indValidFLEX1.4_",k),
              paste0("indValidFLEX2.2_",k), paste0("indValidFLEX2.4_",k),
              paste0("indValidWG_",k)))
}
# Transformer en matrices pour chaque méthode :
AUCval_WHOLE_P    <- do.call(rbind, lapply(resultval_list, function(x) x$P))
AUCval_WHOLE_F1.2 <- do.call(rbind, lapply(resultval_list, function(x) x$F1.2))
AUCval_WHOLE_F1.4 <- do.call(rbind, lapply(resultval_list, function(x) x$F1.4))
AUCval_WHOLE_F2.2 <- do.call(rbind, lapply(resultval_list, function(x) x$F2.2))
AUCval_WHOLE_F2.4 <- do.call(rbind, lapply(resultval_list, function(x) x$F2.4))
AUCval_WHOLE_WG   <- do.call(rbind, lapply(resultval_list, function(x) x$WG))

### pour renommer les colonnes avec les temps de prognostic
for(l in c("P","F1.2","F1.4","F2.2","F2.4","WG")){
  dathold <- get(paste0("AUCval_WHOLE_", l))
  colnames(dathold) <- c("1 year", "3 years", "5 years", "10 years")
  rownames(dathold) <- c(iterations)
  assign(paste0("AUCval_WHOLE_", l),dathold)
}
## moyennes par temps 

for(l in c("P","F1.2","F1.4","F2.2","F2.4","WG")){
  
  assign(paste0("AUCval_means_",l) ,  colMeans( get(paste0("AUCval_WHOLE_",l)) ) )
  
}
################## STRATES

##importation des bases de données / survies individuelles prédites par les modèles 

for(k in iterations){
  ##bases de données 
  for (j in strata_names){
    assign(paste0("DATAValid_",j,"_", k) , read.csv(paste0(path0, "DATAFRAMES/VALID/",j,"/",k,"_dfvalid_",j,".csv"),
                                                    sep = ";") )
    ##survies individuelles 
    #plann
    assign(paste0("indValidPLANN_",j,"_", k) , read.csv(paste0(path0, "VALID/STRATA/PLANN/individual_survival/",k,"_ind_PLANNval_",j,".csv"),
                                                        sep = ";") )
    #flex1
    assign(paste0("indValidFLEX1.2_",j,"_", k) , read.csv(paste0(path0, "VALID/STRATA/FLEX1.2/individual_survival/",k,"_ind_flexval12_",j,".csv"),
                                                          sep = ";") )
    assign(paste0("indValidFLEX1.4_",j,"_", k) , read.csv(paste0(path0, "VALID/STRATA/FLEX1.4/individual_survival/",k,"_ind_flexval14_",j,".csv"),
                                                          sep = ";") )
    #flex2
    assign(paste0("indValidFLEX2.2_",j,"_", k) , read.csv(paste0(path0, "VALID/STRATA/FLEX2.2/individual_survival/",k,"_ind_flexval22_",j,".csv"),
                                                          sep = ";") )
    assign(paste0("indValidFLEX2.4_",j,"_", k) , read.csv(paste0(path0, "VALID/STRATA/FLEX2.4/individual_survival/",k,"_ind_flexval24_",j,".csv"),
                                                          sep = ";") )
    #WG
    assign(paste0("indValidWG_",j,"_", k) , read.csv(paste0(path0, "VALID/STRATA/WG/individual_survival/",k,"_ind_WGval_",j,".csv"),
                                                     sep = ";") )
  }
}

##### calcul de ROC net

AUC_calc_val_S <- function(k){ 
  auc_results <- list()
  for(j in strata_names){
    dataval <- get(paste0("DATAValid_",j,"_", k))
    ind_estimP <- get(paste0("indValidPLANN_",j,"_", k))
    ind_estimF1.2 <- get(paste0("indValidFLEX1.2_",j,"_", k))
    ind_estimF1.4 <- get(paste0("indValidFLEX1.4_",j,"_", k))
    ind_estimF2.2 <- get(paste0("indValidFLEX2.2_",j,"_", k))
    ind_estimF2.4 <- get(paste0("indValidFLEX2.4_",j,"_", k))
    ind_estimWG <- get(paste0("indValidWG_",j,"_", k))
    
    rm(list = c(paste0("DATAValid_",j,"_", k), paste0("indValidPLANN_",j,"_", k), paste0("indValidFLEX1.2_",j,"_", k),paste0("indValidFLEX1.4_",j,"_", k),
                paste0("indValidFLEX2.2_",j,"_", k), paste0("indValidFLEX2.2_",j,"_", k), paste0("indValidWG_",j,"_", k)))
    
    hold_P <- hold_F1.2 <- hold_F1.4 <- hold_F2.2 <- hold_F2.4 <- hold_WG <- c()
    
    for (l in seq_along(newtimes)) {
      hold_P    <- c(hold_P,    roc.net(dataval$times, dataval$status, 1-ind_estimP[,l + 1], dataval$age, dataval$sexchara,
                                        dataval$year, slopop, pro.time = newtimes[l], 
                                        cut.off = unique( quantile(1-ind_estimP[,l+1],probs= seq(0,1, by = 0.01)) ) )$auc)
      hold_F1.2 <- c(hold_F1.2, roc.net(dataval$times, dataval$status, 1-ind_estimF1.2[, l], dataval$age, dataval$sexchara,
                                        dataval$year, slopop, pro.time = newtimes[l],
                                        cut.off = unique( quantile( 1-ind_estimF1.2[, l],  probs= seq(0,1, by = 0.01)) ))$auc)
      hold_F1.4 <- c(hold_F1.4, roc.net(dataval$times, dataval$status, 1-ind_estimF1.4[, l], dataval$age, dataval$sexchara,
                                        dataval$year, slopop, pro.time = newtimes[l],
                                        cut.off = unique( quantile( 1-ind_estimF1.4[, l],  probs= seq(0,1, by = 0.01)) ))$auc)
      hold_F2.2 <- c(hold_F2.2, roc.net(dataval$times, dataval$status, 1-ind_estimF2.2[, l], dataval$age, dataval$sexchara,
                                        dataval$year, slopop, pro.time = newtimes[l], 
                                        cut.off = unique( quantile( 1-ind_estimF2.2[, l],  probs= seq(0,1, by = 0.01)) ))$auc)
      hold_F2.4 <- c(hold_F2.4, roc.net(dataval$times, dataval$status, 1-ind_estimF2.4[, l], dataval$age, dataval$sexchara,
                                        dataval$year, slopop, pro.time = newtimes[l],
                                        cut.off = unique( quantile( 1-ind_estimF2.4[, l],  probs= seq(0,1, by = 0.01)) ))$auc)
      hold_WG   <- c(hold_WG,   roc.net(dataval$times, dataval$status, 1-ind_estimWG[, l],   dataval$age, dataval$sexchara,
                                        dataval$year, slopop, pro.time = newtimes[l],
                                        cut.off = unique( quantile( 1-ind_estimWG[, l],  probs= seq(0,1, by = 0.01)) ))$auc)
    }
    
    auc_results[[j]] <- list(
      k = k,
      P = hold_P,
      F1.2 = hold_F1.2,
      F1.4 = hold_F1.4,
      F2.2 = hold_F2.2,
      F2.4 = hold_F2.4,
      WG = hold_WG
    )
  }
  
  return(auc_results)
}

resval_all <- mclapply(iterations, AUC_calc_val_S, mc.cores = detectCores() - 4)

for(k in iterations){
  for(j in strata_names){
    rm(list = c(
      paste0("DATAValid_", j, "_", k),
      paste0("indValidPLANN_", j, "_", k),
      paste0("indValidFLEX1.2_", j, "_", k),
      paste0("indValidFLEX1.4_", j, "_", k),
      paste0("indValidFLEX2.2_", j, "_", k),
      paste0("indValidFLEX2.4_", j, "_", k),
      paste0("indValidWG_", j, "_", k)
    ))
  }
}

# Transformer en matrices pour chaque méthode :

for(j in strata_names){
  assign(paste0("AUCval_",j,"_P"), do.call(rbind, lapply(resval_all, function(x) x[[j]]$P)) )
  assign(paste0("AUCval_",j,"_F1.2"), do.call(rbind, lapply(resval_all, function(x) x[[j]]$F1.2)) )
  assign(paste0("AUCval_",j,"_F1.4"), do.call(rbind, lapply(resval_all, function(x) x[[j]]$F1.4)) )
  assign(paste0("AUCval_",j,"_F2.2"), do.call(rbind, lapply(resval_all, function(x) x[[j]]$F2.2)) )
  assign(paste0("AUCval_",j,"_F2.4"), do.call(rbind, lapply(resval_all, function(x) x[[j]]$F2.4)) )
  assign(paste0("AUCval_",j,"_WG"), do.call(rbind, lapply(resval_all, function(x) x[[j]]$WG)) )
}
### pour renommer les colonnes avec les temps de prognostic
for(l in c("P","F1.2","F1.4","F2.2","F2.4","WG")){
  for(j in strata_names){
    dathold <- get(paste0("AUCval_", j,"_",l))
    colnames(dathold) <- c("1 year", "3 years", "5 years", "10 years")
    rownames(dathold) <- c(iterations)
    assign(paste0("AUCval_",j,"_", l),dathold)
  }
}
## moyennes par temps 

for(l in c("P","F1.2","F1.4","F2.2","F2.4","WG")){
  for(j in strata_names){
    assign(paste0("AUCval_means_",j,"_",l) ,  colMeans( get(paste0("AUCval_",j,"_",l)) ) )
  }
}


### sous forme de liste pour libérer de l'espace dans l'evt (pour les moyennes)

ROCmeanval_results <- list()

for(l in c("P","F1.2","F1.4","F2.2","F2.4","WG")){
  ROCmeanval_results[['WHOLE']][[l]] <- list()
  ROCmeanval_results[['WHOLE']][[l]] <- get( paste0("AUCval_means_",l) )
  rm(list = paste0("AUCval_means_",l))
  
  for(j in strata_names){
    ROCmeanval_results[[j]][[l]] <- get(paste0("AUCval_means_", j, "_", l)) 
    rm(list = paste0("AUCval_means_", j, "_", l))
  }
  
}
### meme chose pour les valeurs pour chaque itération 


ROCval_results <- list()

for(l in c("P","F1.2","F1.4","F2.2","F2.4","WG")){
  ROCval_results[['WHOLE']][[l]] <- list()
  ROCval_results[['WHOLE']][[l]] <- get( paste0("AUCval_WHOLE_",l) )
  rm(list = paste0("AUCval_WHOLE_",l) )
  for(j in strata_names){
    ROCval_results[[j]][[l]] <- get(paste0("AUCval_", j, "_", l)) 
    rm(list = paste0("AUCval_", j, "_", l))
  }
  
}

print(Sys.time()-start)


save.image("~/Documents/Simulations/Simulations mai 2025/Résultats/1000ite_3000ind_NPH_HC_14_05.Rdata")


rm(list = ls(envir = .GlobalEnv), envir = .GlobalEnv)

##3000IndEND
##################################################################################################
##################################################################################################
##################################################################################################
                                        ##5000##
##################################################################################################
##################################################################################################
##################################################################################################

#5000SiSTA
library(survivalNET)
library(survivalPLANN)
library(parallel)
library(doParallel)


path0 <- "~/Documents/Simulations/Simulations mai 2025/N5000/Output simulations_N5000_NPH_HC/"
path1 <- "~/Documents/Simulations/BASES/"
############

#importation des données
data("fr.ratetable")
data("slopop")
data(colrec)

# data management


colrec$agey <- colrec$age/365.241

colrec$sexchara <- "male"
colrec$sexchara[colrec$sex==2] <- "female"

colrec$colon <- 1*(as.character(colrec$site)=="colon")

colrec <- colrec[colrec$stage!=99,] # suppresion des patients avec valeur manquante
table(colrec$stage)

colrec$stage1 <- 1*(colrec$stage==1)
colrec$stage2 <- 1*(colrec$stage==2)
colrec$stage3 <- 1*(colrec$stage==3)

table(colrec$stat)

colrec$agey10 <- colrec$agey/10
names(colrec)


###################################################################################################
### Fonctions de survie nette, attendue, observée et de simulation de temps
###################################################################################################


Sn <- function(time, sigma, nu, theta, beta, covariates)
{
  exp((1-(1+(time/(exp(sigma)))^exp(nu))^(1/exp(theta))) * exp(sum(covariates*beta)) )
}


simulate_iteration <- function(i){
  
  set.seed(i)
  
  start <- Sys.time()
  
  data <- read.csv(file = paste0(path1,"ind5000NPH_HC/",i,"_df5000.csv"), sep = ";")
  Tsigma <- 12.6 
  Tnu <- -0.5
  Ttheta <- 0
  betaZ <- c(0.9,2.7,0.3,-0.1,-0.1)
  
  TsigmaHC <- 13.0
  TnuHC <- -0.55
  TthetaHC <- 0
  betaZHC <- c(1,2.9,0.2)
  
  #### Hommes Rectum
  TsigmaHR <- 12.75
  TnuHR <- -0.3
  TthetaHR <- 0                      
  betaZHR <- c(1.2,3,0.35)
  
  #### Femmes Colon
  TsigmaFC <- 14
  TnuFC <- -0.55
  TthetaFC <- 0
  betaZFC <- c(0.6,2.6,-0.1)
  
  #### Femmes Rectum
  TsigmaFR <- 11.5
  TnuFR <- -0.4
  TthetaFR <- 0
  betaZFR <- c(0.9,2.7,0.16)
  ############### validation sample
  N <- dim(data)[1]
  N_train <- N/2 #500
  data_train = data[1:N_train,]
  data_valid = data[(N_train+1):dim(data)[1],]
  
  ###############################
  ########## STRATES ############
  ###############################
  
  ### training
  #################### Hommes ####################
  
  ### hommes & colon
  
  data_train_HC <- data_train[data_train$sex == 1 & data_train$colon == 1,]
  
  ### hommes & rectum
  
  data_train_HR <- data_train[data_train$sex == 1 & data_train$colon == 0,]
  
  #################### Femmes ####################
  
  ### femmes & colon
  
  data_train_FC <- data_train[data_train$sex == 2 & data_train$colon == 1,]
  
  ### femmes & rectum
  
  data_train_FR <- data_train[data_train$sex == 2 & data_train$colon == 0,]
  
  
  ########################################
  ### validation
  
  #################### Hommes ####################
  
  ### hommes & colon
  
  data_valid_HC <- data_valid[data_valid$sex == 1 & data_valid$colon == 1,]
  
  ### hommes & rectum
  
  data_valid_HR <- data_valid[data_valid$sex == 1 & data_valid$colon == 0,]
  
  #################### Femmes ####################
  
  ### femmes & colon
  
  data_valid_FC <- data_valid[data_valid$sex == 2 & data_valid$colon == 1,]
  
  ### femmes & rectum
  
  data_valid_FR <- data_valid[data_valid$sex == 2 & data_valid$colon == 0,]
  
  ########################################
  
  ######################valeurs issues du modèle theorique
  
  pro.time <- max(data_train$times[data_train$status==1])
  
  ###temps d'évaluations (t =1, 3, 5, 10 ans)
  
  newtimes <- c(365.241, 365.241*3, 365.241*5, 365.241*10) 
  # @THOMAS : 80% 50% 20% sujets à risque 
  
  Zdata_train <- cbind(data_train$stage2, data_train$stage3, data_train$agey10, data_train$sex, data_train$colon)
  colnames(Zdata_train) <- c("stage2", "stage3", "agey10", "sex", "colon")
  
  funtheo <- function(i){
    if(Zdata_train[i,"sex"] == 1 & Zdata_train[i,"colon"] == 1){
      Tsigma = TsigmaHC
      Tnu = TnuHC
      Ttheta = TthetaHC
      betaZ = betaZHC
    }
    # Hommes rectum
    if(Zdata_train[i,"sex"] == 1 & Zdata_train[i,"colon"] == 0){
      Tsigma = TsigmaHR
      Tnu = TnuHR
      Ttheta = TthetaHR
      betaZ = betaZHR
    }
    # Femmes colon
    if(Zdata_train[i,"sex"] == 2 & Zdata_train[i,"colon"] == 1){
      Tsigma = TsigmaFC
      Tnu = TnuFC
      Ttheta = TthetaFC
      betaZ = betaZFC
    }
    # Femmes rectum
    if(Zdata_train[i,"sex"] == 2 & Zdata_train[i,"colon"] == 0){
      Tsigma = TsigmaFR
      Tnu = TnuFR
      Ttheta = TthetaFR
      betaZ = betaZFR
    }
    Sn(time = newtimes, sigma = Tsigma, nu = Tnu, theta = Ttheta, beta = betaZ, covariates = Zdata_train[i,-c(4,5)]) }
  
  theopred <- t(sapply(1:dim(Zdata_train)[1], FUN = funtheo ))
  
  colnames(theopred) <- newtimes
  
  mean.theo <- apply(theopred, FUN="mean", MARGIN=2)
  
  ###strates
  #####
  
  strata_names <- c("HC", "HR", "FC", "FR")
  
  for(j in strata_names){
    
    data_train_var <- get(paste0("data_train_", j))
    
    Zdata_train <- cbind(data_train_var$stage2, data_train_var$stage3, data_train_var$agey10, data_train_var$sex, data_train_var$colon)
    colnames(Zdata_train) <- c("stage2", "stage3", "agey10", "sex", "colon")
    
    funtheostrat <- function(i){
      if(Zdata_train[i,"sex"] == 1 & Zdata_train[i,"colon"] == 1){
        Tsigma = TsigmaHC
        Tnu = TnuHC
        Ttheta = TthetaHC
        betaZ = betaZHC
      }
      # Hommes rectum
      if(Zdata_train[i,"sex"] == 1 & Zdata_train[i,"colon"] == 0){
        Tsigma = TsigmaHR
        Tnu = TnuHR
        Ttheta = TthetaHR
        betaZ = betaZHR
      }
      # Femmes colon
      if(Zdata_train[i,"sex"] == 2 & Zdata_train[i,"colon"] == 1){
        Tsigma = TsigmaFC
        Tnu = TnuFC
        Ttheta = TthetaFC
        betaZ = betaZFC
      }
      # Femmes rectum
      if(Zdata_train[i,"sex"] == 2 & Zdata_train[i,"colon"] == 0){
        Tsigma = TsigmaFR
        Tnu = TnuFR
        Ttheta = TthetaFR
        betaZ = betaZFR
      }
      Sn(time = newtimes, sigma = Tsigma, nu = Tnu, theta = Ttheta, beta = betaZ, covariates = Zdata_train[i,-c(4,5)]) }
    
    theopredstrat <- t(sapply(1:dim(Zdata_train)[1], FUN = funtheostrat ))
    
    colnames(theopredstrat) <- newtimes
    
    mean.theo_strat <- apply(theopredstrat, FUN = "mean", MARGIN = 2)
    
    assign(paste0("Zdata_train_", j), Zdata_train)
    assign(paste0("theopred_", j), theopredstrat)
    assign(paste0("mean.theo_", j), mean.theo_strat)  
    
  }
  
  #####
  
  ####### sur la base de validation
  
  Zdata_val <- cbind(data_valid$stage2, data_valid$stage3, data_valid$agey10, data_valid$sex, data_valid$colon)
  
  theopredval <- t(sapply(1:dim(Zdata_val)[1], FUN = function(i){
    Sn(time = newtimes, sigma = Tsigma, nu = Tnu, theta = Ttheta, beta = betaZ,  covariates = Zdata_val[i,]) } ) )
  
  colnames(theopredval) <- newtimes
  
  mean.theoval <- apply(theopredval, FUN="mean", MARGIN=2) 
  
  ###strates
  #####
  
  for(j in strata_names){
    
    data_valid_var <- get(paste0("data_valid_", j))
    
    Zdata_valid <- cbind(data_valid_var$stage2, data_valid_var$stage3, data_valid_var$agey10, data_valid_var$sex, data_valid_var$colon)
    
    theopredvalstrat <- t(sapply(1:nrow(Zdata_valid), FUN = function(i) {
      Sn(time = newtimes, sigma = Tsigma, nu = Tnu, theta = Ttheta, beta = betaZ, covariates = Zdata_valid[i,])
    }))
    
    colnames(theopredvalstrat) <- newtimes
    
    mean.theoval_strat <- apply(theopredvalstrat, FUN = "mean", MARGIN = 2)
    
    assign(paste0("Zdata_val_", j), Zdata_valid)
    assign(paste0("theopredval_", j), theopredvalstrat)
    assign(paste0("mean.theoval_", j), mean.theoval_strat)  
    
  }
  
  
  #####
  #######################################################
  ###### choix des hyper-paramètres par cross-validation
  #######################################################
  
  
  # 
  # ## Modèle PLANN
  # # 
  
  startCVplann <- Sys.time()
  
  tune.plann <- cvPLANN(formula = Surv(times, status) ~ stage2 + stage3 + agey10 + sex01 + colon, pro.time = pro.time,
                        data = data_train, cv = 10, inter= 365.241/12, size=c(2, 4, 8, 10),
                        decay=c(0.01, 0.1), maxit=1000, MaxNWts=10000)
  
  TimeCVplann <- Sys.time() - startCVplann
  
  # 
  # ####################
  # #paramètres optimaux
  # ####################
  # 
  
  # inter <- 365.241/12
  # size <- 4
  # decay <- 0.01
  # maxit <- 1000
  # MaxNWts <- 10000
  
  ##### PLANN
  
  inter <- tune.plann$optimal$inter
  
  size <- tune.plann$optimal$size
  
  decay <- tune.plann$optimal$decay
  
  maxit <- tune.plann$optimal$maxit
  
  MaxNWts <- tune.plann$optimal$MaxNWts  
  
  # ##### FLEXNET
  
  #flex1 
  
  m1.2 <- 2
  m1.4 <- 4
  failedindic1.2 <- NULL
  failedindic1.4 <- NULL
  
  #flex2  
  
  m2.2 <- 2
  m2.4 <- 4
  failedindic2.2 <- NULL
  failedindic2.4 <- NULL
  
  
  ##################
  # Modèles fittés
  ##################
  # 
  
  ### plann
  
  plann.model <- sPLANN(formula = Surv(times, status) ~ stage2 + stage3 + agey10 + sex01 + colon,
                        pro.time = pro.time, data = data_train, inter= inter, size=size, 
                        decay=decay, maxit=maxit, MaxNWts=MaxNWts)
  
  ### Weibull Généralisé 
  
  WG.model <- survivalNET(formula = Surv(times, status) ~ stage2 + stage3 + agey10 + sex01 + colon +
                            ratetable(age, year, sexchara), data = data_train,
                          ratetable=slopop, dist = "genweibull")
  
  logcoeffWG <- tail(WG.model$coefficients, 3)
  
  ################## flex1
  ########### 2 noeuds
  m1.2_original <- m1.2
  m1.2_changed <- FALSE
  
  converged <- FALSE
  
  while (m1.2 >= 0 && !converged) {
    tryCatch({
      flex.model1.2 <- survivalFLEXNET(
        formula = Surv(times, status) ~ stage2 + stage3 + agey10 + sex01 + colon +
          ratetable(age, year, sexchara),
        data = data_train,
        ratetable = slopop,
        m = m1.2
      )
      # If the model succeeds, set converged to TRUE and break out of the loop
      converged <- TRUE
    }, error = function(e) {
      # If an error occurs, decrement m1
      m1.2_changed <<- TRUE
      m1.2 <<- m1.2 - 1
    })
  }
  
  if (m1.2_changed & converged) {
    failedindic1.2 <- paste0( " : this value was decreased from ", m1.2_original)
  }
  
  if (!converged) {
    flex.model1.2 <- WG.model
    failedindic1.2 <- "flex.model1.2 didn't converge, even after decreasing the number of knots incrementaly. WG.model was used instead."
  }
  
  ########### 4 noeuds
  m1.4_original <- m1.4
  m1.4_changed <- FALSE
  
  converged <- FALSE
  
  while (m1.4 >= 0 && !converged) {
    tryCatch({
      flex.model1.4 <- survivalFLEXNET(
        formula = Surv(times, status) ~ stage2 + stage3 + agey10 + sex01 + colon +
          ratetable(age, year, sexchara),
        data = data_train,
        ratetable = slopop,
        m = m1.4
      )
      # If the model succeeds, set converged to TRUE and break out of the loop
      converged <- TRUE
    }, error = function(e) {
      # If an error occurs, decrement m1
      m1.4_changed <<- TRUE
      m1.4 <<- m1.4 - 1
    })
  }
  
  if (m1.4_changed & converged) {
    failedindic1.4 <- paste0( " : this value was decreased from ", m1.4_original)
  }
  
  if (!converged) {
    flex.model1.4 <- WG.model
    failedindic1.4 <- "flex.model1.4 didn't converge, even after decreasing the number of knots incrementaly. WG.model was used instead."
  }
  ### flex2
  ## 2 noeuds 
  m2.2_original <- m2.2
  m2.2_changed <- FALSE
  
  converged <- FALSE
  
  while (m2.2 >= 0 && !converged) {
    tryCatch({flex.model2.2 <-survivalFLEXNET(formula = Surv(times, status) ~ stage2 + stage3 + agey10 + strata(sex.organ) +
                                                ratetable(age, year, sexchara), data = data_train,
                                              ratetable=slopop, m = m2.2)
    # If the model succeeds, set converged to TRUE and break out of the loop
    converged <- TRUE
    }, error = function(e) {
      # If an error occurs, decrement m2
      m2.2_changed <<- TRUE
      m2.2 <<- m2.2 - 1
    })
  }
  
  if (m2.2_changed & converged) {
    failedindic2.2 <- paste0(" this value was decreased from ", m2.2_original)
  }
  
  if (!converged) {
    flex.model2.2 <- flex.model1.2
    failedindic2.2 <- "flex.model2.2 didn't converge, even after decreasing the number of knots incrementaly. flex.model1.2 was used instead."
    
  }
  ## 4 noeuds 
  
  m2.4_original <- m2.4
  m2.4_changed <- FALSE
  
  converged <- FALSE
  
  while (m2.4 >= 0 && !converged) {
    tryCatch({flex.model2.4 <-survivalFLEXNET(formula = Surv(times, status) ~ stage2 + stage3 + agey10 + strata(sex.organ) +
                                                ratetable(age, year, sexchara), data = data_train,
                                              ratetable=slopop, m = m2.4)
    # If the model succeeds, set converged to TRUE and break out of the loop
    converged <- TRUE
    }, error = function(e) {
      # If an error occurs, decrement m2
      m2.4_changed <<- TRUE
      m2.4 <<- m2.4 - 1
    })
  }
  
  if (m2.4_changed & converged) {
    failedindic2.4 <- paste0(" this value was decreased from ", m2.4_original)
  }
  
  if (!converged) {
    flex.model2.4 <- flex.model1.4
    failedindic2.4 <- "flex.model2.4 didn't converge, even after decreasing the number of knots incrementaly. flex.model1.4 was used instead."
    
  }
  
  ##############
  #predictions #
  ##############
  
  #bases entières (entraînement et validation)
  ##################################################################################################
  ########### BASE DE TRAIN
  ###plann
  
  plannpred <- predictRS(plann.model, data = data_train, newtimes = newtimes, ratetable=slopop, 
                         age = "age", year = "year", sex = "sexchara")
  
  mean.plann <- plannpred$mpredictions$net_survival 
  
  ###flexnet
  
  ### modèle PH
  ##2 noeuds
  
  flexpred1.2 <- predict(flex.model1.2, newtimes = newtimes)$predictions
  
  mean.flex1.2 <- apply(flexpred1.2, FUN="mean", MARGIN=2)
  
  ##4 noeuds
  
  flexpred1.4 <- predict(flex.model1.4, newtimes = newtimes)$predictions
  
  mean.flex1.4 <- apply(flexpred1.4, FUN="mean", MARGIN=2)
  
  #### modèle NPH
  
  ##2 noeuds
  
  flexpred2.2 <- predict(flex.model2.2 , newtimes = newtimes)$predictions
  
  mean.flex2.2  <- apply(flexpred2.2 , FUN="mean", MARGIN=2)
  
  ##4 noeuds
  
  flexpred2.4 <- predict(flex.model2.4 , newtimes = newtimes)$predictions
  
  mean.flex2.4  <- apply(flexpred2.4 , FUN="mean", MARGIN=2)
  ### genweibull
  
  WGpred <- predict(WG.model, newtimes = newtimes)$predictions
  
  mean.WG <- apply(WGpred, FUN="mean", MARGIN=2)
  
  
  ### Estimateur de Pohar-Perme 
  
  PP <- summary( rs.surv(Surv(times, status) ~ 1, data = data_train,
                         ratetable = slopop, method = "pohar-perme", add.times = newtimes,
                         rmap = list(age = age, sex = sexchara, year = year)), times = newtimes, extend = TRUE) 
  
  poharpred <- PP$surv 
  
  ###### par strates ######
  ######
  #plann
  
  for(j in strata_names){
    
    data_train_var <- get(paste0("data_train_", j))
    
    plannpredS <- predictRS(plann.model, data = data_train_var, newtimes = newtimes, 
                            ratetable = slopop, age = "age", year = "year", sex = "sexchara")
    
    mean.plannS <- plannpredS$mpredictions$net_survival
    
    assign(paste0("plannpred_", j), plannpredS)
    assign(paste0("mean.plann_", j), mean.plannS)
  }
  
  ###flexnet
  
  ### modèle PH
  
  ## 2 noeuds
  
  for(j in strata_names){
    
    data_train_var <- get(paste0("data_train_", j))
    
    flexpred1.2S <- predict(flex.model1.2, newtimes = newtimes, newdata = data_train_var)$predictions
    
    mean.flex1.2S <- apply(flexpred1.2S, FUN="mean", MARGIN=2)
    
    assign(paste0("flexpred1.2_", j), flexpred1.2S)
    assign(paste0("mean.flex1.2_", j), mean.flex1.2S)
  }
  
  ## 4 noeuds
  
  for(j in strata_names){
    
    data_train_var <- get(paste0("data_train_", j))
    
    flexpred1.4S <- predict(flex.model1.4, newtimes = newtimes, newdata = data_train_var)$predictions
    
    mean.flex1.4S <- apply(flexpred1.4S, FUN="mean", MARGIN=2)
    
    assign(paste0("flexpred1.4_", j), flexpred1.4S)
    assign(paste0("mean.flex1.4_", j), mean.flex1.4S)
  }
  
  #### modèle NPH
  
  ##2 noeuds
  
  for(j in strata_names){
    
    data_train_var <- get(paste0("data_train_", j))
    
    flexpred2.2S <- predict(flex.model2.2, newtimes = newtimes, newdata = data_train_var)$predictions
    
    mean.flex2.2S <- apply(flexpred2.2S, FUN="mean", MARGIN=2)
    
    assign(paste0("flexpred2.2_", j), flexpred2.2S)
    assign(paste0("mean.flex2.2_", j), mean.flex2.2S)
  }
  
  ##4 noeuds
  
  for(j in strata_names){
    
    data_train_var <- get(paste0("data_train_", j))
    
    flexpred2.4S <- predict(flex.model2.4, newtimes = newtimes, newdata = data_train_var)$predictions
    
    mean.flex2.4S <- apply(flexpred2.4S, FUN="mean", MARGIN=2)
    
    assign(paste0("flexpred2.4_", j), flexpred2.4S)
    assign(paste0("mean.flex2.4_", j), mean.flex2.4S)
  }
  
  ### genweibull
  
  for(j in strata_names){
    
    data_train_var <- get(paste0("data_train_", j))
    
    WGpredS <- predict(WG.model, newtimes = newtimes, newdata = data_train_var)$predictions
    
    mean.WGS <- apply(WGpredS, FUN="mean", MARGIN=2)
    
    assign(paste0("WGpred_", j), WGpredS)
    assign(paste0("mean.WG_", j), mean.WGS)
  }
  
  
  ##Estimateur de Pohar-Perme 
  
  for(j in strata_names){
    
    data_train_var <- get(paste0("data_train_", j))
    
    PPS <- summary( rs.surv(Surv(times, status) ~ 1, data = data_train_var,
                            ratetable = slopop, method = "pohar-perme", add.times = newtimes,
                            rmap = list(age = age, sex = sexchara, year = year)), times = newtimes, extend = TRUE) 
    
    poharpredS <- PPS$surv
    
    assign(paste0("PP_", j), PPS)
    assign(paste0("poharpred_", j), poharpredS)
  }
  
  ###### 
  
  ############################################
  ########### BASE DE VALIDATION
  ###plann
  
  plannpredval <-  predictRS(plann.model, data = data_valid, newtimes = newtimes, ratetable=slopop, 
                             age = "age", year = "year", sex = "sexchara")
  
  mean.plannval <- plannpredval$mpredictions$net_survival 
  
  ###flexnet
  
  ### modèle PH
  
  ### 2 noeuds
  
  flexpredval1.2 <- predict(flex.model1.2, newtimes = newtimes, newdata = data_valid)$predictions
  
  mean.flexval1.2 <- apply(flexpredval1.2, FUN="mean", MARGIN=2)
  
  ### 4 noeuds
  
  flexpredval1.4 <- predict(flex.model1.4, newtimes = newtimes, newdata = data_valid)$predictions
  
  mean.flexval1.4 <- apply(flexpredval1.4, FUN="mean", MARGIN=2)
  
  #### modèle NPH
  
  ### 2 noeuds
  
  flexpredval2.2 <- predict(flex.model2.2, newtimes = newtimes, newdata = data_valid)$predictions
  
  mean.flexval2.2 <- apply(flexpredval2.2, FUN="mean", MARGIN=2)
  
  ### 4 noeuds
  
  flexpredval2.4 <- predict(flex.model2.4, newtimes = newtimes, newdata = data_valid)$predictions
  
  mean.flexval2.4 <- apply(flexpredval2.4, FUN="mean", MARGIN=2)
  
  ### genweibull
  
  WGpredval <- predict(WG.model, newtimes = newtimes, newdata = data_valid)$predictions
  
  mean.WGval <- apply(WGpredval, FUN="mean", MARGIN=2)
  
  ### Estimateur de Pohar-Perme 
  
  PPval <- summary( rs.surv(Surv(times, status) ~ 1, data = data_valid,
                            ratetable = slopop, method = "pohar-perme", add.times = newtimes,
                            rmap = list(age = age, sex = sexchara, year = year)), times = newtimes, extend = TRUE) 
  
  poharpredval <- PPval$surv 
  
  ###### par strates ######
  
  
  ###plann
  
  for(j in strata_names){
    
    data_valid_var <- get(paste0("data_valid_", j))
    
    plannpredvalS <- predictRS(plann.model, data = data_valid_var, newtimes = newtimes, 
                               ratetable = slopop, age = "age", year = "year", sex = "sexchara")
    
    mean.plannvalS <- plannpredvalS$mpredictions$net_survival
    
    assign(paste0("plannpredval_", j), plannpredvalS)
    assign(paste0("mean.plannval_", j), mean.plannvalS)
  }
  
  ###flexnet
  
  ### modèle PH
  
  ### 2 noeuds 
  
  for(j in strata_names){
    
    data_valid_var <- get(paste0("data_valid_", j))
    
    flexpredval1.2S <- predict(flex.model1.2, newtimes = newtimes, newdata = data_valid_var)$predictions
    
    mean.flexval1.2S <- apply(flexpredval1.2S, FUN="mean", MARGIN=2)
    
    assign(paste0("flexpredval1.2_", j), flexpredval1.2S)
    assign(paste0("mean.flexval1.2_", j), mean.flexval1.2S)
  }
  
  ### 4 noeuds 
  
  for(j in strata_names){
    
    data_valid_var <- get(paste0("data_valid_", j))
    
    flexpredval1.4S <- predict(flex.model1.4, newtimes = newtimes, newdata = data_valid_var)$predictions
    
    mean.flexval1.4S <- apply(flexpredval1.4S, FUN="mean", MARGIN=2)
    
    assign(paste0("flexpredval1.4_", j), flexpredval1.4S)
    assign(paste0("mean.flexval1.4_", j), mean.flexval1.4S)
  }
  #### modèle NPH
  
  ## 2 noeuds
  
  for(j in strata_names){
    
    data_valid_var <- get(paste0("data_valid_", j))
    
    flexpredval2.2S <- predict(flex.model2.2, newtimes = newtimes, newdata = data_valid_var)$predictions
    
    mean.flexval2.2S <- apply(flexpredval2.2S, FUN="mean", MARGIN=2)
    
    assign(paste0("flexpredval2.2_", j), flexpredval2.2S)
    assign(paste0("mean.flexval2.2_", j), mean.flexval2.2S)
  }
  
  ## 4 noeuds
  
  for(j in strata_names){
    
    data_valid_var <- get(paste0("data_valid_", j))
    
    flexpredval2.4S <- predict(flex.model2.4, newtimes = newtimes, newdata = data_valid_var)$predictions
    
    mean.flexval2.4S <- apply(flexpredval2.4S, FUN="mean", MARGIN=2)
    
    assign(paste0("flexpredval2.4_", j), flexpredval2.4S)
    assign(paste0("mean.flexval2.4_", j), mean.flexval2.4S)
  }
  
  ### genweibull
  
  for(j in strata_names){
    
    data_valid_var <- get(paste0("data_valid_", j))
    
    WGpredvalS <- predict(WG.model, newtimes = newtimes, newdata = data_valid_var)$predictions
    
    mean.WGvalS <- apply(WGpredvalS, FUN="mean", MARGIN=2)
    
    assign(paste0("WGpredval_", j), WGpredvalS)
    assign(paste0("mean.WGval_", j), mean.WGvalS)
  }
  
  ##Estimateur de Pohar-Perme 
  
  for(j in strata_names){
    
    data_valid_var <- get(paste0("data_valid_", j))
    
    PPvalS <- summary( rs.surv(Surv(times, status) ~ 1, data = data_valid_var,
                               ratetable = slopop, method = "pohar-perme", add.times = newtimes,
                               rmap = list(age = age, sex = sexchara, year = year)), times = newtimes, extend = TRUE) 
    
    poharpredvalS <- PPvalS$surv
    
    assign(paste0("PPval_", j), PPvalS)
    assign(paste0("poharpredval_", j), poharpredvalS)
  }
  ##################################################################################################
  
  ######
  #time
  ######
  
  time <- Sys.time()-start
  time <- as.numeric(time, units = 'mins')
  ########################
  #       outputs        #
  ########################
  
  #####bases####
  
  ### TRAIN
  write.table(data_train,  paste0(path0, "DATAFRAMES/TRAIN/WHOLE/",i,"_dftrain.csv"), sep = ";", row.names = F, col.names = T)
  ##strates
  for(j in strata_names){
    
    data_train_var <- get(paste0("data_train_", j))
    
    file_path <- paste0(path0, "DATAFRAMES/TRAIN/", j, "/", i, "_dftrain_", j, ".csv")
    
    write.table(data_train_var, file = file_path, sep = ";", row.names = FALSE, col.names = TRUE)
  }
  ### VALID
  write.table(data_valid,  paste0(path0, "DATAFRAMES/VALID/WHOLE/",i,"_dfvalid.csv"), sep = ";", row.names = F, col.names = T)
  ##strates
  
  for(j in strata_names){
    
    data_valid_var <- get(paste0("data_valid_", j))
    
    file_path <- paste0(path0, "DATAFRAMES/VALID/", j, "/", i, "_dfvalid_", j, ".csv")
    
    write.table(data_valid_var, file = file_path, sep = ";", row.names = FALSE, col.names = TRUE)
  }
  
  ############## paramètres estimés par CV ###########
  
  m1.2value <- paste0(m1.2, failedindic1.2)
  m1.4value <- paste0(m1.4, failedindic1.4)
  
  m2.2value <- paste0(m2.2, failedindic2.2)
  m2.4value <- paste0(m2.4, failedindic2.4)
  
  paramsCV <- as.data.frame(c(inter = inter,size = size,decay = decay, maxit = maxit,
                              MaxNWts = MaxNWts, m1.2 = m1.2value, m1.4 = m1.4value, 
                              m2.2 = m2.2value, m2.4 = m2.4value, logcoeffWG = logcoeffWG))
  
  write.table(paramsCV,  paste0(path0, "PARAM/",i,"_param.csv"))
  
  
  ############# temps écoulé total & par CV de chaque modèle et taille de la base #########################
  
  timeComputation <- as.data.frame(c(time = time, plann = TimeCVplann, N = dim(data)[1]))
  
  write.table(timeComputation,  paste0(path0, "TIME/",i,"_time.csv"))
  
  ##############Log(LIKELIHOOD)####################### 
  
  plann.loglik <- plannpred$loglik  
  
  flex1.2.loglik <- unname(flex.model1.2$loglik[1])
  flex1.4.loglik <- unname(flex.model1.4$loglik[1])
  
  flex2.2.loglik <- unname(flex.model2.2$loglik[1])
  flex2.4.loglik <- unname(flex.model2.4$loglik[1])
  
  WG.loglik <- unname(WG.model$loglik[1])
  
  logliks <- as.data.frame(c(plann = plann.loglik, flex1.2 = flex1.2.loglik,
                             flex1.4 = flex1.4.loglik, flex2.2 = flex2.2.loglik,
                             flex2.4 = flex2.4.loglik, WG = WG.loglik))
  
  write.table(logliks,  paste0(path0, "LOGLIK/",i,"_loglik.csv"))
  
  ##################DATA TRAIN
  #################################### MOYENNES
  #######THEO
  
  mean_theo <- t(data.frame(THEO = mean.theo))
  write.table(mean_theo, paste0(path0, "TRAIN/WHOLE/THEO/mean_survival/",i,"_mean_theo.csv"), sep = ";", row.names = F, col.names = T)
  
  
  ######PLANN
  
  mean_PLANN <- t(data.frame(PLANN = mean.plann))
  write.table(mean_PLANN, paste0(path0, "TRAIN/WHOLE/PLANN/mean_survival/",i,"_mean_PLANN.csv"), sep = ";", row.names = F, col.names = T)
  
  #######FLEX1
  
  ##2 noeuds
  mean_flex1.2 <- t(data.frame(flex1.2=mean.flex1.2))
  write.table(mean_flex1.2, paste0(path0, "TRAIN/WHOLE/FLEX1.2/mean_survival/",i,"_mean_flex12.csv"), sep = ";", row.names = F, col.names = T)
  
  ##4 noeuds
  mean_flex1.4 <- t(data.frame(flex1.4=mean.flex1.4))
  write.table(mean_flex1.4, paste0(path0, "TRAIN/WHOLE/FLEX1.4/mean_survival/",i,"_mean_flex14.csv"), sep = ";", row.names = F, col.names = T)
  
  ######FLEX2
  
  ##2 noeuds
  mean_flex2.2 <- t(data.frame(flex2.2=mean.flex2.2))
  write.table(mean_flex2.2, paste0(path0, "TRAIN/WHOLE/FLEX2.2/mean_survival/",i,"_mean_flex22.csv"), sep = ";", row.names = F, col.names = T)
  
  ##4 noeuds
  mean_flex2.4 <- t(data.frame(flex2.4=mean.flex2.4))
  write.table(mean_flex2.4, paste0(path0, "TRAIN/WHOLE/FLEX2.4/mean_survival/",i,"_mean_flex24.csv"), sep = ";", row.names = F, col.names = T)
  
  ######WG
  
  mean_flexWG <- t(data.frame(flexWG=mean.WG))
  write.table(mean_flexWG, paste0(path0, "TRAIN/WHOLE/WG/mean_survival/",i,"_mean_flexWG.csv"), sep = ";",row.names = F, col.names = T)
  
  ######Estimateur de Pohar-Perme 
  
  write.table(poharpred, paste0(path0, "TRAIN/WHOLE/PP/",i,"PP.csv"), sep = ";",row.names = T, col.names = T)
  ###strates 
  ######
  
  ##theo
  
  mean.theo_STRATES = c()
  
  for (j in strata_names) {
    
    mean.theo_STRATES <- rbind(mean.theo_STRATES, get(paste0("mean.theo_", j) ) )
    rownames(mean.theo_STRATES)[nrow(mean.theo_STRATES)] <- paste0("mean.theo_", j)
  }
  
  mean.theo_STRATES <- data.frame(strat = rownames(mean.theo_STRATES), mean.theo_STRATES)
  write.table(mean.theo_STRATES, paste0(path0, "TRAIN/STRATA/THEO/mean_survival/",i,"mean_theo_strates.csv"), sep = ";",row.names = F, col.names = T)
  
  #PLANN
  
  mean.plann_STRATES = c()
  
  for (j in strata_names) {
    
    mean.plann_STRATES <- rbind(mean.plann_STRATES, get(paste0("mean.plann_", j) ) )
    rownames(mean.plann_STRATES)[nrow(mean.plann_STRATES)] <- paste0("mean.plann_", j)
  }
  
  mean.plann_STRATES <- data.frame(strat = rownames(mean.plann_STRATES), mean.plann_STRATES)
  write.table(mean.plann_STRATES, paste0(path0, "TRAIN/STRATA/PLANN/mean_survival/",i,"mean_plann_strates.csv"), sep = ";",row.names = F, col.names = T)
  
  
  #flex1
  
  ##2 noeuds
  mean.flex1.2_STRATES = c()
  
  for (j in strata_names) {
    
    mean.flex1.2_STRATES <- rbind(mean.flex1.2_STRATES, get(paste0("mean.flex1.2_", j) ) )
    rownames(mean.flex1.2_STRATES)[nrow(mean.flex1.2_STRATES)] <- paste0("mean.flex1.2_", j)
  }
  
  mean.flex1.2_STRATES <- data.frame(strat = rownames(mean.flex1.2_STRATES), mean.flex1.2_STRATES)
  write.table(mean.flex1.2_STRATES, paste0(path0, "TRAIN/STRATA/FLEX1.2/mean_survival/",i,"mean_flex12_strates.csv"), sep = ";",row.names = F, col.names = T)
  
  ##4 noeuds
  
  mean.flex1.4_STRATES = c()
  
  for (j in strata_names) {
    
    mean.flex1.4_STRATES <- rbind(mean.flex1.4_STRATES, get(paste0("mean.flex1.4_", j) ) )
    rownames(mean.flex1.4_STRATES)[nrow(mean.flex1.4_STRATES)] <- paste0("mean.flex1.4_", j)
  }
  
  mean.flex1.4_STRATES <- data.frame(strat = rownames(mean.flex1.4_STRATES), mean.flex1.4_STRATES)
  write.table(mean.flex1.4_STRATES, paste0(path0, "TRAIN/STRATA/FLEX1.4/mean_survival/",i,"mean_flex14_strates.csv"), sep = ";",row.names = F, col.names = T)
  
  #flex2
  
  ## 2 noeuds
  
  mean.flex2.2_STRATES = c()
  
  for (j in strata_names) {
    
    mean.flex2.2_STRATES <- rbind(mean.flex2.2_STRATES, get(paste0("mean.flex2.2_", j) ) )
    rownames(mean.flex2.2_STRATES)[nrow(mean.flex2.2_STRATES)] <- paste0("mean.flex2.2_", j)
  }
  
  mean.flex2.2_STRATES <- data.frame(strat = rownames(mean.flex2.2_STRATES), mean.flex2.2_STRATES)
  write.table(mean.flex2.2_STRATES, paste0(path0, "TRAIN/STRATA/FLEX2.2/mean_survival/",i,"mean_flex22_strates.csv"), sep = ";",row.names = F, col.names = T)
  
  ## 4 noeuds
  
  mean.flex2.4_STRATES = c()
  
  for (j in strata_names) {
    
    mean.flex2.4_STRATES <- rbind(mean.flex2.4_STRATES, get(paste0("mean.flex2.4_", j) ) )
    rownames(mean.flex2.4_STRATES)[nrow(mean.flex2.4_STRATES)] <- paste0("mean.flex2.4_", j)
  }
  
  mean.flex2.4_STRATES <- data.frame(strat = rownames(mean.flex2.4_STRATES), mean.flex2.4_STRATES)
  write.table(mean.flex2.4_STRATES, paste0(path0, "TRAIN/STRATA/FLEX2.4/mean_survival/",i,"mean_flex24_strates.csv"), sep = ";",row.names = F, col.names = T)
  
  #WG
  mean.WG_STRATES = c()
  
  for (j in strata_names) {
    
    mean.WG_STRATES <- rbind(mean.WG_STRATES, get(paste0("mean.WG_", j) ) )
    rownames(mean.WG_STRATES)[nrow(mean.WG_STRATES)] <- paste0("mean.WG_", j)
  }
  
  mean.WG_STRATES <- data.frame(strat = rownames(mean.WG_STRATES), mean.WG_STRATES)
  write.table(mean.WG_STRATES, paste0(path0, "TRAIN/STRATA/WG/mean_survival/",i,"mean_WG_strates.csv"), sep = ";",row.names = F, col.names = T)
  
  
  ####Estimateur de Pohar-Perme
  
  poharpred_strates = c()
  
  for (j in strata_names) {
    
    poharpred_strates <- rbind(poharpred_strates, get(paste0("poharpred_", j) ) )
    rownames(poharpred_strates)[nrow(poharpred_strates)] <- paste0("poharpred_", j)
  }
  
  poharpred_strates <- data.frame(strat = rownames(poharpred_strates), poharpred_strates)
  write.table(poharpred_strates, paste0(path0, "TRAIN/STRATA/PP/",i,"PP_strates.csv"), sep = ";",row.names = F, col.names = T)
  
  ######
  
  #################################### SURVIES INDIVIDUELLES
  
  #######THEO
  
  ind_theo <- data.frame(theopred)
  write.table(ind_theo, paste0(path0, "TRAIN/WHOLE/THEO/individual_survival/",
                               i,"_ind_theo.csv"), sep = ";", row.names = F, col.names = T)
  
  ######PLANN
  
  ind_PLANN <- data.frame(plannpred$ipredictions$relative_survival)
  write.table(ind_PLANN, paste0(path0, "TRAIN/WHOLE/PLANN/individual_survival/",
                                i,"_ind_PLANN.csv"), sep = ";", row.names = F, col.names = T)
  
  #######FLEX1
  
  ##2 noeuds
  ind_flex1.2 <- data.frame(flexpred1.2)
  write.table(ind_flex1.2, paste0(path0, "TRAIN/WHOLE/FLEX1.2/individual_survival/"
                                  ,i,"_ind_flex12.csv"), sep = ";", row.names = F, col.names = T)
  
  ##4 noeuds
  ind_flex1.4 <- data.frame(flexpred1.4)
  write.table(ind_flex1.4, paste0(path0, "TRAIN/WHOLE/FLEX1.4/individual_survival/"
                                  ,i,"_ind_flex14.csv"), sep = ";", row.names = F, col.names = T)
  ######FLEX2
  
  ##2 noeuds
  ind_flex2.2 <- data.frame(flexpred2.2)
  write.table(ind_flex2.2, paste0(path0, "TRAIN/WHOLE/FLEX2.2/individual_survival/"
                                  ,i,"_ind_flex22.csv"), sep = ";", row.names = F, col.names = T)
  
  ##4 noeuds
  ind_flex2.4 <- data.frame(flexpred2.4)
  write.table(ind_flex2.4, paste0(path0, "TRAIN/WHOLE/FLEX2.4/individual_survival/"
                                  ,i,"_ind_flex24.csv"), sep = ";", row.names = F, col.names = T)
  
  ######WG
  
  ind_WG <- data.frame(WGpred)
  write.table(ind_WG, paste0(path0, "TRAIN/WHOLE/WG/individual_survival/"
                             ,i,"_ind_WG.csv"), sep = ";",row.names = F, col.names = T)
  
  
  ###strates 
  ######
  
  #theo
  
  for (j in strata_names) {
    
    ind_theo <- data.frame(get(paste0("theopred_", j)))
    
    write.table(ind_theo, paste0(path0, "TRAIN/STRATA/THEO/individual_survival/", i, "_ind_theo_", j, ".csv"),
                sep = ";",row.names = F, col.names = T)
  }
  
  ##plann
  
  for (j in strata_names) {
    
    ind_PLANN <- data.frame(get(paste0("plannpred_", j))$ipredictions$relative_survival)
    
    write.table(ind_PLANN, paste0(path0, "TRAIN/STRATA/PLANN/individual_survival/", i, "_ind_PLANN_", j, ".csv"),
                sep = ";",row.names = F, col.names = T)
  }
  
  #flex1
  
  ## 2 noeuds
  for (j in strata_names) {
    
    ind_flex1.2 <- data.frame(get(paste0("flexpred1.2_", j)))
    
    write.table(ind_flex1.2, paste0(path0, "TRAIN/STRATA/FLEX1.2/individual_survival/", i, "_ind_flex12_", j, ".csv"),
                sep = ";",row.names = F, col.names = T)
  }
  
  ## 4 noeuds
  for (j in strata_names) {
    
    ind_flex1.4 <- data.frame(get(paste0("flexpred1.4_", j)))
    
    write.table(ind_flex1.4, paste0(path0, "TRAIN/STRATA/FLEX1.4/individual_survival/", i, "_ind_flex14_", j, ".csv"),
                sep = ";",row.names = F, col.names = T)
  }
  
  #flex2
  
  ##2 noeuds
  for (j in strata_names) {
    
    ind_flex2.2 <- data.frame(get(paste0("flexpred2.2_", j)))
    
    write.table(ind_flex2.2, paste0(path0, "TRAIN/STRATA/FLEX2.2/individual_survival/", i, "_ind_flex22_", j, ".csv"),
                sep = ";",row.names = F, col.names = T)
  }
  
  ##4 noeuds
  for (j in strata_names) {
    
    ind_flex2.4 <- data.frame(get(paste0("flexpred2.4_", j)))
    
    write.table(ind_flex2.4, paste0(path0, "TRAIN/STRATA/FLEX2.4/individual_survival/", i, "_ind_flex24_", j, ".csv"),
                sep = ";",row.names = F, col.names = T)
  }
  #WG
  
  for (j in strata_names) {
    
    ind_WG <- data.frame(get(paste0("WGpred_", j)))
    
    write.table(ind_WG, paste0(path0, "TRAIN/STRATA/WG/individual_survival/", i, "_ind_WG_", j, ".csv"),
                sep = ";",row.names = F, col.names = T)
  }
  
  ######
  
  ##################DATA VALID
  #################################### MOYENNES
  
  #######THEO
  
  mean_theoval <- t(data.frame(THEO = mean.theoval))
  write.table(mean_theoval, paste0(path0, "VALID/WHOLE/THEO/mean_survival/",
                                   i,"_mean_theoval.csv"), sep = ";", row.names = F, col.names = T)
  
  ######PLANN
  
  mean_PLANNval <- t(data.frame(PLANN = mean.plannval))
  write.table(mean_PLANNval, paste0(path0, "VALID/WHOLE/PLANN/mean_survival/",
                                    i,"_mean_PLANNval.csv"), sep = ";", row.names = F, col.names = T)
  
  #######FLEX1
  
  ##2 noeuds
  mean_flexval1.2 <- t(data.frame(flex1.2=mean.flexval1.2))
  write.table(mean_flexval1.2, paste0(path0, "VALID/WHOLE/FLEX1.2/mean_survival/"
                                      ,i,"_mean_flexval12.csv"), sep = ";", row.names = F, col.names = T)
  ##4 noeuds
  mean_flexval1.4 <- t(data.frame(flex1.4=mean.flexval1.4))
  write.table(mean_flexval1.4, paste0(path0, "VALID/WHOLE/FLEX1.4/mean_survival/"
                                      ,i,"_mean_flexval14.csv"), sep = ";", row.names = F, col.names = T)
  
  ######FLEX2
  
  ## 2 noeuds
  mean_flexval2.2 <- t(data.frame(flex2.2=mean.flexval2.2))
  write.table(mean_flexval2.2, paste0(path0, "VALID/WHOLE/FLEX2.2/mean_survival/"
                                      ,i,"_mean_flexval22.csv"), sep = ";", row.names = F, col.names = T)
  
  ## 4 noeuds
  mean_flexval2.4 <- t(data.frame(flex2.4=mean.flexval2.4))
  write.table(mean_flexval2.4, paste0(path0, "VALID/WHOLE/FLEX2.4/mean_survival/"
                                      ,i,"_mean_flexval24.csv"), sep = ";", row.names = F, col.names = T)
  ######WG
  
  mean_flexWGval <- t(data.frame(flexWG=mean.WGval))
  write.table(mean_flexWGval, paste0(path0, "VALID/WHOLE/WG/mean_survival/"
                                     ,i,"_mean_flexWGval.csv"), sep = ";",row.names = F, col.names = T)
  
  ###### Estimateur de pohar perme
  
  write.table(poharpred, paste0(path0, "VALID/WHOLE/PP/",i,"PPval.csv"), sep = ";",row.names = T, col.names = T)
  
  ###strates 
  ######
  
  ##theo
  
  mean.theoval_STRATES <- c()
  
  for (j in strata_names) {
    
    mean.theoval_STRATES <- rbind(mean.theoval_STRATES, get(paste0("mean.theoval_", j) ) )
    rownames(mean.theoval_STRATES)[nrow(mean.theoval_STRATES)] <- paste0("mean.theoval_", j)
  }
  
  mean.theoval_STRATES <- data.frame(strat = rownames(mean.theo_STRATES), mean.theo_STRATES)
  write.table(mean.theoval_STRATES, paste0(path0, "VALID/STRATA/THEO/mean_survival/",i,"mean_theoval_strates.csv"), sep = ";",row.names = F, col.names = T)
  
  
  #PLANN
  
  mean.plannval_STRATES = c()
  
  for (j in strata_names) {
    
    mean.plannval_STRATES <- rbind(mean.plannval_STRATES, get(paste0("mean.plannval_", j) ) )
    rownames(mean.plannval_STRATES)[nrow(mean.plannval_STRATES)] <- paste0("mean.plannval_", j)
  }
  
  mean.plannval_STRATES <- data.frame(strat = rownames(mean.plannval_STRATES), mean.plannval_STRATES)
  write.table(mean.plannval_STRATES, paste0(path0, "VALID/STRATA/PLANN/mean_survival/",i,"mean_plannval_strates.csv"), sep = ";",row.names = F, col.names = T)
  
  #flex1
  ## 2 noeuds
  mean.flexval1.2_STRATES = c()
  
  for (j in strata_names) {
    
    mean.flexval1.2_STRATES <- rbind(mean.flexval1.2_STRATES, get(paste0("mean.flexval1.2_", j) ) )
    rownames(mean.flexval1.2_STRATES)[nrow(mean.flexval1.2_STRATES)] <- paste0("mean.flexval1.2_", j)
  }
  
  mean.flexval1.2_STRATES <- data.frame(strat = rownames(mean.flexval1.2_STRATES), mean.flexval1.2_STRATES)
  write.table(mean.flexval1.2_STRATES, paste0(path0, "VALID/STRATA/FLEX1.2/mean_survival/",i,"mean_flexval12_strates.csv"), sep = ";",row.names = F, col.names = T)
  
  ##4 noeuds
  mean.flexval1.4_STRATES = c()
  
  for (j in strata_names) {
    
    mean.flexval1.4_STRATES <- rbind(mean.flexval1.4_STRATES, get(paste0("mean.flexval1.4_", j) ) )
    rownames(mean.flexval1.4_STRATES)[nrow(mean.flexval1.4_STRATES)] <- paste0("mean.flexval1.4_", j)
  }
  
  mean.flexval1.4_STRATES <- data.frame(strat = rownames(mean.flexval1.4_STRATES), mean.flexval1.4_STRATES)
  write.table(mean.flexval1.4_STRATES, paste0(path0, "VALID/STRATA/FLEX1.4/mean_survival/",i,"mean_flexval14_strates.csv"), sep = ";",row.names = F, col.names = T)
  
  
  #flex2
  
  ##2 noeuds
  mean.flexval2.2_STRATES = c()
  
  for (j in strata_names) {
    
    mean.flexval2.2_STRATES <- rbind(mean.flexval2.2_STRATES, get(paste0("mean.flexval2.2_", j) ) )
    rownames(mean.flexval2.2_STRATES)[nrow(mean.flexval2.2_STRATES)] <- paste0("mean.flexval2.2_", j)
  }
  
  mean.flexval2.2_STRATES <- data.frame(strat = rownames(mean.flexval2.2_STRATES), mean.flexval2.2_STRATES)
  write.table(mean.flexval2.2_STRATES, paste0(path0, "VALID/STRATA/FLEX2.2/mean_survival/",i,"mean_flexval22_strates.csv"), sep = ";",row.names = F, col.names = T)
  
  ##4 noeuds
  mean.flexval2.4_STRATES = c()
  
  for (j in strata_names) {
    
    mean.flexval2.4_STRATES <- rbind(mean.flexval2.4_STRATES, get(paste0("mean.flexval2.4_", j) ) )
    rownames(mean.flexval2.4_STRATES)[nrow(mean.flexval2.4_STRATES)] <- paste0("mean.flexval2.4_", j)
  }
  
  mean.flexval2.4_STRATES <- data.frame(strat = rownames(mean.flexval2.4_STRATES), mean.flexval2.4_STRATES)
  write.table(mean.flexval2.4_STRATES, paste0(path0, "VALID/STRATA/FLEX2.4/mean_survival/",i,"mean_flexval24_strates.csv"), sep = ";",row.names = F, col.names = T)
  
  #WG
  mean.WGval_STRATES = c()
  
  for (j in strata_names) {
    
    mean.WGval_STRATES <- rbind(mean.WGval_STRATES, get(paste0("mean.WGval_", j) ) )
    rownames(mean.WGval_STRATES)[nrow(mean.WGval_STRATES)] <- paste0("mean.WGval_", j)
  }
  
  mean.WGval_STRATES <- data.frame(strat = rownames(mean.WGval_STRATES), mean.WGval_STRATES)
  write.table(mean.WGval_STRATES, paste0(path0, "VALID/STRATA/WG/mean_survival/",i,"mean_WGval_strates.csv"), sep = ";",row.names = F, col.names = T)
  
  ####Estimateur de Pohar-Perme
  
  poharpredval_strates = c()
  
  for (j in strata_names) {
    
    poharpredval_strates <- rbind(poharpredval_strates, get(paste0("poharpredval_", j) ) )
    rownames(poharpredval_strates)[nrow(poharpredval_strates)] <- paste0("poharpredval_", j)
  }
  
  poharpredval_strates <- data.frame(strat = rownames(poharpredval_strates), poharpredval_strates)
  write.table(poharpredval_strates, paste0(path0, "VALID/STRATA/PP/",i,"PPval_strates.csv"), sep = ";",row.names = F, col.names = T)
  
  ######
  #################################### SURVIES INDIVIDUELLES
  
  #######THEO
  ind_theoval <- data.frame(theopredval)
  write.table(ind_theoval, 
              paste0(path0, "VALID/WHOLE/THEO/individual_survival/",i
                     ,"_ind_theoval.csv"), sep = ";", row.names = F, col.names = T)
  
  ######PLANN
  
  ind_PLANNval <- data.frame(plannpredval$ipredictions$relative_survival)
  write.table(ind_PLANNval, paste0(path0, "VALID/WHOLE/PLANN/individual_survival/",
                                   i,"_ind_PLANNval.csv"), sep = ";", row.names = F, col.names = T)
  
  #######FLEX1
  
  ##2 noeuds
  
  ind_flexval1.2 <- data.frame(flexpredval1.2)
  write.table(ind_flexval1.2, paste0(path0, "VALID/WHOLE/FLEX1.2/individual_survival/"
                                     ,i,"_ind_flexval12.csv"), sep = ";", row.names = F, col.names = T)
  
  ##4 noeuds
  
  ind_flexval1.4 <- data.frame(flexpredval1.4)
  write.table(ind_flexval1.4, paste0(path0, "VALID/WHOLE/FLEX1.4/individual_survival/"
                                     ,i,"_ind_flexval14.csv"), sep = ";", row.names = F, col.names = T)
  ######FLEX2
  
  ##2 noeuds
  ind_flexval2.2 <- data.frame(flexpredval2.2)
  write.table(ind_flexval2.2, paste0(path0, "VALID/WHOLE/FLEX2.2/individual_survival/"
                                     ,i,"_ind_flexval22.csv"), sep = ";", row.names = F, col.names = T)
  
  ##4 noeuds
  ind_flexval2.4 <- data.frame(flexpredval2.4)
  write.table(ind_flexval2.4, paste0(path0, "VALID/WHOLE/FLEX2.4/individual_survival/"
                                     ,i,"_ind_flexval24.csv"), sep = ";", row.names = F, col.names = T)
  ######WG
  
  ind_flexWGval <- data.frame(WGpredval)
  write.table(ind_flexWGval, paste0(path0, "VALID/WHOLE/WG/individual_survival/"
                                    ,i,"_ind_flexWGval.csv"), sep = ";",row.names = F, col.names = T)
  
  
  ###strates ####
  ######
  
  #theo
  
  for (j in strata_names) {
    
    ind_theoval <- data.frame(get(paste0("theopredval_", j)))
    
    write.table(ind_theoval, paste0(path0, "VALID/STRATA/THEO/individual_survival/", i, "_ind_theoval_", j, ".csv"),
                sep = ";",row.names = F, col.names = T)
  }
  
  ##plann
  
  for (j in strata_names) {
    
    ind_PLANNval <- data.frame(get(paste0("plannpredval_", j))$ipredictions$relative_survival)
    
    write.table(ind_PLANNval, paste0(path0, "VALID/STRATA/PLANN/individual_survival/", i, "_ind_PLANNval_", j, ".csv"),
                sep = ";",row.names = F, col.names = T)
  }
  
  #flex1
  
  ##2 noeuds
  for (j in strata_names) {
    
    ind_flexval1.2 <- data.frame(get(paste0("flexpredval1.2_", j)))
    
    write.table(ind_flexval1.2, paste0(path0, "VALID/STRATA/FLEX1.2/individual_survival/", i, "_ind_flexval12_", j, ".csv"),
                sep = ";",row.names = F, col.names = T)
  }
  ##4 noeuds
  for (j in strata_names) {
    
    ind_flexval1.4 <- data.frame(get(paste0("flexpredval1.4_", j)))
    
    write.table(ind_flexval1.4, paste0(path0, "VALID/STRATA/FLEX1.4/individual_survival/", i, "_ind_flexval14_", j, ".csv"),
                sep = ";",row.names = F, col.names = T)
  }
  #flex2
  
  ##2 noeuds
  for (j in strata_names) {
    
    ind_flexval2.2 <- data.frame(get(paste0("flexpredval2.2_", j)))
    
    write.table(ind_flexval2.2, paste0(path0, "VALID/STRATA/FLEX2.2/individual_survival/", i, "_ind_flexval22_", j, ".csv"),
                sep = ";",row.names = F, col.names = T)
  }
  
  ##4 noeuds
  for (j in strata_names) {
    
    ind_flexval2.4 <- data.frame(get(paste0("flexpredval2.4_", j)))
    
    write.table(ind_flexval2.4, paste0(path0, "VALID/STRATA/FLEX2.4/individual_survival/", i, "_ind_flexval24_", j, ".csv"),
                sep = ";",row.names = F, col.names = T)
  }
  
  #WG
  
  for (j in strata_names) {
    
    ind_WGval <- data.frame(get(paste0("WGpredval_", j)))
    
    write.table(ind_WGval, paste0(path0, "VALID/STRATA/WG/individual_survival/", i, "_ind_WGval_", j, ".csv"),
                sep = ";",row.names = F, col.names = T)
  }
  ######
  
}

iterations <- 1:1000

mclapply(iterations, simulate_iteration, mc.cores = detectCores() - 14)

max_attempts <- 3
attempt <- 1
while (TRUE) {
  indic <- c()  
  
  for (i in 1:1000) {
    if (!file.exists(paste0(path0, "DATAFRAMES/TRAIN/WHOLE/", i, "_dftrain.csv"))) {
      indic <- c(indic, i)
    }
  }
  
  if (length(indic) == 0) {
    break
  }
  if (attempt > max_attempts) {
    break
  }
  mclapply(indic, simulate_iteration, mc.cores = detectCores() - 14)
  attempt <- attempt + 1
}

rm(list = ls(envir = .GlobalEnv), envir = .GlobalEnv)

### 5000SiEND
############## INDICATEURS

### 5000IndSTA
library(parallel)
library(doParallel)


path0 <- "~/Documents/Simulations/Simulations mai 2025/Résultats/N5000_results/NPH/Output Simulations_N5000_NPH_HC/"

indic = c()
for(i in 1:1000){
  if(!file.exists(paste0(path0,"DATAFRAMES/TRAIN/WHOLE/",i,"_dftrain.csv"))){ indic <<- c(indic, i)} 
}


if(!is.null(indic)){
  iterations <- c(1:1000)[-indic]
}else{
  iterations <- c(1:1000)
}
strata_names <- c("HC", "HR", "FC", "FR")
newtimes <- c(365.241, 365.241*3, 365.241*5, 365.241*10) 

### DEBUT ROC.NET ~ligne 1108

start <- Sys.time()

### Calculs des indicateurs

### RMSE & BIAIS
############################################ TRAIN #################################################
#WHOLE

###on récupère les valeurs théoriques moyennes Se(t) 
ALL_theo <- data.frame()

for(k in iterations){
  assign(paste0("meanTrainTheo_", k) , read.csv(paste0(path0, "TRAIN/WHOLE/THEO/mean_survival/",k,"_mean_theo.csv"),
                                                sep = ";") )
  
  ALL_theo <-rbind(ALL_theo, get(paste0("meanTrainTheo_", k)))
  rm(list =paste0("meanTrainTheo_", k))
}

#moyenne théorique de la valur de survie
S_theo <- colMeans(ALL_theo)
############################################
#importation des données des différents modèles 
#########
##plann
ALL_PLANN <- data.frame() 
for(k in iterations){
  assign(paste0("meanTrainPLANN_", k) , read.csv(paste0(path0, "TRAIN/WHOLE/PLANN/mean_survival/",k,"_mean_PLANN.csv"), sep = ";") )
  ALL_PLANN <-rbind(ALL_PLANN, get(paste0("meanTrainPLANN_", k)))
  rm(list =paste0("meanTrainPLANN_", k))
}

###flex1
##2 noeuds
ALL_flex1.2 <- data.frame()
for(k in iterations){
  assign(paste0("meanTrainflex1.2_", k) , read.csv(paste0(path0, "TRAIN/WHOLE/FLEX1.2/mean_survival/",k,"_mean_flex12.csv"), sep = ";") )
  ALL_flex1.2 <- rbind(ALL_flex1.2, get(paste0("meanTrainflex1.2_",k)))
  rm(list =paste0("meanTrainflex1.2_", k))
}

##4 noeuds
ALL_flex1.4 <- data.frame()
for(k in iterations){
  assign(paste0("meanTrainflex1.4_", k) , read.csv(paste0(path0, "TRAIN/WHOLE/FLEX1.4/mean_survival/",k,"_mean_flex14.csv"), sep = ";") )
  ALL_flex1.4 <- rbind(ALL_flex1.4, get(paste0("meanTrainflex1.4_",k)))
  rm(list =paste0("meanTrainflex1.4_", k))
}

###flex2

##2 noeuds
ALL_flex2.2 <- data.frame()
for(k in iterations){
  assign(paste0("meanTrainflex2.2_", k) , read.csv(paste0(path0, "TRAIN/WHOLE/FLEX2.2/mean_survival/",k,"_mean_flex22.csv"), sep = ";") )
  ALL_flex2.2 <- rbind(ALL_flex2.2, get(paste0("meanTrainflex2.2_",k)))
  rm(list =paste0("meanTrainflex2.2_", k))
}

##4 noeuds
ALL_flex2.4 <- data.frame()
for(k in iterations){
  assign(paste0("meanTrainflex2.4_", k) , read.csv(paste0(path0, "TRAIN/WHOLE/FLEX2.4/mean_survival/",k,"_mean_flex24.csv"), sep = ";") )
  ALL_flex2.4 <- rbind(ALL_flex2.4, get(paste0("meanTrainflex2.4_",k)))
  rm(list =paste0("meanTrainflex2.4_", k))
}
### WG
ALL_WG <- data.frame()
for(k in iterations){
  assign(paste0("meanTrainWG_", k) , read.csv(paste0(path0, "TRAIN/WHOLE/WG/mean_survival/",k,"_mean_flexWG.csv"), sep = ";") )
  ALL_WG <- rbind(ALL_WG,get(paste0("meanTrainWG_",k))) 
  rm(list =paste0("meanTrainWG_", k))
}

ALL_PP <- data.frame() 
for(k in iterations){
  assign(paste0("meanTrainPP_", k) , read.csv(paste0(path0, "TRAIN/WHOLE/PP/",k,"PP.csv"), sep = ";") )
  ALL_PP <-rbind(ALL_PP, get(paste0("meanTrainPP_", k)))
  rm(list =paste0("meanTrainPP_", k))
}
ALL_PP <- as.data.frame(matrix(ALL_PP$x, ncol = 4, byrow = TRUE))

biais_PLANN <- colMeans(ALL_PLANN[,-1] - S_theo[col( ALL_PLANN[,-1] )] )
RMSE_PLANN <- sqrt( colMeans((ALL_PLANN[,-1] - S_theo[col( ALL_PLANN[,-1] )])^2) )

biais_flex1.2 <- colMeans( ALL_flex1.2 - S_theo[col( ALL_flex1.2)] )
RMSE_flex1.2 <- sqrt( colMeans( (ALL_flex1.2 - S_theo[col( ALL_flex1.2)])^2 ) )

biais_flex1.4 <- colMeans( ALL_flex1.4 - S_theo[col( ALL_flex1.4)] )
RMSE_flex1.4 <- sqrt( colMeans( (ALL_flex1.4 - S_theo[col( ALL_flex1.4)])^2 ) )

biais_flex2.2 <- colMeans( ALL_flex2.2 - S_theo[col( ALL_flex2.2)] )
RMSE_flex2.2 <- sqrt( colMeans( (ALL_flex2.2 - S_theo[col( ALL_flex2.2)])^2 ) )

biais_flex2.4 <- colMeans( ALL_flex2.4 - S_theo[col( ALL_flex2.4)] )
RMSE_flex2.4 <- sqrt( colMeans( (ALL_flex2.4 - S_theo[col( ALL_flex2.4)])^2 ) )

biais_WG <- colMeans( ALL_WG - S_theo[col( ALL_WG)] )
RMSE_WG <- sqrt( colMeans( (ALL_WG - S_theo[col( ALL_WG)])^2 ) )

biais_PP <- colMeans( ALL_PP - S_theo[col( ALL_PP)] )
RMSE_PP <- sqrt( colMeans( (ALL_PP - S_theo[col( ALL_PP)])^2 ) )
#STRATES

###on récupère les valeurs théoriques moyennes Se(t) 
for(j in strata_names){
  assign(paste0("ALL_",j, "_theo"),  data.frame()) 
}

for(k in iterations){
  assign(paste0("meanTrainTheoS_", k) , read.csv(paste0(path0, "TRAIN/STRATA/THEO/mean_survival/",k,"mean_theo_strates.csv"),
                                                 sep = ";", row.names = 1) ) 
  for(j in strata_names){
    
    current_row <- get(paste0("meanTrainTheoS_", k))[paste0("mean.theo_", j), , drop = FALSE]
    rownames(current_row) <- paste0("mean.theo_", j, k)
    assign(paste0("ALL_", j, "_theo"), rbind(get(paste0("ALL_", j, "_theo")), current_row))
    
  }
  
  rm(list =paste0("meanTrainTheoS_", k))
}

#moyenne théorique de la valeur de survie
for(j in strata_names){
  assign(paste0("S_theo", j ),colMeans(get(paste0("ALL_",j,"_theo") ) ) )
}
############################################
#importation des données des différents modèles 
############################################

# PLANN

for(j in strata_names){
  assign(paste0("ALL_PLANN_", j),data.frame()) 
}

for(k in iterations){
  assign(paste0("meanTrainPLANNstr_", k) , read.csv(paste0(path0, "TRAIN/STRATA/PLANN/mean_survival/",k,"mean_plann_strates.csv"), sep = ";", row.names = 1) )
  
  for(j in strata_names){
    assign(paste0("ALL_PLANN_",j) ,rbind(get(paste0("ALL_PLANN_",j)),get(paste0("meanTrainPLANNstr_",k))[paste0("mean.plann_", j),-1]  ) )
  }
  rm(list =paste0("meanTrainPLANNstr_", k))
}

### flex1

##2 noeuds
for(j in strata_names){
  assign(paste0("ALL_FLEX1.2_", j),data.frame()) 
}

for(k in iterations){
  assign(paste0("meanTrainFLEX1.2str_", k) , read.csv(paste0(path0, "TRAIN/STRATA/FLEX1.2/mean_survival/",k,"mean_flex12_strates.csv"), sep = ";", row.names = 1) )
  
  for(j in strata_names){
    assign(paste0("ALL_FLEX1.2_",j) ,rbind(get(paste0("ALL_FLEX1.2_",j)), get(paste0("meanTrainFLEX1.2str_",k))[paste0("mean.flex1.2_", j),]  ) )
  }
  rm(list =paste0("meanTrainFLEX1.2str_", k))
}

##4 noeuds
for(j in strata_names){
  assign(paste0("ALL_FLEX1.4_", j),data.frame()) 
}

for(k in iterations){
  assign(paste0("meanTrainFLEX1.4str_", k) , read.csv(paste0(path0, "TRAIN/STRATA/FLEX1.4/mean_survival/",k,"mean_flex14_strates.csv"), sep = ";", row.names = 1) )
  
  for(j in strata_names){
    assign(paste0("ALL_FLEX1.4_",j) ,rbind(get(paste0("ALL_FLEX1.4_",j)), get(paste0("meanTrainFLEX1.4str_",k))[paste0("mean.flex1.4_", j),]  ) )
  }
  rm(list =paste0("meanTrainFLEX1.4str_", k))
}

### flex2
##2 noeuds
for(j in strata_names){
  assign(paste0("ALL_FLEX2.2_", j),data.frame()) 
}

for(k in iterations){
  assign(paste0("meanTrainFLEX2.2str_", k) , read.csv(paste0(path0, "TRAIN/STRATA/FLEX2.2/mean_survival/",k,"mean_flex22_strates.csv"), sep = ";", row.names = 1) )
  
  for(j in strata_names){
    assign(paste0("ALL_FLEX2.2_",j) ,rbind(get(paste0("ALL_FLEX2.2_",j)), get(paste0("meanTrainFLEX2.2str_",k))[paste0("mean.flex2.2_", j),]  ) )
  }
  rm(list =paste0("meanTrainFLEX2.2str_", k))
}

##4 noeuds
for(j in strata_names){
  assign(paste0("ALL_FLEX2.4_", j),data.frame()) 
}

for(k in iterations){
  assign(paste0("meanTrainFLEX2.4str_", k) , read.csv(paste0(path0, "TRAIN/STRATA/FLEX2.4/mean_survival/",k,"mean_flex24_strates.csv"), sep = ";", row.names = 1) )
  
  for(j in strata_names){
    assign(paste0("ALL_FLEX2.4_",j) ,rbind(get(paste0("ALL_FLEX2.4_",j)), get(paste0("meanTrainFLEX2.4str_",k))[paste0("mean.flex2.4_", j),]  ) )
  }
  rm(list =paste0("meanTrainFLEX2.4str_", k))
}
### WG

for(j in strata_names){
  assign(paste0("ALL_WG_", j),data.frame()) 
}

for(k in iterations){
  assign(paste0("meanTrainWGstr_", k) , read.csv(paste0(path0, "TRAIN/STRATA/WG/mean_survival/",k,"mean_WG_strates.csv"), sep = ";", row.names = 1) )
  
  for(j in strata_names){
    assign(paste0("ALL_WG_",j) ,rbind(get(paste0("ALL_WG_",j)), get(paste0("meanTrainWGstr_",k))[paste0("mean.WG_", j),]  ) )
  }
  rm(list =paste0("meanTrainWGstr_", k))
}

### PP

for(j in strata_names){
  assign(paste0("ALL_PP_", j),data.frame()) 
}

for(k in iterations){
  assign(paste0("meanTrainPPstr_", k) , read.csv(paste0(path0, "TRAIN/STRATA/PP/",k,"PP_strates.csv"), sep = ";", row.names = 1) )
  
  for(j in strata_names){
    assign(paste0("ALL_PP_",j) ,rbind(get(paste0("ALL_PP_",j)), get(paste0("meanTrainPPstr_",k))[paste0("poharpred_", j),]  ) )
  }
  rm(list =paste0("meanTrainPPstr_", k))
}

#### Calcul des indicateurs
##PLANN
biais_PLANN_strates <- data.frame()
RMSE_PLANN_strates <- data.frame()
for(j in strata_names){
  assign(paste0("biais_PLANN_", j), 
         colMeans( get( paste0("ALL_PLANN_", j) ) - get(paste0("S_theo", j) )[col(get( paste0("ALL_PLANN_", j)) )] ) )
  #on concatene tous les résultats dans une dataframe pour ne pas saturer l'evt 
  biais_PLANN_strates <- rbind(biais_PLANN_strates, get(paste0("biais_PLANN_", j)) )
  rownames(biais_PLANN_strates)[dim(biais_PLANN_strates)[1]] <- paste0(j)
  rm(list = paste0("biais_PLANN_", j))
  
  assign(paste0("RMSE_PLANN_",  j) ,
         sqrt( colMeans( ( get( paste0("ALL_PLANN_", j) ) - 
                             get(paste0("S_theo", j) )[col(get( paste0("ALL_PLANN_", j)) )] )^2 ) ) )
  RMSE_PLANN_strates <- rbind(RMSE_PLANN_strates, get(paste0("RMSE_PLANN_", j)) )
  rownames(RMSE_PLANN_strates)[dim(RMSE_PLANN_strates)[1]] <- paste0(j)
  rm(list = paste0("RMSE_PLANN_", j))
}

colnames(biais_PLANN_strates) <- colnames(ALL_flex1.4)
colnames(RMSE_PLANN_strates) <- colnames(ALL_flex1.4)

# FLEX 1
#2 noeuds
biais_FLEX1.2_strates <- data.frame()
RMSE_FLEX1.2_strates <- data.frame()
for(j in strata_names){
  assign(paste0("biais_FLEX1.2_", j), 
         colMeans( get( paste0("ALL_FLEX1.2_", j) ) - get(paste0("S_theo", j) )[col(get( paste0("ALL_FLEX1.2_", j)) )] ) )
  
  #on concatene tous les résultats dans une dataframe pour ne pas saturer l'evt 
  biais_FLEX1.2_strates <- rbind(biais_FLEX1.2_strates, get(paste0("biais_FLEX1.2_", j)) )
  rownames(biais_FLEX1.2_strates)[dim(biais_FLEX1.2_strates)[1]] <- paste0(j)
  rm(list = paste0("biais_FLEX1.2_", j))
  
  assign(paste0("RMSE_FLEX1.2_",  j) ,
         sqrt( colMeans( ( get( paste0("ALL_FLEX1.2_", j) ) - 
                             get(paste0("S_theo", j) )[col(get( paste0("ALL_FLEX1.2_", j)) )] )^2 ) ) )
  
  RMSE_FLEX1.2_strates <- rbind(RMSE_FLEX1.2_strates, get(paste0("RMSE_FLEX1.2_", j)) )
  rownames(RMSE_FLEX1.2_strates)[dim(RMSE_FLEX1.2_strates)[1]] <- paste0(j)
  rm(list = paste0("RMSE_FLEX1.2_", j))
  
}

colnames(biais_FLEX1.2_strates) <- colnames(ALL_flex1.2)
colnames(RMSE_FLEX1.2_strates) <- colnames(ALL_flex1.2)

#4 noeuds
biais_FLEX1.4_strates <- data.frame()
RMSE_FLEX1.4_strates <- data.frame()
for(j in strata_names){
  assign(paste0("biais_FLEX1.4_", j), 
         colMeans( get( paste0("ALL_FLEX1.4_", j) ) - get(paste0("S_theo", j) )[col(get( paste0("ALL_FLEX1.4_", j)) )] ) )
  
  #on concatene tous les résultats dans une dataframe pour ne pas saturer l'evt 
  biais_FLEX1.4_strates <- rbind(biais_FLEX1.4_strates, get(paste0("biais_FLEX1.4_", j)) )
  rownames(biais_FLEX1.4_strates)[dim(biais_FLEX1.4_strates)[1]] <- paste0(j)
  rm(list = paste0("biais_FLEX1.4_", j))
  
  assign(paste0("RMSE_FLEX1.4_",  j) ,
         sqrt( colMeans( ( get( paste0("ALL_FLEX1.4_", j) ) - 
                             get(paste0("S_theo", j) )[col(get( paste0("ALL_FLEX1.4_", j)) )] )^2 ) ) )
  
  RMSE_FLEX1.4_strates <- rbind(RMSE_FLEX1.4_strates, get(paste0("RMSE_FLEX1.4_", j)) )
  rownames(RMSE_FLEX1.4_strates)[dim(RMSE_FLEX1.4_strates)[1]] <- paste0(j)
  rm(list = paste0("RMSE_FLEX1.4_", j))
  
}

colnames(biais_FLEX1.4_strates) <- colnames(ALL_flex1.4)
colnames(RMSE_FLEX1.4_strates) <- colnames(ALL_flex1.4)


# FLEX 2
##2 noeuds
biais_FLEX2.2_strates <- data.frame()
RMSE_FLEX2.2_strates <- data.frame()
for(j in strata_names){
  assign(paste0("biais_FLEX2.2_", j), 
         colMeans( get( paste0("ALL_FLEX2.2_", j) ) - get(paste0("S_theo", j) )[col(get( paste0("ALL_FLEX2.2_", j)) )] ) )
  
  biais_FLEX2.2_strates <- rbind(biais_FLEX2.2_strates, get(paste0("biais_FLEX2.2_", j)) )
  rownames(biais_FLEX2.2_strates)[dim(biais_FLEX2.2_strates)[1]] <- paste0(j)
  rm(list = paste0("biais_FLEX2.2_", j))
  
  assign(paste0("RMSE_FLEX2.2_",  j) ,
         sqrt( colMeans( ( get( paste0("ALL_FLEX2.2_", j) ) - 
                             get(paste0("S_theo", j) )[col(get( paste0("ALL_FLEX2.2_", j)) )] )^2 ) ) )
  
  RMSE_FLEX2.2_strates <- rbind(RMSE_FLEX2.2_strates, get(paste0("RMSE_FLEX2.2_", j)) )
  rownames(RMSE_FLEX2.2_strates)[dim(RMSE_FLEX2.2_strates)[1]] <- paste0(j)
  rm(list = paste0("RMSE_FLEX2.2_", j))
}
colnames(biais_FLEX2.2_strates) <- colnames(ALL_flex2.2)
colnames(RMSE_FLEX2.2_strates) <- colnames(ALL_flex2.2)

##4 noeuds
biais_FLEX2.4_strates <- data.frame()
RMSE_FLEX2.4_strates <- data.frame()
for(j in strata_names){
  assign(paste0("biais_FLEX2.4_", j), 
         colMeans( get( paste0("ALL_FLEX2.4_", j) ) - get(paste0("S_theo", j) )[col(get( paste0("ALL_FLEX2.4_", j)) )] ) )
  
  biais_FLEX2.4_strates <- rbind(biais_FLEX2.4_strates, get(paste0("biais_FLEX2.4_", j)) )
  rownames(biais_FLEX2.4_strates)[dim(biais_FLEX2.4_strates)[1]] <- paste0(j)
  rm(list = paste0("biais_FLEX2.4_", j))
  
  assign(paste0("RMSE_FLEX2.4_",  j) ,
         sqrt( colMeans( ( get( paste0("ALL_FLEX2.4_", j) ) - 
                             get(paste0("S_theo", j) )[col(get( paste0("ALL_FLEX2.4_", j)) )] )^2 ) ) )
  
  RMSE_FLEX2.4_strates <- rbind(RMSE_FLEX2.4_strates, get(paste0("RMSE_FLEX2.4_", j)) )
  rownames(RMSE_FLEX2.4_strates)[dim(RMSE_FLEX2.4_strates)[1]] <- paste0(j)
  rm(list = paste0("RMSE_FLEX2.4_", j))
}
colnames(biais_FLEX2.4_strates) <- colnames(ALL_flex2.4)
colnames(RMSE_FLEX2.4_strates) <- colnames(ALL_flex2.4)

# WG 
biais_WG_strates <- data.frame()
RMSE_WG_strates <- data.frame()
for(j in strata_names){
  assign(paste0("biais_WG_", j), 
         colMeans( get( paste0("ALL_WG_", j) ) - get(paste0("S_theo", j) )[col(get( paste0("ALL_WG_", j)) )] ) )
  
  biais_WG_strates <- rbind(biais_WG_strates, get(paste0("biais_WG_", j)) )
  rownames(biais_WG_strates)[dim(biais_WG_strates)[1]] <- paste0(j)
  rm(list = paste0("biais_WG_", j))
  
  assign(paste0("RMSE_WG_",  j) ,
         sqrt( colMeans( ( get( paste0("ALL_WG_", j) ) - 
                             get(paste0("S_theo", j) )[col(get( paste0("ALL_WG_", j)) )] )^2 ) ) )
  
  RMSE_WG_strates <- rbind(RMSE_WG_strates, get(paste0("RMSE_WG_", j)) )
  rownames(RMSE_WG_strates)[dim(RMSE_WG_strates)[1]] <- paste0(j)
  rm(list = paste0("RMSE_WG_", j))
}

colnames(biais_WG_strates) <- colnames(ALL_WG)
colnames(RMSE_WG_strates) <- colnames(ALL_WG)

# PP
biais_PP_strates <- data.frame()
RMSE_PP_strates <- data.frame()
for(j in strata_names){
  assign(paste0("biais_PP_", j), 
         colMeans( get( paste0("ALL_PP_", j) ) - get(paste0("S_theo", j) )[col(get( paste0("ALL_PP_", j)) )] ) )
  
  biais_PP_strates <- rbind(biais_PP_strates, get(paste0("biais_PP_", j)) )
  rownames(biais_PP_strates)[dim(biais_PP_strates)[1]] <- paste0(j)
  rm(list = paste0("biais_PP_", j))
  
  assign(paste0("RMSE_PP_",  j) ,
         sqrt( colMeans( ( get( paste0("ALL_PP_", j) ) - 
                             get(paste0("S_theo", j) )[col(get( paste0("ALL_PP_", j)) )] )^2 ) ) )
  
  RMSE_PP_strates <- rbind(RMSE_PP_strates, get(paste0("RMSE_PP_", j)) )
  rownames(RMSE_PP_strates)[dim(RMSE_PP_strates)[1]] <- paste0(j)
  rm(list = paste0("RMSE_PP_", j))
}

colnames(biais_PP_strates) <- colnames(ALL_PP)
colnames(RMSE_PP_strates) <- colnames(ALL_PP)

###concaten des S_theo, là aussi pour libérer de la place dans l'evt global
S_theo_strates <- data.frame()
for(j in strata_names){
  S_theo_strates <- rbind(S_theo_strates, get(paste0("S_theo",j)))
  rownames(S_theo_strates)[dim(S_theo_strates)[1]] <- paste0(j)
  rm(list = paste0("S_theo",j))
}
colnames(S_theo_strates) <- colnames(ALL_WG)

############################################ VALID #################################################

# WHOLE
###on récupère les valeurs théoriques moyennes Se(t) 
ALL_theoval <- data.frame()

for(k in iterations){
  assign(paste0("meanValidTheo_", k) , read.csv(paste0(path0, "VALID/WHOLE/THEO/mean_survival/",k,"_mean_theoval.csv"),
                                                sep = ";") )
  
  ALL_theoval <-rbind(ALL_theoval, get(paste0("meanValidTheo_", k)))
  rm(list =paste0("meanValidTheo_", k))
}

#moyenne théorique de la valur de survie
S_theoval <- colMeans(ALL_theoval)

############################################
#importation des données des différents modèles 
#########
##plann
ALL_PLANNval <- data.frame() 
for(k in iterations){
  assign(paste0("meanValidPLANN_", k) , read.csv(paste0(path0, "VALID/WHOLE/PLANN/mean_survival/",k,"_mean_PLANNval.csv"), sep = ";") )
  ALL_PLANNval <-rbind(ALL_PLANNval, get(paste0("meanValidPLANN_", k)))
  rm(list =paste0("meanValidPLANN_", k))
}

###flex1
##2 noeuds
ALL_flex1.2val <- data.frame()
for(k in iterations){
  assign(paste0("meanValidflex1.2_", k) , read.csv(paste0(path0, "VALID/WHOLE/FLEX1.2/mean_survival/",k,"_mean_flexval12.csv"), sep = ";") )
  ALL_flex1.2val <- rbind(ALL_flex1.2val, get(paste0("meanValidflex1.2_",k)))
  rm(list =paste0("meanValidflex1.2_", k))
}
##4 noeuds
ALL_flex1.4val <- data.frame()
for(k in iterations){
  assign(paste0("meanValidflex1.4_", k) , read.csv(paste0(path0, "VALID/WHOLE/FLEX1.4/mean_survival/",k,"_mean_flexval14.csv"), sep = ";") )
  ALL_flex1.4val <- rbind(ALL_flex1.4val, get(paste0("meanValidflex1.4_",k)))
  rm(list =paste0("meanValidflex1.4_", k))
}
### flex2

##2 noeuds
ALL_flex2.2val <- data.frame()
for(k in iterations){
  assign(paste0("meanValidflex2.2_", k) , read.csv(paste0(path0, "VALID/WHOLE/FLEX2.2/mean_survival/",k,"_mean_flexval22.csv"), sep = ";") )
  ALL_flex2.2val <- rbind(ALL_flex2.2val, get(paste0("meanValidflex2.2_",k)))
  rm(list =paste0("meanValidflex2.2_", k))
}

##4 noeuds
ALL_flex2.4val <- data.frame()
for(k in iterations){
  assign(paste0("meanValidflex2.4_", k) , read.csv(paste0(path0, "VALID/WHOLE/FLEX2.4/mean_survival/",k,"_mean_flexval24.csv"), sep = ";") )
  ALL_flex2.4val <- rbind(ALL_flex2.4val, get(paste0("meanValidflex2.4_",k)))
  rm(list =paste0("meanValidflex2.4_", k))
}
### WG
ALL_WGval <- data.frame()
for(k in iterations){
  assign(paste0("meanValidWG_", k) , read.csv(paste0(path0, "VALID/WHOLE/WG/mean_survival/",k,"_mean_flexWGval.csv"), sep = ";") )
  ALL_WGval <- rbind(ALL_WGval, get(paste0("meanValidWG_",k))) 
  rm(list =paste0("meanValidWG_", k))
}

##PP
ALL_PPval <- data.frame()
for(k in iterations){
  assign(paste0("meanValidPP_", k) , read.csv(paste0(path0, "VALID/WHOLE/PP/",k,"PPval.csv"), sep = ";") )
  ALL_PPval <- rbind(ALL_PPval, get(paste0("meanValidPP_",k))) 
  rm(list =paste0("meanValidPP_", k))
}
ALL_PPval <- as.data.frame(matrix(ALL_PPval$x, ncol = 4, byrow = TRUE))

biais_PLANNval <- colMeans(ALL_PLANNval[,-1] - S_theoval[col( ALL_PLANNval[,-1] )] )
RMSE_PLANNval <- sqrt( colMeans((ALL_PLANNval[,-1] - S_theoval[col( ALL_PLANNval[,-1] )])^2) )

biais_flex1.2val <- colMeans( ALL_flex1.2val - S_theoval[col( ALL_flex1.2val)] )
RMSE_flex1.2val <- sqrt( colMeans( (ALL_flex1.2val - S_theoval[col( ALL_flex1.2val)])^2 ) )

biais_flex1.4val <- colMeans( ALL_flex1.4val - S_theoval[col( ALL_flex1.4val)] )
RMSE_flex1.4val <- sqrt( colMeans( (ALL_flex1.4val - S_theoval[col( ALL_flex1.4val)])^2 ) )

biais_flex2.2val <- colMeans( ALL_flex2.2val - S_theoval[col( ALL_flex2.2val)] )
RMSE_flex2.2val <- sqrt( colMeans( (ALL_flex2.2val - S_theoval[col( ALL_flex2.2val)])^2 ) )

biais_flex2.4val <- colMeans( ALL_flex2.4val - S_theoval[col( ALL_flex2.4val)] )
RMSE_flex2.4val <- sqrt( colMeans( (ALL_flex2.4val - S_theoval[col( ALL_flex2.4val)])^2 ) )

biais_WGval <- colMeans( ALL_WGval - S_theoval[col( ALL_WGval)] )
RMSE_WGval <- sqrt( colMeans( (ALL_WGval - S_theoval[col( ALL_WGval)])^2 ) )

biais_PPval <- colMeans( ALL_PP - S_theoval[col( ALL_PPval)] )
RMSE_PPval <- sqrt( colMeans( (ALL_PP - S_theoval[col( ALL_PPval)])^2 ) )
#STRATES

###on récupère les valeurs théoriques moyennes Se(t) 
for(j in strata_names){
  assign(paste0("ALL_",j, "_theoval"),  data.frame()) 
}

for(k in iterations){
  assign(paste0("meanValidTheoS_", k) , read.csv(paste0(path0, "VALID/STRATA/THEO/mean_survival/",k,"mean_theoval_strates.csv"),
                                                 sep = ";", row.names = 1)[,-1] ) ### [,-1] car quand j'ai enregistré les résultats, j'ai dupliué la premiere colonne donc il faut l'enlever
  for(j in strata_names){
    
    current_row <- get(paste0("meanValidTheoS_", k))[paste0("mean.theo_", j), , drop = FALSE]
    rownames(current_row) <- paste0("mean.theoval_", j, k)
    assign(paste0("ALL_", j, "_theoval"), rbind(get(paste0("ALL_", j, "_theoval")), current_row))
    
  }
  
  rm(list =paste0("meanValidTheoS_", k))
}
rm(current_row)
#moyenne théorique de la valeur de survie
for(j in strata_names){
  assign(paste0("S_theoval", j ),colMeans(get(paste0("ALL_",j,"_theoval") ) ) )
}
############################################
#importation des données des différents modèles 
############################################

# PLANN

for(j in strata_names){
  assign(paste0("ALL_PLANNval_", j),data.frame()) 
}

for(k in iterations){
  assign(paste0("meanValidPLANNstr_", k) , read.csv(paste0(path0, "VALID/STRATA/PLANN/mean_survival/",k,"mean_plannval_strates.csv"), sep = ";", row.names = 1) )
  
  for(j in strata_names){
    assign(paste0("ALL_PLANNval_",j) ,rbind(get(paste0("ALL_PLANNval_",j)),get(paste0("meanValidPLANNstr_",k))[paste0("mean.plannval_", j),-1]  ) )
  }
  rm(list =paste0("meanValidPLANNstr_", k))
}

### flex1
##2 noeuds
for(j in strata_names){
  assign(paste0("ALL_FLEX1.2val_", j),data.frame()) 
}

for(k in iterations){
  assign(paste0("meanValidFLEX1.2str_", k) , read.csv(paste0(path0, "VALID/STRATA/FLEX1.2/mean_survival/",k,"mean_flexval12_strates.csv"), sep = ";", row.names = 1) )
  
  for(j in strata_names){
    assign(paste0("ALL_FLEX1.2val_",j) ,rbind(get(paste0("ALL_FLEX1.2val_",j)), get(paste0("meanValidFLEX1.2str_",k))[paste0("mean.flexval1.2_", j),]  ) )
  }
  rm(list =paste0("meanValidFLEX1.2str_", k))
}

##4 noeuds
for(j in strata_names){
  assign(paste0("ALL_FLEX1.4val_", j),data.frame()) 
}

for(k in iterations){
  assign(paste0("meanValidFLEX1.4str_", k) , read.csv(paste0(path0, "VALID/STRATA/FLEX1.4/mean_survival/",k,"mean_flexval14_strates.csv"), sep = ";", row.names = 1) )
  
  for(j in strata_names){
    assign(paste0("ALL_FLEX1.4val_",j) ,rbind(get(paste0("ALL_FLEX1.4val_",j)), get(paste0("meanValidFLEX1.4str_",k))[paste0("mean.flexval1.4_", j),]  ) )
  }
  rm(list =paste0("meanValidFLEX1.4str_", k))
}
### flex2

##2 noeuds
for(j in strata_names){
  assign(paste0("ALL_FLEX2.2val_", j),data.frame()) 
}

for(k in iterations){
  assign(paste0("meanValidFLEX2.2str_", k) , read.csv(paste0(path0, "VALID/STRATA/FLEX2.2/mean_survival/",k,"mean_flexval22_strates.csv"), sep = ";", row.names = 1) )
  
  for(j in strata_names){
    assign(paste0("ALL_FLEX2.2val_",j) ,rbind(get(paste0("ALL_FLEX2.2val_",j)), get(paste0("meanValidFLEX2.2str_",k))[paste0("mean.flexval2.2_", j),]  ) )
  }
  rm(list =paste0("meanValidFLEX2.2str_", k))
}

##4 noeuds
for(j in strata_names){
  assign(paste0("ALL_FLEX2.4val_", j),data.frame()) 
}

for(k in iterations){
  assign(paste0("meanValidFLEX2.4str_", k) , read.csv(paste0(path0, "VALID/STRATA/FLEX2.4/mean_survival/",k,"mean_flexval24_strates.csv"), sep = ";", row.names = 1) )
  
  for(j in strata_names){
    assign(paste0("ALL_FLEX2.4val_",j) ,rbind(get(paste0("ALL_FLEX2.4val_",j)), get(paste0("meanValidFLEX2.4str_",k))[paste0("mean.flexval2.4_", j),]  ) )
  }
  rm(list =paste0("meanValidFLEX2.4str_", k))
}
### WG

for(j in strata_names){
  assign(paste0("ALL_WGval_", j),data.frame()) 
}

for(k in iterations){
  assign(paste0("meanValidWGstr_", k) , read.csv(paste0(path0, "VALID/STRATA/WG/mean_survival/",k,"mean_WGval_strates.csv"), sep = ";", row.names = 1) )
  
  for(j in strata_names){
    assign(paste0("ALL_WGval_",j) ,rbind(get(paste0("ALL_WGval_",j)), get(paste0("meanValidWGstr_",k))[paste0("mean.WGval_", j),]  ) )
  }
  rm(list =paste0("meanValidWGstr_", k))
}

### PP 

for(j in strata_names){
  assign(paste0("ALL_PPval_", j),data.frame()) 
}

for(k in iterations){
  assign(paste0("meanValidPPstr_", k) , read.csv(paste0(path0, "VALID/STRATA/PP/",k,"PPval_strates.csv"), sep = ";", row.names = 1) )
  
  for(j in strata_names){
    assign(paste0("ALL_PPval_",j) ,rbind(get(paste0("ALL_PPval_",j)), get(paste0("meanValidPPstr_",k))[paste0("poharpredval_", j),]  ) )
  }
  rm(list =paste0("meanValidPPstr_", k))
}
#### Calcul des indicateurs
##PLANN
biais_PLANNval_strates <- data.frame()
RMSE_PLANNval_strates <- data.frame()
for(j in strata_names){
  assign(paste0("biais_PLANNval_", j), 
         colMeans( get( paste0("ALL_PLANNval_", j) ) - get(paste0("S_theoval", j) )[col(get( paste0("ALL_PLANNval_", j)) )] ) )
  
  biais_PLANNval_strates <- rbind(biais_PLANNval_strates, get(paste0("biais_PLANNval_", j)) )
  rownames(biais_PLANNval_strates)[dim(biais_PLANNval_strates)[1]] <- paste0(j)
  rm(list = paste0("biais_PLANNval_", j))
  
  assign(paste0("RMSE_PLANNval_",  j) ,
         sqrt( colMeans( ( get( paste0("ALL_PLANNval_", j) ) - 
                             get(paste0("S_theoval", j) )[col(get( paste0("ALL_PLANNval_", j)) )] )^2 ) ) )
  
  RMSE_PLANNval_strates <- rbind(RMSE_PLANNval_strates, get(paste0("RMSE_PLANNval_", j)) )
  rownames(RMSE_PLANNval_strates)[dim(RMSE_PLANNval_strates)[1]] <- paste0(j)
  rm(list = paste0("RMSE_PLANNval_", j))
  
}

colnames(biais_PLANNval_strates) <- colnames(ALL_flex1.2)
colnames(RMSE_PLANNval_strates) <- colnames(ALL_flex1.2)

# FLEX 1

##2 noeuds
biais_FLEX1.2val_strates <- data.frame()
RMSE_FLEX1.2val_strates <- data.frame()
for(j in strata_names){
  assign(paste0("biais_FLEX1.2val_", j), 
         colMeans( get( paste0("ALL_FLEX1.2val_", j) ) - get(paste0("S_theoval", j) )[col(get( paste0("ALL_FLEX1.2val_", j)) )] ) )
  
  biais_FLEX1.2val_strates <- rbind(biais_FLEX1.2val_strates, get(paste0("biais_FLEX1.2val_", j)) )
  rownames(biais_FLEX1.2val_strates)[dim(biais_FLEX1.2val_strates)[1]] <- paste0(j)
  rm(list = paste0("biais_FLEX1.2val_", j))
  
  assign(paste0("RMSE_FLEX1.2val_",  j) ,
         sqrt( colMeans( ( get( paste0("ALL_FLEX1.2val_", j) ) - 
                             get(paste0("S_theoval", j) )[col(get( paste0("ALL_FLEX1.2val_", j)) )] )^2 ) ) )
  
  RMSE_FLEX1.2val_strates <- rbind(RMSE_FLEX1.2val_strates, get(paste0("RMSE_FLEX1.2val_", j)) )
  rownames(RMSE_FLEX1.2val_strates)[dim(RMSE_FLEX1.2val_strates)[1]] <- paste0(j)
  rm(list = paste0("RMSE_FLEX1.2val_", j))
}

colnames(biais_FLEX1.2val_strates) <- colnames(ALL_flex1.2)
colnames(RMSE_FLEX1.2val_strates) <- colnames(ALL_flex1.2)

##4 noeuds
biais_FLEX1.4val_strates <- data.frame()
RMSE_FLEX1.4val_strates <- data.frame()
for(j in strata_names){
  assign(paste0("biais_FLEX1.4val_", j), 
         colMeans( get( paste0("ALL_FLEX1.4val_", j) ) - get(paste0("S_theoval", j) )[col(get( paste0("ALL_FLEX1.4val_", j)) )] ) )
  
  biais_FLEX1.4val_strates <- rbind(biais_FLEX1.4val_strates, get(paste0("biais_FLEX1.4val_", j)) )
  rownames(biais_FLEX1.4val_strates)[dim(biais_FLEX1.4val_strates)[1]] <- paste0(j)
  rm(list = paste0("biais_FLEX1.4val_", j))
  
  assign(paste0("RMSE_FLEX1.4val_",  j) ,
         sqrt( colMeans( ( get( paste0("ALL_FLEX1.4val_", j) ) - 
                             get(paste0("S_theoval", j) )[col(get( paste0("ALL_FLEX1.4val_", j)) )] )^2 ) ) )
  
  RMSE_FLEX1.4val_strates <- rbind(RMSE_FLEX1.4val_strates, get(paste0("RMSE_FLEX1.4val_", j)) )
  rownames(RMSE_FLEX1.4val_strates)[dim(RMSE_FLEX1.4val_strates)[1]] <- paste0(j)
  rm(list = paste0("RMSE_FLEX1.4val_", j))
}

colnames(biais_FLEX1.4val_strates) <- colnames(ALL_flex1.4)
colnames(RMSE_FLEX1.4val_strates) <- colnames(ALL_flex1.4)

# FLEX 2

##2 noeuds
biais_FLEX2.2val_strates <- data.frame()
RMSE_FLEX2.2val_strates <- data.frame()

for(j in strata_names){
  assign(paste0("biais_FLEX2.2val_", j), 
         colMeans( get( paste0("ALL_FLEX2.2val_", j) ) - get(paste0("S_theoval", j) )[col(get( paste0("ALL_FLEX2.2val_", j)) )] ) )
  
  biais_FLEX2.2val_strates <- rbind(biais_FLEX2.2val_strates, get(paste0("biais_FLEX2.2val_", j)) )
  rownames(biais_FLEX2.2val_strates)[dim(biais_FLEX2.2val_strates)[1]] <- paste0(j)
  rm(list = paste0("biais_FLEX2.2val_", j))
  
  assign(paste0("RMSE_FLEX2.2val_",  j) ,
         sqrt( colMeans( ( get( paste0("ALL_FLEX2.2val_", j) ) - 
                             get(paste0("S_theoval", j) )[col(get( paste0("ALL_FLEX2.2val_", j)) )] )^2 ) ) )
  
  RMSE_FLEX2.2val_strates <- rbind(RMSE_FLEX2.2val_strates, get(paste0("RMSE_FLEX2.2val_", j)) )
  rownames(RMSE_FLEX2.2val_strates)[dim(RMSE_FLEX2.2val_strates)[1]] <- paste0(j)
  rm(list = paste0("RMSE_FLEX2.2val_", j))
  
}

colnames(biais_FLEX2.2val_strates) <- colnames(ALL_flex1.2)
colnames(RMSE_FLEX2.2val_strates) <- colnames(ALL_flex1.2)

##4 noeuds
biais_FLEX2.4val_strates <- data.frame()
RMSE_FLEX2.4val_strates <- data.frame()

for(j in strata_names){
  assign(paste0("biais_FLEX2.4val_", j), 
         colMeans( get( paste0("ALL_FLEX2.4val_", j) ) - get(paste0("S_theoval", j) )[col(get( paste0("ALL_FLEX2.4val_", j)) )] ) )
  
  biais_FLEX2.4val_strates <- rbind(biais_FLEX2.4val_strates, get(paste0("biais_FLEX2.4val_", j)) )
  rownames(biais_FLEX2.4val_strates)[dim(biais_FLEX2.4val_strates)[1]] <- paste0(j)
  rm(list = paste0("biais_FLEX2.4val_", j))
  
  assign(paste0("RMSE_FLEX2.4val_",  j) ,
         sqrt( colMeans( ( get( paste0("ALL_FLEX2.4val_", j) ) - 
                             get(paste0("S_theoval", j) )[col(get( paste0("ALL_FLEX2.4val_", j)) )] )^2 ) ) )
  
  RMSE_FLEX2.4val_strates <- rbind(RMSE_FLEX2.4val_strates, get(paste0("RMSE_FLEX2.4val_", j)) )
  rownames(RMSE_FLEX2.4val_strates)[dim(RMSE_FLEX2.4val_strates)[1]] <- paste0(j)
  rm(list = paste0("RMSE_FLEX2.4val_", j))
  
}

colnames(biais_FLEX2.4val_strates) <- colnames(ALL_flex1.2)
colnames(RMSE_FLEX2.4val_strates) <- colnames(ALL_flex1.2)
# WG 

biais_WGval_strates <- data.frame()
RMSE_WGval_strates <- data.frame()

for(j in strata_names){
  assign(paste0("biais_WGval_", j), 
         colMeans( get( paste0("ALL_WGval_", j) ) - get(paste0("S_theoval", j) )[col(get( paste0("ALL_WGval_", j)) )] ) )
  
  biais_WGval_strates <- rbind(biais_WGval_strates, get(paste0("biais_WGval_", j)) )
  rownames(biais_WGval_strates)[dim(biais_WGval_strates)[1]] <- paste0(j)
  rm(list = paste0("biais_WGval_", j))
  
  assign(paste0("RMSE_WGval_",  j) ,
         sqrt( colMeans( ( get( paste0("ALL_WGval_", j) ) - 
                             get(paste0("S_theoval", j) )[col(get( paste0("ALL_WGval_", j)) )] )^2 ) ) )
  
  RMSE_WGval_strates <- rbind(RMSE_WGval_strates, get(paste0("RMSE_WGval_", j)) )
  rownames(RMSE_WGval_strates)[dim(RMSE_WGval_strates)[1]] <- paste0(j)
  rm(list = paste0("RMSE_WGval_", j))
  
}

colnames(biais_WGval_strates) <- colnames(ALL_flex1.2)
colnames(RMSE_WGval_strates) <- colnames(ALL_flex1.2)

# PP ############################################## A CHANGER 
biais_PPval_strates <- data.frame()
RMSE_PPval_strates <- data.frame()
for(j in strata_names){
  assign(paste0("biais_PPval_", j), 
         colMeans( get( paste0("ALL_PP_", j) ) - get(paste0("S_theoval", j) )[col(get( paste0("ALL_PP_", j)) )] ) )
  
  biais_PPval_strates <- rbind(biais_PPval_strates, get(paste0("biais_PPval_", j)) )
  rownames(biais_PPval_strates)[dim(biais_PPval_strates)[1]] <- paste0(j)
  rm(list = paste0("biais_PPval_", j))
  
  assign(paste0("RMSE_PPval_",  j) ,
         sqrt( colMeans( ( get( paste0("ALL_PP_", j) ) - 
                             get(paste0("S_theoval", j) )[col(get( paste0("ALL_PP_", j)) )] )^2 ) ) )
  
  RMSE_PPval_strates <- rbind(RMSE_PPval_strates, get(paste0("RMSE_PPval_", j)) )
  rownames(RMSE_PPval_strates)[dim(RMSE_PPval_strates)[1]] <- paste0(j)
  rm(list = paste0("RMSE_PPval_", j))
}

colnames(biais_PPval_strates) <- colnames(ALL_PPval)
colnames(RMSE_PPval_strates) <- colnames(ALL_PPval)

###concaten des S_theo, là aussi pour libérer de la place dans l'evt global
S_theoval_strates <- data.frame()
for(j in strata_names){
  S_theoval_strates <- rbind(S_theoval_strates, get(paste0("S_theoval",j)))
  rownames(S_theoval_strates)[dim(S_theoval_strates)[1]] <- paste0(j)
  rm(list = paste0("S_theoval",j))
}
colnames(S_theoval_strates) <- colnames(ALL_WG)


## mise en forme de l'evt global 
##############################################
##concatenation STHEO/ STHEOSTRATES
##TRAIN
S_theo <- rbind(S_theo, S_theo_strates)
rownames(S_theo)[1] <- "ALL" 
rm(S_theo_strates)
##VALID

S_theoval <- rbind(S_theoval, S_theoval_strates)
rownames(S_theoval)[1] <- "ALL" 
rm(S_theoval_strates)
## concatenation des résultats sur la base totale avec ceux fait en strates 
##TRAIN
##PLANN
biais_PLANN <- rbind(biais_PLANN, biais_PLANN_strates)
RMSE_PLANN <- rbind(RMSE_PLANN, RMSE_PLANN_strates)
rownames(biais_PLANN)[1] <- "ALL"
rownames(RMSE_PLANN)[1] <- "ALL"
rm(biais_PLANN_strates)
rm(RMSE_PLANN_strates)
##FLEX1
##2noeuds
biais_FLEX1.2 <- rbind(biais_flex1.2, biais_FLEX1.2_strates)
RMSE_FLEX1.2 <- rbind(RMSE_flex1.2, RMSE_FLEX1.2_strates)
rownames(biais_FLEX1.2)[1] <- "ALL"
rownames(RMSE_FLEX1.2)[1] <- "ALL"
rm(biais_flex1.2, biais_FLEX1.2_strates)
rm(RMSE_flex1.2, RMSE_FLEX1.2_strates)

##4noeuds
biais_FLEX1.4 <- rbind(biais_flex1.4, biais_FLEX1.4_strates)
RMSE_FLEX1.4 <- rbind(RMSE_flex1.4, RMSE_FLEX1.4_strates)
rownames(biais_FLEX1.4)[1] <- "ALL"
rownames(RMSE_FLEX1.4)[1] <- "ALL"
rm(biais_flex1.4, biais_FLEX1.4_strates)
rm(RMSE_flex1.4, RMSE_FLEX1.4_strates)

##FLEX2
##2 noeuds
biais_FLEX2.2 <- rbind(biais_flex2.2, biais_FLEX2.2_strates)
RMSE_FLEX2.2 <- rbind(RMSE_flex2.2, RMSE_FLEX2.2_strates)
rownames(biais_FLEX2.2)[1] <- "ALL"
rownames(RMSE_FLEX2.2)[1] <- "ALL"
rm(biais_flex2.2, biais_FLEX2.2_strates)
rm(RMSE_flex2.2, RMSE_FLEX2.2_strates)
##4 noeuds
biais_FLEX2.4 <- rbind(biais_flex2.4, biais_FLEX2.4_strates)
RMSE_FLEX2.4 <- rbind(RMSE_flex2.4, RMSE_FLEX2.4_strates)
rownames(biais_FLEX2.4)[1] <- "ALL"
rownames(RMSE_FLEX2.4)[1] <- "ALL"
rm(biais_flex2.4, biais_FLEX2.4_strates)
rm(RMSE_flex2.4, RMSE_FLEX2.4_strates)
##WG

biais_WG <- rbind(biais_WG, biais_WG_strates)
RMSE_WG <- rbind(RMSE_WG, RMSE_WG_strates)
rownames(biais_WG)[1] <- "ALL"
rownames(RMSE_WG)[1] <- "ALL"
rm(biais_WG_strates)
rm(RMSE_WG_strates)

##PP

biais_PP <- rbind(biais_PP, biais_PP_strates)
RMSE_PP <- rbind(RMSE_PP, RMSE_PP_strates)
rownames(biais_PP)[1] <- "ALL"
rownames(RMSE_PP)[1] <- "ALL"
rm(biais_PP_strates)
rm(RMSE_PP_strates)

######VALID

biais_PLANNval <- rbind(biais_PLANNval, biais_PLANNval_strates)
RMSE_PLANNval <- rbind(RMSE_PLANNval, RMSE_PLANNval_strates)
rownames(biais_PLANNval)[1] <- "ALL"
rownames(RMSE_PLANNval)[1] <- "ALL"
rm(biais_PLANNval_strates)
rm(RMSE_PLANNval_strates)
##FLEX1
##2 noeuds
biais_FLEX1.2val <- rbind(biais_flex1.2val, biais_FLEX1.2val_strates)
RMSE_FLEX1.2val <- rbind(RMSE_flex1.2val, RMSE_FLEX1.2val_strates)
rownames(biais_FLEX1.2val)[1] <- "ALL"
rownames(RMSE_FLEX1.2val)[1] <- "ALL"
rm(biais_flex1.2val, biais_FLEX1.2val_strates)
rm(RMSE_flex1.2val, RMSE_FLEX1.2val_strates)
##4 noeuds
biais_FLEX1.4val <- rbind(biais_flex1.4val, biais_FLEX1.4val_strates)
RMSE_FLEX1.4val <- rbind(RMSE_flex1.4val, RMSE_FLEX1.4val_strates)
rownames(biais_FLEX1.4val)[1] <- "ALL"
rownames(RMSE_FLEX1.4val)[1] <- "ALL"
rm(biais_flex1.4val, biais_FLEX1.4val_strates)
rm(RMSE_flex1.4val, RMSE_FLEX1.4val_strates)
##FLEX2
##2 noeuds
biais_FLEX2.2val <- rbind(biais_flex2.2val, biais_FLEX2.2val_strates)
RMSE_FLEX2.2val <- rbind(RMSE_flex2.2val, RMSE_FLEX2.2val_strates)
rownames(biais_FLEX2.2val)[1] <- "ALL"
rownames(RMSE_FLEX2.2val)[1] <- "ALL"
rm(biais_flex2.2val, biais_FLEX2.2val_strates)
rm(RMSE_flex2.2val, RMSE_FLEX2.2val_strates)
##4 noeuds
biais_FLEX2.4val <- rbind(biais_flex2.4val, biais_FLEX2.4val_strates)
RMSE_FLEX2.4val <- rbind(RMSE_flex2.4val, RMSE_FLEX2.4val_strates)
rownames(biais_FLEX2.4val)[1] <- "ALL"
rownames(RMSE_FLEX2.4val)[1] <- "ALL"
rm(biais_flex2.4val, biais_FLEX2.4val_strates)
rm(RMSE_flex2.4val, RMSE_FLEX2.4val_strates)
##WG

biais_WGval <- rbind(biais_WGval, biais_WGval_strates)
RMSE_WGval <- rbind(RMSE_WGval, RMSE_WGval_strates)
rownames(biais_WGval)[1] <- "ALL"
rownames(RMSE_WGval)[1] <- "ALL"
rm(biais_WGval_strates)
rm(RMSE_WGval_strates)

##PP

biais_PPval <- rbind(biais_PPval, biais_PPval_strates)
RMSE_PPval <- rbind(RMSE_PPval, RMSE_PPval_strates)
rownames(biais_PPval)[1] <- "ALL"
rownames(RMSE_PPval)[1] <- "ALL"
rm(biais_PPval_strates)
rm(RMSE_PPval_strates)

### les dataframes des moyennes sont enregistrées sous forme de liste par rapport au modèle utilisé
########
##TRAIN
THEO_means_names <- c()
for(j in strata_names){
  THEO_means_names <- c(THEO_means_names, paste0("ALL_",j,"_theo"))
}
THEO_means <- c(mget("ALL_theo"), mget(THEO_means_names))
rm(ALL_theo)
rm(THEO_means_names)
for(j in strata_names){
  rm(list = paste0("ALL_",j,"_theo") )
}

##PLANN
PLANN_means_names <- c()
for(j in strata_names){
  PLANN_means_names <- c(PLANN_means_names, paste0("ALL_PLANN_",j))
}
PLANN_means <- c(mget("ALL_PLANN"), mget(PLANN_means_names))
rm(ALL_PLANN)
rm(PLANN_means_names)
for(j in strata_names){
  rm(list = paste0("ALL_PLANN_",j) )
}

##FLEX1
##2 noeuds
FLEX1.2_means_names <- c()
for(j in strata_names){
  FLEX1.2_means_names <- c(FLEX1.2_means_names, paste0("ALL_FLEX1.2_",j))
}
FLEX1.2_means <- c(mget("ALL_flex1.2"), mget(FLEX1.2_means_names))
rm(ALL_flex1.2)
rm(FLEX1.2_means_names)
for(j in strata_names){
  rm(list = paste0("ALL_FLEX1.2_",j) )
}

##4 noeuds
FLEX1.4_means_names <- c()
for(j in strata_names){
  FLEX1.4_means_names <- c(FLEX1.4_means_names, paste0("ALL_FLEX1.4_",j))
}
FLEX1.4_means <- c(mget("ALL_flex1.4"), mget(FLEX1.4_means_names))
rm(ALL_flex1.4)
rm(FLEX1.4_means_names)
for(j in strata_names){
  rm(list = paste0("ALL_FLEX1.4_",j) )
}

##FLEX2
##2 noeuds
FLEX2.2_means_names <- c()
for(j in strata_names){
  FLEX2.2_means_names <- c(FLEX2.2_means_names, paste0("ALL_FLEX2.2_",j))
}
FLEX2.2_means <- c(mget("ALL_flex2.2"), mget(FLEX2.2_means_names))
rm(ALL_flex2.2)
rm(FLEX2.2_means_names)
for(j in strata_names){
  rm(list = paste0("ALL_FLEX2.2_",j) )
}
##4 noeuds
FLEX2.4_means_names <- c()
for(j in strata_names){
  FLEX2.4_means_names <- c(FLEX2.4_means_names, paste0("ALL_FLEX2.4_",j))
}
FLEX2.4_means <- c(mget("ALL_flex2.4"), mget(FLEX2.4_means_names))
rm(ALL_flex2.4)
rm(FLEX2.4_means_names)
for(j in strata_names){
  rm(list = paste0("ALL_FLEX2.4_",j) )
}
##WG
WG_means_names <- c()
for(j in strata_names){
  WG_means_names <- c(WG_means_names, paste0("ALL_WG_",j))
}
WG_means <- c(mget("ALL_WG"), mget(WG_means_names))
rm(ALL_WG)
rm(WG_means_names)
for(j in strata_names){
  rm(list = paste0("ALL_WG_",j) )
}

PP_means_names <- c()
for(j in strata_names){
  PP_means_names <- c(PP_means_names, paste0("ALL_PP_",j))
}
PP_means <- c(mget("ALL_PP"), mget(PP_means_names))
rm(ALL_PP)
rm(PP_means_names)
for(j in strata_names){
  rm(list = paste0("ALL_PP_",j) )
}
############VALID
THEOval_means_names <- c()
for(j in strata_names){
  THEOval_means_names <- c(THEOval_means_names, paste0("ALL_",j,"_theoval"))
}
THEOval_means <- c(mget("ALL_theoval"), mget(THEOval_means_names))
rm(ALL_theoval)
rm(THEOval_means_names)
for(j in strata_names){
  rm(list = paste0("ALL_",j,"_theoval") )
}

##PLANN
PLANNval_means_names <- c()
for(j in strata_names){
  PLANNval_means_names <- c(PLANNval_means_names, paste0("ALL_PLANNval_",j))
}
PLANNval_means <- c(mget("ALL_PLANNval"), mget(PLANNval_means_names))
rm(ALL_PLANNval)
rm(PLANNval_means_names)
for(j in strata_names){
  rm(list = paste0("ALL_PLANNval_",j) )
}

##FLEX1
##2 noeuds
FLEX1.2val_means_names <- c()
for(j in strata_names){
  FLEX1.2val_means_names <- c(FLEX1.2val_means_names, paste0("ALL_FLEX1.2val_",j))
}
FLEX1.2val_means <- c(mget("ALL_flex1.2val"), mget(FLEX1.2val_means_names))
rm(ALL_flex1.2val)
rm(FLEX1.2val_means_names)
for(j in strata_names){
  rm(list = paste0("ALL_FLEX1.2val_",j) )
}

##4 noeuds
FLEX1.4val_means_names <- c()
for(j in strata_names){
  FLEX1.4val_means_names <- c(FLEX1.4val_means_names, paste0("ALL_FLEX1.4val_",j))
}
FLEX1.4val_means <- c(mget("ALL_flex1.4val"), mget(FLEX1.4val_means_names))
rm(ALL_flex1.4val)
rm(FLEX1.4val_means_names)
for(j in strata_names){
  rm(list = paste0("ALL_FLEX1.4val_",j) )
}
##FLEX2
##2 noeuds

FLEX2.2val_means_names <- c()
for(j in strata_names){
  FLEX2.2val_means_names <- c(FLEX2.2val_means_names, paste0("ALL_FLEX2.2val_",j))
}
FLEX2.2val_means <- c(mget("ALL_flex2.2val"), mget(FLEX2.2val_means_names))
rm(ALL_flex2.2val)
rm(FLEX2.2val_means_names)
for(j in strata_names){
  rm(list = paste0("ALL_FLEX2.2val_",j) )
}

##4 noeuds

FLEX2.4val_means_names <- c()
for(j in strata_names){
  FLEX2.4val_means_names <- c(FLEX2.4val_means_names, paste0("ALL_FLEX2.4val_",j))
}
FLEX2.4val_means <- c(mget("ALL_flex2.4val"), mget(FLEX2.4val_means_names))
rm(ALL_flex2.4val)
rm(FLEX2.4val_means_names)
for(j in strata_names){
  rm(list = paste0("ALL_FLEX2.4val_",j) )
}
##WG
WGval_means_names <- c()
for(j in strata_names){
  WGval_means_names <- c(WGval_means_names, paste0("ALL_WGval_",j))
}
WGval_means <- c(mget("ALL_WGval"), mget(WGval_means_names))
rm(ALL_WGval)
rm(WGval_means_names)
for(j in strata_names){
  rm(list = paste0("ALL_WGval_",j) )
}

##PP
PPval_means_names <- c()
for(j in strata_names){
  PPval_means_names <- c(PPval_means_names, paste0("ALL_PPval_",j))
}
PPval_means <- c(mget("ALL_PPval"), mget(PPval_means_names))
rm(ALL_PPval)
rm(PPval_means_names)
for(j in strata_names){
  rm(list = paste0("ALL_PPval_",j) )
}
########
##     5000ROCSTA
################################################### ROC.net ##############################################################

library(RISCA)
##choix de l'espacement des cutoffs 
# c.off <- seq(0,1, by = 0.2)
# sort(unique(1-ind_estimP[, l + 1]))
# unique( quantile( 1-ind_estimP[, l + 1],  probs= seq(0,1, by = 0.01)) ) #pas en dessous de 20 points 0.025
####TRAIN
### WHOLE

##importation des bases de données / survies individuelles prédites par les modèles 
for(k in iterations){
  ##bases de données 
  assign(paste0("DATATrain_", k) , read.csv(paste0(path0, "DATAFRAMES/TRAIN/WHOLE/",k,"_dftrain.csv"),
                                            sep = ";") )
  
  ##survies individuelles 
  #plann
  assign(paste0("indTrainPLANN_", k) , read.csv(paste0(path0, "TRAIN/WHOLE/PLANN/individual_survival/",k,"_ind_PLANN.csv"),
                                                sep = ";") )
  #flex1
  assign(paste0("indTrainFLEX1.2_", k) , read.csv(paste0(path0, "TRAIN/WHOLE/FLEX1.2/individual_survival/",k,"_ind_flex12.csv"),
                                                  sep = ";") )
  assign(paste0("indTrainFLEX1.4_", k) , read.csv(paste0(path0, "TRAIN/WHOLE/FLEX1.4/individual_survival/",k,"_ind_flex14.csv"),
                                                  sep = ";") )
  #flex2
  assign(paste0("indTrainFLEX2.2_", k) , read.csv(paste0(path0, "TRAIN/WHOLE/FLEX2.2/individual_survival/",k,"_ind_flex22.csv"),
                                                  sep = ";") )
  assign(paste0("indTrainFLEX2.4_", k) , read.csv(paste0(path0, "TRAIN/WHOLE/FLEX2.4/individual_survival/",k,"_ind_flex24.csv"),
                                                  sep = ";") )
  #WG
  assign(paste0("indTrainWG_", k) , read.csv(paste0(path0, "TRAIN/WHOLE/WG/individual_survival/",k,"_ind_WG.csv"),
                                             sep = ";") )
}

##### calcul de ROC net

AUC_calc_W <- function(k) {
  datatrain <- get(paste0("DATATrain_", k))
  ind_estimP <- get(paste0("indTrainPLANN_", k))
  ind_estimF1.2 <- get(paste0("indTrainFLEX1.2_", k))
  ind_estimF1.4 <- get(paste0("indTrainFLEX1.4_", k))
  ind_estimF2.2 <- get(paste0("indTrainFLEX2.2_", k))
  ind_estimF2.4 <- get(paste0("indTrainFLEX2.4_", k))
  ind_estimWG <- get(paste0("indTrainWG_", k))
  
  rm(list = c(paste0("DATATrain_",k), paste0("indTrainPLANN_",k),
              paste0("indTrainFLEX1.2_",k), paste0("indTrainFLEX1.4_",k),
              paste0("indTrainFLEX2.2_",k), paste0("indTrainFLEX2.4_",k),
              paste0("indTrainWG_",k)))
  
  # Initialiser les vecteurs de résultats
  hold_P <- hold_F1.2 <- hold_F1.4 <- hold_F2.2 <- hold_F2.4 <- hold_WG <- c()
  # sort(unique(1-ind_estimP[, l + 1]))
  # unique( quantile( 1-ind_estimP[, l + 1],  probs= seq(0,1, by = 0.01)) ) #pas en dessous de 20 points 0.025
  
  for (l in seq_along(newtimes)) {
    hold_P    <- c(hold_P,    roc.net(datatrain$times, datatrain$status, 1-ind_estimP[,l + 1], datatrain$age, datatrain$sexchara,
                                      datatrain$year, slopop, pro.time = newtimes[l], 
                                      cut.off = unique( quantile(1-ind_estimP[,l+1],probs= seq(0,1, by = 0.01)) ) )$auc)
    hold_F1.2 <- c(hold_F1.2, roc.net(datatrain$times, datatrain$status, 1-ind_estimF1.2[, l], datatrain$age, datatrain$sexchara,
                                      datatrain$year, slopop, pro.time = newtimes[l],
                                      cut.off = unique( quantile( 1-ind_estimF1.2[, l],  probs= seq(0,1, by = 0.01)) ))$auc)
    hold_F1.4 <- c(hold_F1.4, roc.net(datatrain$times, datatrain$status, 1-ind_estimF1.4[, l], datatrain$age, datatrain$sexchara,
                                      datatrain$year, slopop, pro.time = newtimes[l],
                                      cut.off = unique( quantile( 1-ind_estimF1.4[, l],  probs= seq(0,1, by = 0.01)) ))$auc)
    hold_F2.2 <- c(hold_F2.2, roc.net(datatrain$times, datatrain$status, 1-ind_estimF2.2[, l], datatrain$age, datatrain$sexchara,
                                      datatrain$year, slopop, pro.time = newtimes[l], 
                                      cut.off = unique( quantile( 1-ind_estimF2.2[, l],  probs= seq(0,1, by = 0.01)) ))$auc)
    hold_F2.4 <- c(hold_F2.4, roc.net(datatrain$times, datatrain$status, 1-ind_estimF2.4[, l], datatrain$age, datatrain$sexchara,
                                      datatrain$year, slopop, pro.time = newtimes[l],
                                      cut.off = unique( quantile( 1-ind_estimF2.4[, l],  probs= seq(0,1, by = 0.01)) ))$auc)
    hold_WG   <- c(hold_WG,   roc.net(datatrain$times, datatrain$status, 1-ind_estimWG[, l],   datatrain$age, datatrain$sexchara,
                                      datatrain$year, slopop, pro.time = newtimes[l],
                                      cut.off = unique( quantile( 1-ind_estimWG[, l],  probs= seq(0,1, by = 0.01)) ))$auc)
  }
  
  return(list(
    P     = hold_P,
    F1.2  = hold_F1.2,
    F1.4  = hold_F1.4,
    F2.2  = hold_F2.2,
    F2.4  = hold_F2.4,
    WG    = hold_WG
  ))
}

results_list <- mclapply(iterations, AUC_calc_W, mc.cores = detectCores() - 4)

#supression des bases inutiles 
for(k in iterations){
  rm(list = c(paste0("DATATrain_",k), paste0("indTrainPLANN_",k),
              paste0("indTrainFLEX1.2_",k), paste0("indTrainFLEX1.4_",k),
              paste0("indTrainFLEX2.2_",k), paste0("indTrainFLEX2.4_",k),
              paste0("indTrainWG_",k)))
}
# Transformer en matrices pour chaque méthode :
AUC_WHOLE_P    <- do.call(rbind, lapply(results_list, function(x) x$P))
AUC_WHOLE_F1.2 <- do.call(rbind, lapply(results_list, function(x) x$F1.2))
AUC_WHOLE_F1.4 <- do.call(rbind, lapply(results_list, function(x) x$F1.4))
AUC_WHOLE_F2.2 <- do.call(rbind, lapply(results_list, function(x) x$F2.2))
AUC_WHOLE_F2.4 <- do.call(rbind, lapply(results_list, function(x) x$F2.4))
AUC_WHOLE_WG   <- do.call(rbind, lapply(results_list, function(x) x$WG))


### pour renommer les colonnes avec les temps de prognostic
for(l in c("P","F1.2","F1.4","F2.2","F2.4","WG")){
  dathold <- get(paste0("AUC_WHOLE_", l))
  colnames(dathold) <- c("1 year", "3 years", "5 years", "10 years")
  rownames(dathold) <- c(iterations)
  assign(paste0("AUC_WHOLE_", l),dathold)
}
## moyennes par temps 

for(l in c("P","F1.2","F1.4","F2.2","F2.4","WG")){
  
  assign(paste0("AUC_means_",l) ,  colMeans( get(paste0("AUC_WHOLE_",l)) ) )
  
}
################## STRATES

##importation des bases de données / survies individuelles prédites par les modèles 

for(k in iterations){
  ##bases de données 
  for (j in strata_names){
    assign(paste0("DATATrain_",j,"_", k) , read.csv(paste0(path0, "DATAFRAMES/TRAIN/",j,"/",k,"_dftrain_",j,".csv"),
                                                    sep = ";") )
    ##survies individuelles 
    #plann
    assign(paste0("indTrainPLANN_",j,"_", k) , read.csv(paste0(path0, "TRAIN/STRATA/PLANN/individual_survival/",k,"_ind_PLANN_",j,".csv"),
                                                        sep = ";") )
    #flex1
    assign(paste0("indTrainFLEX1.2_",j,"_", k) , read.csv(paste0(path0, "TRAIN/STRATA/FLEX1.2/individual_survival/",k,"_ind_flex12_",j,".csv"),
                                                          sep = ";") )
    assign(paste0("indTrainFLEX1.4_",j,"_", k) , read.csv(paste0(path0, "TRAIN/STRATA/FLEX1.4/individual_survival/",k,"_ind_flex14_",j,".csv"),
                                                          sep = ";") )
    #flex2
    assign(paste0("indTrainFLEX2.2_",j,"_", k) , read.csv(paste0(path0, "TRAIN/STRATA/FLEX2.2/individual_survival/",k,"_ind_flex22_",j,".csv"),
                                                          sep = ";") )
    assign(paste0("indTrainFLEX2.4_",j,"_", k) , read.csv(paste0(path0, "TRAIN/STRATA/FLEX2.4/individual_survival/",k,"_ind_flex24_",j,".csv"),
                                                          sep = ";") )
    #WG
    assign(paste0("indTrainWG_",j,"_", k) , read.csv(paste0(path0, "TRAIN/STRATA/WG/individual_survival/",k,"_ind_WG_",j,".csv"),
                                                     sep = ";") )
  }
}

##### calcul de ROC net

AUC_calc_S <- function(k) {
  auc_results <- list()
  
  for (j in strata_names) {
    datatrain <- get(paste0("DATATrain_", j, "_", k))
    ind_estimP <- get(paste0("indTrainPLANN_", j, "_", k))
    ind_estimF1.2 <- get(paste0("indTrainFLEX1.2_", j, "_", k))
    ind_estimF1.4 <- get(paste0("indTrainFLEX1.4_", j, "_", k))
    ind_estimF2.2 <- get(paste0("indTrainFLEX2.2_", j, "_", k))
    ind_estimF2.4 <- get(paste0("indTrainFLEX2.4_", j, "_", k))
    ind_estimWG <- get(paste0("indTrainWG_", j, "_", k))
    
    hold_P <- hold_F1.2 <- hold_F1.4 <- hold_F2.2 <- hold_F2.4 <- hold_WG <- c()
    
    for (l in seq_along(newtimes)) {
      hold_P    <- c(hold_P,    roc.net(datatrain$times, datatrain$status, 1-ind_estimP[,l + 1], datatrain$age, datatrain$sexchara,
                                        datatrain$year, slopop, pro.time = newtimes[l], 
                                        cut.off = unique( quantile(1-ind_estimP[,l+1],probs= seq(0,1, by = 0.01)) ) )$auc)
      hold_F1.2 <- c(hold_F1.2, roc.net(datatrain$times, datatrain$status, 1-ind_estimF1.2[, l], datatrain$age, datatrain$sexchara,
                                        datatrain$year, slopop, pro.time = newtimes[l],
                                        cut.off = unique( quantile( 1-ind_estimF1.2[, l],  probs= seq(0,1, by = 0.01)) ))$auc)
      hold_F1.4 <- c(hold_F1.4, roc.net(datatrain$times, datatrain$status, 1-ind_estimF1.4[, l], datatrain$age, datatrain$sexchara,
                                        datatrain$year, slopop, pro.time = newtimes[l],
                                        cut.off = unique( quantile( 1-ind_estimF1.4[, l],  probs= seq(0,1, by = 0.01)) ))$auc)
      hold_F2.2 <- c(hold_F2.2, roc.net(datatrain$times, datatrain$status, 1-ind_estimF2.2[, l], datatrain$age, datatrain$sexchara,
                                        datatrain$year, slopop, pro.time = newtimes[l], 
                                        cut.off = unique( quantile( 1-ind_estimF2.2[, l],  probs= seq(0,1, by = 0.01)) ))$auc)
      hold_F2.4 <- c(hold_F2.4, roc.net(datatrain$times, datatrain$status, 1-ind_estimF2.4[, l], datatrain$age, datatrain$sexchara,
                                        datatrain$year, slopop, pro.time = newtimes[l],
                                        cut.off = unique( quantile( 1-ind_estimF2.4[, l],  probs= seq(0,1, by = 0.01)) ))$auc)
      hold_WG   <- c(hold_WG,   roc.net(datatrain$times, datatrain$status, 1-ind_estimWG[, l],   datatrain$age, datatrain$sexchara,
                                        datatrain$year, slopop, pro.time = newtimes[l],
                                        cut.off = unique( quantile( 1-ind_estimWG[, l],  probs= seq(0,1, by = 0.01)) ))$auc)
    }
    
    auc_results[[j]] <- list(
      k = k,
      P = hold_P,
      F1.2 = hold_F1.2,
      F1.4 = hold_F1.4,
      F2.2 = hold_F2.2,
      F2.4 = hold_F2.4,
      WG = hold_WG
    )
  }
  
  return(auc_results)
}

res_all <- mclapply(iterations, AUC_calc_S, mc.cores = detectCores() - 4)

##supression des bases inutiles pour faire de la place
for(k in iterations){
  for(j in strata_names){
    rm(list = c(
      paste0("DATATrain_", j, "_", k),
      paste0("indTrainPLANN_", j, "_", k),
      paste0("indTrainFLEX1.2_", j, "_", k),
      paste0("indTrainFLEX1.4_", j, "_", k),
      paste0("indTrainFLEX2.2_", j, "_", k),
      paste0("indTrainFLEX2.4_", j, "_", k),
      paste0("indTrainWG_", j, "_", k)
    ))
  }
}

# Transformer en matrices pour chaque méthode :

for(j in strata_names){
  assign(paste0("AUC_",j,"_P"), do.call(rbind, lapply(res_all, function(x) x[[j]]$P)) )
  assign(paste0("AUC_",j,"_F1.2"), do.call(rbind, lapply(res_all, function(x) x[[j]]$F1.2)) )
  assign(paste0("AUC_",j,"_F1.4"), do.call(rbind, lapply(res_all, function(x) x[[j]]$F1.4)) )
  assign(paste0("AUC_",j,"_F2.2"), do.call(rbind, lapply(res_all, function(x) x[[j]]$F2.2)) )
  assign(paste0("AUC_",j,"_F2.4"), do.call(rbind, lapply(res_all, function(x) x[[j]]$F2.4)) )
  assign(paste0("AUC_",j,"_WG"), do.call(rbind, lapply(res_all, function(x) x[[j]]$WG)) )
}
### pour renommer les colonnes avec les temps de prognostic
for(l in c("P","F1.2","F1.4","F2.2","F2.4","WG")){
  for(j in strata_names){
    dathold <- get(paste0("AUC_", j,"_",l))
    colnames(dathold) <- c("1 year", "3 years", "5 years", "10 years")
    rownames(dathold) <- c(iterations)
    assign(paste0("AUC_",j,"_", l),dathold)
  }
}
## moyennes par temps 

for(l in c("P","F1.2","F1.4","F2.2","F2.4","WG")){
  for(j in strata_names){
    assign(paste0("AUC_means_",j,"_",l) ,  colMeans( get(paste0("AUC_",j,"_",l)) ) )
  }
}


### sous forme de liste pour libérer de l'espace dans l'evt (pour les moyennes)

ROCmean_results <- list()

for(l in c("P","F1.2","F1.4","F2.2","F2.4","WG")){
  ROCmean_results[['WHOLE']][[l]] <- list()
  ROCmean_results[['WHOLE']][[l]] <- get( paste0("AUC_means_",l) )
  rm(list = paste0("AUC_means_",l))
  
  for(j in strata_names){
    ROCmean_results[[j]][[l]] <- get(paste0("AUC_means_", j, "_", l)) 
    rm(list = paste0("AUC_means_", j, "_", l))
  }
  
}
### meme chose pour les valeurs pour chaque itération 


ROC_results <- list()

for(l in c("P","F1.2","F1.4","F2.2","F2.4","WG")){
  ROC_results[['WHOLE']][[l]] <- list()
  ROC_results[['WHOLE']][[l]] <- get( paste0("AUC_WHOLE_",l) )
  rm(list = paste0("AUC_WHOLE_",l) )
  for(j in strata_names){
    ROC_results[[j]][[l]] <- get(paste0("AUC_", j, "_", l)) 
    rm(list = paste0("AUC_", j, "_", l))
  }
  
}
#################### VALIDATION ######################

### WHOLE

##importation des bases de données / survies individuelles prédites par les modèles 
for(k in iterations){
  ##bases de données 
  assign(paste0("DATAValid_", k) , read.csv(paste0(path0, "DATAFRAMES/VALID/WHOLE/",k,"_dfvalid.csv"),
                                            sep = ";") )
  
  ##survies individuelles 
  #plann
  assign(paste0("indValidPLANN_", k) , read.csv(paste0(path0, "VALID/WHOLE/PLANN/individual_survival/",k,"_ind_PLANNval.csv"),
                                                sep = ";") )
  
  #flex1
  assign(paste0("indValidFLEX1.2_",k) , read.csv(paste0(path0, "VALID/WHOLE/FLEX1.2/individual_survival/",k,"_ind_flexval12.csv"),
                                                 sep = ";") )
  assign(paste0("indValidFLEX1.4_",k) , read.csv(paste0(path0, "VALID/WHOLE/FLEX1.4/individual_survival/",k,"_ind_flexval14.csv"),
                                                 sep = ";") )
  #flex2
  assign(paste0("indValidFLEX2.2_",k) , read.csv(paste0(path0, "VALID/WHOLE/FLEX2.2/individual_survival/",k,"_ind_flexval22.csv"),
                                                 sep = ";") )
  assign(paste0("indValidFLEX2.4_",k) , read.csv(paste0(path0, "VALID/WHOLE/FLEX2.4/individual_survival/",k,"_ind_flexval24.csv"),
                                                 sep = ";") )
  #WG
  assign(paste0("indValidWG_", k) , read.csv(paste0(path0, "VALID/WHOLE/WG/individual_survival/",k,"_ind_flexWGval.csv"),
                                             sep = ";") )
}

##### calcul de ROC net

##creation d'objets pour enregistrer les AUC de chaque itération aux différents temps et pour les différents modèles 

AUC_calc_val_W <- function(k){ 
  
  dataval <- get(paste0("DATAValid_",k))
  ind_estimP <- get(paste0("indValidPLANN_",k))
  ind_estimF1.2 <- get(paste0("indValidFLEX1.2_",k))
  ind_estimF1.4 <- get(paste0("indValidFLEX1.4_", k))
  ind_estimF2.2 <- get(paste0("indValidFLEX2.2_", k))
  ind_estimF2.4 <- get(paste0("indValidFLEX2.4_", k))
  ind_estimWG <- get(paste0("indValidWG_",k))
  
  rm(list = c(paste0("DATAValid_",k), paste0("indValidPLANN_",k),paste0("indValidFLEX1.2_",k), paste0("indValidFLEX1.4_", k),
              paste0("indValidFLEX2.2_", k),paste0("indValidFLEX2.4_", k),paste0("indValidWG_",k)))
  
  # Initialiser les vecteurs de résultats
  hold_P <- hold_F1.2 <- hold_F1.4 <- hold_F2.2 <- hold_F2.4 <- hold_WG <- c()
  
  ###PLANN
  for (l in seq_along(newtimes)) {
    hold_P    <- c(hold_P,    roc.net(dataval$times, dataval$status, 1-ind_estimP[,l + 1], dataval$age, dataval$sexchara,
                                      dataval$year, slopop, pro.time = newtimes[l], 
                                      cut.off = unique( quantile(1-ind_estimP[,l+1],probs= seq(0,1, by = 0.01)) ) )$auc)
    hold_F1.2 <- c(hold_F1.2, roc.net(dataval$times, dataval$status, 1-ind_estimF1.2[, l], dataval$age, dataval$sexchara,
                                      dataval$year, slopop, pro.time = newtimes[l],
                                      cut.off = unique( quantile( 1-ind_estimF1.2[, l],  probs= seq(0,1, by = 0.01)) ))$auc)
    hold_F1.4 <- c(hold_F1.4, roc.net(dataval$times, dataval$status, 1-ind_estimF1.4[, l], dataval$age, dataval$sexchara,
                                      dataval$year, slopop, pro.time = newtimes[l],
                                      cut.off = unique( quantile( 1-ind_estimF1.4[, l],  probs= seq(0,1, by = 0.01)) ))$auc)
    hold_F2.2 <- c(hold_F2.2, roc.net(dataval$times, dataval$status, 1-ind_estimF2.2[, l], dataval$age, dataval$sexchara,
                                      dataval$year, slopop, pro.time = newtimes[l], 
                                      cut.off = unique( quantile( 1-ind_estimF2.2[, l],  probs= seq(0,1, by = 0.01)) ))$auc)
    hold_F2.4 <- c(hold_F2.4, roc.net(dataval$times, dataval$status, 1-ind_estimF2.4[, l], dataval$age, dataval$sexchara,
                                      dataval$year, slopop, pro.time = newtimes[l],
                                      cut.off = unique( quantile( 1-ind_estimF2.4[, l],  probs= seq(0,1, by = 0.01)) ))$auc)
    hold_WG   <- c(hold_WG,   roc.net(dataval$times, dataval$status, 1-ind_estimWG[, l],   dataval$age, dataval$sexchara,
                                      dataval$year, slopop, pro.time = newtimes[l],
                                      cut.off = unique( quantile( 1-ind_estimWG[, l],  probs= seq(0,1, by = 0.01)) ))$auc)
  }
  
  return(list(
    P     = hold_P,
    F1.2  = hold_F1.2,
    F1.4  = hold_F1.4,
    F2.2  = hold_F2.2,
    F2.4  = hold_F2.4,
    WG    = hold_WG
  ))
}

resultval_list <- mclapply(iterations, AUC_calc_val_W, mc.cores = detectCores() - 4)

#supression des bases inutiles 
for(k in iterations){
  rm(list = c(paste0("DATAValid_",k), paste0("indValidPLANN_",k),
              paste0("indValidFLEX1.2_",k), paste0("indValidFLEX1.4_",k),
              paste0("indValidFLEX2.2_",k), paste0("indValidFLEX2.4_",k),
              paste0("indValidWG_",k)))
}
# Transformer en matrices pour chaque méthode :
AUCval_WHOLE_P    <- do.call(rbind, lapply(resultval_list, function(x) x$P))
AUCval_WHOLE_F1.2 <- do.call(rbind, lapply(resultval_list, function(x) x$F1.2))
AUCval_WHOLE_F1.4 <- do.call(rbind, lapply(resultval_list, function(x) x$F1.4))
AUCval_WHOLE_F2.2 <- do.call(rbind, lapply(resultval_list, function(x) x$F2.2))
AUCval_WHOLE_F2.4 <- do.call(rbind, lapply(resultval_list, function(x) x$F2.4))
AUCval_WHOLE_WG   <- do.call(rbind, lapply(resultval_list, function(x) x$WG))

### pour renommer les colonnes avec les temps de prognostic
for(l in c("P","F1.2","F1.4","F2.2","F2.4","WG")){
  dathold <- get(paste0("AUCval_WHOLE_", l))
  colnames(dathold) <- c("1 year", "3 years", "5 years", "10 years")
  rownames(dathold) <- c(iterations)
  assign(paste0("AUCval_WHOLE_", l),dathold)
}
## moyennes par temps 

for(l in c("P","F1.2","F1.4","F2.2","F2.4","WG")){
  
  assign(paste0("AUCval_means_",l) ,  colMeans( get(paste0("AUCval_WHOLE_",l)) ) )
  
}
################## STRATES

##importation des bases de données / survies individuelles prédites par les modèles 

for(k in iterations){
  ##bases de données 
  for (j in strata_names){
    assign(paste0("DATAValid_",j,"_", k) , read.csv(paste0(path0, "DATAFRAMES/VALID/",j,"/",k,"_dfvalid_",j,".csv"),
                                                    sep = ";") )
    ##survies individuelles 
    #plann
    assign(paste0("indValidPLANN_",j,"_", k) , read.csv(paste0(path0, "VALID/STRATA/PLANN/individual_survival/",k,"_ind_PLANNval_",j,".csv"),
                                                        sep = ";") )
    #flex1
    assign(paste0("indValidFLEX1.2_",j,"_", k) , read.csv(paste0(path0, "VALID/STRATA/FLEX1.2/individual_survival/",k,"_ind_flexval12_",j,".csv"),
                                                          sep = ";") )
    assign(paste0("indValidFLEX1.4_",j,"_", k) , read.csv(paste0(path0, "VALID/STRATA/FLEX1.4/individual_survival/",k,"_ind_flexval14_",j,".csv"),
                                                          sep = ";") )
    #flex2
    assign(paste0("indValidFLEX2.2_",j,"_", k) , read.csv(paste0(path0, "VALID/STRATA/FLEX2.2/individual_survival/",k,"_ind_flexval22_",j,".csv"),
                                                          sep = ";") )
    assign(paste0("indValidFLEX2.4_",j,"_", k) , read.csv(paste0(path0, "VALID/STRATA/FLEX2.4/individual_survival/",k,"_ind_flexval24_",j,".csv"),
                                                          sep = ";") )
    #WG
    assign(paste0("indValidWG_",j,"_", k) , read.csv(paste0(path0, "VALID/STRATA/WG/individual_survival/",k,"_ind_WGval_",j,".csv"),
                                                     sep = ";") )
  }
}

##### calcul de ROC net

AUC_calc_val_S <- function(k){ 
  auc_results <- list()
  for(j in strata_names){
    dataval <- get(paste0("DATAValid_",j,"_", k))
    ind_estimP <- get(paste0("indValidPLANN_",j,"_", k))
    ind_estimF1.2 <- get(paste0("indValidFLEX1.2_",j,"_", k))
    ind_estimF1.4 <- get(paste0("indValidFLEX1.4_",j,"_", k))
    ind_estimF2.2 <- get(paste0("indValidFLEX2.2_",j,"_", k))
    ind_estimF2.4 <- get(paste0("indValidFLEX2.4_",j,"_", k))
    ind_estimWG <- get(paste0("indValidWG_",j,"_", k))
    
    rm(list = c(paste0("DATAValid_",j,"_", k), paste0("indValidPLANN_",j,"_", k), paste0("indValidFLEX1.2_",j,"_", k),paste0("indValidFLEX1.4_",j,"_", k),
                paste0("indValidFLEX2.2_",j,"_", k), paste0("indValidFLEX2.2_",j,"_", k), paste0("indValidWG_",j,"_", k)))
    
    hold_P <- hold_F1.2 <- hold_F1.4 <- hold_F2.2 <- hold_F2.4 <- hold_WG <- c()
    
    for (l in seq_along(newtimes)) {
      hold_P    <- c(hold_P,    roc.net(dataval$times, dataval$status, 1-ind_estimP[,l + 1], dataval$age, dataval$sexchara,
                                        dataval$year, slopop, pro.time = newtimes[l], 
                                        cut.off = unique( quantile(1-ind_estimP[,l+1],probs= seq(0,1, by = 0.01)) ) )$auc)
      hold_F1.2 <- c(hold_F1.2, roc.net(dataval$times, dataval$status, 1-ind_estimF1.2[, l], dataval$age, dataval$sexchara,
                                        dataval$year, slopop, pro.time = newtimes[l],
                                        cut.off = unique( quantile( 1-ind_estimF1.2[, l],  probs= seq(0,1, by = 0.01)) ))$auc)
      hold_F1.4 <- c(hold_F1.4, roc.net(dataval$times, dataval$status, 1-ind_estimF1.4[, l], dataval$age, dataval$sexchara,
                                        dataval$year, slopop, pro.time = newtimes[l],
                                        cut.off = unique( quantile( 1-ind_estimF1.4[, l],  probs= seq(0,1, by = 0.01)) ))$auc)
      hold_F2.2 <- c(hold_F2.2, roc.net(dataval$times, dataval$status, 1-ind_estimF2.2[, l], dataval$age, dataval$sexchara,
                                        dataval$year, slopop, pro.time = newtimes[l], 
                                        cut.off = unique( quantile( 1-ind_estimF2.2[, l],  probs= seq(0,1, by = 0.01)) ))$auc)
      hold_F2.4 <- c(hold_F2.4, roc.net(dataval$times, dataval$status, 1-ind_estimF2.4[, l], dataval$age, dataval$sexchara,
                                        dataval$year, slopop, pro.time = newtimes[l],
                                        cut.off = unique( quantile( 1-ind_estimF2.4[, l],  probs= seq(0,1, by = 0.01)) ))$auc)
      hold_WG   <- c(hold_WG,   roc.net(dataval$times, dataval$status, 1-ind_estimWG[, l],   dataval$age, dataval$sexchara,
                                        dataval$year, slopop, pro.time = newtimes[l],
                                        cut.off = unique( quantile( 1-ind_estimWG[, l],  probs= seq(0,1, by = 0.01)) ))$auc)
    }
    
    auc_results[[j]] <- list(
      k = k,
      P = hold_P,
      F1.2 = hold_F1.2,
      F1.4 = hold_F1.4,
      F2.2 = hold_F2.2,
      F2.4 = hold_F2.4,
      WG = hold_WG
    )
  }
  
  return(auc_results)
}

resval_all <- mclapply(iterations, AUC_calc_val_S, mc.cores = detectCores() - 4)

for(k in iterations){
  for(j in strata_names){
    rm(list = c(
      paste0("DATAValid_", j, "_", k),
      paste0("indValidPLANN_", j, "_", k),
      paste0("indValidFLEX1.2_", j, "_", k),
      paste0("indValidFLEX1.4_", j, "_", k),
      paste0("indValidFLEX2.2_", j, "_", k),
      paste0("indValidFLEX2.4_", j, "_", k),
      paste0("indValidWG_", j, "_", k)
    ))
  }
}

# Transformer en matrices pour chaque méthode :

for(j in strata_names){
  assign(paste0("AUCval_",j,"_P"), do.call(rbind, lapply(resval_all, function(x) x[[j]]$P)) )
  assign(paste0("AUCval_",j,"_F1.2"), do.call(rbind, lapply(resval_all, function(x) x[[j]]$F1.2)) )
  assign(paste0("AUCval_",j,"_F1.4"), do.call(rbind, lapply(resval_all, function(x) x[[j]]$F1.4)) )
  assign(paste0("AUCval_",j,"_F2.2"), do.call(rbind, lapply(resval_all, function(x) x[[j]]$F2.2)) )
  assign(paste0("AUCval_",j,"_F2.4"), do.call(rbind, lapply(resval_all, function(x) x[[j]]$F2.4)) )
  assign(paste0("AUCval_",j,"_WG"), do.call(rbind, lapply(resval_all, function(x) x[[j]]$WG)) )
}
### pour renommer les colonnes avec les temps de prognostic
for(l in c("P","F1.2","F1.4","F2.2","F2.4","WG")){
  for(j in strata_names){
    dathold <- get(paste0("AUCval_", j,"_",l))
    colnames(dathold) <- c("1 year", "3 years", "5 years", "10 years")
    rownames(dathold) <- c(iterations)
    assign(paste0("AUCval_",j,"_", l),dathold)
  }
}
## moyennes par temps 

for(l in c("P","F1.2","F1.4","F2.2","F2.4","WG")){
  for(j in strata_names){
    assign(paste0("AUCval_means_",j,"_",l) ,  colMeans( get(paste0("AUCval_",j,"_",l)) ) )
  }
}


### sous forme de liste pour libérer de l'espace dans l'evt (pour les moyennes)

ROCmeanval_results <- list()

for(l in c("P","F1.2","F1.4","F2.2","F2.4","WG")){
  ROCmeanval_results[['WHOLE']][[l]] <- list()
  ROCmeanval_results[['WHOLE']][[l]] <- get( paste0("AUCval_means_",l) )
  rm(list = paste0("AUCval_means_",l))
  
  for(j in strata_names){
    ROCmeanval_results[[j]][[l]] <- get(paste0("AUCval_means_", j, "_", l)) 
    rm(list = paste0("AUCval_means_", j, "_", l))
  }
  
}
### meme chose pour les valeurs pour chaque itération 


ROCval_results <- list()

for(l in c("P","F1.2","F1.4","F2.2","F2.4","WG")){
  ROCval_results[['WHOLE']][[l]] <- list()
  ROCval_results[['WHOLE']][[l]] <- get( paste0("AUCval_WHOLE_",l) )
  rm(list = paste0("AUCval_WHOLE_",l) )
  for(j in strata_names){
    ROCval_results[[j]][[l]] <- get(paste0("AUCval_", j, "_", l)) 
    rm(list = paste0("AUCval_", j, "_", l))
  }
  
}

print(Sys.time()-start)

save.image("~/Documents/Simulations/Simulations mai 2025/Résultats/1000ite_5000ind_NPH_HC_14_05.Rdata")
### 5000IndEND
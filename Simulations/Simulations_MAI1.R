# remotes::install_github("ollard-t/survivalNET")
# remotes::install_github("ollard-t/survivalPLANN")

# install.packages("~/Documents/GitHub/survivalNET", repos = NULL, type = "source")
# install.packages("~/Documents/GitHub/survivalPLANN", repos = NULL, type = "source")

##base deja générée pour ne pas a avoir à la refaire (reprendre à la ligne 172 : N/2 data_train)
# data <- read.csv(paste0("~/Documents/Rstudio/Simulations/Simulations février 2025/data_1.csv"),sep = ";")
#il faut juste repréciser N, les valeurs des theta sigma, etc et les fnctions Sn, SP et S_obs
library(survivalNET)
library(survivalPLANN)
library(parallel)
library(doParallel)

#################################################################################################
### path ###

path0 <- "~/Documents/Simulations/Simulations mai 2025/Output simulations/"
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
  
  data <- read.csv(file = paste0(path1,"ind1000/",i,"_df1000.csv"), sep = ";")
  Tsigma <- 12.6 
  Tnu <- -0.5
  Ttheta <- 0
  betaZ <- c(0.9,2.7,0.3,-0.1,-0.1)
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
  
  theopred <- t(sapply(1:dim(Zdata_train)[1], FUN = function(i){
    Sn(time = newtimes, sigma = Tsigma, nu = Tnu, theta = Ttheta, beta = betaZ, covariates = Zdata_train[i,]) } ) )
  
  colnames(theopred) <- newtimes
  
  mean.theo <- apply(theopred, FUN="mean", MARGIN=2)
  
  ###strates
  #####
  
  strata_names <- c("HC", "HR", "FC", "FR")
  
  for(j in strata_names){
    
    data_train_var <- get(paste0("data_train_", j))
    
    Zdata_train <- cbind(data_train_var$stage2, data_train_var$stage3, data_train_var$agey10, data_train_var$sex, data_train_var$colon)
    
    theopredstrat <- t(sapply(1:nrow(Zdata_train), FUN = function(i) {
      Sn(time = newtimes, sigma = Tsigma, nu = Tnu, theta = Ttheta, beta = betaZ, covariates = Zdata_train[i,])
    }))
    
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
    failedindic1.2 <- paste0(m1.2, " : this value was decreased from ", m1.2_original)
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
    failedindic1.4 <- paste0(m1.4, " : this value was decreased from ", m1.4_original)
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

# nb.cluster = detectCores() - 4   
# cl <- makeCluster(nb.cluster, type="FORK") 
# registerDoParallel(cl)
# 
# foreach(i = 1:iterations,  .inorder = FALSE, .verbose = T) %dopar% 
#   { simulate_iteration(i) }
# stopCluster(cl)
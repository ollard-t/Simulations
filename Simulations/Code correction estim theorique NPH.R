cens_val <- "HC"
PH_val <- "NPH"

file.exists("~/Documents/Rstudio/Simulations/BASES/")
for(N in c(1000,3000,5000)){
  path0 <- paste0("~/Documents/Rstudio/Simulations/Simulations mai 2025/Résultats/N",N,"_results/",PH_val,"/Output simulations_N",N,"_",PH_val,"_",cens_val,"/")
  print(file.exists(path0))
}
rm(path0)
rm(N)
#################################################################################################
#################################################################################################
#################################################################################################
##1000##
#################################################################################################
#################################################################################################
#################################################################################################

library(survivalNET)
library(survivalPLANN)
library(parallel)
library(doParallel)
date_launch <- Sys.Date()

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

iterations <- 1:1000

N = 1000
path0 <- paste0("~/Documents/Rstudio/Simulations/Simulations mai 2025/Résultats/N",N,"_results/",PH_val,"/Output simulations_N",N,"_",PH_val,"_",cens_val,"/")
path1 <- "~/Documents/Rstudio/Simulations/BASES/"

correc_theo<- function(i, N){
    
    set.seed(i)
    
    folder_suffix <- ""
    if (PH_val == "PH"  && cens_val == "LC") folder_suffix <- ""
    if (PH_val == "PH"  && cens_val == "HC") folder_suffix <- "_HC"
    if (PH_val == "NPH" && cens_val == "LC") folder_suffix <- "NPH"
    if (PH_val == "NPH" && cens_val == "HC") folder_suffix <- "NPH_HC"
    
    file_suffix <- ""
    if (PH_val == "PH" && cens_val == "HC") file_suffix <- "HC"
    
    
    path_t <- paste0(path1, "ind", N, folder_suffix, "/", i, "_df", N, file_suffix, ".csv")
    data <- read.csv(path_t, sep = ";")
    
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
    betaZFC <- c(0.6, 2.6, -0.1) # c(0.6,2.6,-0.1)
    
    #### Femmes Rectum
    TsigmaFR <- 8
    TnuFR <- 1.5
    TthetaFR <- 0
    betaZFR <- c(0.9,2.7,0.16)
    
    if(PH_val == "PH"){
      TsigmaHC <- TsigmaHR <-TsigmaFC <- TsigmaFR <- Tsigma
      TnuHC <- TnuHR <-TnuFC <- TnuFR <- Tnu
      TthetaHC <- TthetaHR <-TthetaFC <- TthetaFR <- Ttheta
      betaZHC <- betaZHR <-betaZFC <- betaZFR <- betaZ
      
    }
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
      
      covaZ <- Zdata_train
      if(PH_val == "NPH"){
        covaZ <- covaZ[,-c(4,5)]
      }
      
      Sn(time = newtimes, sigma = Tsigma, nu = Tnu, theta = Ttheta, beta = betaZ, covariates = covaZ[i,])
      
    }
    
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
        
        covaZ <- Zdata_train
        if(PH_val == "NPH"){
          covaZ <- covaZ[,-c(4,5)]
        }
        
        Sn(time = newtimes, sigma = Tsigma, nu = Tnu, theta = Ttheta, beta = betaZ, covariates = covaZ[i,]) }
      
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
    colnames(Zdata_val) <- c("stage2", "stage3", "agey10", "sex", "colon")
    
    
    funtheo_val <- function(i){
      if(Zdata_val[i,"sex"] == 1 & Zdata_val[i,"colon"] == 1){
        Tsigma = TsigmaHC
        Tnu = TnuHC
        Ttheta = TthetaHC
        betaZ = betaZHC
      }
      # Hommes rectum
      if(Zdata_val[i,"sex"] == 1 & Zdata_val[i,"colon"] == 0){
        Tsigma = TsigmaHR
        Tnu = TnuHR
        Ttheta = TthetaHR
        betaZ = betaZHR
      }
      # Femmes colon
      if(Zdata_val[i,"sex"] == 2 & Zdata_val[i,"colon"] == 1){
        Tsigma = TsigmaFC
        Tnu = TnuFC
        Ttheta = TthetaFC
        betaZ = betaZFC
      }
      # Femmes rectum
      if(Zdata_val[i,"sex"] == 2 & Zdata_val[i,"colon"] == 0){
        Tsigma = TsigmaFR
        Tnu = TnuFR
        Ttheta = TthetaFR
        betaZ = betaZFR
      }
      
      covaZ <- Zdata_val
      if(PH_val == "NPH"){
        covaZ <- covaZ[,-c(4,5)]
      }
      
      Sn(time = newtimes, sigma = Tsigma, nu = Tnu, theta = Ttheta, beta = betaZ, covariates = covaZ[i,]) }
    
    theopredval <- t(sapply(1:dim(Zdata_val)[1], FUN = funtheo_val ))
    
    colnames(theopredval) <- newtimes
    
    mean.theoval <- apply(theopredval, FUN="mean", MARGIN=2) 
    
    ###strates
    #####
    
    
    for(j in strata_names){
      
      data_val_var <- get(paste0("data_valid_", j))
      
      Zdata_val <- cbind(data_val_var$stage2, data_val_var$stage3, data_val_var$agey10, data_val_var$sex, data_val_var$colon)
      colnames(Zdata_val) <- c("stage2", "stage3", "agey10", "sex", "colon")
      
      funtheovalstrat <- function(i){
        if(Zdata_val[i,"sex"] == 1 & Zdata_val[i,"colon"] == 1){
          Tsigma = TsigmaHC
          Tnu = TnuHC
          Ttheta = TthetaHC
          betaZ = betaZHC
        }
        # Hommes rectum
        if(Zdata_val[i,"sex"] == 1 & Zdata_val[i,"colon"] == 0){
          Tsigma = TsigmaHR
          Tnu = TnuHR
          Ttheta = TthetaHR
          betaZ = betaZHR
        }
        # Femmes colon
        if(Zdata_val[i,"sex"] == 2 & Zdata_val[i,"colon"] == 1){
          Tsigma = TsigmaFC
          Tnu = TnuFC
          Ttheta = TthetaFC
          betaZ = betaZFC
        }
        # Femmes rectum
        if(Zdata_val[i,"sex"] == 2 & Zdata_val[i,"colon"] == 0){
          Tsigma = TsigmaFR
          Tnu = TnuFR
          Ttheta = TthetaFR
          betaZ = betaZFR
        }
        
        covaZ <- Zdata_val
        if(PH_val == "NPH"){
          covaZ <- covaZ[,-c(4,5)]
        }
        
        Sn(time = newtimes, sigma = Tsigma, nu = Tnu, theta = Ttheta, beta = betaZ, covariates = covaZ[i,]) }
      
      theopredvalstrat <- t(sapply(1:dim(Zdata_val)[1], FUN = funtheovalstrat ))
      
      colnames(theopredvalstrat) <- newtimes
      
      mean.theoval_strat <- apply(theopredvalstrat, FUN = "mean", MARGIN = 2)
      
      assign(paste0("Zdata_val_", j), Zdata_val)
      assign(paste0("theopredval_", j), theopredvalstrat)
      assign(paste0("mean.theoval_", j), mean.theoval_strat)  
      
    }
    
    ##################DATA TRAIN
    #################################### MOYENNES
    #######THEO
    
    mean_theo <- t(data.frame(THEO = mean.theo))
    write.table(mean_theo, paste0(path0, "TRAIN/WHOLE/THEO/mean_survival/",i,"_mean_theo.csv"), sep = ";", row.names = F, col.names = T)
    
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
    
    #### validatio,
    
    ##################DATA VALID
    #################################### MOYENNES
    
    #######THEO
    
    mean_theoval <- t(data.frame(THEO = mean.theoval))
    write.table(mean_theoval, paste0(path0, "VALID/WHOLE/THEO/mean_survival/",
                                     i,"_mean_theoval.csv"), sep = ";", row.names = F, col.names = T)
    
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

    #################################### SURVIES INDIVIDUELLES
    
    #######THEO
    
    ind_theo <- data.frame(theopred)
    write.table(ind_theo, paste0(path0, "TRAIN/WHOLE/THEO/individual_survival/",
                                 i,"_ind_theo.csv"), sep = ";", row.names = F, col.names = T)
    
    
    ###strates 
    ######
    
    #theo
    
    for (j in strata_names) {
      
      ind_theo <- data.frame(get(paste0("theopred_", j)))
      
      write.table(ind_theo, paste0(path0, "TRAIN/STRATA/THEO/individual_survival/", i, "_ind_theo_", j, ".csv"),
                  sep = ";",row.names = F, col.names = T)
    }
    
    #################################### SURVIES INDIVIDUELLES
    
    #######THEO
    ind_theoval <- data.frame(theopredval)
    write.table(ind_theoval, 
                paste0(path0, "VALID/WHOLE/THEO/individual_survival/",i
                       ,"_ind_theoval.csv"), sep = ";", row.names = F, col.names = T)
    
    ###strates ####
    ######
    
    #theo
    
    for (j in strata_names) {
      
      ind_theoval <- data.frame(get(paste0("theopredval_", j)))
      
      write.table(ind_theoval, paste0(path0, "VALID/STRATA/THEO/individual_survival/", i, "_ind_theoval_", j, ".csv"),
                  sep = ";",row.names = F, col.names = T)
    }
    
    
    }

mclapply(iterations, correc_theo, N = 1000, mc.cores = detectCores() - 4)

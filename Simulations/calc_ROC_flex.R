cens_val <- "HC"
PH_val <- "PH"

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


calc_indic <- function(N, iterations){
  
  strata_names <- c("HC", "HR", "FC", "FR")
  newtimes <- c(365.241*3, 365.241*5, 365.241*10) 
  
  start <- Sys.time()
  
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

    #flex2
    assign(paste0("indTrainFLEX2.2_", k) , read.csv(paste0(path0, "TRAIN/WHOLE/FLEX2.2/individual_survival/",k,"_ind_flex22.csv"),
                                                    sep = " ") )
    assign(paste0("indTrainFLEX2.4_", k) , read.csv(paste0(path0, "TRAIN/WHOLE/FLEX2.4/individual_survival/",k,"_ind_flex24.csv"),
                                                    sep = " ") )
  }
  
  ##### calcul de ROC net
  
  AUC_calc_W <- function(k) {
    datatrain <- get(paste0("DATATrain_", k))
    ind_estimF2.2 <- get(paste0("indTrainFLEX2.2_", k))
    ind_estimF2.4 <- get(paste0("indTrainFLEX2.4_", k))

    rm(list = c(paste0("DATATrain_",k),
                paste0("indTrainFLEX2.2_",k), paste0("indTrainFLEX2.4_",k)))
    
    # Initialiser les vecteurs de résultats
    hold_F2.2 <- hold_F2.4 <- c()
    # sort(unique(1-ind_estimP[, l + 1]))
    # unique( quantile( 1-ind_estimP[, l + 1],  probs= seq(0,1, by = 0.01)) ) #pas en dessous de 20 points 0.025
    
    for (l in seq_along(newtimes)) {
    
      hold_F2.2 <- c(hold_F2.2, roc.net(datatrain$times, datatrain$status, 1-ind_estimF2.2[, l], datatrain$age, datatrain$sexchara,
                                        datatrain$year, slopop, pro.time = newtimes[l], 
                                        cut.off = unique( quantile( 1-ind_estimF2.2[, l],  probs= seq(0,1, by = 0.01)) ))$auc)
      hold_F2.4 <- c(hold_F2.4, roc.net(datatrain$times, datatrain$status, 1-ind_estimF2.4[, l], datatrain$age, datatrain$sexchara,
                                        datatrain$year, slopop, pro.time = newtimes[l],
                                        cut.off = unique( quantile( 1-ind_estimF2.4[, l],  probs= seq(0,1, by = 0.01)) ))$auc)
    }
    
    return(list(
      F2.2  = hold_F2.2,
      F2.4  = hold_F2.4
    ))
  }
  
  results_list <- mclapply(iterations, AUC_calc_W, mc.cores = detectCores() - 4)
  
  #supression des bases inutiles 
  for(k in iterations){
    rm(list = c(paste0("DATATrain_",k),
                paste0("indTrainFLEX2.2_",k), paste0("indTrainFLEX2.4_",k)))
  }
  # Transformer en matrices pour chaque méthode :
  AUC_WHOLE_F2.2 <- do.call(rbind, lapply(results_list, function(x) x$F2.2))
  AUC_WHOLE_F2.4 <- do.call(rbind, lapply(results_list, function(x) x$F2.4))

  
  ### pour renommer les colonnes avec les temps de prognostic
  for(l in c("F2.2","F2.4")){
    dathold <- get(paste0("AUC_WHOLE_", l))
    colnames(dathold) <- c("3 years", "5 years", "10 years")
    rownames(dathold) <- c(iterations)
    assign(paste0("AUC_WHOLE_", l),dathold)
  }
  ## moyennes par temps 
  
  for(l in c("F2.2","F2.4")){
    
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
    
      #flex2
      assign(paste0("indTrainFLEX2.2_",j,"_", k) , read.csv(paste0(path0, "TRAIN/STRATA/FLEX2.2/individual_survival/",k,"_ind_flex22_",j,".csv"),
                                                            sep = " ") )
      assign(paste0("indTrainFLEX2.4_",j,"_", k) , read.csv(paste0(path0, "TRAIN/STRATA/FLEX2.4/individual_survival/",k,"_ind_flex24_",j,".csv"),
                                                            sep = " ") )
  
    }
  }
  
  ##### calcul de ROC net
  
  AUC_calc_S <- function(k) {
    auc_results <- list()
    
    for (j in strata_names) {
      datatrain <- get(paste0("DATATrain_", j, "_", k))
      ind_estimF2.2 <- get(paste0("indTrainFLEX2.2_", j, "_", k))
      ind_estimF2.4 <- get(paste0("indTrainFLEX2.4_", j, "_", k))
      
      
      hold_F2.2 <- hold_F2.4 <- c()
      
      for (l in seq_along(newtimes)) {
        

        hold_F2.2 <- c(hold_F2.2, roc.net(datatrain$times, datatrain$status, 1-ind_estimF2.2[, l], datatrain$age, datatrain$sexchara,
                                          datatrain$year, slopop, pro.time = newtimes[l], 
                                          cut.off = unique( quantile( 1-ind_estimF2.2[, l],  probs= seq(0,1, by = 0.01)) ))$auc)
        hold_F2.4 <- c(hold_F2.4, roc.net(datatrain$times, datatrain$status, 1-ind_estimF2.4[, l], datatrain$age, datatrain$sexchara,
                                          datatrain$year, slopop, pro.time = newtimes[l],
                                          cut.off = unique( quantile( 1-ind_estimF2.4[, l],  probs= seq(0,1, by = 0.01)) ))$auc)
       
      }
      
      auc_results[[j]] <- list(
        k = k,
        F2.2 = hold_F2.2,
        F2.4 = hold_F2.4
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
        paste0("indTrainFLEX2.2_", j, "_", k),
        paste0("indTrainFLEX2.4_", j, "_", k)
      ))
    }
  }
  
  # Transformer en matrices pour chaque méthode :
  
  for(j in strata_names){
    assign(paste0("AUC_",j,"_F2.2"), do.call(rbind, lapply(res_all, function(x) x[[j]]$F2.2)) )
    assign(paste0("AUC_",j,"_F2.4"), do.call(rbind, lapply(res_all, function(x) x[[j]]$F2.4)) )
  }
  ### pour renommer les colonnes avec les temps de prognostic
  for(l in c("F2.2","F2.4")){
    for(j in strata_names){
      dathold <- get(paste0("AUC_", j,"_",l))
      colnames(dathold) <- c( "3 years", "5 years", "10 years")
      rownames(dathold) <- c(iterations)
      assign(paste0("AUC_",j,"_", l),dathold)
    }
  }
  ## moyennes par temps 
  
  for(l in c("F2.2","F2.4")){
    for(j in strata_names){
      assign(paste0("AUC_means_",j,"_",l) ,  colMeans( get(paste0("AUC_",j,"_",l)) ) )
    }
  }
  
  
  ### sous forme de liste pour libérer de l'espace dans l'evt (pour les moyennes)
  
  ROCmean_results <- list()
  
  for(l in c("F2.2","F2.4")){
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
  
  for(l in c("F2.2","F2.4")){
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

    #flex2
    assign(paste0("indValidFLEX2.2_",k) , read.csv(paste0(path0, "VALID/WHOLE/FLEX2.2/individual_survival/",k,"_ind_flexval22.csv"),
                                                   sep = " ") )
    assign(paste0("indValidFLEX2.4_",k) , read.csv(paste0(path0, "VALID/WHOLE/FLEX2.4/individual_survival/",k,"_ind_flexval24.csv"),
                                                   sep = " ") )
   
  }
  
  ##### calcul de ROC net
  
  ##creation d'objets pour enregistrer les AUC de chaque itération aux différents temps et pour les différents modèles 
  
  AUC_calc_val_W <- function(k){ 
    
    dataval <- get(paste0("DATAValid_",k))
    ind_estimF2.2 <- get(paste0("indValidFLEX2.2_", k))
    ind_estimF2.4 <- get(paste0("indValidFLEX2.4_", k))

    rm(list = c(paste0("DATAValid_",k),
                paste0("indValidFLEX2.2_", k),paste0("indValidFLEX2.4_", k)))
    
    # Initialiser les vecteurs de résultats
    hold_F2.2 <- hold_F2.4 <- c()
    
    ###PLANN
    for (l in seq_along(newtimes)) {
     
      hold_F2.2 <- c(hold_F2.2, roc.net(dataval$times, dataval$status, 1-ind_estimF2.2[, l], dataval$age, dataval$sexchara,
                                        dataval$year, slopop, pro.time = newtimes[l], 
                                        cut.off = unique( quantile( 1-ind_estimF2.2[, l],  probs= seq(0,1, by = 0.01)) ))$auc)
      hold_F2.4 <- c(hold_F2.4, roc.net(dataval$times, dataval$status, 1-ind_estimF2.4[, l], dataval$age, dataval$sexchara,
                                        dataval$year, slopop, pro.time = newtimes[l],
                                        cut.off = unique( quantile( 1-ind_estimF2.4[, l],  probs= seq(0,1, by = 0.01)) ))$auc)
      
    }
    
    return(list(
      F2.2  = hold_F2.2,
      F2.4  = hold_F2.4
    ))
  }
  
  resultval_list <- mclapply(iterations, AUC_calc_val_W, mc.cores = detectCores() - 4)
  
  #supression des bases inutiles 
  for(k in iterations){
    rm(list = c(paste0("DATAValid_",k),
                paste0("indValidFLEX2.2_",k), paste0("indValidFLEX2.4_",k)))
  }
  # Transformer en matrices pour chaque méthode :
  AUCval_WHOLE_F2.2 <- do.call(rbind, lapply(resultval_list, function(x) x$F2.2))
  AUCval_WHOLE_F2.4 <- do.call(rbind, lapply(resultval_list, function(x) x$F2.4))

  ### pour renommer les colonnes avec les temps de prognostic
  for(l in c("F2.2","F2.4")){
    dathold <- get(paste0("AUCval_WHOLE_", l))
    colnames(dathold) <- c("3 years", "5 years", "10 years")
    rownames(dathold) <- c(iterations)
    assign(paste0("AUCval_WHOLE_", l),dathold)
  }
  ## moyennes par temps 
  
  for(l in c("F2.2","F2.4")){
    
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

      #flex2
      assign(paste0("indValidFLEX2.2_",j,"_", k) , read.csv(paste0(path0, "VALID/STRATA/FLEX2.2/individual_survival/",k,"_ind_flexval22_",j,".csv"),
                                                            sep = " ") )
      assign(paste0("indValidFLEX2.4_",j,"_", k) , read.csv(paste0(path0, "VALID/STRATA/FLEX2.4/individual_survival/",k,"_ind_flexval24_",j,".csv"),
                                                            sep = " ") )
      
    }
  }
  
  ##### calcul de ROC net
  
  AUC_calc_val_S <- function(k){ 
    auc_results <- list()
    for(j in strata_names){
      dataval <- get(paste0("DATAValid_",j,"_", k))
      ind_estimF2.2 <- get(paste0("indValidFLEX2.2_",j,"_", k))
      ind_estimF2.4 <- get(paste0("indValidFLEX2.4_",j,"_", k))

      rm(list = c(paste0("DATAValid_",j,"_", k),
                  paste0("indValidFLEX2.2_",j,"_", k), paste0("indValidFLEX2.4_",j,"_", k)))
      
      hold_F2.2 <- hold_F2.4 <- c()
      
      for (l in seq_along(newtimes)) {
        
        hold_F2.2 <- c(hold_F2.2, roc.net(dataval$times, dataval$status, 1-ind_estimF2.2[, l], dataval$age, dataval$sexchara,
                                          dataval$year, slopop, pro.time = newtimes[l], 
                                          cut.off = unique( quantile( 1-ind_estimF2.2[, l],  probs= seq(0,1, by = 0.01)) ))$auc)
        hold_F2.4 <- c(hold_F2.4, roc.net(dataval$times, dataval$status, 1-ind_estimF2.4[, l], dataval$age, dataval$sexchara,
                                          dataval$year, slopop, pro.time = newtimes[l],
                                          cut.off = unique( quantile( 1-ind_estimF2.4[, l],  probs= seq(0,1, by = 0.01)) ))$auc)
      }
      
      auc_results[[j]] <- list(
        k = k,
        F2.2 = hold_F2.2,
        F2.4 = hold_F2.4
      )
    }
    
    return(auc_results)
  }
  
  resval_all <- mclapply(iterations, AUC_calc_val_S, mc.cores = detectCores() - 4)
  
  for(k in iterations){
    for(j in strata_names){
      rm(list = c(
        paste0("DATAValid_", j, "_", k),
      
        paste0("indValidFLEX2.2_", j, "_", k),
        paste0("indValidFLEX2.4_", j, "_", k)
       
      ))
    }
  }
  
  # Transformer en matrices pour chaque méthode :
  
  for(j in strata_names){
    assign(paste0("AUCval_",j,"_F2.2"), do.call(rbind, lapply(resval_all, function(x) x[[j]]$F2.2)) )
    assign(paste0("AUCval_",j,"_F2.4"), do.call(rbind, lapply(resval_all, function(x) x[[j]]$F2.4)) )
  }
  ### pour renommer les colonnes avec les temps de prognostic
  for(l in c("F2.2","F2.4")){
    for(j in strata_names){
      dathold <- get(paste0("AUCval_", j,"_",l))
      colnames(dathold) <- c("3 years", "5 years", "10 years")
      rownames(dathold) <- c(iterations)
      assign(paste0("AUCval_",j,"_", l),dathold)
    }
  }
  ## moyennes par temps 
  
  for(l in c("F2.2","F2.4")){
    for(j in strata_names){
      assign(paste0("AUCval_means_",j,"_",l) ,  colMeans( get(paste0("AUCval_",j,"_",l)) ) )
    }
  }
  
  
  ### sous forme de liste pour libérer de l'espace dans l'evt (pour les moyennes)
  
  ROCmeanval_results <- list()
  
  for(l in c("F2.2","F2.4")){
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
  
  for(l in c("F2.2","F2.4")){
    ROCval_results[['WHOLE']][[l]] <- list()
    ROCval_results[['WHOLE']][[l]] <- get( paste0("AUCval_WHOLE_",l) )
    rm(list = paste0("AUCval_WHOLE_",l) )
    for(j in strata_names){
      ROCval_results[[j]][[l]] <- get(paste0("AUCval_", j, "_", l)) 
      rm(list = paste0("AUCval_", j, "_", l))
    }
    
  }
  
  print(Sys.time()-start)
  
  
  # save.image(paste0("~/Documents/Rstudio/Simulations/Simulations mai 2025/Résultats/",length(iterations),"ite_",N,"ind_",date_launch,".Rdata"))
  save(list = ls(), file = paste0("~/Documents/Rstudio/Simulations/Simulations mai 2025/Résultats/N",N,"_results/",PH_val,"/",length(iterations),"ite_",N,"ind_",PH_val,"_",cens_val,"_",date_launch,".Rdata"))
  
}


iterations <- 1:1000

N = 1000
path0 <- paste0("~/Documents/Rstudio/Simulations/Simulations mai 2025/Résultats/N",N,"_results/",PH_val,"/Output simulations_N",N,"_",PH_val,"_",cens_val,"/")
path1 <- "~/Documents/Rstudio/Simulations/BASES/"

rm(list = setdiff(ls(envir = .GlobalEnv), 
                  c("path0", "date_launch", "fr.ratetable", "slopop", 
                    "colrec", "Sn", "simulate_iteration", "calc_indic", "calc_indic_PP", "cens_val", "PH_val")), envir = .GlobalEnv) ## cleaning de l'environnement excpeté ce qu iva être réutilisé
## 1000SiEND

############################

## 1000IndSTA
# path0 <- "~/Documents/Rstudio/Simulations/Simulations mai 2025/Résultats/N1000_results/NPH/Output Simulations_N1000_NPH_HC/"

indic = c()
for(i in 1:1000){
  if(!file.exists(paste0(path0,"DATAFRAMES/TRAIN/WHOLE/",i,"_dftrain.csv"))){ indic <<- c(indic, i)} 
}

if(!is.null(indic)){
  iterations <- c(1:1000)[-indic]
}else{
  iterations <- c(1:1000)
}


calc_indic(N = 1000, iterations = iterations)

rm(list = setdiff(ls(envir = .GlobalEnv), 
                  c("path0", "date_launch", "fr.ratetable", "slopop", 
                    "colrec", "Sn", "simulate_iteration", "calc_indic", "calc_indic_PP", "iterations", "cens_val", "PH_val")), envir = .GlobalEnv)

##########
  ##3000
##########
iterations <- 1:1000

N = 3000
path0 <- paste0("~/Documents/Rstudio/Simulations/Simulations mai 2025/Résultats/N",N,"_results/",PH_val,"/Output simulations_N",N,"_",PH_val,"_",cens_val,"/")
path1 <- "~/Documents/Rstudio/Simulations/BASES/"


rm(list = setdiff(ls(envir = .GlobalEnv), 
                  c("path0", "date_launch", "fr.ratetable", "slopop", 
                    "colrec", "Sn", "simulate_iteration", "calc_indic", "calc_indic_PP", "cens_val", "PH_val")), envir = .GlobalEnv) ## cleaning de l'environnement excpeté ce qu iva être réutilisé
## 3000SiEND

############################

## 3000IndSTA
# path0 <- "~/Documents/Rstudio/Simulations/Simulations mai 2025/Résultats/N3000_results/NPH/Output Simulations_N3000_NPH_HC/"

indic = c()
for(i in 1:1000){
  if(!file.exists(paste0(path0,"DATAFRAMES/TRAIN/WHOLE/",i,"_dftrain.csv"))){ indic <<- c(indic, i)} 
}

if(!is.null(indic)){
  iterations <- c(1:1000)[-indic]
}else{
  iterations <- c(1:1000)
}


calc_indic(N = 3000, iterations = iterations)

rm(list = setdiff(ls(envir = .GlobalEnv), 
                  c("path0", "date_launch", "fr.ratetable", "slopop", 
                    "colrec", "Sn", "simulate_iteration", "calc_indic", "calc_indic_PP", "iterations", "cens_val", "PH_val")), envir = .GlobalEnv)

### 3000IndEND



####################################################################################################
#                                             5000
####################################################################################################

iterations <- 1:1000

N = 5000
path0 <- paste0("~/Documents/Rstudio/Simulations/Simulations mai 2025/Résultats/N",N,"_results/",PH_val,"/Output simulations_N",N,"_",PH_val,"_",cens_val,"/")
path1 <- "~/Documents/Rstudio/Simulations/BASES/"

## 5000SiEND

############################

## 5000IndSTA
# path0 <- "~/Documents/Rstudio/Simulations/Simulations mai 2025/Résultats/N5000_results/NPH/Output Simulations_N5000_NPH_HC/"

indic = c()
for(i in 1:1000){
  if(!file.exists(paste0(path0,"DATAFRAMES/TRAIN/WHOLE/",i,"_dftrain.csv"))){ indic <<- c(indic, i)} 
}

if(!is.null(indic)){
  iterations <- c(1:1000)[-indic]
}else{
  iterations <- c(1:1000)
}


calc_indic(N = 5000, iterations = iterations)

rm(list = setdiff(ls(envir = .GlobalEnv), 
                  c("path0", "date_launch", "fr.ratetable", "slopop", 
                    "colrec", "Sn", "simulate_iteration", "calc_indic", "calc_indic_PP", "iterations", "cens_val", "PH_val")), envir = .GlobalEnv)

### 5000IndEND


##############################################################################################################################


rm(list = ls()) 


##############################################################################################################################


cens_val <- "LC"
PH_val <- "PH"

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


calc_indic <- function(N, iterations){
  
  strata_names <- c("HC", "HR", "FC", "FR")
  newtimes <- c(365.241*3, 365.241*5, 365.241*10) 
  
  start <- Sys.time()
  
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
    
    #flex2
    assign(paste0("indTrainFLEX2.2_", k) , read.csv(paste0(path0, "TRAIN/WHOLE/FLEX2.2/individual_survival/",k,"_ind_flex22.csv"),
                                                    sep = " ") )
    assign(paste0("indTrainFLEX2.4_", k) , read.csv(paste0(path0, "TRAIN/WHOLE/FLEX2.4/individual_survival/",k,"_ind_flex24.csv"),
                                                    sep = " ") )
  }
  
  ##### calcul de ROC net
  
  AUC_calc_W <- function(k) {
    datatrain <- get(paste0("DATATrain_", k))
    ind_estimF2.2 <- get(paste0("indTrainFLEX2.2_", k))
    ind_estimF2.4 <- get(paste0("indTrainFLEX2.4_", k))
    
    rm(list = c(paste0("DATATrain_",k),
                paste0("indTrainFLEX2.2_",k), paste0("indTrainFLEX2.4_",k)))
    
    # Initialiser les vecteurs de résultats
    hold_F2.2 <- hold_F2.4 <- c()
    # sort(unique(1-ind_estimP[, l + 1]))
    # unique( quantile( 1-ind_estimP[, l + 1],  probs= seq(0,1, by = 0.01)) ) #pas en dessous de 20 points 0.025
    
    for (l in seq_along(newtimes)) {
      
      hold_F2.2 <- c(hold_F2.2, roc.net(datatrain$times, datatrain$status, 1-ind_estimF2.2[, l], datatrain$age, datatrain$sexchara,
                                        datatrain$year, slopop, pro.time = newtimes[l], 
                                        cut.off = unique( quantile( 1-ind_estimF2.2[, l],  probs= seq(0,1, by = 0.01)) ))$auc)
      hold_F2.4 <- c(hold_F2.4, roc.net(datatrain$times, datatrain$status, 1-ind_estimF2.4[, l], datatrain$age, datatrain$sexchara,
                                        datatrain$year, slopop, pro.time = newtimes[l],
                                        cut.off = unique( quantile( 1-ind_estimF2.4[, l],  probs= seq(0,1, by = 0.01)) ))$auc)
    }
    
    return(list(
      F2.2  = hold_F2.2,
      F2.4  = hold_F2.4
    ))
  }
  
  results_list <- mclapply(iterations, AUC_calc_W, mc.cores = detectCores() - 4)
  
  #supression des bases inutiles 
  for(k in iterations){
    rm(list = c(paste0("DATATrain_",k),
                paste0("indTrainFLEX2.2_",k), paste0("indTrainFLEX2.4_",k)))
  }
  # Transformer en matrices pour chaque méthode :
  AUC_WHOLE_F2.2 <- do.call(rbind, lapply(results_list, function(x) x$F2.2))
  AUC_WHOLE_F2.4 <- do.call(rbind, lapply(results_list, function(x) x$F2.4))
  
  
  ### pour renommer les colonnes avec les temps de prognostic
  for(l in c("F2.2","F2.4")){
    dathold <- get(paste0("AUC_WHOLE_", l))
    colnames(dathold) <- c("3 years", "5 years", "10 years")
    rownames(dathold) <- c(iterations)
    assign(paste0("AUC_WHOLE_", l),dathold)
  }
  ## moyennes par temps 
  
  for(l in c("F2.2","F2.4")){
    
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
      
      #flex2
      assign(paste0("indTrainFLEX2.2_",j,"_", k) , read.csv(paste0(path0, "TRAIN/STRATA/FLEX2.2/individual_survival/",k,"_ind_flex22_",j,".csv"),
                                                            sep = " ") )
      assign(paste0("indTrainFLEX2.4_",j,"_", k) , read.csv(paste0(path0, "TRAIN/STRATA/FLEX2.4/individual_survival/",k,"_ind_flex24_",j,".csv"),
                                                            sep = " ") )
      
    }
  }
  
  ##### calcul de ROC net
  
  AUC_calc_S <- function(k) {
    auc_results <- list()
    
    for (j in strata_names) {
      datatrain <- get(paste0("DATATrain_", j, "_", k))
      ind_estimF2.2 <- get(paste0("indTrainFLEX2.2_", j, "_", k))
      ind_estimF2.4 <- get(paste0("indTrainFLEX2.4_", j, "_", k))
      
      
      hold_F2.2 <- hold_F2.4 <- c()
      
      for (l in seq_along(newtimes)) {
        
        
        hold_F2.2 <- c(hold_F2.2, roc.net(datatrain$times, datatrain$status, 1-ind_estimF2.2[, l], datatrain$age, datatrain$sexchara,
                                          datatrain$year, slopop, pro.time = newtimes[l], 
                                          cut.off = unique( quantile( 1-ind_estimF2.2[, l],  probs= seq(0,1, by = 0.01)) ))$auc)
        hold_F2.4 <- c(hold_F2.4, roc.net(datatrain$times, datatrain$status, 1-ind_estimF2.4[, l], datatrain$age, datatrain$sexchara,
                                          datatrain$year, slopop, pro.time = newtimes[l],
                                          cut.off = unique( quantile( 1-ind_estimF2.4[, l],  probs= seq(0,1, by = 0.01)) ))$auc)
        
      }
      
      auc_results[[j]] <- list(
        k = k,
        F2.2 = hold_F2.2,
        F2.4 = hold_F2.4
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
        paste0("indTrainFLEX2.2_", j, "_", k),
        paste0("indTrainFLEX2.4_", j, "_", k)
      ))
    }
  }
  
  # Transformer en matrices pour chaque méthode :
  
  for(j in strata_names){
    assign(paste0("AUC_",j,"_F2.2"), do.call(rbind, lapply(res_all, function(x) x[[j]]$F2.2)) )
    assign(paste0("AUC_",j,"_F2.4"), do.call(rbind, lapply(res_all, function(x) x[[j]]$F2.4)) )
  }
  ### pour renommer les colonnes avec les temps de prognostic
  for(l in c("F2.2","F2.4")){
    for(j in strata_names){
      dathold <- get(paste0("AUC_", j,"_",l))
      colnames(dathold) <- c( "3 years", "5 years", "10 years")
      rownames(dathold) <- c(iterations)
      assign(paste0("AUC_",j,"_", l),dathold)
    }
  }
  ## moyennes par temps 
  
  for(l in c("F2.2","F2.4")){
    for(j in strata_names){
      assign(paste0("AUC_means_",j,"_",l) ,  colMeans( get(paste0("AUC_",j,"_",l)) ) )
    }
  }
  
  
  ### sous forme de liste pour libérer de l'espace dans l'evt (pour les moyennes)
  
  ROCmean_results <- list()
  
  for(l in c("F2.2","F2.4")){
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
  
  for(l in c("F2.2","F2.4")){
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
    
    #flex2
    assign(paste0("indValidFLEX2.2_",k) , read.csv(paste0(path0, "VALID/WHOLE/FLEX2.2/individual_survival/",k,"_ind_flexval22.csv"),
                                                   sep = " ") )
    assign(paste0("indValidFLEX2.4_",k) , read.csv(paste0(path0, "VALID/WHOLE/FLEX2.4/individual_survival/",k,"_ind_flexval24.csv"),
                                                   sep = " ") )
    
  }
  
  ##### calcul de ROC net
  
  ##creation d'objets pour enregistrer les AUC de chaque itération aux différents temps et pour les différents modèles 
  
  AUC_calc_val_W <- function(k){ 
    
    dataval <- get(paste0("DATAValid_",k))
    ind_estimF2.2 <- get(paste0("indValidFLEX2.2_", k))
    ind_estimF2.4 <- get(paste0("indValidFLEX2.4_", k))
    
    rm(list = c(paste0("DATAValid_",k),
                paste0("indValidFLEX2.2_", k),paste0("indValidFLEX2.4_", k)))
    
    # Initialiser les vecteurs de résultats
    hold_F2.2 <- hold_F2.4 <- c()
    
    ###PLANN
    for (l in seq_along(newtimes)) {
      
      hold_F2.2 <- c(hold_F2.2, roc.net(dataval$times, dataval$status, 1-ind_estimF2.2[, l], dataval$age, dataval$sexchara,
                                        dataval$year, slopop, pro.time = newtimes[l], 
                                        cut.off = unique( quantile( 1-ind_estimF2.2[, l],  probs= seq(0,1, by = 0.01)) ))$auc)
      hold_F2.4 <- c(hold_F2.4, roc.net(dataval$times, dataval$status, 1-ind_estimF2.4[, l], dataval$age, dataval$sexchara,
                                        dataval$year, slopop, pro.time = newtimes[l],
                                        cut.off = unique( quantile( 1-ind_estimF2.4[, l],  probs= seq(0,1, by = 0.01)) ))$auc)
      
    }
    
    return(list(
      F2.2  = hold_F2.2,
      F2.4  = hold_F2.4
    ))
  }
  
  resultval_list <- mclapply(iterations, AUC_calc_val_W, mc.cores = detectCores() - 4)
  
  #supression des bases inutiles 
  for(k in iterations){
    rm(list = c(paste0("DATAValid_",k),
                paste0("indValidFLEX2.2_",k), paste0("indValidFLEX2.4_",k)))
  }
  # Transformer en matrices pour chaque méthode :
  AUCval_WHOLE_F2.2 <- do.call(rbind, lapply(resultval_list, function(x) x$F2.2))
  AUCval_WHOLE_F2.4 <- do.call(rbind, lapply(resultval_list, function(x) x$F2.4))
  
  ### pour renommer les colonnes avec les temps de prognostic
  for(l in c("F2.2","F2.4")){
    dathold <- get(paste0("AUCval_WHOLE_", l))
    colnames(dathold) <- c("3 years", "5 years", "10 years")
    rownames(dathold) <- c(iterations)
    assign(paste0("AUCval_WHOLE_", l),dathold)
  }
  ## moyennes par temps 
  
  for(l in c("F2.2","F2.4")){
    
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
      
      #flex2
      assign(paste0("indValidFLEX2.2_",j,"_", k) , read.csv(paste0(path0, "VALID/STRATA/FLEX2.2/individual_survival/",k,"_ind_flexval22_",j,".csv"),
                                                            sep = " ") )
      assign(paste0("indValidFLEX2.4_",j,"_", k) , read.csv(paste0(path0, "VALID/STRATA/FLEX2.4/individual_survival/",k,"_ind_flexval24_",j,".csv"),
                                                            sep = " ") )
      
    }
  }
  
  ##### calcul de ROC net
  
  AUC_calc_val_S <- function(k){ 
    auc_results <- list()
    for(j in strata_names){
      dataval <- get(paste0("DATAValid_",j,"_", k))
      ind_estimF2.2 <- get(paste0("indValidFLEX2.2_",j,"_", k))
      ind_estimF2.4 <- get(paste0("indValidFLEX2.4_",j,"_", k))
      
      rm(list = c(paste0("DATAValid_",j,"_", k),
                  paste0("indValidFLEX2.2_",j,"_", k), paste0("indValidFLEX2.4_",j,"_", k)))
      
      hold_F2.2 <- hold_F2.4 <- c()
      
      for (l in seq_along(newtimes)) {
        
        hold_F2.2 <- c(hold_F2.2, roc.net(dataval$times, dataval$status, 1-ind_estimF2.2[, l], dataval$age, dataval$sexchara,
                                          dataval$year, slopop, pro.time = newtimes[l], 
                                          cut.off = unique( quantile( 1-ind_estimF2.2[, l],  probs= seq(0,1, by = 0.01)) ))$auc)
        hold_F2.4 <- c(hold_F2.4, roc.net(dataval$times, dataval$status, 1-ind_estimF2.4[, l], dataval$age, dataval$sexchara,
                                          dataval$year, slopop, pro.time = newtimes[l],
                                          cut.off = unique( quantile( 1-ind_estimF2.4[, l],  probs= seq(0,1, by = 0.01)) ))$auc)
      }
      
      auc_results[[j]] <- list(
        k = k,
        F2.2 = hold_F2.2,
        F2.4 = hold_F2.4
      )
    }
    
    return(auc_results)
  }
  
  resval_all <- mclapply(iterations, AUC_calc_val_S, mc.cores = detectCores() - 4)
  
  for(k in iterations){
    for(j in strata_names){
      rm(list = c(
        paste0("DATAValid_", j, "_", k),
        
        paste0("indValidFLEX2.2_", j, "_", k),
        paste0("indValidFLEX2.4_", j, "_", k)
        
      ))
    }
  }
  
  # Transformer en matrices pour chaque méthode :
  
  for(j in strata_names){
    assign(paste0("AUCval_",j,"_F2.2"), do.call(rbind, lapply(resval_all, function(x) x[[j]]$F2.2)) )
    assign(paste0("AUCval_",j,"_F2.4"), do.call(rbind, lapply(resval_all, function(x) x[[j]]$F2.4)) )
  }
  ### pour renommer les colonnes avec les temps de prognostic
  for(l in c("F2.2","F2.4")){
    for(j in strata_names){
      dathold <- get(paste0("AUCval_", j,"_",l))
      colnames(dathold) <- c("3 years", "5 years", "10 years")
      rownames(dathold) <- c(iterations)
      assign(paste0("AUCval_",j,"_", l),dathold)
    }
  }
  ## moyennes par temps 
  
  for(l in c("F2.2","F2.4")){
    for(j in strata_names){
      assign(paste0("AUCval_means_",j,"_",l) ,  colMeans( get(paste0("AUCval_",j,"_",l)) ) )
    }
  }
  
  
  ### sous forme de liste pour libérer de l'espace dans l'evt (pour les moyennes)
  
  ROCmeanval_results <- list()
  
  for(l in c("F2.2","F2.4")){
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
  
  for(l in c("F2.2","F2.4")){
    ROCval_results[['WHOLE']][[l]] <- list()
    ROCval_results[['WHOLE']][[l]] <- get( paste0("AUCval_WHOLE_",l) )
    rm(list = paste0("AUCval_WHOLE_",l) )
    for(j in strata_names){
      ROCval_results[[j]][[l]] <- get(paste0("AUCval_", j, "_", l)) 
      rm(list = paste0("AUCval_", j, "_", l))
    }
    
  }
  
  print(Sys.time()-start)
  
  
  # save.image(paste0("~/Documents/Rstudio/Simulations/Simulations mai 2025/Résultats/",length(iterations),"ite_",N,"ind_",date_launch,".Rdata"))
  save(list = ls(), file = paste0("~/Documents/Rstudio/Simulations/Simulations mai 2025/Résultats/N",N,"_results/",PH_val,"/",length(iterations),"ite_",N,"ind_",PH_val,"_",cens_val,"_",date_launch,".Rdata"))
  
}


iterations <- 1:1000

N = 1000
path0 <- paste0("~/Documents/Rstudio/Simulations/Simulations mai 2025/Résultats/N",N,"_results/",PH_val,"/Output simulations_N",N,"_",PH_val,"_",cens_val,"/")
path1 <- "~/Documents/Rstudio/Simulations/BASES/"

rm(list = setdiff(ls(envir = .GlobalEnv), 
                  c("path0", "date_launch", "fr.ratetable", "slopop", 
                    "colrec", "Sn", "simulate_iteration", "calc_indic", "calc_indic_PP", "cens_val", "PH_val")), envir = .GlobalEnv) ## cleaning de l'environnement excpeté ce qu iva être réutilisé
## 1000SiEND

############################

## 1000IndSTA
# path0 <- "~/Documents/Rstudio/Simulations/Simulations mai 2025/Résultats/N1000_results/NPH/Output Simulations_N1000_NPH_HC/"

indic = c()
for(i in 1:1000){
  if(!file.exists(paste0(path0,"DATAFRAMES/TRAIN/WHOLE/",i,"_dftrain.csv"))){ indic <<- c(indic, i)} 
}

if(!is.null(indic)){
  iterations <- c(1:1000)[-indic]
}else{
  iterations <- c(1:1000)
}


calc_indic(N = 1000, iterations = iterations)

rm(list = setdiff(ls(envir = .GlobalEnv), 
                  c("path0", "date_launch", "fr.ratetable", "slopop", 
                    "colrec", "Sn", "simulate_iteration", "calc_indic", "calc_indic_PP", "iterations", "cens_val", "PH_val")), envir = .GlobalEnv)

##########
##3000
##########
iterations <- 1:1000

N = 3000
path0 <- paste0("~/Documents/Rstudio/Simulations/Simulations mai 2025/Résultats/N",N,"_results/",PH_val,"/Output simulations_N",N,"_",PH_val,"_",cens_val,"/")
path1 <- "~/Documents/Rstudio/Simulations/BASES/"


rm(list = setdiff(ls(envir = .GlobalEnv), 
                  c("path0", "date_launch", "fr.ratetable", "slopop", 
                    "colrec", "Sn", "simulate_iteration", "calc_indic", "calc_indic_PP", "cens_val", "PH_val")), envir = .GlobalEnv) ## cleaning de l'environnement excpeté ce qu iva être réutilisé
## 3000SiEND

############################

## 3000IndSTA
# path0 <- "~/Documents/Rstudio/Simulations/Simulations mai 2025/Résultats/N3000_results/NPH/Output Simulations_N3000_NPH_HC/"

indic = c()
for(i in 1:1000){
  if(!file.exists(paste0(path0,"DATAFRAMES/TRAIN/WHOLE/",i,"_dftrain.csv"))){ indic <<- c(indic, i)} 
}

if(!is.null(indic)){
  iterations <- c(1:1000)[-indic]
}else{
  iterations <- c(1:1000)
}


calc_indic(N = 3000, iterations = iterations)

rm(list = setdiff(ls(envir = .GlobalEnv), 
                  c("path0", "date_launch", "fr.ratetable", "slopop", 
                    "colrec", "Sn", "simulate_iteration", "calc_indic", "calc_indic_PP", "iterations", "cens_val", "PH_val")), envir = .GlobalEnv)

### 3000IndEND



####################################################################################################
#                                             5000
####################################################################################################

iterations <- 1:1000

N = 5000
path0 <- paste0("~/Documents/Rstudio/Simulations/Simulations mai 2025/Résultats/N",N,"_results/",PH_val,"/Output simulations_N",N,"_",PH_val,"_",cens_val,"/")
path1 <- "~/Documents/Rstudio/Simulations/BASES/"

## 5000SiEND

############################

## 5000IndSTA
# path0 <- "~/Documents/Rstudio/Simulations/Simulations mai 2025/Résultats/N5000_results/NPH/Output Simulations_N5000_NPH_HC/"

indic = c()
for(i in 1:1000){
  if(!file.exists(paste0(path0,"DATAFRAMES/TRAIN/WHOLE/",i,"_dftrain.csv"))){ indic <<- c(indic, i)} 
}

if(!is.null(indic)){
  iterations <- c(1:1000)[-indic]
}else{
  iterations <- c(1:1000)
}


calc_indic(N = 5000, iterations = iterations)

rm(list = setdiff(ls(envir = .GlobalEnv), 
                  c("path0", "date_launch", "fr.ratetable", "slopop", 
                    "colrec", "Sn", "simulate_iteration", "calc_indic", "calc_indic_PP", "iterations", "cens_val", "PH_val")), envir = .GlobalEnv)

### 5000IndEND


cens_val <- "HC"
PH_val <- "PH"

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


calc_indic <- function(N, iterations){
  
  strata_names <- c("HC", "HR", "FC", "FR")
  newtimes <- c(365.241*3, 365.241*5, 365.241*10) 
  
  start <- Sys.time()
  
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
    
    #flex2
    assign(paste0("indTrainFLEX2.2_", k) , read.csv(paste0(path0, "TRAIN/WHOLE/FLEX2.2/individual_survival/",k,"_ind_flex22.csv"),
                                                    sep = " ") )
    assign(paste0("indTrainFLEX2.4_", k) , read.csv(paste0(path0, "TRAIN/WHOLE/FLEX2.4/individual_survival/",k,"_ind_flex24.csv"),
                                                    sep = " ") )
  }
  
  ##### calcul de ROC net
  
  AUC_calc_W <- function(k) {
    datatrain <- get(paste0("DATATrain_", k))
    ind_estimF2.2 <- get(paste0("indTrainFLEX2.2_", k))
    ind_estimF2.4 <- get(paste0("indTrainFLEX2.4_", k))
    
    rm(list = c(paste0("DATATrain_",k),
                paste0("indTrainFLEX2.2_",k), paste0("indTrainFLEX2.4_",k)))
    
    # Initialiser les vecteurs de résultats
    hold_F2.2 <- hold_F2.4 <- c()
    # sort(unique(1-ind_estimP[, l + 1]))
    # unique( quantile( 1-ind_estimP[, l + 1],  probs= seq(0,1, by = 0.01)) ) #pas en dessous de 20 points 0.025
    
    for (l in seq_along(newtimes)) {
      
      hold_F2.2 <- c(hold_F2.2, roc.net(datatrain$times, datatrain$status, 1-ind_estimF2.2[, l], datatrain$age, datatrain$sexchara,
                                        datatrain$year, slopop, pro.time = newtimes[l], 
                                        cut.off = unique( quantile( 1-ind_estimF2.2[, l],  probs= seq(0,1, by = 0.01)) ))$auc)
      hold_F2.4 <- c(hold_F2.4, roc.net(datatrain$times, datatrain$status, 1-ind_estimF2.4[, l], datatrain$age, datatrain$sexchara,
                                        datatrain$year, slopop, pro.time = newtimes[l],
                                        cut.off = unique( quantile( 1-ind_estimF2.4[, l],  probs= seq(0,1, by = 0.01)) ))$auc)
    }
    
    return(list(
      F2.2  = hold_F2.2,
      F2.4  = hold_F2.4
    ))
  }
  
  results_list <- mclapply(iterations, AUC_calc_W, mc.cores = detectCores() - 4)
  
  #supression des bases inutiles 
  for(k in iterations){
    rm(list = c(paste0("DATATrain_",k),
                paste0("indTrainFLEX2.2_",k), paste0("indTrainFLEX2.4_",k)))
  }
  # Transformer en matrices pour chaque méthode :
  AUC_WHOLE_F2.2 <- do.call(rbind, lapply(results_list, function(x) x$F2.2))
  AUC_WHOLE_F2.4 <- do.call(rbind, lapply(results_list, function(x) x$F2.4))
  
  
  ### pour renommer les colonnes avec les temps de prognostic
  for(l in c("F2.2","F2.4")){
    dathold <- get(paste0("AUC_WHOLE_", l))
    colnames(dathold) <- c("3 years", "5 years", "10 years")
    rownames(dathold) <- c(iterations)
    assign(paste0("AUC_WHOLE_", l),dathold)
  }
  ## moyennes par temps 
  
  for(l in c("F2.2","F2.4")){
    
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
      
      #flex2
      assign(paste0("indTrainFLEX2.2_",j,"_", k) , read.csv(paste0(path0, "TRAIN/STRATA/FLEX2.2/individual_survival/",k,"_ind_flex22_",j,".csv"),
                                                            sep = " ") )
      assign(paste0("indTrainFLEX2.4_",j,"_", k) , read.csv(paste0(path0, "TRAIN/STRATA/FLEX2.4/individual_survival/",k,"_ind_flex24_",j,".csv"),
                                                            sep = " ") )
      
    }
  }
  
  ##### calcul de ROC net
  
  AUC_calc_S <- function(k) {
    auc_results <- list()
    
    for (j in strata_names) {
      datatrain <- get(paste0("DATATrain_", j, "_", k))
      ind_estimF2.2 <- get(paste0("indTrainFLEX2.2_", j, "_", k))
      ind_estimF2.4 <- get(paste0("indTrainFLEX2.4_", j, "_", k))
      
      
      hold_F2.2 <- hold_F2.4 <- c()
      
      for (l in seq_along(newtimes)) {
        
        
        hold_F2.2 <- c(hold_F2.2, roc.net(datatrain$times, datatrain$status, 1-ind_estimF2.2[, l], datatrain$age, datatrain$sexchara,
                                          datatrain$year, slopop, pro.time = newtimes[l], 
                                          cut.off = unique( quantile( 1-ind_estimF2.2[, l],  probs= seq(0,1, by = 0.01)) ))$auc)
        hold_F2.4 <- c(hold_F2.4, roc.net(datatrain$times, datatrain$status, 1-ind_estimF2.4[, l], datatrain$age, datatrain$sexchara,
                                          datatrain$year, slopop, pro.time = newtimes[l],
                                          cut.off = unique( quantile( 1-ind_estimF2.4[, l],  probs= seq(0,1, by = 0.01)) ))$auc)
        
      }
      
      auc_results[[j]] <- list(
        k = k,
        F2.2 = hold_F2.2,
        F2.4 = hold_F2.4
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
        paste0("indTrainFLEX2.2_", j, "_", k),
        paste0("indTrainFLEX2.4_", j, "_", k)
      ))
    }
  }
  
  # Transformer en matrices pour chaque méthode :
  
  for(j in strata_names){
    assign(paste0("AUC_",j,"_F2.2"), do.call(rbind, lapply(res_all, function(x) x[[j]]$F2.2)) )
    assign(paste0("AUC_",j,"_F2.4"), do.call(rbind, lapply(res_all, function(x) x[[j]]$F2.4)) )
  }
  ### pour renommer les colonnes avec les temps de prognostic
  for(l in c("F2.2","F2.4")){
    for(j in strata_names){
      dathold <- get(paste0("AUC_", j,"_",l))
      colnames(dathold) <- c( "3 years", "5 years", "10 years")
      rownames(dathold) <- c(iterations)
      assign(paste0("AUC_",j,"_", l),dathold)
    }
  }
  ## moyennes par temps 
  
  for(l in c("F2.2","F2.4")){
    for(j in strata_names){
      assign(paste0("AUC_means_",j,"_",l) ,  colMeans( get(paste0("AUC_",j,"_",l)) ) )
    }
  }
  
  
  ### sous forme de liste pour libérer de l'espace dans l'evt (pour les moyennes)
  
  ROCmean_results <- list()
  
  for(l in c("F2.2","F2.4")){
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
  
  for(l in c("F2.2","F2.4")){
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
    
    #flex2
    assign(paste0("indValidFLEX2.2_",k) , read.csv(paste0(path0, "VALID/WHOLE/FLEX2.2/individual_survival/",k,"_ind_flexval22.csv"),
                                                   sep = " ") )
    assign(paste0("indValidFLEX2.4_",k) , read.csv(paste0(path0, "VALID/WHOLE/FLEX2.4/individual_survival/",k,"_ind_flexval24.csv"),
                                                   sep = " ") )
    
  }
  
  ##### calcul de ROC net
  
  ##creation d'objets pour enregistrer les AUC de chaque itération aux différents temps et pour les différents modèles 
  
  AUC_calc_val_W <- function(k){ 
    
    dataval <- get(paste0("DATAValid_",k))
    ind_estimF2.2 <- get(paste0("indValidFLEX2.2_", k))
    ind_estimF2.4 <- get(paste0("indValidFLEX2.4_", k))
    
    rm(list = c(paste0("DATAValid_",k),
                paste0("indValidFLEX2.2_", k),paste0("indValidFLEX2.4_", k)))
    
    # Initialiser les vecteurs de résultats
    hold_F2.2 <- hold_F2.4 <- c()
    
    ###PLANN
    for (l in seq_along(newtimes)) {
      
      hold_F2.2 <- c(hold_F2.2, roc.net(dataval$times, dataval$status, 1-ind_estimF2.2[, l], dataval$age, dataval$sexchara,
                                        dataval$year, slopop, pro.time = newtimes[l], 
                                        cut.off = unique( quantile( 1-ind_estimF2.2[, l],  probs= seq(0,1, by = 0.01)) ))$auc)
      hold_F2.4 <- c(hold_F2.4, roc.net(dataval$times, dataval$status, 1-ind_estimF2.4[, l], dataval$age, dataval$sexchara,
                                        dataval$year, slopop, pro.time = newtimes[l],
                                        cut.off = unique( quantile( 1-ind_estimF2.4[, l],  probs= seq(0,1, by = 0.01)) ))$auc)
      
    }
    
    return(list(
      F2.2  = hold_F2.2,
      F2.4  = hold_F2.4
    ))
  }
  
  resultval_list <- mclapply(iterations, AUC_calc_val_W, mc.cores = detectCores() - 4)
  
  #supression des bases inutiles 
  for(k in iterations){
    rm(list = c(paste0("DATAValid_",k),
                paste0("indValidFLEX2.2_",k), paste0("indValidFLEX2.4_",k)))
  }
  # Transformer en matrices pour chaque méthode :
  AUCval_WHOLE_F2.2 <- do.call(rbind, lapply(resultval_list, function(x) x$F2.2))
  AUCval_WHOLE_F2.4 <- do.call(rbind, lapply(resultval_list, function(x) x$F2.4))
  
  ### pour renommer les colonnes avec les temps de prognostic
  for(l in c("F2.2","F2.4")){
    dathold <- get(paste0("AUCval_WHOLE_", l))
    colnames(dathold) <- c("3 years", "5 years", "10 years")
    rownames(dathold) <- c(iterations)
    assign(paste0("AUCval_WHOLE_", l),dathold)
  }
  ## moyennes par temps 
  
  for(l in c("F2.2","F2.4")){
    
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
      
      #flex2
      assign(paste0("indValidFLEX2.2_",j,"_", k) , read.csv(paste0(path0, "VALID/STRATA/FLEX2.2/individual_survival/",k,"_ind_flexval22_",j,".csv"),
                                                            sep = " ") )
      assign(paste0("indValidFLEX2.4_",j,"_", k) , read.csv(paste0(path0, "VALID/STRATA/FLEX2.4/individual_survival/",k,"_ind_flexval24_",j,".csv"),
                                                            sep = " ") )
      
    }
  }
  
  ##### calcul de ROC net
  
  AUC_calc_val_S <- function(k){ 
    auc_results <- list()
    for(j in strata_names){
      dataval <- get(paste0("DATAValid_",j,"_", k))
      ind_estimF2.2 <- get(paste0("indValidFLEX2.2_",j,"_", k))
      ind_estimF2.4 <- get(paste0("indValidFLEX2.4_",j,"_", k))
      
      rm(list = c(paste0("DATAValid_",j,"_", k),
                  paste0("indValidFLEX2.2_",j,"_", k), paste0("indValidFLEX2.4_",j,"_", k)))
      
      hold_F2.2 <- hold_F2.4 <- c()
      
      for (l in seq_along(newtimes)) {
        
        hold_F2.2 <- c(hold_F2.2, roc.net(dataval$times, dataval$status, 1-ind_estimF2.2[, l], dataval$age, dataval$sexchara,
                                          dataval$year, slopop, pro.time = newtimes[l], 
                                          cut.off = unique( quantile( 1-ind_estimF2.2[, l],  probs= seq(0,1, by = 0.01)) ))$auc)
        hold_F2.4 <- c(hold_F2.4, roc.net(dataval$times, dataval$status, 1-ind_estimF2.4[, l], dataval$age, dataval$sexchara,
                                          dataval$year, slopop, pro.time = newtimes[l],
                                          cut.off = unique( quantile( 1-ind_estimF2.4[, l],  probs= seq(0,1, by = 0.01)) ))$auc)
      }
      
      auc_results[[j]] <- list(
        k = k,
        F2.2 = hold_F2.2,
        F2.4 = hold_F2.4
      )
    }
    
    return(auc_results)
  }
  
  resval_all <- mclapply(iterations, AUC_calc_val_S, mc.cores = detectCores() - 4)
  
  for(k in iterations){
    for(j in strata_names){
      rm(list = c(
        paste0("DATAValid_", j, "_", k),
        
        paste0("indValidFLEX2.2_", j, "_", k),
        paste0("indValidFLEX2.4_", j, "_", k)
        
      ))
    }
  }
  
  # Transformer en matrices pour chaque méthode :
  
  for(j in strata_names){
    assign(paste0("AUCval_",j,"_F2.2"), do.call(rbind, lapply(resval_all, function(x) x[[j]]$F2.2)) )
    assign(paste0("AUCval_",j,"_F2.4"), do.call(rbind, lapply(resval_all, function(x) x[[j]]$F2.4)) )
  }
  ### pour renommer les colonnes avec les temps de prognostic
  for(l in c("F2.2","F2.4")){
    for(j in strata_names){
      dathold <- get(paste0("AUCval_", j,"_",l))
      colnames(dathold) <- c("3 years", "5 years", "10 years")
      rownames(dathold) <- c(iterations)
      assign(paste0("AUCval_",j,"_", l),dathold)
    }
  }
  ## moyennes par temps 
  
  for(l in c("F2.2","F2.4")){
    for(j in strata_names){
      assign(paste0("AUCval_means_",j,"_",l) ,  colMeans( get(paste0("AUCval_",j,"_",l)) ) )
    }
  }
  
  
  ### sous forme de liste pour libérer de l'espace dans l'evt (pour les moyennes)
  
  ROCmeanval_results <- list()
  
  for(l in c("F2.2","F2.4")){
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
  
  for(l in c("F2.2","F2.4")){
    ROCval_results[['WHOLE']][[l]] <- list()
    ROCval_results[['WHOLE']][[l]] <- get( paste0("AUCval_WHOLE_",l) )
    rm(list = paste0("AUCval_WHOLE_",l) )
    for(j in strata_names){
      ROCval_results[[j]][[l]] <- get(paste0("AUCval_", j, "_", l)) 
      rm(list = paste0("AUCval_", j, "_", l))
    }
    
  }
  
  print(Sys.time()-start)
  
  
  # save.image(paste0("~/Documents/Rstudio/Simulations/Simulations mai 2025/Résultats/",length(iterations),"ite_",N,"ind_",date_launch,".Rdata"))
  save(list = ls(), file = paste0("~/Documents/Rstudio/Simulations/Simulations mai 2025/Résultats/N",N,"_results/",PH_val,"/",length(iterations),"ite_",N,"ind_",PH_val,"_",cens_val,"_",date_launch,".Rdata"))
  
}


iterations <- 1:1000

N = 1000
path0 <- paste0("~/Documents/Rstudio/Simulations/Simulations mai 2025/Résultats/N",N,"_results/",PH_val,"/Output simulations_N",N,"_",PH_val,"_",cens_val,"/")
path1 <- "~/Documents/Rstudio/Simulations/BASES/"

rm(list = setdiff(ls(envir = .GlobalEnv), 
                  c("path0", "date_launch", "fr.ratetable", "slopop", 
                    "colrec", "Sn", "simulate_iteration", "calc_indic", "calc_indic_PP", "cens_val", "PH_val")), envir = .GlobalEnv) ## cleaning de l'environnement excpeté ce qu iva être réutilisé
## 1000SiEND

############################

## 1000IndSTA
# path0 <- "~/Documents/Rstudio/Simulations/Simulations mai 2025/Résultats/N1000_results/NPH/Output Simulations_N1000_NPH_HC/"

indic = c()
for(i in 1:1000){
  if(!file.exists(paste0(path0,"DATAFRAMES/TRAIN/WHOLE/",i,"_dftrain.csv"))){ indic <<- c(indic, i)} 
}

if(!is.null(indic)){
  iterations <- c(1:1000)[-indic]
}else{
  iterations <- c(1:1000)
}


calc_indic(N = 1000, iterations = iterations)

rm(list = setdiff(ls(envir = .GlobalEnv), 
                  c("path0", "date_launch", "fr.ratetable", "slopop", 
                    "colrec", "Sn", "simulate_iteration", "calc_indic", "calc_indic_PP", "iterations", "cens_val", "PH_val")), envir = .GlobalEnv)


##############################################################################################################################


rm(list = ls()) 


##############################################################################################################################


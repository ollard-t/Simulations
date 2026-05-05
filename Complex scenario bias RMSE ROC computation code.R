path0 <- paste0("/home/thomas/Documents/Rstudio/Simulations/Simulations janvier 2026/Résultatscomplex/Output simulations_N1000_COM/")
path0 <- paste0("/home/thomas/Documents/Rstudio/Simulations/Simulations mars 2026/Résultatscomplex/Output simulations_N1000_COM/")


calc_indic_PP <- function(iterations){
  
  strata_names <- c("H", "F")
  newtimes <- c(365.241, 365.241*3, 365.241*5, 365.241*10) 
  
  ### Calculs des indicateurs
  
  ### RMSE & BIAIS
  ############################################ TRAIN #################################################
  #WHOLE
  ALL <- vector("list", length(iterations))
  
  for (i in seq_along(iterations)) {
    k <- iterations[i]
    ALL[[i]] <- read.csv(
      paste0(path0, "TRAIN/WHOLE/mean/pred_mat_", k, ".csv"),
      sep = ""
    )
  }
  
  diff_list <- lapply(ALL, function(mat) {
    sapply(c("PLANN", "FLEX1", "FLEX2"), function(m) {
      mat[m, ] - mat["PP", ]
    })
  })
  diff_list <- lapply(diff_list, t)
  
  diff_array <- simplify2array(diff_list)
  
  mean_diff <- apply(
    diff_array,
    c(1, 2),
    function(x) mean(as.numeric(x), na.rm = TRUE)
  )  
  
  bias <- mean_diff
  ### RMSE
  
  sq_list <- lapply(ALL, function(mat) {
    sapply(c("PLANN", "FLEX1", "FLEX2"), function(m) {
      (mat[m, ] - mat["PP", ])^2
    })
  })
  sq_list <- lapply(sq_list, t)
  
  sq_array <- simplify2array(sq_list)
  
  sq_mean <- apply(
    sq_array,
    c(1, 2),
    function(x) mean(as.numeric(x), na.rm = TRUE)
  )  
  
  RMSE <- sqrt(sq_mean)
  
  #STRATES
  
  listbias <- vector("list", length(strata_names)+1)
  listRMSE <- vector("list", length(strata_names)+1)
  
  names(listbias) <- c("WHOLE", strata_names)
  names(listRMSE) <- c("WHOLE", strata_names)
  
  listbias[[1]] <- bias
  listRMSE[[1]] <- RMSE
  
  for(j in strata_names){
    
    ALL_hold <- vector("list", length(iterations))
    
    for (i in seq_along(iterations)) {
      k <- iterations[i]
      ALL_hold[[i]] <- read.csv(paste0(path0, "TRAIN/STRATA/mean/pred_mat_",k,"_",j,".csv"), sep = "")
    }
    
    diff_list <- lapply(ALL_hold, function(mat) {
      sapply(c("PLANN", "FLEX1", "FLEX2"), function(m) {
        mat[m, ] - mat["PP", ]
      })
    })
    diff_list <- lapply(diff_list, t)
    
    diff_array <- simplify2array(diff_list)
    
    mean_diff <- apply(
      diff_array,
      c(1, 2),
      function(x) mean(as.numeric(x), na.rm = TRUE)
    )  
    
    biasS <- mean_diff
    ### RMSE
    
    sq_list <- lapply(ALL_hold, function(mat) {
      sapply(c("PLANN", "FLEX1", "FLEX2"), function(m) {
        (mat[m, ] - mat["PP", ])^2
      })
    })
    sq_list <- lapply(sq_list, t)
    
    sq_array <- simplify2array(sq_list)
    
    sq_mean <- apply(
      sq_array,
      c(1, 2),
      function(x) mean(as.numeric(x), na.rm = TRUE)
    )  
    
    RMSES <- sqrt(sq_mean)
    assign(paste0("bias_",j), biasS)
    assign(paste0("RMSE_",j), RMSES)
    
    listbias[[j]] <- biasS
    listRMSE[[j]] <- RMSES
    
  }
  
  ### RMSE & BIAIS
  ############################################ VALID #################################################
  #WHOLE
  ALLval <- vector("list", length(iterations))
  
  for (i in seq_along(iterations)) {
    k <- iterations[i]
    ALLval[[i]] <- read.csv(
      paste0(path0, "VALID/WHOLE/mean/pred_mat_val_", k, ".csv"),
      sep = ""
    )
  }
  
  diff_list <- lapply(ALLval, function(mat) {
    sapply(c("PLANN", "FLEX1", "FLEX2"), function(m) {
      mat[m, ] - mat["PP", ]
    })
  })
  diff_list <- lapply(diff_list, t)
  
  diff_array <- simplify2array(diff_list)
  
  mean_diff <- apply(
    diff_array,
    c(1, 2),
    function(x) mean(as.numeric(x), na.rm = TRUE)
  )  
  
  bias_val <- mean_diff
  ### RMSE
  
  sq_list <- lapply(ALLval, function(mat) {
    sapply(c("PLANN", "FLEX1", "FLEX2"), function(m) {
      (mat[m, ] - mat["PP", ])^2
    })
  })
  sq_list <- lapply(sq_list, t)
  
  sq_array <- simplify2array(sq_list)
  
  sq_mean <- apply(
    sq_array,
    c(1, 2),
    function(x) mean(as.numeric(x), na.rm = TRUE)
  )  
  
  RMSEval <- sqrt(sq_mean)
  
  #STRATES
  
  listbiasval <- vector("list", length(strata_names)+1)
  listRMSEval <- vector("list", length(strata_names)+1)
  
  names(listbiasval) <- c("WHOLE", strata_names)
  names(listRMSEval) <- c("WHOLE", strata_names)
  
  listbiasval[[1]] <- bias_val
  listRMSEval[[1]] <- RMSEval
  
  for(j in strata_names){
    
    ALL_hold <- vector("list", length(iterations))
    
    for (i in seq_along(iterations)) {
      k <- iterations[i]
      ALL_hold[[i]] <- read.csv(paste0(path0, "VALID/STRATA/mean/pred_mat_val_",k,"_",j,".csv"), sep = "")
    }
    
    diff_list <- lapply(ALL_hold, function(mat) {
      sapply(c("PLANN", "FLEX1", "FLEX2"), function(m) {
        mat[m, ] - mat["PP", ]
      })
    })
    diff_list <- lapply(diff_list, t)
    
    diff_array <- simplify2array(diff_list)
    
    mean_diff <- apply(
      diff_array,
      c(1, 2),
      function(x) mean(as.numeric(x), na.rm = TRUE)
    )  
    
    biasS <- mean_diff
    ### RMSE
    
    sq_list <- lapply(ALL_hold, function(mat) {
      sapply(c("PLANN", "FLEX1", "FLEX2"), function(m) {
        (mat[m, ] - mat["PP", ])^2
      })
    })
    sq_list <- lapply(sq_list, t)
    
    sq_array <- simplify2array(sq_list)
    
    sq_mean <- apply(
      sq_array,
      c(1, 2),
      function(x) mean(as.numeric(x), na.rm = TRUE)
    )  
    
    RMSES <- sqrt(sq_mean)
    assign(paste0("bias_",j), biasS)
    assign(paste0("RMSE_",j), RMSES)
    
    listbiasval[[j]] <- biasS
    listRMSEval[[j]] <- RMSES
    
  }
  date_launch <- Sys.Date()
  save(list = ls(), file =paste0("~/Documents/Rstudio/Simulations/Simulations mars 2026/Résultatscomplex/",length(iterations),"ite_complex_",date_launch,"_PP.Rdata"))
  # save(list = ls(), file = paste0("~/Documents/Rstudio/Simulations/Simulations janvier 2026/Résultatscomplex/N",N,"_results/",PH_val,"/",length(iterations),"ite_",N,"ind_",PH_val,"_",cens_val,"_",date_launch,"_PP.Rdata"))
  
}
calc_ROC <- function(iterations){
  
  strata_names <- c("H", "F")
  newtimes <- c(365.241, 365.241*3, 365.241*5, 365.241*10) 
  
  start <- Sys.time()
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
                                              sep = " ") )
    
    ##survies individuelles 
    #plann
    assign(paste0("Pind_",k), read.csv(paste0(path0, "TRAIN/WHOLE/ind/plann_",k,".csv"),
                                       sep = " ")) 
    assign(paste0("F1ind_",k), read.csv(paste0(path0, "TRAIN/WHOLE/ind/flex1_",k,".csv"),
                                        sep = " ") )
    assign(paste0("F2ind_",k),  read.csv(paste0(path0, "TRAIN/WHOLE/ind/flex2_",k,".csv"),
                                         sep = " ") )
  }
  
  ##### calcul de ROC net
  
  AUC_calc_W <- function(k) {
    datatrain <- get(paste0("DATATrain_", k))
    Pind <- get(paste0("Pind_", k))
    F1ind <- get(paste0("F1ind_", k))
    F2ind <- get(paste0("F2ind_", k))
    
    rm(list = c(paste0("DATATrain_",k), paste0("Pind_",k), paste0("F1ind_",k),
                paste0("F2ind_",k)) ) 
    
    # Initialiser les vecteurs de résultats
    hold_P <- hold_F1 <- hold_F2 <- c()
    # sort(unique(1-ind_estimP[, l + 1]))
    # unique( quantile( 1-ind_estimP[, l + 1],  probs= seq(0,1, by = 0.01)) ) #pas en dessous de 20 points 0.025
    
    for (l in seq_along(newtimes)) {
      hold_P    <- c(hold_P,    roc.net(datatrain$times, datatrain$status, 1-Pind[,l], datatrain$age, datatrain$sexchara,
                                        datatrain$year, slopop, pro.time = newtimes[l], 
                                        cut.off = unique( quantile(1-Pind[,l],probs= seq(0,1, by = 0.01)) ) )$auc)
      hold_F1<- c(hold_F1, roc.net(datatrain$times, datatrain$status, 1-F1ind[, l], datatrain$age, datatrain$sexchara,
                                   datatrain$year, slopop, pro.time = newtimes[l],
                                   cut.off = unique( quantile( 1-F1ind[, l],  probs= seq(0,1, by = 0.01)) ))$auc)
      hold_F2 <- c(hold_F2, roc.net(datatrain$times, datatrain$status, 1-F2ind[, l], datatrain$age, datatrain$sexchara,
                                    datatrain$year, slopop, pro.time = newtimes[l], 
                                    cut.off = unique( quantile( 1-F2ind[, l],  probs= seq(0,1, by = 0.01)) ))$auc)
    }
    
    return(list(
      P     = hold_P,
      F1  = hold_F1,
      F2  = hold_F2
    ))
  }
  
  results_list <- mclapply(iterations, AUC_calc_W, mc.cores = detectCores() - 4)
  
  for(k in iterations){
    rm(list = c(
      paste0("DATATrain_", k),
      paste0("Pind_",  k),
      paste0("F1ind_", k),
      paste0("F2ind_", k)
    ),
    envir = .GlobalEnv)
  }
  
  # Transformer en matrices pour chaque méthode :
  AUC_WHOLE_P    <- do.call(rbind, lapply(results_list, function(x) x$P))
  AUC_WHOLE_F1 <- do.call(rbind, lapply(results_list, function(x) x$F1))
  AUC_WHOLE_F2 <- do.call(rbind, lapply(results_list, function(x) x$F2))
  
  
  ### pour renommer les colonnes avec les temps de prognostic
  for(l in c("P","F1", "F2")){
    dathold <- get(paste0("AUC_WHOLE_", l))
    colnames(dathold) <- c("1 year", "3 years", "5 years", "10 years")
    rownames(dathold) <- c(iterations)
    assign(paste0("AUC_WHOLE_", l),dathold)
  }
  ## moyennes par temps 
  
  for(l in c("P","F1","F2")){
    
    assign(paste0("AUC_means_",l) ,  colMeans( get(paste0("AUC_WHOLE_",l)) ) )
    
  }
  ################## STRATES
  
  ##importation des bases de données / survies individuelles prédites par les modèles 
  
  for(k in iterations){
    ##bases de données 
    for (j in strata_names){
      assign(paste0("DATATrain_",j,"_", k) , read.csv(paste0(path0, "DATAFRAMES/TRAIN/",j,"/",k,"_dftrain",j,".csv"),
                                                      sep = " ") )
      ##survies individuelles 
      #plann
      
      assign(paste0("Pind_",j,"_", k), read.csv(paste0(path0, "TRAIN/STRATA/ind/plann_",k,"_",j,".csv"),
                                                sep = ";")) 
      assign(paste0("F1ind_",j,"_", k), read.csv(paste0(path0, "TRAIN/STRATA/ind/flex1_",k,"_",j,".csv"),
                                                 sep = " ") )
      assign(paste0("F2ind_",j,"_", k),  read.csv(paste0(path0, "TRAIN/STRATA/ind/flex2_",k,"_",j,".csv"),
                                                  sep = " ") )
      
    }
  }
  
  ##### calcul de ROC net
  
  AUC_calc_S <- function(k) {
    auc_results <- list()
    
    for (j in strata_names) {
      datatrain <- get(paste0("DATATrain_", j, "_", k))
      Pind <- get(paste0("Pind_",j,"_", k))
      F1ind<- get(paste0("F1ind_",j,"_", k))
      F2ind<- get(paste0("F2ind_",j,"_", k))
      
      rm(list = c(paste0("DATATrain_",k), paste0("Pind_",j,"_", k), paste0("F1ind_",j,"_", k),
                  paste0("F2ind_",j,"_", k)) ) 
      
      # Initialiser les vecteurs de résultats
      hold_P <- hold_F1 <- hold_F2 <- c()
      # sort(unique(1-ind_estimP[, l + 1]))
      # unique( quantile( 1-ind_estimP[, l + 1],  probs= seq(0,1, by = 0.01)) ) #pas en dessous de 20 points 0.025
      
      for (l in seq_along(newtimes)) {
        hold_P    <- c(hold_P,    roc.net(datatrain$times, datatrain$status, 1-Pind[,l], datatrain$age, datatrain$sexchara,
                                          datatrain$year, slopop, pro.time = newtimes[l], 
                                          cut.off = unique( quantile(1-Pind[,l],probs= seq(0,1, by = 0.01)) ) )$auc)
        hold_F1<- c(hold_F1, roc.net(datatrain$times, datatrain$status, 1-F1ind[, l], datatrain$age, datatrain$sexchara,
                                     datatrain$year, slopop, pro.time = newtimes[l],
                                     cut.off = unique( quantile( 1-F1ind[, l],  probs= seq(0,1, by = 0.01)) ))$auc)
        hold_F2 <- c(hold_F2, roc.net(datatrain$times, datatrain$status, 1-F2ind[, l], datatrain$age, datatrain$sexchara,
                                      datatrain$year, slopop, pro.time = newtimes[l], 
                                      cut.off = unique( quantile( 1-F2ind[, l],  probs= seq(0,1, by = 0.01)) ))$auc)
      }
      
      auc_results[[j]] <- list(
        k = k,
        P = hold_P,
        F1 = hold_F1,
        F2 = hold_F2
        
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
        paste0("Pind_", j, "_", k),
        paste0("F1ind_", j, "_", k),
        paste0("F2ind_", j, "_", k)
      ),
      envir = .GlobalEnv)
    }
  }
  
  # Transformer en matrices pour chaque méthode :
  
  for(j in strata_names){
    assign(paste0("AUC_",j,"_P"), do.call(rbind, lapply(res_all, function(x) x[[j]]$P)) )
    assign(paste0("AUC_",j,"_F1"), do.call(rbind, lapply(res_all, function(x) x[[j]]$F1)) )
    assign(paste0("AUC_",j,"_F2"), do.call(rbind, lapply(res_all, function(x) x[[j]]$F2)) )
  }
  ### pour renommer les colonnes avec les temps de prognostic
  for(l in c("P","F1","F2")){
    for(j in strata_names){
      dathold <- get(paste0("AUC_", j,"_",l))
      colnames(dathold) <- c("1 year", "3 years", "5 years", "10 years")
      rownames(dathold) <- c(iterations)
      assign(paste0("AUC_",j,"_", l),dathold)
    }
  }
  ## moyennes par temps 
  
  for(l in c("P","F1","F2")){
    for(j in strata_names){
      assign(paste0("AUC_means_",j,"_",l) ,  colMeans( get(paste0("AUC_",j,"_",l)) ) )
    }
  }
  
  
  ### sous forme de liste pour libérer de l'espace dans l'evt (pour les moyennes)
  
  ROCmean_results <- list()
  
  for(l in c("P","F1","F2")){
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
  
  for(l in c("P","F1","F2")){
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
    assign(paste0("DATATrain_", k) , read.csv(paste0(path0, "DATAFRAMES/VALID/WHOLE/",k,"_dfvalid.csv"),
                                              sep = " ") )
    
    ##survies individuelles 
    #plann
    assign(paste0("Pind_",k), read.csv(paste0(path0, "VALID/WHOLE/ind/plann_val_",k,".csv"),
                                       sep = " ")) 
    assign(paste0("F1ind_",k), read.csv(paste0(path0, "VALID/WHOLE/ind/flex1_val_",k,".csv"),
                                        sep = " ") )
    assign(paste0("F2ind_",k),  read.csv(paste0(path0, "VALID/WHOLE/ind/flex2_val_",k,".csv"),
                                         sep = " ") )
  }
  
  ##### calcul de ROC net
  
  ##creation d'objets pour enregistrer les AUC de chaque itération aux différents temps et pour les différents modèles 
  
  AUC_calc_val_W <- function(k){ 
    datatrain <- get(paste0("DATATrain_", k))
    Pind <- get(paste0("Pind_", k))
    F1ind <- get(paste0("F1ind_", k))
    F2ind <- get(paste0("F2ind_", k))
    
    rm(list = c(paste0("DATATrain_",k), paste0("Pind_",k), paste0("F1ind_",k),
                paste0("F2ind_",k)) ) 
    
    # Initialiser les vecteurs de résultats
    hold_P <- hold_F1 <- hold_F2 <- c()
    # sort(unique(1-ind_estimP[, l + 1]))
    # unique( quantile( 1-ind_estimP[, l + 1],  probs= seq(0,1, by = 0.01)) ) #pas en dessous de 20 points 0.025
    
    for (l in seq_along(newtimes)) {
      hold_P    <- c(hold_P,    roc.net(datatrain$times, datatrain$status, 1-Pind[,l], datatrain$age, datatrain$sexchara,
                                        datatrain$year, slopop, pro.time = newtimes[l], 
                                        cut.off = unique( quantile(1-Pind[,l],probs= seq(0,1, by = 0.01)) ) )$auc)
      hold_F1<- c(hold_F1, roc.net(datatrain$times, datatrain$status, 1-F1ind[, l], datatrain$age, datatrain$sexchara,
                                   datatrain$year, slopop, pro.time = newtimes[l],
                                   cut.off = unique( quantile( 1-F1ind[, l],  probs= seq(0,1, by = 0.01)) ))$auc)
      hold_F2 <- c(hold_F2, roc.net(datatrain$times, datatrain$status, 1-F2ind[, l], datatrain$age, datatrain$sexchara,
                                    datatrain$year, slopop, pro.time = newtimes[l], 
                                    cut.off = unique( quantile( 1-F2ind[, l],  probs= seq(0,1, by = 0.01)) ))$auc)
    }
    
    return(list(
      P     = hold_P,
      F1  = hold_F1,
      F2  = hold_F2
    ))
  }
  
  resultval_list <- mclapply(iterations, AUC_calc_val_W, mc.cores = detectCores() - 4)
  
  for(k in iterations){
    rm(list = c(
      paste0("DATATrain_", k),
      paste0("Pind_",  k),
      paste0("F1ind_", k),
      paste0("F2ind_", k)
    ),
    envir = .GlobalEnv)
  }
  # Transformer en matrices pour chaque méthode :
  AUCval_WHOLE_P    <- do.call(rbind, lapply(resultval_list, function(x) x$P))
  AUCval_WHOLE_F1 <- do.call(rbind, lapply(resultval_list, function(x) x$F1))
  AUCval_WHOLE_F2 <- do.call(rbind, lapply(resultval_list, function(x) x$F2))
  
  
  ### pour renommer les colonnes avec les temps de prognostic
  for(l in c("P","F1", "F2")){
    dathold <- get(paste0("AUCval_WHOLE_", l))
    colnames(dathold) <- c("1 year", "3 years", "5 years", "10 years")
    rownames(dathold) <- c(iterations)
    assign(paste0("AUCval_WHOLE_", l),dathold)
  }
  ## moyennes par temps 
  
  for(l in c("P","F1", "F2")){
    
    assign(paste0("AUCval_means_",l) ,  colMeans( get(paste0("AUCval_WHOLE_",l)) ) )
    
  }
  ################## STRATES
  
  ##importation des bases de données / survies individuelles prédites par les modèles 
  
  for(k in iterations){
    ##bases de données 
    
    for (j in strata_names){
      
      assign(paste0("DATATrain_",j,"_", k) , read.csv(paste0(path0, "DATAFRAMES/VALID/",j,"/",k,"_dfvalid",j,".csv"),
                                                      sep = " ") )
      ##survies individuelles 
      #plann
      
      assign(paste0("Pind_",j,"_", k), read.csv(paste0(path0, "VALID/STRATA/ind/plann_",k,"_",j,".csv"),
                                                sep = ";")) 
      assign(paste0("F1ind_",j,"_", k), read.csv(paste0(path0, "VALID/STRATA/ind/flex1_val_",k,"_",j,".csv"),
                                                 sep = " ") )
      assign(paste0("F2ind_",j,"_", k),  read.csv(paste0(path0, "VALID/STRATA/ind/flex2_val_",k,"_",j,".csv"),
                                                  sep = " ") )
      
    }
  }
  
  
  ##### calcul de ROC net
  
  AUC_calc_val_S <- function(k){ 
    auc_results <- list()
    
    for (j in strata_names) {
      datatrain <- get(paste0("DATATrain_", j, "_", k))
      Pind <- get(paste0("Pind_",j,"_", k))
      F1ind<- get(paste0("F1ind_",j,"_", k))
      F2ind<- get(paste0("F2ind_",j,"_", k))
      
      rm(list = c(paste0("DATATrain_", j, "_", k), paste0("Pind_",j,"_", k), paste0("F1ind_",j,"_", k),
                  paste0("F2ind_",j,"_", k)) ) 
      
      # Initialiser les vecteurs de résultats
      hold_P <- hold_F1 <- hold_F2 <- c()
      # sort(unique(1-ind_estimP[, l + 1]))
      # unique( quantile( 1-ind_estimP[, l + 1],  probs= seq(0,1, by = 0.01)) ) #pas en dessous de 20 points 0.025
      
      for (l in seq_along(newtimes)) {
        hold_P    <- c(hold_P,    roc.net(datatrain$times, datatrain$status, 1-Pind[,l], datatrain$age, datatrain$sexchara,
                                          datatrain$year, slopop, pro.time = newtimes[l], 
                                          cut.off = unique( quantile(1-Pind[,l],probs= seq(0,1, by = 0.01)) ) )$auc)
        hold_F1<- c(hold_F1, roc.net(datatrain$times, datatrain$status, 1-F1ind[, l], datatrain$age, datatrain$sexchara,
                                     datatrain$year, slopop, pro.time = newtimes[l],
                                     cut.off = unique( quantile( 1-F1ind[, l],  probs= seq(0,1, by = 0.01)) ))$auc)
        hold_F2 <- c(hold_F2, roc.net(datatrain$times, datatrain$status, 1-F2ind[, l], datatrain$age, datatrain$sexchara,
                                      datatrain$year, slopop, pro.time = newtimes[l], 
                                      cut.off = unique( quantile( 1-F2ind[, l],  probs= seq(0,1, by = 0.01)) ))$auc)
      }
      
      auc_results[[j]] <- list(
        k = k,
        P = hold_P,
        F1 = hold_F1,
        F2 = hold_F2
        
      )
    }
    
    return(auc_results)
  }
  
  resval_all <- mclapply(iterations, AUC_calc_val_S, mc.cores = detectCores() - 4)
  
  for(k in iterations){
    for(j in strata_names){
      rm(list = c(
        paste0("DATATrain_", j, "_", k),
        paste0("Pind_", j, "_", k),
        paste0("F1ind_", j, "_", k),
        paste0("F2ind_", j, "_", k)
      ),
      envir = .GlobalEnv)
    }
  }
  
  # Transformer en matrices pour chaque méthode :
  
  for(j in strata_names){
    assign(paste0("AUCval_",j,"_P"), do.call(rbind, lapply(resval_all, function(x) x[[j]]$P)) )
    assign(paste0("AUCval_",j,"_F1"), do.call(rbind, lapply(resval_all, function(x) x[[j]]$F1)) )
    assign(paste0("AUCval_",j,"_F2"), do.call(rbind, lapply(resval_all, function(x) x[[j]]$F2)) )
    
  }
  ### pour renommer les colonnes avec les temps de prognostic
  for(l in c("P","F1","F2")){
    for(j in strata_names){
      dathold <- get(paste0("AUCval_", j,"_",l))
      colnames(dathold) <- c("1 year", "3 years", "5 years", "10 years")
      rownames(dathold) <- c(iterations)
      assign(paste0("AUCval_",j,"_", l),dathold)
    }
  }
  ## moyennes par temps 
  
  for(l in c("P","F1","F2")){
    for(j in strata_names){
      assign(paste0("AUCval_means_",j,"_",l) ,  colMeans( get(paste0("AUCval_",j,"_",l)) ) )
    }
  }
  
  
  ### sous forme de liste pour libérer de l'espace dans l'evt (pour les moyennes)
  
  ROCmeanval_results <- list()
  
  for(l in c("P","F1","F2")){
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
  
  for(l in c("P","F1","F2")){
    ROCval_results[['WHOLE']][[l]] <- list()
    ROCval_results[['WHOLE']][[l]] <- get( paste0("AUCval_WHOLE_",l) )
    rm(list = paste0("AUCval_WHOLE_",l) )
    for(j in strata_names){
      ROCval_results[[j]][[l]] <- get(paste0("AUCval_", j, "_", l)) 
      rm(list = paste0("AUCval_", j, "_", l))
    }
    
  }
  
  print(Sys.time()-start)
  
  
  # save.image(paste0("~/Documents/Rstudio/Simulations/Simulations janvier 2026/RésultatsX/",length(iterations),"ite_",N,"ind_",date_launch,".Rdata"))
  save(list = ls(), file = paste0("~/Documents/Rstudio/Simulations/Simulations mars 2026/Résultatscomplex/",length(iterations),"ite_ROC.Rdata"))
  
}



indic = c()
for(i in 1:1000){
  if(!file.exists(paste0(path0,"TRAIN/WHOLE/mean/pred_mat_",i,".csv"))){ indic <<- c(indic, i)} 
}

if(!is.null(indic)){
  iterations <- c(1:1000)[-indic]
}else{
  iterations <- c(1:1000)
}

calc_indic_PP(iterations)

calc_ROC(iterations)
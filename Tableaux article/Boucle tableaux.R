library(xtable)

latex_table <- function(N, strata, file_name, cens_val, PH_val, PP = FALSE){
  path1 <- paste0("~/Documents/Rstudio/Simulations/Simulations mai 2025/Résultats/N",N,"_results/",PH_val,"/",file_name,".Rdata")#nom à changer
  load(path1)
  if(PP == TRUE){
    rm(list = setdiff(ls(), c("ROCmean_results", "ROCmeanval_results", "N", "strata", "file_name", "cens_val", "PH_val", "PP") ) )
    path1 <- paste0("~/Documents/Rstudio/Simulations/Simulations mai 2025/Résultats/N",N,"_results/",PH_val,"/",file_name,"_PP.Rdata")
    load(path1)
  }
  ##pour aller récupérer les logliklihood
  
  path0 <- paste0("~/Documents/Rstudio/Simulations/Simulations mai 2025/Résultats/N",N,"_results/",PH_val,"/Output simulations_N",N,"_",PH_val,"_",cens_val,"/LOGLIK/")
  
  method_names <- c("WG", "FLEX1.2", "FLEX1.4", "FLEX2.2", "FLEX2.4",
                    "PLANN", "PP")
  method_names_val <- c("WGval", "FLEX1.2val", "FLEX1.4val", "FLEX2.2val", "FLEX2.4val",
                        "PLANNval", "PPval")
  article_names <- c("&Weibull", "~ & ~ &Flexnet1.2","~ & ~ &Flexnet1.4","~ & ~ &Flexnet2.2","~ & ~ &Flexnet2.4", "~ & ~ &PLANN", "~ & ~ &Pohar Perme")
  
  if(PP == TRUE){
    method_names <- method_names[-7]
    method_names_val <- method_names_val[-7]
    article_names <- article_names[-7]
  }

  if(PP == FALSE){
    colnames(biais_PP) <- colnames(biais_FLEX1.2)
  
    colnames(biais_PPval) <- colnames(biais_FLEX1.2)
  
    colnames(RMSE_PP) <- colnames(RMSE_FLEX1.2)
    colnames(RMSE_PPval) <- colnames(RMSE_FLEX1.2)
  }
  
  whole_biais <- data.frame()
  for(i in method_names){
    a <- get(paste0("biais_",i))
    whole_biais <- rbind(whole_biais, a[strata,])
  }
  rownames(whole_biais) <- article_names
  
  whole_RMSE <- data.frame()
  for(i in method_names){
    a <- get(paste0("RMSE_",i))
    whole_RMSE <- rbind(whole_RMSE, a[strata,])
  }
  rownames(whole_RMSE) <- article_names
  
  ## roc
  whole_roc <- data.frame()
  for(i in 1:6){
    whole_roc <- rbind(whole_roc, unlist(unname(ROCmean_results[[strata]][c(6,2:5,1)][i])))
  }
  if(PP == FALSE){
    whole_roc <- rbind(whole_roc, c(rep(NA,4)))
  }
  rownames(whole_roc) <- article_names
  colnames(whole_roc) <- colnames(whole_RMSE)
  
  ##loglik
  iterations <- 1:1000
  miss <- c()
  for(i in 1:1000){
    if(!file.exists(paste0(path0,"TRAIN/",i,"_loglik.csv"))){miss <- c(miss, i)}
  }
  
  if(!is.null(miss)){
    iterations <- iterations[-miss]
  }
  
  logll <- unname(read.csv(paste0(path0,"TRAIN/",iterations[1],"_loglik.csv"), sep= " ", quote = "") )[c(6,2:5,1),strata]
  for(i in iterations[-1]){
    # assign(paste0("a_", i) , unname(read.csv(paste0(path0,i,"_loglik.csv")))[c(2:5,1),] )
    a <- unname(read.csv(paste0(path0,"TRAIN/",i,"_loglik.csv"), sep= " ", quote = "") )[c(6,2:5,1),strata]
    logll <- cbind(logll, a)
  }
  colnames(logll) <- NULL
  
  whole_log <- as.data.frame( t(t(apply(logll, mean, MARGIN = 1))) )
  if(PP == FALSE){
    whole_log <-rbind(whole_log, NA)
  }
  rownames(whole_log) <- article_names
  
  whole_tout <- cbind(whole_biais, whole_RMSE, whole_roc, whole_log)
  
  xt <- xtable(whole_tout, digits = c(0,4,4,4,4,4,4,4,4,4,4,4,4,4))
  
  ####### VALID
  
  whole_biais_val <- data.frame()
  for(i in method_names_val){
    a <- get(paste0("biais_",i))
    whole_biais_val <- rbind(whole_biais_val, a[strata,])
  }
  rownames(whole_biais_val) <- article_names
  
  whole_RMSE_val <- data.frame()
  for(i in method_names_val){
    a <- get(paste0("RMSE_",i))
    whole_RMSE_val <- rbind(whole_RMSE_val, a[strata,])
  }
  rownames(whole_RMSE_val) <- article_names
  
  ## roc
  whole_roc_val <- data.frame()
  for(i in 1:6){
    whole_roc_val <- rbind(whole_roc_val, unlist(unname(ROCmeanval_results[[strata]][c(6,2:5,1)][i])))
  }
  
  if(PP == FALSE){
    whole_roc_val <- rbind(whole_roc_val, c(rep(NA,4)))
  }
  rownames(whole_roc_val) <- article_names
  colnames(whole_roc_val) <- colnames(whole_RMSE_val)
  
  ##loglik
  iterations <- 1:1000
  miss <- c()
  for(i in 1:1000){
    if(!file.exists(paste0(path0,"VALID/",i,"_loglikval.csv"))){miss <- c(miss, i)}
  }
  
  if(!is.null(miss)){
    iterations <- iterations[-miss]
  }
  
  logllval <- unname(read.csv(paste0(path0,"VALID/",iterations[1],"_loglikval.csv"), sep= " ", quote = "") )[c(6,2:5,1),strata]
  for(i in iterations[-1]){
    # assign(paste0("a_", i) , unname(read.csv(paste0(path0,i,"_loglik.csv")))[c(2:5,1),] )
    a <- unname(read.csv(paste0(path0,"VALID/",i,"_loglikval.csv"), sep= " ", quote = "") )[c(6,2:5,1),strata]
    logllval <- cbind(logllval, a)
  }
  colnames(logllval) <- NULL
  
  PLANNInfcount <- which(apply(logllval, 2, function(x) any(is.infinite(x))))
  F2.2count <- which(is.na(logllval[4,]))
  F2.4count <- which(is.na(logllval[5,]))
  
  mat <- as.matrix(logllval)
  mat[mat == -Inf] <- NA
  logllval <- as.data.frame(mat)
  
  whole_log_val <- as.data.frame( t(t(apply(logllval, mean, na.rm = TRUE, MARGIN = 1))) )
  
  if(PP == FALSE){
    whole_log_val <-rbind(whole_log_val, NA)
  }
  rownames(whole_log_val) <- article_names
  
  whole_tout_val <- cbind(whole_biais_val, whole_RMSE_val, whole_roc_val, whole_log_val)
  
  xt_val <- xtable(whole_tout_val[-7,], digits = c(0,4,4,4,4,4,4,4,4,4,4,4,4,4))
  
  path3 <- paste0("~/Documents/Rstudio/Simulations/Simulations mai 2025/Résultats/N",N,"_results/",PH_val,"/Output simulations_N",N,"_",PH_val,"_",cens_val,"/")
  
  if(!(strata == 1)){
    if(strata == 2){name <- "HC"}
    if(strata == 3){name <- "HR"}
    if(strata == 4){name <- "FC"}
    if(strata == 5){name <- "FR"}
    prop <- c()
    for(i in iterations){
      a <- read.csv(paste0(path3,"DATAFRAMES/TRAIN/",name,"/",i,"_dftrain_",name,".csv"), sep = ";")
      prop <- c(prop, dim(a)[1])
    }
    print(paste0("Taille sous-groupe : ", mean(prop)/(N/2)))
  }
  print(paste0("Nombre Inf PLANN : ", length(PLANNInfcount),". Détail : ", paste(PLANNInfcount, collapse = ", ")))
  print(paste0("Nombre NA FLEX2.2 : ", length(F2.2count),". Détail : ", paste(F2.2count, collapse = ", ")))
  print(paste0("Nombre NA FLEX2.4 : ", length(F2.4count),". Détail : ", paste(F2.4count, collapse = ", ")))
  
  
  print(xt, type = "latex", sanitize.rownames.function = identity)
  print(xt_val, type = "latex", sanitize.rownames.function = identity)
}

## strata
# 1 = WHOLE
# 2 = HC
# 3 = HR
# 4 = FC
# 5 = FR

file_name = "828ite_1000ind_2025-09-04"
cens_val = "HC"
PH_val = "NPH"
      ###############
      #     1000
      ###############
for(i in 1:5){
  latex_table(1000, i, file_name, cens_val, PH_val)
}


 ###############
#     3000
###############
for(i in 1:5){
  latex_table(3000, i, file_name, cens_val, PH_val, PP =TRUE)
}

###############
#     5000
###############
for(i in 1:5){
  latex_table(5000, i, file_name, cens_val, PH_val)
}



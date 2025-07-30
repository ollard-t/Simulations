library(xtable)

latex_table <- function(N, strata, file_name){
  path1 <- paste0("~/Documents/Rstudio/Simulations/Simulations mai 2025/Résultats/N",N,"_results/PH/",file_name,".Rdata")#nom à changer
  load(path1)
  
  ##pour aller récupérer les logliklihood
  
  path0 <- paste0("~/Documents/Rstudio/Simulations/Simulations mai 2025/N",N,"/Output simulations_N",N,"_PH_HC/LOGLIK/")
  
  method_names <- c("WG", "FLEX1.2", "FLEX1.4", "FLEX2.2", "FLEX2.4",
                    "PLANN", "PP")
  method_names_val <- c("WGval", "FLEX1.2val", "FLEX1.4val", "FLEX2.2val", "FLEX2.4val",
                        "PLANNval", "PPval")
  
  article_names <- c("Weibull", "Flexnet1.2","Flexnet1.4","Flexnet2.2","Flexnet2.4", "PLANN", "Pohar Perme")
  
  colnames(biais_PP) <- colnames(biais_FLEX1.2)
  colnames(biais_PPval) <- colnames(biais_FLEX1.2)
  
  colnames(RMSE_PP) <- colnames(RMSE_FLEX1.2)
  colnames(RMSE_PPval) <- colnames(RMSE_FLEX1.2)
  
  
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
    whole_roc <- rbind(whole_roc, unlist(unname(ROCmean_results[[1]][c(6, 1:5)][i])))
  }
  
  whole_roc <- rbind(whole_roc, c(rep(NA,4)))
  rownames(whole_roc) <- article_names
  colnames(whole_roc) <- colnames(whole_RMSE)
  
  ##loglik
  iterations <- 1:1000
  miss <- c()
  for(i in 1:1000){
    if(!file.exists(paste0(path0,i,"_loglik.csv"))){miss <- c(miss, i)}
  }
  
  if(!is.null(miss)){
    iterations <- iterations[-miss]
  }
  
  logll <- unname(read.csv(paste0(path0,iterations[1],"_loglik.csv"), sep= " ", quote = "") )[c(6,2:5,1),2]
  for(i in iterations[-1]){
    # assign(paste0("a_", i) , unname(read.csv(paste0(path0,i,"_loglik.csv")))[c(2:5,1),] )
    a <- unname(read.csv(paste0(path0,i,"_loglik.csv"), sep= " ", quote = "") )[c(6,2:5,1),2]
    logll <- cbind(logll, a)
  }
  colnames(logll) <- NULL
  
  whole_log <- as.data.frame( t(t(apply(logll, mean, MARGIN = 1))) )
  whole_log <-rbind(whole_log, NA)
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
    whole_roc_val <- rbind(whole_roc_val, unlist(unname(ROCmeanval_results[[1]][c(6,1:5)][i])))
  }
  
  whole_roc_val <- rbind(whole_roc_val, c(rep(NA,4)))
  rownames(whole_roc_val) <- article_names
  colnames(whole_roc_val) <- colnames(whole_RMSE_val)
  
  
  whole_tout_val <- cbind(whole_biais_val, whole_RMSE_val, whole_roc_val, whole_log)
  
  xt_val <- xtable(whole_tout_val, digits = c(0,4,4,4,4,4,4,4,4,4,4,4,4,4))
  
  path3 <- paste0("~/Documents/Rstudio/Simulations/Simulations mai 2025/N",N,"/Output simulations_N",N,"_PH_HC/")
  
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
    print(mean(prop)/(N/2))
  }
  print(xt, type = "latex")
  print(xt_val, type = "latex")
}

      ###############
      #     1000
      ###############
#1 = WHOLE
latex_table(1000, 1, file_name)
#2 = HC
latex_table(1000, 2, file_name)
#3 = HR
latex_table(1000, 3, file_name)
#4 = FC
latex_table(1000, 4, file_name)
#5 = FR
latex_table(1000, 5, file_name)


###############
#     3000
###############
#1 = WHOLE
latex_table(3000, 1, file_name)
#2 = HC
latex_table(3000, 2, file_name)
#3 = HR
latex_table(3000, 3, file_name)
#4 = FC
latex_table(3000, 4, file_name)
#5 = FR
latex_table(3000, 5, file_name)

###############
#     5000
###############
#1 = WHOLE
latex_table(5000, 1, file_name)
#2 = HC
latex_table(5000, 2, file_name)
#3 = HR
latex_table(5000, 3, file_name)
#4 = FC
latex_table(5000, 4, file_name)
#5 = FR
latex_table(5000, 5, file_name)



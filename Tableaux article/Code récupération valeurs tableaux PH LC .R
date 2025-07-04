library(xtable)

##################################################  PH LC
################################## N1000 
path1 <- "~/Documents/Rstudio/Simulations/Simulations mai 2025/Résultats/N1000_results/PH/999ite_1000ind_PH_LC_25_04.Rdata"
load(path1)

##pour aller récupérer les logliklihood

path0 <- "~/Documents/Rstudio/Simulations/Simulations mai 2025/N1000/Output simulations_N1000_PH_LC/LOGLIK/"

method_names <- c("FLEX1.2", "FLEX1.4", "FLEX2.2", "FLEX2.4",
                  "PLANN", "PP")
method_names_val <- c("FLEX1.2val", "FLEX1.4val", "FLEX2.2val", "FLEX2.4val",
                  "PLANNval", "PPval")

article_names <- c("Flexnet1.2","Flexnet1.4","Flexnet2.2","Flexnet2.4", "PLANN", "Pohar Perme")

colnames(biais_PP) <- colnames(biais_FLEX1.2)
colnames(biais_PPval) <- colnames(biais_FLEX1.2)

colnames(RMSE_PP) <- colnames(RMSE_FLEX1.2)
colnames(RMSE_PPval) <- colnames(RMSE_FLEX1.2)

whole_biais <- data.frame()
for(i in method_names){
  a <- get(paste0("biais_",i))
  whole_biais <- rbind(whole_biais, a[1,])
}
rownames(whole_biais) <- article_names

whole_RMSE <- data.frame()
for(i in method_names){
  a <- get(paste0("RMSE_",i))
  whole_RMSE <- rbind(whole_RMSE, a[1,])
}
rownames(whole_RMSE) <- article_names

## roc
whole_roc <- data.frame()
for(i in 1:6){
  whole_roc <- rbind(whole_roc, unlist(unname(ROCmean_results[[1]][c(2:6,1)][i])))
}
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

logll <- unname(read.csv(paste0(path0,iterations[1],"_loglik.csv"), sep= " ", quote = "") )[c(2:5,1),2]
for(i in iterations[-1]){
  # assign(paste0("a_", i) , unname(read.csv(paste0(path0,i,"_loglik.csv")))[c(2:5,1),] )
  a <- unname(read.csv(paste0(path0,i,"_loglik.csv"), sep= " ", quote = "") )[c(2:5,1),2]
  logll <- cbind(logll, a)
  }
colnames(logll) <- NULL

whole_log <- as.data.frame( t(t(apply(logll, mean, MARGIN = 1))) )
whole_log <-rbind(whole_log, NA)
rownames(whole_log) <- article_names

whole_tout <- cbind(whole_biais, whole_RMSE, whole_roc, whole_log)

xt <- xtable(whole_tout, digits = c(0,4,4,4,4,4,4,4,4,4,4,4,4,4))

# Affiche la version LaTeX
print(xt, type = "latex")

####### VALID

whole_biais_val <- data.frame()
for(i in method_names_val){
  a <- get(paste0("biais_",i))
  whole_biais_val <- rbind(whole_biais_val, a[1,])
}
rownames(whole_biais_val) <- article_names

whole_RMSE_val <- data.frame()
for(i in method_names_val){
  a <- get(paste0("RMSE_",i))
  whole_RMSE_val <- rbind(whole_RMSE_val, a[1,])
}
rownames(whole_RMSE_val) <- article_names

## roc
whole_roc_val <- data.frame()
for(i in 1:6){
  whole_roc_val <- rbind(whole_roc_val, unlist(unname(ROCmeanval_results[[1]][c(2:6,1)][i])))
}
rownames(whole_roc_val) <- article_names
colnames(whole_roc_val) <- colnames(whole_RMSE_val)


whole_tout_val <- cbind(whole_biais_val, whole_RMSE_val, whole_roc_val, whole_log)

xt_val <- xtable(whole_tout_val, digits = c(0,4,4,4,4,4,4,4,4,4,4,4,4,4))

# Affiche la version LaTeX
print(xt_val, type = "latex")

################################## N3000 
path1 <- "~/Documents/Rstudio/Simulations/Simulations mai 2025/Résultats/N3000_results/PH/1000ite_3000ind_PH_LC_25_04.Rdata"
load(path1)

##pour aller récupérer les logliklihood

path0 <- "~/Documents/Rstudio/Simulations/Simulations mai 2025/N3000/Output simulations_N3000_PH_LC/LOGLIK/"

method_names <- c("FLEX1.2", "FLEX1.4", "FLEX2.2", "FLEX2.4",
                  "PLANN", "PP")
method_names_val <- c("FLEX1.2val", "FLEX1.4val", "FLEX2.2val", "FLEX2.4val",
                      "PLANNval", "PPval")

article_names <- c("Flexnet1.2","Flexnet1.4","Flexnet2.2","Flexnet2.4", "PLANN", "Pohar Perme")

colnames(biais_PP) <- colnames(biais_FLEX1.2)
colnames(biais_PPval) <- colnames(biais_FLEX1.2)

colnames(RMSE_PP) <- colnames(RMSE_FLEX1.2)
colnames(RMSE_PPval) <- colnames(RMSE_FLEX1.2)

whole_biais <- data.frame()
for(i in method_names){
  a <- get(paste0("biais_",i))
  whole_biais <- rbind(whole_biais, a[1,])
}
rownames(whole_biais) <- article_names

whole_RMSE <- data.frame()
for(i in method_names){
  a <- get(paste0("RMSE_",i))
  whole_RMSE <- rbind(whole_RMSE, a[1,])
}
rownames(whole_RMSE) <- article_names

## roc
whole_roc <- data.frame()
for(i in 1:6){
  whole_roc <- rbind(whole_roc, unlist(unname(ROCmean_results[[1]][c(2:6,1)][i])))
}
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

logll <- unname(read.csv(paste0(path0,iterations[1],"_loglik.csv"), sep= " ", quote = "") )[c(2:5,1),2]
for(i in iterations[-1]){
  # assign(paste0("a_", i) , unname(read.csv(paste0(path0,i,"_loglik.csv")))[c(2:5,1),] )
  a <- unname(read.csv(paste0(path0,i,"_loglik.csv"), sep= " ", quote = "") )[c(2:5,1),2]
  logll <- cbind(logll, a)
}
colnames(logll) <- NULL

whole_log <- as.data.frame( t(t(apply(logll, mean, MARGIN = 1))) )
whole_log <-rbind(whole_log, NA)
rownames(whole_log) <- article_names

whole_tout <- cbind(whole_biais, whole_RMSE, whole_roc, whole_log)

xt <- xtable(whole_tout, digits = c(0,4,4,4,4,4,4,4,4,4,4,4,4,4))

# Affiche la version LaTeX
print(xt, type = "latex")

####### VALID

whole_biais_val <- data.frame()
for(i in method_names_val){
  a <- get(paste0("biais_",i))
  whole_biais_val <- rbind(whole_biais_val, a[1,])
}
rownames(whole_biais_val) <- article_names

whole_RMSE_val <- data.frame()
for(i in method_names_val){
  a <- get(paste0("RMSE_",i))
  whole_RMSE_val <- rbind(whole_RMSE_val, a[1,])
}
rownames(whole_RMSE_val) <- article_names

## roc
whole_roc_val <- data.frame()
for(i in 1:6){
  whole_roc_val <- rbind(whole_roc_val, unlist(unname(ROCmeanval_results[[1]][c(2:6,1)][i])))
}
rownames(whole_roc_val) <- article_names
colnames(whole_roc_val) <- colnames(whole_RMSE_val)


whole_tout_val <- cbind(whole_biais_val, whole_RMSE_val, whole_roc_val, whole_log)

xt_val <- xtable(whole_tout_val, digits = c(0,4,4,4,4,4,4,4,4,4,4,4,4,4))

# Affiche la version LaTeX
print(xt_val, type = "latex")

################################## N5000 
path1 <- "~/Documents/Rstudio/Simulations/Simulations mai 2025/Résultats/N5000_results/PH/1000ite_5000ind_PH_LC_25_04.Rdata"
load(path1)

##pour aller récupérer les logliklihood

path0 <- "~/Documents/Rstudio/Simulations/Simulations mai 2025/N5000/Output simulations_N5000_PH_LC/LOGLIK/"

method_names <- c("FLEX1.2", "FLEX1.4", "FLEX2.2", "FLEX2.4",
                  "PLANN", "PP")
method_names_val <- c("FLEX1.2val", "FLEX1.4val", "FLEX2.2val", "FLEX2.4val",
                      "PLANNval", "PPval")

article_names <- c("Flexnet1.2","Flexnet1.4","Flexnet2.2","Flexnet2.4", "PLANN", "Pohar Perme")

colnames(biais_PP) <- colnames(biais_FLEX1.2)
colnames(biais_PPval) <- colnames(biais_FLEX1.2)

colnames(RMSE_PP) <- colnames(RMSE_FLEX1.2)
colnames(RMSE_PPval) <- colnames(RMSE_FLEX1.2)

whole_biais <- data.frame()
for(i in method_names){
  a <- get(paste0("biais_",i))
  whole_biais <- rbind(whole_biais, a[1,])
}
rownames(whole_biais) <- article_names

whole_RMSE <- data.frame()
for(i in method_names){
  a <- get(paste0("RMSE_",i))
  whole_RMSE <- rbind(whole_RMSE, a[1,])
}
rownames(whole_RMSE) <- article_names

## roc
whole_roc <- data.frame()
for(i in 1:6){
  whole_roc <- rbind(whole_roc, unlist(unname(ROCmean_results[[1]][c(2:6,1)][i])))
}
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

logll <- unname(read.csv(paste0(path0,iterations[1],"_loglik.csv"), sep= " ", quote = "") )[c(2:5,1),2]
for(i in iterations[-1]){
  # assign(paste0("a_", i) , unname(read.csv(paste0(path0,i,"_loglik.csv")))[c(2:5,1),] )
  a <- unname(read.csv(paste0(path0,i,"_loglik.csv"), sep= " ", quote = "") )[c(2:5,1),2]
  logll <- cbind(logll, a)
}
colnames(logll) <- NULL

whole_log <- as.data.frame( t(t(apply(logll, mean, MARGIN = 1))) )
whole_log <-rbind(whole_log, NA)
rownames(whole_log) <- article_names

whole_tout <- cbind(whole_biais, whole_RMSE, whole_roc, whole_log)

xt <- xtable(whole_tout, digits = c(0,4,4,4,4,4,4,4,4,4,4,4,4,4))

# Affiche la version LaTeX
print(xt, type = "latex")

####### VALID

whole_biais_val <- data.frame()
for(i in method_names_val){
  a <- get(paste0("biais_",i))
  whole_biais_val <- rbind(whole_biais_val, a[1,])
}
rownames(whole_biais_val) <- article_names

whole_RMSE_val <- data.frame()
for(i in method_names_val){
  a <- get(paste0("RMSE_",i))
  whole_RMSE_val <- rbind(whole_RMSE_val, a[1,])
}
rownames(whole_RMSE_val) <- article_names

## roc
whole_roc_val <- data.frame()
for(i in 1:6){
  whole_roc_val <- rbind(whole_roc_val, unlist(unname(ROCmeanval_results[[1]][c(2:6,1)][i])))
}
rownames(whole_roc_val) <- article_names
colnames(whole_roc_val) <- colnames(whole_RMSE_val)


whole_tout_val <- cbind(whole_biais_val, whole_RMSE_val, whole_roc_val, whole_log)

xt_val <- xtable(whole_tout_val, digits = c(0,4,4,4,4,4,4,4,4,4,4,4,4,4))

# Affiche la version LaTeX
print(xt_val, type = "latex")


############################################################################################################
######                                    HOMME & COLON SUBSTRATA 
############################################################################################################

library(xtable)

##################################################  PH LC
################################## N1000 
path1 <- "~/Documents/Rstudio/Simulations/Simulations mai 2025/Résultats/N3000_results/PH/1000ite_3000ind_PH_LC_25_04.Rdata"
load(path1)

##pour aller récupérer les logliklihood

path0 <- "~/Documents/Rstudio/Simulations/Simulations mai 2025/N1000/Output simulations_N1000_PH_LC/LOGLIK/"

method_names <- c("FLEX1.2", "FLEX1.4", "FLEX2.2", "FLEX2.4",
                  "PLANN", "PP")
method_names_val <- c("FLEX1.2val", "FLEX1.4val", "FLEX2.2val", "FLEX2.4val",
                      "PLANNval", "PPval")

article_names <- c("Flexnet1.2","Flexnet1.4","Flexnet2.2","Flexnet2.4", "PLANN", "Pohar Perme")

colnames(biais_PP) <- colnames(biais_FLEX1.2)
colnames(biais_PPval) <- colnames(biais_FLEX1.2)

colnames(RMSE_PP) <- colnames(RMSE_FLEX1.2)
colnames(RMSE_PPval) <- colnames(RMSE_FLEX1.2)

whole_biais <- data.frame()
for(i in method_names){
  a <- get(paste0("biais_",i))
  whole_biais <- rbind(whole_biais, a[2,])
}
rownames(whole_biais) <- article_names

whole_RMSE <- data.frame()
for(i in method_names){
  a <- get(paste0("RMSE_",i))
  whole_RMSE <- rbind(whole_RMSE, a[2,])
}
rownames(whole_RMSE) <- article_names

## roc
whole_roc <- data.frame()
for(i in 1:6){
  whole_roc <- rbind(whole_roc, unlist(unname(ROCmean_results[[2]][c(2:6,1)][i])))
}
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

logll <- unname(read.csv(paste0(path0,iterations[1],"_loglik.csv"), sep= " ", quote = "") )[c(2:5,1),2]
for(i in iterations[-1]){
  # assign(paste0("a_", i) , unname(read.csv(paste0(path0,i,"_loglik.csv")))[c(2:5,1),] )
  a <- unname(read.csv(paste0(path0,i,"_loglik.csv"), sep= " ", quote = "") )[c(2:5,1),2]
  logll <- cbind(logll, a)
}
colnames(logll) <- NULL

whole_log <- as.data.frame( t(t(apply(logll, mean, MARGIN = 1))) )
whole_log <-rbind(whole_log, NA)
rownames(whole_log) <- article_names

whole_tout <- cbind(whole_biais, whole_RMSE, whole_roc, whole_log)

xt <- xtable(whole_tout, digits = c(0,4,4,4,4,4,4,4,4,4,4,4,4,4))

# Affiche la version LaTeX
print(xt, type = "latex")

####### VALID

whole_biais_val <- data.frame()
for(i in method_names_val){
  a <- get(paste0("biais_",i))
  whole_biais_val <- rbind(whole_biais_val, a[2,])
}
rownames(whole_biais_val) <- article_names

whole_RMSE_val <- data.frame()
for(i in method_names_val){
  a <- get(paste0("RMSE_",i))
  whole_RMSE_val <- rbind(whole_RMSE_val, a[2,])
}
rownames(whole_RMSE_val) <- article_names

## roc
whole_roc_val <- data.frame()
for(i in 1:6){
  whole_roc_val <- rbind(whole_roc_val, unlist(unname(ROCmeanval_results[[2]][c(2:6,1)][i])))
}
rownames(whole_roc_val) <- article_names
colnames(whole_roc_val) <- colnames(whole_RMSE_val)


whole_tout_val <- cbind(whole_biais_val, whole_RMSE_val, whole_roc_val, whole_log)

xt_val <- xtable(whole_tout_val, digits = c(0,4,4,4,4,4,4,4,4,4,4,4,4,4))

# Affiche la version LaTeX
print(xt_val, type = "latex")

################################## N3000 
path1 <- "~/Documents/Rstudio/Simulations/Simulations mai 2025/Résultats/N3000_results/NPH/1000ite_3000ind_NPH_LC_14_05.Rdata"
load(path1)

##pour aller récupérer les logliklihood

path0 <- "~/Documents/Rstudio/Simulations/Simulations mai 2025/N3000/Output simulations_N3000_NPH_LC/LOGLIK/"

method_names <- c("FLEX1.2", "FLEX1.4", "FLEX2.2", "FLEX2.4",
                  "PLANN", "PP")
method_names_val <- c("FLEX1.2val", "FLEX1.4val", "FLEX2.2val", "FLEX2.4val",
                      "PLANNval", "PPval")

article_names <- c("Flexnet1.2","Flexnet1.4","Flexnet2.2","Flexnet2.4", "PLANN", "Pohar Perme")

colnames(biais_PP) <- colnames(biais_FLEX1.2)
colnames(biais_PPval) <- colnames(biais_FLEX1.2)

colnames(RMSE_PP) <- colnames(RMSE_FLEX1.2)
colnames(RMSE_PPval) <- colnames(RMSE_FLEX1.2)

whole_biais <- data.frame()
for(i in method_names){
  a <- get(paste0("biais_",i))
  whole_biais <- rbind(whole_biais, a[2,])
}
rownames(whole_biais) <- article_names

whole_RMSE <- data.frame()
for(i in method_names){
  a <- get(paste0("RMSE_",i))
  whole_RMSE <- rbind(whole_RMSE, a[2,])
}
rownames(whole_RMSE) <- article_names

## roc
whole_roc <- data.frame()
for(i in 1:6){
  whole_roc <- rbind(whole_roc, unlist(unname(ROCmean_results[[2]][c(2:6,1)][i])))
}
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

logll <- unname(read.csv(paste0(path0,iterations[1],"_loglik.csv"), sep= " ", quote = "") )[c(2:5,1),2]
for(i in iterations[-1]){
  # assign(paste0("a_", i) , unname(read.csv(paste0(path0,i,"_loglik.csv")))[c(2:5,1),] )
  a <- unname(read.csv(paste0(path0,i,"_loglik.csv"), sep= " ", quote = "") )[c(2:5,1),2]
  logll <- cbind(logll, a)
}
colnames(logll) <- NULL

whole_log <- as.data.frame( t(t(apply(logll, mean, MARGIN = 1))) )
whole_log <-rbind(whole_log, NA)
rownames(whole_log) <- article_names

whole_tout <- cbind(whole_biais, whole_RMSE, whole_roc, whole_log)

xt <- xtable(whole_tout, digits = c(0,4,4,4,4,4,4,4,4,4,4,4,4,4))

# Affiche la version LaTeX
print(xt, type = "latex")

####### VALID

whole_biais_val <- data.frame()
for(i in method_names_val){
  a <- get(paste0("biais_",i))
  whole_biais_val <- rbind(whole_biais_val, a[2,])
}
rownames(whole_biais_val) <- article_names

whole_RMSE_val <- data.frame()
for(i in method_names_val){
  a <- get(paste0("RMSE_",i))
  whole_RMSE_val <- rbind(whole_RMSE_val, a[2,])
}
rownames(whole_RMSE_val) <- article_names

## roc
whole_roc_val <- data.frame()
for(i in 1:6){
  whole_roc_val <- rbind(whole_roc_val, unlist(unname(ROCmeanval_results[[2]][c(2:6,1)][i])))
}
rownames(whole_roc_val) <- article_names
colnames(whole_roc_val) <- colnames(whole_RMSE_val)


whole_tout_val <- cbind(whole_biais_val, whole_RMSE_val, whole_roc_val, whole_log)

xt_val <- xtable(whole_tout_val, digits = c(0,4,4,4,4,4,4,4,4,4,4,4,4,4))

# Affiche la version LaTeX
print(xt_val, type = "latex")

################################## N5000 
path1 <- "~/Documents/Rstudio/Simulations/Simulations mai 2025/Résultats/N5000_results/NPH/1000ite_5000ind_NPH_LC_14_05.Rdata"
load(path1)

##pour aller récupérer les logliklihood

path0 <- "~/Documents/Rstudio/Simulations/Simulations mai 2025/N5000/Output simulations_N5000_NPH_LC/LOGLIK/"

method_names <- c("FLEX1.2", "FLEX1.4", "FLEX2.2", "FLEX2.4",
                  "PLANN", "PP")
method_names_val <- c("FLEX1.2val", "FLEX1.4val", "FLEX2.2val", "FLEX2.4val",
                      "PLANNval", "PPval")

article_names <- c("Flexnet1.2","Flexnet1.4","Flexnet2.2","Flexnet2.4", "PLANN", "Pohar Perme")

colnames(biais_PP) <- colnames(biais_FLEX1.2)
colnames(biais_PPval) <- colnames(biais_FLEX1.2)

colnames(RMSE_PP) <- colnames(RMSE_FLEX1.2)
colnames(RMSE_PPval) <- colnames(RMSE_FLEX1.2)

whole_biais <- data.frame()
for(i in method_names){
  a <- get(paste0("biais_",i))
  whole_biais <- rbind(whole_biais, a[2,])
}
rownames(whole_biais) <- article_names

whole_RMSE <- data.frame()
for(i in method_names){
  a <- get(paste0("RMSE_",i))
  whole_RMSE <- rbind(whole_RMSE, a[2,])
}
rownames(whole_RMSE) <- article_names

## roc
whole_roc <- data.frame()
for(i in 1:6){
  whole_roc <- rbind(whole_roc, unlist(unname(ROCmean_results[[2]][c(2:6,1)][i])))
}
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

logll <- unname(read.csv(paste0(path0,iterations[1],"_loglik.csv"), sep= " ", quote = "") )[c(2:5,1),2]
for(i in iterations[-1]){
  # assign(paste0("a_", i) , unname(read.csv(paste0(path0,i,"_loglik.csv")))[c(2:5,1),] )
  a <- unname(read.csv(paste0(path0,i,"_loglik.csv"), sep= " ", quote = "") )[c(2:5,1),2]
  logll <- cbind(logll, a)
}
colnames(logll) <- NULL

whole_log <- as.data.frame( t(t(apply(logll, mean, MARGIN = 1))) )
whole_log <-rbind(whole_log, NA)
rownames(whole_log) <- article_names

whole_tout <- cbind(whole_biais, whole_RMSE, whole_roc, whole_log)

xt <- xtable(whole_tout, digits = c(0,4,4,4,4,4,4,4,4,4,4,4,4,4))

# Affiche la version LaTeX
print(xt, type = "latex")

####### VALID

whole_biais_val <- data.frame()
for(i in method_names_val){
  a <- get(paste0("biais_",i))
  whole_biais_val <- rbind(whole_biais_val, a[2,])
}
rownames(whole_biais_val) <- article_names

whole_RMSE_val <- data.frame()
for(i in method_names_val){
  a <- get(paste0("RMSE_",i))
  whole_RMSE_val <- rbind(whole_RMSE_val, a[2,])
}
rownames(whole_RMSE_val) <- article_names

## roc
whole_roc_val <- data.frame()
for(i in 1:6){
  whole_roc_val <- rbind(whole_roc_val, unlist(unname(ROCmeanval_results[[2]][c(2:6,1)][i])))
}
rownames(whole_roc_val) <- article_names
colnames(whole_roc_val) <- colnames(whole_RMSE_val)


whole_tout_val <- cbind(whole_biais_val, whole_RMSE_val, whole_roc_val, whole_log)

xt_val <- xtable(whole_tout_val, digits = c(0,4,4,4,4,4,4,4,4,4,4,4,4,4))

# Affiche la version LaTeX
print(xt_val, type = "latex")

############################################################################################################
######                                    HOMME & RECTUM SUBSTRATA 
############################################################################################################

library(xtable)

##################################################  PH LC
################################## N1000 
path1 <- "~/Documents/Rstudio/Simulations/Simulations mai 2025/Résultats/N3000_results/PH/1000ite_3000ind_PH_LC_25_04.Rdata"
load(path1)

##pour aller récupérer les logliklihood

path0 <- "~/Documents/Rstudio/Simulations/Simulations mai 2025/N1000/Output simulations_N1000_PH_LC/LOGLIK/"

method_names <- c("FLEX1.2", "FLEX1.4", "FLEX2.2", "FLEX2.4",
                  "PLANN", "PP")
method_names_val <- c("FLEX1.2val", "FLEX1.4val", "FLEX2.2val", "FLEX2.4val",
                      "PLANNval", "PPval")

article_names <- c("Flexnet1.2","Flexnet1.4","Flexnet2.2","Flexnet2.4", "PLANN", "Pohar Perme")

colnames(biais_PP) <- colnames(biais_FLEX1.2)
colnames(biais_PPval) <- colnames(biais_FLEX1.2)

colnames(RMSE_PP) <- colnames(RMSE_FLEX1.2)
colnames(RMSE_PPval) <- colnames(RMSE_FLEX1.2)

whole_biais <- data.frame()
for(i in method_names){
  a <- get(paste0("biais_",i))
  whole_biais <- rbind(whole_biais, a[3,])
}
rownames(whole_biais) <- article_names

whole_RMSE <- data.frame()
for(i in method_names){
  a <- get(paste0("RMSE_",i))
  whole_RMSE <- rbind(whole_RMSE, a[3,])
}
rownames(whole_RMSE) <- article_names

## roc
whole_roc <- data.frame()
for(i in 1:6){
  whole_roc <- rbind(whole_roc, unlist(unname(ROCmean_results[[3]][c(2:6,1)][i])))
}
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

logll <- unname(read.csv(paste0(path0,iterations[1],"_loglik.csv"), sep= " ", quote = "") )[c(2:5,1),2]
for(i in iterations[-1]){
  # assign(paste0("a_", i) , unname(read.csv(paste0(path0,i,"_loglik.csv")))[c(2:5,1),] )
  a <- unname(read.csv(paste0(path0,i,"_loglik.csv"), sep= " ", quote = "") )[c(2:5,1),2]
  logll <- cbind(logll, a)
}
colnames(logll) <- NULL

whole_log <- as.data.frame( t(t(apply(logll, mean, MARGIN = 1))) )
whole_log <-rbind(whole_log, NA)
rownames(whole_log) <- article_names

whole_tout <- cbind(whole_biais, whole_RMSE, whole_roc, whole_log)

xt <- xtable(whole_tout, digits = c(0,4,4,4,4,4,4,4,4,4,4,4,4,4))

# Affiche la version LaTeX
print(xt, type = "latex")

####### VALID

whole_biais_val <- data.frame()
for(i in method_names_val){
  a <- get(paste0("biais_",i))
  whole_biais_val <- rbind(whole_biais_val, a[3,])
}
rownames(whole_biais_val) <- article_names

whole_RMSE_val <- data.frame()
for(i in method_names_val){
  a <- get(paste0("RMSE_",i))
  whole_RMSE_val <- rbind(whole_RMSE_val, a[3,])
}
rownames(whole_RMSE_val) <- article_names

## roc
whole_roc_val <- data.frame()
for(i in 1:6){
  whole_roc_val <- rbind(whole_roc_val, unlist(unname(ROCmeanval_results[[3]][c(2:6,1)][i])))
}
rownames(whole_roc_val) <- article_names
colnames(whole_roc_val) <- colnames(whole_RMSE_val)

whole_log = as.data.frame(t(t(c(rep(NA,6)))))

whole_tout_val <- cbind(whole_biais_val, whole_RMSE_val, whole_roc_val, whole_log)

xt_val <- xtable(whole_tout_val, digits = c(0,4,4,4,4,4,4,4,4,4,4,4,4,4))

# Affiche la version LaTeX
print(xt_val, type = "latex")

################################## N3000 
path1 <- "~/Documents/Rstudio/Simulations/Simulations mai 2025/Résultats/N3000_results/NPH/1000ite_3000ind_NPH_LC_14_05.Rdata"
load(path1)

##pour aller récupérer les logliklihood

path0 <- "~/Documents/Rstudio/Simulations/Simulations mai 2025/N3000/Output simulations_N3000_NPH_LC/LOGLIK/"

method_names <- c("FLEX1.2", "FLEX1.4", "FLEX2.2", "FLEX2.4",
                  "PLANN", "PP")
method_names_val <- c("FLEX1.2val", "FLEX1.4val", "FLEX2.2val", "FLEX2.4val",
                      "PLANNval", "PPval")

article_names <- c("Flexnet1.2","Flexnet1.4","Flexnet2.2","Flexnet2.4", "PLANN", "Pohar Perme")

colnames(biais_PP) <- colnames(biais_FLEX1.2)
colnames(biais_PPval) <- colnames(biais_FLEX1.2)

colnames(RMSE_PP) <- colnames(RMSE_FLEX1.2)
colnames(RMSE_PPval) <- colnames(RMSE_FLEX1.2)

whole_biais <- data.frame()
for(i in method_names){
  a <- get(paste0("biais_",i))
  whole_biais <- rbind(whole_biais, a[3,])
}
rownames(whole_biais) <- article_names

whole_RMSE <- data.frame()
for(i in method_names){
  a <- get(paste0("RMSE_",i))
  whole_RMSE <- rbind(whole_RMSE, a[3,])
}
rownames(whole_RMSE) <- article_names

## roc
whole_roc <- data.frame()
for(i in 1:6){
  whole_roc <- rbind(whole_roc, unlist(unname(ROCmean_results[[3]][c(2:6,1)][i])))
}
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

logll <- unname(read.csv(paste0(path0,iterations[1],"_loglik.csv"), sep= " ", quote = "") )[c(2:5,1),2]
for(i in iterations[-1]){
  # assign(paste0("a_", i) , unname(read.csv(paste0(path0,i,"_loglik.csv")))[c(2:5,1),] )
  a <- unname(read.csv(paste0(path0,i,"_loglik.csv"), sep= " ", quote = "") )[c(2:5,1),2]
  logll <- cbind(logll, a)
}
colnames(logll) <- NULL

whole_log <- as.data.frame( t(t(apply(logll, mean, MARGIN = 1))) )
whole_log <-rbind(whole_log, NA)
rownames(whole_log) <- article_names

whole_tout <- cbind(whole_biais, whole_RMSE, whole_roc, whole_log)

xt <- xtable(whole_tout, digits = c(0,4,4,4,4,4,4,4,4,4,4,4,4,4))

# Affiche la version LaTeX
print(xt, type = "latex")

####### VALID

whole_biais_val <- data.frame()
for(i in method_names_val){
  a <- get(paste0("biais_",i))
  whole_biais_val <- rbind(whole_biais_val, a[3,])
}
rownames(whole_biais_val) <- article_names

whole_RMSE_val <- data.frame()
for(i in method_names_val){
  a <- get(paste0("RMSE_",i))
  whole_RMSE_val <- rbind(whole_RMSE_val, a[3,])
}
rownames(whole_RMSE_val) <- article_names

## roc
whole_roc_val <- data.frame()
for(i in 1:6){
  whole_roc_val <- rbind(whole_roc_val, unlist(unname(ROCmeanval_results[[3]][c(2:6,1)][i])))
}
rownames(whole_roc_val) <- article_names
colnames(whole_roc_val) <- colnames(whole_RMSE_val)

whole_log = as.data.frame(t(t(c(rep(NA,6)))))

whole_tout_val <- cbind(whole_biais_val, whole_RMSE_val, whole_roc_val, whole_log)

xt_val <- xtable(whole_tout_val, digits = c(0,4,4,4,4,4,4,4,4,4,4,4,4,4))

# Affiche la version LaTeX
print(xt_val, type = "latex")

################################## N5000 
path1 <- "~/Documents/Rstudio/Simulations/Simulations mai 2025/Résultats/N5000_results/NPH/1000ite_5000ind_NPH_LC_14_05.Rdata"
load(path1)

##pour aller récupérer les logliklihood

path0 <- "~/Documents/Rstudio/Simulations/Simulations mai 2025/N5000/Output simulations_N5000_NPH_LC/LOGLIK/"

method_names <- c("FLEX1.2", "FLEX1.4", "FLEX2.2", "FLEX2.4",
                  "PLANN", "PP")
method_names_val <- c("FLEX1.2val", "FLEX1.4val", "FLEX2.2val", "FLEX2.4val",
                      "PLANNval", "PPval")

article_names <- c("Flexnet1.2","Flexnet1.4","Flexnet2.2","Flexnet2.4", "PLANN", "Pohar Perme")

colnames(biais_PP) <- colnames(biais_FLEX1.2)
colnames(biais_PPval) <- colnames(biais_FLEX1.2)

colnames(RMSE_PP) <- colnames(RMSE_FLEX1.2)
colnames(RMSE_PPval) <- colnames(RMSE_FLEX1.2)

whole_biais <- data.frame()
for(i in method_names){
  a <- get(paste0("biais_",i))
  whole_biais <- rbind(whole_biais, a[3,])
}
rownames(whole_biais) <- article_names

whole_RMSE <- data.frame()
for(i in method_names){
  a <- get(paste0("RMSE_",i))
  whole_RMSE <- rbind(whole_RMSE, a[3,])
}
rownames(whole_RMSE) <- article_names

## roc
whole_roc <- data.frame()
for(i in 1:6){
  whole_roc <- rbind(whole_roc, unlist(unname(ROCmean_results[[3]][c(2:6,1)][i])))
}
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

logll <- unname(read.csv(paste0(path0,iterations[1],"_loglik.csv"), sep= " ", quote = "") )[c(2:5,1),2]
for(i in iterations[-1]){
  # assign(paste0("a_", i) , unname(read.csv(paste0(path0,i,"_loglik.csv")))[c(2:5,1),] )
  a <- unname(read.csv(paste0(path0,i,"_loglik.csv"), sep= " ", quote = "") )[c(2:5,1),2]
  logll <- cbind(logll, a)
}
colnames(logll) <- NULL

whole_log <- as.data.frame( t(t(apply(logll, mean, MARGIN = 1))) )
whole_log <-rbind(whole_log, NA)
rownames(whole_log) <- article_names

whole_tout <- cbind(whole_biais, whole_RMSE, whole_roc, whole_log)

xt <- xtable(whole_tout, digits = c(0,4,4,4,4,4,4,4,4,4,4,4,4,4))

# Affiche la version LaTeX
print(xt, type = "latex")

####### VALID

whole_biais_val <- data.frame()
for(i in method_names_val){
  a <- get(paste0("biais_",i))
  whole_biais_val <- rbind(whole_biais_val, a[3,])
}
rownames(whole_biais_val) <- article_names

whole_RMSE_val <- data.frame()
for(i in method_names_val){
  a <- get(paste0("RMSE_",i))
  whole_RMSE_val <- rbind(whole_RMSE_val, a[3,])
}
rownames(whole_RMSE_val) <- article_names

## roc
whole_roc_val <- data.frame()
for(i in 1:6){
  whole_roc_val <- rbind(whole_roc_val, unlist(unname(ROCmeanval_results[[3]][c(2:6,1)][i])))
}
rownames(whole_roc_val) <- article_names
colnames(whole_roc_val) <- colnames(whole_RMSE_val)

whole_log = as.data.frame(t(t(c(rep(NA,6)))))

whole_tout_val <- cbind(whole_biais_val, whole_RMSE_val, whole_roc_val, whole_log)

xt_val <- xtable(whole_tout_val, digits = c(0,4,4,4,4,4,4,4,4,4,4,4,4,4))

# Affiche la version LaTeX
print(xt_val, type = "latex")



############################################################################################################
######                                    FEMME & COLON SUBSTRATA 
############################################################################################################

library(xtable)

##################################################  PH LC
################################## N1000 
path1 <- "~/Documents/Rstudio/Simulations/Simulations mai 2025/Résultats/N3000_results/PH/1000ite_3000ind_PH_LC_25_04.Rdata"
load(path1)

##pour aller récupérer les logliklihood

path0 <- "~/Documents/Rstudio/Simulations/Simulations mai 2025/N1000/Output simulations_N1000_PH_LC/LOGLIK/"

method_names <- c("FLEX1.2", "FLEX1.4", "FLEX2.2", "FLEX2.4",
                  "PLANN", "PP")
method_names_val <- c("FLEX1.2val", "FLEX1.4val", "FLEX2.2val", "FLEX2.4val",
                      "PLANNval", "PPval")

article_names <- c("Flexnet1.2","Flexnet1.4","Flexnet2.2","Flexnet2.4", "PLANN", "Pohar Perme")

colnames(biais_PP) <- colnames(biais_FLEX1.2)
colnames(biais_PPval) <- colnames(biais_FLEX1.2)

colnames(RMSE_PP) <- colnames(RMSE_FLEX1.2)
colnames(RMSE_PPval) <- colnames(RMSE_FLEX1.2)

whole_biais <- data.frame()
for(i in method_names){
  a <- get(paste0("biais_",i))
  whole_biais <- rbind(whole_biais, a[4,])
}
rownames(whole_biais) <- article_names

whole_RMSE <- data.frame()
for(i in method_names){
  a <- get(paste0("RMSE_",i))
  whole_RMSE <- rbind(whole_RMSE, a[4,])
}
rownames(whole_RMSE) <- article_names

## roc
whole_roc <- data.frame()
for(i in 1:6){
  whole_roc <- rbind(whole_roc, unlist(unname(ROCmean_results[[4]][c(2:6,1)][i])))
}
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

logll <- unname(read.csv(paste0(path0,iterations[1],"_loglik.csv"), sep= " ", quote = "") )[c(2:5,1),2]
for(i in iterations[-1]){
  # assign(paste0("a_", i) , unname(read.csv(paste0(path0,i,"_loglik.csv")))[c(2:5,1),] )
  a <- unname(read.csv(paste0(path0,i,"_loglik.csv"), sep= " ", quote = "") )[c(2:5,1),2]
  logll <- cbind(logll, a)
}
colnames(logll) <- NULL

whole_log <- as.data.frame( t(t(apply(logll, mean, MARGIN = 1))) )
whole_log <-rbind(whole_log, NA)
rownames(whole_log) <- article_names

whole_tout <- cbind(whole_biais, whole_RMSE, whole_roc, whole_log)

xt <- xtable(whole_tout, digits = c(0,4,4,4,4,4,4,4,4,4,4,4,4,4))

# Affiche la version LaTeX
print(xt, type = "latex")

####### VALID

whole_biais_val <- data.frame()
for(i in method_names_val){
  a <- get(paste0("biais_",i))
  whole_biais_val <- rbind(whole_biais_val, a[4,])
}
rownames(whole_biais_val) <- article_names

whole_RMSE_val <- data.frame()
for(i in method_names_val){
  a <- get(paste0("RMSE_",i))
  whole_RMSE_val <- rbind(whole_RMSE_val, a[4,])
}
rownames(whole_RMSE_val) <- article_names

## roc
whole_roc_val <- data.frame()
for(i in 1:6){
  whole_roc_val <- rbind(whole_roc_val, unlist(unname(ROCmeanval_results[[4]][c(2:6,1)][i])))
}
rownames(whole_roc_val) <- article_names
colnames(whole_roc_val) <- colnames(whole_RMSE_val)

whole_log = as.data.frame(t(t(c(rep(NA,6)))))

whole_tout_val <- cbind(whole_biais_val, whole_RMSE_val, whole_roc_val, whole_log)

xt_val <- xtable(whole_tout_val, digits = c(0,4,4,4,4,4,4,4,4,4,4,4,4,4))

# Affiche la version LaTeX
print(xt_val, type = "latex")

################################## N3000 
path1 <- "~/Documents/Rstudio/Simulations/Simulations mai 2025/Résultats/N3000_results/NPH/1000ite_3000ind_NPH_LC_14_05.Rdata"
load(path1)

##pour aller récupérer les logliklihood

path0 <- "~/Documents/Rstudio/Simulations/Simulations mai 2025/N3000/Output simulations_N3000_NPH_LC/LOGLIK/"

method_names <- c("FLEX1.2", "FLEX1.4", "FLEX2.2", "FLEX2.4",
                  "PLANN", "PP")
method_names_val <- c("FLEX1.2val", "FLEX1.4val", "FLEX2.2val", "FLEX2.4val",
                      "PLANNval", "PPval")

article_names <- c("Flexnet1.2","Flexnet1.4","Flexnet2.2","Flexnet2.4", "PLANN", "Pohar Perme")

colnames(biais_PP) <- colnames(biais_FLEX1.2)
colnames(biais_PPval) <- colnames(biais_FLEX1.2)

colnames(RMSE_PP) <- colnames(RMSE_FLEX1.2)
colnames(RMSE_PPval) <- colnames(RMSE_FLEX1.2)

whole_biais <- data.frame()
for(i in method_names){
  a <- get(paste0("biais_",i))
  whole_biais <- rbind(whole_biais, a[4,])
}
rownames(whole_biais) <- article_names

whole_RMSE <- data.frame()
for(i in method_names){
  a <- get(paste0("RMSE_",i))
  whole_RMSE <- rbind(whole_RMSE, a[4,])
}
rownames(whole_RMSE) <- article_names

## roc
whole_roc <- data.frame()
for(i in 1:6){
  whole_roc <- rbind(whole_roc, unlist(unname(ROCmean_results[[4]][c(2:6,1)][i])))
}
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

logll <- unname(read.csv(paste0(path0,iterations[1],"_loglik.csv"), sep= " ", quote = "") )[c(2:5,1),2]
for(i in iterations[-1]){
  # assign(paste0("a_", i) , unname(read.csv(paste0(path0,i,"_loglik.csv")))[c(2:5,1),] )
  a <- unname(read.csv(paste0(path0,i,"_loglik.csv"), sep= " ", quote = "") )[c(2:5,1),2]
  logll <- cbind(logll, a)
}
colnames(logll) <- NULL

whole_log <- as.data.frame( t(t(apply(logll, mean, MARGIN = 1))) )
whole_log <-rbind(whole_log, NA)
rownames(whole_log) <- article_names

whole_tout <- cbind(whole_biais, whole_RMSE, whole_roc, whole_log)

xt <- xtable(whole_tout, digits = c(0,4,4,4,4,4,4,4,4,4,4,4,4,4))

# Affiche la version LaTeX
print(xt, type = "latex")

####### VALID

whole_biais_val <- data.frame()
for(i in method_names_val){
  a <- get(paste0("biais_",i))
  whole_biais_val <- rbind(whole_biais_val, a[4,])
}
rownames(whole_biais_val) <- article_names

whole_RMSE_val <- data.frame()
for(i in method_names_val){
  a <- get(paste0("RMSE_",i))
  whole_RMSE_val <- rbind(whole_RMSE_val, a[4,])
}
rownames(whole_RMSE_val) <- article_names

## roc
whole_roc_val <- data.frame()
for(i in 1:6){
  whole_roc_val <- rbind(whole_roc_val, unlist(unname(ROCmeanval_results[[4]][c(2:6,1)][i])))
}
rownames(whole_roc_val) <- article_names
colnames(whole_roc_val) <- colnames(whole_RMSE_val)

whole_log = as.data.frame(t(t(c(rep(NA,6)))))

whole_tout_val <- cbind(whole_biais_val, whole_RMSE_val, whole_roc_val, whole_log)

xt_val <- xtable(whole_tout_val, digits = c(0,4,4,4,4,4,4,4,4,4,4,4,4,4))

# Affiche la version LaTeX
print(xt_val, type = "latex")

################################## N5000 
path1 <- "~/Documents/Rstudio/Simulations/Simulations mai 2025/Résultats/N5000_results/NPH/1000ite_5000ind_NPH_LC_14_05.Rdata"
load(path1)

##pour aller récupérer les logliklihood

path0 <- "~/Documents/Rstudio/Simulations/Simulations mai 2025/N5000/Output simulations_N5000_NPH_LC/LOGLIK/"

method_names <- c("FLEX1.2", "FLEX1.4", "FLEX2.2", "FLEX2.4",
                  "PLANN", "PP")
method_names_val <- c("FLEX1.2val", "FLEX1.4val", "FLEX2.2val", "FLEX2.4val",
                      "PLANNval", "PPval")

article_names <- c("Flexnet1.2","Flexnet1.4","Flexnet2.2","Flexnet2.4", "PLANN", "Pohar Perme")

colnames(biais_PP) <- colnames(biais_FLEX1.2)
colnames(biais_PPval) <- colnames(biais_FLEX1.2)

colnames(RMSE_PP) <- colnames(RMSE_FLEX1.2)
colnames(RMSE_PPval) <- colnames(RMSE_FLEX1.2)

whole_biais <- data.frame()
for(i in method_names){
  a <- get(paste0("biais_",i))
  whole_biais <- rbind(whole_biais, a[4,])
}
rownames(whole_biais) <- article_names

whole_RMSE <- data.frame()
for(i in method_names){
  a <- get(paste0("RMSE_",i))
  whole_RMSE <- rbind(whole_RMSE, a[4,])
}
rownames(whole_RMSE) <- article_names

## roc
whole_roc <- data.frame()
for(i in 1:6){
  whole_roc <- rbind(whole_roc, unlist(unname(ROCmean_results[[4]][c(2:6,1)][i])))
}
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

logll <- unname(read.csv(paste0(path0,iterations[1],"_loglik.csv"), sep= " ", quote = "") )[c(2:5,1),2]
for(i in iterations[-1]){
  # assign(paste0("a_", i) , unname(read.csv(paste0(path0,i,"_loglik.csv")))[c(2:5,1),] )
  a <- unname(read.csv(paste0(path0,i,"_loglik.csv"), sep= " ", quote = "") )[c(2:5,1),2]
  logll <- cbind(logll, a)
}
colnames(logll) <- NULL

whole_log <- as.data.frame( t(t(apply(logll, mean, MARGIN = 1))) )
whole_log <-rbind(whole_log, NA)
rownames(whole_log) <- article_names

whole_tout <- cbind(whole_biais, whole_RMSE, whole_roc, whole_log)

xt <- xtable(whole_tout, digits = c(0,4,4,4,4,4,4,4,4,4,4,4,4,4))

# Affiche la version LaTeX
print(xt, type = "latex")

####### VALID

whole_biais_val <- data.frame()
for(i in method_names_val){
  a <- get(paste0("biais_",i))
  whole_biais_val <- rbind(whole_biais_val, a[4,])
}
rownames(whole_biais_val) <- article_names

whole_RMSE_val <- data.frame()
for(i in method_names_val){
  a <- get(paste0("RMSE_",i))
  whole_RMSE_val <- rbind(whole_RMSE_val, a[4,])
}
rownames(whole_RMSE_val) <- article_names

## roc
whole_roc_val <- data.frame()
for(i in 1:6){
  whole_roc_val <- rbind(whole_roc_val, unlist(unname(ROCmeanval_results[[4]][c(2:6,1)][i])))
}
rownames(whole_roc_val) <- article_names
colnames(whole_roc_val) <- colnames(whole_RMSE_val)

whole_log = as.data.frame(t(t(c(rep(NA,6)))))

whole_tout_val <- cbind(whole_biais_val, whole_RMSE_val, whole_roc_val, whole_log)

xt_val <- xtable(whole_tout_val, digits = c(0,4,4,4,4,4,4,4,4,4,4,4,4,4))

# Affiche la version LaTeX
print(xt_val, type = "latex")


############################################################################################################
######                                    FEMME & RECTUM SUBSTRATA 
############################################################################################################

library(xtable)

##################################################  PH LC
################################## N1000 
path1 <- "~/Documents/Rstudio/Simulations/Simulations mai 2025/Résultats/N3000_results/PH/1000ite_3000ind_PH_LC_25_04.Rdata"
load(path1)

##pour aller récupérer les logliklihood

path0 <- "~/Documents/Rstudio/Simulations/Simulations mai 2025/N1000/Output simulations_N1000_PH_LC/LOGLIK/"

method_names <- c("FLEX1.2", "FLEX1.4", "FLEX2.2", "FLEX2.4",
                  "PLANN", "PP")
method_names_val <- c("FLEX1.2val", "FLEX1.4val", "FLEX2.2val", "FLEX2.4val",
                      "PLANNval", "PPval")

article_names <- c("Flexnet1.2","Flexnet1.4","Flexnet2.2","Flexnet2.4", "PLANN", "Pohar Perme")

colnames(biais_PP) <- colnames(biais_FLEX1.2)
colnames(biais_PPval) <- colnames(biais_FLEX1.2)

colnames(RMSE_PP) <- colnames(RMSE_FLEX1.2)
colnames(RMSE_PPval) <- colnames(RMSE_FLEX1.2)

whole_biais <- data.frame()
for(i in method_names){
  a <- get(paste0("biais_",i))
  whole_biais <- rbind(whole_biais, a[5,])
}
rownames(whole_biais) <- article_names

whole_RMSE <- data.frame()
for(i in method_names){
  a <- get(paste0("RMSE_",i))
  whole_RMSE <- rbind(whole_RMSE, a[5,])
}
rownames(whole_RMSE) <- article_names

## roc
whole_roc <- data.frame()
for(i in 1:6){
  whole_roc <- rbind(whole_roc, unlist(unname(ROCmean_results[[5]][c(2:6,1)][i])))
}
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

logll <- unname(read.csv(paste0(path0,iterations[1],"_loglik.csv"), sep= " ", quote = "") )[c(2:5,1),2]
for(i in iterations[-1]){
  # assign(paste0("a_", i) , unname(read.csv(paste0(path0,i,"_loglik.csv")))[c(2:5,1),] )
  a <- unname(read.csv(paste0(path0,i,"_loglik.csv"), sep= " ", quote = "") )[c(2:5,1),2]
  logll <- cbind(logll, a)
}
colnames(logll) <- NULL

whole_log <- as.data.frame( t(t(apply(logll, mean, MARGIN = 1))) )
whole_log <-rbind(whole_log, NA)
rownames(whole_log) <- article_names

whole_tout <- cbind(whole_biais, whole_RMSE, whole_roc, whole_log)

xt <- xtable(whole_tout, digits = c(0,4,4,4,4,4,4,4,4,4,4,4,4,4))

# Affiche la version LaTeX
print(xt, type = "latex")

####### VALID

whole_biais_val <- data.frame()
for(i in method_names_val){
  a <- get(paste0("biais_",i))
  whole_biais_val <- rbind(whole_biais_val, a[5,])
}
rownames(whole_biais_val) <- article_names

whole_RMSE_val <- data.frame()
for(i in method_names_val){
  a <- get(paste0("RMSE_",i))
  whole_RMSE_val <- rbind(whole_RMSE_val, a[5,])
}
rownames(whole_RMSE_val) <- article_names

## roc
whole_roc_val <- data.frame()
for(i in 1:6){
  whole_roc_val <- rbind(whole_roc_val, unlist(unname(ROCmeanval_results[[5]][c(2:6,1)][i])))
}
rownames(whole_roc_val) <- article_names
colnames(whole_roc_val) <- colnames(whole_RMSE_val)

whole_log = as.data.frame(t(t(c(rep(NA,6)))))

whole_tout_val <- cbind(whole_biais_val, whole_RMSE_val, whole_roc_val, whole_log)

xt_val <- xtable(whole_tout_val, digits = c(0,4,4,4,4,4,4,4,4,4,4,4,4,4))

# Affiche la version LaTeX
print(xt_val, type = "latex")

################################## N3000 
path1 <- "~/Documents/Rstudio/Simulations/Simulations mai 2025/Résultats/N3000_results/NPH/1000ite_3000ind_NPH_LC_14_05.Rdata"
load(path1)

##pour aller récupérer les logliklihood

path0 <- "~/Documents/Rstudio/Simulations/Simulations mai 2025/N3000/Output simulations_N3000_NPH_LC/LOGLIK/"

method_names <- c("FLEX1.2", "FLEX1.4", "FLEX2.2", "FLEX2.4",
                  "PLANN", "PP")
method_names_val <- c("FLEX1.2val", "FLEX1.4val", "FLEX2.2val", "FLEX2.4val",
                      "PLANNval", "PPval")

article_names <- c("Flexnet1.2","Flexnet1.4","Flexnet2.2","Flexnet2.4", "PLANN", "Pohar Perme")

colnames(biais_PP) <- colnames(biais_FLEX1.2)
colnames(biais_PPval) <- colnames(biais_FLEX1.2)

colnames(RMSE_PP) <- colnames(RMSE_FLEX1.2)
colnames(RMSE_PPval) <- colnames(RMSE_FLEX1.2)

whole_biais <- data.frame()
for(i in method_names){
  a <- get(paste0("biais_",i))
  whole_biais <- rbind(whole_biais, a[5,])
}
rownames(whole_biais) <- article_names

whole_RMSE <- data.frame()
for(i in method_names){
  a <- get(paste0("RMSE_",i))
  whole_RMSE <- rbind(whole_RMSE, a[5,])
}
rownames(whole_RMSE) <- article_names

## roc
whole_roc <- data.frame()
for(i in 1:6){
  whole_roc <- rbind(whole_roc, unlist(unname(ROCmean_results[[5]][c(2:6,1)][i])))
}
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

logll <- unname(read.csv(paste0(path0,iterations[1],"_loglik.csv"), sep= " ", quote = "") )[c(2:5,1),2]
for(i in iterations[-1]){
  # assign(paste0("a_", i) , unname(read.csv(paste0(path0,i,"_loglik.csv")))[c(2:5,1),] )
  a <- unname(read.csv(paste0(path0,i,"_loglik.csv"), sep= " ", quote = "") )[c(2:5,1),2]
  logll <- cbind(logll, a)
}
colnames(logll) <- NULL

whole_log <- as.data.frame( t(t(apply(logll, mean, MARGIN = 1))) )
whole_log <-rbind(whole_log, NA)
rownames(whole_log) <- article_names

whole_tout <- cbind(whole_biais, whole_RMSE, whole_roc, whole_log)

xt <- xtable(whole_tout, digits = c(0,4,4,4,4,4,4,4,4,4,4,4,4,4))

# Affiche la version LaTeX
print(xt, type = "latex")

####### VALID

whole_biais_val <- data.frame()
for(i in method_names_val){
  a <- get(paste0("biais_",i))
  whole_biais_val <- rbind(whole_biais_val, a[5,])
}
rownames(whole_biais_val) <- article_names

whole_RMSE_val <- data.frame()
for(i in method_names_val){
  a <- get(paste0("RMSE_",i))
  whole_RMSE_val <- rbind(whole_RMSE_val, a[5,])
}
rownames(whole_RMSE_val) <- article_names

## roc
whole_roc_val <- data.frame()
for(i in 1:6){
  whole_roc_val <- rbind(whole_roc_val, unlist(unname(ROCmeanval_results[[5]][c(2:6,1)][i])))
}
rownames(whole_roc_val) <- article_names
colnames(whole_roc_val) <- colnames(whole_RMSE_val)

whole_log = as.data.frame(t(t(c(rep(NA,6)))))

whole_tout_val <- cbind(whole_biais_val, whole_RMSE_val, whole_roc_val, whole_log)

xt_val <- xtable(whole_tout_val, digits = c(0,4,4,4,4,4,4,4,4,4,4,4,4,4))

# Affiche la version LaTeX
print(xt_val, type = "latex")

################################## N5000 
path1 <- "~/Documents/Rstudio/Simulations/Simulations mai 2025/Résultats/N5000_results/NPH/1000ite_5000ind_NPH_LC_14_05.Rdata"
load(path1)

##pour aller récupérer les logliklihood

path0 <- "~/Documents/Rstudio/Simulations/Simulations mai 2025/N5000/Output simulations_N5000_NPH_LC/LOGLIK/"

method_names <- c("FLEX1.2", "FLEX1.4", "FLEX2.2", "FLEX2.4",
                  "PLANN", "PP")
method_names_val <- c("FLEX1.2val", "FLEX1.4val", "FLEX2.2val", "FLEX2.4val",
                      "PLANNval", "PPval")

article_names <- c("Flexnet1.2","Flexnet1.4","Flexnet2.2","Flexnet2.4", "PLANN", "Pohar Perme")

colnames(biais_PP) <- colnames(biais_FLEX1.2)
colnames(biais_PPval) <- colnames(biais_FLEX1.2)

colnames(RMSE_PP) <- colnames(RMSE_FLEX1.2)
colnames(RMSE_PPval) <- colnames(RMSE_FLEX1.2)

whole_biais <- data.frame()
for(i in method_names){
  a <- get(paste0("biais_",i))
  whole_biais <- rbind(whole_biais, a[5,])
}
rownames(whole_biais) <- article_names

whole_RMSE <- data.frame()
for(i in method_names){
  a <- get(paste0("RMSE_",i))
  whole_RMSE <- rbind(whole_RMSE, a[5,])
}
rownames(whole_RMSE) <- article_names

## roc
whole_roc <- data.frame()
for(i in 1:6){
  whole_roc <- rbind(whole_roc, unlist(unname(ROCmean_results[[5]][c(2:6,1)][i])))
}
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

logll <- unname(read.csv(paste0(path0,iterations[1],"_loglik.csv"), sep= " ", quote = "") )[c(2:5,1),2]
for(i in iterations[-1]){
  # assign(paste0("a_", i) , unname(read.csv(paste0(path0,i,"_loglik.csv")))[c(2:5,1),] )
  a <- unname(read.csv(paste0(path0,i,"_loglik.csv"), sep= " ", quote = "") )[c(2:5,1),2]
  logll <- cbind(logll, a)
}
colnames(logll) <- NULL

whole_log <- as.data.frame( t(t(apply(logll, mean, MARGIN = 1))) )
whole_log <-rbind(whole_log, NA)
rownames(whole_log) <- article_names

whole_tout <- cbind(whole_biais, whole_RMSE, whole_roc, whole_log)

xt <- xtable(whole_tout, digits = c(0,4,4,4,4,4,4,4,4,4,4,4,4,4))

# Affiche la version LaTeX
print(xt, type = "latex")

####### VALID

whole_biais_val <- data.frame()
for(i in method_names_val){
  a <- get(paste0("biais_",i))
  whole_biais_val <- rbind(whole_biais_val, a[5,])
}
rownames(whole_biais_val) <- article_names

whole_RMSE_val <- data.frame()
for(i in method_names_val){
  a <- get(paste0("RMSE_",i))
  whole_RMSE_val <- rbind(whole_RMSE_val, a[5,])
}
rownames(whole_RMSE_val) <- article_names

## roc
whole_roc_val <- data.frame()
for(i in 1:6){
  whole_roc_val <- rbind(whole_roc_val, unlist(unname(ROCmeanval_results[[5]][c(2:6,1)][i])))
}
rownames(whole_roc_val) <- article_names
colnames(whole_roc_val) <- colnames(whole_RMSE_val)

whole_log = as.data.frame(t(t(c(rep(NA,6)))))

whole_tout_val <- cbind(whole_biais_val, whole_RMSE_val, whole_roc_val, whole_log)

xt_val <- xtable(whole_tout_val, digits = c(0,4,4,4,4,4,4,4,4,4,4,4,4,4))

# Affiche la version LaTeX
print(xt_val, type = "latex")
#load("~/Documents/Rstudio/Application/Application colrec/AUC_valid_strata.Rdata")

library(survivalNET)
library(survivalPLANN)
library(relsurv)
library(parallel)

library(ggplot2)
library(dplyr)
library(tidyr)
library(future)
library(future.apply)

# plan(multisession, workers = parallel::detectCores() - 8)

# rm(plann.model)
# rm(plannpred)
# rm(plannpred_HC)
# rm(plannpred_HR)
# rm(plannpred_FR)
# rm(plannpred_FC)
# rm(plannpred_s)
# rm(plannpredval)
# rm(plannpredval_HC)
# rm(plannpredval_HR)
# rm(plannpredval_FR)
# rm(plannpredval_FC)


set.seed(10626)

##cross-validation pour hyper-paramètres 
source("/home/thomas/Documents/Rstudio/Application/Application registre/cvFLEXNETpara.R")
source("/home/thomas/Documents/Rstudio/Application/Application registre/cvPLANNpara.R")

data(colrec)
data(slopop)
pro.time = max(colrec$time)
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



cv <- cvPLANNpara(formula = Surv(time, stat) ~ stage2 + stage3 + agey10 + sex01 + colon, pro.time = pro.time,
                  data = colrec, cv = 10, inter= 365.241/12, size=c(6,12,24,32),
                  decay=c(0.01, 0.05, 0.1), maxit=1000, MaxNWts=10000, metric = "ibs",
                  parallel = TRUE, ncores = parallel::detectCores()-20)


cv2 <- cvFLEXNETpara(formula = Surv(time, stat) ~ stage2 + stage3 + agey10 + sex01 + colon+
                       ratetable(age, diag, sexchara), 
                     ratetable=slopop , pro.time = pro.time,
                     data = colrec, cv = 10, m = c(0,1,2,3),
                     parallel = TRUE, ncores = parallel::detectCores()-20)

cv3 <- cvFLEXNETpara(formula = Surv(time, stat) ~ stage2 + stage3 + agey10 + strata(sex.organ)+
                       ratetable(age, diag, sexchara), ratetable=slopop , pro.time = pro.time,
                     data = colrec, cv = 10, m = c(0,1,2,3), m_s = c(0,1),
                     parallel = TRUE, ncores = parallel::detectCores()-20)


tune.flex <- cv2
tune.flex2 <- cv3
tune.plann <- cv

m1 <- tune.flex$optimal$m

m2 <- tune.flex2$optimal$m
m_s <- tune.flex2$optimal$m_s

inter <- tune.plann$optimal$inter

size <- tune.plann$optimal$size

decay <- tune.plann$optimal$decay

maxit <- tune.plann$optimal$maxit

MaxNWts <- tune.plann$optimal$MaxNWts

save(list = ls(), file = paste0("/home/thomas/Documents/Rstudio/Application/Application colrec/CV_ICdone.Rdata"))
#load("~/Documents/Rstudio/Application/Application colrec/CV_ICdone.Rdata")

#models FULL

#plann

plann.model <-  sPLANN(formula = Surv(time, stat) ~ stage2 + stage3 + agey10 + sex01 + colon,
                       pro.time = pro.time, data = colrec, inter= inter, size=size, 
                       decay=decay, maxit=maxit, MaxNWts=MaxNWts)


plannpred <- predictRS(plann.model, data = colrec, newtimes = newtimespred,
                       ratetable = slopop,  age = "age", year = "diag", sex = "sexchara")

mean.plann <- plannpred$mpredictions$net_survival
meanind.plann <- apply(plannpred$ipredictions$relative_survival, mean, MARGIN = 2)

### flex PH

flex.model1 <- survivalFLEXNET(formula = Surv(time, stat) ~ stage2 + stage3 + agey10 + sex01 + colon +
                                 ratetable(age, diag, sexchara), data = colrec, ratetable = slopop, m = m1) 

flexpred1 <- predict(flex.model1, newdata = colrec, newtimes = newtimespred)$predictions

mean.flex1 <- apply(flexpred1, MARGIN = 2, FUN = mean)


## flex NPH

flex.model2 <-survivalFLEXNET(formula = Surv(time, stat) ~ stage2 + stage3 + agey10 + strata(sex.organ) +
                                ratetable(age, diag, sexchara), data = colrec,
                              ratetable=slopop, m = m2, m_s = m_s)

flexpred2 <- predict(flex.model2, newdata = colrec, newtimes = newtimespred)$predictions

mean.flex2 <- apply(flexpred2, MARGIN = 2, FUN = mean)

save(list = ls(), file = paste0("/home/thomas/Documents/Rstudio/Application/Application colrec/modelsestimates.Rdata"))

#load("~/Documents/Rstudio/Application/Application colrec/modelsestimates.Rdata")
library(survivalNET)
library(survivalPLANN)
library(relsurv)
library(parallel)

dirs <- c(
  "~/Documents/Rstudio/Application/Application colrec/boot/params",
  
  "~/Documents/Rstudio/Application/Application colrec/boot/pred/mean/p",
  "~/Documents/Rstudio/Application/Application colrec/boot/pred/mean/f1",
  "~/Documents/Rstudio/Application/Application colrec/boot/pred/mean/f2",
  "~/Documents/Rstudio/Application/Application colrec/boot/pred/mean/PP",
  
  "~/Documents/Rstudio/Application/Application colrec/boot/pred/ind/p",
  "~/Documents/Rstudio/Application/Application colrec/boot/pred/ind/f1",
  "~/Documents/Rstudio/Application/Application colrec/boot/pred/ind/f2",
  
  "~/Documents/Rstudio/Application/Application colrec/boot/models/p",
  "~/Documents/Rstudio/Application/Application colrec/boot/models/f1",
  "~/Documents/Rstudio/Application/Application colrec/boot/models/f2",
  
  "~/Documents/Rstudio/Application/Application colrec/boot/data"
)

for (d in dirs) {
  if (!dir.exists(d)) {
    dir.create(d, recursive = TRUE)
    message("Created: ", d)
  }
}

B = 1000


set.seed(30626)

boot_results <- function(b){
  
  
  set.seed(b)
  
  tryCatch({
    
    idx <- sample(
      1:nrow(colrec),
      replace = TRUE
    )
    
    dboot <- colrec[idx, ]
    dbootval <- colrec[setdiff(seq_len(nrow(colrec)), unique(idx)),]
    
    pro.time <- max(dboot$time)
    
    ### estimations et prédictions
    plann.boot <- sPLANN(formula = Surv(time, stat) ~ stage2 + stage3 + agey10 + sex01 + colon,
                         pro.time = pro.time, data = dboot, inter= inter, size=size, 
                         decay=decay, maxit=maxit, MaxNWts=MaxNWts)
    
    
    plannpredboot <- predictRS(plann.boot, data = dboot, newtimes = newtimespred,
                               ratetable = slopop,  age = "age", year = "diag", sex = "sexchara")
    
    mean.plannboot <- plannpredboot$mpredictions$net_survival
    meanind.plannboot <- apply(plannpredboot$ipredictions$relative_survival, mean, MARGIN = 2)
    
    plannpredbootval <- predictRS(plann.boot, data = dbootval, newtimes = newtimespred,
                                  ratetable = slopop,  age = "age", year = "diag", sex = "sexchara")
    
    mean.plannbootval <- plannpredbootval$mpredictions$net_survival
    meanind.plannbootval <- apply(plannpredbootval$ipredictions$relative_survival, mean, MARGIN = 2)
    ### flex PH
    
    flex.boot1 <- survivalFLEXNET(
      formula = Surv(time, stat) ~ stage2 + stage3 + agey10 + sex01 + colon +
        ratetable(age, diag, sexchara),
      data = dboot,
      ratetable = slopop,
      m = m1
    )
    
    flexpredboot1 <- predict(flex.boot1, newdata = dboot, newtimes = newtimespred)$predictions
    
    mean.flexboot1 <- apply(flexpredboot1, MARGIN = 2, FUN = mean)
    
    flexpredboot1val <- predict(flex.boot1, newdata = dbootval, newtimes = newtimespred)$predictions
    
    mean.flexboot1val <- apply(flexpredboot1val, MARGIN = 2, FUN = mean)
    
    ## flex NPH
    
    flex.boot2 <-survivalFLEXNET(formula = Surv(time, stat) ~ stage2 + stage3 + agey10 + strata(sex.organ) +
                                   ratetable(age, diag, sexchara), data = dboot,
                                 ratetable=slopop, m = m2, m_s = m_s)
    
    flexpredboot2 <- predict(flex.boot2,newdata = dboot, newtimes = newtimespred)$predictions
    
    mean.flexboot2 <- apply(flexpredboot2, MARGIN = 2, FUN = mean)
    
    flexpredboot2val <- predict(flex.boot2,newdata = dbootval, newtimes = newtimespred)$predictions
    
    mean.flexboot2val <- apply(flexpredboot2val, MARGIN = 2, FUN = mean)
    # list(P = mean.plannboot[-1] , F1 = mean.flexboot1, F2 = mean.flexboot2, Pind = meanind.plannboot[-1]
    # )
    ##Pohar-pred sur le train 
    
    PP <- summary( rs.surv(Surv(time, stat) ~ 1, data = dboot,
                           ratetable = slopop, method = "pohar-perme", add.times = newtimespred,
                           rmap = list(age = age, sex = sexchara, year = diag)), times = newtimespred, extend = TRUE) 
    
    poharpred <- PP$surv 
    
    ##Pohar-Perme sur le oob sample
    
    PPval <- summary( rs.surv(Surv(time, stat) ~ 1, data = dbootval,
                              ratetable = slopop, method = "pohar-perme", add.times = newtimespred,
                              rmap = list(age = age, sex = sexchara, year = diag)), times = newtimespred, extend = TRUE) 
    
    poharpredval <- PPval$surv 
    
    
    
    params <- list(P = tune.plann, F1 = tune.flex, F2 = tune.flex2, idxboot = idx)
    
    saveRDS(params, paste0("~/Documents/Rstudio/Application/Application colrec/boot/params/",b,"_params.Rdata"))
    
    write.table(t(as.matrix(meanind.plannboot[-1])), paste0("~/Documents/Rstudio/Application/Application colrec/boot/pred/mean/p/",b,"_p_mean.csv"), sep = ";", row.names = F, col.names = T)
    write.table(t(as.matrix(mean.flexboot1)), paste0("~/Documents/Rstudio/Application/Application colrec/boot/pred/mean/f1/",b,"_f1_mean.csv"), sep = ";", row.names = F, col.names = T)
    write.table(t(as.matrix(mean.flexboot2)), paste0("~/Documents/Rstudio/Application/Application colrec/boot/pred/mean/f2/",b,"_f2_mean.csv"), sep = ";", row.names = F, col.names = T)
    write.table(t(as.matrix(poharpred)), paste0("~/Documents/Rstudio/Application/Application colrec/boot/pred/mean/PP/",b,"_PP_mean.csv"), sep = ";", row.names = F, col.names = T)
    
    write.table(t(as.matrix(meanind.plannbootval[-1])), paste0("~/Documents/Rstudio/Application/Application colrec/boot/pred/mean/p/",b,"_p_meanval.csv"), sep = ";", row.names = F, col.names = T)
    write.table(t(as.matrix(mean.flexboot1val)), paste0("~/Documents/Rstudio/Application/Application colrec/boot/pred/mean/f1/",b,"_f1_meanval.csv"), sep = ";", row.names = F, col.names = T)
    write.table(t(as.matrix(mean.flexboot2val)), paste0("~/Documents/Rstudio/Application/Application colrec/boot/pred/mean/f2/",b,"_f2_meanval.csv"), sep = ";", row.names = F, col.names = T)
    write.table(t(as.matrix(poharpredval)), paste0("~/Documents/Rstudio/Application/Application colrec/boot/pred/mean/PP/",b,"_PP_meanval.csv"), sep = ";", row.names = F, col.names = T)
    
    write.table(as.matrix(plannpredboot$ipredictions$relative_survival), paste0("~/Documents/Rstudio/Application/Application colrec/boot/pred/ind/p/",b,"_p_ind.csv"), sep = ";", row.names = F, col.names = T)
    write.table(as.matrix(flexpredboot1), paste0("~/Documents/Rstudio/Application/Application colrec/boot/pred/ind/f1/",b,"_f1_ind.csv"), sep = ";", row.names = F, col.names = T)
    write.table(as.matrix(flexpredboot2), paste0("~/Documents/Rstudio/Application/Application colrec/boot/pred/ind/f2/",b,"_f2_ind.csv"), sep = ";", row.names = F, col.names = T)
    
    write.table(as.matrix(plannpredbootval$ipredictions$relative_survival), paste0("~/Documents/Rstudio/Application/Application colrec/boot/pred/ind/p/",b,"_p_indval.csv"), sep = ";", row.names = F, col.names = T)
    write.table(as.matrix(flexpredboot1val), paste0("~/Documents/Rstudio/Application/Application colrec/boot/pred/ind/f1/",b,"_f1_indval.csv"), sep = ";", row.names = F, col.names = T)
    write.table(as.matrix(flexpredboot2val), paste0("~/Documents/Rstudio/Application/Application colrec/boot/pred/ind/f2/",b,"_f2_indval.csv"), sep = ";", row.names = F, col.names = T)
    
    saveRDS(plann.boot, paste0("~/Documents/Rstudio/Application/Application colrec/boot/models/p/",b,"_p_mod.Rdata"))
    saveRDS(flex.boot1, paste0("~/Documents/Rstudio/Application/Application colrec/boot/models/f1/",b,"_f1_mod.Rdata"))
    saveRDS(flex.boot2, paste0("~/Documents/Rstudio/Application/Application colrec/boot/models/f2/",b,"_f2_mod.Rdata"))
    
    write.table(as.matrix(dboot), paste0("~/Documents/Rstudio/Application/Application colrec/boot/data/",b,"_data.csv"), sep = ";", row.names = F, col.names = T)
    write.table(as.matrix(dbootval), paste0("~/Documents/Rstudio/Application/Application colrec/boot/data/",b,"_dataval.csv"), sep = ";", row.names = F, col.names = T)
    
  }, error = function(e) {
    stop(sprintf("Error at iteration %d: %s", b, e$message))
  })
}

mclapply(1:B, boot_results, mc.cores = detectCores() - 4)

B = 1000

ind_m <- c()

path0 <- "/media/thomas/Ultra Touch/boot/" #chemin disque dur externe
# path0 <- "/home/thomas/Documents/Rstudio/Application/Application colrec/boot/" chemin perso

ind_missing <- c(1,64,105,836,317,882,279,487,297,922,783,616,843,251,762,659,894)

for(i in 1:1000){
  if(file.exists(paste0(path0, "pred/mean/f1/",i,"_f1_meanval.csv"))){ind_m <- c(ind_m, i)}
}

ind_m <- ind_m[!ind_m %in% ind_missing]

boot_survP <- matrix(NA, nrow=B, ncol=length(newtimespred))
boot_survF1 <- matrix(NA, nrow=B, ncol=length(newtimespred))
boot_survF2 <- matrix(NA, nrow=B, ncol=length(newtimespred))


for(i in ind_m){
  
  boot_survF1[i,] <- as.matrix(read.csv(paste0(path0, "pred/mean/f1/",i,"_f1_meanval.csv"), sep = ";"))
  boot_survF2[i,] <- as.matrix(read.csv(paste0(path0, "pred/mean/f2/",i,"_f2_meanval.csv"), sep = ";"))
  boot_survP[i,] <- as.matrix(read.csv(paste0(path0, "pred/mean/p/",i,"_p_meanval.csv"), sep = ";"))
  
}

lowerP <- apply(boot_survP, 2, quantile, 0.025, na.rm = TRUE)
upperP <- apply(boot_survP, 2, quantile, 0.975, na.rm = TRUE)

lowerF1 <- apply(boot_survF1, 2, quantile, 0.025, na.rm = TRUE)
upperF1 <- apply(boot_survF1, 2, quantile, 0.975, na.rm = TRUE)

lowerF2 <- apply(boot_survF2, 2, quantile, 0.025, na.rm = TRUE)
upperF2 <- apply(boot_survF2, 2, quantile, 0.975, na.rm = TRUE)

mean_boot_Pval <- apply(boot_survP, mean, MARGIN = 2, na.rm = TRUE)
mean_boot_F1val <- apply(boot_survF1, mean, MARGIN = 2, na.rm = TRUE)
mean_boot_F2val <- apply(boot_survF2, mean, MARGIN = 2, na.rm = TRUE)


plot(newtimespred,mean_boot_Pval, type = 'l', xlim = c(0,3500), ylim = c(0,1))
lines(newtimespred, lowerP, lty = 2)
lines(newtimespred, upperP, lty = 2)

lines(newtimespred,mean_boot_F1val, col = 'red')
lines(newtimespred, lowerF1, lty = 2, col = 'red')
lines(newtimespred, upperF1, lty = 2, col = 'red')

lines(newtimespred,mean_boot_F2val, col = 'green')
lines(newtimespred, lowerF2, lty = 2, col = 'green')
lines(newtimespred, upperF2, lty = 2, col = 'green')

#### échantillon de train
boot_survPtrain <- matrix(NA, nrow=B, ncol=length(newtimespred))
boot_survF1train <- matrix(NA, nrow=B, ncol=length(newtimespred))
boot_survF2train <- matrix(NA, nrow=B, ncol=length(newtimespred))

ind_m <- c()
path0 <- "/media/thomas/Ultra Touch/boot/" #chemin disque dur externe
# path0 <- "/home/thomas/Documents/Rstudio/Application/Application colrec/boot/" chemin perso

for(i in 1:1000){
  if(file.exists(paste0(path0, "pred/mean/f1/",i,"_f1_mean.csv"))){ind_m <- c(ind_m, i)}
}

ind_m <- ind_m[!ind_m %in% ind_missing]

for(i in ind_m){
  
  boot_survF1train[i,] <- as.matrix(read.csv(paste0(path0, "pred/mean/f1/",i,"_f1_mean.csv"), sep = ";"))
  boot_survF2train[i,] <- as.matrix(read.csv(paste0(path0, "pred/mean/f2/",i,"_f2_mean.csv"), sep = ";"))
  boot_survPtrain[i,] <- as.matrix(read.csv(paste0(path0, "pred/mean/p/",i,"_p_mean.csv"), sep = ";"))
  
}

lowerPt <- apply(boot_survPtrain, 2, quantile, 0.025, na.rm = TRUE)
upperPt <- apply(boot_survPtrain, 2, quantile, 0.975, na.rm = TRUE)

lowerF1t <- apply(boot_survF1train, 2, quantile, 0.025, na.rm = TRUE)
upperF1t <- apply(boot_survF1train, 2, quantile, 0.975, na.rm = TRUE)

lowerF2t <- apply(boot_survF2train, 2, quantile, 0.025, na.rm = TRUE)
upperF2t <- apply(boot_survF2train, 2, quantile, 0.975, na.rm = TRUE)

mean_boot_Pt <- apply(boot_survPtrain, mean, MARGIN = 2, na.rm = TRUE)
mean_boot_F1t <- apply(boot_survF1train, mean, MARGIN = 2, na.rm = TRUE)
mean_boot_F2t <- apply(boot_survF2train, mean, MARGIN = 2, na.rm = TRUE)

############################## CALCUL INDICATEURS
library(RISCA)
####AUC

ind_estimP <- plannpred$ipredictions$relative_survival[,-1]
ind_estimF1 <- flexpred1
ind_estimF2 <- flexpred2

newtimes = c(3,5,10)*365.241

for(nt in newtimes){
  
  main_name <- main_vec[match(nt, newtimes)]
  mat <- match(newtimes, newtimes)
  correstab <- setNames(mat, newtimes)  
  # Initialiser les vecteurs de résultats
  rocP <- rocF1.2 <- rocF2.2 <- c()
  
  rocP <- roc.net(colrec$time, colrec$stat, 1-ind_estimP[, as.numeric(correstab[as.character(nt)])],
                  colrec$age, colrec$sexchara, colrec$diag, slopop,
                  pro.time = newtimes[ as.numeric(correstab[as.character(nt)])],
                  cut.off = unique(quantile(1-ind_estimP[,  as.numeric(correstab[as.character(nt)])],
                                            probs = seq(0, 1, .01))))
  
  rocF1.2 <- roc.net(colrec$time, colrec$stat, 1-ind_estimF1[,  as.numeric(correstab[as.character(nt)])],
                     colrec$age, colrec$sexchara, colrec$diag,
                     slopop, pro.time = newtimes[ as.numeric(correstab[as.character(nt)])],
                     cut.off = unique(quantile(1-ind_estimF1[,  as.numeric(correstab[as.character(nt)])],
                                               probs = seq(0, 1, .01))))
  
  rocF2.2 <- roc.net(colrec$time, colrec$stat, 1-ind_estimF2[,  as.numeric(correstab[as.character(nt)])],
                     colrec$age, colrec$sexchara, colrec$diag,
                     slopop, pro.time = newtimes[ as.numeric(correstab[as.character(nt)])],
                     cut.off = unique(quantile(1-ind_estimF2[,  as.numeric(correstab[as.character(nt)])],
                                               probs = seq(0, 1, .01))))
  
  assign(paste0("AUC_whole_",nt), list(
    P = rocP, F1 = rocF1.2, F2 = rocF2.2) )
}

# save(list = ls(), file = paste0("/home/thomas/Documents/Rstudio/Application/Application colrec/ROCtraincalc.Rdata"))

#load("~/Documents/Rstudio/Application/Application colrec/ROCtraincalc.Rdata")

##AUC bootstrap pour IC

boot_AUC <- function(i){
  
  newtimes = c(3,5,10)*365.241
  ind_estimP = read.csv(paste0("/home/thomas/Documents/Rstudio/Application/Application colrec/boot/pred/ind/p/",i,"_p_indval.csv"), sep = ";")[,-1]
  ind_estimF1 = read.csv(paste0("/home/thomas/Documents/Rstudio/Application/Application colrec/boot/pred/ind/f1/",i,"_f1_indval.csv"), sep = ";")
  ind_estimF2 = read.csv(paste0("/home/thomas/Documents/Rstudio/Application/Application colrec/boot/pred/ind/f2/",i,"_f2_indval.csv"), sep = ";")
  data_b = read.csv(paste0("/home/thomas/Documents/Rstudio/Application/Application colrec/boot/data/",i,"_dataval.csv"), sep = ";")
  data_b$diag <- as.Date(data_b$diag)
  
  auc_list <- list()
  
  for(nt in newtimes){
    
    main_name <- main_vec[match(nt, newtimes)]
    mat <- match(newtimes, newtimes)
    correstab <- setNames(mat, newtimes)
    
    rocP <- roc.net(data_b$time, data_b$stat, 1-ind_estimP[, as.numeric(correstab[as.character(nt)])],
                    data_b$age, data_b$sexchara, data_b$diag, slopop,
                    pro.time = newtimes[ as.numeric(correstab[as.character(nt)])],
                    cut.off = unique(quantile(1-ind_estimP[,  as.numeric(correstab[as.character(nt)])],
                                              probs = seq(0, 1, .01))))
    
    rocF1.2 <- roc.net(data_b$time, data_b$stat, 1-ind_estimF1[,  as.numeric(correstab[as.character(nt)])],
                       data_b$age, data_b$sexchara, data_b$diag,
                       slopop, pro.time = newtimes[ as.numeric(correstab[as.character(nt)])],
                       cut.off = unique(quantile(1-ind_estimF1[,  as.numeric(correstab[as.character(nt)])],
                                                 probs = seq(0, 1, .01))))
    
    rocF2.2 <- roc.net(data_b$time, data_b$stat, 1-ind_estimF2[,  as.numeric(correstab[as.character(nt)])],
                       data_b$age, data_b$sexchara, data_b$diag,
                       slopop, pro.time = newtimes[ as.numeric(correstab[as.character(nt)])],
                       cut.off = unique(quantile(1-ind_estimF2[,  as.numeric(correstab[as.character(nt)])],
                                                 probs = seq(0, 1, .01))))
    
    auc_list[[as.character(nt)]] <- list(
      aucP = rocP$auc,
      aucF1.2 = rocF1.2$auc,
      aucF2.2 = rocF2.2$auc,
      rocP = rocP,
      rocF1.2 = rocF1.2,
      rocF2.2 = rocF2.2
    )
    
  }
  
  return(auc_list)
  
}


# boot_AUC_values <- mclapply(ind_m, boot_AUC, mc.cores = detectCores() - 4)

# saveRDS(boot_AUC_values, "~/Documents/Rstudio/Application/Application colrec/auc_boot_calculées.Rdata")
# load("~/Documents/Rstudio/Application/Application colrec/ROCtraincalc.Rdata")
boot_AUC_values <- readRDS("~/Documents/Rstudio/Application/Application colrec/auc_boot_calculées.Rdata")


###### récupération des valeurs de quantiles pour IC

boot_AUC3P <- boot_AUC5P <- boot_AUC10P <- boot_AUC3F1 <- boot_AUC5F1 <- boot_AUC10F1 <- boot_AUC3F2 <- boot_AUC5F2 <- boot_AUC10F2 <- c()

ind_m <- 1:1000

for(i in 1:length(ind_m)){
  
  boot_AUC3P <- rbind(boot_AUC3P, boot_AUC_values[[i]]$`1095.723`$aucP)
  boot_AUC3F1 <- rbind(boot_AUC3F1, boot_AUC_values[[i]]$`1095.723`$aucF1.2)
  boot_AUC3F2 <- rbind(boot_AUC3F2, boot_AUC_values[[i]]$`1095.723`$aucF2.2)
  
  boot_AUC5P <- rbind(boot_AUC5P, boot_AUC_values[[i]]$`1826.205`$aucP)
  boot_AUC5F1 <- rbind(boot_AUC5F1, boot_AUC_values[[i]]$`1826.205`$aucF1.2)
  boot_AUC5F2 <- rbind(boot_AUC5F2, boot_AUC_values[[i]]$`1826.205`$aucF2.2)
  
  boot_AUC10P <- rbind(boot_AUC10P, boot_AUC_values[[i]]$`3652.41`$aucP)
  boot_AUC10F1 <- rbind(boot_AUC10F1, boot_AUC_values[[i]]$`3652.41`$aucF1.2)
  boot_AUC10F2 <- rbind(boot_AUC10F2, boot_AUC_values[[i]]$`3652.41`$aucF2.2)
  
}

ind_m <- ind_m[!ind_m %in% ind_missing]

for(p in c("P", "F1", "F2")){
  
  a <- get(paste0("boot_AUC3",p))
  assign(paste0("boot_AUC3",p), as.matrix(a[-ind_missing]) )
  b <- get(paste0("boot_AUC5",p))
  assign(paste0("boot_AUC5",p), as.matrix(b[-ind_missing]) )
  c <- get(paste0("boot_AUC10",p))
  assign(paste0("boot_AUC10",p), as.matrix(c[-ind_missing]) )
  
}


##PLANN
lowerP3AUC <- apply(boot_AUC3P, 2, quantile, 0.025)
upperP3AUC <- apply(boot_AUC3P, 2, quantile, 0.975)

lowerP5AUC <- apply(boot_AUC5P, 2, quantile, 0.025)
upperP5AUC <- apply(boot_AUC5P, 2, quantile, 0.975)

lowerP10AUC <- apply(boot_AUC10P, 2, quantile, 0.025)
upperP10AUC <- apply(boot_AUC10P, 2, quantile, 0.975)

## F1

lowerF13AUC <- apply(boot_AUC3F1, 2, quantile, 0.025)
upperF13AUC <- apply(boot_AUC3F1, 2, quantile, 0.975)

lowerF15AUC <- apply(boot_AUC5F1, 2, quantile, 0.025)
upperF15AUC <- apply(boot_AUC5F1, 2, quantile, 0.975)

lowerF110AUC <- apply(boot_AUC10F1, 2, quantile, 0.025)
upperF110AUC <- apply(boot_AUC10F1, 2, quantile, 0.975)

## F2

lowerF23AUC <- apply(boot_AUC3F2, 2, quantile, 0.025)
upperF23AUC <- apply(boot_AUC3F2, 2, quantile, 0.975)

lowerF25AUC <- apply(boot_AUC5F2, 2, quantile, 0.025)
upperF25AUC <- apply(boot_AUC5F2, 2, quantile, 0.975)

lowerF210AUC <- apply(boot_AUC10F2, 2, quantile, 0.025)
upperF210AUC <- apply(boot_AUC10F2, 2, quantile, 0.975)


IC_AUC <- list(P = list("3y" = list(low = lowerP3AUC, upp = upperP3AUC) ,"5y" = list(low = lowerP5AUC, upp = upperP5AUC), "10y" = list(low = lowerP10AUC, upp = upperP10AUC)),
               F1 = list("3y" = list(low = lowerF13AUC, upp = upperF13AUC) ,"5y" = list(low = lowerF15AUC, upp = upperF15AUC), "10y" = list(low = lowerF110AUC, upp = upperF110AUC)),
               F2 = list("3y" = list(low = lowerF23AUC, upp = upperF23AUC) ,"5y" = list(low = lowerF25AUC, upp = upperF25AUC), "10y" = list(low = lowerF210AUC, upp = upperF210AUC))
)

#mean AUCboot val

boot_AUCP3val <- apply(boot_AUC3P , mean, MARGIN = 2)
boot_AUCP5val <- apply(boot_AUC5P , mean, MARGIN = 2)
boot_AUCP10val <- apply(boot_AUC10P , mean, MARGIN = 2)

boot_AUCF13val <- apply(boot_AUC3F1 , mean, MARGIN = 2)
boot_AUCF15val <- apply(boot_AUC5F1 , mean, MARGIN = 2)
boot_AUCF110val <- apply(boot_AUC10F1 , mean, MARGIN = 2)

boot_AUCF23val <- apply(boot_AUC3F2 , mean, MARGIN = 2)
boot_AUCF25val <- apply(boot_AUC5F2 , mean, MARGIN = 2)
boot_AUCF210val <- apply(boot_AUC10F2 , mean, MARGIN = 2)


##xtable IC_AUC

dfAUCP <- data.frame(boot_AUCP3val, IC_AUC$P$`3y`$low, IC_AUC$P$`3y`$upp)
colnames(dfAUCP) <- c("estimate", "lower", "upper")
dfAUCP <- rbind(dfAUCP, c(boot_AUCP5val, IC_AUC$P$`5y`$low, IC_AUC$P$`5y`$upp))
dfAUCP <- rbind(dfAUCP, c(boot_AUCP10val, IC_AUC$P$`10y`$low, IC_AUC$P$`10y`$upp))
dfAUCP <- round(dfAUCP, digit = 3)

dfAUCF1 <- data.frame(boot_AUCF13val, IC_AUC$F1$`3y`$low, IC_AUC$F1$`3y`$upp)
colnames(dfAUCF1) <- c("estimate", "lower", "upper")
dfAUCF1 <- rbind(dfAUCF1, c(boot_AUCF15val, IC_AUC$F1$`5y`$low, IC_AUC$F1$`5y`$upp))
dfAUCF1 <- rbind(dfAUCF1, c(boot_AUCF110val, IC_AUC$F1$`10y`$low, IC_AUC$F1$`10y`$upp))
dfAUCF1 <- round(dfAUCF1, digit = 3)

dfAUCF2 <- data.frame(boot_AUCF23val, IC_AUC$F2$`3y`$low, IC_AUC$F2$`3y`$upp)
colnames(dfAUCF2) <- c("estimate", "lower", "upper")
dfAUCF2 <- rbind(dfAUCF2, c(boot_AUCF25val, IC_AUC$F2$`5y`$low, IC_AUC$F2$`5y`$upp))
dfAUCF2 <- rbind(dfAUCF2, c(boot_AUCF210val, IC_AUC$F2$`10y`$low, IC_AUC$F2$`10y`$upp))
dfAUCF2 <- round(dfAUCF2, digit = 3)


dfAUC <- t(data.frame(c( paste0(dfAUCP$estimate[1], " [", dfAUCP$lower[1]," : ", dfAUCP$upper[1], "] ")  ,  paste0(dfAUCP$estimate[2], " [", dfAUCP$lower[2]," : ", dfAUCP$upper[2], "] ") ,  paste0(dfAUCP$estimate[3], " [", dfAUCP$lower[3]," : ", dfAUCP$upper[3], "] ") ))   )
colnames(dfAUC) <- c("3 years", "5 years", "10 years")
dfAUC <- rbind(dfAUC, c( paste0(dfAUCF1$estimate[1], " [", dfAUCF1$lower[1]," : ", dfAUCF1$upper[1], "] ")  ,  paste0(dfAUCF1$estimate[2], " [", dfAUCF1$lower[2]," : ", dfAUCF1$upper[2], "] ") ,  paste0(dfAUCF1$estimate[3], " [", dfAUCF1$lower[3]," : ", dfAUCF1$upper[3], "] ") ))
dfAUC <- rbind(dfAUC, c( paste0(dfAUCF2$estimate[1], " [", dfAUCF2$lower[1]," : ", dfAUCF2$upper[1], "] ")  ,  paste0(dfAUCF2$estimate[2], " [", dfAUCF2$lower[2]," : ", dfAUCF2$upper[2], "] ") ,  paste0(dfAUCF2$estimate[3], " [", dfAUCF2$lower[3]," : ", dfAUCF2$upper[3], "] ") ))

rownames(dfAUC) <- c("P", "F1", "F2")


###############################
###############################
###############################

###############################
###############################
###############################

# ### CALIBRATION
#calibration base entière

newtimes = c(3,5,10)*365.241

for(nt in newtimes){
  
  main_name <- main_vec[match(nt, newtimes)]
  mat <- match(newtimes, newtimespred)
  correstab <- setNames(mat, newtimes)
  
  .predF1 <- flexpred1[, as.numeric(correstab[as.character(nt)])]
  .predF2 <- flexpred2[, as.numeric(correstab[as.character(nt)])]
  .predP <- plannpred$ipredictions$relative_survival[,-1][, as.numeric(correstab[as.character(nt)])]
  
  n.groups <- 5
  
  .grpsF1 <- as.numeric(cut(.predF1, breaks = c(-Inf, quantile(.predF1, seq(1/n.groups, 1, 1/n.groups))),
                            labels = 1:n.groups))
  .grpsF2 <- as.numeric(cut(.predF2, breaks = c(-Inf, quantile(.predF2, seq(1/n.groups, 1, 1/n.groups))),
                            labels = 1:n.groups))
  .grpsP <- as.numeric(cut(.predP, breaks = c(-Inf, quantile(.predP, seq(1/n.groups, 1, 1/n.groups))),
                           labels = 1:n.groups))
  
  
  .estF1 <- sapply(1:n.groups, FUN = function(x) { mean(.predF1[.grpsF1==x]) } )
  .estF2 <- sapply(1:n.groups, FUN = function(x) { mean(.predF2[.grpsF2==x]) } )
  .estP <- sapply(1:n.groups, FUN = function(x) { mean(.predP[.grpsP==x]) } )
  
  
  time <- colrec$time
  event <- colrec$stat
  
  .dataF1 <- data.frame(time = time, event = event, age = colrec$age, sexchara = colrec$sexchara, diag = colrec$diag, grps = .grpsF1)
  .dataF2 <- data.frame(time = time, event = event, age = colrec$age, sexchara = colrec$sexchara, diag = colrec$diag, grps = .grpsF2)
  .dataP <- data.frame(time = time, event = event, age = colrec$age, sexchara = colrec$sexchara, diag = colrec$diag, grps = .grpsP)
  
  
  .survfitF1 <- summary( rs.surv(Surv(time, event) ~ .grpsF1, data = .dataF1,
                                 ratetable = slopop, method = "pohar-perme", add.times = nt,
                                 rmap = list(age = age, sex = sexchara, year = diag)), times = nt, extend = TRUE)
  .survfitF2 <- summary( rs.surv(Surv(time, event) ~ .grpsF2, data = .dataF2,
                                 ratetable = slopop, method = "pohar-perme", add.times = nt,
                                 rmap = list(age = age, sex = sexchara, year = diag)), times = nt, extend = TRUE)
  .survfitP <- summary( rs.surv(Surv(time, event) ~ .grpsP, data = .dataP,
                                ratetable = slopop, method = "pohar-perme", add.times = nt,
                                rmap = list(age = age, sex = sexchara, year = diag)), times = nt, extend = TRUE)
  
  .obsF1 <- .survfitF1$surv
  
  .obsF2 <- .survfitF2$surv
  
  .obsP <- .survfitP$surv
  
  
  df_plot <- bind_rows(
    data.frame(Method = "Spline PH", est = .estF1, obs = .obsF1
    ),
    data.frame(Method = "Spline NPH", est = .estF2, obs = .obsF2
    ),
    data.frame(Method = "PLANN", est = .estP, obs = .obsP
    )
  )
  assign(paste0("calib_W_T_",nt) , df_plot ,envir = globalenv())
}

# save(list = ls(), file = paste0("/home/thomas/Documents/Rstudio/Application/Application colrec/ROCetCALItrain.Rdata"))
#load("~/Documents/Rstudio/Application/Application colrec/ROCetCALItrain.Rdata")

##calibration bootstrap pour IC
##Attention, il faut les predictions enregistrées sur le disque dur externe
boot_cali <- function(i){
  path0 <- "/media/thomas/Ultra Touch/boot/" #chemin disque dur externe
  # path0 <- "/home/thomas/Documents/Rstudio/Application/Application colrec/boot/" chemin perso
  newtimes = c(3,5,10)*365.241
  pest = read.csv(paste0(path0, "pred/ind/p/",i,"_p_indval.csv"), sep = ";")
  f1est = read.csv(paste0(path0, "pred/ind/f1/",i,"_f1_indval.csv"), sep = ";")
  f2est = read.csv(paste0(path0, "pred/ind/f2/",i,"_f2_indval.csv"), sep = ";")
  data_b = read.csv(paste0(path0, "data/",i,"_dataval.csv"), sep = ";")
  data_b$diag <- as.Date(data_b$diag)
  
  time_list <- list()
  
  for(nt in newtimes){
    
    main_name <- main_vec[match(nt, newtimes)]
    mat <- match(newtimes, newtimespred)
    correstab <- setNames(mat, newtimes)
    
    .predF1 <- f1est[, as.numeric(correstab[as.character(nt)])]
    .predF2 <- f2est[, as.numeric(correstab[as.character(nt)])]
    .predP <- pest[, as.numeric(correstab[as.character(nt)])]
    
    n.groups <- 5
    
    .grpsF1 <- as.numeric(cut(.predF1, breaks = c(-Inf, quantile(.predF1, seq(1/n.groups, 1, 1/n.groups))),
                              labels = 1:n.groups))
    .grpsF2 <- as.numeric(cut(.predF2, breaks = c(-Inf, quantile(.predF2, seq(1/n.groups, 1, 1/n.groups))),
                              labels = 1:n.groups))
    .grpsP <- as.numeric(cut(.predP, breaks = c(-Inf, quantile(.predP, seq(1/n.groups, 1, 1/n.groups))),
                             labels = 1:n.groups))
    
    
    .estF1 <- sapply(1:n.groups, FUN = function(x) { mean(.predF1[.grpsF1==x]) } )
    .estF2 <- sapply(1:n.groups, FUN = function(x) { mean(.predF2[.grpsF2==x]) } )
    .estP <- sapply(1:n.groups, FUN = function(x) { mean(.predP[.grpsP==x]) } )
    
    
    time <- data_b$time
    event <- data_b$stat
    
    .dataF1 <- data.frame(time = time, event = event, age = data_b$age, sexchara = data_b$sexchara, diag = data_b$diag, grps = .grpsF1)
    .dataF2 <- data.frame(time = time, event = event, age = data_b$age, sexchara = data_b$sexchara, diag = data_b$diag, grps = .grpsF2)
    .dataP <- data.frame(time = time, event = event, age = data_b$age, sexchara = data_b$sexchara, diag = data_b$diag, grps = .grpsP)
    
    
    .survfitF1 <- summary( rs.surv(Surv(time, event) ~ .grpsF1, data = .dataF1,
                                   ratetable = slopop, method = "pohar-perme", add.times = nt,
                                   rmap = list(age = age, sex = sexchara, year = diag)), times = nt, extend = TRUE)
    .survfitF2 <- summary( rs.surv(Surv(time, event) ~ .grpsF2, data = .dataF2,
                                   ratetable = slopop, method = "pohar-perme", add.times = nt,
                                   rmap = list(age = age, sex = sexchara, year = diag)), times = nt, extend = TRUE)
    .survfitP <- summary( rs.surv(Surv(time, event) ~ .grpsP, data = .dataP,
                                  ratetable = slopop, method = "pohar-perme", add.times = nt,
                                  rmap = list(age = age, sex = sexchara, year = diag)), times = nt, extend = TRUE)
    
    .obsF1 <- .survfitF1$surv
    
    .obsF2 <- .survfitF2$surv
    
    .obsP <- .survfitP$surv
    
    res <- list(p = list(ep = .estP, op = .obsP), f1 = list(ef1 = .estF1, of1 = .obsF1), f2 = list(ef2 = .estF2, of2 = .obsF2))
    time_list[[as.character(nt)]] <- res
  }
  
  return(time_list)
  
}


boot_cali_values <- mclapply(ind_m, boot_cali, mc.cores = detectCores() - 4)

saveRDS(boot_cali_values,"~/Documents/Rstudio/Application/Application colrec/cali_boot_calculées.Rdata")

# load("~/Documents/Rstudio/Application/Application colrec/ROCetCALItrain.Rdata")

boot_AUC_values <- readRDS("~/Documents/Rstudio/Application/Application colrec/auc_boot_calculées.Rdata")
boot_cali_values <- readRDS("~/Documents/Rstudio/Application/Application colrec/cali_boot_calculées.Rdata")

###### récupération des valeurs de quantiles pour IC

boot_cali3P <- boot_cali5P <- boot_cali10P <- boot_cali3F1 <- boot_cali5F1 <- boot_cali10F1 <- boot_cali3F2 <- boot_cali5F2 <- boot_cali10F2 <- c()


for(i in 1:length(ind_m)){
  
  boot_cali3P <- rbind(boot_cali3P, boot_cali_values[[i]]$`1095.723`$p$ep)
  boot_cali3F1 <- rbind(boot_cali3F1, boot_cali_values[[i]]$`1095.723`$f1$ef1)
  boot_cali3F2 <- rbind(boot_cali3F2, boot_cali_values[[i]]$`1095.723`$f2$ef2)
  
  boot_cali5P <- rbind(boot_cali5P, boot_cali_values[[i]]$`1826.205`$p$ep)
  boot_cali5F1 <- rbind(boot_cali5F1, boot_cali_values[[i]]$`1826.205`$f1$ef1)
  boot_cali5F2 <- rbind(boot_cali5F2, boot_cali_values[[i]]$`1826.205`$f2$ef2)
  
  boot_cali10P <- rbind(boot_cali10P, boot_cali_values[[i]]$`3652.41`$p$ep)
  boot_cali10F1 <- rbind(boot_cali10F1, boot_cali_values[[i]]$`3652.41`$f1$ef1)
  boot_cali10F2 <- rbind(boot_cali10F2, boot_cali_values[[i]]$`3652.41`$f2$ef2)
  
}


for(p in c("P", "F1", "F2")){
  
  a <- get(paste0("boot_cali3",p))
  assign(paste0("boot_cali3",p), as.matrix(a[-ind_missing,]) )
  b <- get(paste0("boot_cali5",p))
  assign(paste0("boot_cali5",p), as.matrix(b[-ind_missing,]) )
  c <- get(paste0("boot_cali10",p))
  assign(paste0("boot_cali10",p), as.matrix(c[-ind_missing,]) )
  
}


##PLANN
lowerP3cali <- apply(boot_cali3P, 2, quantile, 0.025)
upperP3cali <- apply(boot_cali3P, 2, quantile, 0.975)

lowerP5cali <- apply(boot_cali5P, 2, quantile, 0.025)
upperP5cali <- apply(boot_cali5P, 2, quantile, 0.975)

lowerP10cali <- apply(boot_cali10P, 2, quantile, 0.025)
upperP10cali <- apply(boot_cali10P, 2, quantile, 0.975)

## F1

lowerF13cali <- apply(boot_cali3F1, 2, quantile, 0.025)
upperF13cali <- apply(boot_cali3F1, 2, quantile, 0.975)

lowerF15cali <- apply(boot_cali5F1, 2, quantile, 0.025)
upperF15cali <- apply(boot_cali5F1, 2, quantile, 0.975)

lowerF110cali <- apply(boot_cali10F1, 2, quantile, 0.025)
upperF110cali <- apply(boot_cali10F1, 2, quantile, 0.975)

## F2

lowerF23cali <- apply(boot_cali3F2, 2, quantile, 0.025)
upperF23cali <- apply(boot_cali3F2, 2, quantile, 0.975)

lowerF25cali <- apply(boot_cali5F2, 2, quantile, 0.025)
upperF25cali <- apply(boot_cali5F2, 2, quantile, 0.975)

lowerF210cali <- apply(boot_cali10F2, 2, quantile, 0.025)
upperF210cali <- apply(boot_cali10F2, 2, quantile, 0.975)


IC_cali <- list(P = list("3y" = list(low = lowerP3cali, upp = upperP3cali) ,"5y" = list(low = lowerP5cali, upp = upperP5cali), "10y" = list(low = lowerP10cali, upp = upperP10cali)),
                F1 = list("3y" = list(low = lowerF13cali, upp = upperF13cali) ,"5y" = list(low = lowerF15cali, upp = upperF15cali), "10y" = list(low = lowerF110cali, upp = upperF110cali)),
                F2 = list("3y" = list(low = lowerF23cali, upp = upperF23cali) ,"5y" = list(low = lowerF25cali, upp = upperF25cali), "10y" = list(low = lowerF210cali, upp = upperF210cali))
)

boot_caliP3val <- apply(boot_cali3P , mean, MARGIN = 2)
boot_caliP5val <- apply(boot_cali5P , mean, MARGIN = 2)
boot_caliP10val <- apply(boot_cali10P , mean, MARGIN = 2)

boot_caliF13val <- apply(boot_cali3F1 , mean, MARGIN = 2)
boot_caliF15val <- apply(boot_cali5F1 , mean, MARGIN = 2)
boot_caliF110val <- apply(boot_cali10F1 , mean, MARGIN = 2)

boot_caliF23val <- apply(boot_cali3F2 , mean, MARGIN = 2)
boot_caliF25val <- apply(boot_cali5F2 , mean, MARGIN = 2)
boot_caliF210val <- apply(boot_cali10F2 , mean, MARGIN = 2)

##xtable IC_cali

dfcaliP3 <- data.frame(boot_caliP3val, IC_cali$P$`3y`$low, IC_cali$P$`3y`$upp)
colnames(dfcaliP3) <- c("estimate", "lower", "upper")
dfcaliP5 <-  data.frame(boot_caliP5val, IC_cali$P$`5y`$low, IC_cali$P$`5y`$upp)
colnames(dfcaliP5) <- c("estimate", "lower", "upper")
dfcaliP10 <- data.frame(boot_caliP10val, IC_cali$P$`10y`$low, IC_cali$P$`10y`$upp)
colnames(dfcaliP10) <- c("estimate", "lower", "upper")

dfcaliP <- cbind(dfcaliP3, dfcaliP5, dfcaliP10)
dfcaliP <- round(dfcaliP, digit = 3)
colnames(dfcaliP) <- c("3y estimate", "3y lower", "3y upper", "5y estimate", "5y lower", "5y upper", "10y estimate", "10y lower", "10y upper")

dfcaliF13 <- data.frame(boot_caliF13val, IC_cali$F1$`3y`$low, IC_cali$F1$`3y`$upp)
colnames(dfcaliF13) <- c("estimate", "lower", "upper")
dfcaliF15 <- data.frame(boot_caliF15val, IC_cali$F1$`5y`$low, IC_cali$F1$`5y`$upp)
colnames(dfcaliF15) <- c("estimate", "lower", "upper")
dfcaliF110 <- data.frame(boot_caliF110val, IC_cali$F1$`10y`$low, IC_cali$F1$`10y`$upp)
colnames(dfcaliF110) <- c("estimate", "lower", "upper")

dfcaliF1 <- cbind(dfcaliF13, dfcaliF15, dfcaliF110)
dfcaliF1 <- round(dfcaliF1, digit = 3)
colnames(dfcaliF1) <- c("3y estimate", "3y lower", "3y upper", "5y estimate", "5y lower", "5y upper", "10y estimate", "10y lower", "10y upper")

dfcaliF23 <- data.frame(boot_caliF23val, IC_cali$F2$`3y`$low, IC_cali$F2$`3y`$upp)
colnames(dfcaliF23) <- c("estimate", "lower", "upper")
dfcaliF25 <- data.frame(boot_caliF25val, IC_cali$F2$`5y`$low, IC_cali$F2$`5y`$upp)
colnames(dfcaliF25) <- c("estimate", "lower", "upper")
dfcaliF210 <- data.frame(boot_caliF210val, IC_cali$F2$`10y`$low, IC_cali$F2$`10y`$upp)
colnames(dfcaliF210) <- c("estimate", "lower", "upper")


dfcaliF2 <- cbind(dfcaliF23, dfcaliF25, dfcaliF210)
dfcaliF2 <- round(dfcaliF2, digit = 3)
colnames(dfcaliF2) <- c("3y estimate", "3y lower", "3y upper", "5y estimate", "5y lower", "5y upper", "10y estimate", "10y lower", "10y upper")

dfcali <- rbind(dfcaliP, dfcaliF1, dfcaliF2)
dfcali <- cbind(c(rep("PLANN", 5), rep("Spline PH", 5), rep("Spline NPH", 5)), dfcali)
colnames(dfcali)[1] <- "Method"

tmpP <- data.frame(
  Y3  = paste0(dfcaliP$`3y estimate`,  " [", dfcaliP$`3y lower`,  " : ", dfcaliP$`3y upper`,  "]"),
  Y5  = paste0(dfcaliP$`5y estimate`,  " [", dfcaliP$`5y lower`,  " : ", dfcaliP$`5y upper`,  "]"),
  Y10 = paste0(dfcaliP$`10y estimate`, " [", dfcaliP$`10y lower`, " : ", dfcaliP$`10y upper`, "]")
)

tmpF1 <- data.frame(
  Y3  = paste0(dfcaliF1$`3y estimate`,  " [", dfcaliF1$`3y lower`,  " : ", dfcaliF1$`3y upper`,  "]"),
  Y5  = paste0(dfcaliF1$`5y estimate`,  " [", dfcaliF1$`5y lower`,  " : ", dfcaliF1$`5y upper`,  "]"),
  Y10 = paste0(dfcaliF1$`10y estimate`, " [", dfcaliF1$`10y lower`, " : ", dfcaliF1$`10y upper`, "]")
)

tmpF2 <- data.frame(
  Y3  = paste0(dfcaliF2$`3y estimate`,  " [", dfcaliF2$`3y lower`,  " : ", dfcaliF2$`3y upper`,  "]"),
  Y5  = paste0(dfcaliF2$`5y estimate`,  " [", dfcaliF2$`5y lower`,  " : ", dfcaliF2$`5y upper`,  "]"),
  Y10 = paste0(dfcaliF2$`10y estimate`, " [", dfcaliF2$`10y lower`, " : ", dfcaliF2$`10y upper`, "]")
)

dfcali_sim <- rbind(tmpP, tmpF1, tmpF2)
dfcali_sim <- cbind(c(rep("PLANN", 5), rep("Spline PH", 5), rep("Spline NPH", 5)) , dfcali_sim)
colnames(dfcali_sim)[1] <- "Method"

# save(list = ls(), file = paste0("/home/thomas/Documents/Rstudio/Application/Application colrec/tableAUCetCALI.Rdata"))
load("~/Documents/Rstudio/Application/Application colrec/tableAUCetCALI.Rdata")


### table moyenne a 3 5 et 10 ans

##train
mat <- match(newtimes, newtimespred)
correstab <- setNames(mat, newtimes)

dfMEANt <- data.frame(rbind( paste0(round( mean_boot_Pt[correstab], digit = 4 ) , " [", round(lowerPt[correstab], digit = 4 ),  " : ", round(upperPt[correstab], digit = 4 ),  "]"), 
                             paste0(round (mean_boot_F1t[correstab], digit = 4 ), " [", round(lowerF1t[correstab], digit = 4 ),  " : ", round(upperF1t[correstab], digit = 4 ),  "]"), 
                             paste0(round( mean_boot_F2t[correstab], digit = 4 ), " [", round(lowerF2t[correstab], digit = 4 ),  " : ", round(upperF2t[correstab], digit = 4 ),  "]")))

##valid

dfMEANv <- data.frame(rbind( paste0(round( mean_boot_Pval[correstab], digit = 4 ) , " [", round(lowerP[correstab], digit = 4 ),  " : ", round(upperP[correstab], digit = 4 ),  "]"), 
                             paste0(round (mean_boot_F1val[correstab], digit = 4 ), " [", round(lowerF1[correstab], digit = 4 ),  " : ", round(upperF1[correstab], digit = 4 ),  "]"), 
                             paste0(round( mean_boot_F2val[correstab], digit = 4 ), " [", round(lowerF2[correstab], digit = 4 ),  " : ", round(upperF2[correstab], digit = 4 ),  "]")))


# save(list = ls(), file = paste0("/home/thomas/Documents/Rstudio/Application/Application colrec/tableAUCetCALIetMEAN.Rdata"))
# load("~/Documents/Rstudio/Application/Application colrec/tableAUCetCALIetMEAN.Rdata") 

View(dfMEANt)
View(dfMEANv)
View(dfAUC)
View(dfcali_sim)

############CALCUL des survies nettes par sous groupes

library(survivalNET)
library(survivalPLANN)
library(relsurv)
library(parallel)

library(ggplot2)
library(dplyr)
library(tidyr)
library(future)
library(future.apply)

strata_calc <- function(i){
  
  pmod <- readRDS(paste0(path0,"models/p/",i,"_p_mod.Rdata"))
  f1mod <- readRDS(paste0(path0,"models/f1/",i,"_f1_mod.Rdata"))
  f2mod <- readRDS(paste0(path0,"models/f2/",i,"_f2_mod.Rdata"))
  
  dt <- read.csv(paste0(path0, "data/",i,"_data.csv"), sep = ";")
  dv <- read.csv(paste0(path0, "data/",i,"_dataval.csv"), sep = ";")
  dt$diag <- as.Date(dt$diag)
  dv$diag <- as.Date(dv$diag)
  
  ##sous strates preds
  
  strata_names <- c("HC", "HR", "FC", "FR")
  
  dtHC <- dt[dt$sex.organ == 2,]
  dtHR <- dt[dt$sex.organ == 1,]
  dtFC <- dt[dt$sex.organ == 4,]
  dtFR <- dt[dt$sex.organ == 3,]
  #val
  dvHC <- dv[dv$sex.organ == 2,]
  dvHR <- dv[dv$sex.organ == 1,]
  dvFC <- dv[dv$sex.organ == 4,]
  dvFR <- dv[dv$sex.organ == 3,]
  
  for(j in strata_names){
    
    dtj <- get(paste0("dt",j))
    dvj <- get(paste0("dv",j))
    
    plannpredS <- predictRS(pmod, data = dtj, newtimes = newtimespred,
                            ratetable = slopop, age = "age", year = "diag", sex = "sexchara")
    
    mean.plannS <- plannpredS$mpredictions$net_survival
    meanind.plannS <- apply(plannpredS$ipredictions$relative_survival, mean, MARGIN = 2)
    
    plannpredvalS <- predictRS(pmod, data = dvj, newtimes = newtimespred,
                               ratetable = slopop, age = "age", year = "diag", sex = "sexchara")
    
    mean.plannvalS <- plannpredvalS$mpredictions$net_survival
    meanind.plannvalS <- apply(plannpredvalS$ipredictions$relative_survival, mean, MARGIN = 2)
    
    flexpred1S <- predict(f1mod, newdata = dtj, newtimes = newtimespred)$predictions
    flexpred1valS <- predict(f1mod, newdata = dvj, newtimes = newtimespred)$predictions
    
    mean.flex1S <- apply(flexpred1S, MARGIN = 2, FUN = mean)
    mean.flex1valS <- apply(flexpred1valS, MARGIN = 2, FUN = mean)
    
    flexpred2S <- predict(f2mod, newdata = dtj, newtimes = newtimespred)$predictions
    flexpred2valS <- predict(f2mod, newdata = dvj, newtimes = newtimespred)$predictions
    
    mean.flex2S <- apply(flexpred2S, MARGIN = 2, FUN = mean)
    mean.flex2valS <- apply(flexpred2valS, MARGIN = 2, FUN = mean)
    
    PPS <- summary( rs.surv(Surv(time, stat) ~ 1, data = dtj, ratetable = slopop, method = "pohar-perme", add.times = newtimespred,
                            rmap = list(age = age, sex = sexchara, year = diag)), times = newtimespred, extend = TRUE)
    poharpredS <- PPS$surv
    
    PPvalS <- summary( rs.surv(Surv(time, stat) ~ 1, data = dvj, ratetable = slopop, method = "pohar-perme", add.times = newtimespred,
                               rmap = list(age = age, sex = sexchara, year = diag)), times = newtimespred, extend = TRUE)
    
    poharpredvalS <- PPvalS$surv
    
    write.table(t(as.matrix(meanind.plannS[-1])), paste0(path0, "pred_strata/mean/p/",i,"_p_mean_",j,".csv"), sep = ";", row.names = F, col.names = T)
    write.table(t(as.matrix(mean.flex1S)), paste0(path0, "pred_strata/mean/f1/",i,"_f1_mean_",j,".csv"), sep = ";", row.names = F, col.names = T)
    write.table(t(as.matrix(mean.flex2S)), paste0(path0, "pred_strata/mean/f2/",i,"_f2_mean_",j,".csv"), sep = ";", row.names = F, col.names = T)
    write.table(t(as.matrix(poharpredS)), paste0(path0, "pred_strata/mean/PP/",i,"_PP_mean_",j,".csv"), sep = ";", row.names = F, col.names = T)
    
    write.table(t(as.matrix(meanind.plannvalS[-1])), paste0(path0, "pred_strata/mean/p/",i,"_p_meanval_",j,".csv"), sep = ";", row.names = F, col.names = T)
    write.table(t(as.matrix(mean.flex1valS)), paste0(path0, "pred_strata/mean/f1/",i,"_f1_meanval_",j,".csv"), sep = ";", row.names = F, col.names = T)
    write.table(t(as.matrix(mean.flex2valS)), paste0(path0, "pred_strata/mean/f2/",i,"_f2_meanval_",j,".csv"), sep = ";", row.names = F, col.names = T)
    write.table(t(as.matrix(poharpredvalS)), paste0(path0, "pred_strata/mean/PP/",i,"_PP_meanval_",j,".csv"), sep = ";", row.names = F, col.names = T)
    
    write.table(as.matrix(plannpredS$ipredictions$relative_survival), paste0(path0, "pred_strata/ind/p/",i,"_p_ind_",j,".csv"), sep = ";", row.names = F, col.names = T)
    write.table(as.matrix(flexpred1S), paste0(path0, "pred_strata/ind/f1/",i,"_f1_ind_",j,".csv"), sep = ";", row.names = F, col.names = T)
    write.table(as.matrix(flexpred2S), paste0(path0, "pred_strata/ind/f2/",i,"_f2_ind_",j,".csv"), sep = ";", row.names = F, col.names = T)
    
    write.table(as.matrix(plannpredvalS$ipredictions$relative_survival), paste0(path0, "pred_strata/ind/p/",i,"_p_indval_",j,".csv"), sep = ";", row.names = F, col.names = T)
    write.table(as.matrix(flexpred1valS), paste0(path0, "pred_strata/ind/f1/",i,"_f1_indval_",j,".csv"), sep = ";", row.names = F, col.names = T)
    write.table(as.matrix(flexpred2valS), paste0(path0, "pred_strata/ind/f2/",i,"_f2_indval_",j,".csv"), sep = ";", row.names = F, col.names = T)
    
  }
  
  
}

path0 <- "/media/thomas/Ultra Touch/boot/" #chemin disque dur externe

ind_m_s <- c()
for(i in 1:1000){
  if(!file.exists(paste0(path0,"pred_strata/mean/p/",i,"_p_mean_","HC",".csv"))){ind_m_s <- c(ind_m_s, i)}
}

mclapply(ind_m_s, strata_calc, mc.cores = detectCores() - 10)


##### calcul moyenne et ic
##
strata_names <- c("HC", "HR", "FC", "FR")

for(j in strata_names){
  
  assign(paste0("boot_survF1_",j), matrix(NA, nrow=B, ncol=length(newtimespred)))
  assign(paste0("boot_survF2_",j), matrix(NA, nrow=B, ncol=length(newtimespred)))
  assign(paste0("boot_survP_",j), matrix(NA, nrow=B, ncol=length(newtimespred)))
  
}


ind_m <- ind_m[!ind_m %in% ind_missing]

for(i in ind_m){
  for(j in strata_names){
    
    a <- get(paste0("boot_survF1_",j))
    b <- get(paste0("boot_survF2_",j))
    c <- get(paste0("boot_survP_",j))
    
    a[i,] <- as.matrix(read.csv(paste0(path0, "pred_strata/mean/f1/",i,"_f1_meanval_",j,".csv"), sep = ";"))
    b[i,] <- as.matrix(read.csv(paste0(path0, "pred_strata/mean/f2/",i,"_f2_meanval_",j,".csv"), sep = ";"))
    c[i,] <- as.matrix(read.csv(paste0(path0, "pred_strata/mean/p/",i,"_p_meanval_",j,".csv"), sep = ";"))
    
    assign(paste0("boot_survF1_",j), a)
    assign(paste0("boot_survF2_",j), b)
    assign(paste0("boot_survP_",j), c)
    
  }
}



for(j in strata_names){
  assign(paste0("lowerP_",j), apply(get(paste0("boot_survP_",j)), 2, quantile, 0.025, na.rm = TRUE))
  assign(paste0("upperP_",j), apply(get(paste0("boot_survP_",j)), 2, quantile, 0.975, na.rm = TRUE))
  
  assign(paste0("lowerF1_",j), apply(get(paste0("boot_survF1_",j)), 2, quantile, 0.025, na.rm = TRUE))
  assign(paste0("upperF1_",j), apply(get(paste0("boot_survF1_",j)), 2, quantile, 0.975, na.rm = TRUE))
  
  assign(paste0("lowerF2_",j), apply(get(paste0("boot_survF2_",j)), 2, quantile, 0.025, na.rm = TRUE))
  assign(paste0("upperF2_",j), apply(get(paste0("boot_survF2_",j)), 2, quantile, 0.975, na.rm = TRUE))
  
  assign(paste0("mean_boot_Pval_",j), apply(get(paste0("boot_survP_",j)), mean, MARGIN = 2, na.rm = TRUE))
  assign(paste0("mean_boot_F1val_",j), apply(get(paste0("boot_survF1_",j)), mean, MARGIN = 2, na.rm = TRUE))
  assign(paste0("mean_boot_F2val_",j), apply(get(paste0("boot_survF2_",j)), mean, MARGIN = 2, na.rm = TRUE))
}

#### échantillon de train
for(j in strata_names){
  
  assign(paste0("boot_survF1train_",j), matrix(NA, nrow=B, ncol=length(newtimespred)))
  assign(paste0("boot_survF2train_",j), matrix(NA, nrow=B, ncol=length(newtimespred)))
  assign(paste0("boot_survPtrain_",j), matrix(NA, nrow=B, ncol=length(newtimespred)))
  
}


for(i in ind_m){
  for(j in strata_names){
    
    a <- get(paste0("boot_survF1train_",j))
    b <- get(paste0("boot_survF2train_",j))
    c <- get(paste0("boot_survPtrain_",j))
    
    a[i,] <- as.matrix(read.csv(paste0(path0, "pred_strata/mean/f1/",i,"_f1_mean_",j,".csv"), sep = ";"))
    b[i,] <- as.matrix(read.csv(paste0(path0, "pred_strata/mean/f2/",i,"_f2_mean_",j,".csv"), sep = ";"))
    c[i,] <- as.matrix(read.csv(paste0(path0, "pred_strata/mean/p/",i,"_p_mean_",j,".csv"), sep = ";"))
    
    assign(paste0("boot_survF1train_",j), a)
    assign(paste0("boot_survF2train_",j), b)
    assign(paste0("boot_survPtrain_",j), c)
    
  }
}


for(j in strata_names){
  assign(paste0("lowerPt_",j), apply(get(paste0("boot_survPtrain_",j)), 2, quantile, 0.025, na.rm = TRUE))
  assign(paste0("upperPt_",j), apply(get(paste0("boot_survPtrain_",j)), 2, quantile, 0.975, na.rm = TRUE))
  
  assign(paste0("lowerF1t_",j), apply(get(paste0("boot_survF1train_",j)), 2, quantile, 0.025, na.rm = TRUE))
  assign(paste0("upperF1t_",j), apply(get(paste0("boot_survF1train_",j)), 2, quantile, 0.975, na.rm = TRUE))
  
  assign(paste0("lowerF2t_",j), apply(get(paste0("boot_survF2train_",j)), 2, quantile, 0.025, na.rm = TRUE))
  assign(paste0("upperF2t_",j), apply(get(paste0("boot_survF2train_",j)), 2, quantile, 0.975, na.rm = TRUE))
  
  assign(paste0("mean_boot_Pt_",j), apply(get(paste0("boot_survPtrain_",j)), mean, MARGIN = 2, na.rm = TRUE))
  assign(paste0("mean_boot_F1t_",j), apply(get(paste0("boot_survF1train_",j)), mean, MARGIN = 2, na.rm = TRUE))
  assign(paste0("mean_boot_F2t_",j), apply(get(paste0("boot_survF2train_",j)), mean, MARGIN = 2, na.rm = TRUE))
}


# save(list = ls(), file = paste0("/home/thomas/Documents/Rstudio/Application/Application colrec/valstrates.Rdata"))
load("~/Documents/Rstudio/Application/Application colrec/valstrates.Rdata")


###############################################################################################################
###############################################################################################################
##################" récupération des valeurs de Pohar Perme pour les forests plots de calibration 
###############################################################################################################
###############################################################################################################


boot_caliPP3P <- boot_caliPP5P <- boot_caliPP10P <- boot_caliPP3F1 <- boot_caliPP5F1 <- boot_caliPP10F1 <- boot_caliPP3F2 <- boot_caliPP5F2 <- boot_caliPP10F2 <- c()


for(i in 1:length(ind_m)){
  
  boot_caliPP3P <- rbind(boot_caliPP3P, boot_cali_values[[i]]$`1095.723`$p$op)
  boot_caliPP3F1 <- rbind(boot_caliPP3F1, boot_cali_values[[i]]$`1095.723`$f1$of1)
  boot_caliPP3F2 <- rbind(boot_caliPP3F2, boot_cali_values[[i]]$`1095.723`$f2$of2)
  
  boot_caliPP5P <- rbind(boot_caliPP5P, boot_cali_values[[i]]$`1826.205`$p$op)
  boot_caliPP5F1 <- rbind(boot_caliPP5F1, boot_cali_values[[i]]$`1826.205`$f1$of1)
  boot_caliPP5F2 <- rbind(boot_caliPP5F2, boot_cali_values[[i]]$`1826.205`$f2$of2)
  
  boot_caliPP10P <- rbind(boot_caliPP10P, boot_cali_values[[i]]$`3652.41`$p$op)
  boot_caliPP10F1 <- rbind(boot_caliPP10F1, boot_cali_values[[i]]$`3652.41`$f1$of1)
  boot_caliPP10F2 <- rbind(boot_caliPP10F2, boot_cali_values[[i]]$`3652.41`$f2$of2)
  
}


for(p in c("P", "F1", "F2")){
  
  a <- get(paste0("boot_caliPP3",p))
  assign(paste0("boot_caliPP3",p), as.matrix(a[-ind_missing,]) )
  b <- get(paste0("boot_caliPP5",p))
  assign(paste0("boot_caliPP5",p), as.matrix(b[-ind_missing,]) )
  c <- get(paste0("boot_caliPP10",p))
  assign(paste0("boot_caliPP10",p), as.matrix(c[-ind_missing,]) )
  
}

##différence estimé-observé 

hold_names <- c("3P", "5P", "10P", "3F1", "5F1", "10F1", "3F2", "5F2", "10F2")

for(j in hold_names){
  
  assign(paste0("boot_cali_diff",j), (get(paste0("boot_cali", j ) ) - get(paste0("boot_caliPP", j) ) )  )  
  
}




# ##PLANN
lowerP3caliPP <- apply(boot_cali_diff3P, 2, quantile, 0.025)
upperP3caliPP <- apply(boot_cali_diff3P, 2, quantile, 0.975)

lowerP5caliPP <- apply(boot_cali_diff5P, 2, quantile, 0.025)
upperP5caliPP <- apply(boot_cali_diff5P, 2, quantile, 0.975)

lowerP10caliPP <- apply(boot_cali_diff10P, 2, quantile, 0.025)
upperP10caliPP <- apply(boot_cali_diff10P, 2, quantile, 0.975)

## F1

lowerF13caliPP <- apply(boot_cali_diff3F1, 2, quantile, 0.025)
upperF13caliPP <- apply(boot_cali_diff3F1, 2, quantile, 0.975)

lowerF15caliPP <- apply(boot_cali_diff5F1, 2, quantile, 0.025)
upperF15caliPP <- apply(boot_cali_diff5F1, 2, quantile, 0.975)

lowerF110caliPP <- apply(boot_cali_diff10F1, 2, quantile, 0.025)
upperF110caliPP <- apply(boot_cali_diff10F1, 2, quantile, 0.975)

## F2

lowerF23caliPP <- apply(boot_cali_diff3F2, 2, quantile, 0.025)
upperF23caliPP <- apply(boot_cali_diff3F2, 2, quantile, 0.975)

lowerF25caliPP <- apply(boot_cali_diff5F2, 2, quantile, 0.025)
upperF25caliPP <- apply(boot_cali_diff5F2, 2, quantile, 0.975)

lowerF210caliPP <- apply(boot_cali_diff10F2, 2, quantile, 0.025)
upperF210caliPP <- apply(boot_cali_diff10F2, 2, quantile, 0.975)


IC_caliPP <- list(P = list("3y" = list(low = lowerP3caliPP, upp = upperP3caliPP) ,"5y" = list(low = lowerP5caliPP, upp = upperP5caliPP), "10y" = list(low = lowerP10caliPP, upp = upperP10caliPP)),
                  F1 = list("3y" = list(low = lowerF13caliPP, upp = upperF13caliPP) ,"5y" = list(low = lowerF15caliPP, upp = upperF15caliPP), "10y" = list(low = lowerF110caliPP, upp = upperF110caliPP)),
                  F2 = list("3y" = list(low = lowerF23caliPP, upp = upperF23caliPP) ,"5y" = list(low = lowerF25caliPP, upp = upperF25caliPP), "10y" = list(low = lowerF210caliPP, upp = upperF210caliPP))
)

boot_cali_diffP3val <- apply(boot_cali_diff3P , mean, MARGIN = 2)
boot_cali_diffP5val <- apply(boot_cali_diff5P , mean, MARGIN = 2)
boot_cali_diffP10val <- apply(boot_cali_diff10P , mean, MARGIN = 2)

boot_cali_diffF13val <- apply(boot_cali_diff3F1 , mean, MARGIN = 2)
boot_cali_diffF15val <- apply(boot_cali_diff5F1 , mean, MARGIN = 2)
boot_cali_diffF110val <- apply(boot_cali_diff10F1 , mean, MARGIN = 2)

boot_cali_diffF23val <- apply(boot_cali_diff3F2 , mean, MARGIN = 2)
boot_cali_diffF25val <- apply(boot_cali_diff5F2 , mean, MARGIN = 2)
boot_cali_diffF210val <- apply(boot_cali_diff10F2 , mean, MARGIN = 2)

##xtable IC_cali

dfcaliPPP3 <- data.frame(boot_cali_diffP3val, IC_caliPP$P$`3y`$low, IC_caliPP$P$`3y`$upp)
colnames(dfcaliPPP3) <- c("estimate", "lower", "upper")
dfcaliPPP5 <-  data.frame(boot_cali_diffP5val, IC_caliPP$P$`5y`$low, IC_caliPP$P$`5y`$upp)
colnames(dfcaliPPP5) <- c("estimate", "lower", "upper")
dfcaliPPP10 <- data.frame(boot_cali_diffP10val, IC_caliPP$P$`10y`$low, IC_caliPP$P$`10y`$upp)
colnames(dfcaliPPP10) <- c("estimate", "lower", "upper")

dfcaliPPP <- cbind(dfcaliPPP3, dfcaliPPP5, dfcaliPPP10)
dfcaliPPP <- round(dfcaliPPP, digit = 3)
colnames(dfcaliPPP) <- c("3y estimate", "3y lower", "3y upper", "5y estimate", "5y lower", "5y upper", "10y estimate", "10y lower", "10y upper")

dfcaliPPF13 <- data.frame(boot_cali_diffF13val, IC_caliPP$F1$`3y`$low, IC_caliPP$F1$`3y`$upp)
colnames(dfcaliPPF13) <- c("estimate", "lower", "upper")
dfcaliPPF15 <- data.frame(boot_cali_diffF15val, IC_caliPP$F1$`5y`$low, IC_caliPP$F1$`5y`$upp)
colnames(dfcaliPPF15) <- c("estimate", "lower", "upper")
dfcaliPPF110 <- data.frame(boot_cali_diffF110val, IC_caliPP$F1$`10y`$low, IC_caliPP$F1$`10y`$upp)
colnames(dfcaliPPF110) <- c("estimate", "lower", "upper")

dfcaliPPF1 <- cbind(dfcaliPPF13, dfcaliPPF15, dfcaliPPF110)
dfcaliPPF1 <- round(dfcaliPPF1, digit = 3)
colnames(dfcaliPPF1) <- c("3y estimate", "3y lower", "3y upper", "5y estimate", "5y lower", "5y upper", "10y estimate", "10y lower", "10y upper")

dfcaliPPF23 <- data.frame(boot_cali_diffF23val, IC_caliPP$F2$`3y`$low, IC_caliPP$F2$`3y`$upp)
colnames(dfcaliPPF23) <- c("estimate", "lower", "upper")
dfcaliPPF25 <- data.frame(boot_cali_diffF25val, IC_caliPP$F2$`5y`$low, IC_caliPP$F2$`5y`$upp)
colnames(dfcaliPPF25) <- c("estimate", "lower", "upper")
dfcaliPPF210 <- data.frame(boot_cali_diffF210val, IC_caliPP$F2$`10y`$low, IC_caliPP$F2$`10y`$upp)
colnames(dfcaliPPF210) <- c("estimate", "lower", "upper")


dfcaliPPF2 <- cbind(dfcaliPPF23, dfcaliPPF25, dfcaliPPF210)
dfcaliPPF2 <- round(dfcaliPPF2, digit = 3)
colnames(dfcaliPPF2) <- c("3y estimate", "3y lower", "3y upper", "5y estimate", "5y lower", "5y upper", "10y estimate", "10y lower", "10y upper")

dfcaliPP <- rbind(dfcaliPPP, dfcaliPPF1, dfcaliPPF2)
dfcaliPP <- cbind(c(rep("PLANN", 5), rep("Spline PH", 5), rep("Spline NPH", 5)), dfcaliPP)
colnames(dfcaliPP)[1] <- "Method"

tmpP <- data.frame(
  Y3  = paste0(dfcaliPPP$`3y estimate`,  " [", dfcaliPPP$`3y lower`,  " : ", dfcaliPPP$`3y upper`,  "]"),
  Y5  = paste0(dfcaliPPP$`5y estimate`,  " [", dfcaliPPP$`5y lower`,  " : ", dfcaliPPP$`5y upper`,  "]"),
  Y10 = paste0(dfcaliPPP$`10y estimate`, " [", dfcaliPPP$`10y lower`, " : ", dfcaliPPP$`10y upper`, "]")
)

tmpF1 <- data.frame(
  Y3  = paste0(dfcaliPPF1$`3y estimate`,  " [", dfcaliPPF1$`3y lower`,  " : ", dfcaliPPF1$`3y upper`,  "]"),
  Y5  = paste0(dfcaliPPF1$`5y estimate`,  " [", dfcaliPPF1$`5y lower`,  " : ", dfcaliPPF1$`5y upper`,  "]"),
  Y10 = paste0(dfcaliPPF1$`10y estimate`, " [", dfcaliPPF1$`10y lower`, " : ", dfcaliPPF1$`10y upper`, "]")
)

tmpF2 <- data.frame(
  Y3  = paste0(dfcaliPPF2$`3y estimate`,  " [", dfcaliPPF2$`3y lower`,  " : ", dfcaliPPF2$`3y upper`,  "]"),
  Y5  = paste0(dfcaliPPF2$`5y estimate`,  " [", dfcaliPPF2$`5y lower`,  " : ", dfcaliPPF2$`5y upper`,  "]"),
  Y10 = paste0(dfcaliPPF2$`10y estimate`, " [", dfcaliPPF2$`10y lower`, " : ", dfcaliPPF2$`10y upper`, "]")
)

dfcaliPP_sim <- rbind(tmpP, tmpF1, tmpF2)
dfcaliPP_sim <- cbind(c(rep("PLANN", 5), rep("Spline PH", 5), rep("Spline NPH", 5)) , dfcaliPP_sim)
colnames(dfcaliPP_sim)[1] <- "Method"

################################################################

###########   PLOT

###############################################################
##plot mean
library(ggplot2)
library(patchwork)

#train

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

##################### plann

### récupération des valeurs de Pohar Perme
PPholdt <- c()

for(i in ind_m){
  
  a <- read.csv(paste0(path0, "pred/mean/PP/",i,"_PP_mean.csv"), sep = ";")
  PPholdt <- rbind(PPholdt,a)
}

PPmeant <- apply(PPholdt , mean, MARGIN = 2)
PPmeantALL <-PPmeant 
lowerPPtALL <- apply(PPholdt, 2, quantile, 0.025)
upperPPtALL <- apply(PPholdt, 2, quantile, 0.975)

PPmeant <- PPmeant[correstab]
lowerPPt <- lowerPPtALL[correstab]
upperPPt <- upperPPtALL[correstab]

for(j in strata_names){
  
  PPhold <- c()
  
  for(i in ind_m){
    a <- read.csv(paste0(path0, "pred_strata/mean/PP/",i,"_PP_mean_",j,".csv"), sep = ";")
    PPhold <- rbind(PPhold,a)
  }
  
  z <- apply(PPhold , mean, MARGIN = 2)
  assign(paste0("PPholdt_",j),z)
  z <- z[correstab]
  assign(paste0("PPmeant_",j), z)
  
  x <- apply(PPhold,  2, quantile, 0.025)
  y <- apply(PPhold,  2, quantile, 0.975)
  assign(paste0("lowerPPtALL_",j),x)
  x <- x[correstab]
  assign(paste0("lowerPPt_",j), x)
  assign(paste0("upperPPtALL_",j),y)
  y <- y[correstab]
  assign(paste0("upperPPt_",j), y)
}

# save(list = c("PPholdt", "PPholdt_HC", "PPholdt_HR", "PPholdt_FC", "PPholdt_FR","PPmeantALL", "PPmeant", "PPmeant_HC", "PPmeant_HR", "PPmeant_FC", "PPmeant_FR"), file = paste0("/home/thomas/Documents/Rstudio/Application/Application colrec/PPloaded.Rdata"))
load("~/Documents/Rstudio/Application/Application colrec/PPloaded.Rdata")

#valid

PPholdv <- c()

for(i in ind_m){
  
  a <- read.csv(paste0(path0, "pred/mean/PP/",i,"_PP_meanval.csv"), sep = ";")
  PPholdv <- rbind(PPholdv,a)
}

PPmeanv <- apply(PPholdv , mean, MARGIN = 2)
PPmeanvALL <-PPmeanv
lowerPPvALL <- apply(PPholdv, 2, quantile, 0.025)
upperPPvALL <- apply(PPholdv, 2, quantile, 0.975)

PPmeanv <- PPmeanv[correstab]
lowerPPv <- lowerPPvALL[correstab]
upperPPv <- upperPPvALL[correstab]

for(j in strata_names){
  
  PPhold <- c()
  
  for(i in ind_m){
    
    a <- read.csv(paste0(path0, "pred_strata/mean/PP/",i,"_PP_meanval_",j,".csv"), sep = ";")
    PPhold <- rbind(PPhold,a)
  }
  
  z <- apply(PPhold , mean, MARGIN = 2)
  assign(paste0("PPholdv_",j), z)
  z <- z[correstab]
  assign(paste0("PPmeanv_",j), z)
  
  x <- apply(PPhold,  2, quantile, 0.025)
  y <- apply(PPhold,  2, quantile, 0.975)
  assign(paste0("lowerPPvALL_",j),x)
  x <- x[correstab]
  assign(paste0("lowerPPv_",j), x)
  assign(paste0("upperPPvALL_",j),y)
  y <- y[correstab]
  assign(paste0("upperPPv_",j), y)
}

# save(list = c("PPholdv", "PPholdv_HR", "PPholdv_HC", "PPholdv_FR", "PPholdv_FC" ,"PPmeanv", "PPmeanv_HC", "PPmeanv_HR", "PPmeanv_FC", "PPmeanv_FR"), file = paste0("/home/thomas/Documents/Rstudio/Application/Application colrec/PPvalloaded.Rdata"))
load("~/Documents/Rstudio/Application/Application colrec/PPvalloaded.Rdata")
PPmeanvALL <- apply(PPholdv , mean, MARGIN = 2)

######

# save(list = ls(), file = paste0("/home/thomas/Documents/Rstudio/Application/Application colrec/plotsready.Rdata"))
load("/home/thomas/Documents/Rstudio/Application/Application colrec/plotsready.Rdata")


common_scale <- scale_colour_manual(name = "Method", breaks = c("PLANN", "Spline PH", "Spline NPH", "Pohar-Perme"),
                                    values = c("PLANN" = "#F8766D", "Spline PH" = "#00BA38", "Spline NPH" = "#619CFF", "Pohar-Perme" = "black"))
common_linetype <- scale_linetype_manual(name = "Method", values = c(
  "PLANN" = "solid","Spline PH" = "solid","Spline NPH" = "solid","Pohar-Perme" = "dashed"),breaks = c("PLANN", "Spline PH", "Spline NPH", "Pohar-Perme"))
### mean train
df_plot <- data.frame(time  = c(0,newtimespred/365.241), mean  = c(1,mean_boot_Pt),
                      lower = c(1,lowerPt), upper = c(1,upperPt), pp = c(1,PPmeantALL))

meanWPt <-  ggplot(df_plot, aes(x = time)) + 
  geom_line(aes(y = mean, colour = "PLANN", linetype = "PLANN"), linewidth = 1) +
  geom_line(aes(y = lower), colour = "#F8766D", linetype = "dashed") +
  geom_line(aes(y = upper), colour = "#F8766D", linetype = "dashed") +
  geom_ribbon(aes(ymin = lower, ymax = upper), fill = "#F8766D", alpha = 0.15) +
  geom_line(aes(y = pp, colour = "Pohar-Perme", linetype = "Pohar-Perme"), linewidth = 1) +
  # scale_colour_manual(name = "Method", values = c("PLANN" = "#F8766D","Pohar-Perme" = "black"))+
  common_scale + 
  common_linetype + 
  coord_cartesian(ylim = c(0, 1)) +
  labs(title = "PLANN", x = "Time (years)", y = "Net Survival") +
  theme_minimal(base_family = "CMU Serif",base_size = 14) +
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
        plot.title = element_text(hjust = 0.5),
        axis.ticks = element_line(color = "black", size = 0.5),  # tick line
        axis.ticks.length = unit(0.25, "cm"),
        legend.position =  "none", #c(0.05, 0.05),
        legend.justification = c("left", "bottom"),
        legend.background = element_rect(fill = "white", color = "black")
  ) + 
  scale_x_continuous(limits = c(0, 10.5), expand = expansion(mult = 0), labels = labels) +
  scale_y_continuous(limits = c(0, 1), expand = expansion(mult = 0), breaks = function(x) {
    b <- pretty(x)
    b[b != 0]           
  }, labels = labels) 

plot(meanWPt)
plot(meanWPt + theme(legend.position = c(0.05, 0.05) ))

ggsave(filename = paste0("/home/thomas/Documents/Articles/Article survivalPLANN/Figures v2/Figures 4 strates/MEAN_BOOT_P_W_T.png"),
       plot = meanWPt, width = 8, height = 6, dpi = 800)

ggsave(filename = paste0("/home/thomas/Documents/Articles/Article survivalPLANN/Figures v2/Figures 4 strates/MEAN_BOOT_P_W_T_leg.png"),
       plot = meanWPt + theme(legend.position = c(0.05, 0.05) ), width = 8, height = 6, dpi = 800)


##################### spline ph

df_plot <- data.frame(time  = c(0,newtimespred/365.241), mean  = c(1,mean_boot_F1t),
                      lower = c(1,lowerF1t), upper = c(1,upperF1t), pp = c(1,PPmeantALL))

meanWF1t <-  ggplot(df_plot, aes(x = time)) + 
  geom_line(aes(y = mean, colour = "Spline PH", linetype = "Spline PH"), linewidth = 1) +
  geom_line(aes(y = lower), colour = "#00BA38", linetype = "dashed") +
  geom_line(aes(y = upper), colour = "#00BA38", linetype = "dashed") +
  geom_ribbon(aes(ymin = lower, ymax = upper), fill = "#00BA38", alpha = 0.15) +
  geom_line(aes(y = pp, colour = "Pohar-Perme", linetype = "Pohar-Perme"), linewidth = 1) +
  # scale_colour_manual(name = "Method", values = c("Spline PH" = "#00BA38","Pohar-Perme" = "black"))+
  common_scale + 
  common_linetype + 
  coord_cartesian(ylim = c(0, 1)) +
  labs(title = "Spline PH", x = "Time (years)", y = "Net Survival"
  ) + theme_minimal(base_family = "CMU Serif",base_size = 14) +
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
        plot.title = element_text(hjust = 0.5),
        axis.ticks = element_line(color = "black", size = 0.5),  # tick line
        axis.ticks.length = unit(0.25, "cm"),
        legend.position =  "none", #c(0.98, 0.10)
        legend.justification = c("left", "bottom"),
        legend.background = element_rect(fill = "white", color = "black")
  ) + 
  scale_x_continuous(limits = c(0, 10.5), expand = expansion(mult = 0), labels = labels) +
  scale_y_continuous(limits = c(0, 1), expand = expansion(mult = 0), breaks = function(x) {
    b <- pretty(x)
    b[b != 0]           
  }, labels = labels) 


ggsave(filename = paste0("/home/thomas/Documents/Articles/Article survivalPLANN/Figures v2/Figures 4 strates/MEAN_BOOT_F1_W_T.png"),
       plot = meanWF1t, width = 8, height = 6, dpi = 800)

ggsave(filename = paste0("/home/thomas/Documents/Articles/Article survivalPLANN/Figures v2/Figures 4 strates/MEAN_BOOT_F1_W_T_leg.png"),
       plot = meanWF1t + theme(legend.position = c(0.05, 0.05) ), width = 8, height = 6, dpi = 800)


##################### spline nph

df_plot <- data.frame(time  = c(0,newtimespred/365.241), mean  = c(1,mean_boot_F2t),
                      lower = c(1,lowerF2t), upper = c(1,upperF2t), pp = c(1,PPmeantALL))

meanWF2t <-  ggplot(df_plot, aes(x = time)) + 
  geom_line(aes(y = mean, colour = "Spline NPH", linetype = "Spline NPH"), linewidth = 1) +
  geom_line(aes(y = lower), colour = "#619CFF", linetype = "dashed") +
  geom_line(aes(y = upper), colour = "#619CFF", linetype = "dashed") +
  geom_ribbon(aes(ymin = lower, ymax = upper), fill = "#619CFF", alpha = 0.15) +
  geom_line(aes(y = pp, colour = "Pohar-Perme", linetype = "Pohar-Perme"), linewidth = 1) +
  # scale_colour_manual(name = "Method",values = c("Spline NPH" = "#619CFF","Pohar-Perme" = "black"))+
  common_scale + 
  common_linetype + 
  coord_cartesian(ylim = c(0, 1)) +
  labs(title = "Spline NPH", x = "Time (years)", y = "Net Survival"
  ) + theme_minimal(base_family = "CMU Serif",base_size = 14) +
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
        plot.title = element_text(hjust = 0.5),
        axis.ticks = element_line(color = "black", size = 0.5),  # tick line
        axis.ticks.length = unit(0.25, "cm"),
        legend.position =  "none", #c(0.98, 0.10)
        legend.justification = c("left", "bottom"),
        legend.background = element_rect(fill = "white", color = "black")
  ) + 
  scale_x_continuous(limits = c(0, 10.5), expand = expansion(mult = 0), labels = labels) +
  scale_y_continuous(limits = c(0, 1), expand = expansion(mult = 0), breaks = function(x) {
    b <- pretty(x)
    b[b != 0]           
  }, labels = labels) 

ggsave(filename = paste0("/home/thomas/Documents/Articles/Article survivalPLANN/Figures v2/Figures 4 strates/MEAN_BOOT_F2_W_T.png"),
       plot = meanWF2t, width = 8, height = 6, dpi = 800)

ggsave(filename = paste0("/home/thomas/Documents/Articles/Article survivalPLANN/Figures v2/Figures 4 strates/MEAN_BOOT_F2_W_T_leg.png"),
       plot = meanWF2t + theme(legend.position = c(0.05, 0.05) ), width = 8, height = 6, dpi = 800)

### agrégation
mean_plot_t_agg <- (meanWPt + meanWF1t + meanWF2t) + plot_layout(nrow = 1, guides = "collect") 

ggsave(filename = paste0("/home/thomas/Documents/Articles/Article survivalPLANN/Figures v2/Figures 4 strates/MEAN_BOOT_ALL_W_T.png"),
       plot = mean_plot_t_agg, width = 12, height = 4, dpi = 800)

a = meanWPt + theme(legend.position = "bottom")
b = meanWF1t + theme(legend.position = "bottom")
c = meanWF2t + theme(legend.position = "bottom")

mean_plot_t_agg_leg <- (a+b+c)

ggsave(filename = paste0("/home/thomas/Documents/Articles/Article survivalPLANN/Figures v2/Figures 4 strates/MEAN_BOOT_ALL_W_T_leg.png"),
       plot = mean_plot_t_agg_leg, width = 12, height = 4, dpi = 800)

#############################  validation


##################### plann

df_plot <- data.frame(time  = c(0,newtimespred/365.241), mean  = c(1,mean_boot_Pval),
                      lower = c(1,lowerP), upper = c(1,upperP), pp = c(1,PPmeanvALL))

meanWPv <-  ggplot(df_plot, aes(x = time)) + 
  geom_line(aes(y = mean, colour = "PLANN", linetype = "PLANN"), linewidth = 1) +
  geom_line(aes(y = lower), colour = "#F8766D", linetype = "dashed") +
  geom_line(aes(y = upper), colour = "#F8766D", linetype = "dashed") +
  geom_ribbon(aes(ymin = lower, ymax = upper), fill = "#F8766D", alpha = 0.15) +
  geom_line(aes(y = pp, colour = "Pohar-Perme", linetype = "Pohar-Perme"), linewidth = 1) +
  # scale_colour_manual(values = c("PLANN" = "#F8766D","Pohar-Perme" = "black"))+
  common_scale + 
  common_linetype + 
  coord_cartesian(ylim = c(0, 1)) +
  labs(title = "PLANN", x = "Time (years)", y = "Net Survival"
  ) + theme_minimal(base_family = "CMU Serif",base_size = 14) +
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
        plot.title = element_text(hjust = 0.5),
        axis.ticks = element_line(color = "black", size = 0.5),  # tick line
        axis.ticks.length = unit(0.25, "cm"),
        legend.position =  "none", #c(0.98, 0.10)
        legend.justification = c("left", "bottom"),
        legend.background = element_rect(fill = "white", color = "black")
  ) + 
  scale_x_continuous(limits = c(0, 10.5), expand = expansion(mult = 0), labels = labels) +
  scale_y_continuous(limits = c(0, 1), expand = expansion(mult = 0), breaks = function(x) {
    b <- pretty(x)
    b[b != 0]           
  }, labels = labels) 

##################### spline ph

df_plot <- data.frame(time  = c(0,newtimespred/365.241), mean  = c(1,mean_boot_F1val),
                      lower = c(1,lowerF1), upper = c(1,upperF1), pp = c(1,PPmeanvALL))

meanWF1v <-  ggplot(df_plot, aes(x = time)) + 
  geom_line(aes(y = mean, colour = "Spline PH", linetype = "Spline PH"), linewidth = 1) +
  geom_line(aes(y = lower), colour = "#00BA38", linetype = "dashed") +
  geom_line(aes(y = upper), colour = "#00BA38", linetype = "dashed") +
  geom_ribbon(aes(ymin = lower, ymax = upper), fill = "#00BA38", alpha = 0.15) +
  geom_line(aes(y = pp, colour = "Pohar-Perme", linetype = "Pohar-Perme"), linewidth = 1) +
  # scale_colour_manual(values = c("Spline PH" = "#00BA38","Pohar-Perme" = "black"))+
  common_scale +  
  common_linetype + 
  coord_cartesian(ylim = c(0, 1)) +
  labs(title = "Spline PH", x = "Time (years)", y = "Net Survival"
  ) + theme_minimal(base_family = "CMU Serif",base_size = 14) +
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
        plot.title = element_text(hjust = 0.5),
        axis.ticks = element_line(color = "black", size = 0.5),  # tick line
        axis.ticks.length = unit(0.25, "cm"),
        legend.position =  "none", #c(0.98, 0.10)
        legend.justification = c("left", "bottom"),
        legend.background = element_rect(fill = "white", color = "black")
  ) + 
  scale_x_continuous(limits = c(0, 10.5), expand = expansion(mult = 0), labels = labels) +
  scale_y_continuous(limits = c(0, 1), expand = expansion(mult = 0), breaks = function(x) {
    b <- pretty(x)
    b[b != 0]           
  }, labels = labels) 

##################### spline nph
df_plot <- data.frame(time  = c(0,newtimespred/365.241), mean  = c(1,mean_boot_F2val),
                      lower = c(1,lowerF2), upper = c(1,upperF2), pp = c(1,PPmeanvALL))

meanWF2v <-  ggplot(df_plot, aes(x = time)) + 
  geom_line(aes(y = mean, colour = "Spline NPH", linetype = "Spline NPH"), linewidth = 1) +
  geom_line(aes(y = lower), colour = "#619CFF", linetype = "dashed") +
  geom_line(aes(y = upper), colour = "#619CFF", linetype = "dashed") +
  geom_ribbon(aes(ymin = lower, ymax = upper), fill = "#619CFF", alpha = 0.15) +
  geom_line(aes(y = pp, colour = "Pohar-Perme", linetype = "Pohar-Perme"), linewidth = 1) +
  # scale_colour_manual(values = c("Spline NPH" = "#619CFF","Pohar-Perme" = "black"))+
  common_scale +  
  common_linetype + 
  coord_cartesian(ylim = c(0, 1)) +
  labs(title = "Spline NPH", x = "Time (years)", y = "Net Survival"
  ) + theme_minimal(base_family = "CMU Serif",base_size = 14) +
  theme(panel.border = element_blank(), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
        plot.title = element_text(hjust = 0.5),
        axis.ticks = element_line(color = "black", size = 0.5),  # tick line
        axis.ticks.length = unit(0.25, "cm"),
        legend.position =  "none", #c(0.98, 0.10)
        legend.justification = c("left", "bottom"),
        legend.background = element_rect(fill = "white", color = "black")
  ) + 
  scale_x_continuous(limits = c(0, 10.5), expand = expansion(mult = 0), labels = labels) +
  scale_y_continuous(limits = c(0, 1), expand = expansion(mult = 0), breaks = function(x) {
    b <- pretty(x)
    b[b != 0]           
  }, labels = labels) 


#plann save
ggsave(filename = paste0("/home/thomas/Documents/Articles/Article survivalPLANN/Figures v2/Figures 4 strates/MEAN_BOOT_P_W_V.png"),
       plot = meanWPv, width = 8, height = 6, dpi = 800)

ggsave(filename = paste0("/home/thomas/Documents/Articles/Article survivalPLANN/Figures v2/Figures 4 strates/MEAN_BOOT_P_W_V_leg.png"),
       plot = meanWPv + theme(legend.position = c(0.05, 0.05) ), width = 8, height = 6, dpi = 800)

#flex1 save
ggsave(filename = paste0("/home/thomas/Documents/Articles/Article survivalPLANN/Figures v2/Figures 4 strates/MEAN_BOOT_F1_W_V.png"),
       plot = meanWF1v, width = 8, height = 6, dpi = 800)

ggsave(filename = paste0("/home/thomas/Documents/Articles/Article survivalPLANN/Figures v2/Figures 4 strates/MEAN_BOOT_F1_W_V_leg.png"),
       plot = meanWF1v + theme(legend.position = c(0.05, 0.05) ), width = 8, height = 6, dpi = 800)

#flex2 save
ggsave(filename = paste0("/home/thomas/Documents/Articles/Article survivalPLANN/Figures v2/Figures 4 strates/MEAN_BOOT_F2_W_V.png"),
       plot = meanWF2v, width = 8, height = 6, dpi = 800)

ggsave(filename = paste0("/home/thomas/Documents/Articles/Article survivalPLANN/Figures v2/Figures 4 strates/MEAN_BOOT_F2_W_V_leg.png"),
       plot = meanWF2v + theme(legend.position = c(0.05, 0.05) ), width = 8, height = 6, dpi = 800)


###ALL
mean_plot_v_agg <- (meanWPv + meanWF1v + meanWF2v) + plot_layout(nrow = 1, guides = "collect") 

ggsave(filename = paste0("/home/thomas/Documents/Articles/Article survivalPLANN/Figures v2/Figures 4 strates/MEAN_BOOT_ALL_W_V.png"),
       plot = mean_plot_v_agg, width = 12, height = 4, dpi = 800)

a = meanWPv + theme(legend.position = "bottom")
b = meanWF1v + theme(legend.position = "bottom")
c = meanWF2v + theme(legend.position = "bottom")

mean_plot_v_agg_leg <- (a+b+c)

ggsave(filename = paste0("/home/thomas/Documents/Articles/Article survivalPLANN/Figures v2/Figures 4 strates/MEAN_BOOT_ALL_W_V_leg.png"),
       plot = mean_plot_v_agg_leg, width = 12, height = 4, dpi = 800)


###################################" PAR STRATES

for(j in strata_names){
  
  for(q in c("t", "v")){
    
    if(q == "t"){
      #PLANN
      mean_bootPj <- get(paste0("mean_boot_P", q, "_",j))
      lowerPj <- get(paste0("lowerP",q,"_",j))
      upperPj <- get(paste0("upperP",q,"_",j))
      PPmeanALLj <- get(paste0("PPholdt_",j))
      #f1
      mean_bootF1j <- get(paste0("mean_boot_F1", q, "_",j))
      lowerF1j <- get(paste0("lowerF1",q,"_",j))
      upperF1j <- get(paste0("upperF1",q,"_",j))
      #F2
      mean_bootF2j <- get(paste0("mean_boot_F2", q, "_",j))
      lowerF2j <- get(paste0("lowerF2",q,"_",j))
      upperF2j <- get(paste0("upperF2",q,"_",j))
    }else{
      #plann
      mean_bootPj <- get(paste0("mean_boot_Pval_",j))
      lowerPj <- get(paste0("lowerP_",j))
      upperPj <- get(paste0("upperP_",j))
      PPmeanALLj <- get(paste0("PPholdv_",j))
      #f1
      mean_bootF1j <- get(paste0("mean_boot_F1val_",j))
      lowerF1j <- get(paste0("lowerF1_",j))
      upperF1j <- get(paste0("upperF1_",j))  
      #F2
      mean_bootF2j <- get(paste0("mean_boot_F2val_",j))
      lowerF2j <- get(paste0("lowerF2_",j))
      upperF2j <- get(paste0("upperF2_",j))  
    }
    
    df_plot <- data.frame(time  = c(0,newtimespred/365.241), mean  = c(1,mean_bootPj),
                          lower = c(1,lowerPj), upper = c(1,upperPj), pp = c(1,PPmeanALLj))
    
    meanWPt <-  ggplot(df_plot, aes(x = time)) + 
      geom_line(aes(y = mean, colour = "PLANN", linetype = "PLANN"), linewidth = 1) +
      geom_line(aes(y = lower), colour = "#F8766D", linetype = "dashed") +
      geom_line(aes(y = upper), colour = "#F8766D", linetype = "dashed") +
      geom_ribbon(aes(ymin = lower, ymax = upper), fill = "#F8766D", alpha = 0.15) +
      geom_line(aes(y = pp, colour = "Pohar-Perme", linetype = "Pohar-Perme"), linewidth = 1) +
      # scale_colour_manual(name = "Method", values = c("PLANN" = "#F8766D","Pohar-Perme" = "black"))+
      common_scale +  
      common_linetype + 
      coord_cartesian(ylim = c(0, 1)) +
      labs(title = "PLANN", x = "Time (years)", y = "Net Survival") +
      theme_minimal(base_family = "CMU Serif",base_size = 14) +
      theme(panel.border = element_blank(), panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
            plot.title = element_text(hjust = 0.5),
            axis.ticks = element_line(color = "black", size = 0.5),  # tick line
            axis.ticks.length = unit(0.25, "cm"),
            legend.position =  "none", #c(0.05, 0.05),
            legend.justification = c("left", "bottom"),
            legend.background = element_rect(fill = "white", color = "black")
      ) + 
      scale_x_continuous(limits = c(0, 10.5), expand = expansion(mult = 0), labels = labels) +
      scale_y_continuous(limits = c(0, 1), expand = expansion(mult = 0), breaks = function(x) {
        b <- pretty(x)
        b[b != 0]           
      }, labels = labels) 
    
    plot(meanWPt)
    plot(meanWPt + theme(legend.position = c(0.05, 0.05) ))
    
    ggsave(filename = paste0("/home/thomas/Documents/Articles/Article survivalPLANN/Figures v2/Figures 4 strates/Strates/MEAN_BOOT_P_",j,"_",q,".png"),
           plot = meanWPt, width = 8, height = 6, dpi = 800)
    
    ggsave(filename = paste0("/home/thomas/Documents/Articles/Article survivalPLANN/Figures v2/Figures 4 strates/Strates/MEAN_BOOT_P_",j,"_",q,"_leg.png"),
           plot = meanWPt + theme(legend.position = c(0.05, 0.05) ), width = 8, height = 6, dpi = 800)
    
    
    ##################### spline ph
    df_plot <- data.frame(time  = c(0,newtimespred/365.241), mean  = c(1,mean_bootF1j),
                          lower = c(1,lowerF1j), upper = c(1,upperF1j), pp = c(1,PPmeanALLj))
    
    meanWF1t <-  ggplot(df_plot, aes(x = time)) + 
      geom_line(aes(y = mean, colour = "Spline PH", linetype = "Spline PH"), linewidth = 1) +
      geom_line(aes(y = lower), colour = "#00BA38", linetype = "dashed") +
      geom_line(aes(y = upper), colour = "#00BA38", linetype = "dashed") +
      geom_ribbon(aes(ymin = lower, ymax = upper), fill = "#00BA38", alpha = 0.15) +
      geom_line(aes(y = pp, colour = "Pohar-Perme", linetype = "Pohar-Perme"), linewidth = 1) +
      # scale_colour_manual(name = "Method", values = c("Spline PH" = "#00BA38","Pohar-Perme" = "black"))+
      common_scale +  
      common_linetype + 
      coord_cartesian(ylim = c(0, 1)) +
      labs(title = "Spline PH", x = "Time (years)", y = "Net Survival"
      ) + theme_minimal(base_family = "CMU Serif",base_size = 14) +
      theme(panel.border = element_blank(), panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
            plot.title = element_text(hjust = 0.5),
            axis.ticks = element_line(color = "black", size = 0.5),  # tick line
            axis.ticks.length = unit(0.25, "cm"),
            legend.position =  "none", #c(0.98, 0.10)
            legend.justification = c("left", "bottom"),
            legend.background = element_rect(fill = "white", color = "black")
      ) + 
      scale_x_continuous(limits = c(0, 10.5), expand = expansion(mult = 0), labels = labels) +
      scale_y_continuous(limits = c(0, 1), expand = expansion(mult = 0), breaks = function(x) {
        b <- pretty(x)
        b[b != 0]           
      }, labels = labels) 
    
    
    ggsave(filename = paste0("/home/thomas/Documents/Articles/Article survivalPLANN/Figures v2/Figures 4 strates/Strates/MEAN_BOOT_F1_",j,"_",q,".png"),
           plot = meanWF1t, width = 8, height = 6, dpi = 800)
    
    ggsave(filename = paste0("/home/thomas/Documents/Articles/Article survivalPLANN/Figures v2/Figures 4 strates/Strates/MEAN_BOOT_F1_",j,"_",q,"_leg.png"),
           plot = meanWF1t + theme(legend.position = c(0.05, 0.05) ), width = 8, height = 6, dpi = 800)
    
    
    ##################### spline nph
    
    df_plot <- data.frame(time  = c(0,newtimespred/365.241), mean  = c(1,mean_bootF2j),
                          lower = c(1,lowerF2j), upper = c(1,upperF2j), pp = c(1,PPmeanALLj))
    
    meanWF2t <-  ggplot(df_plot, aes(x = time)) + 
      geom_line(aes(y = mean, colour = "Spline NPH", linetype = "Spline NPH"), linewidth = 1) +
      geom_line(aes(y = lower), colour = "#619CFF", linetype = "dashed") +
      geom_line(aes(y = upper), colour = "#619CFF", linetype = "dashed") +
      geom_ribbon(aes(ymin = lower, ymax = upper), fill = "#619CFF", alpha = 0.15) +
      geom_line(aes(y = pp, colour = "Pohar-Perme", linetype = "Pohar-Perme"), linewidth = 1) +
      # scale_colour_manual(name = "Method",values = c("Spline NPH" = "#619CFF","Pohar-Perme" = "black"))+
      common_scale +  
      common_linetype + 
      coord_cartesian(ylim = c(0, 1)) +
      labs(title = "Spline NPH", x = "Time (years)", y = "Net Survival"
      ) + theme_minimal(base_family = "CMU Serif",base_size = 14) +
      theme(panel.border = element_blank(), panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
            plot.title = element_text(hjust = 0.5),
            axis.ticks = element_line(color = "black", size = 0.5),  # tick line
            axis.ticks.length = unit(0.25, "cm"),
            legend.position =  "none", #c(0.98, 0.10)
            legend.justification = c("left", "bottom"),
            legend.background = element_rect(fill = "white", color = "black")
      ) + 
      scale_x_continuous(limits = c(0, 10.5), expand = expansion(mult = 0), labels = labels) +
      scale_y_continuous(limits = c(0, 1), expand = expansion(mult = 0), breaks = function(x) {
        b <- pretty(x)
        b[b != 0]           
      }, labels = labels) 
    
    ggsave(filename = paste0("/home/thomas/Documents/Articles/Article survivalPLANN/Figures v2/Figures 4 strates/Strates/MEAN_BOOT_F2_",j,"_",q,".png"),
           plot = meanWF2t, width = 8, height = 6, dpi = 800)
    
    ggsave(filename = paste0("/home/thomas/Documents/Articles/Article survivalPLANN/Figures v2/Figures 4 strates/Strates/MEAN_BOOT_F2_",j,"_",q,"_leg.png"),
           plot = meanWF2t + theme(legend.position = c(0.05, 0.05) ), width = 8, height = 6, dpi = 800)
    
    ### agrégation
    mean_plot_t_agg <- (meanWPt + meanWF1t + meanWF2t) + plot_layout(nrow = 1, guides = "collect") 
    
    ggsave(filename = paste0("/home/thomas/Documents/Articles/Article survivalPLANN/Figures v2/Figures 4 strates/Strates/MEAN_BOOT_ALL_",j,"_",q,".png"),
           plot = mean_plot_t_agg, width = 12, height = 4, dpi = 800)
    
    a = meanWPt + theme(legend.position = "bottom")
    b = meanWF1t + theme(legend.position = "bottom")
    c = meanWF2t + theme(legend.position = "bottom")
    
    mean_plot_t_agg_leg <- (a+b+c)
    
    ggsave(filename = paste0("/home/thomas/Documents/Articles/Article survivalPLANN/Figures v2/Figures 4 strates/Strates/MEAN_BOOT_ALL_",j,"_",q,"_leg.png"),
           plot = mean_plot_t_agg_leg, width = 12, height = 4, dpi = 800)
    
    
  }
  
}


#################################################################################
#forest plot pour les moyennes 

###################### TRAIN ##############################################


dfmeanPt <- data.frame(mean_boot_Pt[correstab], lowerPt[correstab], upperPt[correstab])
colnames(dfmeanPt) <- c("estimate", "lower", "upper")
dfmeanPt <- round(dfmeanPt, digit = 3)

dfmeanF1t <- data.frame(mean_boot_F1t[correstab], lowerF1t[correstab], upperF1t[correstab])
colnames(dfmeanF1t) <- c("estimate", "lower", "upper")
dfmeanF1t <- round(dfmeanF1t, digit = 3)

dfmeanF2t <- data.frame(mean_boot_F2t[correstab], lowerF2t[correstab], upperF2t[correstab])
colnames(dfmeanF2t) <- c("estimate", "lower", "upper")
dfmeanF2t <- round(dfmeanF2t, digit = 3)

dfmeanPPt <- data.frame(PPmeantALL[correstab], lowerPPtALL[correstab], upperPPtALL[correstab])
colnames(dfmeanPPt) <- c("estimate", "lower", "upper")
dfmeanPPt <- round(dfmeanPPt, digit = 3)


dfmeant <- rbind(dfmeanPt, dfmeanF1t, dfmeanF2t, dfmeanPPt)
dfmeant <- cbind(c(rep("PLANN", 3), rep("Spline PH", 3), rep("Spline NPH", 3), rep("PP",3)), dfmeant)
colnames(dfmeant)[1] <- "Method"

for(j in strata_names){
  
  Phold <- data.frame(get(paste0("mean_boot_Pt_",j))[correstab], get(paste0("lowerPt_",j))[correstab], get(paste0("upperPt_",j))[correstab])
  colnames(Phold) <- c("estimate", "lower", "upper")
  Phold <- round(Phold, digit = 3)
  assign(paste0("dfmeanPt_",j), Phold )
  
  F1hold <- data.frame(get(paste0("mean_boot_F1t_",j))[correstab], get(paste0("lowerF1t_",j))[correstab], get(paste0("upperF1t_",j))[correstab])
  colnames(F1hold) <- c("estimate", "lower", "upper")
  F1hold <- round(F1hold, digit = 3)
  assign(paste0("dfmeanF1t_",j), F1hold )
  
  F2hold <- data.frame(get(paste0("mean_boot_F2t_",j))[correstab], get(paste0("lowerF2t_",j))[correstab], get(paste0("upperF2t_",j))[correstab])
  colnames(F2hold) <- c("estimate", "lower", "upper")
  F2hold <- round(F2hold, digit = 3)
  assign(paste0("dfmeanF2t_",j), F2hold )
  
  PPhold <- data.frame(get(paste0("PPmeant_",j)), get(paste0("lowerPPt_",j)), get(paste0("upperPPt_",j)))
  colnames(PPhold) <- c("estimate", "lower", "upper")
  PPhold <- round(PPhold, digit = 3)
  assign(paste0("dfmeanPPt_",j), PPhold )
  
  dfthold <- rbind(Phold, F1hold, F2hold, PPhold)
  dfthold <- cbind(c(rep("PLANN", 3), rep("Spline PH", 3), rep("Spline NPH", 3), rep("PP", 3)), dfthold)
  colnames(dfthold)[1] <- "Method"
  
  assign(paste0("dfmeant_",j), dfthold )
  
}



pp_df <- rbind(data.frame(Group = "Whole", Time = c("3 years","5 years","10 years"), PP = as.numeric(PPmeant)),
               data.frame(Group = "HC", Time = c("3 years","5 years","10 years"), PP = as.numeric(PPmeant_HC)),
               data.frame(Group = "HR", Time = c("3 years","5 years","10 years"), PP = as.numeric(PPmeant_HR)),
               data.frame(Group = "FC", Time = c("3 years","5 years","10 years"), PP = as.numeric(PPmeant_FC)),
               data.frame(Group = "FR", Time = c("3 years","5 years","10 years"), PP = as.numeric(PPmeant_FR))
)

################# data frame pour le plotting
df_long <- rbind(dfmeant, dfmeant_HC, dfmeant_HR, dfmeant_FC, dfmeant_FR)
colnames(df_long)[2:4] <- c("Estimate", "Lower", "Upper")

df_long$Group <- c(rep("Whole",12), rep("HC",12), rep("HR",12), rep("FC",12), rep("FR",12))

df_long$Time <- rep(c("3 years", "5 years", "10 years"),20)


df_long$y_base <- as.numeric(factor(df_long$Group, levels = c("FR", "FC", "HR", "HC", "Whole")))

df_long$y_pos <- df_long$y_base +
  ifelse(df_long$Method == "PLANN",      0.30,
         ifelse(df_long$Method == "Spline PH",  0.10,
                ifelse(df_long$Method == "Spline NPH", -0.10,
                       -0.30)))

pp_df$y_baseP <- c(FR = 1, FC = 2, HR = 3, HC = 4, Whole = 5)[pp_df$Group]
pp_df$Time <- factor(
  pp_df$Time,
  levels = c("3 years", "5 years", "10 years")
)
df_long <- merge(df_long, pp_df,
                 by = c("Group", "Time"))

df_long$Time <- factor(
  df_long$Time,
  levels = c("3 years", "5 years", "10 years")
)
df_long$Method <- factor(
  df_long$Method,
  levels = c("PLANN", "Spline PH", "Spline NPH", "PP")
)

df_long$Group <- factor(
  df_long$Group,
  levels = c("Whole", "HC", "HR", "FC", "FR")
)

ForestMEANt <- ggplot(df_long, aes(x = Estimate, xmin = Lower, xmax = Upper, colour = Method)) +
  geom_errorbarh(aes(y = y_pos), width = 0.15) +
  geom_point(aes(y = y_pos,shape = Method), size = 3) +
  facet_wrap(~Time, nrow = 1) +
  scale_y_continuous(breaks = 1:5, labels = c("Female & Rectum", "Female & Colon", "Male & Rectum", "Male & Colon", "Whole"),
                     limits = c(0.5, 5.5),
                     expand = c(0, 0)) +
  theme_bw() +
  theme_minimal(base_family = "CMU Serif",base_size = 14) +
  labs(x = "Estimated survival",y = NULL)+
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_rect(fill = "white"),
    strip.background = element_rect(fill = "white"),
    # legend.position = "none",
    # panel.spacing.x = unit(0, "cm"),
    # panel.border = element_rect(colour = "black", fill = NA)
  )+ scale_shape_discrete(name = "Model") +
  scale_colour_manual(name = "Model",
                      values = c("PLANN" = "#F8766D","Spline PH" = "#00BA38",
                                 "Spline NPH" = "#619CFF","PP" = "black"),
                      breaks = c("PLANN","Spline PH","Spline NPH","PP"))+
  geom_hline(yintercept = c(1.5, 2.5, 3.5, 4.5), colour = "black",linewidth = 0.3)+
  geom_segment(
    data = pp_df,
    aes( x = PP, xend = PP,  y = y_baseP - 0.49, yend = y_baseP + 0.49
    ), inherit.aes = FALSE, colour = "black", linewidth = 0.35, linetype = "dashed"
  )

plot(ForestMEANt)

ggsave(filename = paste0("/home/thomas/Documents/Articles/Article survivalPLANN/Figures v2/Figures 4 strates/MEAN_BOOT_FOREST_t_leg.png"),
       plot = ForestMEANt, width = 8, height = 6, dpi = 800)


###################### VALIDATION ##############################################


dfmeanPv <- data.frame(mean_boot_Pval[correstab], lowerP[correstab], upperP[correstab])
colnames(dfmeanPv) <- c("estimate", "lower", "upper")
dfmeanPv <- round(dfmeanPv, digit = 3)

dfmeanF1v <- data.frame(mean_boot_F1val[correstab], lowerF1[correstab], upperF1[correstab])
colnames(dfmeanF1v) <- c("estimate", "lower", "upper")
dfmeanF1v <- round(dfmeanF1v, digit = 3)

dfmeanF2v <- data.frame(mean_boot_F2val[correstab], lowerF2[correstab], upperF2[correstab])
colnames(dfmeanF2v) <- c("estimate", "lower", "upper")
dfmeanF2v <- round(dfmeanF2v, digit = 3)

dfmeanPPv <- data.frame(PPmeanvALL[correstab], lowerPPvALL[correstab], upperPPvALL[correstab])
colnames(dfmeanPPv) <- c("estimate", "lower", "upper")
dfmeanPPv <- round(dfmeanPPv, digit = 3)

dfmeanv <- rbind(dfmeanPv, dfmeanF1v, dfmeanF2v, dfmeanPPv)
dfmeanv <- cbind(c(rep("PLANN", 3), rep("Spline PH", 3), rep("Spline NPH", 3), rep("PP",3)), dfmeanv)
colnames(dfmeanv)[1] <- "Method"

for(j in strata_names){
  
  Phold <- data.frame(get(paste0("mean_boot_Pval_",j))[correstab], get(paste0("lowerP_",j))[correstab], get(paste0("upperP_",j))[correstab])
  colnames(Phold) <- c("estimate", "lower", "upper")
  Phold <- round(Phold, digit = 3)
  assign(paste0("dfmeanPval_",j), Phold )
  
  F1hold <- data.frame(get(paste0("mean_boot_F1val_",j))[correstab], get(paste0("lowerF1_",j))[correstab], get(paste0("upperF1_",j))[correstab])
  colnames(F1hold) <- c("estimate", "lower", "upper")
  F1hold <- round(F1hold, digit = 3)
  assign(paste0("dfmeanF1val_",j), F1hold )
  
  F2hold <- data.frame(get(paste0("mean_boot_F2val_",j))[correstab], get(paste0("lowerF2_",j))[correstab], get(paste0("upperF2_",j))[correstab])
  colnames(F2hold) <- c("estimate", "lower", "upper")
  F2hold <- round(F2hold, digit = 3)
  assign(paste0("dfmeanF2val_",j), F2hold )
  
  PPhold <- data.frame(get(paste0("PPmeanv_",j)), get(paste0("lowerPPv_",j)), get(paste0("upperPPv_",j)))
  colnames(PPhold) <- c("estimate", "lower", "upper")
  PPhold <- round(PPhold, digit = 3)
  assign(paste0("dfmeanPPv_",j), PPhold )
  
  dfthold <- rbind(Phold, F1hold, F2hold, PPhold)
  dfthold <- cbind(c(rep("PLANN", 3), rep("Spline PH", 3), rep("Spline NPH", 3), rep("PP", 3)), dfthold)
  colnames(dfthold)[1] <- "Method"
  
  assign(paste0("dfmeanv_",j), dfthold )
  
}

### récupération des valeurs de Pohar Perme

ppv_df <- rbind(data.frame(Group = "Whole", Time = c("3 years","5 years","10 years"), PP = as.numeric(PPmeanv)),
                data.frame(Group = "HC", Time = c("3 years","5 years","10 years"), PP = as.numeric(PPmeanv_HC)),
                data.frame(Group = "HR", Time = c("3 years","5 years","10 years"), PP = as.numeric(PPmeanv_HR)),
                data.frame(Group = "FC", Time = c("3 years","5 years","10 years"), PP = as.numeric(PPmeanv_FC)),
                data.frame(Group = "FR", Time = c("3 years","5 years","10 years"), PP = as.numeric(PPmeanv_FR))
)

################# data frame pour le plotting
df_long <- rbind(dfmeanv, dfmeanv_HC, dfmeanv_HR, dfmeanv_FC, dfmeanv_FR)
colnames(df_long)[2:4] <- c("Estimate", "Lower", "Upper")

df_long$Group <- c(rep("Whole",12), rep("HC",12), rep("HR",12), rep("FC",12), rep("FR",12))

df_long$Time <- rep(c("3 years", "5 years", "10 years"),20)


df_long$y_base <- as.numeric(factor(df_long$Group, levels = c("FR", "FC", "HR", "HC", "Whole")))

df_long$y_pos <- df_long$y_base +
  ifelse(df_long$Method == "PLANN",      0.30,
         ifelse(df_long$Method == "Spline PH",  0.10,
                ifelse(df_long$Method == "Spline NPH", -0.10,
                       -0.30)))

ppv_df$y_baseP <- c(FR = 1, FC = 2, HR = 3, HC = 4, Whole = 5)[ppv_df$Group]
ppv_df$Time <- factor(
  ppv_df$Time,
  levels = c("3 years", "5 years", "10 years")
)
df_long <- merge(df_long, ppv_df,
                 by = c("Group", "Time"))

df_long$Time <- factor(
  df_long$Time,
  levels = c("3 years", "5 years", "10 years")
)
df_long$Method <- factor(
  df_long$Method,
  levels = c("PLANN", "Spline PH", "Spline NPH", "PP")
)

df_long$Group <- factor(
  df_long$Group,
  levels = c("Whole", "HC", "HR", "FC", "FR")
)


ForestMEANv <- ggplot(df_long, aes(x = Estimate, xmin = Lower, xmax = Upper, colour = Method)) +
  geom_errorbarh(aes(y = y_pos), width = 0.15) +
  geom_point(aes(y = y_pos,shape = Method), size = 3) +
  facet_wrap(~Time, nrow = 1) +
  scale_y_continuous(breaks = 1:5, labels = c("Female & Rectum", "Female & Colon", "Male & Rectum", "Male & Colon", "Whole"),
                     limits = c(0.5, 5.5),
                     expand = c(0, 0)) +
  theme_bw() +
  theme_minimal(base_family = "CMU Serif",base_size = 14) +
  labs(x = "Estimated survival",y = NULL)+
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_rect(fill = "white"),
    strip.background = element_rect(fill = "white"),
    # legend.position = "none",
    # panel.spacing.x = unit(0, "cm"),
    # panel.border = element_rect(colour = "black", fill = NA)
  )+ scale_shape_discrete(name = "Model") +
  scale_colour_manual(name = "Model",
                      values = c("PLANN" = "#F8766D","Spline PH" = "#00BA38",
                                 "Spline NPH" = "#619CFF","PP" = "black"),
                      breaks = c("PLANN","Spline PH","Spline NPH","PP"))+
  geom_hline(yintercept = c(1.5, 2.5, 3.5, 4.5), colour = "black",linewidth = 0.3)+
  geom_segment(
    data = ppv_df,
    aes( x = PP, xend = PP,  y = y_baseP - 0.49, yend = y_baseP + 0.49
    ), inherit.aes = FALSE, colour = "black", linewidth = 0.35, linetype = "dashed"
  )


ggsave(filename = paste0("/home/thomas/Documents/Articles/Article survivalPLANN/Figures v2/Figures 4 strates/MEAN_BOOT_FOREST_v_leg.png"),
       plot = ForestMEANv, width = 8, height = 6, dpi = 800)



#################################################################################
#### plot AUC
# 

##################### PLOT des AUC SEUELEMNT SUR LA BASE DE TRAIN (forest plot sur les deux plus bas)

## JE NAI PAS SAVE LES VALEURS DE Se ET Sp nécessaire pour faire une moyenne et construire de là, je peux les re faire tourner si nécessaire mais sinon bof
for(p in c(1,2,3)){
  
  nametime <- c(3,5,10)
  main_name <- main_vec[p]
  df <- rbind(
    data.frame(FPR = 1 - get(paste0("AUC_whole_", newtimes[p] ))$P$table$sp,TPR = get(paste0("AUC_whole_", newtimes[p]))$P$table$se,
               Model = "PLANN", lower = get(paste0("lowerP",nametime[p],"AUC")), upper = get(paste0("upperP",nametime[p],"AUC"))),
    data.frame(FPR   = 1 - get(paste0("AUC_whole_", newtimes[p]))$F1$table$sp,TPR = get(paste0("AUC_whole_", newtimes[p]))$F1$table$se,
               Model = "Spline PH", lower = get(paste0("lowerF1",nametime[p],"AUC")), upper = get(paste0("upperF1",nametime[p],"AUC"))),
    data.frame(FPR = 1 - get(paste0("AUC_whole_", newtimes[p]))$F2$table$sp, TPR = get(paste0("AUC_whole_", newtimes[p]))$F2$table$se,
               Model = "Spline NPH", lower = get(paste0("lowerF2",nametime[p],"AUC")), upper = get(paste0("upperF2",nametime[p],"AUC")))
  )
  
  plot_to_print <- ggplot(df, aes(x = FPR, y = TPR, color = Model, linetype = Model)) +
    
    geom_line(size = 1.2) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
    scale_color_manual(limits = c("PLANN", "Spline PH", "Spline NPH"),
                       values = c("Spline PH" = "#00BA38", "Spline NPH" = "#619CFF", "PLANN" = "#F8766D"),labels = c(
                         PLANN = paste0("PLANN            (AUC : ", sprintf("%.4f", get(paste0("AUC_whole_", newtimes[p] ))$P$auc ), ")"),
                         'Spline PH'    = paste0("Spline PH          (AUC : ", sprintf("%.4f", get(paste0("AUC_whole_", newtimes[p] ))$F1$auc), ")"),
                         'Spline NPH'   = paste0("Spline NPH        (AUC : ", sprintf("%.4f", get(paste0("AUC_whole_", newtimes[p] ))$F2$auc), ")"))
    )+
    scale_linetype_manual(limits = c("PLANN", "Spline PH", "Spline NPH"), values=c("Spline PH" ="dotted", "Spline NPH" ="longdash", "PLANN" ="twodash"), labels = c(
      PLANN = paste0("PLANN            (AUC : ", sprintf("%.4f", get(paste0("AUC_whole_", newtimes[p] ))$P$auc), ")"),
      'Spline PH'    = paste0("Spline PH          (AUC : ", sprintf("%.4f", get(paste0("AUC_whole_", newtimes[p] ))$F1$auc), ")"),
      'Spline NPH'   = paste0("Spline NPH        (AUC : ", sprintf("%.4f", get(paste0("AUC_whole_", newtimes[p] ))$F2$auc), ")"))) +
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
  
  ggsave(filename = paste0("/home/thomas/Documents/Articles/Article survivalPLANN/Figures v2/Figures 4 strates/AUC_", nametime[p] ,"_W_v.png"),
         plot = plot_to_print, width = 8, height = 6, dpi = 800)
  
}

################### A LA PLACE
# FOREST PLOT POUR VALEUR AUC 

dfAUC_pl <- data.frame(rbind(c(boot_AUCP3val, lowerP3AUC, upperP3AUC, boot_AUCP5val, lowerP5AUC, upperP5AUC, boot_AUCP10val, lowerP10AUC, upperP10AUC),
                             c(boot_AUCF13val, lowerF13AUC, upperF13AUC, boot_AUCF15val, lowerF15AUC, upperF15AUC, boot_AUCF110val, lowerF110AUC, upperF110AUC),
                             c(boot_AUCF23val, lowerF23AUC, upperF23AUC, boot_AUCF25val, lowerF25AUC, upperF25AUC, boot_AUCF210val, lowerF210AUC, upperF210AUC)
))
dfAUC_pl <- cbind(c("PLANN", "Spline PH", "Spline NPH"), dfAUC_pl)
colnames(dfAUC_pl) <- c("Method","3y estimate", "3y lower", "3y upper","5y estimate" ,"5y lower","5y upper","10y estimate", "10y lower", "10y upper")

df_long <- bind_rows(
  
  dfAUC_pl %>%
    transmute(
      Method, Time = "3 years", Estimate = `3y estimate`, Lower = `3y lower`, Upper = `3y upper`),
  
  dfAUC_pl %>%
    transmute(
      Method,  Time = "5 years", Estimate = `5y estimate`, Lower = `5y lower`, Upper = `5y upper`),
  
  dfAUC_pl %>% transmute(
    Method, Time = "10 years", Estimate = `10y estimate`, Lower = `10y lower`, Upper = `10y upper`)
  
)

df_long$Time <- factor(
  df_long$Time,
  levels = c("3 years", "5 years", "10 years")
)
df_long$Method <- factor(
  df_long$Method,
  levels = c("PLANN", "Spline PH", "Spline NPH")
)
df_long$y <- 1
df_long$Group <- rep(1, dim(df_long)[1])

df_long$y_base <- as.numeric(factor(df_long$Group, levels = 1))

df_long$y_pos <- df_long$y_base +
  ifelse(df_long$Method == "PLANN", 0.2,
         ifelse(df_long$Method == "Spline PH", 0,
                -0.2))


forestAUCv <- ggplot(df_long, aes( x = Estimate,
                                   xmin = Lower, xmax = Upper, colour = Method)) +
  geom_errorbarh(aes(y = y_pos), width = 0.05) +
  geom_point(aes(y = y_pos, shape = Method), size = 3) +
  facet_wrap(~Time, nrow = 1) +
  scale_linetype_manual(limits = c("PLANN", "Spline PH", "Spline NPH"), values=c("Spline PH" ="dotted", "Spline NPH" ="longdash", "PLANN" ="twodash"))+
  labs(x = "AUC values", y = "") +
  theme_bw() +
  theme_minimal(base_family = "CMU Serif",base_size = 14) +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_rect(fill = "white"),
    strip.background = element_rect(fill = "white"),
    axis.ticks.y = element_blank(),
    axis.text.y = element_blank()
  )+
  scale_y_continuous(limits = c(0.7, 1.3), expand = c(0, 0))+
  scale_colour_discrete(name = "Model")+
  scale_shape_discrete(name = "Model")


ggsave(filename = paste0("/home/thomas/Documents/Articles/Article survivalPLANN/Figures v2/Figures 4 strates/AUC_BOOT_FOREST_v_leg.png"),
       plot = forestAUCv, width = 8, height = 3, dpi = 800)


############### FOREST PLOT CALIBRATION

#dfcaliPP_sim
#dfcaliPP

dfcaliPP$Group <- rep(1:5, times = 3)

df_long <- bind_rows(
  
  dfcaliPP %>%
    transmute(
      Method, Group, Time = "3 years", Estimate = `3y estimate`, Lower = `3y lower`, Upper = `3y upper`),
  
  dfcaliPP %>%
    transmute(
      Method, Group, Time = "5 years", Estimate = `5y estimate`, Lower = `5y lower`, Upper = `5y upper`),
  
  dfcaliPP %>% transmute(
    Method, Group, Time = "10 years", Estimate = `10y estimate`, Lower = `10y lower`, Upper = `10y upper`)
  
)

df_long$Time <- factor(
  df_long$Time,
  levels = c("3 years", "5 years", "10 years")
)
df_long$Method <- factor(
  df_long$Method,
  levels = c("PLANN", "Spline PH", "Spline NPH")
)
df_long$Group <- 6 - df_long$Group

df_long$y_base <- as.numeric(factor(df_long$Group, levels = 5:1))

df_long$y_pos <- df_long$y_base +
  ifelse(df_long$Method == "PLANN", 0.25,
         ifelse(df_long$Method == "Spline PH", 0,
                -0.25))
forestCALIv <- ggplot(df_long, aes( x = Estimate,
                                    xmin = Lower, xmax = Upper, colour = Method)) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  geom_errorbarh(aes(y = y_pos), width = 0.05) +
  geom_point(aes(y = y_pos, shape = Method), size = 3) +
  facet_wrap(~Time, nrow = 1) +
  scale_linetype_manual(limits = c("PLANN", "Spline PH", "Spline NPH"), values=c("Spline PH" ="dotted", "Spline NPH" ="longdash", "PLANN" ="twodash"))+
  scale_x_continuous(limits = c(-0.2, 0.2)) + 
  labs(x = "Model-based estimations minus PP estimations", y = "Risk group by using quantiles of indiviudal predictions") +
  theme_bw() +
  theme_minimal(base_family = "CMU Serif",base_size = 14) +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_rect(fill = "white"),
    strip.background = element_rect(fill = "white"),
    # legend.position = "none",
    # panel.spacing.x = unit(0, "cm"),
    # panel.border = element_rect(colour = "black", fill = NA)
  )+
  geom_hline(
    yintercept = c(1.5, 2.5, 3.5, 4.5),
    colour = "black",
    linewidth = 0.3
  )+scale_colour_discrete(name = "Model") +
  scale_shape_discrete(name = "Model")

ggsave(filename = paste0("/home/thomas/Documents/Articles/Article survivalPLANN/Figures v2/Figures 4 strates/CALI_BOOT_FOREST_v_leg.png"),
       plot = forestCALIv, width = 8, height = 6, dpi = 800)

##########
# code pour concatener les 12 plots de survie par sous strates 

for(j in strata_names){
  
  for(q in c("t", "v")){
    
    if(q == "t"){
      #PLANN
      mean_bootPj <- get(paste0("mean_boot_P", q, "_",j))
      lowerPj <- get(paste0("lowerP",q,"_",j))
      upperPj <- get(paste0("upperP",q,"_",j))
      PPmeanALLj <- get(paste0("PPholdt_",j))
      #f1
      mean_bootF1j <- get(paste0("mean_boot_F1", q, "_",j))
      lowerF1j <- get(paste0("lowerF1",q,"_",j))
      upperF1j <- get(paste0("upperF1",q,"_",j))
      #F2
      mean_bootF2j <- get(paste0("mean_boot_F2", q, "_",j))
      lowerF2j <- get(paste0("lowerF2",q,"_",j))
      upperF2j <- get(paste0("upperF2",q,"_",j))
    }else{
      #plann
      mean_bootPj <- get(paste0("mean_boot_Pval_",j))
      lowerPj <- get(paste0("lowerP_",j))
      upperPj <- get(paste0("upperP_",j))
      PPmeanALLj <- get(paste0("PPholdv_",j))
      #f1
      mean_bootF1j <- get(paste0("mean_boot_F1val_",j))
      lowerF1j <- get(paste0("lowerF1_",j))
      upperF1j <- get(paste0("upperF1_",j))  
      #F2
      mean_bootF2j <- get(paste0("mean_boot_F2val_",j))
      lowerF2j <- get(paste0("lowerF2_",j))
      upperF2j <- get(paste0("upperF2_",j))  
    }
    
    df_plot <- data.frame(time  = c(0,newtimespred/365.241), mean  = c(1,mean_bootPj),
                          lower = c(1,lowerPj), upper = c(1,upperPj), pp = c(1,PPmeanALLj))
    
    meanWPt <-  ggplot(df_plot, aes(x = time)) + 
      geom_line(aes(y = mean, colour = "PLANN", linetype = "PLANN"), linewidth = 1) +
      geom_line(aes(y = lower), colour = "#F8766D", linetype = "dashed") +
      geom_line(aes(y = upper), colour = "#F8766D", linetype = "dashed") +
      geom_ribbon(aes(ymin = lower, ymax = upper), fill = "#F8766D", alpha = 0.15) +
      geom_line(aes(y = pp, colour = "Pohar-Perme", linetype = "Pohar-Perme"), linewidth = 1) +
      # scale_colour_manual(name = "Method", values = c("PLANN" = "#F8766D","Pohar-Perme" = "black"))+
      common_scale +  
      common_linetype + 
      coord_cartesian(ylim = c(0, 1)) +
      labs(title = "PLANN", x = "Time (years)", y = "Net Survival") +
      theme_minimal(base_family = "CMU Serif",base_size = 14) +
      theme(panel.border = element_blank(), panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
            plot.title = element_text(hjust = 0.5),
            axis.ticks = element_line(color = "black", size = 0.5),  # tick line
            axis.ticks.length = unit(0.25, "cm"),
            legend.position =  "none", #c(0.05, 0.05),
            legend.justification = c("left", "bottom"),
            legend.background = element_rect(fill = "white", color = "black")
      ) + 
      scale_x_continuous(limits = c(0, 10.5), expand = expansion(mult = 0), labels = labels) +
      scale_y_continuous(limits = c(0, 1), expand = expansion(mult = 0), breaks = function(x) {
        b <- pretty(x)
        b[b != 0]           
      }, labels = labels) 
    
    
    
    ##################### spline ph
    df_plot <- data.frame(time  = c(0,newtimespred/365.241), mean  = c(1,mean_bootF1j),
                          lower = c(1,lowerF1j), upper = c(1,upperF1j), pp = c(1,PPmeanALLj))
    
    meanWF1t <-  ggplot(df_plot, aes(x = time)) + 
      geom_line(aes(y = mean, colour = "Spline PH", linetype = "Spline PH"), linewidth = 1) +
      geom_line(aes(y = lower), colour = "#00BA38", linetype = "dashed") +
      geom_line(aes(y = upper), colour = "#00BA38", linetype = "dashed") +
      geom_ribbon(aes(ymin = lower, ymax = upper), fill = "#00BA38", alpha = 0.15) +
      geom_line(aes(y = pp, colour = "Pohar-Perme", linetype = "Pohar-Perme"), linewidth = 1) +
      # scale_colour_manual(name = "Method", values = c("Spline PH" = "#00BA38","Pohar-Perme" = "black"))+
      common_scale +  
      common_linetype + 
      coord_cartesian(ylim = c(0, 1)) +
      labs(title = "Spline PH", x = "Time (years)", y = "Net Survival"
      ) + theme_minimal(base_family = "CMU Serif",base_size = 14) +
      theme(panel.border = element_blank(), panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
            plot.title = element_text(hjust = 0.5),
            axis.ticks = element_line(color = "black", size = 0.5),  # tick line
            axis.ticks.length = unit(0.25, "cm"),
            legend.position =  "none", #c(0.98, 0.10)
            legend.justification = c("left", "bottom"),
            legend.background = element_rect(fill = "white", color = "black")
      ) + 
      scale_x_continuous(limits = c(0, 10.5), expand = expansion(mult = 0), labels = labels) +
      scale_y_continuous(limits = c(0, 1), expand = expansion(mult = 0), breaks = function(x) {
        b <- pretty(x)
        b[b != 0]           
      }, labels = labels) 
    
    ##################### spline nph
    
    df_plot <- data.frame(time  = c(0,newtimespred/365.241), mean  = c(1,mean_bootF2j),
                          lower = c(1,lowerF2j), upper = c(1,upperF2j), pp = c(1,PPmeanALLj))
    
    meanWF2t <-  ggplot(df_plot, aes(x = time)) + 
      geom_line(aes(y = mean, colour = "Spline NPH", linetype = "Spline NPH"), linewidth = 1) +
      geom_line(aes(y = lower), colour = "#619CFF", linetype = "dashed") +
      geom_line(aes(y = upper), colour = "#619CFF", linetype = "dashed") +
      geom_ribbon(aes(ymin = lower, ymax = upper), fill = "#619CFF", alpha = 0.15) +
      geom_line(aes(y = pp, colour = "Pohar-Perme", linetype = "Pohar-Perme"), linewidth = 1) +
      # scale_colour_manual(name = "Method",values = c("Spline NPH" = "#619CFF","Pohar-Perme" = "black"))+
      common_scale +  
      common_linetype + 
      coord_cartesian(ylim = c(0, 1)) +
      labs(title = "Spline NPH", x = "Time (years)", y = "Net Survival"
      ) + theme_minimal(base_family = "CMU Serif",base_size = 14) +
      theme(panel.border = element_blank(), panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
            plot.title = element_text(hjust = 0.5),
            axis.ticks = element_line(color = "black", size = 0.5),  # tick line
            axis.ticks.length = unit(0.25, "cm"),
            legend.position =  "none", #c(0.98, 0.10)
            legend.justification = c("left", "bottom"),
            legend.background = element_rect(fill = "white", color = "black")
      ) + 
      scale_x_continuous(limits = c(0, 10.5), expand = expansion(mult = 0), labels = labels) +
      scale_y_continuous(limits = c(0, 1), expand = expansion(mult = 0), breaks = function(x) {
        b <- pretty(x)
        b[b != 0]           
      }, labels = labels) 
    
    ### agrégation
    mean_plot_t_agg <- (meanWPt + meanWF1t + meanWF2t) + plot_layout(nrow = 1, guides = "collect") 
    
    a = meanWPt + theme(legend.position = "bottom")
    b = meanWF1t + theme(legend.position = "bottom")
    c = meanWF2t + theme(legend.position = "bottom")
    
    mean_plot_t_agg_leg <- (a+b+c)
    
    assign(paste0("mean_",j,"_P",q) , meanWPt, envir = globalenv())
    assign(paste0("mean_",j,"_F1",q) , meanWF1t, envir = globalenv())
    assign(paste0("mean_",j,"_F2",q) , meanWF2t, envir = globalenv())
    
  }
  
}

###########
mean_HC_Pt_plot <- mean_HC_Pt + labs(y = NULL)
mean_HR_Pt_plot <- mean_HR_Pt + labs(y = NULL)
mean_FC_Pt_plot <- mean_FC_Pt + labs(y = NULL)
mean_FR_Pt_plot <- mean_FR_Pt + labs(y = NULL)


mean_HR_Pt_plot <- mean_HR_Pt_plot + labs(title = NULL)
mean_HR_F1t_plot <- mean_HR_F1t + labs(title = NULL)
mean_HR_F2t_plot <- mean_HR_F2t + labs(title = NULL) 

mean_FC_Pt_plot <- mean_FC_Pt_plot + labs(title = NULL)
mean_FC_F1t_plot <- mean_FC_F1t+ labs(title = NULL)
mean_FC_F2t_plot <- mean_FC_F2t + labs(title = NULL)

mean_FR_Pt_plot <- mean_FR_Pt_plot + labs(title = NULL)+ labs(x = NULL)
mean_FR_F1t_plot <- mean_FR_F1t + labs(title = NULL)+ labs(x = NULL)
mean_FR_F2t_plot <- mean_FR_F2t + labs(title = NULL) + labs(x = NULL)


no_y <- theme(axis.title.y = element_blank(), axis.text.y  = element_blank(),axis.ticks.y = element_blank())

mean_HC_F1t_plot <- mean_HC_F1t + no_y
mean_HC_F2t_plot <- mean_HC_F2t + no_y 

mean_HR_F1t_plot <- mean_HR_F1t_plot + no_y
mean_HR_F2t_plot <- mean_HR_F2t_plot + no_y

mean_FC_F1t_plot <- mean_FC_F1t_plot + no_y
mean_FC_F2t_plot <- mean_FC_F2t_plot + no_y

mean_FR_F1t_plot <- mean_FR_F1t_plot + no_y
mean_FR_F2t_plot <- mean_FR_F2t_plot + no_y

no_x <- theme(axis.title.x = element_blank(), axis.text.x  = element_blank(),axis.ticks.x = element_blank())

mean_HC_Pt_plot <- mean_HC_Pt_plot + no_x
mean_HC_F1t_plot <- mean_HC_F1t_plot + no_x
mean_HC_F2t_plot <- mean_HC_F2t_plot + no_x

mean_HR_Pt_plot <- mean_HR_Pt_plot + no_x
mean_HR_F1t_plot <- mean_HR_F1t_plot + no_x
mean_HR_F2t_plot <- mean_HR_F2t_plot + no_x

mean_FC_Pt_plot <- mean_FC_Pt_plot + no_x
mean_FC_F1t_plot <- mean_FC_F1t_plot + no_x
mean_FC_F2t_plot <- mean_FC_F2t_plot + no_x

glob_plot_t_s <- ((mean_HC_Pt_plot + mean_HC_F1t_plot + mean_HC_F2t_plot) / (mean_HR_Pt_plot + mean_HR_F1t_plot + mean_HR_F2t_plot) / (mean_FC_Pt_plot + mean_FC_F1t_plot + mean_FC_F2t_plot) / (mean_FR_Pt_plot + mean_FR_F1t_plot + mean_FR_F2t_plot))

ggsave(filename = paste0("/home/thomas/Documents/Articles/Article survivalPLANN/Figures v2/Figures 4 strates/MEAN_BOOT_T_S_CONCA.png"),
       plot = glob_plot_t_s, width = 8, height = 6, dpi = 800)

########################## VALIDATION

mean_HC_Pv_plot <- mean_HC_Pv + labs(y = NULL)
mean_HR_Pv_plot <- mean_HR_Pv + labs(y = NULL)
mean_FC_Pv_plot <- mean_FC_Pv + labs(y = NULL)
mean_FR_Pv_plot <- mean_FR_Pv + labs(y = NULL)


mean_HR_Pv_plot <- mean_HR_Pv_plot + labs(title = NULL)
mean_HR_F1v_plot <- mean_HR_F1v + labs(title = NULL)
mean_HR_F2v_plot <- mean_HR_F2v + labs(title = NULL) 

mean_FC_Pv_plot <- mean_FC_Pv_plot + labs(title = NULL)
mean_FC_F1v_plot <- mean_FC_F1v+ labs(title = NULL)
mean_FC_F2v_plot <- mean_FC_F2v + labs(title = NULL)

mean_FR_Pv_plot <- mean_FR_Pv_plot + labs(title = NULL)+ labs(x = NULL)
mean_FR_F1v_plot <- mean_FR_F1v + labs(title = NULL)+ labs(x = NULL)
mean_FR_F2v_plot <- mean_FR_F2v + labs(title = NULL) + labs(x = NULL)

mean_HC_F1v_plot <- mean_HC_F1v + no_y
mean_HC_F2v_plot <- mean_HC_F2v + no_y 

mean_HR_F1v_plot <- mean_HR_F1v_plot + no_y
mean_HR_F2v_plot <- mean_HR_F2v_plot + no_y

mean_FC_F1v_plot <- mean_FC_F1v_plot + no_y
mean_FC_F2v_plot <- mean_FC_F2v_plot + no_y

mean_FR_F1v_plot <- mean_FR_F1v_plot + no_y
mean_FR_F2v_plot <- mean_FR_F2v_plot + no_y

mean_HC_Pv_plot <- mean_HC_Pv_plot + no_x
mean_HC_F1v_plot <- mean_HC_F1v_plot + no_x
mean_HC_F2v_plot <- mean_HC_F2v_plot + no_x

mean_HR_Pv_plot <- mean_HR_Pv_plot + no_x
mean_HR_F1v_plot <- mean_HR_F1v_plot + no_x
mean_HR_F2v_plot <- mean_HR_F2v_plot + no_x

mean_FC_Pv_plot <- mean_FC_Pv_plot + no_x
mean_FC_F1v_plot <- mean_FC_F1v_plot + no_x
mean_FC_F2v_plot <- mean_FC_F2v_plot + no_x

glob_plot_v_s <- ((mean_HC_Pv_plot + mean_HC_F1v_plot + mean_HC_F2v_plot) / (mean_HR_Pv_plot + mean_HR_F1v_plot + mean_HR_F2v_plot) / (mean_FC_Pv_plot + mean_FC_F1v_plot + mean_FC_F2v_plot) / (mean_FR_Pv_plot + mean_FR_F1v_plot + mean_FR_F2v_plot))

ggsave(filename = paste0("/home/thomas/Documents/Articles/Article survivalPLANN/Figures v2/Figures 4 strates/MEAN_BOOT_V_S_CONCA.png"),
       plot = glob_plot_v_s, width = 8, height = 6, dpi = 800)

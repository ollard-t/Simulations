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

path0 <- paste0("/home/thomas/Documents/Rstudio/Simulations/Simulations janvier 2026/Résultatscomplex/Output simulations_N1000_COM/")
path1 <- paste0("/home/thomas/Documents/Rstudio/Simulations/BASES_COM/LC/")

###################################################################################################
### Fonctions de survie nette, attendue, observée et de simulation de temps
###################################################################################################


Sn <- function(time, sigma, nu, theta, beta, covariates)
{
  exp((1-(1+(time/(exp(sigma)))^exp(nu))^(1/exp(theta))) * exp(sum(covariates*beta)) )
}


estimation_iterations <- function(i){
  
  set.seed(i)
  
  data <- read.csv(paste0(path1, i, "_df1000.csv"), sep = ";")
  
  data_train <- data[1:500,]
  data_valid <- data[501:1000,]
  
  newtimes = c(1,3,5,10)*365.241
  
  pro.time <- max(data_train$times[data_train$status==1])
  strata_names = c("H", "F")
  data_H <- data_train[data_train$sexchara == "male" ,]
  
  data_F <- data_train[data_train$sexchara == "female" ,]
  
  data_val_H <- data_valid[data_valid$sexchara == "male" ,]
  
  data_val_F <- data_valid[data_valid$sexchara == "female" ,]
  
  ##########
  
  tune.flex <- cvFLEXNET(formula = Surv(times,status)~ stage2 + stage3 +agey10+ sex 
                         + colon + stage2_colon + stage3_colon + colon_agey10 + ratetable(age, year, sexchara), 
                         ratetable=slopop , pro.time = pro.time,
                         data = data_train, cv = 10, m = c(0,1,2,3))
  
  tune.flex2 <- cvFLEXNET(formula = Surv(times,status)~ stage2 + stage3 + agey10+  colon + stage2_colon +
                            stage3_colon + colon_agey10 + strata(sex) + ratetable(age, year, sexchara),
                          ratetable=slopop , pro.time = pro.time,
                          data = data_train, cv = 10, m = c(0,1,2,3), m_s = c(0,1))
  
  tune.plann <- cvPLANN(formula = Surv(times, status) ~ stage2 + stage3 + agey10 + sex01 + colon
                        + stage2_colon + stage3_colon + colon_agey10, pro.time = pro.time,
                        data = data_train, cv = 10, inter= 365.241/12, size=c(6,12,24,32),
                        decay=c(0.01, 0.05, 0.1), maxit=1000, MaxNWts=10000, metric = "ibs")
  
  inter <- tune.plann$optimal$inter
  
  size <- tune.plann$optimal$size
  
  decay <- tune.plann$optimal$decay
  
  maxit <- tune.plann$optimal$maxit
  
  MaxNWts <- tune.plann$optimal$MaxNWts  
  
  m1 <- tune.flex$optimal$m
  
  m2 <- tune.flex2$optimal$m
  
  m2s <- tune.flex2$optimal$m_s
  
  
  flex.model1 <- survivalFLEXNET(formula = Surv(times,status)~ stage2 + stage3 +agey10+ sex 
                                 + colon + stage2_colon + stage3_colon + colon_agey10 + ratetable(age, year, sexchara),
                                 data = data_train, ratetable = slopop, m=m1)
  
  flex.model2 <- survivalFLEXNET(formula = Surv(times,status)~ stage2 + stage3 + agey10+  colon + stage2_colon +
                                   stage3_colon + colon_agey10 + strata(sex)  + ratetable(age, year, sexchara),
                                 data = data_train, ratetable = slopop, m=m2, m_s = m2s)
  
  plann.model <- sPLANN(formula = Surv(times,status)~ stage2 + stage3 + agey10 + sex01 + colon
                        + stage2_colon + stage3_colon + colon_agey10, pro.time = pro.time, data = data_train,
                        inter = inter, size = size, decay = decay,
                        maxit = maxit, MaxNWts = MaxNWts)
  
  
  pred.flex1 <- predict(flex.model1, newtimes = newtimes)
  
  pred.flex2 <- predict(flex.model2, newtimes = newtimes)
  
  pred.plann <- predictRS(plann.model, data_train, newtimes, slopop, "age", "year", "sexchara")
  
  PP <- summary( rs.surv(Surv(times, status) ~ 1, data = data_train,
                         ratetable = slopop, method = "pohar-perme", add.times = newtimes,
                         rmap = list(age = age, sex = sexchara, year = year)), times = newtimes, extend = TRUE) 
  
  mean.flex1 <- apply(pred.flex1$predictions, mean, MARGIN =2)
  mean.flex2 <- apply(pred.flex2$predictions, mean, MARGIN =2)
  mean.plann <- apply(pred.plann$ipredictions$relative_survival, mean, MARGIN =2)[-1]
  
  
  for(j in strata_names){
    
    data_train_var <- get(paste0("data_", j))
    
    plannpredS <- predictRS(plann.model, data = data_train_var, newtimes = newtimes, 
                            ratetable = slopop, age = "age", year = "year", sex = "sexchara")
    
    mean.plannS <- plannpredS$mpredictions$net_survival
    
    assign(paste0("plannpred_", j), plannpredS)
    assign(paste0("mean.plann_", j), mean.plannS)
  }
  
  
  for(j in strata_names){
    
    data_train_var <- get(paste0("data_", j))
    
    flexpred1S <- predict(flex.model1, newtimes = newtimes, newdata = data_train_var)$predictions
    
    mean.flex1S <- apply(flexpred1S, FUN="mean", MARGIN=2)
    
    assign(paste0("flexpred1_", j), flexpred1S)
    assign(paste0("mean.flex1_", j), mean.flex1S)
  }
  
  
  for(j in strata_names){
    
    data_train_var <- get(paste0("data_", j))
    
    flexpred2S <- predict(flex.model2, newtimes = newtimes, newdata = data_train_var)$predictions
    
    mean.flex2S <- apply(flexpred2S, FUN="mean", MARGIN=2)
    
    assign(paste0("flexpred2_", j), flexpred2S)
    assign(paste0("mean.flex2_", j), mean.flex2S)
  }
  
  for(j in strata_names){
    
    data_train_var <- get(paste0("data_", j))
    
    PPS <- summary( rs.surv(Surv(times, status) ~ 1, data = data_train_var,
                            ratetable = slopop, method = "pohar-perme", add.times = newtimes,
                            rmap = list(age = age, sex = sexchara, year = year)), times = newtimes, extend = TRUE) 
    
    poharpredS <- PPS$surv
    
    assign(paste0("PP_", j), PPS)
    assign(paste0("poharpred_", j), poharpredS)
  }
  
  pred_mat <-  matrix(data = NA, nrow = 4, ncol = 4)
  pred_mat[1,] <- mean.plann
  pred_mat[2,] <- mean.flex1
  pred_mat[3,] <- mean.flex2
  pred_mat[4,] <- PP$surv
  
  colnames(pred_mat) <- newtimes
  rownames(pred_mat) <- c("PLANN", "FLEX1", "FLEX2", "PP")
  
  stra_pred_H <- matrix(data = NA, nrow = 4, ncol = 4)
  stra_pred_H[1,] <- c(apply(plannpred_H$ipredictions$relative_survival,mean,MARGIN = 2)[-1])
  stra_pred_H[2,] <- mean.flex1_H
  stra_pred_H[3,] <- mean.flex2_H
  stra_pred_H[4,] <- poharpred_H
  
  colnames(stra_pred_H) <- newtimes
  rownames(stra_pred_H) <- c("PLANN", "FLEX1", "FLEX2", "PP")
  
  stra_pred_F <- matrix(data = NA, nrow = 4, ncol = 4)
  stra_pred_F[1,] <- c(apply(plannpred_F$ipredictions$relative_survival,mean,MARGIN = 2)[-1])
  stra_pred_F[2,] <- mean.flex1_F
  stra_pred_F[3,] <- mean.flex2_F
  stra_pred_F[4,] <- poharpred_F
  
  colnames(stra_pred_F) <- newtimes
  rownames(stra_pred_F) <- c("PLANN", "FLEX1", "FLEX2", "PP")
  
  
  ################ validation 
  
  pred.flex1_val <- predict(flex.model1, newdata = data_valid, newtimes = newtimes)
  
  pred.flex2_val <- predict(flex.model2, newdata = data_valid, newtimes = newtimes)
  
  pred.plann_val <- predictRS(plann.model, data_valid, newtimes, slopop, "age", "year", "sexchara")
  
  PP_val <- summary( rs.surv(Surv(times, status) ~ 1, data = data_valid,
                             ratetable = slopop, method = "pohar-perme", add.times = newtimes,
                             rmap = list(age = age, sex = sexchara, year = year)), times = newtimes, extend = TRUE) 
  
  mean.flex1_val <- apply(pred.flex1$predictions, mean, MARGIN =2)
  mean.flex2_val <- apply(pred.flex2$predictions, mean, MARGIN =2)
  mean.plann_val <- apply(pred.plann$ipredictions$relative_survival, mean, MARGIN =2)[-1]
  
  
  
  strata_names = c("H", "F")
  for(j in strata_names){
    
    data_train_var <- get(paste0("data_val_", j))
    
    plannpredS <- predictRS(plann.model, data = data_train_var, newtimes = newtimes, 
                            ratetable = slopop, age = "age", year = "year", sex = "sexchara")
    
    mean.plannS <- plannpredS$mpredictions$net_survival
    
    assign(paste0("plannpred_val_", j), plannpredS)
    assign(paste0("mean.plann_val_", j), mean.plannS)
  }
  
  
  for(j in strata_names){
    
    data_train_var <- get(paste0("data_val_", j))
    
    flexpred1S <- predict(flex.model1, newtimes = newtimes, newdata = data_train_var)$predictions
    
    mean.flex1S <- apply(flexpred1S, FUN="mean", MARGIN=2)
    
    assign(paste0("flexpred1_val_", j), flexpred1S)
    assign(paste0("mean.flex1_val_", j), mean.flex1S)
  }
  
  
  for(j in strata_names){
    
    data_train_var <- get(paste0("data_val_", j))
    
    flexpred2S <- predict(flex.model2, newtimes = newtimes, newdata = data_train_var)$predictions
    
    mean.flex2S <- apply(flexpred2S, FUN="mean", MARGIN=2)
    
    assign(paste0("flexpred2_val_", j), flexpred2S)
    assign(paste0("mean.flex2_val_", j), mean.flex2S)
  }
  
  for(j in strata_names){
    
    data_train_var <- get(paste0("data_val_", j))
    
    PPS <- summary( rs.surv(Surv(times, status) ~ 1, data = data_train_var,
                            ratetable = slopop, method = "pohar-perme", add.times = newtimes,
                            rmap = list(age = age, sex = sexchara, year = year)), times = newtimes, extend = TRUE) 
    
    poharpredS <- PPS$surv
    
    assign(paste0("PP_val_", j), PPS)
    assign(paste0("poharpred_val_", j), poharpredS)
  }
  
  pred_mat_val <-  matrix(data = NA, nrow = 4, ncol = 4)
  pred_mat_val[1,] <- mean.plann_val
  pred_mat_val[2,] <- mean.flex1_val
  pred_mat_val[3,] <- mean.flex2_val
  pred_mat_val[4,] <- PP_val$surv
  
  colnames(pred_mat_val) <- newtimes
  rownames(pred_mat_val) <- c("PLANN", "FLEX1", "FLEX2", "PP")
  
  stra_pred_val_H <- matrix(data = NA, nrow = 4, ncol = 4)
  stra_pred_val_H[1,] <- c(apply(plannpred_val_H$ipredictions$relative_survival,mean,MARGIN = 2)[-1])
  stra_pred_val_H[2,] <- mean.flex1_val_H
  stra_pred_val_H[3,] <- mean.flex2_val_H
  stra_pred_val_H[4,] <- poharpred_val_H
  
  colnames(stra_pred_val_H) <- newtimes
  rownames(stra_pred_val_H) <- c("PLANN", "FLEX1", "FLEX2", "PP")
  
  stra_pred_val_F <- matrix(data = NA, nrow = 4, ncol = 4)
  stra_pred_val_F[1,] <- c(apply(plannpred_val_F$ipredictions$relative_survival,mean,MARGIN = 2)[-1])
  stra_pred_val_F[2,] <- mean.flex1_val_F
  stra_pred_val_F[3,] <- mean.flex2_val_F
  stra_pred_val_F[4,] <- poharpred_val_F
  
  colnames(stra_pred_val_F) <- newtimes
  rownames(stra_pred_val_F) <- c("PLANN", "FLEX1", "FLEX2", "PP")
  
  
  ############### saving
  
  ##data
  
  write.table(data_train, paste0(path0, "DATAFRAMES/TRAIN/WHOLE/",i,"_dftrain.csv"))
  write.table(data_valid, paste0(path0, "DATAFRAMES/VALID/WHOLE/",i,"_dfvalid.csv"))
  
  for(j in strata_names){
    write.table(get(paste0("data_",j)), paste0(path0, "DATAFRAMES/TRAIN/",j,"/",i,"_dftrain",j,".csv") )
    write.table(get(paste0("data_val_",j)), paste0(path0, "DATAFRAMES/VALID/",j,"/",i,"_dfvalid",j,".csv") )
    
  }
  
  ##param
  
  flex1.coeff <- flex.model1$coefficients
  flex2.coeff <- flex.model2$coefficients
  
  paramsCV <- as.data.frame(c(inter = inter,size = size,decay = decay, maxit = maxit,
                              MaxNWts = MaxNWts, m1 = m1, coeff1 = flex1.coeff,
                              m2 = m2, m2s = m2s, coeff2 = flex2.coeff ))
  write.table(paramsCV, paste0(path0, "PARAM/",i,"_params.csv"))
  
  ##logliks
  ##train
  plann.loglik <- pred.plann$loglik  
  
  flex1.loglik <- unname(flex.model1$loglik[1])
  
  flex2.loglik <- unname(flex.model2$loglik[1])
  
  logliks <- as.data.frame(c(plann = plann.loglik, flex1 = flex1.loglik,
                             flex2 = flex2.loglik))
  ########### STRATAS
  ind_strat <- 2
  for(j in strata_names){
    
    a <- get(paste0(paste0("data_", j)))
    event_s <- a$status 
    time_s <- a$times 
    
    hP_s <- c()
    
    for(d in 1:dim(a)[1] ){
      hP_s <- c(hP_s, expectedhaz(slopop, age=a[d, "age"], sex=a[d, "sexchara"],
                                  year=a[d, "year"], time=time_s[d]) )
    }
    ### plann
    
    assign("loglik_s_plann", get(paste0("plannpred_",j))$loglik)
    
    
    ### flex 1
    #1.2
    beta_est1.2 <- flex1.coeff[1:(length(flex1.coeff)-(flex.model1$m+2))]
    gamma_est1.2 <- tail(flex1.coeff, flex.model1$m+2)
    cova_s <- as.matrix(a[,names(beta_est1.2)])
    
    
    loglik_s_1.2 <- sum( event_s * log(hP_s + (1/time_s)*splinecubeP(time_s, gamma_est1.2, flex.model1$m, flex.model1$mpos)$spln *
                                         exp(splinecube(time_s, gamma_est1.2, flex.model1$m, flex.model1$mpos)$spln + cova_s %*% beta_est1.2) ) -
                           exp(splinecube(time_s, gamma_est1.2, flex.model1$m, flex.model1$mpos)$spln 
                               + cova_s %*% beta_est1.2) )
    
    ### flex 2
    #2.2
    beta_est2.2 <- unname( flex.model2$coefficients[(1:(dim(flex.model2$x)[2]) )] ) 
    gamma_est2.2 <- unname( flex.model2$coefficients[((dim(flex.model2$x)[2]+1): (length(flex.model2$coefficients)))] ) 
    
    value <- c()
    K <- sort(unique(a$sex))
    cova_s2 <- cova_s[,attr(terms(flex.model2$formula), "term.labels")[!grepl("^(strata|ratetable)\\(", attr(terms(flex.model2$formula), "term.labels"))]]
    
    gamma_base2.2 <- gamma_est2.2[1:(flex.model2$m+2)]
    
    for(k in K){
      betak <- beta_est2.2
      gammak <- gamma_est2.2[((flex.model2$m+2)+1+(k-1)*(flex.model2$m_s+2)):((flex.model2$m+2)+(k)*(flex.model2$m_s+2))]
      idx <- a$sex == k
      
      timek <- time_s[idx]
      eventk <- event_s[idx]
      hPk <- hP_s[idx]
      covak <- cova_s2[idx, , drop = FALSE]
      # wk <- w[idx]
      splk_base <- splinecube(timek, gamma_base2.2, flex.model2$m, flex.model2$mpos)$spln
      splkP_base <- splinecubeP(timek, gamma_base2.2, flex.model2$m, flex.model2$mpos)$spln
      
      Kref <- flex.model2$Kref
      
      if (k != Kref) {
        splk  <- splinecube(timek, gammak, flex.model2$m_s, flex.model2$mpos_s[[k]])$spln #ici ça marche car pour [[k]], k et l'ordre de la list ont les meme indices mais attention si le Kref est un autre, il y aurait un décalage.
        splkP <- splinecubeP(timek, gammak, flex.model2$m_s, flex.model2$mpos_s[[k]])$spln
      } else {
        splk <- 0
        splkP <- 0
      }
      linpred <- splk_base + splk + covak %*% betak
      
      calc <- hPk + (1/timek) * (splkP_base + splkP) * exp(linpred)
      
      ##[calc >= 0] pour éviter des NaN quand les parametres estimés sur une autre base pourraient créer un spln < 0 (surtout sur base de validation)
      value_k <- sum( (
        eventk[calc >= 0] * log( calc[calc >= 0] )) -
          exp(linpred[calc >= 0]))
      
      value <- c(value, value_k)
    }
    
    loglik_s_2.2 <- sum(value)
    
    logliks[,ind_strat] <- c(loglik_s_plann, loglik_s_1.2, loglik_s_2.2)
    ind_strat <- ind_strat+1
    
  }
  
  colnames(logliks) <- c("WHOLE", strata_names)
  
  write.table(logliks,  paste0(path0, "LOGLIK/TRAIN/",i,"_loglik.csv"))
  
  ##logliks
  ##validation
  event_valid <- data_valid$status 
  time_valid <- data_valid$times 
  
  hP_valid <- c()
  
  for(d in 1:dim(data_valid)[1] ){
    hP_valid <- c(hP_valid, expectedhaz(slopop, age=data_valid[d, "age"], sex=data_valid[d, "sexchara"],
                                        year=data_valid[d, "year"], time=time_valid[d]) )
  }
  plannval.loglik <- pred.plann_val$loglik  
  
  
  cova_valid <- as.matrix(data_valid[,names(beta_est1.2)])
  
  loglikval_1.2 <- sum( event_valid * log(hP_valid + (1/time_valid)*splinecubeP(time_valid, gamma_est1.2, flex.model1$m, flex.model1$mpos)$spln *
                                            exp(splinecube(time_valid, gamma_est1.2, flex.model1$m, flex.model1$mpos)$spln + cova_valid %*% beta_est1.2) ) -
                          exp(splinecube(time_valid, gamma_est1.2, flex.model1$m, flex.model1$mpos)$spln 
                              + cova_valid %*% beta_est1.2) )
  
  
  value <- c()
  K <- sort(unique(data_valid$sex))
  cova_valid2 <- cova_valid[,attr(terms(flex.model2$formula), "term.labels")[!grepl("^(strata|ratetable)\\(", attr(terms(flex.model2$formula), "term.labels"))]]
  
  for(k in K){
    betak <- beta_est2.2
    gammak <- gamma_est2.2[((flex.model2$m+2)+1+(k-1)*(flex.model2$m_s+2)):((flex.model2$m+2)+(k)*(flex.model2$m_s+2))]
    idx <- data_valid$sex == k
    
    timek <- time_valid[idx]
    eventk <- event_valid[idx]
    hPk <- hP_valid[idx]
    covak <- cova_valid2[idx, , drop = FALSE]
    # wk <- w[idx]
    splk_base <- splinecube(timek, gamma_base2.2, flex.model2$m, flex.model2$mpos)$spln
    splkP_base <- splinecubeP(timek, gamma_base2.2, flex.model2$m, flex.model2$mpos)$spln
    
    Kref <- flex.model2$Kref
    
    if (k != Kref) {
      splk  <- splinecube(timek, gammak, flex.model2$m_s, flex.model2$mpos_s[[k]])$spln
      splkP <- splinecubeP(timek, gammak, flex.model2$m_s, flex.model2$mpos_s[[k]])$spln
    } else {
      splk <- 0
      splkP <- 0
    }
    linpred <- splk_base + splk + covak %*% betak
    
    calc <- hPk + (1/timek) * (splkP_base + splkP) * exp(linpred)
    
    ##[calc >= 0] pour éviter des NaN quand les parametres estimés sur une autre base pourraient créer un spln < 0 (surtout sur base de validation)
    value_k <- sum( (
      eventk[calc >= 0] * log( calc[calc >= 0] )) -
        exp(linpred[calc >= 0]))
    
    value <- c(value, value_k)
  }
  
  loglikval_2.2 <- sum(value)
  
  logliksval <- as.data.frame(c(plann = plannval.loglik, flex1 = loglikval_1.2,
                                flex2 = loglikval_2.2))
  ########### STRATAS
  ind_strat <- 2
  for(j in strata_names){
    
    a <- get(paste0(paste0("data_val_", j)))
    event_s <- a$status 
    time_s <- a$times 
    
    hP_s <- c()
    
    for(d in 1:dim(a)[1] ){
      hP_s <- c(hP_s, expectedhaz(slopop, age=a[d, "age"], sex=a[d, "sexchara"],
                                  year=a[d, "year"], time=time_s[d]) )
    }
    ### plann
    
    assign("loglik_s_plann", get(paste0("plannpred_val_",j))$loglik)
    
    
    ### flex 1
    #1.2
    beta_est1.2 <- flex1.coeff[1:(length(flex1.coeff)-(flex.model1$m+2))]
    gamma_est1.2 <- tail(flex1.coeff, flex.model1$m+2)
    cova_s <- as.matrix(a[,names(beta_est1.2)])
    
    
    loglik_s_1.2 <- sum( event_s * log(hP_s + (1/time_s)*splinecubeP(time_s, gamma_est1.2, flex.model1$m, flex.model1$mpos)$spln *
                                         exp(splinecube(time_s, gamma_est1.2, flex.model1$m, flex.model1$mpos)$spln + cova_s %*% beta_est1.2) ) -
                           exp(splinecube(time_s, gamma_est1.2, flex.model1$m, flex.model1$mpos)$spln 
                               + cova_s %*% beta_est1.2) )
    
    ### flex 2
    #2.2
    beta_est2.2 <- unname( flex.model2$coefficients[(1:(dim(flex.model2$x)[2]) )] ) 
    gamma_est2.2 <- unname( flex.model2$coefficients[((dim(flex.model2$x)[2]+1): (length(flex.model2$coefficients)))] ) 
    
    value <- c()
    K <- sort(unique(a$sex))
    cova_s2 <- cova_s[,attr(terms(flex.model2$formula), "term.labels")[!grepl("^(strata|ratetable)\\(", attr(terms(flex.model2$formula), "term.labels"))]]
    
    gamma_base2.2 <- gamma_est2.2[1:(flex.model2$m+2)]
    
    for(k in K){
      betak <- beta_est2.2
      gammak <- gamma_est2.2[((flex.model2$m+2)+1+(k-1)*(flex.model2$m_s+2)):((flex.model2$m+2)+(k)*(flex.model2$m_s+2))]
      idx <- a$sex == k
      
      timek <- time_s[idx]
      eventk <- event_s[idx]
      hPk <- hP_s[idx]
      covak <- cova_s2[idx, , drop = FALSE]
      # wk <- w[idx]
      splk_base <- splinecube(timek, gamma_base2.2, flex.model2$m, flex.model2$mpos)$spln
      splkP_base <- splinecubeP(timek, gamma_base2.2, flex.model2$m, flex.model2$mpos)$spln
      
      Kref <- flex.model2$Kref
      
      if (k != Kref) {
        splk  <- splinecube(timek, gammak, flex.model2$m_s, flex.model2$mpos_s[[k]])$spln #ici ça marche car pour [[k]], k et l'ordre de la list ont les meme indices mais attention si le Kref est un autre, il y aurait un décalage.
        splkP <- splinecubeP(timek, gammak, flex.model2$m_s, flex.model2$mpos_s[[k]])$spln
      } else {
        splk <- 0
        splkP <- 0
      }
      linpred <- splk_base + splk + covak %*% betak
      
      calc <- hPk + (1/timek) * (splkP_base + splkP) * exp(linpred)
      
      ##[calc >= 0] pour éviter des NaN quand les parametres estimés sur une autre base pourraient créer un spln < 0 (surtout sur base de validation)
      value_k <- sum( (
        eventk[calc >= 0] * log( calc[calc >= 0] )) -
          exp(linpred[calc >= 0]))
      
      value <- c(value, value_k)
    }
    
    loglik_s_2.2 <- sum(value)
    
    logliksval[,ind_strat] <- c(loglik_s_plann, loglik_s_1.2, loglik_s_2.2)
    ind_strat <- ind_strat+1
    
  }
  
  colnames(logliksval) <- c("WHOLE", strata_names)
  
  write.table(logliksval,  paste0(path0, "LOGLIK/VALID/",i,"_loglikval.csv"))
  
  
  ##train
  
  write.table(pred.plann$ipredictions$relative_survival[,-1], paste0(path0, "TRAIN/WHOLE/ind/plann_",i,".csv"))
  write.table(pred.flex1$predictions, paste0(path0, "TRAIN/WHOLE/ind/flex1_",i,".csv"))
  write.table(pred.flex2$predictions, paste0(path0, "TRAIN/WHOLE/ind/flex2_",i,".csv"))
  
  write.table(pred_mat, paste0(path0, "TRAIN/WHOLE/mean/pred_mat_",i,".csv"))
  
  for(j in strata_names){
    write.table(get(paste0("plannpred_",j))$ipredictions$relative_survival[,-1], paste0(path0, "TRAIN/STRATA/ind/plann_",i,"_",j,".csv"), sep = ";", col.names = T)
    write.table(get(paste0("flexpred1_",j)), paste0(path0, "TRAIN/STRATA/ind/flex1_",i,"_",j,".csv"))
    write.table(get(paste0("flexpred2_",j)), paste0(path0, "TRAIN/STRATA/ind/flex2_",i,"_",j,".csv"))
    
    write.table(get(paste0("stra_pred_",j)), paste0(path0, "TRAIN/STRATA/mean/pred_mat_",i,"_",j,".csv"))
  }
  
  #### validation
  
  write.table(pred.plann_val$ipredictions$relative_survival[,-1], paste0(path0, "VALID/WHOLE/ind/plann_val_",i,".csv"))
  write.table(pred.flex1_val$predictions, paste0(path0, "VALID/WHOLE/ind/flex1_val_",i,".csv"))
  write.table(pred.flex2_val$predictions, paste0(path0, "VALID/WHOLE/ind/flex2_val_",i,".csv"))
  
  write.table(pred_mat_val, paste0(path0, "VALID/WHOLE/mean/pred_mat_val_",i,".csv"))
  
  for(j in strata_names){
    write.table(get(paste0("plannpred_val_",j))$ipredictions$relative_survival[,-1], paste0(path0, "VALID/STRATA/ind/plann_",i,"_",j,".csv"), sep = ";", col.names = T)
    write.table(get(paste0("flexpred1_val_",j)), paste0(path0, "VALID/STRATA/ind/flex1_val_",i,"_",j,".csv"))
    write.table(get(paste0("flexpred2_val_",j)), paste0(path0, "VALID/STRATA/ind/flex2_val_",i,"_",j,".csv"))
    
    write.table(get(paste0("stra_pred_val_",j)), paste0(path0, "VALID/STRATA/mean/pred_mat_val_",i,"_",j,".csv"))
  }
}

##utiliser code Indics_complexes plutot Github/Simulations/Indics_complexes

iterations = 1:1000

mclapply(iterations, estimation_iterations, mc.cores = detectCores() - 5)


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
  mclapply(indic, estimation_iterations, N = 1000, mc.cores = detectCores() - 4)
  attempt <- attempt + 1
}

rm(list = setdiff(ls(envir = .GlobalEnv), 
                  c("path0", "date_launch", "fr.ratetable", "slopop", 
                    "colrec", "Sn", "estimation_iterations", "calc_indic")), envir = .GlobalEnv) ## cleaning de l'environnement excpeté ce qu iva être réutilisé
## 1000SiEND

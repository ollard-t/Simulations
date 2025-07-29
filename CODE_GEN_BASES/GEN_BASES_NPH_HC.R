library(survivalNET)
library(survivalPLANN)
library(parallel)
library(doParallel)

#################################################################################################
### path ###

path0 <- "~/Documents/Simulations/BASES/ind1000NPH_HC/"
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

colHC <- colrec[colrec$sex == 1 & colrec$colon == 1,]
colHR <- colrec[colrec$sex == 1 & colrec$colon == 0,]
colFC <- colrec[colrec$sex == 2 & colrec$colon == 1,]
colFR <- colrec[colrec$sex == 2 & colrec$colon == 0,]

rsHC <- survivalNET(formula = Surv(time, stat) ~ stage2 + stage3 + agey10 +   
                      ratetable(age, diag, sexchara), data = colHC,
                    ratetable=slopop, dist="weibull", weights=NULL) 


rsHR <- survivalNET(formula = Surv(time, stat) ~ stage2 + stage3 + agey10 + 
                      ratetable(age, diag, sexchara), data = colHR,
                    ratetable=slopop, dist="weibull", weights=NULL) 

rsFC <- survivalNET(formula = Surv(time, stat) ~ stage2 + stage3 + agey10 + 
                      ratetable(age, diag, sexchara), data = colFC,
                    ratetable=slopop, dist="weibull", weights=NULL) 

rsFR <- survivalNET(formula = Surv(time, stat) ~ stage2 + stage3 + agey10 +  
                      ratetable(age, diag, sexchara), data = colFR,
                    ratetable=slopop, dist="weibull", weights=NULL) 

###################################################################################################
### Fonctions de survie nette, attendue, observée et de simulation de temps
###################################################################################################

Sp <- function(time, ratetable, age, sex, year)
{
  return(exp(-expectedcumhaz(ratetable, age, year, sex, time, method = "table")))
}

Sn <- function(time, sigma, nu, theta, beta, covariates)
{
  exp((1-(1+(time/(exp(sigma)))^exp(nu))^(1/exp(theta))) * exp(sum(covariates*beta)) )
}

Sobs <- function(time, sigma, nu, theta, beta, covariates, ratetable, age, sex, year)
{
  Sp(time = time, ratetable = ratetable, age = age, year = year, sex = sex) * 
    Sn(time = time, sigma = sigma, nu = nu, theta = theta, beta = beta, covariates = covariates)
}

Tsim <- function(sigma, nu, theta, covariates, beta, ratetable, age, sex, year)
{
  stimes <- lapply(1:length(age), function(i)
  {
    u <- runif(1)
    time <- try(uniroot(function(x) 1- Sobs(time=x, sigma=sigma, nu=nu, theta=theta, beta=beta, covariates=covariates[i,], 
                                            ratetable=ratetable, age=age[i], sex=sex[i], year=year[i])-u,
                        interval = c(.Machine$double.eps, 130000))$root, silent = TRUE) 
    return(time)
  })
}
###################################################################################################
##############
# simulation#
#############

simulate_iteration <- function(i){
  
  set.seed(i)
  
  start <- Sys.time()
  
  N <- 1000
  
  ###paramètres arrondis
  #### Hommes Colon
  TsigmaHC <- 13.0
  TnuHC <- -0.55
  TthetaHC <- 0
  betaZHC <- c(1,2.9,0.2)
  names(betaZHC) <- names(rsHC$coefficients[1:3])
  
  #### Hommes Rectum
  TsigmaHR <- 12.75
  TnuHR <- -0.3
  TthetaHR <- 0
  betaZHR <- c(1.2,3,0.35)
  names(betaZHR) <- names(rsHR$coefficients[1:3])
  
  #### Femmes Colon
  TsigmaFC <- 14 
  TnuFC <- -0.55
  TthetaFC <- 0
  betaZFC <- c(0.6,2.6,-0.1)
  names(betaZFC) <- names(rsFC$coefficients[1:3])
  
  #### Femmes Rectum
  TsigmaFR <- 11.5
  TnuFR <- -0.4
  TthetaFR <- 0
  betaZFR <- c(0.9,2.7,0.16)
  names(betaZFR) <- names(rsFR$coefficients[1:3])
  
  ######
  
  data <- data.frame(sex = 1 + rbinom(N, 1, 0.5))
  
  data$sexchara <- "male"
  data$sexchara[data$sex==2] <- "female"
  
  data$sex01 <- (data$sex)-1
  
  data$agey <- floor(rnorm(N, mean=65, sd=10))
  data$agey10 <- data$agey /10
  data$age <- data$agey*365.241
  
  # glm.colon <- glm(colon ~ sex + agey10, data = colrec, family = "binomial")
  
  lp1 <- -0.4 + 0.2 * data$sex + 0.07 * data$agey10
  
  data$colon <- rbinom(N, 1, prob =  exp(lp1 ) / (1 + exp(lp1)) )
  
  # glm.stage1 <- glm(stage1 ~ sex + agey10 + colon, data = colrec, family = "binomial")
  
  lp2 <- -1.9 + 0.1 * data$sex + 0.06 * data$agey10 - 0.7 * data$colon
  
  data$stage1 <- rbinom(N, 1, prob =  exp(lp2) / (1 + exp(lp2)) )
  
  # glm.stage2 <- glm(stage2 ~ sex + agey10 + colon, data = colrec[colrec$stage1==0,], family = "binomial")
  
  lp3 <- 0.7 + 0.2 * data$sex + 0.02 * data$agey10 - 0.2 * data$colon
  
  data$stage2 <- 0
  
  data$stage2[data$stage1==0] <- rbinom(sum(data$stage1==0), 1,
                                        prob =  exp(lp3[data$stage1==0]) / (1 + exp(lp3[data$stage1==0])) )
  
  data$stage3 <- 1 - data$stage1 - data$stage2
  
  data$stage <- 1 + 1*(data$stage2==1) + 2*(data$stage3==1)
  
  data$year <- as.numeric(sample(seq(
    as.Date('1994/01/01', origin = "1960-01-01"),
    as.Date('2000/12/30', origin = "1960-01-01"),
    by="day"), N, replace = TRUE) - as.Date("1960-01-01"))
  
  Z <- cbind(data$stage2, data$stage3, data$agey10, data$sex, data$colon)
  colnames(Z) <- c("stage2", "stage3", "agey10", "sex", "colon")
  
  ## Hommes colon
  dataHC <- data[data$sex ==1 & data$colon == 1,]
  ZHC <- cbind(dataHC$stage2, dataHC$stage3, dataHC$agey10)
  colnames(ZHC) <- c("stage2", "stage3", "agey10")
  
  data[data$sex ==1 & data$colon == 1, "timesT"] <- unlist(Tsim(sigma = TsigmaHC, nu = TnuHC, theta = TthetaHC,
                                                                beta = betaZHC, covariates = ZHC,  
                                                                ratetable = slopop,
                                                                age = dataHC$age, sex = dataHC$sexchara, year = dataHC$year))
  
  ## Hommes rectum
  dataHR <- data[data$sex ==1 & data$colon == 0,]
  ZHR <- cbind(dataHR$stage2, dataHR$stage3, dataHR$agey10)
  colnames(ZHR) <- c("stage2", "stage3", "agey10")
  
  data[data$sex ==1 & data$colon == 0, "timesT"] <- unlist(Tsim(sigma = TsigmaHR, nu = TnuHR, theta = TthetaHR,
                                                                beta = betaZHR, covariates = ZHR,  
                                                                ratetable = slopop,
                                                                age = dataHR$age, sex = dataHR$sexchara, year = dataHR$year) )
  
  ## Femmes Colon 
  dataFC <- data[data$sex ==2 & data$colon == 1,]
  ZFC <- cbind(dataFC$stage2, dataFC$stage3, dataFC$agey10) 
  colnames(ZFC) <- c("stage2", "stage3", "agey10")
  
  data[data$sex == 2 & data$colon == 1, "timesT"] <- unlist(Tsim(sigma = TsigmaFC, nu = TnuFC, theta = TthetaFC,
                                                                 beta = betaZFC, covariates = ZFC,  
                                                                 ratetable = slopop,
                                                                 age = dataFC$age, sex = dataFC$sexchara, year = dataFC$year))
  
  ## Femmes Rectum 
  dataFR <- data[data$sex ==2 & data$colon == 0,]
  ZFR <- cbind(dataFR$stage2, dataFR$stage3, dataFR$agey10)
  colnames(ZFR) <- c("stage2", "stage3", "agey10")  
  data[data$sex ==2 & data$colon == 0, "timesT"] <- unlist(Tsim(sigma = TsigmaFR, nu = TnuFR, theta = TthetaFR,
                                                                beta = betaZFR, covariates = ZFR,  
                                                                ratetable = slopop,
                                                                age = dataFR$age, sex = dataFR$sexchara, year = dataFR$year))
  
  
  data$timesT <- as.numeric(data$timesT)
  
  
  if (anyNA(data$timesT)) {
    while(anyNA(data$timesT) != FALSE ){
      na_indices <- which(is.na(data$timesT))
      
      for (p in na_indices) {
        # Hommes colon
        if(data[p,]$sex == 1 & data[p,]$colon == 1){
          Tsigma = TsigmaHC
          Tnu = TnuHC
          Ttheta = TthetaHC
          betaZ = betaZHC
        }
        # Hommes rectum
        if(data[p,]$sex == 1 & data[p,]$colon == 0){
          Tsigma = TsigmaHR
          Tnu = TnuHR
          Ttheta = TthetaHR
          betaZ = betaZHR
        }
        # Femmes colon
        if(data[p,]$sex == 2 & data[p,]$colon == 1){
          Tsigma = TsigmaFC
          Tnu = TnuFC
          Ttheta = TthetaFC
          betaZ = betaZFC
        }
        # Femmes rectum
        if(data[p,]$sex == 2 & data[p,]$colon == 0){
          Tsigma = TsigmaFR
          Tnu = TnuFR
          Ttheta = TthetaFR
          betaZ = betaZFR
        }
        data$timesT[p] <- Tsim(sigma = Tsigma, nu = Tnu, theta = Ttheta,
                               beta = betaZ, covariates = Z[p, -c(4,5), drop = FALSE],  
                               ratetable = slopop,
                               age = data$age[p], sex = data$sexchara[p], year = data$year[p]) 
      }
      data$timesT <- as.numeric(data$timesT)
      
    }
  }
  data$timesT <- as.numeric(data$timesT)
  
  #data <- na.omit(data)
  #N <- dim(data)[1] 
  ################## génération des temps de censure
  ZC <- Z[,c("agey10", "sex")] 
  
  Csim <- function(u, sigma, betaC, ZC) { -1*log(1-u)*(sigma*365.241)/exp(ZC%*%betaC) }
  
  sigmaC <- 13.5 ### 75 -> 7.5% ; 55 -> 10% ;  13.5-> 30%; 6 <- 35%, <10% d'evt au dernier temps;  5 -> ~40% ; 2.5 -> ~50% avec pas assez d'ind à risque à 10 ans 
  betaC <- c(-0.01,-0.3)
  
  data$timesC <- Csim(u = runif(n = N, min=0, max=1), sigma = sigmaC, betaC = betaC, ZC = ZC) 
  
  # boxplot(data$timesC ~ data$sex, main = "Censoring Times by Sex")
  # plot(data$age, data$timesC, main = "Censoring Times vs Age", xlab = "Age", ylab = "Censoring Time")
  
  data$times <- pmin(data$timesT, data$timesC)
  
  data$status <- 1*(data$timesT <= data$timesC)
  
  data$sex.organ <- paste0(data$sexchara, data$colon)
  data$sex.organ <- as.numeric(factor(data$sex.organ, levels = c("male0", "male1", "female0", "female1")))
  
  data$times[data$times < 10e-4] <- 10e-4
  write.table(data, paste0(path0,i,"_df",N,".csv"), sep = ';', row.names = F, col.names = T)
}


iterations <- 1:1000

mclapply(iterations, simulate_iteration, mc.cores = detectCores() - 10)


#################################################################################################

#####################################      3000           #######################################

#################################################################################################
### path ###

path0 <- "~/Documents/Simulations/BASES/ind3000NPH_HC/"
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

colHC <- colrec[colrec$sex == 1 & colrec$colon == 1,]
colHR <- colrec[colrec$sex == 1 & colrec$colon == 0,]
colFC <- colrec[colrec$sex == 2 & colrec$colon == 1,]
colFR <- colrec[colrec$sex == 2 & colrec$colon == 0,]

rsHC <- survivalNET(formula = Surv(time, stat) ~ stage2 + stage3 + agey10 +   
                      ratetable(age, diag, sexchara), data = colHC,
                    ratetable=slopop, dist="weibull", weights=NULL) 


rsHR <- survivalNET(formula = Surv(time, stat) ~ stage2 + stage3 + agey10 + 
                      ratetable(age, diag, sexchara), data = colHR,
                    ratetable=slopop, dist="weibull", weights=NULL) 

rsFC <- survivalNET(formula = Surv(time, stat) ~ stage2 + stage3 + agey10 + 
                      ratetable(age, diag, sexchara), data = colFC,
                    ratetable=slopop, dist="weibull", weights=NULL) 

rsFR <- survivalNET(formula = Surv(time, stat) ~ stage2 + stage3 + agey10 +  
                      ratetable(age, diag, sexchara), data = colFR,
                    ratetable=slopop, dist="weibull", weights=NULL) 

###################################################################################################
### Fonctions de survie nette, attendue, observée et de simulation de temps
###################################################################################################

Sp <- function(time, ratetable, age, sex, year)
{
  return(exp(-expectedcumhaz(ratetable, age, year, sex, time, method = "table")))
}

Sn <- function(time, sigma, nu, theta, beta, covariates)
{
  exp((1-(1+(time/(exp(sigma)))^exp(nu))^(1/exp(theta))) * exp(sum(covariates*beta)) )
}

Sobs <- function(time, sigma, nu, theta, beta, covariates, ratetable, age, sex, year)
{
  Sp(time = time, ratetable = ratetable, age = age, year = year, sex = sex) * 
    Sn(time = time, sigma = sigma, nu = nu, theta = theta, beta = beta, covariates = covariates)
}

Tsim <- function(sigma, nu, theta, covariates, beta, ratetable, age, sex, year)
{
  stimes <- lapply(1:length(age), function(i)
  {
    u <- runif(1)
    time <- try(uniroot(function(x) 1- Sobs(time=x, sigma=sigma, nu=nu, theta=theta, beta=beta, covariates=covariates[i,], 
                                            ratetable=ratetable, age=age[i], sex=sex[i], year=year[i])-u,
                        interval = c(.Machine$double.eps, 130000))$root, silent = TRUE) 
    return(time)
  })
}
###################################################################################################
##############
# simulation#
#############

simulate_iteration <- function(i){
  
  set.seed(i)
  
  start <- Sys.time()
  
  N <- 3000
  
  ###paramètres arrondis
  #### Hommes Colon
  TsigmaHC <- 13.0
  TnuHC <- -0.55
  TthetaHC <- 0
  betaZHC <- c(1,2.9,0.2)
  names(betaZHC) <- names(rsHC$coefficients[1:3])
  
  #### Hommes Rectum
  TsigmaHR <- 12.75
  TnuHR <- -0.3
  TthetaHR <- 0
  betaZHR <- c(1.2,3,0.35)
  names(betaZHR) <- names(rsHR$coefficients[1:3])
  
  #### Femmes Colon
  TsigmaFC <- 14 
  TnuFC <- -0.55
  TthetaFC <- 0
  betaZFC <- c(0.6,2.6,-0.1)
  names(betaZFC) <- names(rsFC$coefficients[1:3])
  
  #### Femmes Rectum
  TsigmaFR <- 11.5
  TnuFR <- -0.4
  TthetaFR <- 0
  betaZFR <- c(0.9,2.7,0.16)
  names(betaZFR) <- names(rsFR$coefficients[1:3])
  
  ######
  
  data <- data.frame(sex = 1 + rbinom(N, 1, 0.5))
  
  data$sexchara <- "male"
  data$sexchara[data$sex==2] <- "female"
  
  data$sex01 <- (data$sex)-1
  
  data$agey <- floor(rnorm(N, mean=65, sd=10))
  data$agey10 <- data$agey /10
  data$age <- data$agey*365.241
  
  # glm.colon <- glm(colon ~ sex + agey10, data = colrec, family = "binomial")
  
  lp1 <- -0.4 + 0.2 * data$sex + 0.07 * data$agey10
  
  data$colon <- rbinom(N, 1, prob =  exp(lp1 ) / (1 + exp(lp1)) )
  
  # glm.stage1 <- glm(stage1 ~ sex + agey10 + colon, data = colrec, family = "binomial")
  
  lp2 <- -1.9 + 0.1 * data$sex + 0.06 * data$agey10 - 0.7 * data$colon
  
  data$stage1 <- rbinom(N, 1, prob =  exp(lp2) / (1 + exp(lp2)) )
  
  # glm.stage2 <- glm(stage2 ~ sex + agey10 + colon, data = colrec[colrec$stage1==0,], family = "binomial")
  
  lp3 <- 0.7 + 0.2 * data$sex + 0.02 * data$agey10 - 0.2 * data$colon
  
  data$stage2 <- 0
  
  data$stage2[data$stage1==0] <- rbinom(sum(data$stage1==0), 1,
                                        prob =  exp(lp3[data$stage1==0]) / (1 + exp(lp3[data$stage1==0])) )
  
  data$stage3 <- 1 - data$stage1 - data$stage2
  
  data$stage <- 1 + 1*(data$stage2==1) + 2*(data$stage3==1)
  
  data$year <- as.numeric(sample(seq(
    as.Date('1994/01/01', origin = "1960-01-01"),
    as.Date('2000/12/30', origin = "1960-01-01"),
    by="day"), N, replace = TRUE) - as.Date("1960-01-01"))
  
  Z <- cbind(data$stage2, data$stage3, data$agey10, data$sex, data$colon)
  colnames(Z) <- c("stage2", "stage3", "agey10", "sex", "colon")
  
  ## Hommes colon
  dataHC <- data[data$sex ==1 & data$colon == 1,]
  ZHC <- cbind(dataHC$stage2, dataHC$stage3, dataHC$agey10)
  colnames(ZHC) <- c("stage2", "stage3", "agey10")
  
  data[data$sex ==1 & data$colon == 1, "timesT"] <- unlist(Tsim(sigma = TsigmaHC, nu = TnuHC, theta = TthetaHC,
                                                                beta = betaZHC, covariates = ZHC,  
                                                                ratetable = slopop,
                                                                age = dataHC$age, sex = dataHC$sexchara, year = dataHC$year))
  
  ## Hommes rectum
  dataHR <- data[data$sex ==1 & data$colon == 0,]
  ZHR <- cbind(dataHR$stage2, dataHR$stage3, dataHR$agey10)
  colnames(ZHR) <- c("stage2", "stage3", "agey10")
  
  data[data$sex ==1 & data$colon == 0, "timesT"] <- unlist(Tsim(sigma = TsigmaHR, nu = TnuHR, theta = TthetaHR,
                                                                beta = betaZHR, covariates = ZHR,  
                                                                ratetable = slopop,
                                                                age = dataHR$age, sex = dataHR$sexchara, year = dataHR$year) )
  
  ## Femmes Colon 
  dataFC <- data[data$sex ==2 & data$colon == 1,]
  ZFC <- cbind(dataFC$stage2, dataFC$stage3, dataFC$agey10) 
  colnames(ZFC) <- c("stage2", "stage3", "agey10")  
  data[data$sex == 2 & data$colon == 1, "timesT"] <- unlist(Tsim(sigma = TsigmaFC, nu = TnuFC, theta = TthetaFC,
                                                                 beta = betaZFC, covariates = ZFC,  
                                                                 ratetable = slopop,
                                                                 age = dataFC$age, sex = dataFC$sexchara, year = dataFC$year))
  
  ## Femmes Rectum 
  dataFR <- data[data$sex ==2 & data$colon == 0,]
  ZFR <- cbind(dataFR$stage2, dataFR$stage3, dataFR$agey10)
  colnames(ZFR) <- c("stage2", "stage3", "agey10")  
  
  data[data$sex ==2 & data$colon == 0, "timesT"] <- unlist(Tsim(sigma = TsigmaFR, nu = TnuFR, theta = TthetaFR,
                                                                beta = betaZFR, covariates = ZFR,  
                                                                ratetable = slopop,
                                                                age = dataFR$age, sex = dataFR$sexchara, year = dataFR$year))
  
  
  data$timesT <- as.numeric(data$timesT)
  
  
  if (anyNA(data$timesT)) {
    while(anyNA(data$timesT) != FALSE ){
      na_indices <- which(is.na(data$timesT))
      
      for (p in na_indices) {
        # Hommes colon
        if(data[p,]$sex == 1 & data[p,]$colon == 1){
          Tsigma = TsigmaHC
          Tnu = TnuHC
          Ttheta = TthetaHC
          betaZ = betaZHC
        }
        # Hommes rectum
        if(data[p,]$sex == 1 & data[p,]$colon == 0){
          Tsigma = TsigmaHR
          Tnu = TnuHR
          Ttheta = TthetaHR
          betaZ = betaZHR
        }
        # Femmes colon
        if(data[p,]$sex == 2 & data[p,]$colon == 1){
          Tsigma = TsigmaFC
          Tnu = TnuFC
          Ttheta = TthetaFC
          betaZ = betaZFC
        }
        # Femmes rectum
        if(data[p,]$sex == 2 & data[p,]$colon == 0){
          Tsigma = TsigmaFR
          Tnu = TnuFR
          Ttheta = TthetaFR
          betaZ = betaZFR
        }
        data$timesT[p] <- Tsim(sigma = Tsigma, nu = Tnu, theta = Ttheta,
                               beta = betaZ, covariates = Z[p, -c(4,5), drop = FALSE],  
                               ratetable = slopop,
                               age = data$age[p], sex = data$sexchara[p], year = data$year[p]) 
      }
      data$timesT <- as.numeric(data$timesT)
      
    }
  }
  data$timesT <- as.numeric(data$timesT)
  
  #data <- na.omit(data)
  #N <- dim(data)[1] 
  ################## génération des temps de censure
  ZC <- Z[,c("agey10", "sex")] 
  
  Csim <- function(u, sigma, betaC, ZC) { -1*log(1-u)*(sigma*365.241)/exp(ZC%*%betaC) }
  
  sigmaC <- 13.5 ### 75 -> 7.5% ; 55 -> 10% ;  13.5-> 30%; 6 <- 35%, <10% d'evt au dernier temps;  5 -> ~40% ; 2.5 -> ~50% avec pas assez d'ind à risque à 10 ans 
  betaC <- c(-0.01,-0.3)
  
  data$timesC <- Csim(u = runif(n = N, min=0, max=1), sigma = sigmaC, betaC = betaC, ZC = ZC) 
  
  # boxplot(data$timesC ~ data$sex, main = "Censoring Times by Sex")
  # plot(data$age, data$timesC, main = "Censoring Times vs Age", xlab = "Age", ylab = "Censoring Time")
  
  data$times <- pmin(data$timesT, data$timesC)
  
  data$status <- 1*(data$timesT <= data$timesC)
  
  data$sex.organ <- paste0(data$sexchara, data$colon)
  data$sex.organ <- as.numeric(factor(data$sex.organ, levels = c("male0", "male1", "female0", "female1")))
  
  data$times[data$times < 10e-4] <- 10e-4
  write.table(data, paste0(path0,i,"_df",N,".csv"), sep = ';', row.names = F, col.names = T)
}


iterations <- 1:1000

mclapply(iterations, simulate_iteration, mc.cores = detectCores() - 10)



#################################################################################################

#####################################      5000           #######################################

#################################################################################################
### path ###

path0 <- "~/Documents/Simulations/BASES/ind5000NPH_HC/"
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

colHC <- colrec[colrec$sex == 1 & colrec$colon == 1,]
colHR <- colrec[colrec$sex == 1 & colrec$colon == 0,]
colFC <- colrec[colrec$sex == 2 & colrec$colon == 1,]
colFR <- colrec[colrec$sex == 2 & colrec$colon == 0,]

rsHC <- survivalNET(formula = Surv(time, stat) ~ stage2 + stage3 + agey10 +   
                      ratetable(age, diag, sexchara), data = colHC,
                    ratetable=slopop, dist="weibull", weights=NULL) 


rsHR <- survivalNET(formula = Surv(time, stat) ~ stage2 + stage3 + agey10 + 
                      ratetable(age, diag, sexchara), data = colHR,
                    ratetable=slopop, dist="weibull", weights=NULL) 

rsFC <- survivalNET(formula = Surv(time, stat) ~ stage2 + stage3 + agey10 + 
                      ratetable(age, diag, sexchara), data = colFC,
                    ratetable=slopop, dist="weibull", weights=NULL) 

rsFR <- survivalNET(formula = Surv(time, stat) ~ stage2 + stage3 + agey10 +  
                      ratetable(age, diag, sexchara), data = colFR,
                    ratetable=slopop, dist="weibull", weights=NULL) 

###################################################################################################
### Fonctions de survie nette, attendue, observée et de simulation de temps
###################################################################################################

Sp <- function(time, ratetable, age, sex, year)
{
  return(exp(-expectedcumhaz(ratetable, age, year, sex, time, method = "table")))
}

Sn <- function(time, sigma, nu, theta, beta, covariates)
{
  exp((1-(1+(time/(exp(sigma)))^exp(nu))^(1/exp(theta))) * exp(sum(covariates*beta)) )
}

Sobs <- function(time, sigma, nu, theta, beta, covariates, ratetable, age, sex, year)
{
  Sp(time = time, ratetable = ratetable, age = age, year = year, sex = sex) * 
    Sn(time = time, sigma = sigma, nu = nu, theta = theta, beta = beta, covariates = covariates)
}

Tsim <- function(sigma, nu, theta, covariates, beta, ratetable, age, sex, year)
{
  stimes <- lapply(1:length(age), function(i)
  {
    u <- runif(1)
    time <- try(uniroot(function(x) 1- Sobs(time=x, sigma=sigma, nu=nu, theta=theta, beta=beta, covariates=covariates[i,], 
                                            ratetable=ratetable, age=age[i], sex=sex[i], year=year[i])-u,
                        interval = c(.Machine$double.eps, 130000))$root, silent = TRUE) 
    return(time)
  })
}
###################################################################################################
##############
# simulation#
#############

simulate_iteration <- function(i){
  
  set.seed(i)
  
  start <- Sys.time()
  
  N <- 5000
  
  ###paramètres arrondis
  #### Hommes Colon
  TsigmaHC <- 13.0
  TnuHC <- -0.55
  TthetaHC <- 0
  betaZHC <- c(1,2.9,0.2)
  names(betaZHC) <- names(rsHC$coefficients[1:3])
  
  #### Hommes Rectum
  TsigmaHR <- 12.75
  TnuHR <- -0.3
  TthetaHR <- 0
  betaZHR <- c(1.2,3,0.35)
  names(betaZHR) <- names(rsHR$coefficients[1:3])
  
  #### Femmes Colon
  TsigmaFC <- 14 
  TnuFC <- -0.55
  TthetaFC <- 0
  betaZFC <- c(0.6,2.6,-0.1)
  names(betaZFC) <- names(rsFC$coefficients[1:3])
  
  #### Femmes Rectum
  TsigmaFR <- 11.5
  TnuFR <- -0.4
  TthetaFR <- 0
  betaZFR <- c(0.9,2.7,0.16)
  names(betaZFR) <- names(rsFR$coefficients[1:3])
  
  ######
  
  data <- data.frame(sex = 1 + rbinom(N, 1, 0.5))
  
  data$sexchara <- "male"
  data$sexchara[data$sex==2] <- "female"
  
  data$sex01 <- (data$sex)-1
  
  data$agey <- floor(rnorm(N, mean=65, sd=10))
  data$agey10 <- data$agey /10
  data$age <- data$agey*365.241
  
  # glm.colon <- glm(colon ~ sex + agey10, data = colrec, family = "binomial")
  
  lp1 <- -0.4 + 0.2 * data$sex + 0.07 * data$agey10
  
  data$colon <- rbinom(N, 1, prob =  exp(lp1 ) / (1 + exp(lp1)) )
  
  # glm.stage1 <- glm(stage1 ~ sex + agey10 + colon, data = colrec, family = "binomial")
  
  lp2 <- -1.9 + 0.1 * data$sex + 0.06 * data$agey10 - 0.7 * data$colon
  
  data$stage1 <- rbinom(N, 1, prob =  exp(lp2) / (1 + exp(lp2)) )
  
  # glm.stage2 <- glm(stage2 ~ sex + agey10 + colon, data = colrec[colrec$stage1==0,], family = "binomial")
  
  lp3 <- 0.7 + 0.2 * data$sex + 0.02 * data$agey10 - 0.2 * data$colon
  
  data$stage2 <- 0
  
  data$stage2[data$stage1==0] <- rbinom(sum(data$stage1==0), 1,
                                        prob =  exp(lp3[data$stage1==0]) / (1 + exp(lp3[data$stage1==0])) )
  
  data$stage3 <- 1 - data$stage1 - data$stage2
  
  data$stage <- 1 + 1*(data$stage2==1) + 2*(data$stage3==1)
  
  data$year <- as.numeric(sample(seq(
    as.Date('1994/01/01', origin = "1960-01-01"),
    as.Date('2000/12/30', origin = "1960-01-01"),
    by="day"), N, replace = TRUE) - as.Date("1960-01-01"))
  
  Z <- cbind(data$stage2, data$stage3, data$agey10, data$sex, data$colon)
  colnames(Z) <- c("stage2", "stage3", "agey10", "sex", "colon")
  
  ## Hommes colon
  dataHC <- data[data$sex ==1 & data$colon == 1,]
  ZHC <- cbind(dataHC$stage2, dataHC$stage3, dataHC$agey10)
  colnames(ZHC) <- c("stage2", "stage3", "agey10")
  
  data[data$sex ==1 & data$colon == 1, "timesT"] <- unlist(Tsim(sigma = TsigmaHC, nu = TnuHC, theta = TthetaHC,
                                                                beta = betaZHC, covariates = ZHC,  
                                                                ratetable = slopop,
                                                                age = dataHC$age, sex = dataHC$sexchara, year = dataHC$year))
  
  ## Hommes rectum
  dataHR <- data[data$sex ==1 & data$colon == 0,]
  ZHR <- cbind(dataHR$stage2, dataHR$stage3, dataHR$agey10)
  colnames(ZHR) <- c("stage2", "stage3", "agey10")
  
  data[data$sex ==1 & data$colon == 0, "timesT"] <- unlist(Tsim(sigma = TsigmaHR, nu = TnuHR, theta = TthetaHR,
                                                                beta = betaZHR, covariates = ZHR,  
                                                                ratetable = slopop,
                                                                age = dataHR$age, sex = dataHR$sexchara, year = dataHR$year) )
  
  ## Femmes Colon 
  dataFC <- data[data$sex ==2 & data$colon == 1,]
  ZFC <- cbind(dataFC$stage2, dataFC$stage3, dataFC$agey10)
  colnames(ZFC) <- c("stage2", "stage3", "agey10")
  
  data[data$sex == 2 & data$colon == 1, "timesT"] <- unlist(Tsim(sigma = TsigmaFC, nu = TnuFC, theta = TthetaFC,
                                                                 beta = betaZFC, covariates = ZFC,  
                                                                 ratetable = slopop,
                                                                 age = dataFC$age, sex = dataFC$sexchara, year = dataFC$year))
  
  ## Femmes Rectum 
  dataFR <- data[data$sex ==2 & data$colon == 0,]
  ZFR <- cbind(dataFR$stage2, dataFR$stage3, dataFR$agey10)
  colnames(ZFR) <- c("stage2", "stage3", "agey10")
  
  data[data$sex ==2 & data$colon == 0, "timesT"] <- unlist(Tsim(sigma = TsigmaFR, nu = TnuFR, theta = TthetaFR,
                                                                beta = betaZFR, covariates = ZFR,  
                                                                ratetable = slopop,
                                                                age = dataFR$age, sex = dataFR$sexchara, year = dataFR$year))
  
  
  data$timesT <- as.numeric(data$timesT)
  
  
  if (anyNA(data$timesT)) {
    while(anyNA(data$timesT) != FALSE ){
      na_indices <- which(is.na(data$timesT))
      
      for (p in na_indices) {
        # Hommes colon
        if(data[p,]$sex == 1 & data[p,]$colon == 1){
          Tsigma = TsigmaHC
          Tnu = TnuHC
          Ttheta = TthetaHC
          betaZ = betaZHC
        }
        # Hommes rectum
        if(data[p,]$sex == 1 & data[p,]$colon == 0){
          Tsigma = TsigmaHR
          Tnu = TnuHR
          Ttheta = TthetaHR
          betaZ = betaZHR
        }
        # Femmes colon
        if(data[p,]$sex == 2 & data[p,]$colon == 1){
          Tsigma = TsigmaFC
          Tnu = TnuFC
          Ttheta = TthetaFC
          betaZ = betaZFC
        }
        # Femmes rectum
        if(data[p,]$sex == 2 & data[p,]$colon == 0){
          Tsigma = TsigmaFR
          Tnu = TnuFR
          Ttheta = TthetaFR
          betaZ = betaZFR
        }
        data$timesT[p] <- Tsim(sigma = Tsigma, nu = Tnu, theta = Ttheta,
                               beta = betaZ, covariates = Z[p, -c(4,5), drop = FALSE],  
                               ratetable = slopop,
                               age = data$age[p], sex = data$sexchara[p], year = data$year[p]) 
      }
      data$timesT <- as.numeric(data$timesT)
      
    }
  }
  data$timesT <- as.numeric(data$timesT)
  
  #data <- na.omit(data)
  #N <- dim(data)[1] 
  ################## génération des temps de censure
  ZC <- Z[,c("agey10", "sex")] 
  
  Csim <- function(u, sigma, betaC, ZC) { -1*log(1-u)*(sigma*365.241)/exp(ZC%*%betaC) }
  #newchange
  sigmaC <- 13.5 ### 75 -> 7.5% ; 55 -> 10% ;  13.5-> 30%; 6 <- 35%, <10% d'evt au dernier temps;  5 -> ~40% ; 2.5 -> ~50% avec pas assez d'ind à risque à 10 ans 
  betaC <- c(-0.01,-0.3)
  
  data$timesC <- Csim(u = runif(n = N, min=0, max=1), sigma = sigmaC, betaC = betaC, ZC = ZC) 
  
  # boxplot(data$timesC ~ data$sex, main = "Censoring Times by Sex")
  # plot(data$age, data$timesC, main = "Censoring Times vs Age", xlab = "Age", ylab = "Censoring Time")
  
  data$times <- pmin(data$timesT, data$timesC)
  
  data$status <- 1*(data$timesT <= data$timesC)
  
  data$sex.organ <- paste0(data$sexchara, data$colon)
  data$sex.organ <- as.numeric(factor(data$sex.organ, levels = c("male0", "male1", "female0", "female1")))
  
  data$times[data$times < 10e-4] <- 10e-4
  write.table(data, paste0(path0,i,"_df",N,".csv"), sep = ';', row.names = F, col.names = T)
}


iterations <- 1:1000

mclapply(iterations, simulate_iteration, mc.cores = detectCores() - 10)
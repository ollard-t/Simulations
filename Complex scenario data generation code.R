library(survivalNET)
library(survivalPLANN)
library(parallel)
library(doParallel)


#################################################################################################
#                                        5000
#################################################################################################
### path ###

path0 <- "/home/thomas/Documents/Rstudio/Simulations/BASES_COM/"
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


rs0 <- survivalNET(formula = Surv(time, stat) ~ stage2 + stage3 + agey10 + sex + colon + 
                     ratetable(age, diag, sexchara), data = colrec,
                   ratetable=slopop, dist="weibull", weights=NULL) 

rsH <- survivalNET(formula = Surv(time, stat) ~ stage2 + stage3 + agey10 + colon + 
                     ratetable(age, diag, sexchara), data = colrec[colrec$sexchara == "male",],
                   ratetable=slopop, dist="weibull", weights=NULL) 

rsF <- survivalNET(formula = Surv(time, stat) ~ stage2 + stage3 + agey10 + colon + 
                     ratetable(age, diag, sexchara), data = colrec[colrec$sexchara == "female",],
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

simulate_iteration <- function(i,N){
  
  set.seed(i)
  
  start <- Sys.time()
  
  ###paramètres arrondis
  #### Hommes 
  TsigmaH <- 14 
  TnuH <- -1.5 
  TthetaH <- 0 
  
  #### Femmes 
  TsigmaF <- 8.5 
  TnuF <- 3.2 
  TthetaF <- 0
  
  data <- data.frame(sex = 1 + rbinom(N, 1, 0.5))
  
  data$sexchara <- "male"
  data$sexchara[data$sex==2] <- "female"
  
  data$sex01 <- (data$sex)-1
  
  data$agey <- floor(rnorm(N, mean=65, sd=10))
  data$agey10 <- data$agey /10
  data$age <- data$agey*365.241
  
  #param pour l'age (en années) c(18:85)*0.09+(0.003)*c(18:85)^2) pour après
  
  # glm.colon <- glm(colon ~ sex + agey10, data = colrec, family = "binomial")
  
  lp1 <- -0.4 + 0.2 * data$sex + 0.07 * data$agey10
  
  data$colon <- rbinom(N, 1, prob =  exp(lp1 ) / (1 + exp(lp1)) )
  
  p_colon   <- c(stage1 = 0.07, stage2 = 0.80, stage3 = 0.13)
  
  # rectum = 0
  p_rectum  <- c(stage1 = 0.32, stage2 = 0.50, stage3 = 0.18)
  
  data$stage <- NA
  
  data$stage[which(data$colon == 1)] <- sample(
    x = c(1, 2, 3),
    size = length(which(data$colon == 1)),
    replace = TRUE,
    prob = p_colon
  )
  
  data$stage[which(data$colon == 0)] <- sample(
    x = c(1, 2, 3),
    size = length(which(data$colon == 0)),
    replace = TRUE,
    prob = p_rectum
  )
  
  table(data$colon, data$stage)
  
  data$stage1 <- 1* (data$stage == 1) 
  data$stage2 <- 1* (data$stage == 2) 
  data$stage3 <- 1* (data$stage == 3) 
  
  data$year <- as.numeric(sample(seq(
    as.Date('1994/01/01', origin = "1960-01-01"),
    as.Date('2000/12/30', origin = "1960-01-01"),
    by="day"), N, replace = TRUE) - as.Date("1960-01-01"))
  
  ##interaction localisation/stade
  
  data$stage2_colon <- data$stage2 * data$colon
  data$stage3_colon <- data$stage3 * data$colon
  
  data$colon_agey10 <- data$colon * data$agey10
  
  # beta0.x = -0.4
  # beta1.x = log(2)
  # 
  # data$x7 <- rnorm(N, beta0.x + beta1.x * data$sex, 1)
  # data$x8 <- rnorm(N, beta0.x - beta1.x * data$sex - beta1.x * data$x7, 1)
  # data$x9 <- 1 * (rnorm(N, beta0.x - beta1.x * data$x8, 1) > (-0.40))
  # data$x10 <- rnorm(N, beta0.x + beta1.x * data$agey10, 1)
  # data$x11 <- 1 * (rnorm(N, beta0.x + beta1.x * data$x9, 1) > (-0.80))
  # data$x12 <- rnorm(N, beta0.x + beta1.x * data$x10, 1)
  # data$x13 <- 1 * (rnorm(N, beta0.x + beta1.x * data$x11, 1) > (0.84))
  # data$x14 <- rnorm(N, beta0.x - beta1.x * data$x12 - beta1.x * data$x11, 1)
  # data$x15 <- rnorm(N, beta0.x - beta1.x * data$x12, 1)
  # data$x16 <- 1 * (rnorm(N, 0, 1) > qnorm(0.75))
  # data$x17 <- 1 * (rnorm(N, 0, 1) > (0.66))
  # data$x18 <- rnorm(N, beta0.x - beta1.x * data$x9 - beta1.x * data$x17, 1)
  # data$x19 <- rnorm(N, beta1.x * data$x18 - beta1.x * data$x15, 1)
  # data$x20 <- 1 * (rnorm(N, beta0.x + beta1.x * data$stage3 - beta1.x*data$x16, 1) > (0.84))
  
  Z <- cbind(data$stage2, data$stage3, data$agey10, data$sex, data$colon, data$stage2_colon, data$stage3_colon, data$colon_agey10)
  # # Z <- cbind(data$stage2, data$stage3, data$agey10, data$sex, data$colon, data$stage2_colon, data$stage3_colon, data[15:28])
  colnames(Z) <- c("stage2", "stage3", "agey10", "sex", "colon", "stage2_colon", "stage3_colon", "colon_agey10")
  # # colnames(Z) <- c("stage2", "stage3", "agey10", "sex", "colon", "stage2_colon", "stage3_colon", colnames(data[15:28]))
  betaZ <- c(0.7,1.3,0.05,-0.3, 0.9, -1.2, -1.6)
  
  
  ## Hommes 
  dataH <- data[data$sex ==1,]
  ZH <- cbind(dataH$stage2, dataH$stage3, dataH$agey10, dataH$colon, dataH$stage2_colon, dataH$stage3_colon, dataH$colon_agey10)
  colnames(ZH) <- c("stage2", "stage3", "agey10", "colon", "stage2_colon", "stage3_colon", "colon_agey10")  
  # betaZH <- c(1,3,0.25, -0.1, -1.8, -1.6,0.5)
  betaZH <- c(0.7,3,0.3, 0.9, -1.2, 1.3,-0.2)
  # betaZH <- c(0,0,0, 0, 0, 0,0)
  
  data[data$sex ==1, "timesT"] <- unlist(Tsim(sigma = TsigmaH, nu = TnuH, theta = TthetaH,
                                              beta = betaZH, covariates = ZH,  
                                              ratetable = slopop,
                                              age = dataH$age, sex = dataH$sexchara, year = dataH$year) )
  
  ## Femmes 
  dataF <- data[data$sex ==2,]
  ZF <- cbind(dataF$stage2, dataF$stage3, dataF$agey10, dataF$colon, dataF$stage2_colon, dataF$stage3_colon, dataF$colon_agey10)
  colnames(ZF) <- c("stage2", "stage3", "agey10", "colon", "stage2_colon", "stage3_colon", "colon_agey10")
  betaZF <- c(0.7,2.6,0.01, 1.6, -0.5, 2,-0.01)# betaZF <- c(0,0,0, 0, 0, 0,0)
  
  data[data$sex == 2, "timesT"] <- unlist(Tsim(sigma = TsigmaF, nu = TnuF, theta = TthetaF,
                                               beta = betaZF, covariates = ZF,  
                                               ratetable = slopop,
                                               age = dataF$age, sex = dataF$sexchara, year = dataF$year))
  
  
  data$timesT <- as.numeric(data$timesT)
  
  
  if (anyNA(data$timesT)) {
    while(anyNA(data$timesT) != FALSE ){
      na_indices <- which(is.na(data$timesT))
      
      for (p in na_indices) {
        
        if(data[p,]$sex == 1){
          Tsigma = TsigmaH
          Tnu = TnuH
          Ttheta = TthetaH
          betaZ = betaZH
        }
        # Hommes rectum
        if(data[p,]$sex == 2 ){
          Tsigma = TsigmaF
          Tnu = TnuF
          Ttheta = TthetaF
          betaZ = betaZF
        }
        
        
        data$timesT[p] <- Tsim(sigma = Tsigma, nu = Tnu, theta = Ttheta,
                               beta = betaZ, covariates = Z[p,-4, drop = FALSE],  # -4 pour enelver la colonne sex dont on a besoin plus tard par contre pour la censure
                               ratetable = slopop,
                               age = data$age[p], sex = data$sexchara[p], year = data$year[p]) 
      }
      data$timesT <- as.numeric(data$timesT)
      
    }
  }
  data$timesT <- as.numeric(data$timesT)
  
  data$stage.organ <- paste0(data$stage, data$colon)
  data$stage.organ <- as.numeric(factor(data$stage.organ, levels = c("10","11","20", "21", "30", "31")))
  #data <- na.omit(data)
  #N <- dim(data)[1] 
  ################## génération des temps de censure
  dataC <- data
  
  ZC <- as.matrix(Z[,c("agey10", "sex")]) 
  
  Csim <- function(u, sigma, betaC, ZC) { -1*log(1-u)*(sigma*365.241)/exp(ZC%*%betaC) }
  #newchange
  sigmaC <- 32 ### 75 -> 7.5% ; 55 -> 10% ;  13.5-> 30%; 6 <- 35%, <10% d'evt au dernier temps;  5 -> ~40% ; 2.5 -> ~50% avec pas assez d'ind à risque à 10 ans 
  betaC <- c(-0.01,-0.3)
  
  data$timesC <- Csim(u = runif(n = N, min=0, max=1), sigma = sigmaC, betaC = betaC, ZC = ZC) 
  
  # boxplot(data$timesC ~ data$sex, main = "Censoring Times by Sex")
  # plot(data$age, data$timesC, main = "Censoring Times vs Age", xlab = "Age", ylab = "Censoring Time")
  
  data$times <- pmin(data$timesT, data$timesC)
  
  data$status <- 1*(data$timesT <= data$timesC)
  # 
  # data$sex.organ <- paste0(data$sexchara, data$colon)
  # data$sex.organ <- as.numeric(factor(data$sex.organ, levels = c("male0", "male1", "female0", "female1")))
  
  data$times[data$times < 10e-4] <- 10e-4
  
  #### High censoring
  sigmaC <- 6
  dataC$timesC <- Csim(u = runif(n = N, min=0, max=1), sigma = sigmaC, betaC = betaC, ZC = ZC) 
  
  dataC$times <- pmin(dataC$timesT, dataC$timesC)
  
  dataC$status <- 1*(dataC$timesT <= dataC$timesC)
  
  dataC$sex.organ <- paste0(dataC$sexchara, dataC$colon)
  dataC$sex.organ <- as.numeric(factor(dataC$sex.organ, levels = c("male0", "male1", "female0", "female1")))
  
  dataC$times[dataC$times < 10e-4] <- 10e-4
  # mean(data$times)
  # mean(dataC$times)
  
  write.table(data, paste0(path0,"/LC/",i,"_df",N,".csv"), sep = ';', row.names = F, col.names = T)
  write.table(dataC, paste0(path0,"/HC/",i,"_df",N,".csv"), sep = ';', row.names = F, col.names = T)
}

iterations <- 1:1000
# for(i in iterations){
#   simulate_iteration(i, N)
# }
mclapply(iterations, simulate_iteration, N = 1000, mc.cores = detectCores() - 4)

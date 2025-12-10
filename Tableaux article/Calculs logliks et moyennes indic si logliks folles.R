cens_val = "LC"
PH_val = "NPH"
strata = 1
PP = FALSE
yearoff = TRUE
percent = TRUE
N =3000

dir_path <-  paste0("~/Documents/Rstudio/Simulations/Simulations mai 2025/Résultats/N",N,"_results/",PH_val)

pattern <- paste0( N, "ind_", PH_val, "_", cens_val)

matching_files <- list.files(path = dir_path, pattern = pattern, full.names = TRUE)
matching_files <- matching_files[!grepl("PP", matching_files)]
file_name <- tools::file_path_sans_ext(basename(matching_files))

path1 <- paste0("~/Documents/Rstudio/Simulations/Simulations mai 2025/Résultats/N",N,"_results/",PH_val,"/",file_name,".Rdata")#nom à changer
load(path1)
if(PP == TRUE){
  rm(list = setdiff(ls(), c("ROCmean_results", "ROCmeanval_results", "N", "strata", "file_name", "cens_val", "PH_val", "PP", "yearoff", "percent") ) )
  path1 <- paste0("~/Documents/Rstudio/Simulations/Simulations mai 2025/Résultats/N",N,"_results/",PH_val,"/",file_name,"_PP.Rdata")
  load(path1)
}
##pour aller récupérer les logliklihood

path0 <- paste0("~/Documents/Rstudio/Simulations/Simulations mai 2025/Résultats/N",N,"_results/",PH_val,"/Output simulations_N",N,"_",PH_val,"_",cens_val,"/LOGLIK/")

method_names <- c("WG", "FLEX1.2", "FLEX1.4", "FLEX2.2", "FLEX2.4",
                  "PLANN", "PP")
method_names_val <- c("WGval", "FLEX1.2val", "FLEX1.4val", "FLEX2.2val", "FLEX2.4val",
                      "PLANNval", "PPval")
article_names <- c("&Weibull", "~ & ~ &Flexible 2","~ & ~ &Flexible 4","~ & ~ &Flexible 2 (TD)","~ & ~ &Flexible 4 (TD)", "~ & ~ &PLANN", "~ & ~ &Pohar Perme")


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

median(as.numeric(logllval[1,]))

which(as.numeric(logllval[1,]) < median(as.numeric(logllval[1,]))+0.5*median(as.numeric(logllval[1,])))

iterations <- 1:1000

path0 <- paste0("~/Documents/Rstudio/Simulations/Simulations mai 2025/Résultats/N",N,"_results/",PH_val,"/Output simulations_N",N,"_",PH_val,"_",cens_val,"/")
path1 <- "~/Documents/Rstudio/Simulations/BASES/"
rm(list = setdiff(ls(envir = .GlobalEnv), 
                  c("path0", "date_launch", "fr.ratetable", "slopop", 
                    "colrec", "Sn", "simulate_iteration", "calc_indic", "calc_indic_PP", "cens_val", "PH_val", "logllval")), envir = .GlobalEnv) ## cleaning de l'environnement excpeté ce qu iva être réutilisé
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

iterations = iterations[-which(as.numeric(logllval[1,]) < median(as.numeric(logllval[1,]))+0.5*median(as.numeric(logllval[1,]))  )]


strata_names <- c("HC", "HR", "FC", "FR")
newtimes <- c(365.241, 365.241*3, 365.241*5, 365.241*10) 


### Calculs des indicateurs

### RMSE & BIAIS
############################################ TRAIN #################################################
#WHOLE

###on récupère les valeurs théoriques moyennes Se(t) 
ALL_theo <- data.frame()

for(k in iterations){
  assign(paste0("meanTrainTheo_", k) , read.csv(paste0(path0, "TRAIN/WHOLE/THEO/mean_survival/",k,"_mean_theo.csv"),
                                                sep = " ") )
  
  ALL_theo <-rbind(ALL_theo, get(paste0("meanTrainTheo_", k)))
  rm(list =paste0("meanTrainTheo_", k))
}

#moyenne théorique de la valur de survie
S_theo <- colMeans(ALL_theo)
############################################
#importation des données des différents modèles 
#########
##plann
# ALL_PLANN <- data.frame() 
# for(k in iterations){
#   assign(paste0("meanTrainPLANN_", k) , read.csv(paste0(path0, "TRAIN/WHOLE/PLANN/mean_survival/",k,"_mean_PLANN.csv"), sep = " ") )
#   ALL_PLANN <-rbind(ALL_PLANN, get(paste0("meanTrainPLANN_", k)))
#   rm(list =paste0("meanTrainPLANN_", k))
# }

ALL_PLANN <- data.frame() 
for(k in iterations){
  assign(paste0("indTrainPLANN_", k) , read.csv(paste0(path0, "TRAIN/WHOLE/PLANN/individual_survival/",k,"_ind_PLANN.csv"), sep = ";") )
  assign(paste0("indTrainPLANN_", k), apply(get(paste0("indTrainPLANN_", k)), mean, MARGIN = 2) ) 
  ALL_PLANN<-rbind(ALL_PLANN, get(paste0("indTrainPLANN_", k)))
  rm(list =paste0("indTrainPLANN_", k))
}

###flex1
##2 noeuds
ALL_flex1.2 <- data.frame()
for(k in iterations){
  assign(paste0("meanTrainflex1.2_", k) , read.csv(paste0(path0, "TRAIN/WHOLE/FLEX1.2/mean_survival/",k,"_mean_flex12.csv"), sep = " ") )
  ALL_flex1.2 <- rbind(ALL_flex1.2, get(paste0("meanTrainflex1.2_",k)))
  rm(list =paste0("meanTrainflex1.2_", k))
}

##4 noeuds
ALL_flex1.4 <- data.frame()
for(k in iterations){
  assign(paste0("meanTrainflex1.4_", k) , read.csv(paste0(path0, "TRAIN/WHOLE/FLEX1.4/mean_survival/",k,"_mean_flex14.csv"), sep = " ") )
  ALL_flex1.4 <- rbind(ALL_flex1.4, get(paste0("meanTrainflex1.4_",k)))
  rm(list =paste0("meanTrainflex1.4_", k))
}

###flex2

##2 noeuds
ALL_flex2.2 <- data.frame()
for(k in iterations){
  assign(paste0("meanTrainflex2.2_", k) , read.csv(paste0(path0, "TRAIN/WHOLE/FLEX2.2/mean_survival/",k,"_mean_flex22.csv"), sep = " ") )
  ALL_flex2.2 <- rbind(ALL_flex2.2, get(paste0("meanTrainflex2.2_",k)))
  rm(list =paste0("meanTrainflex2.2_", k))
}

##4 noeuds
ALL_flex2.4 <- data.frame()
for(k in iterations){
  assign(paste0("meanTrainflex2.4_", k) , read.csv(paste0(path0, "TRAIN/WHOLE/FLEX2.4/mean_survival/",k,"_mean_flex24.csv"), sep = " ") )
  ALL_flex2.4 <- rbind(ALL_flex2.4, get(paste0("meanTrainflex2.4_",k)))
  rm(list =paste0("meanTrainflex2.4_", k))
}
### WG
ALL_WG <- data.frame()
for(k in iterations){
  assign(paste0("meanTrainWG_", k) , read.csv(paste0(path0, "TRAIN/WHOLE/WG/mean_survival/",k,"_mean_flexWG.csv"), sep = " ") )
  ALL_WG <- rbind(ALL_WG,get(paste0("meanTrainWG_",k))) 
  rm(list =paste0("meanTrainWG_", k))
}

ALL_PP <- data.frame() 
for(k in iterations){
  assign(paste0("meanTrainPP_", k) , read.csv(paste0(path0, "TRAIN/WHOLE/PP/",k,"PP.csv"), sep = " ") )
  ALL_PP <-rbind(ALL_PP, get(paste0("meanTrainPP_", k)))
  rm(list =paste0("meanTrainPP_", k))
}
ALL_PP <- as.data.frame(matrix(ALL_PP$x, ncol = 4, byrow = TRUE))

biais_PLANN <- colMeans(ALL_PLANN[,-1] - S_theo[col( ALL_PLANN[,-1] )] )
RMSE_PLANN <- sqrt( colMeans((ALL_PLANN[,-1] - S_theo[col( ALL_PLANN[,-1] )])^2) )

biais_flex1.2 <- colMeans( ALL_flex1.2 - S_theo[col( ALL_flex1.2)] )
RMSE_flex1.2 <- sqrt( colMeans( (ALL_flex1.2 - S_theo[col( ALL_flex1.2)])^2 ) )

biais_flex1.4 <- colMeans( ALL_flex1.4 - S_theo[col( ALL_flex1.4)] )
RMSE_flex1.4 <- sqrt( colMeans( (ALL_flex1.4 - S_theo[col( ALL_flex1.4)])^2 ) )

biais_flex2.2 <- colMeans( ALL_flex2.2 - S_theo[col( ALL_flex2.2)] )
RMSE_flex2.2 <- sqrt( colMeans( (ALL_flex2.2 - S_theo[col( ALL_flex2.2)])^2 ) )

biais_flex2.4 <- colMeans( ALL_flex2.4 - S_theo[col( ALL_flex2.4)] )
RMSE_flex2.4 <- sqrt( colMeans( (ALL_flex2.4 - S_theo[col( ALL_flex2.4)])^2 ) )

biais_WG <- colMeans( ALL_WG - S_theo[col( ALL_WG)] )
RMSE_WG <- sqrt( colMeans( (ALL_WG - S_theo[col( ALL_WG)])^2 ) )

biais_PP <- colMeans( ALL_PP - S_theo[col( ALL_PP)] )
RMSE_PP <- sqrt( colMeans( (ALL_PP - S_theo[col( ALL_PP)])^2 ) )
#STRATES

###on récupère les valeurs théoriques moyennes Se(t) 
for(j in strata_names){
  assign(paste0("ALL_",j, "_theo"),  data.frame()) 
}

for(k in iterations){
  assign(paste0("meanTrainTheoS_", k) , read.csv(paste0(path0, "TRAIN/STRATA/THEO/mean_survival/",k,"mean_theo_strates.csv"),
                                                 sep = " ", row.names = 2)[,-1] ) 
  for(j in strata_names){
    
    current_row <- get(paste0("meanTrainTheoS_", k))[paste0("mean.theo_", j), , drop = FALSE]
    rownames(current_row) <- paste0("mean.theo_", j, k)
    assign(paste0("ALL_", j, "_theo"), rbind(get(paste0("ALL_", j, "_theo")), current_row))
    
  }
  
  rm(list =paste0("meanTrainTheoS_", k))
}

#moyenne théorique de la valeur de survie
for(j in strata_names){
  assign(paste0("S_theo", j ),colMeans(get(paste0("ALL_",j,"_theo") ) ) )
}
############################################
#importation des données des différents modèles 
############################################

# PLANN
# 
for(j in strata_names){
  assign(paste0("ALL_PLANN_", j),data.frame())
}

# for(k in iterations){
#   assign(paste0("meanTrainPLANNstr_", k) , read.csv(paste0(path0, "TRAIN/STRATA/PLANN/mean_survival/",k,"mean_plann_strates.csv"), sep = " ", row.names = 1) )
# 
#   for(j in strata_names){
#     assign(paste0("ALL_PLANN_",j) ,rbind(get(paste0("ALL_PLANN_",j)),get(paste0("meanTrainPLANNstr_",k))[paste0("mean.plann_", j),-1]  ) )
#   }
#   rm(list =paste0("meanTrainPLANNstr_", k))
# }

for(k in iterations){
  HC <- read.csv(paste0(path0, "TRAIN/STRATA/PLANN/individual_survival/",k,"_ind_PLANN_HC.csv"), sep = ";")
  HC <- apply(HC, mean, MARGIN = 2)[-1]
  
  HR <- read.csv(paste0(path0, "TRAIN/STRATA/PLANN/individual_survival/",k,"_ind_PLANN_HR.csv"), sep = ";")
  HR <- apply(HR, mean, MARGIN = 2)[-1]
  
  FC <- read.csv(paste0(path0, "TRAIN/STRATA/PLANN/individual_survival/",k,"_ind_PLANN_FC.csv"), sep = ";")
  FC <- apply(FC, mean, MARGIN = 2)[-1]
  
  FR <- read.csv(paste0(path0, "TRAIN/STRATA/PLANN/individual_survival/",k,"_ind_PLANN_FR.csv"), sep = ";")
  FR <- apply(FR, mean, MARGIN = 2)[-1]
  
  for(j in strata_names){
    assign(paste0("ALL_PLANN_",j) ,rbind(get(paste0("ALL_PLANN_",j)), get(j) ) )
  }
  rm(list =c("HC", "HR", "FC", "FR"))
}

### flex1

##2 noeuds
for(j in strata_names){
  assign(paste0("ALL_FLEX1.2_", j),data.frame()) 
}

for(k in iterations){
  assign(paste0("meanTrainFLEX1.2str_", k) , read.csv(paste0(path0, "TRAIN/STRATA/FLEX1.2/mean_survival/",k,"mean_flex12_strates.csv"), sep = " ", row.names = 2)[,-1] )
  
  for(j in strata_names){
    assign(paste0("ALL_FLEX1.2_",j) ,rbind(get(paste0("ALL_FLEX1.2_",j)), get(paste0("meanTrainFLEX1.2str_",k))[paste0("mean.flex1.2_", j),]  ) )
  }
  rm(list =paste0("meanTrainFLEX1.2str_", k))
}

##4 noeuds
for(j in strata_names){
  assign(paste0("ALL_FLEX1.4_", j),data.frame()) 
}

for(k in iterations){
  assign(paste0("meanTrainFLEX1.4str_", k) , read.csv(paste0(path0, "TRAIN/STRATA/FLEX1.4/mean_survival/",k,"mean_flex14_strates.csv"), sep = " ", row.names = 2)[,-1] )
  
  for(j in strata_names){
    assign(paste0("ALL_FLEX1.4_",j) ,rbind(get(paste0("ALL_FLEX1.4_",j)), get(paste0("meanTrainFLEX1.4str_",k))[paste0("mean.flex1.4_", j),]  ) )
  }
  rm(list =paste0("meanTrainFLEX1.4str_", k))
}

### flex2
##2 noeuds
for(j in strata_names){
  assign(paste0("ALL_FLEX2.2_", j),data.frame()) 
}

for(k in iterations){
  assign(paste0("meanTrainFLEX2.2str_", k) , read.csv(paste0(path0, "TRAIN/STRATA/FLEX2.2/mean_survival/",k,"mean_flex22_strates.csv"), sep = " ", row.names = 2)[,-1] )
  
  for(j in strata_names){
    assign(paste0("ALL_FLEX2.2_",j) ,rbind(get(paste0("ALL_FLEX2.2_",j)), get(paste0("meanTrainFLEX2.2str_",k))[paste0("mean.flex2.2_", j),]  ) )
  }
  rm(list =paste0("meanTrainFLEX2.2str_", k))
}

##4 noeuds
for(j in strata_names){
  assign(paste0("ALL_FLEX2.4_", j),data.frame()) 
}

for(k in iterations){
  assign(paste0("meanTrainFLEX2.4str_", k) , read.csv(paste0(path0, "TRAIN/STRATA/FLEX2.4/mean_survival/",k,"mean_flex24_strates.csv"), sep = " ", row.names = 2)[,-1] )
  
  for(j in strata_names){
    assign(paste0("ALL_FLEX2.4_",j) ,rbind(get(paste0("ALL_FLEX2.4_",j)), get(paste0("meanTrainFLEX2.4str_",k))[paste0("mean.flex2.4_", j),]  ) )
  }
  rm(list =paste0("meanTrainFLEX2.4str_", k))
}
### WG

for(j in strata_names){
  assign(paste0("ALL_WG_", j),data.frame()) 
}

for(k in iterations){
  assign(paste0("meanTrainWGstr_", k) , read.csv(paste0(path0, "TRAIN/STRATA/WG/mean_survival/",k,"mean_WG_strates.csv"), sep = " ", row.names = 2)[,-1] )
  
  for(j in strata_names){
    assign(paste0("ALL_WG_",j) ,rbind(get(paste0("ALL_WG_",j)), get(paste0("meanTrainWGstr_",k))[paste0("mean.WG_", j),]  ) )
  }
  rm(list =paste0("meanTrainWGstr_", k))
}

### PP

for(j in strata_names){
  assign(paste0("ALL_PP_", j),data.frame()) 
}

for(k in iterations){
  assign(paste0("meanTrainPPstr_", k) , read.csv(paste0(path0, "TRAIN/STRATA/PP/",k,"PP_strates.csv"), sep = " ", row.names = 2)[,-1] )
  
  for(j in strata_names){
    assign(paste0("ALL_PP_",j) ,rbind(get(paste0("ALL_PP_",j)), get(paste0("meanTrainPPstr_",k))[paste0("poharpred_", j),]  ) )
  }
  rm(list =paste0("meanTrainPPstr_", k))
}

#### Calcul des indicateurs
##PLANN
biais_PLANN_strates <- data.frame()
RMSE_PLANN_strates <- data.frame()
for(j in strata_names){
  assign(paste0("biais_PLANN_", j), 
         colMeans( get( paste0("ALL_PLANN_", j) ) - get(paste0("S_theo", j) )[col(get( paste0("ALL_PLANN_", j)) )] ) )
  #on concatene tous les résultats dans une dataframe pour ne pas saturer l'evt 
  biais_PLANN_strates <- rbind(biais_PLANN_strates, get(paste0("biais_PLANN_", j)) )
  rownames(biais_PLANN_strates)[dim(biais_PLANN_strates)[1]] <- paste0(j)
  rm(list = paste0("biais_PLANN_", j))
  
  assign(paste0("RMSE_PLANN_",  j) ,
         sqrt( colMeans( ( get( paste0("ALL_PLANN_", j) ) - 
                             get(paste0("S_theo", j) )[col(get( paste0("ALL_PLANN_", j)) )] )^2 ) ) )
  RMSE_PLANN_strates <- rbind(RMSE_PLANN_strates, get(paste0("RMSE_PLANN_", j)) )
  rownames(RMSE_PLANN_strates)[dim(RMSE_PLANN_strates)[1]] <- paste0(j)
  rm(list = paste0("RMSE_PLANN_", j))
}

colnames(biais_PLANN_strates) <- colnames(ALL_flex1.4)
colnames(RMSE_PLANN_strates) <- colnames(ALL_flex1.4)

# FLEX 1
#2 noeuds
biais_FLEX1.2_strates <- data.frame()
RMSE_FLEX1.2_strates <- data.frame()
for(j in strata_names){
  assign(paste0("biais_FLEX1.2_", j), 
         colMeans( get( paste0("ALL_FLEX1.2_", j) ) - get(paste0("S_theo", j) )[col(get( paste0("ALL_FLEX1.2_", j)) )] ) )
  
  #on concatene tous les résultats dans une dataframe pour ne pas saturer l'evt 
  biais_FLEX1.2_strates <- rbind(biais_FLEX1.2_strates, get(paste0("biais_FLEX1.2_", j)) )
  rownames(biais_FLEX1.2_strates)[dim(biais_FLEX1.2_strates)[1]] <- paste0(j)
  rm(list = paste0("biais_FLEX1.2_", j))
  
  assign(paste0("RMSE_FLEX1.2_",  j) ,
         sqrt( colMeans( ( get( paste0("ALL_FLEX1.2_", j) ) - 
                             get(paste0("S_theo", j) )[col(get( paste0("ALL_FLEX1.2_", j)) )] )^2 ) ) )
  
  RMSE_FLEX1.2_strates <- rbind(RMSE_FLEX1.2_strates, get(paste0("RMSE_FLEX1.2_", j)) )
  rownames(RMSE_FLEX1.2_strates)[dim(RMSE_FLEX1.2_strates)[1]] <- paste0(j)
  rm(list = paste0("RMSE_FLEX1.2_", j))
  
}

colnames(biais_FLEX1.2_strates) <- colnames(ALL_flex1.2)
colnames(RMSE_FLEX1.2_strates) <- colnames(ALL_flex1.2)

#4 noeuds
biais_FLEX1.4_strates <- data.frame()
RMSE_FLEX1.4_strates <- data.frame()
for(j in strata_names){
  assign(paste0("biais_FLEX1.4_", j), 
         colMeans( get( paste0("ALL_FLEX1.4_", j) ) - get(paste0("S_theo", j) )[col(get( paste0("ALL_FLEX1.4_", j)) )] ) )
  
  #on concatene tous les résultats dans une dataframe pour ne pas saturer l'evt 
  biais_FLEX1.4_strates <- rbind(biais_FLEX1.4_strates, get(paste0("biais_FLEX1.4_", j)) )
  rownames(biais_FLEX1.4_strates)[dim(biais_FLEX1.4_strates)[1]] <- paste0(j)
  rm(list = paste0("biais_FLEX1.4_", j))
  
  assign(paste0("RMSE_FLEX1.4_",  j) ,
         sqrt( colMeans( ( get( paste0("ALL_FLEX1.4_", j) ) - 
                             get(paste0("S_theo", j) )[col(get( paste0("ALL_FLEX1.4_", j)) )] )^2 ) ) )
  
  RMSE_FLEX1.4_strates <- rbind(RMSE_FLEX1.4_strates, get(paste0("RMSE_FLEX1.4_", j)) )
  rownames(RMSE_FLEX1.4_strates)[dim(RMSE_FLEX1.4_strates)[1]] <- paste0(j)
  rm(list = paste0("RMSE_FLEX1.4_", j))
  
}

colnames(biais_FLEX1.4_strates) <- colnames(ALL_flex1.4)
colnames(RMSE_FLEX1.4_strates) <- colnames(ALL_flex1.4)


# FLEX 2
##2 noeuds
biais_FLEX2.2_strates <- data.frame()
RMSE_FLEX2.2_strates <- data.frame()
for(j in strata_names){
  assign(paste0("biais_FLEX2.2_", j), 
         colMeans( get( paste0("ALL_FLEX2.2_", j) ) - get(paste0("S_theo", j) )[col(get( paste0("ALL_FLEX2.2_", j)) )] ) )
  
  biais_FLEX2.2_strates <- rbind(biais_FLEX2.2_strates, get(paste0("biais_FLEX2.2_", j)) )
  rownames(biais_FLEX2.2_strates)[dim(biais_FLEX2.2_strates)[1]] <- paste0(j)
  rm(list = paste0("biais_FLEX2.2_", j))
  
  assign(paste0("RMSE_FLEX2.2_",  j) ,
         sqrt( colMeans( ( get( paste0("ALL_FLEX2.2_", j) ) - 
                             get(paste0("S_theo", j) )[col(get( paste0("ALL_FLEX2.2_", j)) )] )^2 ) ) )
  
  RMSE_FLEX2.2_strates <- rbind(RMSE_FLEX2.2_strates, get(paste0("RMSE_FLEX2.2_", j)) )
  rownames(RMSE_FLEX2.2_strates)[dim(RMSE_FLEX2.2_strates)[1]] <- paste0(j)
  rm(list = paste0("RMSE_FLEX2.2_", j))
}
colnames(biais_FLEX2.2_strates) <- colnames(ALL_flex2.2)
colnames(RMSE_FLEX2.2_strates) <- colnames(ALL_flex2.2)

##4 noeuds
biais_FLEX2.4_strates <- data.frame()
RMSE_FLEX2.4_strates <- data.frame()
for(j in strata_names){
  assign(paste0("biais_FLEX2.4_", j), 
         colMeans( get( paste0("ALL_FLEX2.4_", j) ) - get(paste0("S_theo", j) )[col(get( paste0("ALL_FLEX2.4_", j)) )] ) )
  
  biais_FLEX2.4_strates <- rbind(biais_FLEX2.4_strates, get(paste0("biais_FLEX2.4_", j)) )
  rownames(biais_FLEX2.4_strates)[dim(biais_FLEX2.4_strates)[1]] <- paste0(j)
  rm(list = paste0("biais_FLEX2.4_", j))
  
  assign(paste0("RMSE_FLEX2.4_",  j) ,
         sqrt( colMeans( ( get( paste0("ALL_FLEX2.4_", j) ) - 
                             get(paste0("S_theo", j) )[col(get( paste0("ALL_FLEX2.4_", j)) )] )^2 ) ) )
  
  RMSE_FLEX2.4_strates <- rbind(RMSE_FLEX2.4_strates, get(paste0("RMSE_FLEX2.4_", j)) )
  rownames(RMSE_FLEX2.4_strates)[dim(RMSE_FLEX2.4_strates)[1]] <- paste0(j)
  rm(list = paste0("RMSE_FLEX2.4_", j))
}
colnames(biais_FLEX2.4_strates) <- colnames(ALL_flex2.4)
colnames(RMSE_FLEX2.4_strates) <- colnames(ALL_flex2.4)

# WG 
biais_WG_strates <- data.frame()
RMSE_WG_strates <- data.frame()
for(j in strata_names){
  assign(paste0("biais_WG_", j), 
         colMeans( get( paste0("ALL_WG_", j) ) - get(paste0("S_theo", j) )[col(get( paste0("ALL_WG_", j)) )] ) )
  
  biais_WG_strates <- rbind(biais_WG_strates, get(paste0("biais_WG_", j)) )
  rownames(biais_WG_strates)[dim(biais_WG_strates)[1]] <- paste0(j)
  rm(list = paste0("biais_WG_", j))
  
  assign(paste0("RMSE_WG_",  j) ,
         sqrt( colMeans( ( get( paste0("ALL_WG_", j) ) - 
                             get(paste0("S_theo", j) )[col(get( paste0("ALL_WG_", j)) )] )^2 ) ) )
  
  RMSE_WG_strates <- rbind(RMSE_WG_strates, get(paste0("RMSE_WG_", j)) )
  rownames(RMSE_WG_strates)[dim(RMSE_WG_strates)[1]] <- paste0(j)
  rm(list = paste0("RMSE_WG_", j))
}

colnames(biais_WG_strates) <- colnames(ALL_WG)
colnames(RMSE_WG_strates) <- colnames(ALL_WG)

# PP
biais_PP_strates <- data.frame()
RMSE_PP_strates <- data.frame()
for(j in strata_names){
  assign(paste0("biais_PP_", j), 
         colMeans( get( paste0("ALL_PP_", j) ) - get(paste0("S_theo", j) )[col(get( paste0("ALL_PP_", j)) )] ) )
  
  biais_PP_strates <- rbind(biais_PP_strates, get(paste0("biais_PP_", j)) )
  rownames(biais_PP_strates)[dim(biais_PP_strates)[1]] <- paste0(j)
  rm(list = paste0("biais_PP_", j))
  
  assign(paste0("RMSE_PP_",  j) ,
         sqrt( colMeans( ( get( paste0("ALL_PP_", j) ) - 
                             get(paste0("S_theo", j) )[col(get( paste0("ALL_PP_", j)) )] )^2 ) ) )
  
  RMSE_PP_strates <- rbind(RMSE_PP_strates, get(paste0("RMSE_PP_", j)) )
  rownames(RMSE_PP_strates)[dim(RMSE_PP_strates)[1]] <- paste0(j)
  rm(list = paste0("RMSE_PP_", j))
}

colnames(biais_PP_strates) <- colnames(ALL_PP)
colnames(RMSE_PP_strates) <- colnames(ALL_PP)

###concaten des S_theo, là aussi pour libérer de la place dans l'evt global
S_theo_strates <- data.frame()
for(j in strata_names){
  S_theo_strates <- rbind(S_theo_strates, get(paste0("S_theo",j)))
  rownames(S_theo_strates)[dim(S_theo_strates)[1]] <- paste0(j)
  rm(list = paste0("S_theo",j))
}
colnames(S_theo_strates) <- colnames(ALL_WG)

############################################ VALID #################################################

# WHOLE
###on récupère les valeurs théoriques moyennes Se(t) 
ALL_theoval <- data.frame()

for(k in iterations){
  assign(paste0("meanValidTheo_", k) , read.csv(paste0(path0, "VALID/WHOLE/THEO/mean_survival/",k,"_mean_theoval.csv"),
                                                sep = " ") )
  
  ALL_theoval <-rbind(ALL_theoval, get(paste0("meanValidTheo_", k)))
  rm(list =paste0("meanValidTheo_", k))
}

#moyenne théorique de la valur de survie
S_theoval <- colMeans(ALL_theoval)

############################################
#importation des données des différents modèles 
#########
##plann
# ALL_PLANNval <- data.frame() 
# for(k in iterations){
#   assign(paste0("meanValidPLANN_", k) , read.csv(paste0(path0, "VALID/WHOLE/PLANN/mean_survival/",k,"_mean_PLANNval.csv"), sep = " ") )
#   ALL_PLANNval <-rbind(ALL_PLANNval, get(paste0("meanValidPLANN_", k)))
#   rm(list =paste0("meanValidPLANN_", k))
# }

ALL_PLANNval <- data.frame() 
for(k in iterations){
  assign(paste0("indValidPLANN_", k) , read.csv(paste0(path0, "VALID/WHOLE/PLANN/individual_survival/",k,"_ind_PLANNval.csv"), sep = ";") )
  assign(paste0("indValidPLANN_", k), apply(get(paste0("indValidPLANN_", k)), mean, MARGIN = 2) ) 
  ALL_PLANNval <-rbind(ALL_PLANNval, get(paste0("indValidPLANN_", k)))
  rm(list =paste0("indValidPLANN_", k))
}
###flex1
##2 noeuds
ALL_flex1.2val <- data.frame()
for(k in iterations){
  assign(paste0("meanValidflex1.2_", k) , read.csv(paste0(path0, "VALID/WHOLE/FLEX1.2/mean_survival/",k,"_mean_flexval12.csv"), sep = " ") )
  ALL_flex1.2val <- rbind(ALL_flex1.2val, get(paste0("meanValidflex1.2_",k)))
  rm(list =paste0("meanValidflex1.2_", k))
}
##4 noeuds
ALL_flex1.4val <- data.frame()
for(k in iterations){
  assign(paste0("meanValidflex1.4_", k) , read.csv(paste0(path0, "VALID/WHOLE/FLEX1.4/mean_survival/",k,"_mean_flexval14.csv"), sep = " ") )
  ALL_flex1.4val <- rbind(ALL_flex1.4val, get(paste0("meanValidflex1.4_",k)))
  rm(list =paste0("meanValidflex1.4_", k))
}
### flex2

##2 noeuds
ALL_flex2.2val <- data.frame()
for(k in iterations){
  assign(paste0("meanValidflex2.2_", k) , read.csv(paste0(path0, "VALID/WHOLE/FLEX2.2/mean_survival/",k,"_mean_flexval22.csv"), sep = " ") )
  ALL_flex2.2val <- rbind(ALL_flex2.2val, get(paste0("meanValidflex2.2_",k)))
  rm(list =paste0("meanValidflex2.2_", k))
}

##4 noeuds
ALL_flex2.4val <- data.frame()
for(k in iterations){
  assign(paste0("meanValidflex2.4_", k) , read.csv(paste0(path0, "VALID/WHOLE/FLEX2.4/mean_survival/",k,"_mean_flexval24.csv"), sep = " ") )
  ALL_flex2.4val <- rbind(ALL_flex2.4val, get(paste0("meanValidflex2.4_",k)))
  rm(list =paste0("meanValidflex2.4_", k))
}
### WG
ALL_WGval <- data.frame()
for(k in iterations){
  assign(paste0("meanValidWG_", k) , read.csv(paste0(path0, "VALID/WHOLE/WG/mean_survival/",k,"_mean_flexWGval.csv"), sep = " ") )
  ALL_WGval <- rbind(ALL_WGval, get(paste0("meanValidWG_",k))) 
  rm(list =paste0("meanValidWG_", k))
}

##PP
ALL_PPval <- data.frame()
for(k in iterations){
  assign(paste0("meanValidPP_", k) , read.csv(paste0(path0, "VALID/WHOLE/PP/",k,"PPval.csv"), sep = " ") )
  ALL_PPval <- rbind(ALL_PPval, get(paste0("meanValidPP_",k))) 
  rm(list =paste0("meanValidPP_", k))
}
ALL_PPval <- as.data.frame(matrix(ALL_PPval$x, ncol = 4, byrow = TRUE))

biais_PLANNval <- colMeans(ALL_PLANNval[,-1] - S_theoval[col( ALL_PLANNval[,-1] )] )
RMSE_PLANNval <- sqrt( colMeans((ALL_PLANNval[,-1] - S_theoval[col( ALL_PLANNval[,-1] )])^2) )

biais_flex1.2val <- colMeans( ALL_flex1.2val - S_theoval[col( ALL_flex1.2val)] )
RMSE_flex1.2val <- sqrt( colMeans( (ALL_flex1.2val - S_theoval[col( ALL_flex1.2val)])^2 ) )

biais_flex1.4val <- colMeans( ALL_flex1.4val - S_theoval[col( ALL_flex1.4val)] )
RMSE_flex1.4val <- sqrt( colMeans( (ALL_flex1.4val - S_theoval[col( ALL_flex1.4val)])^2 ) )

biais_flex2.2val <- colMeans( ALL_flex2.2val - S_theoval[col( ALL_flex2.2val)] )
RMSE_flex2.2val <- sqrt( colMeans( (ALL_flex2.2val - S_theoval[col( ALL_flex2.2val)])^2 ) )

biais_flex2.4val <- colMeans( ALL_flex2.4val - S_theoval[col( ALL_flex2.4val)] )
RMSE_flex2.4val <- sqrt( colMeans( (ALL_flex2.4val - S_theoval[col( ALL_flex2.4val)])^2 ) )

biais_WGval <- colMeans( ALL_WGval - S_theoval[col( ALL_WGval)] )
RMSE_WGval <- sqrt( colMeans( (ALL_WGval - S_theoval[col( ALL_WGval)])^2 ) )

biais_PPval <- colMeans( ALL_PPval - S_theoval[col( ALL_PPval)] )
RMSE_PPval <- sqrt( colMeans( (ALL_PPval - S_theoval[col( ALL_PPval)])^2 ) )
#STRATES

###on récupère les valeurs théoriques moyennes Se(t) 
for(j in strata_names){
  assign(paste0("ALL_",j, "_theoval"),  data.frame()) 
}

for(k in iterations){
  assign(paste0("meanValidTheoS_", k) , read.csv(paste0(path0, "VALID/STRATA/THEO/mean_survival/",k,"mean_theoval_strates.csv"),
                                                 sep = " ", row.names = 2)[,-c(1,2)] ) ### [,-1] car quand j'ai enregistré les résultats, j'ai dupliué la premiere colonne donc il faut l'enlever
  for(j in strata_names){
    
    current_row <- get(paste0("meanValidTheoS_", k))[paste0("mean.theo_", j), , drop = FALSE]
    rownames(current_row) <- paste0("mean.theoval_", j, k)
    assign(paste0("ALL_", j, "_theoval"), rbind(get(paste0("ALL_", j, "_theoval")), current_row))
    
  }
  
  rm(list =paste0("meanValidTheoS_", k))
}
rm(current_row)
#moyenne théorique de la valeur de survie
for(j in strata_names){
  assign(paste0("S_theoval", j ),colMeans(get(paste0("ALL_",j,"_theoval") ) ) )
}
############################################
#importation des données des différents modèles 
############################################

# PLANN

for(j in strata_names){
  assign(paste0("ALL_PLANNval_", j),data.frame()) 
}

# for(k in iterations){
#   assign(paste0("meanValidPLANNstr_", k) , read.csv(paste0(path0, "VALID/STRATA/PLANN/mean_survival/",k,"mean_plannval_strates.csv"), sep = " ", row.names = 1) )
#   
#   for(j in strata_names){
#     assign(paste0("ALL_PLANNval_",j) ,rbind(get(paste0("ALL_PLANNval_",j)),get(paste0("meanValidPLANNstr_",k))[paste0("mean.plannval_", j),-1]  ) )
#   }
#   rm(list =paste0("meanValidPLANNstr_", k))
# }

for(k in iterations){
  HC <- read.csv(paste0(path0, "VALID/STRATA/PLANN/individual_survival/",k,"_ind_PLANNval_HC.csv"), sep = ";")
  HC <- apply(HC, mean, MARGIN = 2)[-1]
  
  HR <- read.csv(paste0(path0, "VALID/STRATA/PLANN/individual_survival/",k,"_ind_PLANNval_HR.csv"), sep = ";")
  HR <- apply(HR, mean, MARGIN = 2)[-1]
  
  FC <- read.csv(paste0(path0, "VALID/STRATA/PLANN/individual_survival/",k,"_ind_PLANNval_FC.csv"), sep = ";")
  FC <- apply(FC, mean, MARGIN = 2)[-1]
  
  FR <- read.csv(paste0(path0, "VALID/STRATA/PLANN/individual_survival/",k,"_ind_PLANNval_FR.csv"), sep = ";")
  FR <- apply(FR, mean, MARGIN = 2)[-1]
  
  for(j in strata_names){
    assign(paste0("ALL_PLANNval_",j) ,rbind(get(paste0("ALL_PLANNval_",j)), get(j) ) )
  }
  rm(list =c("HC", "HR", "FC", "FR"))
}

### flex1
##2 noeuds
for(j in strata_names){
  assign(paste0("ALL_FLEX1.2val_", j),data.frame()) 
}

for(k in iterations){
  assign(paste0("meanValidFLEX1.2str_", k) , read.csv(paste0(path0, "VALID/STRATA/FLEX1.2/mean_survival/",k,"mean_flexval12_strates.csv"), sep = " ", row.names = 2)[,-1] )
  
  for(j in strata_names){
    assign(paste0("ALL_FLEX1.2val_",j) ,rbind(get(paste0("ALL_FLEX1.2val_",j)), get(paste0("meanValidFLEX1.2str_",k))[paste0("mean.flexval1.2_", j),]  ) )
  }
  rm(list =paste0("meanValidFLEX1.2str_", k))
}

##4 noeuds
for(j in strata_names){
  assign(paste0("ALL_FLEX1.4val_", j),data.frame()) 
}

for(k in iterations){
  assign(paste0("meanValidFLEX1.4str_", k) , read.csv(paste0(path0, "VALID/STRATA/FLEX1.4/mean_survival/",k,"mean_flexval14_strates.csv"), sep = " ", row.names = 2)[,-1] )
  
  for(j in strata_names){
    assign(paste0("ALL_FLEX1.4val_",j) ,rbind(get(paste0("ALL_FLEX1.4val_",j)), get(paste0("meanValidFLEX1.4str_",k))[paste0("mean.flexval1.4_", j),]  ) )
  }
  rm(list =paste0("meanValidFLEX1.4str_", k))
}
### flex2

##2 noeuds
for(j in strata_names){
  assign(paste0("ALL_FLEX2.2val_", j),data.frame()) 
}

for(k in iterations){
  assign(paste0("meanValidFLEX2.2str_", k) , read.csv(paste0(path0, "VALID/STRATA/FLEX2.2/mean_survival/",k,"mean_flexval22_strates.csv"), sep = " ", row.names = 2)[,-1] )
  
  for(j in strata_names){
    assign(paste0("ALL_FLEX2.2val_",j) ,rbind(get(paste0("ALL_FLEX2.2val_",j)), get(paste0("meanValidFLEX2.2str_",k))[paste0("mean.flexval2.2_", j),]  ) )
  }
  rm(list =paste0("meanValidFLEX2.2str_", k))
}

##4 noeuds
for(j in strata_names){
  assign(paste0("ALL_FLEX2.4val_", j),data.frame()) 
}

for(k in iterations){
  assign(paste0("meanValidFLEX2.4str_", k) , read.csv(paste0(path0, "VALID/STRATA/FLEX2.4/mean_survival/",k,"mean_flexval24_strates.csv"), sep = " ", row.names = 2)[,-1] )
  
  for(j in strata_names){
    assign(paste0("ALL_FLEX2.4val_",j) ,rbind(get(paste0("ALL_FLEX2.4val_",j)), get(paste0("meanValidFLEX2.4str_",k))[paste0("mean.flexval2.4_", j),]  ) )
  }
  rm(list =paste0("meanValidFLEX2.4str_", k))
}
### WG

for(j in strata_names){
  assign(paste0("ALL_WGval_", j),data.frame()) 
}

for(k in iterations){
  assign(paste0("meanValidWGstr_", k) , read.csv(paste0(path0, "VALID/STRATA/WG/mean_survival/",k,"mean_WGval_strates.csv"), sep = " ", row.names = 2)[,-1] )
  
  for(j in strata_names){
    assign(paste0("ALL_WGval_",j) ,rbind(get(paste0("ALL_WGval_",j)), get(paste0("meanValidWGstr_",k))[paste0("mean.WGval_", j),]  ) )
  }
  rm(list =paste0("meanValidWGstr_", k))
}

### PP 

for(j in strata_names){
  assign(paste0("ALL_PPval_", j),data.frame()) 
}

for(k in iterations){
  assign(paste0("meanValidPPstr_", k) , read.csv(paste0(path0, "VALID/STRATA/PP/",k,"PPval_strates.csv"), sep = " ", row.names = 2)[,-1] )
  
  for(j in strata_names){
    assign(paste0("ALL_PPval_",j) ,rbind(get(paste0("ALL_PPval_",j)), get(paste0("meanValidPPstr_",k))[paste0("poharpredval_", j),]  ) )
  }
  rm(list =paste0("meanValidPPstr_", k))
}
#### Calcul des indicateurs
##PLANN
biais_PLANNval_strates <- data.frame()
RMSE_PLANNval_strates <- data.frame()
for(j in strata_names){
  assign(paste0("biais_PLANNval_", j), 
         colMeans( get( paste0("ALL_PLANNval_", j) ) - get(paste0("S_theoval", j) )[col(get( paste0("ALL_PLANNval_", j)) )] ) )
  
  biais_PLANNval_strates <- rbind(biais_PLANNval_strates, get(paste0("biais_PLANNval_", j)) )
  rownames(biais_PLANNval_strates)[dim(biais_PLANNval_strates)[1]] <- paste0(j)
  rm(list = paste0("biais_PLANNval_", j))
  
  assign(paste0("RMSE_PLANNval_",  j) ,
         sqrt( colMeans( ( get( paste0("ALL_PLANNval_", j) ) - 
                             get(paste0("S_theoval", j) )[col(get( paste0("ALL_PLANNval_", j)) )] )^2 ) ) )
  
  RMSE_PLANNval_strates <- rbind(RMSE_PLANNval_strates, get(paste0("RMSE_PLANNval_", j)) )
  rownames(RMSE_PLANNval_strates)[dim(RMSE_PLANNval_strates)[1]] <- paste0(j)
  rm(list = paste0("RMSE_PLANNval_", j))
  
}

colnames(biais_PLANNval_strates) <- colnames(ALL_flex1.2)
colnames(RMSE_PLANNval_strates) <- colnames(ALL_flex1.2)

# FLEX 1

##2 noeuds
biais_FLEX1.2val_strates <- data.frame()
RMSE_FLEX1.2val_strates <- data.frame()
for(j in strata_names){
  assign(paste0("biais_FLEX1.2val_", j), 
         colMeans( get( paste0("ALL_FLEX1.2val_", j) ) - get(paste0("S_theoval", j) )[col(get( paste0("ALL_FLEX1.2val_", j)) )] ) )
  
  biais_FLEX1.2val_strates <- rbind(biais_FLEX1.2val_strates, get(paste0("biais_FLEX1.2val_", j)) )
  rownames(biais_FLEX1.2val_strates)[dim(biais_FLEX1.2val_strates)[1]] <- paste0(j)
  rm(list = paste0("biais_FLEX1.2val_", j))
  
  assign(paste0("RMSE_FLEX1.2val_",  j) ,
         sqrt( colMeans( ( get( paste0("ALL_FLEX1.2val_", j) ) - 
                             get(paste0("S_theoval", j) )[col(get( paste0("ALL_FLEX1.2val_", j)) )] )^2 ) ) )
  
  RMSE_FLEX1.2val_strates <- rbind(RMSE_FLEX1.2val_strates, get(paste0("RMSE_FLEX1.2val_", j)) )
  rownames(RMSE_FLEX1.2val_strates)[dim(RMSE_FLEX1.2val_strates)[1]] <- paste0(j)
  rm(list = paste0("RMSE_FLEX1.2val_", j))
}

colnames(biais_FLEX1.2val_strates) <- colnames(ALL_flex1.2)
colnames(RMSE_FLEX1.2val_strates) <- colnames(ALL_flex1.2)

##4 noeuds
biais_FLEX1.4val_strates <- data.frame()
RMSE_FLEX1.4val_strates <- data.frame()
for(j in strata_names){
  assign(paste0("biais_FLEX1.4val_", j), 
         colMeans( get( paste0("ALL_FLEX1.4val_", j) ) - get(paste0("S_theoval", j) )[col(get( paste0("ALL_FLEX1.4val_", j)) )] ) )
  
  biais_FLEX1.4val_strates <- rbind(biais_FLEX1.4val_strates, get(paste0("biais_FLEX1.4val_", j)) )
  rownames(biais_FLEX1.4val_strates)[dim(biais_FLEX1.4val_strates)[1]] <- paste0(j)
  rm(list = paste0("biais_FLEX1.4val_", j))
  
  assign(paste0("RMSE_FLEX1.4val_",  j) ,
         sqrt( colMeans( ( get( paste0("ALL_FLEX1.4val_", j) ) - 
                             get(paste0("S_theoval", j) )[col(get( paste0("ALL_FLEX1.4val_", j)) )] )^2 ) ) )
  
  RMSE_FLEX1.4val_strates <- rbind(RMSE_FLEX1.4val_strates, get(paste0("RMSE_FLEX1.4val_", j)) )
  rownames(RMSE_FLEX1.4val_strates)[dim(RMSE_FLEX1.4val_strates)[1]] <- paste0(j)
  rm(list = paste0("RMSE_FLEX1.4val_", j))
}

colnames(biais_FLEX1.4val_strates) <- colnames(ALL_flex1.4)
colnames(RMSE_FLEX1.4val_strates) <- colnames(ALL_flex1.4)

# FLEX 2

##2 noeuds
biais_FLEX2.2val_strates <- data.frame()
RMSE_FLEX2.2val_strates <- data.frame()

for(j in strata_names){
  assign(paste0("biais_FLEX2.2val_", j), 
         colMeans( get( paste0("ALL_FLEX2.2val_", j) ) - get(paste0("S_theoval", j) )[col(get( paste0("ALL_FLEX2.2val_", j)) )] ) )
  
  biais_FLEX2.2val_strates <- rbind(biais_FLEX2.2val_strates, get(paste0("biais_FLEX2.2val_", j)) )
  rownames(biais_FLEX2.2val_strates)[dim(biais_FLEX2.2val_strates)[1]] <- paste0(j)
  rm(list = paste0("biais_FLEX2.2val_", j))
  
  assign(paste0("RMSE_FLEX2.2val_",  j) ,
         sqrt( colMeans( ( get( paste0("ALL_FLEX2.2val_", j) ) - 
                             get(paste0("S_theoval", j) )[col(get( paste0("ALL_FLEX2.2val_", j)) )] )^2 ) ) )
  
  RMSE_FLEX2.2val_strates <- rbind(RMSE_FLEX2.2val_strates, get(paste0("RMSE_FLEX2.2val_", j)) )
  rownames(RMSE_FLEX2.2val_strates)[dim(RMSE_FLEX2.2val_strates)[1]] <- paste0(j)
  rm(list = paste0("RMSE_FLEX2.2val_", j))
  
}

colnames(biais_FLEX2.2val_strates) <- colnames(ALL_flex1.2)
colnames(RMSE_FLEX2.2val_strates) <- colnames(ALL_flex1.2)

##4 noeuds
biais_FLEX2.4val_strates <- data.frame()
RMSE_FLEX2.4val_strates <- data.frame()

for(j in strata_names){
  assign(paste0("biais_FLEX2.4val_", j), 
         colMeans( get( paste0("ALL_FLEX2.4val_", j) ) - get(paste0("S_theoval", j) )[col(get( paste0("ALL_FLEX2.4val_", j)) )] ) )
  
  biais_FLEX2.4val_strates <- rbind(biais_FLEX2.4val_strates, get(paste0("biais_FLEX2.4val_", j)) )
  rownames(biais_FLEX2.4val_strates)[dim(biais_FLEX2.4val_strates)[1]] <- paste0(j)
  rm(list = paste0("biais_FLEX2.4val_", j))
  
  assign(paste0("RMSE_FLEX2.4val_",  j) ,
         sqrt( colMeans( ( get( paste0("ALL_FLEX2.4val_", j) ) - 
                             get(paste0("S_theoval", j) )[col(get( paste0("ALL_FLEX2.4val_", j)) )] )^2 ) ) )
  
  RMSE_FLEX2.4val_strates <- rbind(RMSE_FLEX2.4val_strates, get(paste0("RMSE_FLEX2.4val_", j)) )
  rownames(RMSE_FLEX2.4val_strates)[dim(RMSE_FLEX2.4val_strates)[1]] <- paste0(j)
  rm(list = paste0("RMSE_FLEX2.4val_", j))
  
}

colnames(biais_FLEX2.4val_strates) <- colnames(ALL_flex1.2)
colnames(RMSE_FLEX2.4val_strates) <- colnames(ALL_flex1.2)
# WG 

biais_WGval_strates <- data.frame()
RMSE_WGval_strates <- data.frame()

for(j in strata_names){
  assign(paste0("biais_WGval_", j), 
         colMeans( get( paste0("ALL_WGval_", j) ) - get(paste0("S_theoval", j) )[col(get( paste0("ALL_WGval_", j)) )] ) )
  
  biais_WGval_strates <- rbind(biais_WGval_strates, get(paste0("biais_WGval_", j)) )
  rownames(biais_WGval_strates)[dim(biais_WGval_strates)[1]] <- paste0(j)
  rm(list = paste0("biais_WGval_", j))
  
  assign(paste0("RMSE_WGval_",  j) ,
         sqrt( colMeans( ( get( paste0("ALL_WGval_", j) ) - 
                             get(paste0("S_theoval", j) )[col(get( paste0("ALL_WGval_", j)) )] )^2 ) ) )
  
  RMSE_WGval_strates <- rbind(RMSE_WGval_strates, get(paste0("RMSE_WGval_", j)) )
  rownames(RMSE_WGval_strates)[dim(RMSE_WGval_strates)[1]] <- paste0(j)
  rm(list = paste0("RMSE_WGval_", j))
  
}

colnames(biais_WGval_strates) <- colnames(ALL_flex1.2)
colnames(RMSE_WGval_strates) <- colnames(ALL_flex1.2)

# PP ############################################## A CHANGER 
biais_PPval_strates <- data.frame()
RMSE_PPval_strates <- data.frame()
for(j in strata_names){
  assign(paste0("biais_PPval_", j), 
         colMeans( get( paste0("ALL_PP_", j) ) - get(paste0("S_theoval", j) )[col(get( paste0("ALL_PP_", j)) )] ) )
  
  biais_PPval_strates <- rbind(biais_PPval_strates, get(paste0("biais_PPval_", j)) )
  rownames(biais_PPval_strates)[dim(biais_PPval_strates)[1]] <- paste0(j)
  rm(list = paste0("biais_PPval_", j))
  
  assign(paste0("RMSE_PPval_",  j) ,
         sqrt( colMeans( ( get( paste0("ALL_PP_", j) ) - 
                             get(paste0("S_theoval", j) )[col(get( paste0("ALL_PP_", j)) )] )^2 ) ) )
  
  RMSE_PPval_strates <- rbind(RMSE_PPval_strates, get(paste0("RMSE_PPval_", j)) )
  rownames(RMSE_PPval_strates)[dim(RMSE_PPval_strates)[1]] <- paste0(j)
  rm(list = paste0("RMSE_PPval_", j))
}

colnames(biais_PPval_strates) <- colnames(ALL_PPval)
colnames(RMSE_PPval_strates) <- colnames(ALL_PPval)

###concaten des S_theo, là aussi pour libérer de la place dans l'evt global
S_theoval_strates <- data.frame()
for(j in strata_names){
  S_theoval_strates <- rbind(S_theoval_strates, get(paste0("S_theoval",j)))
  rownames(S_theoval_strates)[dim(S_theoval_strates)[1]] <- paste0(j)
  rm(list = paste0("S_theoval",j))
}
colnames(S_theoval_strates) <- colnames(ALL_WG)


## mise en forme de l'evt global 
##############################################
##concatenation STHEO/ STHEOSTRATES
##TRAIN
S_theo <- rbind(S_theo, S_theo_strates)
rownames(S_theo)[1] <- "ALL" 
rm(S_theo_strates)
##VALID

S_theoval <- rbind(S_theoval, S_theoval_strates)
rownames(S_theoval)[1] <- "ALL" 
rm(S_theoval_strates)
## concatenation des résultats sur la base totale avec ceux fait en strates 
##TRAIN
##PLANN
biais_PLANN <- rbind(biais_PLANN, biais_PLANN_strates)
RMSE_PLANN <- rbind(RMSE_PLANN, RMSE_PLANN_strates)
rownames(biais_PLANN)[1] <- "ALL"
rownames(RMSE_PLANN)[1] <- "ALL"
rm(biais_PLANN_strates)
rm(RMSE_PLANN_strates)
##FLEX1
##2noeuds
biais_FLEX1.2 <- rbind(biais_flex1.2, biais_FLEX1.2_strates)
RMSE_FLEX1.2 <- rbind(RMSE_flex1.2, RMSE_FLEX1.2_strates)
rownames(biais_FLEX1.2)[1] <- "ALL"
rownames(RMSE_FLEX1.2)[1] <- "ALL"
rm(biais_flex1.2, biais_FLEX1.2_strates)
rm(RMSE_flex1.2, RMSE_FLEX1.2_strates)

##4noeuds
biais_FLEX1.4 <- rbind(biais_flex1.4, biais_FLEX1.4_strates)
RMSE_FLEX1.4 <- rbind(RMSE_flex1.4, RMSE_FLEX1.4_strates)
rownames(biais_FLEX1.4)[1] <- "ALL"
rownames(RMSE_FLEX1.4)[1] <- "ALL"
rm(biais_flex1.4, biais_FLEX1.4_strates)
rm(RMSE_flex1.4, RMSE_FLEX1.4_strates)

##FLEX2
##2 noeuds
biais_FLEX2.2 <- rbind(biais_flex2.2, biais_FLEX2.2_strates)
RMSE_FLEX2.2 <- rbind(RMSE_flex2.2, RMSE_FLEX2.2_strates)
rownames(biais_FLEX2.2)[1] <- "ALL"
rownames(RMSE_FLEX2.2)[1] <- "ALL"
rm(biais_flex2.2, biais_FLEX2.2_strates)
rm(RMSE_flex2.2, RMSE_FLEX2.2_strates)
##4 noeuds
biais_FLEX2.4 <- rbind(biais_flex2.4, biais_FLEX2.4_strates)
RMSE_FLEX2.4 <- rbind(RMSE_flex2.4, RMSE_FLEX2.4_strates)
rownames(biais_FLEX2.4)[1] <- "ALL"
rownames(RMSE_FLEX2.4)[1] <- "ALL"
rm(biais_flex2.4, biais_FLEX2.4_strates)
rm(RMSE_flex2.4, RMSE_FLEX2.4_strates)
##WG

biais_WG <- rbind(biais_WG, biais_WG_strates)
RMSE_WG <- rbind(RMSE_WG, RMSE_WG_strates)
rownames(biais_WG)[1] <- "ALL"
rownames(RMSE_WG)[1] <- "ALL"
rm(biais_WG_strates)
rm(RMSE_WG_strates)

##PP

biais_PP <- rbind(biais_PP, biais_PP_strates)
RMSE_PP <- rbind(RMSE_PP, RMSE_PP_strates)
rownames(biais_PP)[1] <- "ALL"
rownames(RMSE_PP)[1] <- "ALL"
rm(biais_PP_strates)
rm(RMSE_PP_strates)

######VALID

biais_PLANNval <- rbind(biais_PLANNval, biais_PLANNval_strates)
RMSE_PLANNval <- rbind(RMSE_PLANNval, RMSE_PLANNval_strates)
rownames(biais_PLANNval)[1] <- "ALL"
rownames(RMSE_PLANNval)[1] <- "ALL"
rm(biais_PLANNval_strates)
rm(RMSE_PLANNval_strates)
##FLEX1
##2 noeuds
biais_FLEX1.2val <- rbind(biais_flex1.2val, biais_FLEX1.2val_strates)
RMSE_FLEX1.2val <- rbind(RMSE_flex1.2val, RMSE_FLEX1.2val_strates)
rownames(biais_FLEX1.2val)[1] <- "ALL"
rownames(RMSE_FLEX1.2val)[1] <- "ALL"
rm(biais_flex1.2val, biais_FLEX1.2val_strates)
rm(RMSE_flex1.2val, RMSE_FLEX1.2val_strates)
##4 noeuds
biais_FLEX1.4val <- rbind(biais_flex1.4val, biais_FLEX1.4val_strates)
RMSE_FLEX1.4val <- rbind(RMSE_flex1.4val, RMSE_FLEX1.4val_strates)
rownames(biais_FLEX1.4val)[1] <- "ALL"
rownames(RMSE_FLEX1.4val)[1] <- "ALL"
rm(biais_flex1.4val, biais_FLEX1.4val_strates)
rm(RMSE_flex1.4val, RMSE_FLEX1.4val_strates)
##FLEX2
##2 noeuds
biais_FLEX2.2val <- rbind(biais_flex2.2val, biais_FLEX2.2val_strates)
RMSE_FLEX2.2val <- rbind(RMSE_flex2.2val, RMSE_FLEX2.2val_strates)
rownames(biais_FLEX2.2val)[1] <- "ALL"
rownames(RMSE_FLEX2.2val)[1] <- "ALL"
rm(biais_flex2.2val, biais_FLEX2.2val_strates)
rm(RMSE_flex2.2val, RMSE_FLEX2.2val_strates)
##4 noeuds
biais_FLEX2.4val <- rbind(biais_flex2.4val, biais_FLEX2.4val_strates)
RMSE_FLEX2.4val <- rbind(RMSE_flex2.4val, RMSE_FLEX2.4val_strates)
rownames(biais_FLEX2.4val)[1] <- "ALL"
rownames(RMSE_FLEX2.4val)[1] <- "ALL"
rm(biais_flex2.4val, biais_FLEX2.4val_strates)
rm(RMSE_flex2.4val, RMSE_FLEX2.4val_strates)
##WG

biais_WGval <- rbind(biais_WGval, biais_WGval_strates)
RMSE_WGval <- rbind(RMSE_WGval, RMSE_WGval_strates)
rownames(biais_WGval)[1] <- "ALL"
rownames(RMSE_WGval)[1] <- "ALL"
rm(biais_WGval_strates)
rm(RMSE_WGval_strates)

##PP

biais_PPval <- rbind(biais_PPval, biais_PPval_strates)
RMSE_PPval <- rbind(RMSE_PPval, RMSE_PPval_strates)
rownames(biais_PPval)[1] <- "ALL"
rownames(RMSE_PPval)[1] <- "ALL"
rm(biais_PPval_strates)
rm(RMSE_PPval_strates)

### les dataframes des moyennes sont enregistrées sous forme de liste par rapport au modèle utilisé
########
##TRAIN
THEO_means_names <- c()
for(j in strata_names){
  THEO_means_names <- c(THEO_means_names, paste0("ALL_",j,"_theo"))
}
THEO_means <- c(mget("ALL_theo"), mget(THEO_means_names))
rm(ALL_theo)
rm(THEO_means_names)
for(j in strata_names){
  rm(list = paste0("ALL_",j,"_theo") )
}

##PLANN
PLANN_means_names <- c()
for(j in strata_names){
  PLANN_means_names <- c(PLANN_means_names, paste0("ALL_PLANN_",j))
}
PLANN_means <- c(mget("ALL_PLANN"), mget(PLANN_means_names))
rm(ALL_PLANN)
rm(PLANN_means_names)
for(j in strata_names){
  rm(list = paste0("ALL_PLANN_",j) )
}

##FLEX1
##2 noeuds
FLEX1.2_means_names <- c()
for(j in strata_names){
  FLEX1.2_means_names <- c(FLEX1.2_means_names, paste0("ALL_FLEX1.2_",j))
}
FLEX1.2_means <- c(mget("ALL_flex1.2"), mget(FLEX1.2_means_names))
rm(ALL_flex1.2)
rm(FLEX1.2_means_names)
for(j in strata_names){
  rm(list = paste0("ALL_FLEX1.2_",j) )
}

##4 noeuds
FLEX1.4_means_names <- c()
for(j in strata_names){
  FLEX1.4_means_names <- c(FLEX1.4_means_names, paste0("ALL_FLEX1.4_",j))
}
FLEX1.4_means <- c(mget("ALL_flex1.4"), mget(FLEX1.4_means_names))
rm(ALL_flex1.4)
rm(FLEX1.4_means_names)
for(j in strata_names){
  rm(list = paste0("ALL_FLEX1.4_",j) )
}

##FLEX2
##2 noeuds
FLEX2.2_means_names <- c()
for(j in strata_names){
  FLEX2.2_means_names <- c(FLEX2.2_means_names, paste0("ALL_FLEX2.2_",j))
}
FLEX2.2_means <- c(mget("ALL_flex2.2"), mget(FLEX2.2_means_names))
rm(ALL_flex2.2)
rm(FLEX2.2_means_names)
for(j in strata_names){
  rm(list = paste0("ALL_FLEX2.2_",j) )
}
##4 noeuds
FLEX2.4_means_names <- c()
for(j in strata_names){
  FLEX2.4_means_names <- c(FLEX2.4_means_names, paste0("ALL_FLEX2.4_",j))
}
FLEX2.4_means <- c(mget("ALL_flex2.4"), mget(FLEX2.4_means_names))
rm(ALL_flex2.4)
rm(FLEX2.4_means_names)
for(j in strata_names){
  rm(list = paste0("ALL_FLEX2.4_",j) )
}
##WG
WG_means_names <- c()
for(j in strata_names){
  WG_means_names <- c(WG_means_names, paste0("ALL_WG_",j))
}
WG_means <- c(mget("ALL_WG"), mget(WG_means_names))
rm(ALL_WG)
rm(WG_means_names)
for(j in strata_names){
  rm(list = paste0("ALL_WG_",j) )
}

PP_means_names <- c()
for(j in strata_names){
  PP_means_names <- c(PP_means_names, paste0("ALL_PP_",j))
}
PP_means <- c(mget("ALL_PP"), mget(PP_means_names))
rm(ALL_PP)
rm(PP_means_names)
for(j in strata_names){
  rm(list = paste0("ALL_PP_",j) )
}
############VALID
THEOval_means_names <- c()
for(j in strata_names){
  THEOval_means_names <- c(THEOval_means_names, paste0("ALL_",j,"_theoval"))
}
THEOval_means <- c(mget("ALL_theoval"), mget(THEOval_means_names))
rm(ALL_theoval)
rm(THEOval_means_names)
for(j in strata_names){
  rm(list = paste0("ALL_",j,"_theoval") )
}

##PLANN
PLANNval_means_names <- c()
for(j in strata_names){
  PLANNval_means_names <- c(PLANNval_means_names, paste0("ALL_PLANNval_",j))
}
PLANNval_means <- c(mget("ALL_PLANNval"), mget(PLANNval_means_names))
rm(ALL_PLANNval)
rm(PLANNval_means_names)
for(j in strata_names){
  rm(list = paste0("ALL_PLANNval_",j) )
}

##FLEX1
##2 noeuds
FLEX1.2val_means_names <- c()
for(j in strata_names){
  FLEX1.2val_means_names <- c(FLEX1.2val_means_names, paste0("ALL_FLEX1.2val_",j))
}
FLEX1.2val_means <- c(mget("ALL_flex1.2val"), mget(FLEX1.2val_means_names))
rm(ALL_flex1.2val)
rm(FLEX1.2val_means_names)
for(j in strata_names){
  rm(list = paste0("ALL_FLEX1.2val_",j) )
}

##4 noeuds
FLEX1.4val_means_names <- c()
for(j in strata_names){
  FLEX1.4val_means_names <- c(FLEX1.4val_means_names, paste0("ALL_FLEX1.4val_",j))
}
FLEX1.4val_means <- c(mget("ALL_flex1.4val"), mget(FLEX1.4val_means_names))
rm(ALL_flex1.4val)
rm(FLEX1.4val_means_names)
for(j in strata_names){
  rm(list = paste0("ALL_FLEX1.4val_",j) )
}
##FLEX2
##2 noeuds

FLEX2.2val_means_names <- c()
for(j in strata_names){
  FLEX2.2val_means_names <- c(FLEX2.2val_means_names, paste0("ALL_FLEX2.2val_",j))
}
FLEX2.2val_means <- c(mget("ALL_flex2.2val"), mget(FLEX2.2val_means_names))
rm(ALL_flex2.2val)
rm(FLEX2.2val_means_names)
for(j in strata_names){
  rm(list = paste0("ALL_FLEX2.2val_",j) )
}

##4 noeuds

FLEX2.4val_means_names <- c()
for(j in strata_names){
  FLEX2.4val_means_names <- c(FLEX2.4val_means_names, paste0("ALL_FLEX2.4val_",j))
}
FLEX2.4val_means <- c(mget("ALL_flex2.4val"), mget(FLEX2.4val_means_names))
rm(ALL_flex2.4val)
rm(FLEX2.4val_means_names)
for(j in strata_names){
  rm(list = paste0("ALL_FLEX2.4val_",j) )
}
##WG
WGval_means_names <- c()
for(j in strata_names){
  WGval_means_names <- c(WGval_means_names, paste0("ALL_WGval_",j))
}
WGval_means <- c(mget("ALL_WGval"), mget(WGval_means_names))
rm(ALL_WGval)
rm(WGval_means_names)
for(j in strata_names){
  rm(list = paste0("ALL_WGval_",j) )
}

##PP
PPval_means_names <- c()
for(j in strata_names){
  PPval_means_names <- c(PPval_means_names, paste0("ALL_PPval_",j))
}
PPval_means <- c(mget("ALL_PPval"), mget(PPval_means_names))
rm(ALL_PPval)
rm(PPval_means_names)
for(j in strata_names){
  rm(list = paste0("ALL_PPval_",j) )
}


packages <- c("PLNmodels", "ggplot2","corrplot","future.apply")
for (package in packages){
  if(! package %in% installed.packages()){
    install.packages(package, dependencies = TRUE)
  }
}

library(PLNmodels)
library(ggplot2)
library(corrplot)
library(future.apply)
setwd("/home/bernardt/Documents/StageTanguy/results/")

options(future.globals.maxSize = +Inf)
source("functions_jsdm.R")

COV <- read.csv("covariable_r.csv", header = F)
list_n <- c("abundance_r_0.csv", 
            "abundance_r_1.csv",
            "abundance_r_2.csv",
            "abundance_r_3.csv",
            "abundance_r_4.csv",
            "abundance_r_5.csv",
            "abundance_r_6.csv",
            "abundance_r_7.csv",
            "abundance_r_8.csv",
            "abundance_r_9.csv",
            "abundance_r_10.csv", 
            "abundance_r_11.csv",
            "abundance_r_12.csv",
            "abundance_r_13.csv",
            "abundance_r_14.csv",
            "abundance_r_15.csv",
            "abundance_r_16.csv",
            "abundance_r_17.csv",
            "abundance_r_18.csv",
            "abundance_r_19.csv",
            "abundance_r_20.csv", 
            "abundance_r_21.csv",
            "abundance_r_22.csv",
            "abundance_r_23.csv",
            "abundance_r_24.csv",
            "abundance_r_25.csv",
            "abundance_r_26.csv"
            )

list_r <- c(1.4, 
  1.5,
  1.6,
  1.7,
  1.8,
  1.9,
  2,
  2.1,
  2.2,
  2.3,
  2.4,
  2.5,
  2.6,
  2.7,
  2.8,
  2.9,
  3,
  3.1,
  3.2,
  3.3,
  3.4,
  3.5,
  3.6,
  3.7,
  3.8,
  3.9,
  4
)

ITER <- list()
for(i in 0:29){
  for(j in 0:39){
    ITER[[length(ITER)+1]] <- c(i,j)
  }
}

plan(multisession, workers = 27)

RESULTS_JSDM <- future_lapply(list_n, MAIN_PROG)


MEAN_PDI_GE_1 <- list()
MEAN_PDI_GE_2 <- list()
MEAN_PDI_SP1_1 <- list()
MEAN_PDI_SP1_2 <- list()
MEAN_PDI_SP2_1 <- list()
MEAN_PDI_SP2_2 <- list()

SD_PDI_GE_1 <- list()
SD_PDI_GE_2 <- list()
SD_PDI_SP1_1 <- list()
SD_PDI_SP1_2 <- list()
SD_PDI_SP2_1 <- list()
SD_PDI_SP2_2 <- list()

MEAN_SIG_GE_GE_0 <- list()
MEAN_SIG_GE_GE_1 <- list()
MEAN_SIG_GE_GE_2 <- list()
MEAN_SIG_GE_SP1_0 <- list()
MEAN_SIG_GE_SP1_1 <- list()
MEAN_SIG_GE_SP1_2 <- list()
MEAN_SIG_GE_SP2_0 <- list()
MEAN_SIG_GE_SP2_1 <- list()
MEAN_SIG_GE_SP2_2 <- list()
MEAN_SIG_SP1_SP1_0 <- list()
MEAN_SIG_SP1_SP1_1 <- list()
MEAN_SIG_SP1_SP1_2 <- list()
MEAN_SIG_SP1_SP2_0 <- list()
MEAN_SIG_SP1_SP2_1 <- list()
MEAN_SIG_SP1_SP2_2 <- list()

SD_SIG_GE_GE_0 <- list()
SD_SIG_GE_GE_1 <- list()
SD_SIG_GE_GE_2 <- list()
SD_SIG_GE_SP1_0 <- list()
SD_SIG_GE_SP1_1 <- list()
SD_SIG_GE_SP1_2 <- list()
SD_SIG_GE_SP2_0 <- list()
SD_SIG_GE_SP2_1 <- list()
SD_SIG_GE_SP2_2 <- list()
SD_SIG_SP1_SP1_0 <- list()
SD_SIG_SP1_SP1_1 <- list()
SD_SIG_SP1_SP1_2 <- list()
SD_SIG_SP1_SP2_0 <- list()
SD_SIG_SP1_SP2_1 <- list()
SD_SIG_SP1_SP2_2 <- list()

for(i in 1:length(RESULTS_JSDM)){
  for(j in 1:4){
    RESULTS_JSDM[[i]][[j]] <- cbind(RESULTS_JSDM[[i]][[j]], r = list_r[i])
  }
  
  MEAN_PDI_GE_1[[i]] <- cbind(subset(RESULTS_JSDM[[i]][[1]][c(1,4:6)], Model == "Res"), Sp = "GE")
  MEAN_PDI_GE_2[[i]] <- cbind(subset(RESULTS_JSDM[[i]][[1]][c(1,4:6)], Model == "Res+t"), Sp = "GE")
  MEAN_PDI_SP1_1[[i]] <- cbind(subset(RESULTS_JSDM[[i]][[1]][c(2,4:6)], Model == "Res"), Sp = "SP1")
  MEAN_PDI_SP1_2[[i]] <- cbind(subset(RESULTS_JSDM[[i]][[1]][c(2,4:6)], Model == "Res+t"), Sp = "SP1")
  MEAN_PDI_SP2_1[[i]] <- cbind(subset(RESULTS_JSDM[[i]][[1]][c(3,4:6)], Model == "Res"), Sp = "SP2")
  MEAN_PDI_SP2_2[[i]] <- cbind(subset(RESULTS_JSDM[[i]][[1]][c(3,4:6)], Model == "Res+t"), Sp = "SP2")
  
  SD_PDI_GE_1[[i]] <-  cbind(subset(RESULTS_JSDM[[i]][[2]][c(1,4:6)], Model == "Res"), Sp = "GE")
  SD_PDI_GE_2[[i]] <-  cbind(subset(RESULTS_JSDM[[i]][[2]][c(1,4:6)], Model == "Res+t"), Sp = "GE")
  SD_PDI_SP1_1[[i]] <-  cbind(subset(RESULTS_JSDM[[i]][[2]][c(2,4:6)], Model == "Res"), Sp = "SP1")
  SD_PDI_SP1_2[[i]] <-  cbind(subset(RESULTS_JSDM[[i]][[2]][c(2,4:6)], Model == "Res+t"), Sp = "SP1")
  SD_PDI_SP2_1[[i]] <-  cbind(subset(RESULTS_JSDM[[i]][[2]][c(3,4:6)], Model == "Res"), Sp = "SP2")
  SD_PDI_SP2_2[[i]] <-  cbind(subset(RESULTS_JSDM[[i]][[2]][c(3,4:6)], Model == "Res+t"), Sp = "SP2")
  
  MEAN_SIG_GE_GE_0[[i]] <- cbind(subset(RESULTS_JSDM[[i]][[3]][c(1,6:8)], Model == "Null"), Int = "GE:GE")
  MEAN_SIG_GE_GE_1[[i]] <- cbind(subset(RESULTS_JSDM[[i]][[3]][c(1,6:8)], Model == "Res"), Int = "GE:GE")
  MEAN_SIG_GE_GE_2[[i]] <- cbind(subset(RESULTS_JSDM[[i]][[3]][c(1,6:8)], Model == "Res+t"), Int = "GE:GE")
  MEAN_SIG_GE_SP1_0[[i]] <- cbind(subset(RESULTS_JSDM[[i]][[3]][c(2,6:8)], Model == "Null"), Int = "GE:SP1")
  MEAN_SIG_GE_SP1_1[[i]] <- cbind(subset(RESULTS_JSDM[[i]][[3]][c(2,6:8)], Model == "Res"), Int = "GE:SP1")
  MEAN_SIG_GE_SP1_2[[i]] <- cbind(subset(RESULTS_JSDM[[i]][[3]][c(2,6:8)], Model == "Res+t"), Int = "GE:SP1")
  MEAN_SIG_GE_SP2_0[[i]] <- cbind(subset(RESULTS_JSDM[[i]][[3]][c(3,6:8)], Model == "Null"), Int = "GE:SP2")
  MEAN_SIG_GE_SP2_1[[i]] <- cbind(subset(RESULTS_JSDM[[i]][[3]][c(3,6:8)], Model == "Res"), Int = "GE:SP2")
  MEAN_SIG_GE_SP2_2[[i]] <- cbind(subset(RESULTS_JSDM[[i]][[3]][c(3,6:8)], Model == "Res+t"), Int = "GE:SP2")
  MEAN_SIG_SP1_SP1_0[[i]] <- cbind(subset(RESULTS_JSDM[[i]][[3]][c(4,6:8)], Model == "Null"), Int = "SP1:SP1")
  MEAN_SIG_SP1_SP1_1[[i]] <- cbind(subset(RESULTS_JSDM[[i]][[3]][c(4,6:8)], Model == "Res"), Int = "SP1:SP1")
  MEAN_SIG_SP1_SP1_2[[i]] <- cbind(subset(RESULTS_JSDM[[i]][[3]][c(4,6:8)], Model == "Res+t"), Int = "SP1:SP1")
  MEAN_SIG_SP1_SP2_0[[i]] <- cbind(subset(RESULTS_JSDM[[i]][[3]][c(5,6:8)], Model == "Null"), Int = "SP1:SP2")
  MEAN_SIG_SP1_SP2_1[[i]] <- cbind(subset(RESULTS_JSDM[[i]][[3]][c(5,6:8)], Model == "Res"), Int = "SP1:SP2")
  MEAN_SIG_SP1_SP2_2[[i]] <- cbind(subset(RESULTS_JSDM[[i]][[3]][c(5,6:8)], Model == "Res+t"), Int = "SP1:SP2")
  
  SD_SIG_GE_GE_0[[i]] <- cbind(subset(RESULTS_JSDM[[i]][[4]][c(1,6:8)], Model == "Null"), Int = "GE:GE")
  SD_SIG_GE_GE_1[[i]] <- cbind(subset(RESULTS_JSDM[[i]][[4]][c(1,6:8)], Model == "Res"), Int = "GE:GE")
  SD_SIG_GE_GE_2[[i]] <- cbind(subset(RESULTS_JSDM[[i]][[4]][c(1,6:8)], Model == "Res+t"), Int = "GE:GE")
  SD_SIG_GE_SP1_0[[i]] <- cbind(subset(RESULTS_JSDM[[i]][[4]][c(2,6:8)], Model == "Null"), Int = "GE:SP1")
  SD_SIG_GE_SP1_1[[i]] <- cbind(subset(RESULTS_JSDM[[i]][[4]][c(2,6:8)], Model == "Res"), Int = "GE:SP1")
  SD_SIG_GE_SP1_2[[i]] <- cbind(subset(RESULTS_JSDM[[i]][[4]][c(2,6:8)], Model == "Res+t"), Int = "GE:SP1")
  SD_SIG_GE_SP2_0[[i]] <- cbind(subset(RESULTS_JSDM[[i]][[4]][c(3,6:8)], Model == "Null"), Int = "GE:SP2")
  SD_SIG_GE_SP2_1[[i]] <- cbind(subset(RESULTS_JSDM[[i]][[4]][c(3,6:8)], Model == "Res"), Int = "GE:SP2")
  SD_SIG_GE_SP2_2[[i]] <- cbind(subset(RESULTS_JSDM[[i]][[4]][c(3,6:8)], Model == "Res+t"), Int = "GE:SP2")
  SD_SIG_SP1_SP1_0[[i]] <- cbind(subset(RESULTS_JSDM[[i]][[4]][c(4,6:8)], Model == "Null"), Int = "SP1:SP1")
  SD_SIG_SP1_SP1_1[[i]] <- cbind(subset(RESULTS_JSDM[[i]][[4]][c(4,6:8)], Model == "Res"), Int = "SP1:SP1")
  SD_SIG_SP1_SP1_2[[i]] <- cbind(subset(RESULTS_JSDM[[i]][[4]][c(4,6:8)], Model == "Res+t"), Int = "SP1:SP1")
  SD_SIG_SP1_SP2_0[[i]] <- cbind(subset(RESULTS_JSDM[[i]][[4]][c(5,6:8)], Model == "Null"), Int = "SP1:SP2")
  SD_SIG_SP1_SP2_1[[i]] <- cbind(subset(RESULTS_JSDM[[i]][[4]][c(5,6:8)], Model == "Res"), Int = "SP1:SP2")
  SD_SIG_SP1_SP2_2[[i]] <- cbind(subset(RESULTS_JSDM[[i]][[4]][c(5,6:8)], Model == "Res+t"), Int = "SP1:SP2")
}

MEAN_PDI_GE_1 <- do.call(rbind.data.frame, MEAN_PDI_GE_1)
MEAN_PDI_GE_2 <- do.call(rbind.data.frame, MEAN_PDI_GE_2)
MEAN_PDI_SP1_1 <- do.call(rbind.data.frame, MEAN_PDI_SP1_1)
MEAN_PDI_SP1_2 <- do.call(rbind.data.frame, MEAN_PDI_SP1_2)
MEAN_PDI_SP2_1 <- do.call(rbind.data.frame, MEAN_PDI_SP2_1)
MEAN_PDI_SP2_2 <- do.call(rbind.data.frame, MEAN_PDI_SP2_2)

SD_PDI_GE_1 <- do.call(rbind.data.frame, SD_PDI_GE_1)
SD_PDI_GE_2 <- do.call(rbind.data.frame, SD_PDI_GE_2)
SD_PDI_SP1_1 <- do.call(rbind.data.frame, SD_PDI_SP1_1)
SD_PDI_SP1_2 <- do.call(rbind.data.frame, SD_PDI_SP1_2)
SD_PDI_SP2_1 <- do.call(rbind.data.frame, SD_PDI_SP2_1)
SD_PDI_SP2_2 <- do.call(rbind.data.frame, SD_PDI_SP2_2)

MEAN_SIG_GE_GE_0 <- do.call(rbind.data.frame, MEAN_SIG_GE_GE_0)
MEAN_SIG_GE_GE_1 <- do.call(rbind.data.frame, MEAN_SIG_GE_GE_1)
MEAN_SIG_GE_GE_2 <- do.call(rbind.data.frame, MEAN_SIG_GE_GE_2)
MEAN_SIG_GE_SP1_0 <- do.call(rbind.data.frame, MEAN_SIG_GE_SP1_0)
MEAN_SIG_GE_SP1_1 <- do.call(rbind.data.frame, MEAN_SIG_GE_SP1_1)
MEAN_SIG_GE_SP1_2 <- do.call(rbind.data.frame, MEAN_SIG_GE_SP1_2)
MEAN_SIG_GE_SP2_0 <- do.call(rbind.data.frame, MEAN_SIG_GE_SP2_0)
MEAN_SIG_GE_SP2_1 <- do.call(rbind.data.frame, MEAN_SIG_GE_SP2_1)
MEAN_SIG_GE_SP2_2 <- do.call(rbind.data.frame, MEAN_SIG_GE_SP2_2)
MEAN_SIG_SP1_SP1_0 <- do.call(rbind.data.frame, MEAN_SIG_SP1_SP1_0)
MEAN_SIG_SP1_SP1_1 <- do.call(rbind.data.frame, MEAN_SIG_SP1_SP1_1)
MEAN_SIG_SP1_SP1_2 <- do.call(rbind.data.frame, MEAN_SIG_SP1_SP1_2)
MEAN_SIG_SP1_SP2_0 <- do.call(rbind.data.frame, MEAN_SIG_SP1_SP2_0)
MEAN_SIG_SP1_SP2_1 <- do.call(rbind.data.frame, MEAN_SIG_SP1_SP2_1)
MEAN_SIG_SP1_SP2_2 <- do.call(rbind.data.frame, MEAN_SIG_SP1_SP2_2)

SD_SIG_GE_GE_0 <- do.call(rbind.data.frame, SD_SIG_GE_GE_0)
SD_SIG_GE_GE_1 <- do.call(rbind.data.frame, SD_SIG_GE_GE_1)
SD_SIG_GE_GE_2 <- do.call(rbind.data.frame, SD_SIG_GE_GE_2)
SD_SIG_GE_SP1_0 <- do.call(rbind.data.frame, SD_SIG_GE_SP1_0)
SD_SIG_GE_SP1_1 <- do.call(rbind.data.frame, SD_SIG_GE_SP1_1)
SD_SIG_GE_SP1_2 <- do.call(rbind.data.frame, SD_SIG_GE_SP1_2)
SD_SIG_GE_SP2_0 <- do.call(rbind.data.frame, SD_SIG_GE_SP2_0)
SD_SIG_GE_SP2_1 <- do.call(rbind.data.frame, SD_SIG_GE_SP2_1)
SD_SIG_GE_SP2_2 <- do.call(rbind.data.frame, SD_SIG_GE_SP2_2)
SD_SIG_SP1_SP1_0 <- do.call(rbind.data.frame, SD_SIG_SP1_SP1_0)
SD_SIG_SP1_SP1_1 <- do.call(rbind.data.frame, SD_SIG_SP1_SP1_1)
SD_SIG_SP1_SP1_2 <- do.call(rbind.data.frame, SD_SIG_SP1_SP1_2)
SD_SIG_SP1_SP2_0 <- do.call(rbind.data.frame, SD_SIG_SP1_SP2_0)
SD_SIG_SP1_SP2_1 <- do.call(rbind.data.frame, SD_SIG_SP1_SP2_1)
SD_SIG_SP1_SP2_2 <- do.call(rbind.data.frame, SD_SIG_SP1_SP2_2)

MEAN_PDI <- do.call(rbind.data.frame,lapply(list(MEAN_PDI_GE_1, MEAN_PDI_GE_2, MEAN_PDI_SP1_1, MEAN_PDI_SP1_2, MEAN_PDI_SP2_1, MEAN_PDI_SP2_2),
       setNames, c("Mean_PDI", "Model", "Age", "r", "Guild")))

SD_PDI <- do.call(rbind.data.frame,lapply(list(SD_PDI_GE_1, SD_PDI_GE_2, SD_PDI_SP1_1, SD_PDI_SP1_2, SD_PDI_SP2_1, SD_PDI_SP2_2),
       setNames, c("SD_PDI", "Model", "Age", "r", "Guild")))

MEAN_SIGMA <- do.call(rbind.data.frame,lapply(list(MEAN_SIG_GE_GE_0, MEAN_SIG_GE_GE_1, MEAN_SIG_GE_GE_2, MEAN_SIG_GE_SP1_0, MEAN_SIG_GE_SP1_1, MEAN_SIG_GE_SP1_2,
            MEAN_SIG_GE_SP2_0, MEAN_SIG_GE_SP2_1, MEAN_SIG_GE_SP2_2, MEAN_SIG_SP1_SP1_0, MEAN_SIG_SP1_SP1_1, MEAN_SIG_SP1_SP1_2,
            MEAN_SIG_SP1_SP2_0, MEAN_SIG_SP1_SP2_1, MEAN_SIG_SP1_SP2_2),
       setNames, c("Mean_Sigma", "Model", "Age", "r", "Interaction")))

SD_SIGMA <- do.call(rbind.data.frame,lapply(list(SD_SIG_GE_GE_0, SD_SIG_GE_GE_1, SD_SIG_GE_GE_2, SD_SIG_GE_SP1_0, SD_SIG_GE_SP1_1, SD_SIG_GE_SP1_2,
            SD_SIG_GE_SP2_0, SD_SIG_GE_SP2_1, SD_SIG_GE_SP2_2, SD_SIG_SP1_SP1_0, SD_SIG_SP1_SP1_1, SD_SIG_SP1_SP1_2,
            SD_SIG_SP1_SP2_0, SD_SIG_SP1_SP2_1, SD_SIG_SP1_SP2_2),
       setNames, c("SD_Sigma", "Model", "Age", "r", "Interaction")))

MEAN_PDI$r <- as.factor(MEAN_PDI$r)
SD_PDI$r <- as.factor(SD_PDI$r)
MEAN_SIGMA$r <- as.factor(MEAN_SIGMA$r)
SD_SIGMA$r <- as.factor(SD_SIGMA$r)

write.csv(MEAN_PDI, "mean_PDI.csv", row.names=FALSE)
write.csv(SD_PDI, "sd_PDI.csv", row.names=FALSE)
write.csv(MEAN_SIGMA, "mean_sigma.csv", row.names=FALSE)
write.csv(SD_SIGMA, "sd_sigma.csv", row.names=FALSE)


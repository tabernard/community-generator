PDI <- function(v){
  v[] <- v/max(v)
  v[] <- max(v) - v
  s <- sum(v)/2
  return(s)
}

LAP_ITER <- function(i, n, cov){
  sub_cov = subset(cov, 
                   Iter==i[1]
                   & Plan==i[2]
  )
  sub_n = subset(n, 
                 Iter==i[1]
                 & Plan==i[2]
  )
  
  n_pln = sub_n[,1:8]
  
  pln <- prepare_data(n_pln, sub_cov)
  PLN0 <- PLN(Abundance ~ 1, data = pln)
  PLN1 <- PLN(Abundance ~ 0+Res, data = pln)
  PLN2 <- PLN(Abundance ~ 0+Res+Time, data = pln)
  
  a <- as.matrix(coef(PLN1))
  a[] <- exp(a)
  
  b <- as.matrix(coef(PLN2)[1:3,])
  b[] <- exp(b)
  
  pdi_a <- apply(a,2,PDI)
  pdi_b <- apply(b,2,PDI)
  
  pdi <- rbind(pdi_a,pdi_b)
  rownames(pdi) <- c("Res","Res+t")
  
  sig0 = cov2cor(sigma(PLN0))
  siga = cov2cor(sigma(PLN1))
  sigb = cov2cor(sigma(PLN2))
  
  return(list(pdi=pdi,sig0=sig0,siga=siga,sigb=sigb))
}

LAP_PDI <- function(i,n,cov){
  cov = subset(cov,
               Age==i
  )
  n = subset(n, 
             Age==i
  )
  iter <- ITER
  results <- lapply(iter,LAP_ITER,cov=cov,n=n)
  res_pdi <- lapply(results, "[[", 1)
  res_sig0 <- lapply(results, "[[", 2)
  res_siga <- lapply(results, "[[", 3)
  res_sigb <- lapply(results, "[[", 4)
  
  mean_pdi <- apply(simplify2array(res_pdi), 1:2, mean)
  sd_pdi <- apply(simplify2array(res_pdi), 1:2, sd)
  
  mean_pdi_gp <- data.frame(matrix(ncol = 5, nrow = 2))
  mean_pdi_gp[,1] <- mean(mean_pdi[,1:4])
  mean_pdi_gp[,2] <- mean(mean_pdi[,5:7])
  mean_pdi_gp[,3] <- mean_pdi[,8]
  mean_pdi_gp[,4] <- rownames(mean_pdi)
  mean_pdi_gp[,5] <- i
  colnames(mean_pdi_gp) <- c("GE", "SP1", "SP2","Model","Age")
  
  sd_pdi_gp <- data.frame(matrix(ncol = 5, nrow = 2))
  sd_pdi_gp[,1] <- mean(sd_pdi[,1:4])
  sd_pdi_gp[,2] <- mean(sd_pdi[,5:7])
  sd_pdi_gp[,3] <- sd_pdi[,8]
  sd_pdi_gp[,4] <- rownames(sd_pdi)
  sd_pdi_gp[,5] <- i
  colnames(sd_pdi_gp) <- c("GE", "SP1", "SP2","Model","Age")
  
  
  mean_sig0 <- apply(simplify2array(res_sig0), 1:2, mean)
  sd_sig0 <- apply(simplify2array(res_sig0), 1:2, sd)
  
  mean_siga <- apply(simplify2array(res_siga), 1:2, mean)
  sd_siga <- apply(simplify2array(res_siga), 1:2, sd)
  
  mean_sigb <- apply(simplify2array(res_sigb), 1:2, mean)
  sd_sigb <- apply(simplify2array(res_sigb), 1:2, sd)
  
  mean_sig_gp <- data.frame(matrix(ncol = 7, nrow = 3))
  mean_sig_gp[1,1] <- mean(c(mean_sig0[1,2:4],mean_sig0[2,3:4],mean_sig0[3,4]))
  mean_sig_gp[1,2] <- mean(c(mean_sig0[1,5:7],mean_sig0[2,5:7],mean_sig0[3,5:7],mean_sig0[4,5:7]))
  mean_sig_gp[1,3] <- mean(mean_sig0[1:4,8])
  mean_sig_gp[1,4] <- mean(c(mean_sig0[5,6:7],mean_sig0[6,7]))
  mean_sig_gp[1,5] <- mean(mean_sig0[5:7,8])
  mean_sig_gp[1,6] <- "Null"
  mean_sig_gp[1,7] <- i
  mean_sig_gp[2,1] <- mean(c(mean_siga[1,2:4],mean_siga[2,3:4],mean_siga[3,4]))
  mean_sig_gp[2,2] <- mean(c(mean_siga[1,5:7],mean_siga[2,5:7],mean_siga[3,5:7],mean_siga[4,5:7]))
  mean_sig_gp[2,3] <- mean(mean_siga[1:4,8])
  mean_sig_gp[2,4] <- mean(c(mean_siga[5,6:7],mean_siga[6,7]))
  mean_sig_gp[2,5] <- mean(mean_siga[5:7,8])
  mean_sig_gp[2,6] <- "Res"
  mean_sig_gp[2,7] <- i
  mean_sig_gp[3,1] <- mean(c(mean_sigb[1,2:4],mean_sigb[2,3:4],mean_sigb[3,4]))
  mean_sig_gp[3,2] <- mean(c(mean_sigb[1,5:7],mean_sigb[2,5:7],mean_sigb[3,5:7],mean_sigb[4,5:7]))
  mean_sig_gp[3,3] <- mean(mean_sigb[1:4,8])
  mean_sig_gp[3,4] <- mean(c(mean_sigb[5,6:7],mean_sigb[6,7]))
  mean_sig_gp[3,5] <- mean(mean_sigb[5:7,8])
  mean_sig_gp[3,6] <- "Res+t"
  mean_sig_gp[3,7] <- i
  colnames(mean_sig_gp) <- c("GE:GE", "GE:SP1", "GE:SP2","SP1:SP1","SP1:SP2","Model","Age")
  
  sd_sig_gp <- data.frame(matrix(ncol = 7, nrow = 3))
  sd_sig_gp[1,1] <- mean(c(sd_sig0[1,2:4],sd_sig0[2,3:4],sd_sig0[3,4]))
  sd_sig_gp[1,2] <- mean(c(sd_sig0[1,5:7],sd_sig0[2,5:7],sd_sig0[3,5:7],sd_sig0[4,5:7]))
  sd_sig_gp[1,3] <- mean(sd_sig0[1:4,8])
  sd_sig_gp[1,4] <- mean(c(sd_sig0[5,6:7],sd_sig0[6,7]))
  sd_sig_gp[1,5] <- mean(sd_sig0[5:7,8])
  sd_sig_gp[1,6] <- "Null"
  sd_sig_gp[1,7] <- i
  sd_sig_gp[2,1] <- mean(c(sd_siga[1,2:4],sd_siga[2,3:4],sd_siga[3,4]))
  sd_sig_gp[2,2] <- mean(c(sd_siga[1,5:7],sd_siga[2,5:7],sd_siga[3,5:7],sd_siga[4,5:7]))
  sd_sig_gp[2,3] <- mean(sd_siga[1:4,8])
  sd_sig_gp[2,4] <- mean(c(sd_siga[5,6:7],sd_siga[6,7]))
  sd_sig_gp[2,5] <- mean(sd_siga[5:7,8])
  sd_sig_gp[2,6] <- "Res"
  sd_sig_gp[2,7] <- i
  sd_sig_gp[3,1] <- mean(c(sd_sigb[1,2:4],sd_sigb[2,3:4],sd_sigb[3,4]))
  sd_sig_gp[3,2] <- mean(c(sd_sigb[1,5:7],sd_sigb[2,5:7],sd_sigb[3,5:7],sd_sigb[4,5:7]))
  sd_sig_gp[3,3] <- mean(sd_sigb[1:4,8])
  sd_sig_gp[3,4] <- mean(c(sd_sigb[5,6:7],sd_sigb[6,7]))
  sd_sig_gp[3,5] <- mean(sd_sigb[5:7,8])
  sd_sig_gp[3,6] <- "Res+t"
  sd_sig_gp[3,7] <- i
  colnames(sd_sig_gp) <- c("GE:GE", "GE:SP1", "GE:SP2","SP1:SP1","SP1:SP2","Model","Age")
  
  return(list(mean_pdi_gp,sd_pdi_gp, mean_sig_gp, sd_sig_gp))
}


MAIN_PROG <- function(name){
  test <- read.csv(name, header = F)
  test_cov <- COV

  test_cov[,1] = as.factor(test_cov[,1])
  test_cov[,3] = as.factor(test_cov[,3])
  test_cov[,5] = as.factor(test_cov[,5])
  test_cov[,6] = as.factor(test_cov[,6])

  test[,9] = test_cov[,4]
  test[,10] = test_cov[,6]
  test[,11] = test_cov[,2]
  test[,12] = test_cov[,7]

  colnames(test) <- c("GE_1","GE_2","GE_3","GE_4","SP1_1","SP1_2","SP1_3","SP2","Time","Iter","Age","Plan")
  colnames(test_cov) <- c("Disp","Age","Res","Time","Site","Iter","Plan")


  cov = subset(test_cov, 
             Time < 100 
             & Time > 20
  )
  n = subset(test, 
           Time < 100 
           & Time > 20
  )

  age <- 1:15
  RESULTS <- lapply(age,LAP_PDI,cov=cov,n=n)

  RES_PDI_GP <- lapply(RESULTS, "[[", 1)
  RES_PDI_GP_SD <- lapply(RESULTS, "[[", 2)
  MEAN_SIG_GP <- lapply(RESULTS, "[[", 3)
  SD_SIG_GP <- lapply(RESULTS, "[[", 4)

  return(list(do.call(rbind.data.frame, RES_PDI_GP),do.call(rbind.data.frame, RES_PDI_GP_SD), do.call(rbind.data.frame, MEAN_SIG_GP), do.call(rbind.data.frame, SD_SIG_GP)))
}
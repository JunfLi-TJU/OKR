setwd("D:/experiment/Conference Paper/ICML/ICML2023/code")
rm(list = ls())
library(MASS)

dpath          <- file.path("D:/experiment/online learning dataset/regression/") 

d_index <- 4
Dataset       <- c("elevators_all","bank_all", "Year_test","ailerons_all","calhousing","N-cpusmall",
                   "N-parkinsons","N-TomsHardware")              
savepath      <- paste0("D:/experiment/Conference Paper/ICML/ICML2023/result/",
                        paste0("NOGD-sq-all-",Dataset[d_index],".txt"))

savepath1       <- paste0("D:/experiment/Conference Paper/ICML/ICML2023/result/",
                          paste0("NOGD-sq-",Dataset[d_index],".txt"))

traindatapath  <- file.path(dpath, paste0(Dataset[d_index], ".train"))

traindatamatrix <- as.matrix(read.table(traindatapath))
trdata     <- traindatamatrix[ ,-1]
ylabel     <- traindatamatrix[ ,1]                                             

length_tr  <- nrow(trdata)                                               
feature_tr <- ncol(trdata)              

# -4 -3 -2 -1 0 1 2 3 4

sigma    <- 16
#B        <- 10*round(sqrt(length_tr)*log(log(length_tr))/feature_tr)
B        <- 600
r        <- 1*B

reptimes <- 10
runtime  <- c(rep(0, reptimes))
errorrate<- c(rep(0, reptimes))
eta      <- 10/sqrt(length_tr)

for(re in 1:reptimes)
{
  order    <- sample(1:length_tr,length_tr,replace = F)   #dis
  svpara   <- array(0,1)    # store the parameter for each support vector
  svmat    <- matrix(0,nrow = feature_tr,ncol=1)
  Gram     <- matrix(1,nrow = B,ncol=B)  
  
  t1       <- proc.time()  #proc.time()
  
  ### the first instance
  error      <- 1
  svmat[,1]  <- trdata[order[1], ]
  svpara[1]  <- 2*eta*ylabel[order[1]]
  k          <- 1
  
  for (t in 2:length_tr)
  {
    err <- 0
    if(k<B)
    {
      diff_S_i <- svmat[,1:k] - trdata[order[t], ]
      if(k>1)
      {
        tem   <- colSums(diff_S_i*diff_S_i)
      }else
      {
        tem   <- crossprod(diff_S_i,diff_S_i)[1,1]
      }
      kt      <- exp(tem/(-2*(sigma^2)))
      fx      <- crossprod(svpara[1:k],kt)[1,1]
      err     <- (fx-ylabel[order[t]])^2
      svmat       <- cbind(svmat,trdata[order[t],])
      
      svpara[k+1] <- -2*eta*(fx-ylabel[order[t]])
      k           <- k+1
      Gram[k,1:(k-1)] <- kt
      Gram[1:(k-1),k] <- kt
    }
    if(k==B)
    {
      T <- svd(Gram)
      U <- T$u
      D <- T$d  ##3 vector
      V <- T$v

      D <- D^{-0.5}
      
      tem <- D*V
      tem <- solve(tem)
      w_0 <- svpara%*%tem
      w_0 <- t(w_0)
      w_t <- w_0[1:r]
      
      Diag <- diag(D[1:r])
      V    <- t(V[,1:r])
      
#      w_t<-c(rep(0,r))
    }
    if(k>=B)
    {
      k   <- k+1
      
      diff_St <- svmat - trdata[order[t], ]
      tem     <- colSums(diff_St*diff_St)
      kt      <- exp(tem/(-2*(sigma^2)))
      tem     <- V%*%kt
      ph_t    <- D[1:r]*tem
      
      fx      <- crossprod(w_t,ph_t)[1,1]
      err     <- (fx-ylabel[order[t]])^2
      w_t   <- w_t-2*eta*(fx-ylabel[order[t]])*ph_t
    }
    error <- error + err   #record the err of selected superarm
  }
  t2 <- proc.time()
  runtime[re] <- (t2 - t1)[3]
  errorrate[re] <- error/length_tr
}

save_result <- list(
  note     = c(" the next term are:alg_name--dataname--sam_num--sigma--sv_num--gamma--eta--run_time--err_num--tot_run_time--ave_run_time--ave_err_rate--sd_time--sd_error"),
  alg_name = c("NOGD-sq"),
  dataname = paste0(Dataset[d_index], ".train"),
  sam_num  = length_tr,
  sv_num   = B,
  run_time = as.character(runtime),
  err_num  = as.character(errorrate), 
  tot_run_time = sum(runtime),
  ave_run_time = sum(runtime)/reptimes,
  ave_err_rate = sum(errorrate)/reptimes,
  sd_time  <- sd(runtime),
  sd_err    <-sd(errorrate)
)

write.table(save_result,file=savepath1,row.names =TRUE, col.names =FALSE, quote = F) 

sprintf("the candidate kernel parameter is %f", sigma)
sprintf("the number of sample is %d", length_tr)
sprintf("total training time is %.4f in dataset", sum(runtime))
sprintf("the average running time is %.5f in dataset", sum(runtime)/reptimes)
sprintf("the average error rate is %f", sum(errorrate)/reptimes)
sprintf("standard deviation of run_time is %.5f in dataset", sd(runtime))
sprintf("standard deviation of average error rate is %.5f in dataset", sd(errorrate))

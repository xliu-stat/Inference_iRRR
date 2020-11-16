require("MASS")
require("Matrix")
source("functions_inference.R")


proc.no <- 0 # setting 1: use the random seed sequence 0-19, each with 5 replications (100 runs in total)  
             # setting 2: use the random seed sequence 0-99, each with 1 replication (100 runs in total)
Rep <- 5
model <- 1 # 1 or 2
a <- 0.4 # SNR from (0.1, 0.2, 0.4)
rho_x <- 0.5 # rho in (0, 0.5)
lam1_XA_wiRRR <- 1 # for obtaining score matrix


# p < n 
if (model == 1){
  n <- 500
  q <- 5 
  p.vec <- c(10,10,10,10,10)
  r.vec <- c(2,0,0,0,0)
  view.num <- 5
  set_ind <- c(1,2,3)
}

# p > n
if (model == 2){
  n <- 200
  q <- 10 
  p.vec <- rep(20,20)
  r.vec <- c(1,rep(0,19))
  view.num <- 20
  set_ind <- c(1,2,3)
}

# tuning sequence
if (model == 2){
  if (rho_x == 0){
    lam1_scaled_wiRRR <- 10^seq(-0.75, 0.45, 0.05)
  } else {
    lam1_scaled_wiRRR <- 10^seq(-0.9, 0.3, 0.05)
  }
} else {
  lam1_scaled_wiRRR <- 10^seq(-1.3, 0.7, 0.05) 
}

# fix a random seed to get one uniform beta among all tasks
set.seed(1) 
# generate true coefficient matrix
# X: use function newdata2 to get correlation among groups
# X: use function newdata3 to get correlation only within each group
mydata <- newdata2(n,view.num,newB=TRUE,
                   SNR=a,p.vec,r.vec,q,
                   rho_X=rho_x,Sigma=CorrAR)
B <- mydata$B
cB <- mydata$cB



# estimation
err.sigma_seq <- rep(0, Rep)
uni_sigma_hat <- matrix(0, nrow = Rep, ncol = q)
uni_sigma_dif <- matrix(0, nrow = Rep, ncol = q)
uni_lam_est <- matrix(0, nrow = Rep, ncol = q)
uni_rankfit <- array(0, dim = c(Rep, q, view.num))
uni_flag_scaled_seq <- matrix(0, nrow = Rep, ncol = q)
uni_flag_B_seq <- matrix(0, nrow = Rep, ncol = q)
uni_flag_scaled_mat <- array(0, dim = c(Rep, q, length(lam1_scaled_wiRRR), 5))
uni_flag_B_mat <- array(0, dim = c(Rep, q, length(lam1_scaled_wiRRR), 5))
# inference
test_stat_mat <- array(0, dim = c(Rep, q, view.num))
p_val_mat <- array(0, dim = c(Rep, q, view.num))
cover_stat_mat <- array(0, dim = c(Rep, q, view.num))
cover_ind_mat <- array(0, dim = c(Rep, q, view.num))
# score matrix
d1 <- matrix(0, Rep, view.num)
flag_XA <- matrix(0, Rep, view.num)
P_rank <- matrix(0, Rep, view.num)




set.seed(proc.no)
fail <- 0
repnum <- 1


while(repnum <= Rep){
  tryCatch({
    ########################### Generate Data ##############################
    # generate design matrix, random error and response
    # identity error covariance, with sigma decided by SNR
    # X: use function newdata2 to get correlation among groups
    # X: use function newdata3 to get correlation only within each group
    mydata <- newdata2(n,view.num,newB=FALSE,
                       SNR=a,p.vec,r.vec,q,
                       rho_X=rho_x,Sigma=CorrAR)
    X <- mydata$X          
    cX <- mydata$cX
    SigmaX <- mydata$SigmaX
    err.sigma_seq[repnum] <- err.sigma <- mydata$SNR_sigma
    E <- matrix(rnorm(n*q, mean = 0, sd = err.sigma), nrow = n) 
    Y <- cX%*%cB + scale(E, center = TRUE, scale = FALSE)
    
    # score matrix estimation
    P_score <- array(0, dim = c(view.num, n, n))
    for (k in set_ind){
      Y_k <- X[[k]]
      X_k <- X[-k]
      cX_k <- X_k[[1]]
      if (view.num > 2){
        for (j in 2:(view.num-1)){
          cX_k <- cbind(cX_k,X_k[[j]])
        }
      }
      weight.k <- rep(0, view.num-1)
      for (i in 1:(view.num-1)){
        weight.k[i] <- (sqrt(1*ncol(X_k[[i]])) + sqrt(2*log(view.num)))/sqrt(1)/nrow(X_k[[i]])
      }
      fit_k <- iRRR_normal_XA(Y_k, X_k, weight=weight.k,
                              lam1 = 1, Tol = 1e-3,
                              Niter = 5000, varyrho = 0, rho = 0.1, lam0 = 0,
                              maxrho = 5, randomstart = 0)
      Z <- Y_k - fit_k$cXA
      P_z <- Z%*%solve(t(Z)%*%Z + 0.000001*diag(ncol(Z)))%*%t(Z)
      P_score[k,,] <- P_z
      P_rank[repnum,k] <- sum(svd(Z)$d > 0.01)
      d1[repnum,k] <- max(svd(P_z%*%(diag(n) - Y_k%*%solve(t(Y_k)%*%Y_k + 0.000001*diag(ncol(Y_k)) )%*%t(Y_k) ))$d)
      flag_XA[repnum,k] <- fit_k$flag_ADMM_XA
    }
    
    
    Y_all <- Y
    for (c in 1:q){
      Y <- as.matrix(Y_all[,c])
      fit3 <- iRRR_normal_CV_scaled_ADMM(Y, X, cX, p.vec, weight=NULL, 
                                                 lam.seq = lam1_scaled_wiRRR, nfold = 5,
                                                 norder = NULL, 
                                                 Tol = 1e-5, 
                                                 Niter = 5000, varyrho = 0, rho = 0.1, lam0 = 0,
                                                 maxrho = 5, randomstart = 0, new = 1)
      uni_sigma_hat[repnum, c] <- fit3$sigma_hat 
      uni_sigma_dif[repnum, c] <- fit3$sigma_hat/err.sigma
      uni_lam_est[repnum, c] <- fit3$lam.est 
      uni_rankfit[repnum, c,] <- fit3$r_hat_vec 
      
      
      uni_flag_scaled_seq[repnum, c] <- fit3$total_scaled_flag
      uni_flag_B_seq[repnum, c] <- fit3$total_B_flag
      uni_flag_scaled_mat[repnum, c, ,] <- fit3$flag_scaled_mat
      uni_flag_B_mat[repnum, c, ,] <- fit3$flag_B_mat
      
      B_init <- fit3$B_hat_list
      sigma_init <- fit3$sigma_hat
      
      
      # obtain de-biased estimator and test statistic
      B_debiased <- B_init
      test_stat <- rep(0, view.num)
      p_val <- rep(0, view.num)
      cover_stat <- rep(0, view.num)
      cover_ind <- rep(0, view.num)
      
      for (i in set_ind){
        B_debiased[[i]] <- B_init[[i]] - 
          ginv(P_score[i,,]%*%X[[i]])%*%P_score[i,,]%*%(Y-rep(1,nrow(Y))%*%fit3$mu_hat-cX%*%fit3$B_hat)
        # test statistic
        test_stat[i] <- sum((P_score[i,,]%*%(Y-rep(1,nrow(Y))%*%fit3$mu_hat
                                             -cX%*%fit3$B_hat+X[[i]]%*%B_init[[i]]))^2)/sigma_init^2
        # P-value (use chi-square distribution)
        p_val[i] <- 1-pchisq(test_stat[i], P_rank[repnum,i]*1)
        
        # verify the asymptotic distribution
        cover_stat[i] <- sum((P_score[i,,]%*%(Y-cX%*%fit3$B_hat+X[[i]]%*%B_init[[i]] -
                                                       X[[i]]%*%B[[i]][,c]))^2)/sigma_init^2
        cover_ind[i] <- (cover_stat[i] < qchisq(0.95, P_rank[repnum,i]*1))
      }
      test_stat_mat[repnum,c,] <- test_stat
      p_val_mat[repnum,c,] <- p_val
      cover_stat_mat[repnum,c,] <- cover_stat
      cover_ind_mat[repnum,c,] <- cover_ind
    }    
    
    repnum <- repnum + 1
  }, error=function(e) {
    fail <- fail + 1
    cat("At repnum:", repnum, "Error :",conditionMessage(e),"\n")})
}

save(n, q, set_ind, p.vec, r.vec, err.sigma_seq, view.num, Rep, a, rho_x,
     lam1_XA_wiRRR,
     uni_sigma_hat, uni_lam_est,
     uni_flag_scaled_seq,
     uni_flag_B_seq,
     uni_flag_scaled_mat,
     uni_flag_B_mat,
     uni_sigma_dif,
     uni_rankfit,
     p_val_mat, 
     P_rank, d1, flag_XA,
     test_stat_mat, cover_stat_mat, cover_ind_mat, 
     fail, file=paste("SNR_XALam",lam1_XA_wiRRR,
                      "_proc",proc.no,
                      "_model",model,"_rho",rho_x,"_a",a,"_5CV_allCorr_scaledADMM_uni.rda",sep=""))






require("MASS")
require("Matrix")
source("functions_inference.R")


proc.no <- 0 # for setting 1: use the random seed sequence 0-19, each with 5 replications (100 runs in total)
             # for setting 2: use the random seed sequence 0-99, each with 1 replication (100 runs in total) 
Rep <- 5
model <- 1 # 1 or 2
a <- 0.1 # SNR from (0.1, 0.2, 0.4)
rho_x <- 0 # rho in (0, 0.5)
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
    lam1_scaled_wiRRR <- 10^seq(-0.7, 0.5, 0.05)
  } else {
    lam1_scaled_wiRRR <- 10^seq(-1.2, 0, 0.05) 
  }
} else {
  lam1_scaled_wiRRR <- 10^seq(-2, 0, 0.04)
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



nmethod <- 2

err.sigma_seq <- rep(0, Rep)
sigma_est <- rep(0, Rep)
sigma_dif <- rep(0, Rep)
lambda_seq <- rep(0, Rep)

flag_scaled_seq <- rep(0, Rep)
flag_iRRR_seq <- rep(0, Rep)
flag_B_seq <- rep(0, Rep)
flag_whole_scaled_mat <- array(0, dim = c(Rep, length(lam1_scaled_wiRRR), 5))
flag_whole_iRRR_mat <- array(0, dim = c(Rep, length(lam1_scaled_wiRRR), 5))
flag_whole_B_mat <- array(0, dim = c(Rep, length(lam1_scaled_wiRRR), 5))


# for scaled iRRR, de-biased scaled iRRR
B_est <- array(0, dim = c(Rep, nmethod, sum(p.vec), q))
Esterror <- matrix(0, nrow = Rep, ncol = nmethod)
Prederror <- matrix(0, nrow = Rep, ncol = nmethod)
Rankfit <- array(0, dim = c(nmethod, view.num, Rep))
Timefit <- matrix(0, nrow = Rep, ncol = nmethod)

# for score matrix
P_rank <- matrix(0, nrow = Rep, ncol = view.num)
d1 <- matrix(0, nrow = Rep, ncol = view.num)
flag_XA <- matrix(0, nrow = Rep, ncol = view.num)


# for inference
test_stat <- matrix(0, nrow = Rep, ncol = view.num)
p_val <- matrix(0, nrow = Rep, ncol = view.num)

cover_stat <- matrix(0, nrow = Rep, ncol = view.num)
cover_ind <- matrix(0, nrow = Rep, ncol = view.num)



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
    
    
    #################### obtain score matrix ###################
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
      ### obtain weights for computing score matrix 
      weight.k <- rep(0, view.num-1)
      for (i in 1:(view.num-1)){ 
        weight.k[i] <- (sqrt(q*ncol(X_k[[i]])) + sqrt(2*log(view.num)))/sqrt(q)/nrow(X_k[[i]])
      }
      fit_k <- iRRR_normal_XA(Y_k, X_k, weight=weight.k,
                              lam1 = lam1_XA_wiRRR, Tol = 1e-3, 
                              Niter = 5000, varyrho = 0, rho = 0.1, lam0 = 0,
                              maxrho = 5, randomstart = 0)
      Z <- Y_k - fit_k$cXA
      P_z <- Z%*%solve(t(Z)%*%Z)%*%t(Z)
      P_score[k,,] <- P_z
      P_rank[repnum,k] <- sum(svd(Z)$d > 0.01)
      d1[repnum,k] <- max(svd(P_z%*%(diag(n) - Y_k%*%solve(t(Y_k)%*%Y_k)%*%t(Y_k) ))$d)
      flag_XA[repnum,k] <- fit_k$flag_ADMM_XA
    }
    
    #################### scaled iRRR ###################
    ptm <- system.time(
      fit2 <- iRRR_normal_CV_scaled(Y, X, cX, p.vec, weight, 
                                 lam.seq = lam1_scaled_wiRRR, nfold = 5,
                                 norder = NULL, 
                                 Tol = 1e-5, 
                                 Niter = 5000, varyrho = 0, rho = 0.1, lam0 = 0,
                                 maxrho = 5, randomstart = 0, new = 1))
    Timefit[repnum,1] <- ptm[1]
    
    B_est[repnum,1,,] <- fit2$B_hat
    Esterror[repnum,1] <- sqrt(sum((fit2$B_hat-cB)^2))
    Prederror[repnum,1] <- sum(diag(t(fit2$B_hat - cB)%*%SigmaX%*%(fit2$B_hat - cB)))
    
    sigma_est[repnum] <- fit2$sigma_hat
    sigma_dif[repnum] <- fit2$sigma_hat/err.sigma
    lambda_seq[repnum] <- fit2$lam.est
    Rankfit[1,,repnum] <- fit2$r_hat_vec
    
    flag_scaled_seq[repnum] <- fit2$total_scaled_flag
    flag_iRRR_seq[repnum] <- fit2$total_iRRR_flag
    flag_B_seq[repnum] <- fit2$total_B_flag
    
    flag_whole_scaled_mat[repnum,,] <- fit2$flag_scaled_mat
    flag_whole_iRRR_mat[repnum,,] <- fit2$flag_iRRR_mat
    flag_whole_B_mat[repnum,,] <- fit2$flag_B_mat
    
    
    ############################################################################
    ################################ inference #################################
    ############################################################################
    B_init <- fit2$B_hat_list
    sigma_init <- fit2$sigma_hat
    
    
    
    #################### obtain de-biased estimator #####################
    B_debiased <- B_init
    
    for (i in set_ind){
      B_debiased[[i]] <- B_init[[i]] - 
        ginv(P_score[i,,]%*%X[[i]])%*%P_score[i,,]%*%(Y-rep(1,nrow(Y))%*%fit2$mu_hat-cX%*%fit2$B_hat)
      # rank estimation
      Rankfit[2,i,repnum] <- sum(svd(B_debiased[[i]])$d>1e-2)
      # test statistic
      test_stat[repnum,i] <- sum((P_score[i,,]%*%(Y-rep(1,nrow(Y))%*%fit2$mu_hat
                                                  -cX%*%fit2$B_hat+X[[i]]%*%B_init[[i]]))^2)/sigma_init^2
      
      # P-value (use chi-square distribution)
      p_val[repnum,i] <- 1-pchisq(test_stat[repnum,i], P_rank[repnum,i]*q)
      
      # verify the asymptotic distribution
      cover_stat[repnum,i] <- sum((P_score[i,,]%*%(Y-cX%*%fit2$B_hat+X[[i]]%*%B_init[[i]] - 
                                                     X[[i]]%*%B[[i]]))^2)/sigma_init^2
      cover_ind[repnum,i] <- (cover_stat[repnum,i] < qchisq(0.95, P_rank[repnum,i]*q))
    }
    
    fit3_B <- B_debiased[[1]]
    if (view.num >= 2){
      for (j in 2:(view.num)){
        fit3_B <- rbind(fit3_B,B_debiased[[j]])
      }
    }
    
    B_est[repnum,2,,] <- fit3_B
    Esterror[repnum,2] <- sqrt(sum((fit3_B-cB)^2))
    Prederror[repnum,2] <- sum(diag(t(fit3_B - cB)%*%SigmaX%*%(fit3_B - cB)))
    
    
    repnum <- repnum + 1
  }, error=function(e) {
    fail <- fail + 1
    cat("At repnum:", repnum, "Error :",conditionMessage(e),"\n")})
}

save(n, q, set_ind, p.vec, r.vec, err.sigma_seq, view.num, Rep, a, rho_x,
     lam1_XA_wiRRR,lambda_seq,
     flag_scaled_seq,
     flag_iRRR_seq,
     flag_B_seq,
     flag_whole_scaled_mat,
     flag_whole_iRRR_mat,
     flag_whole_B_mat,
     sigma_est, sigma_dif, p_val, P_rank, d1, flag_XA,
     test_stat, cover_stat, cover_ind, B_est,
     Esterror, Prederror, 
     Rankfit, Timefit,
     fail, file=paste("SNR_XALam",lam1_XA_wiRRR,
                      "_proc",proc.no,
                      "_model",model,"_rho",rho_x,"_a",a,"_5CV_allCorr.rda",sep=""))


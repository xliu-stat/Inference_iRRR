#############################################################
#####                                                   #####
#####                                                   #####
#####                    Data Generation                #####
#####                                                   #####
#####                                                   #####
#############################################################



########### covariance matrix generation ###########
### This function generates covariance matrix    ###
### to generate predictors.                      ### 


# autoregressive structure                     
CorrAR <- function(p,rho){  
  Sigma <- matrix(nrow=p,ncol=p,NA)
  for(i in 1:p){
    for(j in 1:p){
      Sigma[i,j] <- rho^(abs(i-j))
    }
  }
  Sigma
}





# compound symmetry structure                  
CorrCS <- function(p,rho){  
  Sigma <- matrix(nrow=p,ncol=p,rho)
  diag(Sigma) <- 1
  Sigma
}
####################################################









############# simulation data generation #############
### This function generates normally distributed   ###
### predictors for simulation part 1:              ###



# correlation among all covariates               
newdata2 <- function(n=300,view.num=2,newB=FALSE,
                     SNR=1,
                     p.vec=c(50,50),r.vec=c(10,10),
                     q=50,rho_X=0.5,Sigma=CorrAR){
  require(MASS)
  
  # input:
  # n: sample size
  # q: number of responses
  # view.num: view number
  # p.vec: a vector to specify group size
  # r.vec: a vector to specify group-specific rank
  # rho_X: correlation strength to specify covariance matrix
  # Sigma: covariance structure CorrAR: autoregressive and CorrCS: compound symmetry
  # newB: TRUE: to generate new coefficient matrix, FALSE: do not generate new coefficient matrix
  # SNR: signal to noise ratio
  
  # output:
  # B: a list, the generated true coefficient matrix
  # cB: a matrix, the generated true coefficient matrix
  # X: a list, the generated multi-view data 
  # cX: a matrix, the generated integrated design matrix
  # SNR_sigma: a scalar, the standard error of the random noise
  # weight: a vector, the weights used in scaled iRRR
  # XhfSigma: the square root of the covariance matrix of X
  
  
  # whether or not generate new B
  if (newB){
    B <- list()
    for (i in 1:view.num){
      if (r.vec[i] == 0){
        B[[i]] <- matrix(0, nrow = p.vec[i], ncol = q)
      }else{
        L <- matrix(rnorm(p.vec[i]*r.vec[i]), nrow = p.vec[i], ncol = r.vec[i])
        R <- matrix(rnorm(q*r.vec[i]), nrow = q, ncol = r.vec[i])
        B[[i]] <- L%*%t(R)
      }
    }
    
    cB <- B[[1]]
    for (i in 2:view.num){
      cB <- rbind(cB,B[[i]])
    }
  }
  
  #  generate design matrix
  SigmaX <- Sigma(sum(p.vec), rho_X)
  dec <- svd(SigmaX)
  XhfSigma <- dec$u%*%diag(sqrt(dec$d))%*%t(dec$u)
  cX <- matrix(rnorm(n*sum(p.vec)), nrow = n)%*%XhfSigma;
  cX <- scale(cX, center = TRUE, scale = TRUE)*sqrt(n/(n-1))
  X <- list()
  K <- view.num
  tempp1 <- 1+cumsum(c(0,p.vec[-K]))
  tempp2 <- cumsum(p.vec)
  for (i in 1:K){
    X[[i]] <- cX[,c(tempp1[i]:tempp2[i])]
    X[[i]] <- sqrt(n)*X[[i]]/max(svd(X[[i]])$d)
  }
  cX <- X[[1]]
  for (i in 2:K){
    cX <- cbind(cX, X[[i]])
  }
  
  B_max <- max(abs(cB))
  cB <- cB/B_max
  B <- list()
  for (i in 1:view.num){
    B[[i]] <- cB[c(tempp1[i]:tempp2[i]),]
  }
  
  temp <- cX%*%cB
  SNR_sigma <- sd(as.vector(temp))/SNR
  
  # weight for scaled iRRR
  weight <- rep(0, K)
  for (i in 1:K){ 
    weight[i] <- max(svd(X[[i]])$d)*(sqrt(q*p.vec[i]) +  sqrt(2*log(K)))/n
  }
  
  return(list(B=B, cB=cB, X=X, cX=cX,
              SNR_sigma=SNR_sigma,
              weight=weight, 
              SigmaX=SigmaX,
              XhfSigma=XhfSigma))
}





# correlation only within each group             
newdata3 <- function(n=300,view.num=2,newB=FALSE,
                     SNR=1,
                     p.vec=c(50,50),r.vec=c(10,10),
                     q=50,rho_X=0.5,Sigma=CorrAR){
  require(MASS)
  
  # input:
  # n: sample size
  # q: number of responses
  # view.num: view number
  # p.vec: a vector to specify group size
  # r.vec: a vector to specify group-specific rank
  # rho_X: correlation strength to specify covariance matrix
  # Sigma: covariance structure CorrAR: autoregressive and CorrCS: compound symmetry
  # newB: TRUE: to generate new coefficient matrix, FALSE: do not generate new coefficient matrix
  # SNR: signal to noise ratio
  
  # output:
  # B: a list, the generated true coefficient matrix
  # cB: a matrix, the generated true coefficient matrix
  # X: a list, the generated multi-view data 
  # cX: a matrix, the generated integrated design matrix
  # SNR_sigma: a scalar, the standard error of the random noise
  # weight: a vector, the weights used in scaled iRRR
  # SigmaX: the covariance matrix of X
  # XhfSigma: the square root of the covariance matrix of X
  
  
  # whether or not generate new B
  if (newB){
    B <- list()
    for (i in 1:view.num){
      if (r.vec[i] == 0){
        B[[i]] <- matrix(0, nrow = p.vec[i], ncol = q)
      }else{
        L <- matrix(rnorm(p.vec[i]*r.vec[i]), nrow = p.vec[i], ncol = r.vec[i])
        R <- matrix(rnorm(q*r.vec[i]), nrow = q, ncol = r.vec[i])
        B[[i]] <- L%*%t(R)
      }
    }
    
    cB <- B[[1]]
    for (i in 2:view.num){
      cB <- rbind(cB,B[[i]])
    }
  }
  
  
  #  generate design matrix
  SigmaX <- kronecker(diag(view.num), Sigma(p.vec[1], rho_X)) 
  dec <- svd(SigmaX)
  XhfSigma <- dec$u%*%diag(sqrt(dec$d))%*%t(dec$u)
  cX <- matrix(rnorm(n*sum(p.vec)), nrow = n)%*%XhfSigma;
  cX <- scale(cX, center = TRUE, scale = TRUE)*sqrt(n/(n-1))
  X <- list()
  K <- view.num
  tempp1 <- 1+cumsum(c(0,p.vec[-K]))
  tempp2 <- cumsum(p.vec)
  for (i in 1:K){
    X[[i]] <- cX[,c(tempp1[i]:tempp2[i])]
    X[[i]] <- sqrt(n)*X[[i]]/max(svd(X[[i]])$d)
  }
  cX <- X[[1]]
  for (i in 2:K){
    cX <- cbind(cX, X[[i]])
  }
  
  B_max <- max(abs(cB))
  cB <- cB/B_max
  B <- list()
  for (i in 1:view.num){
    B[[i]] <- cB[c(tempp1[i]:tempp2[i]),]
  }
  
  temp <- cX%*%cB
  SNR_sigma <- sd(as.vector(temp))/SNR
  
  # weight for scaled iRRR
  weight <- rep(0, K)
  for (i in 1:K){ 
    weight[i] <- max(svd(X[[i]])$d)*(sqrt(q*p.vec[i]) +  sqrt(2*log(K)))/n
  }
  
  return(list(B=B, cB=cB, X=X, cX=cX, 
              SNR_sigma=SNR_sigma,
              weight=weight, 
              SigmaX=SigmaX,
              XhfSigma=XhfSigma))
}

####################################################








############# simulation data generation #############
### This function generates compositional          ###
### predictors for simulation part 2.              ###

compdata <- function(n=300,view.num=2,newB=FALSE,use.SNR=TRUE,fix.sigma=1,
                     SNR=1, meanX=rep(1,p),
                     p.vec=c(50,50),r.vec=c(10,10),
                     q=50,rho_X=0.5,Sigma=CorrAR){
  require(MASS)
  require(mvtnorm)
  
  # input:
  # n: sample size
  # q: number of responses
  # view.num: view number
  # p.vec: a vector to specify group size
  # r.vec: a vector to specify group-specific rank
  # meanX: mean vector in multivariate normal distribution
  # rho_X: correlation strength to specify covariance matrix of multivariate normal distribution
  # Sigma: covariance structure CorrAR: autoregressive and CorrCS: compound symmetry
  # newB: TRUE: to generate new coefficient matrix, FALSE: do not generate new coefficient matrix
  # use.SNR: TRUE: use SNR to find sigma, FALSE: use fixed sigma
  # SNR: signal to noise ratio
  
  # output:
  # B: a list, the generated true coefficient matrix
  # cB: a matrix, the generated true coefficient matrix
  # X: a list, the generated multi-view data 
  # cX: a matrix, the generated integrated design matrix
  # temp: a matrix, the true response
  # Y: a matrix, the noise polluted response
  # matW: a matrix, the generated count data
  # listW: a list, the generated count data
  # SNR_sigma: a scalar, the standard error of the random noise
  # weight: a vector, the weights used in scaled iRRR

  
  K <- view.num
  tempp1 <- 1+cumsum(c(0,p.vec[-K]))
  tempp2 <- cumsum(p.vec)
  
  # whether or not generate new B
  if (newB){
    B <- list()
    for (i in 1:view.num){
      if (r.vec[i] == 0){
        B[[i]] <- matrix(0, nrow = p.vec[i], ncol = q)
      }else{
        L <- matrix(rnorm(p.vec[i]*r.vec[i], mean = 0), nrow = p.vec[i], ncol = r.vec[i])
        R <- matrix(rnorm(q*r.vec[i], mean = 0), nrow = q, ncol = r.vec[i])
        B[[i]] <- L%*%t(R)
      }
    }
    
    cB <- B[[1]]
    for (i in 2:view.num){
      cB <- rbind(cB,B[[i]])
    }
  }
  
  # generate counts from LogNormal distribution 
  SigmaW <- Sigma(sum(p.vec), rho_X)
  matW <- rmvnorm(n, mean = meanX, sigma = SigmaW)
  matW <- exp(matW)
  
  # the list of count
  listW <- list()
  # log-subcomp list
  listZ <- list()
  for (i in 1:K){
    listW[[i]] <- matW[,c(tempp1[i]:tempp2[i])]
    # each log-subcomp
    listZ[[i]] <- log(listW[[i]]/(apply(listW[[i]], 1, sum)%*%t(rep(1,ncol(listW[[i]])))))
  }
  # log-subcomp mat
  cZ <- listZ[[1]]
  for (i in 2:K){
    cZ <- cbind(cZ, listZ[[i]])
  }
  
  # linear transformed design
  X <- list()
  for (i in 1:K){
    p1 <- p.vec[i]
    X[[i]] <- listZ[[i]]%*%(diag(p1) - rep(1,p1)%*%t(rep(1,p1))/p1)
  }
  cX <- X[[1]]
  for (i in 2:K){
    cX <- cbind(cX, X[[i]])
  }
  
  # true Y
  temp <- cX%*%cB
  # noise level
  if(use.SNR==TRUE){
    SNR_sigma <- sd(as.vector(temp))/SNR
    E <- matrix(rnorm(n*q, mean = 0, sd = SNR_sigma), nrow = n)
    Y <- temp + scale(E, center = TRUE, scale = FALSE)
  }else{
    SNR_sigma <- fix.sigma
    E <- matrix(rnorm(n*q, mean = 0, sd = SNR_sigma), nrow = n)
    Y <- temp + scale(E, center = TRUE, scale = FALSE)
  }
  
  
  # weight for scaled iRRR
  weight <- rep(0, K)
  for (i in 1:K){ 
    weight[i] <- max(svd(X[[i]])$d)*(sqrt(q*p.vec[i]) +  sqrt(2*log(K)))/n
  }
  
  return(list(B=B, cB=cB, X=X, cX=cX, temp=temp, Y=Y,
              matW=matW, listW=listW,
              SNR_sigma=SNR_sigma,
              weight=weight))
}

####################################################








#############################################################
#####                                                   #####
#####                                                   #####
#####                    Data Estimation                #####
#####                                                   #####
#####                                                   #####
#############################################################




#####################   iRRR estimation    ##################
### This function solves the iRRR problem with penalty:   ###
### \sum_k w_k \|B_k\|_* by ADMM algorithm.               ###


# ADMM 
iRRR_normal <- function(Y, X, weight = NULL, lam1, Tol = 1e-3, 
                        Niter = 500, varyrho = 0, rho = 0.1, lam0 = 0,
                        maxrho = 5, randomstart = 0){
  
  # input:
  # X: a list of multi-view data
  # Y: a matrix of response
  # weight: a vector of weights used in iRRR, if NULL then generate it 
  # lam1: the tuning parameter in iRRR
  # Tol: the threshold value of ADMM
  # Niter: the maximum number of iteration in ADMM
  # veryrho: 1 increase rho gradually, 0 do not change rho
  # rho: step size in ADMM
  # lam0: the tuning parameter for ridge penalty
  # maxrho: the largest possible rho if varyrho=1
  # randomstart: whether randomly generate initial values of ADMM
  
  # output:
  # C: the estimated coefficient matrix 
  # mu: the estimated intercept vector
  # A: a list, each list item is an estimated sub-coefficient matrix
  # flag_ADMM: 0 ADMM converges, 1 ADMM does not converge after Niter many iterations
  # flag_B: 0 the estimated matrix is not a zero matrix, 1 it is a zero matrix
  # rec_obj_all: a vector of full objective function value (with penalties)
  # rec_obj_ls_all: a vector of partial objective function value (no penalties)
  # rec_primal: a vector of primal residuals
  # rec_dual: a vector of dual residuals
  # rec_diff: a vector of max(primal residual, dual residual)

  K <- length(X) # a list
  p <- rep(0, K) # a vector of size of each view
  Y <- as.matrix(Y)
  n <- nrow(Y)
  q <- ncol(Y)
  
  for (i in 1:K){
    p[i] <- ncol(X[[i]])
  }
  
  # default weights
  if (is.null(weight)){
    weight <- rep(0,K)
    for (i in 1:K){
      weight[i] <- max(svd(X[[i]])$d)*(sqrt(q) + sqrt(sum(svd(X[[i]])$d>0)))/n
    }
  }
  
  
  # horizonatally concatenated X, also column centered and weighted
  cX <- matrix(0, nrow = n, ncol = sum(p))
  meanX <- vector()
  
  for (i in 1:K){
    meanX <- c(meanX, apply(X[[i]], 2, mean)) # collect all means
    # first, do column centering
    X[[i]] <- scale(X[[i]], center = TRUE, scale = FALSE) 
    # second, normalize centered X[[i]]'s by wi's
    X[[i]] <- X[[i]]/weight[i] 
    cX[,c((sum(p[1:(i-1)])*(i != 1)+1):sum(p[1:i]))] <- X[[i]] # column centered X
  }
  
  # initial parameter estimates
  mu <- apply(Y, 2, mean, na.rm = TRUE) # a q*1 vector 
  # majorize Y to get a working Y
  wY <- Y
  temp <- rep(1, n)%*%t(mu)
  wY[is.na(wY)] <- temp[is.na(wY)] # wY should be a complete matrix
  mu <- apply(wY, 2, mean) # new est of mu, b/c cX is col centered
  wY1 <- scale(wY, center = TRUE, scale = FALSE) # column centered wY
  
  B <- list()    
  Theta <- list()   # Lagrange params for B
  cB <- matrix(0, nrow = sum(p), ncol = q) # vertically concatenated B
  for (i in 1:K){
    if (randomstart){
      B[[i]] <- matrix(rnorm(p[i]*q), nrow = p[i]) # initialize the algorithm
    } else {
      B[[i]] <- ginv(t(X[[i]])%*%X[[i]])%*%t(X[[i]])%*%wY1 
    }
    Theta[[i]] <- matrix(0, nrow = p[i], ncol = q)
    cB[c((sum(p[1:(i-1)])*(i != 1) + 1):sum(p[1:i])),] <- B[[i]]
  }
  A <- B # low-rank alias
  cA <- cB
  cTheta <- matrix(0, nrow = sum(p), ncol = q) 
  
  decomp <- svd(cX/sqrt(n), nv = min(n, sum(p)))
  if (varyrho == 0){
    # compute inv(1/nX'X+(lam0+rho)I)
    DeltaMat <- decomp$v%*%diag(1/(decomp$d^2 + lam0 + rho))%*%t(decomp$v) + 
      (diag(rep(1,sum(p)))-decomp$v%*%t(decomp$v))/(lam0 + rho)
  }
  
  # check obj value
  # full objective function (with penalties) on observed data
  obj <- ObjValue1(Y, X, mu, A, lam0, lam1) 
  # only the least square part on observed data
  obj_ls <- ObjValue1(Y, X, mu, A, 0, 0)
  
  
  #################
  # ADMM
  flag_ADMM <- 0
  niter <- 0
  diff <- 1e+10
  rec_obj_all <- obj  # record obj value
  rec_obj_ls_all <- obj_ls # record obj_ls value
  rec_Theta <- vector() # record the Fro norm of Theta
  rec_primal <- vector() # record primal residual
  rec_dual <- vector() # record dual residual
  rec_diff <- vector() # record diff
  
  while ((niter < Niter)&(abs(diff) > Tol)){
    niter <- niter + 1
    cB_old <- cB
    
    ########### majorization #########################
    Eta <- rep(1, n)%*%t(mu) + cX%*%cB  # current linear predictor
    wY <- Y
    wY[is.na(wY)] <- Eta[is.na(wY)]
    mu <- apply(wY, 2, mean)
    wY1 <- scale(wY, center = TRUE, scale = FALSE)
    
    # est concatenated B
    if (varyrho) {
      DeltaMat <- decomp$v%*%diag(1/(decomp$d^2 + lam0 + rho))%*%t(decomp$v) + 
        (diag(rep(1,sum(p)))-decomp$v%*%t(decomp$v))/(lam0 + rho)
    }
    cB <- DeltaMat%*%(t(cX)%*%wY1/n+rho*cA+cTheta)
    for (j in 1:K){
      B[[j]] <- cB[c((sum(p[1:(j-1)])*(j != 1)+1):sum(p[1:j])),]
    }
    
    # est A_i 
    # update Theta_i right after est A_i
    for (k in 1:K){
      # est A
      temp <- B[[k]] - Theta[[k]]/rho
      decomp1 <- svd(temp, nu = min(q, p[k]), nv = min(q, p[k]))
      if (length(decomp1$d)>1){
        A[[k]] <- decomp1$u%*%diag(SoftThres(decomp1$d, lam1/rho))%*%t(decomp1$v)
      }else{
        A[[k]] <- decomp1$u%*%max(decomp1$d-lam1/rho,0)%*%t(decomp1$v)
      }
      
      # update Theta
      Theta[[k]] <- Theta[[k]] + rho*(A[[k]]-B[[k]])
      
      # update cA and cTheta
      cA[c((sum(p[1:(k-1)])*(k != 1)+1):sum(p[1:k])),] <- A[[k]]
      cTheta[c((sum(p[1:(k-1)])*(k != 1)+1):sum(p[1:k])),] <- Theta[[k]]
    }
    
    # update rho
    if (varyrho) { rho <- min(maxrho, 1.01*rho) } # steadily increasing rho
    
    # stopping rule
    # primal and dual residuals
    primal <- sum((cA-cB)^2)
    primal <- sqrt(primal)
    rec_primal <- c(rec_primal, primal)
    dual <- sum((cB-cB_old)^2)
    dual <- rho*sqrt(dual)
    rec_dual <- c(rec_dual, dual)
    # Fro norm of Theta
    rec_Theta <- c(rec_Theta, sqrt(sum(cTheta^2)))
    
    # objective function value
    obj <- ObjValue1(Y,X,mu,A,lam0,lam1)
    obj_ls <- ObjValue1(Y,X,mu,A,0,0)
    rec_obj_all <- c(rec_obj_all, obj)
    rec_obj_ls_all <- c(rec_obj_ls_all, obj_ls)
    
    # stopping rule
    diff <- max(primal, dual)
    rec_diff <- c(rec_diff,diff)
  }
  
  if (niter == Niter){
    flag_ADMM <- 1
    cat("iRRR does NOT converge after_", niter, "_iterations!", "\n", sep = "")
  } else {
    cat("iRRR converges after_", niter, "_iterations.", "\n", sep = "")
  }
  
  
  # output
  # rescale parameter estimate, and add back mean
  C <- matrix(0, nrow = sum(p), ncol = q)
  for (i in 1:K){
    A[[i]] <- A[[i]]/weight[i]
    B[[i]] <- B[[i]]/weight[i]
    C[c((sum(p[1:(i-1)])*(i != 1)+1):sum(p[1:i])),] <- A[[i]]
  }
  mu <- mu-(t(meanX)%*%C)
  
  flag_B <- 0
  if (sum(C^2)==0) flag_B <- 1
  
  return(list(C=C, mu=mu, A=A, B=B, Theta=Theta, 
              flag_ADMM=flag_ADMM, flag_B=flag_B,
              rec_obj_all=rec_obj_all, rec_obj_ls_all=rec_obj_ls_all,
              rec_primal=rec_primal, rec_dual=rec_dual,
              rec_Theta=rec_Theta, rec_diff=rec_diff))
}

#############################################################





# functions used in ADMM
SoftThres <- function(Din, lam){
  # this function soft thresholds the diagonal values of Din
  # Din is a diagonal matrix
  # lam is a positive threshold
  # output is also a diagonal matrix
  d <- diag(Din)
  d[d>0] <- pmax(d[d>0]-lam,0)
  d[d<0] <- pmin(d[d<0]+lam,0)
  Dout <- diag(d)
  return(Dout)
}

ObjValue1 <- function(Y, X, mu, B, lam0, lam1){
  # Calc 1/(2n)|Y-1*mu'-sum(Xi*Bi)|^2 + lam0/2*sum(|Bi|_F^2) + lam1*sum(|Bi|_*) 
  # with column centered Xi's and (potentially non-centered and missing) Y
  n <- nrow(Y)
  K <- length(X)
  obj <- 0
  pred <- rep(1, n)%*%t(mu)
  for (i in 1:K){
    pred <- pred + X[[i]]%*%B[[i]]
    obj <- obj + lam0/2*sum((B[[i]])^2) + lam1*sum(svd(B[[i]])$d)
  }
  obj <- obj + 1/2/n*sum((Y-pred)^2, na.rm = TRUE)
  return(obj)
}

#############################################################





# tuning by cross validation (the weights are for scaled iRRR)
iRRR_normal_CV <- function(Y, X, cX, p.vec, weight=NULL, lam.seq, nfold, norder = NULL, Tol = 1e-3, 
                           Niter = 500, varyrho = 0, rho = 0.1, lam0 = 0,
                           maxrho = 5, randomstart = 0){
  
  # input:
  # X: a list of multi-view data
  # cX: the integrated design matrix
  # Y: a matrix of response
  # p.vec: a vector of group sizes
  # weight: a vector of weights used in iRRR, if NULL then generate it 
  # lam.seq: the tuning sequence
  # nfold: used in cross validation
  # norder: used in cross validation, if NULL then randomly generate it
  # Tol: the threshold value of ADMM
  # Niter: the maximum number of iteration in ADMM
  # veryrho: 1 increase rho gradually, 0 do not change rho
  # rho: step size in ADMM
  # lam0: the tuning parameter for ridge penalty
  # maxrho: the largest possible rho if varyrho=1
  # randomstart: whether randomly generate initial values of ADMM
  
  # output:
  # B_hat: the estimated coefficient matrix 
  # mu_hat: the estimated intercept vector
  # B_hat_list: a list, each list item is an estimated sub-coefficient matrix
  # r_hat_vec: a vector of estimated ranks
  # lam.est: the selected tuning parameter 
  # norder, lam_path: details in cross validation 
  # flag_iRRR_mat: a matrix of the same dimension as lam_path, 
  #                whether or not iRRR converges in each cross validation fit
  # total_iRRR_flag: >0: some iRRR do not converge, 0: all iRRR converge
  # flag_B_mat: a matrix of the same dimension as lam_path,
  #             whether or not the estimated B is zero in each cross validation fit
  # total_B_flag: >0: some fitted B are zero matrices, 0: all fitted B are not zero matrices
  
  
  Y <- as.matrix(Y)
  n <- nrow(Y)
  p <- ncol(cX)
  q <- ncol(Y)
  K <- length(X) # number of views
  
  if (is.null(weight)){
    weight <- rep(0,K)
    for (i in 1:K){
      weight[i] <- max(svd(X[[i]])$d)*(sqrt(q) + sqrt(sum(svd(X[[i]])$d>0)))/n
    }
  }
  
  if (is.null(norder))
    norder <- sample(seq_len(n),n)
  
  lam_path <- matrix(ncol = nfold, nrow = length(lam.seq), NA)
  # iRRR convergence
  flag_iRRR_mat <- matrix(ncol = nfold, nrow = length(lam.seq), NA)
  # whether B is zero or not
  flag_B_mat <- matrix(ncol = nfold, nrow = length(lam.seq), NA)
  
  ndel <- round(n/nfold)
  for (f in seq_len(nfold)){
    if (f != nfold) {
      iddel <- norder[(1 + ndel * (f - 1)):(ndel * f)]
    }
    else {
      iddel <- norder[(1 + ndel * (f - 1)):n]
    }
    ndel <- length(iddel)
    nf <- n - ndel
    idkeep <- (seq_len(n))[-iddel]
    
    Xf <- cX[-iddel, ]
    Xfdel <- cX[iddel, ]
    Yf <- as.matrix(Y[-iddel, ])
    Yfdel <- as.matrix(Y[iddel, ])
    
    Xf.list <- list()
    tempp1 <- 1+cumsum(c(0,p.vec[-K]))
    tempp2 <- cumsum(p.vec)
    for (i in 1:K){
      Xf.list[[i]] <- Xf[,c(tempp1[i]:tempp2[i])]
    }
    
    weight.f <- rep(0,K)
    for (i in 1:K){
      weight.f[i] <- max(svd(Xf.list[[i]])$d)*(sqrt(q) + 
                                                 sqrt(sum(svd(Xf.list[[i]])$d>0)))/nrow(Xf.list[[i]])
    }
    
    
    for (i in 1:length(lam.seq)){
      fit <- iRRR_normal(Yf, Xf.list, weight = weight.f, lam1 = lam.seq[i], Tol, 
                         Niter, varyrho, rho, lam0,
                         maxrho, randomstart)
      lam_path[i,f] <- sum((Yfdel-rep(1,nrow(Yfdel))%*%fit$mu-Xfdel%*%fit$C)^2)
      # iRRR convergence 
      flag_iRRR_mat[i,f] <- fit$flag_ADMM
      # whether B is zero or not
      flag_B_mat[i,f] <- fit$flag_B
    }
  }
  index <- order(colSums(lam_path))
  crerr <- rowSums(lam_path[, index])/length(index) * nfold
  lam.est <- lam.seq[which.min(crerr)]
  
  # use all data to get estimate of B with the selected lambda
  fit <- iRRR_normal(Y, X, weight = weight, lam1 = lam.est, Tol, 
                     Niter, varyrho, rho, lam0,
                     maxrho, randomstart)
  B <- fit$C # estimated matrix
  A <- fit$A # estimated list
  mu <- fit$mu
  r_hat_vec <- rep(0,K)
  for (i in 1:K){
    r_hat_vec[i] <- sum(svd(A[[i]])$d>1e-2)
  }
  
  total_iRRR_flag <- sum(flag_iRRR_mat)
  total_B_flag <- sum(flag_B_mat)
  
  return(list(B_hat=B, B_hat_list=A, mu_hat=mu, r_hat_vec=r_hat_vec, 
              norder=norder, lam.est=lam.est, lam_path=lam_path,
              flag_iRRR_mat=flag_iRRR_mat, flag_B_mat=flag_B_mat,
              total_iRRR_flag=total_iRRR_flag,
              total_B_flag=total_B_flag))
}

#############################################################







################## scaled iRRR estimation #################
### This function solves the scaled iRRR problem        ###
### by Blockwise coordinate descent algorithm.          ###


# Blockwise coordinate descent algorithm
iRRR_normal_scaled <- function(Y, X, cX, p.vec, weight=NULL, lam1, Tol = 1e-3, 
                               Niter = 500, varyrho = 0, rho = 0.1, lam0 = 0,
                               maxrho = 5, randomstart = 0){
  
  
  # input:
  # X: a list of multi-view data
  # cX: an integrated design matrix
  # Y: a matrix of response
  # p.vec: a vector of group sizes
  # weight: a vector of weights used in scaled iRRR, if NULL then generate it 
  # lam1: the tuning parameter in iRRR
  # Tol: the threshold value of ADMM
  # Niter: the maximum number of iteration in ADMM
  # veryrho: 1 increase rho gradually, 0 do not change rho
  # rho: step size in ADMM
  # lam0: the tuning parameter for ridge penalty
  # maxrho: the largest possible rho if varyrho=1
  # randomstart: whether randomly generate initial values of ADMM
  
  # output:
  # B_hat: the estimated coefficient matrix 
  # sigma_hat: the estimated noise level
  # mu_hat: the estimated intercept vector
  # B_hat_list: a list, each list item is an estimated sub-coefficient matrix
  # r_hat_vec: a vector of estimated ranks
  # iRRR_flag: >0: some iRRR do not converge, 0: all iRRR converge
  # conv_flag: 0 the iterative algorithm for scaled iRRR converges, 1 does not converge
  # flag_B: 0 the estimated matrix is not a zero matrix, 1 it is a zero matrix
  # obj: a vector of full objective function value (with penalties)
  
  
  Y <- as.matrix(Y)
  n <- nrow(Y)
  p <- ncol(cX)
  q <- ncol(Y)
  K <- length(X)         # number of views
  
  # compute weight
  if(is.null(weight)){
    weight <- rep(0, K)
    for (i in 1:K){ 
      weight[i] <- max(svd(X[[i]])$d)*(sqrt(q*p.vec[i]) +  sqrt(2*log(K)))/n
    }
  }
  
  # the iterative algorithm
  max_iter <- 100
  sigma_tol_val <- 1e-4  
  iter_num <- 0
  
  # initialization
  iRRR_flag <- 0 # >0: some iRRR do not converge, 0: all iRRR converge
  flag <- 0
  sigma_init <- 0.1
  sigma_new <- 1
  sigma_dif <- sigma_new - sigma_init
  obj1 <- vector()
  sigma_new_list <- vector()
  sumnu <- vector()
  
  while ((iter_num < max_iter)&(abs(sigma_dif) > sigma_tol_val)){
    iter_num <- iter_num + 1
    sigma_init <- sigma_new
    weight1 <- weight*sigma_init
    fit <- iRRR_normal(Y, X, weight = weight1, lam1, Tol, 
                       Niter, varyrho, rho, lam0,
                       maxrho, randomstart)
    iRRR_flag <- iRRR_flag + fit$flag_ADMM
    new_B <- fit$C
    sum_nu <- 0
    for (g in 1:K){
      sum_nu <- sum_nu + weight[g]*sum(svd(fit$B[[g]])$d)
    }
    sumnu <- c(sumnu,sum_nu)
    new_mu <- fit$mu
    sigma_new <- sqrt(sum( (Y - rep(1,nrow(Y))%*%new_mu - cX%*%new_B)^2 ) )/sqrt(n*q)
    sigma_new_list <- c(sigma_new_list,sigma_new)
    sigma_dif <- sigma_new - sigma_init
    obj1_new <- sigma_new^2/2/sigma_init+sigma_init/2+lam1*sum_nu/q
    obj1 <- c(obj1, obj1_new)
  }
  
  if (iter_num == max_iter){
    flag <- 1
    cat("scaled iRRR does NOT converge after_", iter_num, "_iterations!", "\n", sep = "")
  } else {
    cat("scaled iRRR converges after_", iter_num, "_iterations.", "\n", sep = "")
  }
  
  
  mu <- fit$mu
  A <- fit$A
  r_hat_vec <- rep(0,K)
  for (i in 1:K){
    r_hat_vec[i] <- sum(svd(A[[i]])$d>1e-2)
  }
  
  flag_B <- 0
  if (sum(new_B^2)==0) flag_B <- 1
  
  return(list(B_hat=new_B, sigma_hat=sigma_new, B_hat_list=A, mu_hat=mu, flag_B=flag_B,
              r_hat_vec=r_hat_vec, conv.flag=flag, iRRR_flag=iRRR_flag, obj=obj1,
              sigma_new_list=sigma_new_list,
              sumnu_list=sumnu))
}

#############################################################





# tuning by cross validation
iRRR_normal_CV_scaled <- function(Y, X, cX, p.vec, weight=NULL, lam.seq, nfold, norder = NULL, Tol = 1e-3, 
                                  Niter = 500, varyrho = 0, rho = 0.1, lam0 = 0,
                                  maxrho = 5, randomstart = 0, new = 0){
  
  # input:
  # X: a list of multi-view data
  # cX: the integrated design matrix
  # Y: a matrix of response
  # p.vec: a vector of group sizes
  # weight: a vector of weights used in scaled iRRR, if NULL then generate it 
  # lam.seq: the tuning sequence
  # nfold: used in cross validation
  # norder: used in cross validation, if NULL then randomly generate it
  # Tol: the threshold value of ADMM
  # Niter: the maximum number of iteration in ADMM
  # veryrho: 1 increase rho gradually, 0 do not change rho
  # rho: step size in ADMM
  # lam0: the tuning parameter for ridge penalty
  # maxrho: the largest possible rho if varyrho=1
  # randomstart: whether randomly generate initial values of ADMM
  # new: whether or not use the negative log-likelihood function as an error measurement in cross validation
  #      0 do not use, 1 use
  
  # output:
  # B_hat: the estimated coefficient matrix
  # sigma_hat: the estimated noise level
  # mu_hat: the estimated intercept vector
  # B_hat_list: a list, each list item is an estimated sub-coefficient matrix
  # r_hat_vec: a vector of estimated ranks
  # lam.est: the selected tuning parameter 
  # norder, lam_path: details in cross validation
  # flag_scaled_mat: a matrix of the same dimension as lam_path, 
  #                whether or not scaled iRRR converges in each cross validation fit
  # flag_iRRR_mat: a matrix of the same dimension as lam_path, 
  #                whether or not all the related iRRR converge in each scaled iRRR fit
  # total_scaled_flag: >0: some scaled iRRR do not converge, 0: all scaled iRRR converge
  # total_iRRR_flag: >0: some iRRR do not converge, 0: all iRRR converge
  # flag_B_mat: a matrix of the same dimension as lam_path,
  #             whether or not the estimated B is zero in each cross validation fit
  # total_B_flag: >0: some fitted B are zero matrices, 0: all fitted B are not zero matrices
  
  Y <- as.matrix(Y)
  n <- nrow(Y)
  p <- ncol(cX)
  q <- ncol(Y)
  K <- length(X) # number of views
  
  if (is.null(norder))
    norder <- sample(seq_len(n),n)
  
  lam_path <- matrix(ncol = nfold, nrow = length(lam.seq), NA)
  # scaled iRRR convergence
  flag_scaled_mat <- matrix(ncol = nfold, nrow = length(lam.seq), NA)
  # iRRR convergence within scaled iRRR
  flag_iRRR_mat <- matrix(ncol = nfold, nrow = length(lam.seq), NA)
  # whether B is zero or not
  flag_B_mat <- matrix(ncol = nfold, nrow = length(lam.seq), NA)
  
  
  ndel <- round(n/nfold)
  for (f in seq_len(nfold)){
    if (f != nfold) {
      iddel <- norder[(1 + ndel * (f - 1)):(ndel * f)]
    }
    else {
      iddel <- norder[(1 + ndel * (f - 1)):n]
    }
    ndel <- length(iddel)
    nf <- n - ndel
    idkeep <- (seq_len(n))[-iddel]
    
    Xf <- cX[-iddel, ]
    Xfdel <- cX[iddel, ]
    Yf <- as.matrix(Y[-iddel, ])
    Yfdel <- as.matrix(Y[iddel, ])
    
    Xf.list <- list()
    tempp1 <- 1+cumsum(c(0,p.vec[-K]))
    tempp2 <- cumsum(p.vec)
    for (i in 1:K){
      Xf.list[[i]] <- Xf[,c(tempp1[i]:tempp2[i])]
    }
    
    for (i in 1:length(lam.seq)){
      fit <- iRRR_normal_scaled(Yf, Xf.list, Xf, p.vec, weight=NULL, lam1 = lam.seq[i], Tol, 
                                Niter, varyrho, rho, lam0,
                                maxrho, randomstart)
      new_B <- fit$B_hat
      new_mu <- fit$mu_hat
      # if new=0: use old version; new=1: use new version 
      lam_path[i,f] <- (1-new)*sum((Yfdel-rep(1,nrow(Yfdel))%*%new_mu-Xfdel%*%new_B)^2) + 
                             new*(sum((Yfdel-rep(1,nrow(Yfdel))%*%new_mu-Xfdel%*%new_B)^2)/fit$sigma_hat^2 + 
                                     ndel*q*log(fit$sigma_hat^2))
      # scaled iRRR convergence
      flag_scaled_mat[i,f] <- fit$conv.flag
      # iRRR convergence within scaled iRRR
      flag_iRRR_mat[i,f] <- fit$iRRR_flag
      # whether B is zero or not
      flag_B_mat[i,f] <- fit$flag_B
    }
  }
  index <- order(colSums(lam_path))
  crerr <- rowSums(lam_path[, index])/length(index) * nfold
  lam.est <- lam.seq[which.min(crerr)]
  
  
  fit <- iRRR_normal_scaled(Y, X, cX, p.vec, weight=NULL, lam1 = lam.est, Tol, 
                            Niter, varyrho, rho, lam0,
                            maxrho, randomstart)
  new_B <- fit$B_hat
  new_mu <- fit$mu_hat
  new_sigma <- fit$sigma_hat
  A <- fit$B_hat_list
  r_hat_vec <- fit$r_hat_vec
  
  total_scaled_flag <- sum(flag_scaled_mat)
  total_iRRR_flag <- sum(flag_iRRR_mat)
  total_B_flag <- sum(flag_B_mat)
  
  return(list(B_hat=new_B, sigma_hat=new_sigma, B_hat_list=A, mu_hat=new_mu,
              lam.est=lam.est, r_hat_vec=r_hat_vec,
              lam_path=lam_path,
              flag_scaled_mat=flag_scaled_mat,
              flag_iRRR_mat=flag_iRRR_mat,
              flag_B_mat=flag_B_mat,
              total_scaled_flag=total_scaled_flag,
              total_iRRR_flag=total_iRRR_flag,
              total_B_flag=total_B_flag,
              norder=norder))
}

#############################################################









#####################   iRRR estimation    ##################
### This function solves the iRRR problem with penalty:   ###
### \sum_k w_k \|X_kB_k\|_* by ADMM algorithm.            ###
### The score matrix estimation relies on this function.  ###


# ADMM
iRRR_normal_XA <- function(Y, X, weight = NULL, lam1, Tol = 1e-3, 
                           Niter = 500, varyrho = 0, rho = 0.1, lam0 = 0,
                           maxrho = 5, randomstart = 0){
  
  # input:
  # X: a list of multi-view data
  # Y: a matrix of response
  # weight: a vector of weights used in iRRR, if NULL then generate it 
  # lam1: the tuning parameter in iRRR
  # Tol: the threshold value of ADMM
  # Niter: the maximum number of iteration in ADMM
  # veryrho: 1 increase rho gradually, 0 do not change rho
  # rho: step size in ADMM
  # lam0: the tuning parameter for ridge penalty
  # maxrho: the largest possible rho if varyrho=1
  # randomstart: whether randomly generate initial values of ADMM
  
  # output:
  # C: the estimated coefficient matrix 
  # mu: the estimated intercept vector
  # B: a list, each list item is an estimated sub-coefficient matrix
  # XA: a list, each list item is the estimated X_kA_k
  # cXA: a matrix of estimated response, i.e., cXA=1%*%mu + \sum_k X_kA_k
  # flag_ADMM_XA: 0 ADMM converges, 1 ADMM does not converge after Niter many iterations
  # obj: a vector of full objective function value (with penalties)
  # obj_ls: a vector of partial objective function value (no penalties)
  # primal: a vector of primal residuals
  # dual: a vector of dual residuals
  
  Y <- as.matrix(Y)
  K <- length(X) # X should be a list
  p <- rep(0, K) # size of each view
  for (i in 1:K){
    p[i] <- ncol(X[[i]])
  }
  n <- nrow(Y)
  q <- ncol(Y)
  
  
  # calculate weights
  if (is.null(weight)){
    weight <- rep(0, K)
    for (i in 1:K){
      # get rank
      rk <- sum(svd(X[[i]])$d > 0.00001)
      weight[i] <- (sqrt(q)+sqrt(rk))/n
    }
  }
  
  
  # horizonatally concatenated X, also column centered and weighted
  cX <- matrix(0, nrow = n, ncol = sum(p))
  meanX <- vector()
  X2 <- list()
  
  for (i in 1:K){
    meanX <- c(meanX, apply(X[[i]], 2, mean)) # collect all means
    # do column centering
    X[[i]] <- scale(X[[i]], center = TRUE, scale = FALSE) 
    X2[[i]] <- t(X[[i]])%*%X[[i]]
    cX[,c((sum(p[1:(i-1)])*(i != 1)+1):sum(p[1:i]))] <- X[[i]] # column centered X
  }
  
  # initial parameter estimates
  mu <- apply(Y, 2, mean, na.rm = TRUE) # a q*1 vector 
  # majorize Y to get a working Y
  wY <- Y
  temp <- rep(1, n)%*%t(mu)
  wY[is.na(wY)] <- temp[is.na(wY)] # wY should be a complete matrix
  mu <- apply(wY, 2, mean) # new est of mu, b/c cX is col centered
  wY1 <- scale(wY, center = TRUE, scale = FALSE) # column centered wY
  
  B <- list()    
  cB <- matrix(0, nrow = sum(p), ncol = q) # vertically concatenated B
  A <- list()
  XA <- list()
  vXA <- matrix(0, nrow = n*K, ncol = q)
  Theta <- list()   # Lagrange params for XB
  cTheta <- matrix(0, nrow = n*K, ncol = q)
  
  for (i in 1:K){
    if (randomstart){
      B[[i]] <- matrix(rnorm(p[i]*q), nrow = p[i]) # initialize the algorithm
    } else {
      B[[i]] <- ginv(t(X[[i]])%*%X[[i]])%*%t(X[[i]])%*%wY1 
    }
    Theta[[i]] <- matrix(0, nrow = n, ncol = q)
    cB[c((sum(p[1:(i-1)])*(i != 1) + 1):sum(p[1:i])),] <- B[[i]]
    A[[i]] <- B[[i]]
    XA[[i]] <- X[[i]]%*%A[[i]]
    vXA[c(((i-1)*n+1):(i*n)),] <- XA[[i]]
  }
  
  
  
  
  if (varyrho == 0){
    # compute inv(1/nX'X+rho*X_d'X_d)
    DeltaMat <- solve( t(cX)%*%cX/n + rho*bdiag(X2) + 0.000001*diag(ncol(cX)))
  }
  
  # check obj value
  # full objective function (with penalties) on observed data
  obj <- ObjValue3(Y, X, mu, XA, B, weight, lam1) 
  # only the least square part on observed data
  obj_ls <- ObjValue3(Y, X, mu, XA, B, weight, 0)
  
  
  #################
  # ADMM
  flag_ADMM_XA <- 0
  niter <- 0
  diff <- 1e+10
  rec_obj_all <- obj  # record obj value
  rec_obj_ls_all <- obj_ls # record obj_ls value
  rec_Theta <- vector() # record the Fro norm of Theta
  rec_primal <- vector() # record total primal residual
  rec_dual <- vector() # record total dual residual
  
  while ((niter < Niter)&(abs(diff) > Tol)){
    niter <- niter + 1
    cB_old <- cB
    B_old <- B
    XA_old <- XA
    
    ########### majorization #########################
    Eta <- rep(1, n)%*%t(mu) + cX%*%cB  # current linear predictor
    wY <- Y
    wY[is.na(wY)] <- Eta[is.na(wY)]
    mu <- apply(wY, 2, mean)
    wY1 <- scale(wY, center = TRUE, scale = FALSE)
    
    # est concatenated B
    if (varyrho) {
      DeltaMat <- solve( t(cX)%*%cX/n + rho*bdiag(X2)+ 0.000001*diag(ncol(cX)))
    }
    
    cB <- DeltaMat%*%(t(cX)%*%wY1/n+rho*t(bdiag(X))%*%(vXA+cTheta))
    for (j in 1:K){
      B[[j]] <- cB[c((sum(p[1:(j-1)])*(j != 1)+1):sum(p[1:j])),]
    }
    
    # est X_iA_i 
    # update Theta_i right after est X_iA_i
    cXA <- matrix(0, nrow = n, ncol = q) 
    for (k in 1:K){
      # est X_iA_i
      temp <- X[[k]]%*%B[[k]] - Theta[[k]]/rho  
      decomp1 <- svd(temp, nu = min(q, n), nv = min(q, n))
      if (length(decomp1$d)>1){
        XA[[k]] <- decomp1$u%*%diag(SoftThres(decomp1$d, lam1*weight[k]/rho))%*%t(decomp1$v)
      }else{
        XA[[k]] <- decomp1$u%*%max(decomp1$d-weight[k]*lam1/rho,0)%*%t(decomp1$v)
      }
      
      
      
      # update Theta
      Theta[[k]] <- Theta[[k]] + rho*(XA[[k]]-X[[k]]%*%B[[k]])
      
      # update cXA and cTheta
      cXA <- cXA + XA[[k]]
      vXA[c(((k-1)*n+1):(k*n)),] <- XA[[k]]
      cTheta[c(((k-1)*n+1):(k*n)),] <- as.matrix(Theta[[k]])
    }
    
    # update rho
    if (varyrho) { rho <- min(maxrho, 1.001*rho) } # steadily increasing rho
    
    # stopping rule
    # primal and dual residuals
    primal <- 0
    dual <- 0
    # to compute relative difference
    rela.A <- 0
    rela.B <- 0
    rela.theta <- 0
    for (i in 1:K){
      primal <- primal + sum((XA[[i]]-X[[i]]%*%B[[i]])^2)
      dual <- dual + sum((rho*t(X[[i]])%*%(XA_old[[i]]-XA[[i]]))^2)
      rela.A <- rela.A + sum(XA[[i]]^2)
      rela.B <- rela.B + sum((X[[i]]%*%B[[i]])^2)
      rela.theta <- rela.theta + sum((t(X[[i]])%*%Theta[[i]])^2)
    }
    primal <- sqrt(primal)
    dual <- sqrt(dual)
    rela.A <- sqrt(rela.A)
    rela.B <- sqrt(rela.B)
    rela.theta <- sqrt(rela.theta)
    
    rec_primal <- c(rec_primal, primal)
    rec_dual <- c(rec_dual, dual)
    rec_Theta <- c(rec_Theta,sum((cTheta)^2))
    
    
    # objective function value
    obj <- ObjValue3(Y,X,mu,XA,B,weight,lam1)
    obj_ls <- ObjValue3(Y,X,mu,XA,B,weight,0)
    rec_obj_all <- c(rec_obj_all, obj)
    rec_obj_ls_all <- c(rec_obj_ls_all, obj_ls)
    
    # stopping rule
    #diff <- max(dual/rela.theta,primal/max(rela.A,rela.B))
    diff <- max(dual,primal)
  }
  
  
  if (niter == Niter){
    flag_ADMM_XA <- 1
    cat("iRRR does NOT converge after_", niter, "_iterations!", "\n", sep = "")
  } else {
    cat("iRRR converges after_", niter, "_iterations.", "\n", sep = "")
  }
  
  
  # output
  # add back mean
  cXA <- cXA+rep(1,n)%*%t(mu)
  mu <- mu-(t(meanX)%*%cB)
  
  return(list(C=cB, mu=mu, B=B, Theta=Theta, primal=rec_primal, dual=rec_dual,
              obj=rec_obj_all, obj_ls=rec_obj_ls_all, XA=XA, cXA=cXA,
              flag_ADMM_XA=flag_ADMM_XA))
}

#############################################################





# tuning by cross validation
iRRR_normal_CV_XA <- function(Y, X, cX, p.vec, weight, lam.seq, nfold, norder = NULL, Tol = 1e-3, 
                           Niter = 500, varyrho = 0, rho = 0.1, lam0 = 0,
                           maxrho = 5, randomstart = 0){
  
  # input:
  # X: a list of multi-view data
  # cX: the integrated design matrix
  # Y: a matrix of response
  # p.vec: a vector of group sizes
  # weight: a vector of weights used in iRRR, if NULL then generate it 
  # lam.seq: the tuning sequence
  # nfold: used in cross validation
  # norder: used in cross validation, if NULL then randomly generate it
  # Tol: the threshold value of ADMM
  # Niter: the maximum number of iteration in ADMM
  # veryrho: 1 increase rho gradually, 0 do not change rho
  # rho: step size in ADMM
  # lam0: the tuning parameter for ridge penalty
  # maxrho: the largest possible rho if varyrho=1
  # randomstart: whether randomly generate initial values of ADMM
  
  # output:
  # B_hat: the estimated coefficient matrix 
  # mu_hat: the estimated intercept vector
  # B_hat_list: a list, each list item is an estimated sub-coefficient matrix
  # r_hat_vec: a vector of estimated ranks
  # lam.est: the selected tuning parameter 
  # norder, lam_path: details in cross validation 
  # flag_iRRR_mat: a matrix of the same dimension as lam_path, 
  #                whether or not iRRR converges in each cross validation fit
  # total_iRRR_flag: >0: some iRRR do not converge, 0: all iRRR converge
  # flag_B_mat: a matrix of the same dimension as lam_path,
  #             whether or not the estimated B is zero in each cross validation fit
  # total_B_flag: >0: some fitted B are zero matrices, 0: all fitted B are not zero matrices
  
  Y <- as.matrix(Y)
  n <- nrow(Y)
  p <- ncol(cX)
  q <- ncol(Y)
  K <- length(X) # number of views
  
  if (is.null(norder))
    norder <- sample(seq_len(n),n)
  
  lam_path <- matrix(ncol = nfold, nrow = length(lam.seq), NA)
  
  ndel <- round(n/nfold)
  for (f in seq_len(nfold)){
    if (f != nfold) {
      iddel <- norder[(1 + ndel * (f - 1)):(ndel * f)]
    }
    else {
      iddel <- norder[(1 + ndel * (f - 1)):n]
    }
    ndel <- length(iddel)
    nf <- n - ndel
    idkeep <- (seq_len(n))[-iddel]
    
    Xf <- cX[-iddel, ]
    Xfdel <- cX[iddel, ]
    Yf <- as.matrix(Y[-iddel, ])
    Yfdel <- as.matrix(Y[iddel, ])
    
    Xf.list <- list()
    weight.f <- rep(0, K)
    tempp1 <- 1+cumsum(c(0,p.vec[-K]))
    tempp2 <- cumsum(p.vec)
    for (i in 1:K){
      Xf.list[[i]] <- Xf[,c(tempp1[i]:tempp2[i])]
      #################################
      #################################
      ####### change K and n #########
      #################################
      #################################
      weight.f[i] <- (sqrt(q*p.vec[i]) + sqrt(2*log(K)))/sqrt(q)/nrow(Xf.list[[i]])
    }
    
    for (i in 1:length(lam.seq)){
      fit <- iRRR_normal_XA(Yf, Xf.list, weight = weight.f, lam1 = lam.seq[i], Tol, 
                           Niter, varyrho, rho, lam0,
                           maxrho, randomstart)
      lam_path[i,f] <- sum((Yfdel- rep(1,nrow(Yfdel))%*%fit$mu - Xfdel%*%fit$C)^2)
    }
  }
  index <- order(colSums(lam_path))
  crerr <- rowSums(lam_path[, index])/length(index) * nfold
  lam.est <- lam.seq[which.min(crerr)]
  
  # use all data to get score matrix with the selected lambda
  fit <- iRRR_normal_XA(Y, X, weight = weight, lam.est, Tol, 
                       Niter, varyrho, rho, lam0,
                       maxrho, randomstart)
  Z <- Y - fit$cXA
  P_z <- Z%*%solve(t(Z)%*%Z)%*%t(Z)
  
  return(list(Z=Z, P_z=P_z))
}

#############################################################






# function used in ADMM
ObjValue3 <- function(Y, X, mu, XA, B, weight, lam1){
  # Calc 1/(2n)|Y-1*mu'-sum(Xi*Bi)|^2 + lam1*sum(|XiBi|_*) 
  # with column centered Xi's and (potentially non-centered and missing) Y
  n <- nrow(Y)
  K <- length(X)
  obj <- 0
  pred <- rep(1, n)%*%t(mu)
  for (i in 1:K){
    pred <- pred + X[[i]]%*%B[[i]]
    obj <- obj + lam1*weight[i]*sum(svd(X[[i]]%*%B[[i]])$d)
  }
  obj <- obj + 1/2/n*sum((Y-pred)^2, na.rm = TRUE)
  return(obj)
}

#############################################################








################## scaled iRRR estimation #################
### This function solves the scaled iRRR problem        ###
### by ADMM algorithm.                                  ###



# ADMM
scalediRRR_ADMM <- function(Y, X, weight = NULL, lam1, Tol = 1e-3, 
                            Niter = 500, varyrho = 0, rho = 0.1, lam0 = 0,
                            maxrho = 5, randomstart = 0){
  
  # input:
  # X: a list of multi-view data
  # Y: a matrix of response
  # weight: a vector of weights used in scaled iRRR, if NULL then generate it 
  # lam1: the tuning parameter in scaled iRRR
  # Tol: the threshold value of ADMM
  # Niter: the maximum number of iteration in ADMM
  # veryrho: 1 increase rho gradually, 0 do not change rho
  # rho: step size in ADMM
  # lam0: the tuning parameter for ridge penalty
  # maxrho: the largest possible rho if varyrho=1
  # randomstart: whether randomly generate initial values of ADMM
  
  # output:
  # C: the estimated coefficient matrix 
  # mu: the estimated intercept vector
  # A: a list, each list item is an estimated sub-coefficient matrix
  # flag_ADMM: 0 ADMM converges, 1 ADMM does not converge after Niter many iterations
  # flag_C: 0 the estimated matrix is not a zero matrix, 1 it is a zero matrix
  # r_hat_vec: a vector of estimated ranks
  # rec_obj_all: a vector of full objective function value (with penalties)
  # rec_obj_ls_all: a vector of partial objective function value (no penalties)
  # rec_primal: a vector of primal residuals
  # rec_dual: a vector of dual residuals
  # rec_diff: a vector of max(primal residual, dual residual)
  
  
  Y <- as.matrix(Y)
  K <- length(X) # X should be a list
  p <- rep(0, K) # size of each view
  n <- nrow(Y)
  q <- ncol(Y)
  
  for (i in 1:K){
    p[i] <- ncol(X[[i]])
  }
  
  # default weights
  if(is.null(weight)){
    weight <- rep(0, K)
    for (i in 1:K){ 
      weight[i] <- max(svd(X[[i]])$d)*(sqrt(q*p.vec[i]) +  sqrt(2*log(K)))/n/q
    }
  }
  
  
  # horizonatally concatenated X, also column centered and weighted
  cX <- matrix(0, nrow = n, ncol = sum(p))
  meanX <- vector()
  
  for (i in 1:K){
    meanX <- c(meanX, apply(X[[i]], 2, mean)) # collect all means
    # first, do column centering
    X[[i]] <- scale(X[[i]], center = TRUE, scale = FALSE) 
    # second, normalize centered X[[i]]'s by wi's
    X[[i]] <- X[[i]]/weight[i] 
    cX[,c((sum(p[1:(i-1)])*(i != 1)+1):sum(p[1:i]))] <- X[[i]] # column centered X
  }
  
  # initial parameter estimates
  mu <- apply(Y, 2, mean, na.rm = TRUE) # a q*1 vector 
  # majorize Y to get a working Y
  wY <- Y
  temp <- rep(1, n)%*%t(mu)
  wY[is.na(wY)] <- temp[is.na(wY)] # wY should be a complete matrix
  mu <- apply(wY, 2, mean) # new est of mu, b/c cX is col centered
  wY1 <- scale(wY, center = TRUE, scale = FALSE) # column centered wY
  
  B <- list()    
  Theta <- list()   # Lagrange params for B
  cB <- matrix(0, nrow = sum(p), ncol = q) # vertically concatenated B
  for (i in 1:K){
    if (randomstart){
      B[[i]] <- matrix(rnorm(p[i]*q), nrow = p[i]) # initialize the algorithm
    } else {
      B[[i]] <- ginv(t(X[[i]])%*%X[[i]])%*%t(X[[i]])%*%wY1 
    }
    Theta[[i]] <- matrix(0, nrow = p[i], ncol = q)
    cB[c((sum(p[1:(i-1)])*(i != 1) + 1):sum(p[1:i])),] <- B[[i]]
  }
  A <- B # low-rank alias
  cA <- cB
  cTheta <- matrix(0, nrow = sum(p), ncol = q) 
  
  
  # initialize noise level estimation
  noise_est <- max(sqrt(2*ObjValue1(Y, X, mu, B, 0, 0)/q), 0.0001)
  
  
  # check obj value
  # full objective function (with penalties) on observed data
  obj <- ObjValue2(Y, X, mu, A, lam1, noise_est) 
  # only the least square part on observed data
  obj_ls <- ObjValue1(Y, X, mu, A, 0, 0)
  
  
  
  #################
  # ADMM
  flag_ADMM <- 0
  niter <- 0
  diff <- 1e+10
  rec_obj_all <- obj  # record obj value
  rec_obj_ls_all <- obj_ls # record obj_ls value
  rec_Theta <- vector() # record the Fro norm of Theta
  rec_primal <- vector() # record primal residual
  rec_dual <- vector() # record dual residual
  rec_diff <- vector() # record diff
  
  while ((niter < Niter)&(abs(diff) > Tol)){
    niter <- niter + 1
    cA_old <- cA
    
    
    ########### majorization #########################
    Eta <- rep(1, n)%*%t(mu) + cX%*%cB  # current linear predictor
    wY <- Y
    wY[is.na(wY)] <- Eta[is.na(wY)]
    mu <- apply(wY, 2, mean)
    wY1 <- scale(wY, center = TRUE, scale = FALSE)
    
    # est concatenated B
    decomp <- svd(cX/sqrt(n*q*noise_est), nv = min(n, sum(p)))
    DeltaMat <- decomp$v%*%diag(1/(decomp$d^2 + lam0 + rho))%*%t(decomp$v) + 
      (diag(rep(1,sum(p)))-decomp$v%*%t(decomp$v))/(lam0 + rho)
    cB <- DeltaMat%*%(t(cX)%*%wY1/n/q/noise_est+rho*cA-cTheta)
    for (j in 1:K){
      B[[j]] <- cB[c((sum(p[1:(j-1)])*(j != 1)+1):sum(p[1:j])),]
    }
    
    # est sigma
    noise_est <- sqrt(2*ObjValue1(Y, X, mu, B, 0, 0)/q)
    
    # est A_i 
    # update Theta_i right after est A_i
    for (k in 1:K){
      # est A
      temp <- B[[k]] + Theta[[k]]/rho
      decomp1 <- svd(temp, nu = min(q, p[k]), nv = min(q, p[k]))
      if (length(decomp1$d)>1){
        A[[k]] <- decomp1$u%*%diag(SoftThres(decomp1$d, lam1/rho))%*%t(decomp1$v)
      }else{
        A[[k]] <- decomp1$u%*%max(decomp1$d-lam1/rho,0)%*%t(decomp1$v)
      }
      
      # update Theta
      Theta[[k]] <- Theta[[k]] + rho*(B[[k]]-A[[k]])
      
      # update cA and cTheta
      cA[c((sum(p[1:(k-1)])*(k != 1)+1):sum(p[1:k])),] <- A[[k]]
      cTheta[c((sum(p[1:(k-1)])*(k != 1)+1):sum(p[1:k])),] <- Theta[[k]]
    }
    
    
    # update rho
    if (varyrho) { rho <- min(maxrho, 1.01*rho) } # steadily increasing rho
    
    
    
    # stopping rule
    # primal and dual residuals
    primal <- sum((cA-cB)^2)
    primal <- sqrt(primal)
    rec_primal <- c(rec_primal, primal)
    
    dual <- sum((cA-cA_old)^2)
    dual <- rho*sqrt(dual)
    rec_dual <- c(rec_dual, dual)
    
    # Fro norm of Theta
    rec_Theta <- c(rec_Theta, sqrt(sum(cTheta^2)))
    
    
    
    # objective function value
    obj <- ObjValue2(Y,X,mu,A,lam1,noise_est)
    obj_ls <- ObjValue1(Y,X,mu,A,0,0)
    rec_obj_all <- c(rec_obj_all, obj)
    rec_obj_ls_all <- c(rec_obj_ls_all, obj_ls)
    
    # stopping rule
    diff <- max(primal, dual)
    rec_diff <- c(rec_diff,diff)
  }
  
  if (niter == Niter){
    flag_ADMM <- 1
    cat("iRRR does NOT converge after_", niter, "_iterations!", "\n", sep = "")
  } else {
    cat("iRRR converges after_", niter, "_iterations.", "\n", sep = "")
  }
  
  
  # output
  # rescale parameter estimate, and add back mean
  C <- matrix(0, nrow = sum(p), ncol = q)
  r_hat_vec <- rep(0,K)
  for (i in 1:K){
    A[[i]] <- A[[i]]/weight[i]
    B[[i]] <- B[[i]]/weight[i]
    C[c((sum(p[1:(i-1)])*(i != 1)+1):sum(p[1:i])),] <- A[[i]]
    r_hat_vec[i] <- sum(svd(A[[i]])$d>1e-2)
  }
  mu <- mu-(t(meanX)%*%C)
  
  flag_C <- 0
  if (sum(C^2)==0) flag_C <- 1
  
  
  return(list(C=C, noise_est=noise_est, mu=mu, A=A, B=B, Theta=Theta, 
              flag_ADMM=flag_ADMM, flag_C=flag_C,
              r_hat_vec=r_hat_vec,
              rec_obj_all=rec_obj_all, rec_obj_ls_all=rec_obj_ls_all,
              rec_primal=rec_primal, rec_dual=rec_dual,
              rec_Theta=rec_Theta, rec_diff=rec_diff))
}

#############################################################






# function used in ADMM
ObjValue2 <- function(Y, X, mu, B, lam1, sigma_est){
  # Calc 1/(2nq\sigma)|Y-1*mu'-sum(Xi*Bi)|^2 + sigma/2 + lam1*sum(|Bi|_*) 
  # with column centered Xi's and (potentially non-centered and missing) Y
  n <- nrow(Y)
  K <- length(X)
  obj <- 0
  pred <- rep(1, n)%*%t(mu)
  for (i in 1:K){
    pred <- pred + X[[i]]%*%B[[i]]
    obj <- obj + lam1*sum(svd(B[[i]])$d)
  }
  obj <- obj + sum((Y-pred)^2, na.rm = TRUE)/2/n/q/sigma_est + sigma_est/2
  return(obj)
}

#############################################################






# tuning by cross validation
iRRR_normal_CV_scaled_ADMM <- function(Y, X, cX, p.vec, weight=NULL, lam.seq, nfold, norder = NULL, Tol = 1e-3, 
                                       Niter = 500, varyrho = 0, rho = 0.1, lam0 = 0,
                                       maxrho = 5, randomstart = 0, new = 0){
  
  # input:
  # X: a list of multi-view data
  # cX: the integrated design matrix
  # Y: a matrix of response
  # p.vec: a vector of group sizes
  # weight: a vector of weights used in scaled iRRR, if NULL then generate it 
  # lam.seq: the tuning sequence
  # nfold: used in cross validation
  # norder: used in cross validation, if NULL then randomly generate it
  # Tol: the threshold value of ADMM
  # Niter: the maximum number of iteration in ADMM
  # veryrho: 1 increase rho gradually, 0 do not change rho
  # rho: step size in ADMM
  # lam0: the tuning parameter for ridge penalty
  # maxrho: the largest possible rho if varyrho=1
  # randomstart: whether randomly generate initial values of ADMM
  # new: whether or not use the negative log-likelihood function as an error measurement in cross validation
  #      0 do not use, 1 use
  
  # output:
  # B_hat: the estimated coefficient matrix
  # sigma_hat: the estimated noise level
  # mu_hat: the estimated intercept vector
  # B_hat_list: a list, each list item is an estimated sub-coefficient matrix
  # r_hat_vec: a vector of estimated ranks
  # lam.est: the selected tuning parameter 
  # norder, lam_path: details in cross validation
  # flag_scaled_mat: a matrix of the same dimension as lam_path, 
  #                whether or not scaled iRRR converges in each cross validation fit
  # total_scaled_flag: >0: some scaled iRRR do not converge, 0: all scaled iRRR converge
  # flag_B_mat: a matrix of the same dimension as lam_path,
  #             whether or not the estimated B is zero in each cross validation fit
  # total_B_flag: >0: some fitted B are zero matrices, 0: all fitted B are not zero matrices
  

  Y <- as.matrix(Y)
  n <- nrow(Y)
  p <- ncol(cX)
  q <- ncol(Y)
  K <- length(X) # number of views
  
  if (is.null(norder))
    norder <- sample(seq_len(n),n)
  
  lam_path <- matrix(ncol = nfold, nrow = length(lam.seq), NA)
  # iRRR convergence within scaled iRRR
  flag_scaled_mat <- matrix(ncol = nfold, nrow = length(lam.seq), NA)
  # whether B is zero or not
  flag_B_mat <- matrix(ncol = nfold, nrow = length(lam.seq), NA)
  
  
  ndel <- round(n/nfold)
  for (f in seq_len(nfold)){
    if (f != nfold) {
      iddel <- norder[(1 + ndel * (f - 1)):(ndel * f)]
    }
    else {
      iddel <- norder[(1 + ndel * (f - 1)):n]
    }
    ndel <- length(iddel)
    nf <- n - ndel
    idkeep <- (seq_len(n))[-iddel]
    
    Xf <- cX[-iddel, ]
    Xfdel <- cX[iddel, ]
    Yf <- as.matrix(Y[-iddel, ])
    Yfdel <- as.matrix(Y[iddel, ])
    
    Xf.list <- list()
    tempp1 <- 1+cumsum(c(0,p.vec[-K]))
    tempp2 <- cumsum(p.vec)
    for (i in 1:K){
      Xf.list[[i]] <- Xf[,c(tempp1[i]:tempp2[i])]
    }
    
    for (i in 1:length(lam.seq)){
      fit <- scalediRRR_ADMM(Yf, Xf.list, weight=NULL, lam1 = lam.seq[i], Tol, 
                             Niter, varyrho, rho, lam0,
                             maxrho, randomstart)
      new_B <- fit$C
      new_mu <- fit$mu
      # if new=0: use old version; new=1: use new version 
      lam_path[i,f] <- (1-new)*sum((Yfdel-rep(1,nrow(Yfdel))%*%new_mu-Xfdel%*%new_B)^2) + 
        new*(sum((Yfdel-rep(1,nrow(Yfdel))%*%new_mu-Xfdel%*%new_B)^2)/fit$noise_est^2 + 
               ndel*q*log(fit$noise_est^2))
      # convergence of scaled iRRR
      flag_scaled_mat[i,f] <- fit$flag_ADMM
      # whether B is zero or not
      flag_B_mat[i,f] <- fit$flag_C
    }
  }
  index <- order(colSums(lam_path))
  crerr <- rowSums(lam_path[, index])/length(index) * nfold
  lam.est <- lam.seq[which.min(crerr)]
  
  fit <- scalediRRR_ADMM(Y, X, weight=NULL, lam1 = lam.est, Tol, 
                         Niter, varyrho, rho, lam0,
                         maxrho, randomstart)
  new_B <- fit$C
  new_mu <- fit$mu
  new_sigma <- fit$noise_est
  A <- fit$A
  r_hat_vec <- fit$r_hat_vec
  
  total_scaled_flag <- sum(flag_scaled_mat)
  total_B_flag <- sum(flag_B_mat)
  
  return(list(B_hat=new_B, sigma_hat=new_sigma, B_hat_list=A, mu_hat=new_mu,
              lam.est=lam.est, r_hat_vec=r_hat_vec,
              lam_path=lam_path,
              flag_scaled_mat=flag_scaled_mat,
              flag_B_mat=flag_B_mat,
              total_scaled_flag=total_scaled_flag,
              total_B_flag=total_B_flag,
              norder=norder))
}

#############################################################





################## scaled iRRR estimation #################
### This function solves the scaled iRRR problem        ###
### by ADMM algorithm.                                  ###
### Control variables are added into the model.         ###


# ADMM
scalediRRR_ADMM_general <- function(Y, X, weight = NULL, lam1, Tol = 1e-3, 
                                    Niter = 500, varyrho = 0, rho = 0.1, lam0 = 0,
                                    maxrho = 5, randomstart = 0){
  
  # input:
  # X: a list of multi-view data
  # Y: a matrix of response
  # weight: a vector of weights used in scaled iRRR, if NULL then generate it 
  # lam1: the tuning parameter in scaled iRRR
  # Tol: the threshold value of ADMM
  # Niter: the maximum number of iteration in ADMM
  # veryrho: 1 increase rho gradually, 0 do not change rho
  # rho: step size in ADMM
  # lam0: the tuning parameter for ridge penalty
  # maxrho: the largest possible rho if varyrho=1
  # randomstart: whether randomly generate initial values of ADMM
  
  # output:
  # C: the estimated coefficient matrix 
  # noise_est: the estimated noise level
  # mu: the estimated intercept vector
  # A: a list, each list item is an estimated sub-coefficient matrix
  # flag_ADMM: 0 ADMM converges, 1 ADMM does not converge after Niter many iterations
  # flag_C: 0 the estimated matrix is not a zero matrix, 1 it is a zero matrix
  # r_hat_vec: a vector of estimated ranks
  # rec_obj_all: a vector of full objective function value (with penalties)
  # rec_obj_ls_all: a vector of partial objective function value (no penalties)
  # rec_primal: a vector of primal residuals
  # rec_dual: a vector of dual residuals
  # rec_diff: a vector of max(primal residual, dual residual)
  
  
  Y <- as.matrix(Y)
  K <- length(X) # X should be a list
  p <- rep(0, K) # size of each view
  n <- nrow(Y)
  q <- ncol(Y)
  
  for (i in 1:K){
    p[i] <- ncol(X[[i]])
  }
  
  # default weights
  if(is.null(weight)){
    weight <- rep(0, K)
    for (i in 1:K){ 
      weight[i] <- max(svd(X[[i]])$d)*(sqrt(q*p.vec[i]) +  sqrt(2*log(K)))/n/q
    }
  }
  
  
  # horizonatally concatenated X, also column centered
  cX <- matrix(0, nrow = n, ncol = sum(p))
  meanX <- vector()
  
  for (i in 1:K){
    meanX <- c(meanX, apply(X[[i]], 2, mean)) # collect all means
    # do column centering
    X[[i]] <- scale(X[[i]], center = TRUE, scale = FALSE) 
    cX[,c((sum(p[1:(i-1)])*(i != 1)+1):sum(p[1:i]))] <- X[[i]] # column centered X
  }
  
  # initial parameter estimates
  mu <- apply(Y, 2, mean, na.rm = TRUE) # a q*1 vector 
  # majorize Y to get a working Y
  wY <- Y
  temp <- rep(1, n)%*%t(mu)
  wY[is.na(wY)] <- temp[is.na(wY)] # wY should be a complete matrix
  mu <- apply(wY, 2, mean) # new est of mu, b/c cX is col centered
  wY1 <- scale(wY, center = TRUE, scale = FALSE) # column centered wY
  
  B <- list()    
  Theta <- list()   # Lagrange params for B
  cB <- matrix(0, nrow = sum(p), ncol = q) # vertically concatenated B
  for (i in 1:K){
    if (randomstart){
      B[[i]] <- matrix(rnorm(p[i]*q), nrow = p[i]) # initialize the algorithm
    } else {
      B[[i]] <- ginv(t(X[[i]])%*%X[[i]])%*%t(X[[i]])%*%wY1 
    }
    Theta[[i]] <- matrix(0, nrow = p[i], ncol = q)
    cB[c((sum(p[1:(i-1)])*(i != 1) + 1):sum(p[1:i])),] <- B[[i]]
  }
  A <- B # low-rank alias
  cA <- cB
  cTheta <- matrix(0, nrow = sum(p), ncol = q) 
  
  
  # initialize noise level estimation
  noise_est <- max(sqrt(2*ObjValue1(Y, X, mu, B, 0, 0)/q), 0.0001)
  
  
  # check obj value
  # full objective function (with penalties) on observed data
  obj <- ObjValue2_general(Y, X, mu, A, lam1, noise_est, weight) 
  # only the least square part on observed data
  obj_ls <- ObjValue1(Y, X, mu, A, 0, 0)
  
  
  
  #################
  # ADMM
  flag_ADMM <- 0
  niter <- 0
  diff <- 1e+10
  rec_obj_all <- obj  # record obj value
  rec_obj_ls_all <- obj_ls # record obj_ls value
  rec_Theta <- vector() # record the Fro norm of Theta
  rec_primal <- vector() # record total primal residual
  rec_dual <- vector() # record total dual residual
  rec_diff <- vector() # record total diff
  
  while ((niter < Niter)&(abs(diff) > Tol)){
    niter <- niter + 1
    cA_old <- cA
    
    
    ########### majorization #########################
    Eta <- rep(1, n)%*%t(mu) + cX%*%cB  # current linear predictor
    wY <- Y
    wY[is.na(wY)] <- Eta[is.na(wY)]
    mu <- apply(wY, 2, mean)
    wY1 <- scale(wY, center = TRUE, scale = FALSE)
    
    # est concatenated B
    decomp <- svd(cX/sqrt(n*q*noise_est), nv = min(n, sum(p)))
    DeltaMat <- decomp$v%*%diag(1/(decomp$d^2 + lam0 + rho))%*%t(decomp$v) + 
      (diag(rep(1,sum(p)))-decomp$v%*%t(decomp$v))/(lam0 + rho)
    cB <- DeltaMat%*%(t(cX)%*%wY1/n/q/noise_est+rho*cA-cTheta)
    for (j in 1:K){
      B[[j]] <- cB[c((sum(p[1:(j-1)])*(j != 1)+1):sum(p[1:j])),]
    }
    
    # est sigma
    noise_est <- sqrt(2*ObjValue1(Y, X, mu, B, 0, 0)/q)
    
    # est A_i 
    # update Theta_i right after est A_i
    for (k in 1:K){
      # est A
      temp <- B[[k]] + Theta[[k]]/rho
      decomp1 <- svd(temp, nu = min(q, p[k]), nv = min(q, p[k]))
      if (length(decomp1$d)>1){
        A[[k]] <- decomp1$u%*%diag(SoftThres(decomp1$d, weight[k]*lam1/rho))%*%t(decomp1$v)
      }else{
        A[[k]] <- decomp1$u%*%max(decomp1$d-weight[k]*lam1/rho,0)%*%t(decomp1$v)
      }
      
      # update Theta
      Theta[[k]] <- Theta[[k]] + rho*(B[[k]]-A[[k]])
      
      # update cA and cTheta
      cA[c((sum(p[1:(k-1)])*(k != 1)+1):sum(p[1:k])),] <- A[[k]]
      cTheta[c((sum(p[1:(k-1)])*(k != 1)+1):sum(p[1:k])),] <- Theta[[k]]
    }
    
    
    # update rho
    if (varyrho) { rho <- min(maxrho, 1.01*rho) } # steadily increasing rho
    
    
    
    # stopping rule
    # primal and dual residuals
    primal <- sum((cA-cB)^2)
    primal <- sqrt(primal)
    rec_primal <- c(rec_primal, primal)
    
    dual <- sum((cA-cA_old)^2)
    dual <- rho*sqrt(dual)
    rec_dual <- c(rec_dual, dual)
    
    # Fro norm of Theta
    rec_Theta <- c(rec_Theta, sqrt(sum(cTheta^2)))
    
    
    
    # objective function value
    obj <- ObjValue2_general(Y,X,mu,A,lam1,noise_est,weight)
    obj_ls <- ObjValue1(Y,X,mu,A,0,0)
    rec_obj_all <- c(rec_obj_all, obj)
    rec_obj_ls_all <- c(rec_obj_ls_all, obj_ls)
    
    # stopping rule
    # diff <- max(primal/max(sqrt(sum(cB^2)),sqrt(sum(cA^2))), dual/sqrt(sum(cTheta^2)) )
    diff <- max(primal, dual)
    rec_diff <- c(rec_diff,diff)
  }
  
  if (niter == Niter){
    flag_ADMM <- 1
    cat("iRRR does NOT converge after_", niter, "_iterations!", "\n", sep = "")
  } else {
    cat("iRRR converges after_", niter, "_iterations.", "\n", sep = "")
  }
  
  
  # output
  C <- matrix(0, nrow = sum(p), ncol = q)
  r_hat_vec <- rep(0,K)
  for (i in 1:K){
    C[c((sum(p[1:(i-1)])*(i != 1)+1):sum(p[1:i])),] <- A[[i]]
    r_hat_vec[i] <- sum(svd(A[[i]])$d>1e-2)
  }
  mu <- mu-(t(meanX)%*%C)
  
  flag_C <- 0
  if (sum(C^2)==0) flag_C <- 1
  
  
  return(list(C=C, noise_est=noise_est, mu=mu, A=A, B=B, Theta=Theta, 
              flag_ADMM=flag_ADMM, flag_C=flag_C,
              r_hat_vec=r_hat_vec,
              rec_obj_all=rec_obj_all, rec_obj_ls_all=rec_obj_ls_all,
              rec_primal=rec_primal, rec_dual=rec_dual,
              rec_Theta=rec_Theta, rec_diff=rec_diff))
}

#############################################################








# function used in ADMM
ObjValue2_general <- function(Y, X, mu, B, lam1, sigma_est, weight){
  # Calc 1/(2nq\sigma)|Y-1*mu'-sum(Xi*Bi)|^2 + sigma/2 + lam1 \sum_{i=1}^K weight_i*sum(|Bi|_*) 
  # with column centered Xi's and (potentially non-centered and missing) Y
  n <- nrow(Y)
  K <- length(X)
  obj <- 0
  pred <- rep(1, n)%*%t(mu)
  for (i in 1:K){
    pred <- pred + X[[i]]%*%B[[i]]
    obj <- obj + lam1*sum(svd(B[[i]])$d)*weight[i]
  }
  obj <- obj + sum((Y-pred)^2, na.rm = TRUE)/2/n/q/sigma_est + sigma_est/2
  return(obj)
}

#############################################################







# tuning by cross validation
iRRR_normal_CV_scaled_ADMM_general <- function(Y, X, cX, p.vec, weight=NULL, lam.seq, nfold, 
                                               norder = NULL, Tol = 1e-3, 
                                               Niter = 500, varyrho = 0, rho = 0.1, lam0 = 0,
                                               maxrho = 5, randomstart = 0, new = 0){
  
  # input:
  # X: a list of multi-view data, let the last view contain the control variables
  # cX: the integrated design matrix
  # Y: a matrix of response
  # p.vec: a vector of group sizes
  # weight: a vector of weights used in scaled iRRR, if NULL then generate it 
  # lam.seq: the tuning sequence
  # nfold: used in cross validation
  # norder: used in cross validation, if NULL then randomly generate it
  # Tol: the threshold value of ADMM
  # Niter: the maximum number of iteration in ADMM
  # veryrho: 1 increase rho gradually, 0 do not change rho
  # rho: step size in ADMM
  # lam0: the tuning parameter for ridge penalty
  # maxrho: the largest possible rho if varyrho=1
  # randomstart: whether randomly generate initial values of ADMM
  # new: whether or not use the negative log-likelihood function as an error measurement in cross validation
  #      0 do not use, 1 use
  
  # output:
  # B_hat: the estimated coefficient matrix
  # sigma_hat: the estimated noise level
  # mu_hat: the estimated intercept vector
  # B_hat_list: a list, each list item is an estimated sub-coefficient matrix
  # r_hat_vec: a vector of estimated ranks
  # lam.est: the selected tuning parameter 
  # norder, lam_path: details in cross validation
  # flag_scaled_mat: a matrix of the same dimension as lam_path, 
  #                whether or not scaled iRRR converges in each cross validation fit
  # total_scaled_flag: >0: some scaled iRRR do not converge, 0: all scaled iRRR converge
  # flag_B_mat: a matrix of the same dimension as lam_path,
  #             whether or not the estimated B is zero in each cross validation fit
  # total_B_flag: >0: some fitted B are zero matrices, 0: all fitted B are not zero matrices
  
  
  
  Y <- as.matrix(Y)
  n <- nrow(Y)
  p <- ncol(cX)
  q <- ncol(Y)
  K <- length(X) # number of views
  
  if (is.null(norder))
    norder <- sample(seq_len(n),n)
  
  lam_path <- matrix(ncol = nfold, nrow = length(lam.seq), NA)
  # iRRR convergence within scaled iRRR
  flag_scaled_mat <- matrix(ncol = nfold, nrow = length(lam.seq), NA)
  # whether B is zero or not
  flag_B_mat <- matrix(ncol = nfold, nrow = length(lam.seq), NA)
  
  
  ndel <- round(n/nfold)
  for (f in seq_len(nfold)){
    if (f != nfold) {
      iddel <- norder[(1 + ndel * (f - 1)):(ndel * f)]
    }
    else {
      iddel <- norder[(1 + ndel * (f - 1)):n]
    }
    ndel <- length(iddel)
    nf <- n - ndel
    idkeep <- (seq_len(n))[-iddel]
    
    Xf <- cX[-iddel, ]
    Xfdel <- cX[iddel, ]
    Yf <- as.matrix(Y[-iddel, ])
    Yfdel <- as.matrix(Y[iddel, ])
    
    
    Xf.list <- list()
    tempp1 <- 1+cumsum(c(0,p.vec[-K]))
    tempp2 <- cumsum(p.vec)
    for (i in 1:K){
      Xf.list[[i]] <- as.matrix(Xf[,c(tempp1[i]:tempp2[i])])
    }
    
    if (is.null(weight)){
      weight <- rep(0, K)
      for (i in 1:(K-1)){ 
        weight[i] <- max(svd(Xf.list[[i]])$d)*(sqrt(q*p.vec[i]) +  sqrt(2*log(K)))/nrow(Xf.list[[i]])/q
      }
      weight[K] <- 0  # no penalty on controling variables
    }
    
    
    for (i in 1:length(lam.seq)){
      fit <- scalediRRR_ADMM_general(Yf, Xf.list, weight=weight, lam1 = lam.seq[i], Tol, 
                                     Niter, varyrho, rho, lam0,
                                     maxrho, randomstart)
      new_B <- fit$C
      new_mu <- fit$mu
      # if new=0: use old version; new=1: use new version 
      lam_path[i,f] <- (1-new)*sum((Yfdel-rep(1,nrow(Yfdel))%*%new_mu-Xfdel%*%new_B)^2) + 
        new*(sum((Yfdel-rep(1,nrow(Yfdel))%*%new_mu-Xfdel%*%new_B)^2)/fit$noise_est^2 + 
               ndel*q*log(fit$noise_est^2))
      # convergence of scaled iRRR
      flag_scaled_mat[i,f] <- fit$flag_ADMM
      # whether B is zero or not
      flag_B_mat[i,f] <- fit$flag_C
    }
  }
  index <- order(colSums(lam_path))
  crerr <- rowSums(lam_path[, index])/length(index) * nfold
  lam.est <- lam.seq[which.min(crerr)]
  
  weight.all <- rep(0, K)
  for (i in 1:(K-1)){ 
    weight.all[i] <- max(svd(X[[i]])$d)*(sqrt(q*p.vec[i]) +  sqrt(2*log(K)))/n/q
  }
  weight.all[K] <- 0  # no penalty on controling variables
  
  
  fit <- scalediRRR_ADMM_general(Y, X, weight=weight.all, lam1 = lam.est, Tol, 
                                 Niter, varyrho, rho, lam0,
                                 maxrho, randomstart)
  new_B <- fit$C
  new_mu <- fit$mu
  new_sigma <- fit$noise_est
  A <- fit$A
  r_hat_vec <- fit$r_hat_vec
  
  total_scaled_flag <- sum(flag_scaled_mat)
  total_B_flag <- sum(flag_B_mat)
  
  return(list(B_hat=new_B, sigma_hat=new_sigma, B_hat_list=A, mu_hat=new_mu,
              lam.est=lam.est, r_hat_vec=r_hat_vec,
              lam_path=lam_path,
              flag_scaled_mat=flag_scaled_mat,
              flag_B_mat=flag_B_mat,
              total_scaled_flag=total_scaled_flag,
              total_B_flag=total_B_flag,
              norder=norder))
}

#############################################################


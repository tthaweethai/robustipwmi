# Robust inference when combining inverse-probability weighting and multiple imputation 
# to address missing data with application to an electronic health records-based study of bariatric surgery

# Code for implementing robust variance estimator

# Tanayott (Tony) Thaweethai
# February 27, 2020

library(numDeriv)
library(rootSolve)

  get.ipwmi.se <- function(dset,
                         M,
                         analysis.model.formula,
                         ipw.model.formula,
                         imp.model.formula,
                         S.theta.i,
                         S.alpha.i,
                         imp.method,
                         pi.i){
    
    # Split data into 
    # Rule (User-specified rule; indicator of inclusion. IPW model outcome variable)
    # H (IPW model predictors)
    # Dm  (Variable to be imputed)
    # Dod (Imputation model predictors: D_i^o and D_i^\dagger)
    # Complete (Indicator of D being completely observed)
    # Y (Analysis model outcome)
    # X (Analysis model predictors)
    
    N <- nrow(dset)
    Dm.var <- all.vars(imp.model.formula)[1] 
    Rule.var <- all.vars(ipw.model.formula)[1]

    # Store data in matrices 
    Rule.N <- matrix(dset[,Rule.var],ncol=1)
    H.N <- model.matrix.lm(ipw.model.formula,data=dset, na.action='na.pass')
    Dm.N <- matrix(dset[,Dm.var], ncol=1)
    Dod.N <- model.matrix.lm(imp.model.formula,data=dset, na.action='na.pass')  
    Complete.N <- matrix(!is.na(Dm.N), ncol=1)
    Y.N <- matrix(model.response(model.frame(analysis.model.formula,data=dset,na.action=NULL)), ncol=1)
    X.N <- model.matrix.lm(analysis.model.formula, data=dset, na.action='na.pass')
    
    if (sum(is.na(H.N)) > 0) {stop("IPW model predictors are not fully observed")}
    if (sum(is.na(Dod.N[which(Rule.N==1),])) > 0) {stop("Imputation model predictors are not fully observed among individuals with Rule=1")}
    if (sum(is.na(Y.N[which(Complete.N==1),])) > 0) {stop("Analysis model outcome is not fully observed among individuals among whom the imputation model is fit, 
                                                    and it should be, because those people should have complete data")}
    if (sum(is.na(X.N[which(Complete.N==1),])) > 0) {stop("Analysis model predictors are not fully observed among individuals among whom the imputation model is fit, 
                                                    and it should be, because those people should have complete data")}
    
    # Extract the data for the ith individual 
    get.Di <- function(i){
      list(X = matrix(X.N[i, ], ncol = 1),
           Y = Y.N[i],
           H = matrix(H.N[i, ], ncol = 1),
           Dod = matrix(Dod.N[i, ], ncol = 1),
           Rule = Rule.N[i],
           Complete = Complete.N[i],
           Dm = Dm.N[i])
    }
    
    ### FIRST ANALYTICAL APPROACH: Complete case analysis
    
    # Solve for theta, which indexes the analysis model
    S.theta.cc.i <- function(Di, theta.hat){
      if (Di$Complete == 1){
        return(S.theta.i(Di, theta.hat))
      } else {
        return(matrix(data = 0, nrow = nrow(Di$X), ncol = 1))
      }
    }
    
    n <- sum(Complete.N) # number of complete cases
    
    S.theta.cc.multiroot <- function(theta.hat){
      rowMeans(do.call("cbind",lapply(1:N,function(i){
          Di <- get.Di(i)
          return(S.theta.cc.i(Di, theta.hat))
      })))
    }
    
    theta.cc.hat <- matrix(data = multiroot(S.theta.cc.multiroot,rep(0,nrow(get.Di(1)$X)))$root, ncol=1)
    
    S.theta.cc <-do.call("cbind",lapply(1:N,function(i){
        S.theta.cc.i(get.Di(i), theta.cc.hat)
      }))
    
    # Use a sandwich estimator to estimate standard errors
    
    tau.cc.hat <- 1/n * Reduce("+",lapply(1:N,function(i){
        S.theta.cc <- function(theta.cc.hat){S.theta.cc.i(get.Di(i),theta.cc.hat)}
        return(-jacobian(S.theta.cc,theta.cc.hat))
    }))
    
    se.cc <-sqrt(diag(1/n * solve(tau.cc.hat) %*% S.theta.cc %*% (t(S.theta.cc)/n) %*% t(solve(tau.cc.hat))))
    
    ### SECOND ANALYTICAL APPROACH: IPW/MI
    
    # Estimate alpha hat (IPW model)
    S.alpha.multiroot <- function(alpha.hat){
      rowMeans(do.call("cbind",lapply(1:N,function(i){
        Di <- get.Di(i)
        return(S.alpha.i(Di, alpha.hat))
      })))
    }
    
    alpha.hat <- matrix(data = multiroot(S.alpha.multiroot,rep(0,nrow(get.Di(1)$H)))$root, ncol=1)
    
    # Partition indiviuals who meet the initial inclusion rule into those whose Dm will be imputed and those who will not
    d.to.impute <- dset[which(Rule.N==1 & Complete.N==0),]
    d.not.imputed <- dset[which(!(Rule.N==1 & Complete.N==0)),]
    
    if (imp.method == "linear"){
      # Linear imputation model
      
      # Estimating equation for beta and sigma (collectively psi) which index the imputation model
      S.psi.i <- function(Di, psi.hat){
        beta.hat <- psi.hat[1:(length(psi.hat)-1)]
        sigma.hat <- psi.hat[length(psi.hat)]
        return(rbind(1/sigma.hat^2 * Di$Dod %*% (Di$Dm - t(beta.hat) %*% Di$Dod),
                     -1/sigma.hat+((Di$Dm - t(beta.hat) %*% Di$Dod)^2/sigma.hat^3)))
      }
      
      # Solve for psi hat (beta hat and sigma hat)
      imp.fit <- glm(imp.model.formula, data = dset[which(Rule.N==1),])
      beta.hat <- matrix(coef(imp.fit), ncol = 1)
      sigma.hat <- sigma(imp.fit)
      psi.hat <- matrix(c(beta.hat, sigma.hat), ncol=1)
      
      # Generate imputations
      imputed <- do.call("rbind",lapply(1:M,function(j){
        d.imputed <- d.to.impute
        d.imputed[, Dm.var] <- as.numeric(model.matrix.lm(imp.model.formula, data=d.imputed,na.action='na.pass') %*% beta.hat) + rnorm(nrow(d.imputed),0,sigma.hat)
        stacked <- rbind(d.imputed, d.not.imputed)
        stacked$j <- j
        return(stacked[order(stacked$i,stacked$j),])
      }))
    } else if (imp.method == "logistic") {
      # Logistic imputation model
      
      # Estimating equation for psi (imputation model)
      S.psi.i <- function(Di, psi.hat){
        return(Di$Dod %*% (Di$Dm - plogis(as.numeric(t(psi.hat) %*% Di$Dod))))
      }
      
      # Solve for psi.hat
      imp.fit <- glm(imp.model.formula, data = dset[which(Rule.N==1),], family="binomial")
      psi.hat <- matrix(coef(imp.fit), ncol=1)
      
      # Generate imputations
      imputed <- do.call("rbind",lapply(1:M,function(j){
        d.imputed <- d.to.impute
        d.imputed[, Dm.var] <- rbinom(nrow(d.imputed), 1, predict(imp.fit, d.imputed, type = "response"))
        stacked <- rbind(d.imputed, d.not.imputed)
        stacked$j <- j
        return(stacked[order(stacked$i,stacked$j),])
      }))
    } else {
      print("Must specify imp.method as linear or logistic")
      stop()
    }
    
    # Store the data of the imputed datasets
    
    Rule.NM <- matrix(imputed[,Rule.var],ncol=1)
    H.NM <- model.matrix.lm(ipw.model.formula,data=imputed, na.action='na.pass')
    Dm.NM <- matrix(imputed[,Dm.var], ncol=1)
    Dod.NM <- model.matrix.lm(imp.model.formula,data=imputed, na.action='na.pass')  
    Complete.NM <- matrix(rep(Complete.N,M),ncol=1)
    Y.NM <- matrix(model.response(model.frame(analysis.model.formula, data=imputed,na.action=NULL)), ncol=1)
    X.NM <- model.matrix.lm(analysis.model.formula, data=imputed, na.action='na.pass')
    
    get.Dij <- function(i,j){
      ij <- (j-1)*N + i
      list(X = matrix(X.NM[ij,], ncol=1),
           Y = Y.NM[ij],
           H = matrix(H.NM[ij,],ncol=1),
           Dod = matrix(Dod.NM[ij,],ncol=1),
           Rule = Rule.NM[ij],
           Complete = Complete.NM[ij],
           Dm = Dm.NM[ij])
    }
    
    ### INFERENCE FOR IPW/MI
    
    # IPW/MI estimating equation for theta
    S.theta.ipwmi.i <- function(Di, theta.hat, alpha.hat){
      Wi <- Di$Rule / as.numeric(pi.i(Di, alpha.hat))
      if (Wi == 0){
        return(matrix(data = 0, nrow = nrow(Di$X), ncol = 1))
      } else {
        return(Wi*S.theta.i(Di, theta.hat))
      }
    }
    
    # Estimate theta hat, using IPW/MI
    S.theta.multiroot <- function(theta.hat){
      rowMeans(do.call("cbind",lapply(1:N,function(i){
        do.call("cbind",lapply(1:M,function(j){
          Dij <- get.Dij(i,j)
          return(S.theta.ipwmi.i(Dij, theta.hat, alpha.hat))
        }))
      })))
    }
    
    theta.hat <- matrix(data = multiroot(S.theta.multiroot,rep(0,nrow(get.Di(1)$X)))$root, ncol=1)
    
    # Perform robust inference
    S.theta <- do.call("cbind",lapply(1:M,function(j){
      do.call("cbind",lapply(1:N,function(i){
        S.theta.ipwmi.i(get.Dij(i,j), theta.hat, alpha.hat)
      }))
    }))
    
    S.psi.obs.i <- function(Di, psi.hat){
      if (Di$Complete == 0){
        return(matrix(data = 0, nrow = length(psi.hat), ncol = 1))
      } else {
        return(S.psi.i(Di, psi.hat))
      }
    }
    
    S.psi.mis.i <- function(Di, psi.hat){
      if (Di$Rule == 1 & Di$Complete == 0){ 
        return(S.psi.i(Di, psi.hat))
      } else {
        return(matrix(data = 0, nrow = length(psi.hat), ncol = 1))
      }
    }
    
    S.psi.obs <- do.call("cbind",lapply(1:M,function(j){
      do.call("cbind",lapply(1:N,function(i){
        S.psi.obs.i(get.Dij(i,j), psi.hat)
      }))
    }))
    
    S.psi.mis <- do.call("cbind",lapply(1:M,function(j){
      do.call("cbind",lapply(1:N,function(i){
        S.psi.mis.i(get.Dij(i,j), psi.hat)
      }))
    }))
  
    S.alpha <- do.call("cbind",lapply(1:M,function(j){
      do.call("cbind",lapply(1:N,function(i){
        return(S.alpha.i(get.Dij(i,j),alpha.hat))
      }))
    }))
    
    I.alpha.hat <- 1/N * Reduce("+",lapply(1:N,function(i){
      S.alpha <- function(alpha.hat){S.alpha.i(get.Di(i),alpha.hat)}
      return(-jacobian(S.alpha,alpha.hat))
    }))
    
    I.psi.hat <- 1/N * Reduce("+",lapply(1:N,function(i){
      S.psi <- function(psi.hat){S.psi.obs.i(get.Di(i),psi.hat)}
      return(-jacobian(S.psi,psi.hat))
    }))
    
    kappa.hat <- -1/(N*M) * S.theta %*% t(S.psi.mis)
    
    delta.hat <- 1/(N*M) * Reduce("+",lapply(1:N,function(i){
      Reduce("+", lapply(1:M,function(j){
        S.theta.ipwmi <- function(alpha.hat){S.theta.ipwmi.i(get.Dij(i,j),theta.hat,alpha.hat)}
        return(-jacobian(S.theta.ipwmi,alpha.hat))
      }))
    }))
    
    tau.hat <- 1/(N*M) * Reduce("+",lapply(1:N,function(i){
      Reduce("+", lapply(1:M, function(j){
        S.theta.ipwmi <- function(theta.hat){S.theta.ipwmi.i(get.Dij(i,j),theta.hat,alpha.hat)}
        return(-jacobian(S.theta.ipwmi,theta.hat))
      }))
    }))
    
    I.alpha.hat.inv <- solve(I.alpha.hat)
    I.psi.hat.inv <- solve(I.psi.hat)
    
    S.theta.bar.N<-do.call("cbind",lapply(1:N,function(i){
      apply(S.theta[,i+(seq(0,(M-1)*N,N))],1,mean)
    }))
    
    S.psi.obs.N<-do.call("cbind",lapply(1:N,function(i){
      apply(S.psi.obs[,i+(seq(0,(M-1)*N,N))],1,mean)
    }))
    
    S.alpha.N<-S.alpha[,1:N]
    
    vec<-do.call("cbind",lapply(1:N,function(i){
      return(S.theta.bar.N[, i] - kappa.hat %*% I.psi.hat.inv %*% S.psi.obs.N[, i] - delta.hat %*% I.alpha.hat.inv %*% S.alpha.N[, i])
    }))
    
    se.prop <-sqrt(diag(1/N * solve(tau.hat) %*% vec %*% (t(vec)/N) %*% t(solve(tau.hat))))

    ### PERFORM INFERENCE USING RUBIN'S RULES
    
    # Point estimates for each of the M imputed datasets
    theta.hats <- lapply(1:M, function(j){
      S.theta.j.multiroot <- function(theta.hat){
        rowMeans(do.call("cbind",lapply(1:N,function(i){
            Dij <- get.Dij(i,j)
            return(S.theta.ipwmi.i(Dij, theta.hat, alpha.hat))
        })))
      }
      return(matrix(data = multiroot(S.theta.j.multiroot,rep(0,nrow(get.Di(1)$X)))$root, ncol=1))
    })
    
    # Estimates of the variance of theta hat for each of the M imputed datasets, 
    # using sandwich estimators incorporating estimation of the weights
    V.hats <- lapply(1:M,function(j){
      
      S.theta.j <- S.theta[,1:N+(j-1)*N]
      S.alpha.j <- S.alpha[,1:N+(j-1)*N]
      
      delta.hat.j <- 1/N * Reduce("+",lapply(1:N,function(i){
        S.theta.ipwmi <- function(alpha.hat){S.theta.ipwmi.i(get.Dij(i,j),theta.hat,alpha.hat)}
          return(-jacobian(S.theta.ipwmi,alpha.hat))
      }))
      
      tau.hat.j <- 1/N * Reduce("+",lapply(1:N,function(i){
        S.theta.ipwmi <- function(theta.hat){S.theta.ipwmi.i(get.Dij(i,j),theta.hat,alpha.hat)}
          return(-jacobian(S.theta.ipwmi,theta.hat))
      }))
      
      vec.j <- do.call("cbind",lapply(1:N,function(i){
        return(S.theta.j[, i] - delta.hat.j %*% I.alpha.hat.inv %*% S.alpha.j[, i])
      }))
      
      return(1/N * solve(tau.hat.j) %*% vec.j %*% (t(vec.j)/N) %*% t(solve(tau.hat.j)))
    })
    
    # Apply Rubin's rules
    se.rr <- sqrt(diag(Reduce("+", V.hats)/M + (1+1/M)*1/(M-1)*Reduce("+",lapply(theta.hats,function(x){
      return((x-theta.hat)%*%t(x-theta.hat))
    }))))
    
    result <- data.frame(theta.hat = theta.hat, 
                         se.prop = se.prop,
                         se.rr = se.rr,
                         theta.cc.hat = theta.cc.hat,
                         se.cc = se.cc)
    rownames(result) <- colnames(X.NM)
    return(result)
}



# Vignette
  
  N <- 10^3
  set.seed(1000)
  d <- data.frame(matrix(nrow=N,ncol=0))
  d$i <- 1:N
  d$X1 <- rbinom(N,1,0.5)
  d$X2 <- rnorm(N)
  d$X3 <- rnorm(N)
  d$X4 <- rnorm(N)
  d$X5 <- rnorm(N,d$X2*d$X3)
  d$Y1 <- rnorm(N,-3 + d$X1*d$X2 + d$X1*d$X3 + 0.5*d$X2*d$X3 + d$X4 + 0.5*d$X5,ifelse(d$X1==1,2,1))
  
  d$R1 <- rbinom(N,1,0.8-0.6*d$X1)
  d$R2 <- ifelse(d$R1==1,rbinom(N,1,1/(1+exp(-1.5+0.6*d$X2*d$X4))),rep(0,N))
  d.obs <- d
  d.obs$X2 <- ifelse(d$R1, d$X2, NA)
  d.obs$X3 <- ifelse(d$R1, d$X3, NA)
  d.obs$X4 <- ifelse(d$R1, d$X4, NA)
  d.obs$X5 <- ifelse(d$R1, d$X5, NA)
  d.obs$Y1 <- ifelse(d$R2, d$Y1, NA)
  
  get.ipwmi.se(dset=d.obs,
               M = 10,
               analysis.model.formula = Y1 ~ X2*X3,
               ipw.model.formula = R1 ~ X1,
               imp.model.formula = Y1 ~ X1*X2*X3 + X4 + X5,
               S.theta.i = function(Di,theta.hat){Di$X %*% (Di$Y - t(theta.hat) %*% Di$X)},
               S.alpha.i = function(Di,alpha.hat){Di$H %*% (Di$Rule - t(alpha.hat) %*% Di$H)},
               pi.i = function(Di,alpha.hat){t(alpha.hat) %*% Di$H},
               imp.method = "linear")
  

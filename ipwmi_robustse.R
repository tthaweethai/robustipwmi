# Combine IPW with MI for missing data
# Robust variance estimator

# Tanayott (Tony) Thaweethai
# Monday January 6, 2020

rm(list=ls())
library(numDeriv)
library(rootSolve)

N <- 10^3
M <- 10
bias <- function(true,est){(est-true)/true}

set.seed(1000)
het <- 1
d<-data.frame(matrix(nrow=N,ncol=0))
d$i <- 1:N
d$X1<-rbinom(N,1,0.5)
d$X2<-rnorm(N)
d$X3<-rnorm(N)
d$X4<-rnorm(N)
d$X5<-rnorm(N,d$X2*d$X3)
d$Y<-rnorm(N,-3+d$X1*d$X2+d$X1*d$X3+0.5*d$X2*d$X3+d$X4+0.5*d$X5,ifelse(d$X1==1,ifelse(het==1,2,1),1))
d$R1<-rbinom(N,1,0.8-0.6*d$X1)
d$R2<-ifelse(d$R1==1,rbinom(N,1,1/(1+exp(-1.5+0.6*d$X2*d$X4))),rep(0,N))

d.obs <- d
d.obs$X2 <- ifelse(d$R1, d$X2, NA)
d.obs$X3 <- ifelse(d$R1, d$X3, NA)
d.obs$X4 <- ifelse(d$R1, d$X4, NA)
d.obs$X5 <- ifelse(d$R1, d$X5, NA)
d.obs$Y <- ifelse(d$R2, d$Y, NA)

  get.ipwmi.se <- function(dset,
                         analysis.model.formula,
                         ipw.model.formula,
                         imp.model.formula,
                         S.theta.i,
                         S.alpha.i,
                         imp.method,
                         W.i){
    
    # Split data into 
    # R1 (IPW model outcome), 
    # Z1 (IPW model predictors),
    # A  (Variable to be imputed)
    # Z2 (Imputation model predictors)
    # Y (Analysis model outcome)
    # X (Analysis model predictors)
    
    N <- nrow(dset)
    A.var <- all.vars(imp.model.formula)[1]
    R1.var <- all.vars(ipw.model.formula)[1]
    
    X.N <- model.matrix.lm(analysis.model.formula, data=dset, na.action='na.pass')
    Y.N <- matrix(model.response(model.frame(analysis.model.formula,data=dset,na.action=NULL)), ncol=1)
    Z1.N <- model.matrix.lm(ipw.model.formula,data=dset, na.action='na.pass')
    R1.N <- matrix(dset[,R1.var],ncol=1)
    Z2.N <- model.matrix.lm(imp.model.formula,data=dset, na.action='na.pass')  
    A.N <- matrix(dset[,A.var], ncol=1)
    R2.N <- matrix(!is.na(A.N), ncol=1)
    
    get.Di <- function(i){
      list(X = matrix(X.N[i,],nrow=1),
           Y = Y.N[i],
           Z1 = matrix(Z1.N[i,],nrow=1),
           Z2 = matrix(Z2.N[i,],nrow=1),
           R1 = R1.N[i],
           R2 = R2.N[i],
           A = A.N[i])
    }
    
    # Complete case analysis
    S.theta.cc.i <- function(Di, theta.hat){
      if (Di$R2 == 1){
        return(S.theta.i(Di, theta.hat))
      } else {
        return(matrix(data = 0, nrow = ncol(Di$X), ncol = 1))
      }
    }
    
    # Estimate theta hat, using IPW/MI
    n <- sum(dset$R2)
    S.theta.cc.multiroot <- function(theta.hat){
      rowMeans(do.call("cbind",lapply(1:N,function(i){
          Di <- get.Di(i)
          return(S.theta.cc.i(Di, theta.hat))
      })))
    }
    
    theta.cc.hat <- matrix(data = multiroot(S.theta.cc.multiroot,rep(0,ncol(get.Di(1)$X)))$root, ncol=1)
    
    S.theta.cc <-do.call("cbind",lapply(1:N,function(i){
        S.theta.cc.i(get.Di(i), theta.cc.hat)
      }))
    
    tau.cc.hat <- 1/n * Reduce("+",lapply(1:N,function(i){
        S.theta.cc <- function(theta.cc.hat){S.theta.cc.i(get.Di(i),theta.cc.hat)}
        return(-jacobian(S.theta.cc,theta.cc.hat))
    }))
    

    # X <- X.N[d$R2==1,]
    # Y <- matrix(Y.N[d$R2==1,],ncol=1)
    # theta.hat.cc <- solve(t(X)%*%X)%*%t(X)%*%Y
    # sigma.hat.cc <- sqrt(sum((Y-X%*%theta.hat.cc)^2)/(nrow(X)-4))
    # sigma(glm(analysis.model.formula,data=dset))
    # sqrt(diag(sigma.hat.cc^2*solve(t(X)%*%X)))
    # summary(glm(analysis.model.formula,data=dset))$coef[,2]
    
    se.cc <-sqrt(diag(1/n * solve(tau.cc.hat) %*% S.theta.cc %*% (t(S.theta.cc)/n) %*% t(solve(tau.cc.hat))))
    
    # Estimate alpha hat
    S.alpha.multiroot <- function(alpha.hat){
      rowMeans(do.call("cbind",lapply(1:N,function(i){
        Di <- get.Di(i)
        return(S.alpha.i(Di, alpha.hat))
    })))}
    
    alpha.hat <- matrix(data = multiroot(S.alpha.multiroot,rep(0,ncol(get.Di(1)$Z1)))$root, ncol=1)
    
    d.to.impute <- dset[dset$R1==1 & dset$R2==0,]
    d.not.imputed <- dset[!(dset$R1==1 & dset$R2==0),]
    
    # Linear imputation model
    
    if (imp.method == "linear"){
      # Estimating equation for beta and sigma (collectively psi) (imputation model)
      S.psi.i <- function(Di, psi.hat){
        beta.hat <- psi.hat[1:(length(psi.hat)-1)]
        sigma.hat <- psi.hat[length(psi.hat)]
        return(rbind(1/sigma.hat^2 * t(Di$Z2) %*% (Di$A - Di$Z2 %*% beta.hat),
                     -1/sigma.hat+((Di$A - Di$Z2 %*% beta.hat)^2/sigma.hat^3)))
      }
      
      # Solve for psi hat (beta hat and sigma hat)
      imp.fit <- glm(imp.model.formula, data = dset)
      beta.hat <- matrix(coef(imp.fit), ncol = 1)
      sigma.hat <- sigma(imp.fit)
      psi.hat <- matrix(c(beta.hat, sigma.hat), ncol=1)
      
      imputed <- do.call("rbind",lapply(1:M,function(j){
        d.imputed <- d.to.impute
        d.imputed[, A.var] <- as.numeric(model.matrix.lm(imp.model.formula, data=d.imputed,na.action='na.pass') %*% beta.hat) + rnorm(nrow(d.imputed),0,sigma.hat)
        stacked <- rbind(d.imputed, d.not.imputed)
        stacked$j <- j
        return(stacked[order(stacked$i,stacked$j),])
      }))
    } else if (imp.method == "logistic") {
      # Logistic imputation model
      
      # Estimating equation for beta and sigma (collectively psi) (imputation model)
      S.psi.i <- function(Di, psi.hat){
        return(t(Di$Z2) %*% (Di$A - Di$Z2 %*% beta.hat))
      }
      
      # Solve for psi.hat
      imp.fit <- glm(imp.model.formula, data = dset, family="binomial")
      psi.hat <- matrix(coef(imp.fit), ncol=1)
      
      imputed <- do.call("rbind",lapply(1:M,function(j){
        d.imputed <- d.to.impute
        d.imputed[, A.var] <- rbinom(nrow(d.imputed), 1, predict(imp.fit, d.imputed, type = "response"))
        stacked <- rbind(d.imputed, d.not.imputed)
        stacked$j <- j
        return(stacked[order(stacked$i,stacked$j),])
      }))
    } else {
      print("Must specify imp.method as linear or logistic")
      stop()
    }
    
    X.NM <- model.matrix.lm(analysis.model.formula, data=imputed, na.action='na.pass')
    Y.NM <- matrix(model.response(model.frame(analysis.model.formula, data=imputed,na.action=NULL)), ncol=1)
    Z1.NM <- model.matrix.lm(ipw.model.formula,data=imputed, na.action='na.pass')
    R1.NM <- matrix(imputed[,R1.var],ncol=1)
    Z2.NM <- model.matrix.lm(imp.model.formula,data=imputed, na.action='na.pass')  
    A.NM <- matrix(imputed[,A.var], ncol=1)
    R2.NM <- matrix(imputed[,"R2"], ncol=1)
    
    get.Dij <- function(i,j){
      ij <- (j-1)*N + i
      list(X = matrix(X.NM[ij,],nrow=1),
           Y = Y.NM[ij],
           Z1 = matrix(Z1.NM[ij,],nrow=1),
           Z2 = matrix(Z2.NM[ij,],nrow=1),
           R1 = R1.NM[ij],
           R2 = R2.NM[ij],
           A = A.NM[ij])
    }
    
    # IPW/MI estimating equation for theta
    S.theta.ipwmi.i <- function(Di, theta.hat, alpha.hat){
      Wi <- W.i(Di, alpha.hat)
      if (Wi == 0){
        return(matrix(data = 0, nrow = ncol(Di$X), ncol = 1))
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
    
    theta.hat <- matrix(data = multiroot(S.theta.multiroot,rep(0,ncol(get.Di(1)$X)))$root, ncol=1)
    
    # Perform inference
    S.theta <- do.call("cbind",lapply(1:M,function(j){
      do.call("cbind",lapply(1:N,function(i){
        S.theta.ipwmi.i(get.Dij(i,j), theta.hat, alpha.hat)
      }))
    }))
    
    S.psi.obs.i <- function(Di, psi.hat){
      if (Di$R2 == 0){
        return(matrix(data = 0, nrow = ncol(Di$Z2)+1, ncol = 1))
      } else {
        return(S.psi.i(Di, psi.hat))
      }
    }
    
    S.psi.mis.i <- function(Di, psi.hat){
      if (Di$R1 == 1 & Di$R2 == 0){ 
        return(S.psi.i(Di, psi.hat))
      } else {
        return(matrix(data = 0, nrow = ncol(Di$Z2) + 1, ncol = 1))
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


    # Rubin's rules
    
    theta.hats <- lapply(1:M, function(j){
      S.theta.j.multiroot <- function(theta.hat){
        rowMeans(do.call("cbind",lapply(1:N,function(i){
            Dij <- get.Dij(i,j)
            return(S.theta.ipwmi.i(Dij, theta.hat, alpha.hat))
        })))
      }
      
      return(matrix(data = multiroot(S.theta.j.multiroot,rep(0,ncol(get.Di(1)$X)))$root, ncol=1))
      
    })
    
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

get.ipwmi.se(dset=d.obs,
              analysis.model.formula = (Y-X4)/X4 ~ X2*X3,
              ipw.model.formula = R1 ~ X1,
              imp.model.formula = Y ~ X1*X2*X3 + X4 + X5,
              S.theta.i = function(Di,theta.hat){t(Di$X)%*%(Di$Y-Di$X%*%theta.hat)},
              S.alpha.i = function(Di,alpha.hat){t(Di$Z1)%*%(Di$R1-Di$Z1%*%alpha.hat)},
              W.i = function(Di,alpha.hat){Di$R1/as.numeric(Di$Z1%*%alpha.hat)},
              imp.method = "linear")

dset<-d.obs
analysis.model.formula <- (Y-X4)/X4 ~ X2*X3
ipw.model.formula <- R1 ~ X1
imp.model.formula <- Y ~ X1*X2*X3 + X4 + X5
S.theta.i <- function(Di,theta.hat){t(Di$X)%*%(Di$Y-Di$X%*%theta.hat)}
S.alpha.i <- function(Di,alpha.hat){t(Di$Z1)%*%(Di$R1-Di$Z1%*%alpha.hat)}
W.i <- function(Di,alpha.hat){Di$R1/as.numeric(Di$Z1%*%alpha.hat)}
imp.method <- "linear"



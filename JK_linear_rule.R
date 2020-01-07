list.of.packages <- c("ranger", "DTRreg", "parallel", "gam","rgenoud")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)

library(ranger)
library(parallel)
library(gam)
library(rgenoud)

fit<-function(x, y, A, psm=2, Qm=2){
  x<-as.matrix(x)
  n<-length(y)
  loc.a1<-which(A==1)
  loc.a0<-which(A==0)
  dat <- data.frame(cbind(y, x))
  dat2 <- data.frame(cbind(y, x, A,  A*x))
  
  ## use random forest to est Q
  lm.xnew<- x
  lm.y<- y[loc.a1]
  lm.x<- x[loc.a1,]
  lm.x<- as.matrix(lm.x)
  # Gam.object <- gam(y ~ s(x1) + s(x2) + s(x3) + s(x4), data=dat[loc.a1, ])
  # reg.mu1.gam = predict(Gam.object, type="response", newdata=dat)
  fit1 <- ranger(y~., data=dat[loc.a1,])
  reg.mu1.rf <- predict(fit1, data = dat)$predictions
  
  lm.y<- y[loc.a0]
  lm.x<- x[loc.a0,]
  lm.x<- as.matrix(lm.x)
  # Gam.object <- gam(y ~ s(x1) + s(x2) + s(x3) + s(x4), data=dat[loc.a0, ])
  # reg.mu0.gam = predict(Gam.object, type="response", newdata=dat)
  fit0 <- ranger(y~., data=dat[loc.a0,])
  reg.mu0.rf <- predict(fit0, data = dat)$predictions
  
  ## aipw
  if (psm==1) {
    phat <- mean(A)
  } else if (psm==2) {
    ## use logistic regression to est PS
    alpha.est<-glm(A~x,family="binomial")$coeff
    lpshat<-cbind(1,x)%*%alpha.est
    phat<-1/(1+exp(-lpshat ))
  } else if (psm==3) {
    ## use gam to est PS
    Gam.object <- gam(A ~ s(x1) + s(x2) + s(x3) + s(x4), family=binomial, data=dat2)
    phat = predict(Gam.object, type="response", newdata=dat2)
  }
  
  lin <- NA
  if (Qm==1) {
    reg.mu1 <- reg.mu1.rf
    reg.mu0 <- reg.mu0.rf
  } else if (Qm==2) {
    ## use linear model to est Q in aipw
    lin <- lm(y ~ ., data=dat2)
    newdf <- data.frame(x, rep(0, n), x*0)
    names(newdf) <- names(dat2)[-1]
    reg.mu0.lin <- predict(lin, newdata = newdf)
    newdf <-  data.frame(x, rep(1, n), x)
    names(newdf) <- names(dat2)[-1]
    reg.mu1.lin <- predict(lin, newdata = newdf)
    
    reg.mu1 <- reg.mu1.lin
    reg.mu0 <- reg.mu0.lin
  } else if (Qm==3) {
    ## use GAM to est Q in aipw
    reg.mu1 <- reg.mu1.gam
    reg.mu0 <- reg.mu0.gam
  }
  psi.aug <- y*A/phat-y*(1-A)/(1-phat)+(reg.mu1*(1-A/phat)-reg.mu0*(1-(1-A)/(1-phat)) )
  aipw<-mean(psi.aug)
  ps_rule <- phat*A + (1-phat)*(1-A)
  return(list(aipw = aipw, lin = lin, ps_rule = ps_rule, phat=phat,
              reg.mu0 = reg.mu0.rf, reg.mu1 = reg.mu1.rf))
}

main <- function(seed, n, psm=2, prob=0.2) {
  set.seed(seed)
  
  x1<-rnorm(n)
  x2<-rnorm(n)
  x3<-rnorm(n)
  x4<-rnorm(n)
  # true.contrast <-  1 + cos(x1) - 2*cos(x2) + cos(x3) - 2*cos(x4) 
  true.contrast <- (x1 > 0.5) * (x2 < -0.5) + 2*abs(x3) - 2 + (cos(x4) < 0.5)
  x<-as.matrix(cbind(x1,x2,x3,x4))
  
  ## generate A
  if (psm == 1) {
    prob <- prob
  } else if(psm==2) {
    lps <- x1/2 + x2/2-x3/4 +x4/4
    prob <- exp(lps)/(1 + exp(lps))
  } else if(psm==3) {
    lps <- (x1+x2)^2/2 -(x3+x4)^2/2
    prob <- exp(lps)/(1 + exp(lps))
  }
  A <- rbinom(n, 1, prob)
  
  y1<- (x1+x2+x3)/2+rnorm(n,0,1) + true.contrast
  y0<- (x1+x2+x3)/2+rnorm(n,0,1)
  y<-y1*A+y0*(1-A)
  tau <- mean(y1-y0)
  method <- fit(x,y,A)
  aipw <- method$aipw
  
  ## compare method 1: E(Y|X,A=1)-E(Y|X,A=0) as the outcome to learn the linear contrast function
  diffs <- method$reg.mu1 - method$reg.mu0
  coef1 <- lm(diffs~x1+x2)$coef
  
  # ## compare method 2: A-learning via G-estimation
  # # blip model
  # blip.mod <- list(~x1+x2)
  # # treatment model (correctly specified)
  # treat.mod <- list(A1~1)
  # # treatment-free model (incorrectly specified)
  # tf.mod <- list(~x1+x2+x3+x4)
  # # perform G-estimation
  # mod1 <- DTRreg(Y, blip.mod, treat.mod, tf.mod, method = "gest")
  # coef2 <- mod1$psi[[1]]
  
  ## compare method 3: value search using genetic algorithm
  lin <- method$lin
  ps_rule <- method$ps_rule
  aipw_d <- function(betas) {
    rule <- as.integer(cbind(rep(1,n), x[,1:(length(betas)-1)]) %*% betas > 0)
    Cd = (rule==A)
    Q.hat <- predict(lin, newdata = data.frame(x, A=rule, rule*x))
    Vd <- mean(Cd*y/ps_rule - (Cd/ps_rule-1)*Q.hat)
    return(Vd)
  }
  coef3 <- genoud(aipw_d, nvars=3, max=TRUE)$par
  
  ## compare method 4: adjusted-outcome contrast estimator
  phat <- method$phat
  adj_out <- (A/phat - (1-A)/(1-phat))*y
  coef4 <- lm(adj_out~x1+x2)$coef
  
  ## proposed method
  tau_jk <- rep(0, n)
  for (i in 1:n) {
    out <- fit(x[-i,],y[-i],A[-i])
    tau_jk[i] <- n*aipw-(n-1)* out$aipw
  }
  coef5 <- lm(tau_jk~x1+x2)$coef
  
  coef <- lm(true.contrast~x1+x2)$coef
  ## results
  coefs <- as.matrix(cbind(coef1,coef3, coef4, coef5, coef))
  return(list(coefs=coefs))
}

quality <- function(N, res, X, true.contrast, true.opt.action, n) {
  ## to calcualte the true value of the learned rules and the classification accuracy 
  accuracies <- matrix(NA, length(ranseed), 5)
  non_baselines <- matrix(NA, length(ranseed), 6)
  for (rep in 1:length(res)) {
    coefs <- res[[rep]]$coefs
    contrast <- cbind(rep(1,N), X[,1:(dim(coefs)[1]-1)]) %*% coefs 
    action <- (contrast > 0)
    non_baseline <- colMeans(true.contrast*cbind(action, true.opt.action))
    accuracy <- c()
    for (j in 1:dim(action)[2]) {
      accuracy <- c(accuracy, mean(action[,j]==true.opt.action))
    }
    accuracies[rep,] <- accuracy
    non_baselines[rep,] <- non_baseline
  }
  print(paste0("n=", n))
  print(paste("accrcy1","accrcyValueSearch", "accrcyAdjOutcome", "accrcyJK", "accrcyTrue"))
  print(apply(accuracies, 2, mean))
  print(paste("value1", "valueValueSearch", "valueAdjOutcome", "vlaueJK", "value_true_cont", "opt_value"))
  print(apply(non_baselines, 2, mean))
  return(list(accuracies=accuracies, non_baselines=non_baselines))
}


numWorkers <- 56
ranseed <- 1:(4*56)
## n=100
n=200
res <- mclapply(ranseed, main, n=n, mc.cores=numWorkers)

N=5000000
set.seed(1)
x1<-rnorm(N)
x2<-rnorm(N)
x3<-rnorm(N)
x4<-rnorm(N)
# true.contrast <-  1 + cos(x1) - 2*cos(x2) + cos(x3) - 2*cos(x4) 
true.contrast <- (x1>0.5)*(x2<-0.5) + 2*abs(x3) - 2 + (cos(x4)<0.5)
true.opt.action <- (true.contrast > 0 )
X<-as.matrix(cbind(x1,x2,x3,x4))
rm(x1);rm(x2);rm(x3);rm(x4)
results <- quality(N=N, res=res, X=X, true.contrast, true.opt.action, n=n)


# taus <- matrix(NA, length(ranseed), 3)
# accuracies <- matrix(NA, length(ranseed), 3)
# non_baselines <- matrix(NA, length(ranseed), 4)
# for (rep in ranseed) {
#   taus[rep,] <- res[[rep]]$taus
#   accuracies[rep,] <- res[[rep]]$accuracy
#   non_baselines[rep,] <- res[[rep]]$non_baseline
# }
# print(paste0("n=", n))
# print(paste("est.tau1","est.tauJK", "true.tau"))
# print(apply(taus, 2, mean))
# print(paste("accrcy1","accrcyJK", "accrcyTrue"))
# print(apply(accuracies, 2, mean))
# print(paste("value1","vlaueJK", "value_true_cont", "opt_value"))
# print(apply(non_baselines, 2, mean))

# par(mfrow=c(1,2))
# taus <- data.frame(taus)
# names(taus) <- c("true.tau","est.tau1","est.tauJK")
# boxplot(taus)
# 
# accuracies <- data.frame(accuracies)
# names(accuracies) <- c("accrcy1","accrcyJK", "accrcyTrue")
# boxplot(accuracies)
par(mfrow=c(2,2))
for (n in c(200, 500, 1000, 2000)) {
  filename <- paste0("/Users/lwu9/Documents/M2C/n",n,".RData")
  load(filename)
  df <- data.frame(results[[1]])
  names(df) <- c("Q","VS", "AO", "JK", "TC")
  boxplot(df, ylab="Accuracy", main=paste0("n=",n))
  df <- data.frame(results[[2]])
  names(df) <- c("Q","VS", "AO", "JK", "TC", "TO")
  boxplot(results[[2]], ylab="Value", main=paste0("n=",n))
}
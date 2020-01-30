#source('S:\\Documents\\lab2018\\transfer learning of dtr\\functions.R')
source("/home/lwu9/transfer/functions.R")
# path <- "/Users/lwu9/Documents/transfer_learner_dtr"; setwd(path); source(paste0(path,'/functions.R'))

#### benchmark: use true w.tilde and true contrast values
benchmark <- function(m, alphas, seed, N, mu, p, method, Qm, y.hat.meth=1, nfold=30, specify=F) {
  set.seed(seed)
  # x <- cbind(rep(1, N), matrix(runif(p*N, min = -1, max = 1), N, p))
  x <- cbind(rep(1, N), rmvnorm(n=N, mean=rep(mu, p), sigma=diag(rep(1, p))))
  # x <- cbind(rep(1, N), rmvnorm(n=N, mean=rep(mu, p), sigma=cov_matrix))
  prob.one = 0.5 # P(A=1)
  a <- rbinom(N, 1, prob.one)
  if (specify==T) {
    # correct specification: true trt effect model is linear
    baseline <- x %*% betas1
    contrast <- x %*% betas2
  } else {
    # misspecification 
    baseline <- x %*% betas1
    # contrast <- x**3 %*% betas2
    # contrast <- x**2 %*% betas2 
    contrast <- (cos(x)-0.5) %*% betas2 ## dist are different in rwe and rct
    # contrast <- betas2[1]*(cos(x[,1])-0.5) + betas2[2]*(cos(x[,2])-0.5) +betas2[3]*(cos(x[,3])-0.5)
    # contrast <- atan(exp(1+x[,2])-3*x[,3]-5) ## dist are similar in rwe and rct
    # contrast <- 3 * ((x[,2]<=1) & (x[,3]>-0.6)) - 1
  }
  y1.star <-  baseline + contrast + rnorm(N, sd=sigma)
  y0.star <-  baseline + rnorm(N, sd=sigma)
  opt.actions <- (sign(contrast)+1)/2
  opt.value <-  baseline + opt.actions*(contrast)
  y.star <- y1.star - y0.star
  # rct
  delta.tilde <- rbinom(N, 1, delta.tilde.prob(x, alphas))
  (n <- sum(delta.tilde)) # sample size in rct
  # print(paste(mean(opt.actions), n))
  delta.tilde.is.one.ind <- which(delta.tilde==1)
  w.tilde <- (1/delta.tilde.prob(x, alphas))[delta.tilde.is.one.ind]
  x.in.rct <- x[delta.tilde.is.one.ind, ]
  y <-  rep(0, N); y[which(a==1)] <- y1.star[which(a==1)]; y[which(a==0)] <- y0.star[which(a==0)] # observed outcomes in RCT
  y.in.rct <- y[delta.tilde.is.one.ind]
  a.in.rct <- a[delta.tilde.is.one.ind]
  # y.star.hat <- est.contrast(a, x, y, delta.tilde.is.one.ind, meth=y.hat.meth, psm=1,Qm=Qm, regr.direct=F)
  # y.star.hat.in.rct <- y.star.hat[delta.tilde.is.one.ind]
  action <- learn_rule_given_w.tilde(y.in.rct, x.in.rct, a.in.rct, y.star.hat.in.rct=contrast[delta.tilde.is.one.ind], 
                                     x=x, method=method, prob.one=NA, w.tilde=w.tilde)
  (evaluation <- quality(action, opt.actions, opt.value, baseline, contrast))
  true.values.mse.rule <- evaluation[2]
  # # rwe
  # delta.is.one.ind <- sample(N, size=m)
  # x.in.rwe <- x[delta.is.one.ind, ]
  # ## visualization for covariates distributions in RCT population, and RWE
  # # par(mfrow=c(p+1, 3))
  # # for (hi in 0:p) {
  # #   hist(x.in.rct[,hi+1])
  # #   hist(x.in.rwe[,hi+1])
  # #   hist(x.in.rct[,hi+1]*w.tilde)
  # #   print(c(mean(x.in.rct[,hi+1]), mean(x.in.rwe[,hi+1]), mean(x.in.rct[,hi+1]*w.tilde)))
  # # }
  # true.values.mse.rule <- c()
  # for (w.method in 4:1) {
  #   action <- learn_rule(y.in.rct, x.in.rct, a.in.rct, x.in.rwe, x, N, y.star.hat.in.rct, 
  #                        w.method, method, prob.one=prob.one, only.weight=F)$action 
  #   (evaluation <- quality(action, opt.actions, opt.value, baseline, contrast))
  #   # print(paste(method, w.method, ": ", round(evaluation[1],3), round(evaluation[2],3)))
  #   true.values.mse.rule <- c(true.values.mse.rule, evaluation[2])
  # }
  return(list(true.values.mse.rule=true.values.mse.rule, n=n))
}


########################################### Simulations #####################################################
N=1000000
p <- 2; 

sigma=0.5; 
# alphas <- as.matrix(c(-8, 1), rep(0, p-1), p+1, 1)
betas1 <- betas2 <- betas <- as.matrix(1:(p+1), p+1, 1)
# alphas <- as.matrix(c(-10, 1, 1, -2, rep(0, p-3)), p+1, 1)
alphas <- as.matrix(c(-10, 1,-2, rep(0, p-2)), p+1, 1)
# betas <- betas1 <- betas2 <- as.matrix(c(-1/2,1), p+1, 1)
# betas1 <- c(2, 1, 0, 1, 0, rep(0, p-4)); betas2 <- c(0, 1, -1, 1, -1, rep(0, p-4))
seed=1; mu=1
mus <- c(1)
nfold = 5; methods=c("linear.VS")

for (method in methods) {
  for (alpha1 in c(-10,-8))   {
    Results <- list()
    if (alpha1==-10) rs <- 1:(56*5) else
      rs <- 1:(56*8)
    for (m in c(5000)) {
      for (num.mu in 1:length(mus)) {
        mu <- mus[num.mu]
        alphas[1] <- alpha1
        registerDoParallel(56) 
        results <- foreach(seed=rs, .combine=rbind, .errorhandling=c('remove'),
                           .packages=c('mvtnorm','ranger', 'rpart', 'randomForest','CVXR','rgenoud')) %dopar% {
                             benchmark(m, alphas, p=p, seed = seed, method=method, Qm=NA, y.hat.meth=NA, N=1000000, mu=mu, nfold=NA, specify = F)
                           }
        stopImplicitCluster()
        print(paste("benchmark",mean(unlist(results[,1])), "; method",method, "; alpha=",alpha1))
        Results[paste0("alpha",alpha1,"yhatmeth0Qm0")] <- list(results)
      }
    }
    save(Results, file=paste0("results_cv_nfoldNA", "_method",method,"_alpha",abs(alpha1),"_diff_benchmark.Rdata"))
  }
}
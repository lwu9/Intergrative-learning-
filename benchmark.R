source("/home/lwu9/transfer/functions3.R")
# path <- "/Users/lwu9/Documents/transfer_learner_dtr"; setwd(path); source(paste0(path,'/functions3.R'))

#### benchmark: use true w.tilde and true contrast values
benchmark <- function(m, alphas, seed, N, mu, p, specify=F, dbn) {
  set.seed(seed)
  x <- cbind(rep(1, N), rmvnorm(n=N, mean=rep(mu, p), sigma=diag(rep(1, p))))
  prob.one = 0.5 # P(A=1)
  a <- rbinom(N, 1, prob.one)
  if (specify==T) {
    # correct specification: true trt effect model is linear
    baseline <- x %*% betas1
    contrast <- x %*% betas2
  } else {
    # misspecification 
    baseline <- x %*% betas1
    if (dbn == "diff") {
      contrast <- (cos(x)-0.5) %*% betas2 ## dist are different in rwe and rct
    } else {
      contrast <- atan(exp(1+x[,2])-3*x[,3]-5) ## dist are similar in rwe and rct
    }
  }
  y1.star <-  baseline + contrast + rnorm(N, sd=sigma)
  y0.star <-  baseline + rnorm(N, sd=sigma)
  opt.actions <- (sign(contrast)+1)/2
  opt.value <-  baseline + opt.actions*(contrast)
  y.star <- y1.star - y0.star
  # rct
  delta.tilde <- rbinom(N, 1, delta.tilde.prob(x, alphas))
  (n <- sum(delta.tilde)) # sample size in rct
  delta.tilde.is.one.ind <- which(delta.tilde==1)
  w.tilde <- (1/delta.tilde.prob(x, alphas))[delta.tilde.is.one.ind]
  x.in.rct <- x[delta.tilde.is.one.ind, ]
  y <-  rep(0, N); y[which(a==1)] <- y1.star[which(a==1)]; y[which(a==0)] <- y0.star[which(a==0)] # observed outcomes in RCT
  y.in.rct <- y[delta.tilde.is.one.ind]
  a.in.rct <- a[delta.tilde.is.one.ind]
  beta.hat.star <- get_eta(x.in.rct, w.tilde=w.tilde/sum(w.tilde), y.star.hat.in.rct=contrast[delta.tilde.is.one.ind])
  action <- as.numeric(x %*% beta.hat.star > 0)
  # action <- learn_rule_given_w.tilde(y.in.rct, x.in.rct, a.in.rct, y.star.hat.in.rct=contrast[delta.tilde.is.one.ind], 
  #                                    x=x, method=method, prob.one=NA, w.tilde=w.tilde)$action
  (evaluation <- quality(action, opt.actions, opt.value, baseline, contrast))
  true.values.mse.rule <- evaluation[2]
  return(list(true.values.mse.rule=true.values.mse.rule, n=n))
}


########################################### Simulations #####################################################
N=1000000
p <- 2; 

sigma=0.5; 
betas1 <- betas2 <- betas <- as.matrix(1:(p+1), p+1, 1)
alphas <- as.matrix(c(-10, 1,-2, rep(0, p-2)), p+1, 1)
mus <- c(1); rs <- 1:200
Results <- list()
for (dbn in c("similar", "diff")) {
  for (alpha1 in c(-10,-9,-8))   {
    for (m in c(5000)) {
      for (num.mu in 1:length(mus)) {
        mu <- mus[num.mu]
        alphas[1] <- alpha1
        registerDoParallel(67) 
        results <- foreach(seed=rs, .combine=rbind, .errorhandling=c('remove'),
                           .packages=c('mvtnorm','ranger', 'rpart', 'randomForest','CVXR','rgenoud')) %dopar% {
                             benchmark(m=m, alphas=alphas, seed, N=N, mu=mu, p=p, specify=F, dbn=dbn)
                           }
        stopImplicitCluster()
        print(paste("benchmark",mean(unlist(results[,1])), "; n", round(mean(unlist(results[,2]))),"; dbn:", dbn, "; alpha=",alpha1))
        Results[paste0("alpha",alpha1,"dbn",dbn)] <- list(results)
      }
    }
  }
}
save(Results, file=paste0("results_benchmark.Rdata"))


#source('S:\\Documents\\lab2018\\transfer learning of dtr\\functions.R')
source("/home/lwu9/transfer/functions.R")

cv.compare <- function(m, alphas, seed, N, mu, p, nfold=50, specify=F) {
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
    # contrast <- (cos(x)-0.5) %*% betas2
    # contrast <- betas2[1]*(cos(x[,1])-0.5) + betas2[2]*(cos(x[,2])-0.5) +betas2[3]*(cos(x[,3])-0.5)
    contrast <- atan(exp(1+x[,2])-3*x[,3]-5)
    # contrast <- 3 * ((x[,2]<=1) & (x[,3]>-0.6)) - 1
  }
  y1.star <-  baseline + contrast + rnorm(N, sd=sigma)
  y0.star <-  baseline + rnorm(N, sd=sigma)
  opt.actions <- (sign(contrast)+1)/2
  opt.value <-  baseline + opt.actions*(contrast)
  y.star <- y1.star - y0.star
  # opt <- optim(rep(2,p+1), cost3, y.hat.in.rct=y.star, x.in.rct=x, w.tilde=rep(1, N), control = list(maxit=10000))
  # beta.N.star <- opt$par
  # print(beta.N.star)
  (beta.N.star <- solve(t(x) %*% x) %*% as.matrix(apply(x * c(y.star), 2, sum), p+1, 1))
  # rct
  delta.tilde <- rbinom(N, 1, delta.tilde.prob(x, alphas))
  (n <- sum(delta.tilde)) # sample size in rct
  # print(paste(mean(opt.actions), n))
  delta.tilde.is.one.ind <- which(delta.tilde==1)
  x.in.rct <- x[delta.tilde.is.one.ind, ]
  y <-  rep(0, N); y[which(a==1)] <- y1.star[which(a==1)]; y[which(a==0)] <- y0.star[which(a==0)] # observed outcomes in RCT
  y.in.rct <- y[delta.tilde.is.one.ind]
  a.in.rct <- a[delta.tilde.is.one.ind]
  y.star.hat <- a*y/prob.one - (1-a)*y/(1-prob.one) # outcome adjusted
  y.star.hat.in.rct <- y.star.hat[delta.tilde.is.one.ind]
  # rwe
  delta.is.one.ind <- sample(N, size=m)
  x.in.rwe <- x[delta.is.one.ind, ]
  ## visualization for covariates distributions in RCT population, and RWE
  # par(mfrow=c(p+1, 3))
  # for (hi in 0:p) {
  #   hist(x.in.rct[,hi+1])
  #   hist(x.in.rwe[,hi+1])
  #   hist(x.in.rct[,hi+1]*w.tilde)
  #   print(c(mean(x.in.rct[,hi+1]), mean(x.in.rwe[,hi+1]), mean(x.in.rct[,hi+1]*w.tilde)))
  # }
  true.values.mse.rule <- c()
  for (w.method in 1:2) {
    for (method in c("linear", "tree")) {
      action <- learn_rule(y.in.rct, x.in.rct, a.in.rct, x.in.rwe, x, N, y.star.hat.in.rct, 
                          w.method, method, prob.one=prob.one, only.weight=F)$action 
      (evaluation <- quality(action, opt.actions, opt.value, baseline, contrast))
      # print(paste(method, w.method, ": ", round(evaluation[1],3), round(evaluation[2],3)))
      true.values.mse.rule <- c(true.values.mse.rule, evaluation[2])
    } 
  }
  (cv.results <- cv(y.in.rct, x.in.rct, a.in.rct, x.in.rwe, N, nfold=nfold, y.star.hat.in.rct,
                    prob.one, ipw=F, baseline, contrast, x, calculate.true.val=F))
  return(list(true.values.mse.rule=true.values.mse.rule, cv.value.rule=cv.results$values, n=n))
}

########################################### Simulations #####################################################
N=1000000
m=8000 # sample size in rwe
p <- 3; 
 
sigma=0.5; 
# alphas <- as.matrix(c(-8, 1), rep(0, p-1), p+1, 1)
betas1 <- betas2 <- betas <- as.matrix(1:(p+1), p+1, 1)
alphas <- as.matrix(c(-10, 1, 1, -2, rep(0, p-3)), p+1, 1)
# betas <- betas1 <- betas2 <- as.matrix(c(-1/2,1), p+1, 1)
# betas1 <- c(2, 1, 0, 1, 0, rep(0, p-4)); betas2 <- c(0, 1, -1, 1, -1, rep(0, p-4))
seed=1; mu=1
mus <- c(1)
Results <- list()
system.time(for (alpha1 in c(-10,-9,-8,-7,-6))   {
  for (m in c(8000)) {
    for (num.mu in 1:length(mus)) {
      mu <- mus[num.mu]
      alphas[1] <- alpha1
      r <- 48*4; 
      rs <-  1:r
      registerDoParallel(48) 
      results <- foreach(seed=rs, .combine=rbind, .errorhandling=c('remove'),
                         .packages=c('mvtnorm','ranger', 'rpart', 'randomForest')) %dopar% {
                           cv.compare(m, alphas, p=p, seed = seed,N=1000000, mu=mu, specify = F)
                         }
      stopImplicitCluster()
      # print(head(Reduce(rbind, results[,1])))
      # print(learned.rule <- apply(Reduce(rbind, results[,1]), 1, which.min))
      # orders <- apply(Reduce(rbind, results[,1]), 1, order)
      # print(cv.rule <- apply(Reduce(rbind, results[,2]), 1, which.max))
      learned.rule.mse <- matrix(unlist(results[,1]), byrow = T, ncol=4)
      cv.rule.value.est <- matrix(unlist(results[,2]), byrow = T, ncol=4)
      orders <- apply(learned.rule.mse, 1, order)
      cv.rule <- apply(cv.rule.value.est, 1, which.max)
      print(mean(unlist(results[,3])))
      print(dim(cv.rule.value.est))
      print(paste("################ correct rate using cv:", mean(cv.rule==orders[1,])))
      print(paste("################ correct rate using cv 2nd:", mean(cv.rule==orders[2,])))
      print(paste("################ correct rate using cv 3rd:", mean(cv.rule==orders[3,])))
      print(paste("################ correct rate using cv 4th:", mean(cv.rule==orders[4,])))
      Results[paste0("alpha",alpha1)] <- list(results)
     }
  }
})
save(Results, file="results_cv.Rdata")

source('S:\\Documents\\lab2018\\transfer learning of dtr\\functions.R')
N=1000000
m=8000 # sample size in rwe
p <- 3; 
# cov_matrix <- matrix(0, nrow = p, ncol = p)
# for (j1 in 1:p) {
#   for (j2 in 1:p) {
#     cov_matrix[j1, j2] <- 0.2**abs(j1-j2)*4
#   }
# }

sigma=0.5; mu=0
alphas=as.matrix(c(-9, 1, -2, 1), p+1, 1)
betas <- as.matrix(1:(p+1), p+1, 1)

beta_alpha <- function(m, alphas, seed, N, mu, specify=F) {
  set.seed(seed)
  x <- cbind(rep(1, N), rmvnorm(n=N, mean=rep(mu, p), sigma=diag(rep(1, p))))
  prob.one = 0.5 # P(A=1)
  a <- rbinom(N, 1, prob.one)
  if (specify==T) {
    # correct specification: true trt effect model is linear
    y1.star <-  x %*% betas + x %*% betas + rnorm(N, sd=sigma)
    y0.star <-  x %*% betas + rnorm(N, sd=sigma)
    opt.actions <- (sign(x%*%betas)+1)/2
  } else {
    # misspecification 
    # y1.star <-  x %*% betas + x**3 %*% betas + rnorm(N, sd=sigma)
    y1.star <- x %*% betas + 3 * ((x[,2]<=1) & (x[,3]>-0.6)) - 1 + rnorm(N, sd=sigma)
    y0.star <-  x %*% betas + rnorm(N, sd=sigma)
    # opt.actions <- (sign(x**3 %*% betas)+1)/2
    opt.actions <- (sign( 3 * ((x[,2]<=1) & (x[,3]>-0.6)) - 1)+1)/2
  }
  # opt.value <-  x %*% betas + opt.actions * (x**3 %*% betas)
  opt.value <- x %*% betas + opt.actions* (3 * ((x[,2]<=1) & (x[,3]>-0.6)) - 1)
  y.star <- y1.star - y0.star
  beta.N.star <- solve(t(x) %*% x) %*% as.matrix(apply(x * c(y.star), 2, sum), p+1, 1)
  # rct
  delta.tilde <- rbinom(N, 1, delta.tilde.prob(x, alphas))
  (n <- sum(delta.tilde)) # sample size in rct
  print(paste(mean(opt.actions), n))
  delta.tilde.is.one.ind <- which(delta.tilde==1)
  x.in.rct <- x[delta.tilde.is.one.ind, ]
  y <-  rep(0, N); y[which(a==1)] <- y1.star[which(a==1)]; y[which(a==0)] <- y0.star[which(a==0)] # observed outcomes in RCT
  y.star.hat <- a*y/prob.one - (1-a)*y/(1-prob.one) # outcome adjusted
 
    # rwe
    delta.is.one.ind <- sample(N, size=m)
    x.in.rwe <- x[delta.is.one.ind, ]
    left.hand <- apply(x[delta.is.one.ind, ], 2, mean)
    
    # the second method to reweight; refer to "method(2)" in the draft
    #print(cost2(alphas, x.in.rct, x.in.rwe, N))
    opt2 <- optim(rep(0,p+1), cost2, x.in.rct=x.in.rct, x.in.rwe=x.in.rwe, N=N, method = "BFGS", control = list(maxit=10000))
    alphas.hat2 <- opt2$par
    w.tilde <- 1 + 1/(exp(x.in.rct %*% as.matrix(alphas.hat2, length(alphas.hat2), 1)))
    x.tilde <- x.in.rct * c(sqrt(w.tilde))
    beta.hat.star2 <- solve(t(x.tilde) %*% x.tilde) %*% as.matrix(apply(x.in.rct * y.star.hat[delta.tilde.is.one.ind]*c(w.tilde), 2, sum), n, 1)
    action2 <- (sign(x%*%beta.hat.star2)+1)/2
    rate2 <- mean(action2==opt.actions)
    y.temp <- rep(0,N)
    y.temp[action2==1] <- y1.star[action2==1]
    y.temp[action2==0] <- y0.star[action2==0]
    v2.mse <- mean((y.temp - opt.value)**2)
    
    # refer to "method(3)" in the draft
    # cost(alphas, x.in.rct, left.hand)
    opt <- optim(c(0,0,0,0), cost, x.in.rct=x.in.rct, left.hand=left.hand,N=N, method = "BFGS", control = list(maxit=10000))
    alphas.hat3 <- opt$par
    w.tilde <- 1 + 1/(exp(x.in.rct %*% as.matrix(alphas.hat3, length(alphas.hat3), 1)))
    x.tilde <- x.in.rct * c(sqrt(w.tilde))
    beta.hat.star3 <- solve(t(x.tilde) %*% x.tilde) %*% as.matrix(apply(x.in.rct * y.star.hat[delta.tilde.is.one.ind]*c(w.tilde), 2, sum), n, 1)
    action3 <- (sign(x%*%beta.hat.star3)+1)/2
    rate3 <- mean(action3==opt.actions)
    y.temp <- rep(0,N)
    y.temp[action3==1] <- y1.star[action3==1]
    y.temp[action3==0] <- y0.star[action3==0]
    v3.mse <- mean((y.temp - opt.value)**2)
  
    # results
    beta_diff2 <- mean(abs(beta.hat.star2-beta.N.star))
    alpha_diff2 <- mean(abs(alphas.hat2-alphas))
    beta_diff3 <- mean(abs(beta.hat.star3-beta.N.star))
    alpha_diff3 <- mean(abs(alphas.hat3-alphas))
    # unweighted
    w.tilde = rep(1, n)
    x.tilde <- x.in.rct * c(sqrt(w.tilde))
    beta.hat.star <- solve(t(x.tilde) %*% x.tilde) %*% as.matrix(apply(x.in.rct * y.star.hat[delta.tilde.is.one.ind]*c(w.tilde), 2, sum), n, 1)
    beta_diff <- mean(abs(beta.hat.star-beta.N.star))
    action <- (sign(x%*%beta.hat.star)+1)/2
    rate <- mean(action==opt.actions)
    y.temp <- rep(0,N)
    y.temp[action==1] <- y1.star[action==1]
    y.temp[action==0] <- y0.star[action==0]
    v.mse <- mean((y.temp - opt.value)**2)
    
    # the approximate optimal linear policy
    approx.opt.action <- (sign(x %*% beta.N.star) + 1)/2
    approx.opt.rate <- mean(approx.opt.action==opt.actions)
    y.temp <- rep(0,N)
    y.temp[approx.opt.action==1] <- y1.star[approx.opt.action==1]
    y.temp[approx.opt.action==0] <- y0.star[approx.opt.action==0]
    approx.opt.v.mse <- mean((y.temp - opt.value)**2)
    return(c(rate, rate2, rate3, approx.opt.rate, v.mse, v2.mse,v3.mse,approx.opt.v.mse, n, 
             beta_diff,beta_diff2,beta_diff3,alpha_diff2, alpha_diff3, mean(opt.actions)))
}
par(mfrow=c(2,2))
mus <- c(-1, -0.5, 0, 0.5, 1)
mus <- c(0)
system.time(for (alpha1 in c(-11, -10))   {
  for (m in c(8000)) {
    for (num.mu in 1:length(mus)) {
      mu <- mus[num.mu]
      alphas[1] <- alpha1
      r <- 90; 
      rs <- 1:r
      registerDoParallel(11) 
      results <- foreach(seed=rs, .combine=rbind,
                         .packages=c('mvtnorm')) %dopar% {
                           beta_alpha(m, alphas, seed = seed,N=1000000, mu=mu, specify = F)
                         }
      stopImplicitCluster()
      print(paste0("m=",m,", alpha1=",alpha1))
      print(rbind(round(apply(results, 2, mean), 2),
                  round(apply(results, 2, sd)/sqrt(length(rs)), 2)))
      data <- data.frame(value=c(results[,5], results[,6],results[,7] ,results[,8]),
                         accuracy=c(results[,1], results[,2], results[,3],results[,4]),
                         method <- c(rep("Unweighted", r), rep("Weighted(1)", r), rep("Weighted(2)", r),rep("Optimal Linear", r)))
      # name1 <- paste0("S:\\Documents\\lab2018\\transfer learning of dtr\\Alearning_alpha1_",abs(alpha1), "_mu", num.mu)
      # save(results, file=paste0(name1, ".RData"))
      boxplot(value~method, data=data, col="lightgray",
              main=paste0("Value MSE"," (m=",m, ", n=",round(mean(results[,9])) ,", mu=", mu,  ")"))
      boxplot(accuracy~method, data=data, col="lightgray",
              main=paste0("Optimal action accuracy"," (m=",m, ", n=",round(mean(results[,9])) ,", mu=", mu, ")"))
    } 
  }
})





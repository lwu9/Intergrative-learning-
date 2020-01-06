#source('S:\\Documents\\lab2018\\transfer learning of dtr\\functions.R')
source("/home/lwu9/transfer/functions.R")
# source('/Users/lwu9/Documents/transfer_learner_dtr/functions.R')

N=1000000
m=5000 # sample size in rwe
p <- 3; 
# cov_matrix <- matrix(0, nrow = p, ncol = p)
# for (j1 in 1:p) {
#   for (j2 in 1:p) {
#     cov_matrix[j1, j2] <- 0.2**abs(j1-j2)*4
#   }
# }

sigma=0.5; 
# alphas <- as.matrix(c(-8, 1), rep(0, p-1), p+1, 1)
betas1 <- betas2 <- betas <- as.matrix(1:(p+1), p+1, 1)
alphas <- as.matrix(c(-10, 1, 1, -2, rep(0, p-3)), p+1, 1)
# betas <- betas1 <- betas2 <- as.matrix(c(-1/2,1), p+1, 1)
# betas1 <- c(2, 1, 0, 1, 0, rep(0, p-4)); betas2 <- c(0, 1, -1, 1, -1, rep(0, p-4))
seed=1; mu=1

linear_tree <- function(m, alphas, seed, N, mu, specify=F) {
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
    # contrast <- betas2[1]*cos(x[,1]-0.5) + betas2[2]*cos(x[,2]-0.5) + betas2[3]*cos(x[,3]-0.5)
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
  delta.tilde <- rbinom(N, 1, delta.tilde.prob(x, alphas)) ## logistic correctly specify
  # delta.tilde <- rbinom(N, 1, delta.tilde.prob(cos(x-0.4), alphas)) ## logistic NOT correctly specify
  (n <- sum(delta.tilde)) # sample size in rct
  print(paste(mean(opt.actions), n))
  delta.tilde.is.one.ind <- which(delta.tilde==1)
  x.in.rct <- x[delta.tilde.is.one.ind, ]
  y <-  rep(0, N); y[which(a==1)] <- y1.star[which(a==1)]; y[which(a==0)] <- y0.star[which(a==0)] # observed outcomes in RCT
  # y.star.hat <- a*y/prob.one - (1-a)*y/(1-prob.one) # outcome adjusted
  y.star.hat <- est.contrast(a, x, y, delta.tilde.is.one.ind, prob.one, meth=2)
  # rwe
  delta.is.one.ind <- sample(N, size=m)
  x.in.rwe <- x[delta.is.one.ind, ]
  left.hand <- apply(x[delta.is.one.ind, ], 2, mean)
  ## visualization for covariates distributions in RCT population, and RWE
  # par(mfrow=c(p+1, 3))
  # for (hi in 0:p) {
  #   hist(x.in.rct[,hi+1])
  #   hist(x.in.rwe[,hi+1])
  #   hist(x.in.rct[,hi+1]*w.tilde)
  #   print(c(mean(x.in.rct[,hi+1]), mean(x.in.rwe[,hi+1]), mean(x.in.rct[,hi+1]*w.tilde)))
  # }
  
  # the second method to reweight; refer to "method(2)" in the draft
  #print(cost2(alphas, x.in.rct, x.in.rwe, N))
  opt2 <- optim(rep(0,p+1), cost2, x.in.rct=x.in.rct, x.in.rwe=x.in.rwe, N=N, method = "BFGS", control = list(maxit=10000))
  (alphas.hat2 <- opt2$par)
  w.tilde <- 1 + 1/(exp(x.in.rct %*% as.matrix(alphas.hat2, length(alphas.hat2), 1)))
  w.tilde2 <- w.tilde
  x.tilde <- x.in.rct * c(sqrt(w.tilde))
  # beta.hat.star2 <- optim(rep(2,p+1), cost3, y.hat.in.rct=y.star[delta.tilde.is.one.ind],
  # x.in.rct=x.in.rct, w.tilde=w.tilde, control = list(maxit=10000))$par
  # print(beta.hat.star2)
  ## Already use loop to checked "solve(t(x.tilde) %*% x.tilde)" is correct
  (beta.hat.star2 <- solve(t(x.tilde) %*% x.tilde) %*% as.matrix(apply(x.in.rct * y.star.hat[delta.tilde.is.one.ind]*c(w.tilde), 2, sum), p+1, 1))
  action2 <- (sign(x%*%beta.hat.star2)+1)/2
  (quality2.linear <- quality(action2, opt.actions, opt.value, baseline, contrast))
  
  # refer to "method(3)" in the draft
  # cost(alphas, x.in.rct, left.hand)
  opt <- optim(rep(0, p+1), cost, x.in.rct=x.in.rct, left.hand=left.hand,N=N, method = "BFGS", control = list(maxit=10000))
  (alphas.hat3 <- opt$par)
  w.tilde <- 1 + 1/(exp(x.in.rct %*% as.matrix(alphas.hat3, length(alphas.hat3), 1)))
  w.tilde3 <- w.tilde
  x.tilde <- x.in.rct * c(sqrt(w.tilde))
  # beta.hat.star3 <- optim(rep(2,p+1), cost3, y.hat.in.rct=y.star[delta.tilde.is.one.ind],
  #                         x.in.rct=x.in.rct, w.tilde=w.tilde)$par
  # print(beta.hat.star3)
  (beta.hat.star3 <- solve(t(x.tilde) %*% x.tilde) %*% as.matrix(apply(x.in.rct * y.star.hat[delta.tilde.is.one.ind]*c(w.tilde), 2, sum), p+1, 1))
  action3 <- (sign(x%*%beta.hat.star3)+1)/2
  (quality3.linear <- quality(action3, opt.actions, opt.value, baseline, contrast))
  
  ## method 4: cov balancing (refer to: https://github.com/yixinwang/minimal-weights-public )
  w.tilde4 <- weight_tuned_entr(tols=c(0.000, 0.001, 0.002, 0.005, 0.010, 0.020, 0.050, 0.100, 0.200, 0.500, 1.000), 
                                x.in.rct, x.in.rwe)
  x.tilde <- x.in.rct * c(sqrt(w.tilde4))
  (beta.hat.star4 <- solve(t(x.tilde) %*% x.tilde) %*% as.matrix(apply(x.in.rct * y.star.hat[delta.tilde.is.one.ind]*c(w.tilde4), 2, sum), p+1, 1))
  action4 <- (sign(x%*%beta.hat.star4)+1)/2
  (quality4.linear <- quality(action4, opt.actions, opt.value, baseline, contrast))
  
  ## method 5: kernel-based cov balancing (refer to: https://github.com/raymondkww/ATE.ncb)
  # w.tilde5 <- kern_cov_balance(x.in.rct, x.in.rwe)
  # x.tilde <- x.in.rct * c(sqrt(w.tilde5))
  # (beta.hat.star5 <- solve(t(x.tilde) %*% x.tilde) %*% as.matrix(apply(x.in.rct * y.star.hat[delta.tilde.is.one.ind]*c(w.tilde5), 2, sum), n, 1))
  # action5 <- (sign(x%*%beta.hat.star5)+1)/2
  # (quality5.linear <- quality(action5, opt.actions, opt.value, baseline, contrast))
  
  # unweighted
  w.tilde = rep(1, n)
  x.tilde <- x.in.rct * c(sqrt(w.tilde))
  # beta.hat.star <- optim(rep(2,p+1), cost3, y.hat.in.rct=y.star[delta.tilde.is.one.ind],
  #                         x.in.rct=x.in.rct, w.tilde=w.tilde, control = list(maxit=10000))$par
  # print(beta.hat.star)
  (beta.hat.star <- solve(t(x.tilde) %*% x.tilde) %*% as.matrix(apply(x.in.rct * y.star.hat[delta.tilde.is.one.ind]*c(w.tilde), 2, sum), p+1, 1))
  action <- (sign(x%*%beta.hat.star)+1)/2
  (quality.linear <- quality(action, opt.actions, opt.value, baseline, contrast))
  
  # use true alphas to estimate optimal linear policy
  w.tilde.true <- 1 + 1/(exp(x.in.rct %*% as.matrix(alphas, length(alphas), 1)))
  x.tilde <- x.in.rct * c(sqrt(w.tilde.true))
  (beta.hat.star.true.alpha <- solve(t(x.tilde) %*% x.tilde) %*% as.matrix(apply(x.in.rct * y.star.hat[delta.tilde.is.one.ind]*c(w.tilde.true), 2, sum), p+1, 1))
  action <- (sign(x%*%beta.hat.star.true.alpha)+1)/2
  (quality.linear.true.alpha <- quality(action, opt.actions, opt.value, baseline, contrast))
  ## results of absolute errors of beta.hat and alpha.hat
  beta_diff2 <- mean(abs(beta.hat.star2-beta.N.star))
  alpha_diff2 <- mean(abs(alphas.hat2-alphas))
  beta_diff3 <- mean(abs(beta.hat.star3-beta.N.star))
  alpha_diff3 <- mean(abs(alphas.hat3-alphas))
  beta_diff <- mean(abs(beta.hat.star-beta.N.star))
  beta_diff4 <- mean(abs(beta.hat.star4-beta.N.star))
  beta_diff_true.alpha <- mean(abs(beta.hat.star.true.alpha-beta.N.star))
  # beta_diff5 <- mean(abs(beta.hat.star5-beta.N.star))
  
  # the approximate optimal linear policy
  approx.opt.action <- (sign(x %*% beta.N.star) + 1)/2
  (quality.approx.opt.linear <- quality(approx.opt.action, opt.actions, opt.value, baseline, contrast))
  print(c(quality2.linear, quality3.linear, quality4.linear, quality.linear, quality.linear.true.alpha, quality.approx.opt.linear))
  
  
  ## olearning
  y.in.rct <- y[delta.tilde.is.one.ind] # observed outcomes in RCT
  # obtain the weights and labels needed in o-learning
  # R.hat <- randomForest(x.in.rct[,2:(p+1)], y.in.rct)$predicted
  R.hat <- ranger(y.in.rct ~ ., data=data.frame(cbind(y.in.rct, x.in.rct[,2:(p+1)])))$predictions
  R.tilde <- y.in.rct-R.hat
  label <- factor(((2*a[delta.tilde.is.one.ind]-1)*sign(R.tilde)+1)/2)
  data.classification <- data.frame(label, x.in.rct[,2:(p+1)], y.in.rct = y.in.rct,
                                    A=as.factor(a[delta.tilde.is.one.ind])) 
  ## weighted method (2)
  rule <- rpart(label ~ ., data=data.classification[,1:(p+1)], method = "class",
                weights=w.tilde2*abs(R.tilde/prob.one))
  action2 <- predict(rule, newdata=data.frame(x[, 2:(p+1)]), type="class")
  (quality2.tr <- quality(action2, opt.actions, opt.value, baseline, contrast))  
  ## weighted method (3)
  rule <- rpart(label ~ ., data=data.classification[,1:(p+1)], method = "class",
                weights=w.tilde3*abs(R.tilde/prob.one))
  action3 <- predict(rule, newdata=data.frame(x[, 2:(p+1)]), type="class")
  (quality3.tr <- quality(action3, opt.actions, opt.value, baseline, contrast))
  ## weighted method (4): cov balancing
  rule <- rpart(label ~ ., data=data.classification[,1:(p+1)], method = "class",
                weights=w.tilde4*abs(R.tilde/prob.one))
  action4 <- predict(rule, newdata=data.frame(x[, 2:(p+1)]), type="class")
  (quality4.tr <- quality(action4, opt.actions, opt.value, baseline, contrast))  
  ## weighted method (5): kernel-based cov balancing
  # rule <- rpart(label ~ ., data=data.classification[,1:(p+1)], method = "class",
  #               weights=w.tilde5*abs(R.tilde/prob.one))
  # action5 <- predict(rule, newdata=data.frame(x[, 2:(p+1)]), type="class")
  # (quality5.tr <- quality(action5, opt.actions, opt.value, baseline, contrast))  
    
  ## unweighted
  rule <- rpart(label ~ ., data=data.classification[,1:(p+1)], method = "class",
                weights = abs(R.tilde/prob.one))
  action <- predict(rule, newdata=data.frame(x[, 2:(p+1)]), type="class")
  (quality.tr <- quality(action, opt.actions, opt.value, baseline, contrast))
  ## use true alphas
  w.tilde <- 1 + 1/(exp(x.in.rct %*% as.matrix(alphas, length(alphas), 1)))
  rule <-  rpart(label ~ ., data=data.classification[,1:(p+1)], method = "class",
                 weights =w.tilde * abs(R.tilde/prob.one))
  action <- predict(rule, newdata=data.frame(x[, 2:(p+1)]), type="class")
  (quality.true.alphas.tr <- quality(action, opt.actions, opt.value, baseline, contrast))
  return(list(evaluation_linear=c(quality.linear, quality.linear.true.alpha, quality2.linear, quality3.linear, quality4.linear, quality.approx.opt.linear,
                                  beta_diff,beta_diff_true.alpha,beta_diff2,beta_diff3, beta_diff4, 
                                  alpha_diff2, alpha_diff3), 
              evaluation_tree=c(quality.tr, quality2.tr, quality3.tr, quality4.tr,quality.true.alphas.tr), 
              n=n, mean.opt.rate=mean(opt.actions), delta.tilde.is.one.ind=delta.tilde.is.one.ind, 
              w.tilde2=w.tilde2, w.tilde3=w.tilde3, w.tilde4=w.tilde4,
              mu.y.error=c(abs(mean(y[delta.tilde.is.one.ind])-mean(y)),
                           abs(sum(y[delta.tilde.is.one.ind]*w.tilde.true)/sum(w.tilde.true)-mean(y)),
                           abs(sum(y[delta.tilde.is.one.ind]*w.tilde2)/sum(w.tilde2)-mean(y)),
                           abs(sum(y[delta.tilde.is.one.ind]*w.tilde3)/sum(w.tilde3)-mean(y)),
                           abs(sum(y[delta.tilde.is.one.ind]*w.tilde4)/sum(w.tilde4)-mean(y)) )))
}



par(mfrow=c(5,2))
mus <- c(-1, -0.5, 0, 0.5, 1)
mus <- c(1)
Results <- list()
system.time(for (alpha1 in c(-5))   {
  for (m in c(5000)) {
    for (num.mu in 1:length(mus)) {
      mu <- mus[num.mu]
      alphas[1] <- alpha1
      r <- 300; 
      rs <- 1:r
      registerDoParallel(96) 
      results <- foreach(seed=rs, .combine=rbind, .errorhandling = 'remove',
                         .packages=c('mvtnorm','ranger', 'rpart', 'randomForest','CVXR')) %dopar% {
                           linear_tree(m, alphas, seed = seed,N=1000000, mu=mu, specify = F)
                         }
      stopImplicitCluster()
      print(paste0("mu=",mu,", alpha1=",alpha1, ", num of repeats: ", dim(results)[1]))
      print(round(apply(Reduce(rbind, results[,1]), 2, mean), 3))
      print(round(apply(Reduce(rbind, results[,2]), 2, mean), 3))
      print(round(apply(Reduce(rbind, results[,9]), 2, mean), 3))
      print(c(round(mean(Reduce(rbind, results[,3]))), round(mean(Reduce(rbind, results[,4])), 3)))
      Results[paste0("alpha",alpha1)] <- list(results)
      # print(rbind(round(apply(results$evaluation_linear, 2, mean), 2),
      # round(apply(results$evaluation_linear, 2, sd)/sqrt(length(rs)), 2)))
      # print(rbind(round(apply(results$evaluation_tree, 2, mean), 2),
      # round(apply(results$evaluation_tree, 2, sd)/sqrt(length(rs)), 2)))
      # data <- data.frame(value=c(results[,5], results[,6],results[,7] ,results[,8]),
      #                    accuracy=c(results[,1], results[,2], results[,3],results[,4]),
      #                    method <- c(rep("Unweighted", r), rep("Weighted(1)", r), rep("Weighted(2)", r),rep("Optimal Linear", r)))
      # name1 <- paste0("S:\\Documents\\lab2018\\transfer learning of dtr\\Alearning_alpha1_",abs(alpha1), "_mu", num.mu)
      # save(results, file=paste0(name1, ".RData"))
      # boxplot(value~method, data=data, col="lightgray",
      #         main=paste0("Value MSE"," (m=",m, ", n=",round(mean(results[,9])) ,", mu=", mu,  ")"))
      # boxplot(accuracy~method, data=data, col="lightgray",
      #         main=paste0("Optimal action accuracy"," (m=",m, ", n=",round(mean(results[,9])) ,", mu=", mu, ")"))
    }
  }
})
save(Results, file="results_est_contrast_rf.Rdata")

# ## Correctly specify
# >       print(paste0("m=",m,", alpha1=",alpha1))
# [1] "m=8000, alpha1=-10"
# >       print(round(apply(Reduce(rbind, results[,1]), 2, mean), 3))
# [1]  0.659  1.393  0.607  2.076  0.565  2.644  0.625  1.690  0.934  0.030  1.594  4.577 10.655  2.807  0.181  5.040
# >       print(round(apply(Reduce(rbind, results[,2]), 2, mean), 3))
# [1] 0.747 0.829 0.709 1.158 0.656 1.647 0.730 0.937 0.711 1.142
# >       print(c(round(mean(Reduce(rbind, results[,3]))), round(mean(Reduce(rbind, results[,4])), 3)))
# [1] 318.000   0.659
# 
# 
# 
# ## Not correctly specify
# [1] "m=8000, alpha1=-11"
# [1] 0.743 0.743 0.665 1.604 0.645 1.857 0.730 0.906 0.934 0.030 1.375 2.959 3.302 1.723 1.086 0.989
# [1] 0.788 0.619 0.758 0.827 0.752 0.898 0.763 0.847 0.703 1.253
# [1] 369.000   0.659
# [1] "m=8000, alpha1=-10"
# [1] 0.782 0.504 0.714 1.093 0.709 1.190 0.765 0.625 0.934 0.030 0.923 2.073 2.104 1.263 1.045 1.021
# [1] 0.831 0.360 0.806 0.516 0.807 0.533 0.828 0.411 0.733 1.128
# [1] 844.000   0.659
# [1] "m=8000, alpha1=-8"
# [1] 0.858 0.166 0.801 0.442 0.799 0.474 0.848 0.198 0.934 0.030 0.408 1.066 1.104 0.558 0.950 0.959
# [1] 0.911 0.080 0.874 0.183 0.876 0.196 0.898 0.102 0.765 0.852
# [1] 4311.000    0.659

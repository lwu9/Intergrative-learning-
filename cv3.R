#################################################
### add the weighted estimators besides value ###
#################################################
#source('S:\\Documents\\lab2018\\transfer learning of dtr\\functions.R')
source("/home/lwu9/transfer/functions2.R")
# path <- "/Users/lwu9/Documents/transfer_learner_dtr"; setwd(path); source(paste0(path,'/functions2.R'))

cv.compare <- function(m, alphas, seed, N, mu, p, method, Qm, y.hat.meth, nfold=30, specify=F, misspecify, dbn, cv=T, weighted.para) {
  ## y.hat.meth: specify which method to est contrast
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
    if (dbn == "diff") {
      contrast <- (cos(x)-0.5) %*% betas2 ## dist are different in rwe and rct
    } else {
      contrast <- atan(exp(1+x[,2])-3*x[,3]-5) ## dist are similar in rwe and rct
    }
    # contrast <- betas2[1]*(cos(x[,1])-0.5) + betas2[2]*(cos(x[,2])-0.5) +betas2[3]*(cos(x[,3])-0.5)
    # contrast <- 3 * ((x[,2]<=1) & (x[,3]>-0.6)) - 1
  }
  y1.star <-  baseline + contrast + rnorm(N, sd=sigma)
  y0.star <-  baseline + rnorm(N, sd=sigma)
  opt.actions <- (sign(contrast)+1)/2
  opt.value <-  baseline + opt.actions*(contrast)
  y.star <- y1.star - y0.star
  y <-  rep(0, N); y[which(a==1)] <- y1.star[which(a==1)]; y[which(a==0)] <- y0.star[which(a==0)] # observed outcomes in RCT
  # ############### optimal linear rule in the target population ##################
  # w.tilde <- learn_rule(y, x, a, x.in.rwe=NA, x=NA, N=NA, y.star.hat.in.rct=NA, 
  #                       w.method=1, method=NA, prob.one=NA, only.weight=T, misspecify = NA)
  # y.star.hat <- est.contrast(a, x, y, meth=y.hat.meth, psm=1, Qm=Qm, 
  #                            regr.direct=F, weighted.para=F, w.tilde=w.tilde)
  # opt.actions.linear <- learn_rule_given_w.tilde(y, x, a, y.star, x, method="linear.VS", prob.one = NA, w.tilde) 
  # opt.value.linear <- baseline + opt.actions.linear*(contrast)
  # print(quality(opt.actions.linear, opt.actions, opt.value, baseline, contrast))
  ######## generate rwe #########
  delta.is.one.ind <- sample(N, size=m)
  x.in.rwe <- x[delta.is.one.ind, ]
  ######## generate rct ########
  delta.tilde <- rbinom(N, 1, delta.tilde.prob(x, alphas)) ## correctly specify
  # delta.tilde <- rbinom(N, 1, delta.tilde.prob(cos(x)+0.5, alphas)) ## misspecify
  (n <- sum(delta.tilde)) # sample size in rct
  # print(paste(mean(opt.actions), n))
  delta.tilde.is.one.ind <- which(delta.tilde==1)
  x.in.rct <- x[delta.tilde.is.one.ind, ]
  y.in.rct <- y[delta.tilde.is.one.ind]
  a.in.rct <- a[delta.tilde.is.one.ind]
  
  ## visualization for covariates distributions in RCT population, and RWE
  # par(mfrow=c(p+1, 3))
  # for (hi in 0:p) {
  #   hist(x.in.rct[,hi+1])
  #   hist(x.in.rwe[,hi+1])
  #   hist(x.in.rct[,hi+1]*w.tilde)
  #   print(c(mean(x.in.rct[,hi+1]), mean(x.in.rwe[,hi+1]), mean(x.in.rct[,hi+1]*w.tilde)))
  # }
  true.values.mse.rule <- true.values.mse.linear <- c()
  for (w.method in c(4, 3, 2, 1)) {
    ####################################
    w.tilde <- learn_rule(y.in.rct, x.in.rct, a.in.rct, x.in.rwe, x, N, y.star.hat.in.rct=NA, 
                          w.method, method, prob.one=NA, only.weight=T, misspecify = misspecify)
    ####################################
    y.star.hat.in.rct <- est.contrast(a.in.rct, x.in.rct, y.in.rct, meth=y.hat.meth, psm=1, Qm=Qm, 
                                      regr.direct=F, weighted.para=weighted.para, w.tilde=w.tilde)
    ## Calculate contrast estimator MSE in RCT
    (ctr.mse <- mean((y.star.hat.in.rct-contrast[delta.tilde.is.one.ind])**2)) 
    print(paste(w.method, ctr.mse))
    action <- learn_rule_given_w.tilde(y.in.rct, x.in.rct, a.in.rct, y.star.hat.in.rct, 
                                       x, method, prob.one = NA, w.tilde) 
    (evaluation <- quality(action, opt.actions, opt.value, baseline, contrast))
    # print(paste(method, w.method, ": ", round(evaluation[1],3), round(evaluation[2],3)))
    true.values.mse.rule <- c(true.values.mse.rule, evaluation[2])
    # (evaluation.linear <- quality(action, opt.actions.linear, opt.value.linear, baseline, contrast))
    # true.values.mse.linear <- c(true.values.mse.linear, evaluation.linear[2])
    # print(true.values.mse.rule)
    # print(true.values.mse.linear)
    # print(paste("In the above, weighted.para is", weighted.para, "; n is", n, "#############"))
  }
  if (cv) {
    #### need to fix the cv code below ###
    (cv.results <- cv(y.in.rct, x.in.rct, a.in.rct, x.in.rwe, N, nfold=nfold, y.star.hat.in.rct,
                      prob.one, ipw=F, method=method, baseline, contrast, x, calculate.true.val=F, misspecify=misspecify))
    return(list(true.values.mse.rule=true.values.mse.rule, cv.value.rule=cv.results$values, 
                cv.chosen.rule.value.mse=true.values.mse.rule[cv.results$num], n=n, 
                ctr.mse=ctr.mse))
  } else {
    return(list(true.values.mse.rule=true.values.mse.rule, n=n, ctr.mse=ctr.mse))
  }

}

##########################################################################################

N=1000000
p <- 2; 

sigma=0.5; 
# alphas <- as.matrix(c(-8, 1), rep(0, p-1), p+1, 1)
betas1 <- betas2 <- betas <- as.matrix(1:(p+1), p+1, 1)
# alphas <- as.matrix(c(-10, 1, -2, 1, 1,rep(0, p-6)), p+1, 1)
alphas <- as.matrix(c(-10, 1,-2, rep(0, p-2)), p+1, 1)
# betas <- betas1 <- betas2 <- as.matrix(c(-1/2,1), p+1, 1)
# betas1 <- c(2, 1, 0, 1, 0, rep(0, p-4)); betas2 <- c(0, 1, -1, 1, -1, rep(0, p-4))
seed=1; mu=1
mus <- c(1)
nfold = 5; methods=c("linear.VS")
misspecify = F
dbn <- "similar"
cv <- FALSE
for (method in methods) {
  for (alpha1 in c(-10, -8)) {
    Results <- list()
    if (alpha1==-10) rs <- 1:(56*5) else
      rs <- 1:(56*7)
    for (m in c(5000)) {
      for (num.mu in 1:length(mus)) {
        mu <- mus[num.mu]
        alphas[1] <- alpha1
        for (y.hat.meth in 1:3) {
          if (y.hat.meth==1) Qms <- c(0) else
            Qms <- 1:2
          for (Qm in Qms) {
            for (weighted.para in c(T, F)) {
              registerDoParallel(2) 
              results <- foreach(seed=rs, .combine=rbind, .errorhandling=c('remove'),
                                 .packages=c('mvtnorm','ranger', 'rpart', 'randomForest','CVXR','rgenoud')) %dopar% {
                                   cv.compare(m, alphas, seed, N, mu, p, method, Qm, y.hat.meth, nfold=30, specify=F, 
                                              misspecify, dbn, cv=cv, weighted.para)
                                 }
              stopImplicitCluster()
              # print(head(Reduce(rbind, results[,1])))
              # print(learned.rule <- apply(Reduce(rbind, results[,1]), 1, which.min))
              # orders <- apply(Reduce(rbind, results[,1]), 1, order)
              # print(cv.rule <- apply(Reduce(rbind, results[,2]), 1, which.max))
              learned.rule.mse <- matrix(unlist(results[,1]), byrow = T, ncol=4)
              rep <- dim(learned.rule.mse)[1]
              # cv.rule.value.est <- matrix(unlist(results[,2]), byrow = T, ncol=4)
              # orders <- apply(learned.rule.mse, 1, order)
              # cv.rule <- apply(cv.rule.value.est, 1, which.max)
              print(paste(mean(unlist(results[,2])), "; method", method, "; y.hat.meth is", y.hat.meth, "; Qm=",
                          Qm, "; alpha=", alpha1, "; rep=", rep, "; weighted.para", weighted.para, "###########"))
              # print(dim(cv.rule.value.est))
              # print(paste("################ correct rate using cv:", mean(cv.rule==orders[1,])))
              # print(paste("################ correct rate using cv 2nd:", mean(cv.rule==orders[2,])))
              # print(paste("################ correct rate using cv 3rd:", mean(cv.rule==orders[3,])))
              # print(paste("################ correct rate using cv 4th:", mean(cv.rule==orders[4,])))
              Results[paste0("alpha",alpha1,"yhatmeth",y.hat.meth,"Qm",Qm)] <- list(results)
              print("MSE of w_np, w2, w1, unweight and its SE:")
              print(apply(learned.rule.mse, 2, mean))
              print(apply(learned.rule.mse, 2, sd)) / sqrt(rep)
              # print(paste("CV rule MSE:",mean(unlist(results[,3])) ))   
              print(paste("contrast estimator MSE:",mean(unlist(results[,3])) )) 
            }
          }
        }
      }
    }
    save(Results, file=paste0("results_cv_nfold", nfold, "_method",method,"_alpha",abs(alpha1),
                              "_diff_dist_all_except_NN_mis", misspecify, "_CV", cv, ".Rdata"))
  }
}


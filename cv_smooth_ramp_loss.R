source("/home/lwu9/transfer/functions3.R")
# path <- "/Users/lwu9/Documents/transfer_learner_dtr"; setwd(path); source(paste0(path,'/functions3.R'))

cv.compare.sramp <- function(m, alphas, seed, N, mu, p, method, Qm, y.hat.meth, nfold, specify=F, 
                       misspecify, dbn, cv.flag, weighted.para=T) {
  ## y.hat.meth: specify which method to est contrast
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
  y <-  rep(0, N); y[which(a==1)] <- y1.star[which(a==1)]; y[which(a==0)] <- y0.star[which(a==0)] # observed outcomes in RCT
  # ############### optimal linear rule in the target population 
  # w.tilde <- learn_rule(y, x, a, x.in.rwe=NA, x=NA, N=NA, y.star.hat.in.rct=NA, 
  #                       w.method=1, method=NA, prob.one=NA, only.weight=T, misspecify = NA)
  # y.star.hat <- est.contrast(a, x, y, meth=y.hat.meth, psm=1, Qm=Qm, regr.direct=F, weighted.para=F, w.tilde=w.tilde)
  # opt.actions.linear <- learn_rule_given_w.tilde(y, x, a, y.star, x, method="linear.VS", prob.one = NA, w.tilde) 
  # opt.value.linear <- baseline + opt.actions.linear$action*(contrast)
  # print(quality(opt.actions.linear$action, opt.actions, opt.value, baseline, contrast)) ## measure how close the opt linear rule to the opt rule
  ######## Generate rwe 
  delta.is.one.ind <- sample(N, size=m)
  x.in.rwe <- x[delta.is.one.ind, ]
  ######## Generate rct 
  delta.tilde <- rbinom(N, 1, delta.tilde.prob(x, alphas)) 
  (n <- sum(delta.tilde)) # sample size in rct
  # print(paste(mean(opt.actions), n))
  delta.tilde.is.one.ind <- which(delta.tilde==1)
  x.in.rct <- x[delta.tilde.is.one.ind, ]
  y.in.rct <- y[delta.tilde.is.one.ind]
  a.in.rct <- a[delta.tilde.is.one.ind]
  
  true.values.mse.rule <-  true.values.mse.rule.genound <- c()
  for (w.method in c(4, 3, 2, 1)) {
    w.tilde <- learn_rule(y.in.rct, x.in.rct, a.in.rct, x.in.rwe, x=NA, N, y.star.hat.in.rct=NA, 
                          w.method, method, prob.one=NA, only.weight=T, misspecify = misspecify)
    w.tilde <- w.tilde/sum(w.tilde) ## normailize the weights s.t. their sum equal to 1
    y.star.hat.in.rct <- est.contrast(a.in.rct, x.in.rct, y.in.rct, meth=y.hat.meth, psm=1, Qm=Qm, 
                                      regr.direct=F, weighted.para=weighted.para, w.tilde=w.tilde)
    
    beta.hat.star <- get_eta(x.in.rct, w.tilde, y.star.hat.in.rct)
    # cat( beta.hat.star/abs(beta.hat.star[1]), "\n")
    action <- as.numeric(x %*% beta.hat.star > 0)
    (evaluation <- quality(action, opt.actions, opt.value, baseline, contrast)) ## value mse to the opt rule
    true.values.mse.rule <- c(true.values.mse.rule, evaluation[2])
    
    # action <- learn_rule_given_w.tilde(y.in.rct, x.in.rct, a.in.rct, y.star.hat.in.rct,
    #           x, method, prob.one = NA, w.tilde, iter=10)$action
    # (evaluation <- quality(action, opt.actions, opt.value, baseline, contrast)) ## value mse to the opt rule
    # true.values.mse.rule.genound <- c(true.values.mse.rule.genound, evaluation[2])
  }
  if (cv.flag) {
    #### need to fix the cv code below ###
    (cv.results <- cv.new.sramp(y.in.rct, x.in.rct, a.in.rct, x.in.rwe, N, nfold=nfold, prob.one, y.hat.meth, Qm, ipw=F, method=method, 
                                baseline, contrast, x=NA, calculate.true.val=F, misspecify=misspecify, weighted.para=weighted.para))
    return(list(true.values.mse.rule=true.values.mse.rule, n=n, cv.value.rule=cv.results$values, 
                cv.chosen.rule.value.mse=true.values.mse.rule[cv.results$num]))
  } else {
    return(list(true.values.mse.rule=true.values.mse.rule, n=n))#,  
    # true.values.mse.rule.genound=true.values.mse.rule.genound))
  }
}

##########################################################################################

N=1000000
p <- 2

sigma=0.5
# alphas <- as.matrix(c(-8, 1), rep(0, p-1), p+1, 1)
betas1 <- betas2 <- betas <- as.matrix(1:(p+1), p+1, 1)
# alphas <- as.matrix(c(-10, 1, -2, 1, 1,rep(0, p-6)), p+1, 1)
alphas <- as.matrix(c(-10, 1,-2, rep(0, p-2)), p+1, 1)
# betas <- betas1 <- betas2 <- as.matrix(c(-1/2,1), p+1, 1)
# betas1 <- c(2, 1, 0, 1, 0, rep(0, p-4)); betas2 <- c(0, 1, -1, 1, -1, rep(0, p-4))
# seed=1; mu=1
mus <- c(1)
nfold = 10; methods=c("linear.VS")
misspecify = F
# dbn <- "similar"
cv.flag <- T
rs <- 1:200
for (dbn in c("similar", "diff")) {
  for (method in methods) {
    for (alpha1 in c(-10, -8)) {
      Results <- list()
      for (m in c(5000)) {
        for (num.mu in 1:length(mus)) {
          mu <- mus[num.mu]
          alphas[1] <- alpha1
          for (y.hat.meth in 1:3) {
            if (y.hat.meth==1) Qms <- c(0) else
              Qms <- 1:2
            for (Qm in Qms) {
              for (weighted.para in c(T)) {
                registerDoParallel(67) 
                results <- foreach(seed=rs, .combine=rbind, .errorhandling=c('remove'),
                                   .packages=c('mvtnorm','ranger', 'rpart', 'randomForest','CVXR','rgenoud')) %dopar% {
                                     cv.compare.sramp(m, alphas, seed, N, mu, p, method, Qm, y.hat.meth, nfold=nfold, specify=F, 
                                                misspecify, dbn, cv.flag, weighted.para)
                                   }
                stopImplicitCluster()
                learned.rule.mse <- matrix(unlist(results[,1]), byrow = T, ncol=4)
                rep <- dim(learned.rule.mse)[1]
                print(paste0("n=", round(mean(unlist(results[, 2]))), "; y.hat.meth=", y.hat.meth, "; Qm=", Qm, "; alpha=", alpha1, "; rep=", 
                             rep, "; dbn is ", dbn, "; misspecify=", misspecify, " ###########"))
                Results[paste0("alpha",alpha1,"yhatmeth",y.hat.meth,"Qm",Qm, "weighted.para", weighted.para)] <- list(results)
                if (cv.flag) {
                  print("MSE of w_np, w2, w1, unweight, cv, and its SE:")
                  print(c(apply(learned.rule.mse, 2, mean), mean(unlist(results[, 4]))))
                  print(c(apply(learned.rule.mse, 2, sd), sd(unlist(results[, 4]))) / sqrt(rep))
                } else {
                  print("MSE of w_np, w2, w1, unweight and its SE:")
                  print(apply(learned.rule.mse, 2, mean))
                  print(apply(learned.rule.mse, 2, sd) / sqrt(rep))
                }
              }
            }
          }
        }
      }
      save(Results, file=paste0("sramp_results_cv_nfold", nfold, "_method", method, "_alpha", abs(alpha1),
                                "_", dbn, "dist_mis", misspecify, "_CV", cv.flag, ".Rdata"))
    }
  }
}


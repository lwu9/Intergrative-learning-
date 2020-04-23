source("/home/lwu9/transfer/functions3.R")
list.of.packages <- c("FindIt", "BART")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)
library(FindIt)
library(BART)

imai.compare <- function(m, alphas, seed, N, mu, p, specify=F, dbn) {
  ## Results from Imai paper, as a comparison
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
  # y.star.hat <- est.contrast(a, x, y, meth=y.hat.meth, psm=1, Qm=Qm, regr.direct=F, weighted.para=F, w.tilde=w.tilde)
  # opt.actions.linear <- learn_rule_given_w.tilde(y, x, a, y.star, x, method="linear.VS", prob.one = NA, w.tilde) 
  # opt.value.linear <- baseline + opt.actions.linear$action*(contrast)
  # print(quality(opt.actions.linear$action, opt.actions, opt.value, baseline, contrast)) ## measure how close the opt linear rule to the opt rule
  ######## Generate rwe #########
  delta.is.one.ind <- sample(N, size=m)
  x.in.rwe <- x[delta.is.one.ind, ]
  ######## Generate rct ########
  delta.tilde <- rbinom(N, 1, delta.tilde.prob(x, alphas)) 
  (n <- sum(delta.tilde)) # sample size in rct
  # print(paste(mean(opt.actions), n))
  delta.tilde.is.one.ind <- which(delta.tilde==1)
  x.in.rct <- x[delta.tilde.is.one.ind, ]
  y.in.rct <- y[delta.tilde.is.one.ind]
  a.in.rct <- a[delta.tilde.is.one.ind]
  
  train <- rbind(x.in.rct[,-1], x.in.rwe[,-1])
  y.train <- c(rep(1, n), rep(0, m))
  post <- pbart(train, y.train)
  post.prob.mean <- apply(pnorm(post$yhat.train), 2, mean)
  prob <- post.prob.mean[1:n]
  wts <- 1/prob / mean(1/prob)
  dat <- data.frame(y=y.in.rct, x=x.in.rct[,-1], a=a.in.rct)
  dat$wts <- wts
  F1 <-FindIt(model.treat= y ~ a,
              model.main= ~ x.1+x.2,
              model.int= ~ x.1+x.2,
              data = dat, wts=wts, 
              type="continuous",
              treat.type="single")
  pred <- predict(F1, newdata=data.frame(x=x[,-1], a=a))
  trt.effect <- pred$data[order(as.numeric(rownames(pred$data))), 1]
  action <- as.numeric(trt.effect > 0)
  (evaluation <- quality(action, opt.actions, opt.value, baseline, contrast)) 
  return(list(value.mse=evaluation[2], n=n))
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
methods=c("linear.VS")
for (dbn in c("similar", "diff")) {
  for (method in methods) {
    for (alpha1 in c(-10, -8)) {
      Results <- list()
      rs <- 1:200
      for (m in c(5000)) {
        for (num.mu in 1:length(mus)) {
          mu <- mus[num.mu]
          alphas[1] <- alpha1
          registerDoParallel(50)
          results <- foreach(seed=rs, .combine=rbind, .errorhandling=c('remove'),
                             .packages=c('FindIt', 'BART')) %dopar% {
                               imai.compare(m, alphas, seed, N, mu, p, specify=F, dbn)
                             }
          stopImplicitCluster()
                # print(head(Reduce(rbind, results[,1])))
                # print(learned.rule <- apply(Reduce(rbind, results[,1]), 1, which.min))
                # orders <- apply(Reduce(rbind, results[,1]), 1, order)
                # print(cv.rule <- apply(Reduce(rbind, results[,2]), 1, which.max))
          learned.rule.mse <- unlist(results[,1])
          rep <- length(learned.rule.mse)
          # cv.rule.value.est <- matrix(unlist(results[,2]), byrow = T, ncol=4)
          # orders <- apply(learned.rule.mse, 1, order)
          # cv.rule <- apply(cv.rule.value.est, 1, which.max)
          print(paste0("n=", round(mean(unlist(results[, 2]))), "; MES=", round(mean(learned.rule.mse), 5), 
                       "; alpha=", alpha1, "; dbn=", dbn, "; rep=", rep, " ###########"))
                # print(dim(cv.rule.value.est))
                # print(paste("################ correct rate using cv:", mean(cv.rule==orders[1,])))
                # print(paste("################ correct rate using cv 2nd:", mean(cv.rule==orders[2,])))
                # print(paste("################ correct rate using cv 3rd:", mean(cv.rule==orders[3,])))
                # print(paste("################ correct rate using cv 4th:", mean(cv.rule==orders[4,])))
          Results[paste0("alpha",alpha1, "_dbn_", dbn)] <- list(results)
        }
        save(Results, file=paste0("results", "_method", method, "_alpha", abs(alpha1),
                                  "_", dbn, "dist_mis", misspecify, "_imai.Rdata"))
      }
    }
  }
}




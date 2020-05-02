# source("/home/lwu9/transfer/functions.R")
######## FindIt
rm(list=ls())
set.seed(1)
library(FindIt)
library(BART)
data("LaLonde")
devtools::install_github("jjchern/lalonde")

rwe <- lalonde::cps_controls
rct <- lalonde::nsw_dw
p <- 8

# rwe <- lalonde::psid_controls; rwe$re74 <- NULL
# rct <- lalonde::nsw
# p <- 7

m <- dim(rwe)[1]
n <- dim(rct)[1]
x.in.rwe <- matrix(NA, m, p)
x.in.rct <- matrix(NA, n, p)

for (j in 1:p) {
  x.in.rwe[,j] <- rwe[[j+2]]
  x.in.rct[,j] <- rct[[j+2]]
  if (j > 6) {
    x.in.rwe[,j] <- log(1+rwe[[j+2]])
    x.in.rct[,j] <- log(1+rct[[j+2]])
  }
}

x.in.rwe <- cbind(rep(1, m), x.in.rwe)
x.in.rct <- cbind(rep(1, n), x.in.rct)

a.in.rct <- rct$treat
y.in.rct <- log(1+rct$re78)
colMeans(x.in.rwe)
colMeans(x.in.rct)
dat <- data.frame(outcome=y.in.rct, treat=a.in.rct, age=x.in.rct[,2], educ=x.in.rct[,3], black=x.in.rct[,4], hisp=x.in.rct[,5], 
                  marr=x.in.rct[,6], nodegr=x.in.rct[,7], log.re74=x.in.rct[,8], log.re75=x.in.rct[,9])
train <- rbind(x.in.rct[,-1], x.in.rwe[,-1])
m <- dim(x.in.rwe)[1]
y.train <- c(rep(1, n), rep(0, m))
post <- pbart(train, y.train)
post.prob.mean <- apply(pnorm(post$yhat.train), 2, mean)
prob <- post.prob.mean[1:n]
wts <- 1/prob / mean(1/prob)
dat$wts <- wts

F1 <-FindIt(model.treat= outcome ~ treat,
            model.main= ~ age+educ+black+hisp+
              marr+nodegr+log.re74+log.re75,
            model.int= ~ age+educ+black+hisp+
              marr+nodegr+log.re74+log.re75,
            data = dat, wts=wts, 
            type="continuous",
            treat.type="single")
summary(F1)
## Returns all the estimated treatment effects.
pred1 <- predict(F1)
head(pred1$data, n=10)
tail(pred1$data, n=10)


library(optmatch)
path <- "/Users/lwu9/Documents/transfer_learner_dtr"; setwd(path); source(paste0(path,'/functions3.R'))

ladonde <- function(y.hat.meth, Qm, iter, cv.flag=T, nfold, reg.mu0.rf, reg.mu1.rf, w.tilde, 
                    a.in.rct, x.in.rct, y.in.rct, x.in.rwe) {
  results <- list()
  ### nonpara weighted
  y.star.hat.in.rct <- est.contrast(a.in.rct, x.in.rct, y.in.rct, meth=y.hat.meth, psm=1, Qm=Qm, 
                                    regr.direct=F, weighted.para=T, w.tilde=w.tilde)
  ####### smooth ramp loss
  weight.beta.hat.star <- get_eta(x.in.rct, w.tilde, y.star.hat.in.rct)
  action.weight <- as.numeric(x.in.rwe %*% weight.beta.hat.star > 0)
  ####### 0-1 loss
  # weighted.learn <- learn_rule_given_w.tilde(y.in.rct, x.in.rct, a.in.rct, y.star.hat.in.rct,
  #                                            x=x.in.rwe, method, prob.one = NA, w.tilde, iter=iter)
  # action.weight <- weighted.learn$action
  # weight.beta.hat.star <- weighted.learn$beta.hat.star
  
  (weight.beta.hat.star <- weight.beta.hat.star / abs(weight.beta.hat.star[1]))
  ###### unweighted
  n <- length(y.in.rct)
  y.star.hat.in.rct <- est.contrast(a.in.rct, x.in.rct, y.in.rct, meth=y.hat.meth, psm=1, Qm=Qm, 
                                    regr.direct=F, weighted.para=T, w.tilde=rep(1,n)/n)
  ####### smooth ramp loss
  unweight.beta.hat.star <- get_eta(x.in.rct, w.tilde=rep(1,n)/n, y.star.hat.in.rct)
  action.unweight <- as.numeric(x.in.rwe %*% unweight.beta.hat.star > 0)
  ####### 0-1 loss
  # unweighted.learn <- learn_rule_given_w.tilde(y.in.rct, x.in.rct, a.in.rct, y.star.hat.in.rct,
  #                                              x=x.in.rwe, method="linear.VS", prob.one = NA, rep(1,n)/n, iter=iter)
  # action.unweight <- unweighted.learn$action
  # unweight.beta.hat.star <- unweighted.learn$beta.hat.star
  
  (unweight.beta.hat.star <- unweight.beta.hat.star / abs(unweight.beta.hat.star[1]))
  #### estimated values ####
  weighted.value <- reg.mu0.rf; weighted.value[action.weight==1] <- reg.mu1.rf[action.weight==1]
  unweighted.value <- reg.mu0.rf; unweighted.value[action.unweight==1] <- reg.mu1.rf[action.unweight==1]
  (val.est.rf <- c(mean(weighted.value), mean(unweighted.value)))
  
  if (cv.flag) {
    ### consider cross validation
    (cv.results <- cv.real.data.new(y.in.rct, x.in.rct, a.in.rct, x.in.rwe, nfold=nfold, y.hat.meth=y.hat.meth, Qm=Qm))
    val.est.rf <- c(val.est.rf, val.est.rf[cv.results$num])
  }
  results[[paste0("y.hat.meth", y.hat.meth, "Qm", Qm)]] <- list(val.ests=val.est.rf,
                                                                beta.hat=rbind(weight.beta.hat.star, unweight.beta.hat.star))
  return(results)
}

test <- function(res, x.in.rwe, x.in.rct, a.in.rct, rwe, testm=100, Rep=100) {
  action.weight <- as.numeric(x.in.rwe %*% res[[1]]$beta.hat[1,] > 0)
  action.unweight <- as.numeric(x.in.rwe %*% res[[1]]$beta.hat[2,] > 0)
  m <- dim(rwe)[1]
  imai <- predict(F1, newdata=data.frame(treat=rep(0, m), age=x.in.rwe[,2], educ=x.in.rwe[,3], black=x.in.rwe[,4], hisp=x.in.rwe[,5], 
                                                marr=x.in.rwe[,6], nodegr=x.in.rwe[,7], log.re74=x.in.rwe[,8], log.re75=x.in.rwe[,9]))
  action.imai <- as.numeric(imai$data[,1] > 0)
  value.ests <- c()
  for (rep in 1:Rep) {
    set.seed(rep)
    m <- dim(x.in.rwe)[1]
    m.test <- sample(1:m, testm)
    x.in.rwe.test <- data.frame(x.in.rwe[m.test, ])
    names(x.in.rwe.test)[1] <- "rwe" ## 1: RWE samples; 0: RCT samples
    x.in.rct.test <- data.frame(x.in.rct); x.in.rct.test[,1] <- 0; names(x.in.rct.test)[1] <- "rwe"
    test.data <- data.frame(rbind(x.in.rwe.test, x.in.rct.test[a.in.rct==1, ])) ## only match the treated samples in RCT
    pm1 <- pairmatch(rwe ~ ., data = test.data)
    table(test.data[which(is.na(pm1)==F),]$rwe)
    y.in.rwe <- log(1 + rwe$re78)
    y.test <- c(y.in.rwe, y.in.rct[a.in.rct==1])[which(is.na(pm1)==F)]
    matched.data <- test.data[which(is.na(pm1)==F),]
    ## random forest
    loc.a1<-which(matched.data$rwe==0)
    loc.a0<-which(matched.data$rwe==1) ## RWE samples are in control samples
    dat <- data.frame(cbind(y.test, matched.data[,-1]))
    dat.rwe <- data.frame(x.in.rwe[,-1]); names(dat.rwe) <- names(dat)[-1]
    fit1 <- ranger(y.test ~ ., data=dat[loc.a1,])
    reg.mu1.rf <- predict(fit1, data = dat.rwe)$predictions
    fit0 <- ranger(y.test ~ ., data=dat[loc.a0,])
    reg.mu0.rf <- predict(fit0, data = dat.rwe)$predictions
    
    weighted.value <- reg.mu0.rf; weighted.value[action.weight==1] <- reg.mu1.rf[action.weight==1]
    unweighted.value <- reg.mu0.rf; unweighted.value[action.unweight==1] <- reg.mu1.rf[action.unweight==1]
    imai.value <- reg.mu0.rf; imai.value[action.imai==1] <- reg.mu1.rf[action.imai==1]
    val.est.rf <- c(mean(weighted.value), mean(unweighted.value), mean(imai.value))
    value.ests <- rbind(value.ests, val.est.rf)
  }
  # cat(apply(value.ests, 2, mean), apply(value.ests, 2, sd), "\n")
  return(list(means=apply(value.ests, 2, mean), ses=apply(value.ests, 2, sd)/sqrt(Rep)))
}

head(x.in.rct)
head(rct)
summary(lm(y.in.rct~., data=data.frame(cbind(y.in.rct, x.in.rct[,-1], a.in.rct*x.in.rct))))

method <- "linear.VS"
w.tilde <- learn_rule(y.in.rct, x.in.rct, a.in.rct, x.in.rwe, x=NA, N=NA, y.star.hat.in.rct=NA, 
                      w.method=4, method, prob.one=NA, only.weight=T, misspecify = F)

loc.a1<-which(a.in.rct==1)
loc.a0<-which(a.in.rct==0)
dat <- data.frame(cbind(y.in.rct, x.in.rct))
dat.rwe <- data.frame(x.in.rwe); names(dat.rwe) <- names(dat)[-1]
fit1 <- ranger(y.in.rct~., data=dat[loc.a1,], case.weights = w.tilde[loc.a1])
reg.mu1.rf <- predict(fit1, data = dat.rwe)$predictions
fit0 <- ranger(y.in.rct~., data=dat[loc.a0,], case.weights = w.tilde[loc.a0])
reg.mu0.rf <- predict(fit0, data = dat.rwe)$predictions
opt.value <- apply(rbind(reg.mu1.rf, reg.mu0.rf), 2, max)
print(mean(opt.value))

### Choose y.hat.meth = 3, Qm=2
res <- ladonde(y.hat.meth = 3, Qm = 2, iter=NA, cv.flag=T, nfold=10, reg.mu0.rf = reg.mu0.rf, 
               reg.mu1.rf = reg.mu1.rf, w.tilde=w.tilde, a.in.rct, x.in.rct, y.in.rct, x.in.rwe)
(eval <- test(res, x.in.rwe, x.in.rct, a.in.rct, rwe, testm=100))
print(res)
x.test.rwe.pos <- Reduce(cbind, pred1$data[1:10,4:11])
x.test.rwe.neg <- Reduce(cbind, tail(pred1$data, n=10)[, 4:11])
t(round(cbind(rep(1, 10),x.test.rwe.pos) %*% t(res[[1]]$beta.hat), 2))
t(round(cbind(rep(1, 10),x.test.rwe.neg)%*% t(res[[1]]$beta.hat), 2))

                                         
library(xtable)
xtable(rbind(colMeans(x.in.rct), colMeans(x.in.rwe)))
xtable(res[[1]][[2]])     
paste0(t(round(res[[1]][[2]], 2)), collapse = ",")

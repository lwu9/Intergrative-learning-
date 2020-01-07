############# Transfer learner projet: functions ############

list.of.packages <- c("ggplot2", "mvtnorm", "doParallel", "parallel", "rpart", "rpart.plot", "ranger", "MASS",
                      "randomForest", "wSVM", "DTRlearn2","CVXR","rgenoud")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)

# devtools::install_github("raymondkww/ATE.ncb")
# library(ATE.ncb)
library(mvtnorm)
library(foreach)
library(doParallel)
library(parallel)
library(rpart)
library(rpart.plot)
library(ranger)
library(MASS)
library(randomForest)
library(wSVM)
# library(DynTxRegime)
library(DTRlearn2)
library(CVXR)
library(rgenoud)

delta.tilde.prob <- function(x, alphas) {
  # x: covariate vecter or matrix including intercep 1
  # return: the probability of \widetilde{\delta}=1 given x
  x.alpha <- x %*% alphas
  prob <- exp(x.alpha)/(1 + exp(x.alpha))
  return(prob)
}

cost3 <- function(betas, y.hat.in.rct, x.in.rct, w.tilde) {
  # to obtrain beta_hat using loss function which targets optimal action rate
  action.contrast <- (sign(y.hat.in.rct) + 1)/2
  action.linear <- (sign(x.in.rct %*% betas) + 1)/2
  loss <- sum(w.tilde * (action.contrast - action.linear)**2)
  return(loss)
}

cost2 <- function(alphas, x.in.rct, x.in.rwe, N) {
  # to obtain weights of method(2)
  part1 <- sum(x.in.rct %*% as.matrix(alphas, length(alphas), 1))
  part2 <- mean(log(1 + exp(x.in.rwe %*% as.matrix(alphas, length(alphas), 1))))*N
  return(part2-part1)
} 

cost <- function(alphas, x.in.rct, left.hand, N) {
  # to obtain weights of method(3)
  right.hand <- apply(x.in.rct / c(delta.tilde.prob(x.in.rct, alphas)), 2, sum)/N
  diff <- right.hand - left.hand
  return(sqrt(mean(diff**2)))
}

linear_policy <- function(betas, x.in.rct, weight, label) {
  policy <- (x.in.rct %*% as.matrix(betas, length(betas), 1))>0
  cost <- mean(weight*(policy!=as.numeric(as.character(label))))
  return(cost)
}

policy <- function(weight, cls.meth, n, p, N, x.in.rct, x.in.rwe, data.classification, 
                   R.tilde, prob.one) {
  # weight=1: unweighted; weight=2: method(2): weight=3: method(3)
  # cls.meth: classification method
  # return: actions of target population assigned by the estimated policy
  if (weight==1) {
    w.tilde <- rep(1, n)
  } else if (weight==2) {
    opt <- optim(rep(0, p+1), cost2, x.in.rct=x.in.rct, x.in.rwe=x.in.rwe, N=N, 
                 method = "BFGS", control = list(maxit=10000))
    alphas.hat <- opt$par
    w.tilde <- 1/delta.tilde.prob(x.in.rct, alphas.hat)
  } else {
    opt <- optim(rep(0, p+1), cost, x.in.rct=x.in.rct, left.hand=left.hand, N=N,
                 method = "BFGS", control = list(maxit=10000))
    alphas.hat <- opt$par
    w.tilde <- 1/delta.tilde.prob(x.in.rct, alphas.hat)
  }
  if (cls.meth=="RF") {
    rule <- ranger(label ~ ., data=data.classification,
                   case.weights=w.tilde*abs(R.tilde/prob.one))
    action <- predict(rule, data=data.frame(x[, 2:(p+1)]))$predictions
  } else if (cls.meth=="tree") {
    rule <- rpart(label ~ ., data=data.classification, method = "class",
                  weights=w.tilde*abs(R.tilde/prob.one))
    action <- predict(rule, newdata=data.frame(x[, 2:(p+1)]), type="class")
  }
  return(action)
}

quality <- function(action, opt.action, opt.value, baseline, contrast) {
  # calculate accuracy of estimated optimal actions to true optimal acitons
  (rate <- mean(action==opt.action))
  # calculate true value of the estimated policy
  y.temp <- baseline + as.numeric(as.character(action))*contrast; 
  # calculate MSE between value of the estimated optimal policy and true optimal value
  (v.mse <- mean((y.temp-opt.value)^2))
  return(c(rate, v.mse))
}

cost_robust_weights <- function(lambda, x.in.rct, x.in.rwe, weights.raw) {
  p <- dim(x.in.rct)[2] 
  lambda <- as.matrix(lambda, p, 1)
  l <- weights.raw/exp(x.in.rct %*% lambda + 1)
  diff <- c(apply(c(l) * x.in.rct[, 2:p], 2, sum), sum(l)) - 
    c(apply(x.in.rwe, 2, mean)[2:p], 1)
  return(mean(diff**2))
}

robust_weights <- function(lambda, x.in.rct, weights.raw) {
  p <- dim(x.in.rct)[2] 
  lambda <- as.matrix(lambda, p, 1)
  robust_weights <- weights.raw/exp(x.in.rct %*% lambda+1)
  return(robust_weights)
}

cost_linear_rule <- function(beta, x.in.rct, label, weights) {
  linear_rule <- (x.in.rct %*% as.matrix(beta, length(beta), 1) > 0)
  return(sum(weights*(as.integer(as.character(label)) != linear_rule)))
}

learn_rule <- function(y.in.rct, x.in.rct, a.in.rct, x.in.rwe, x, N, y.star.hat.in.rct, 
                       w.method, method, prob.one, only.weight=F) {
  n <- dim(x.in.rct)[1]
  p <- dim(x.in.rct)[2] - 1
  if (w.method == 1) {
    ## unweighted
    w.tilde <- rep(1, n)
  } else if (w.method == 2) {
    ## weight (2)
    opt2 <- optim(rep(0,p+1), cost2, x.in.rct=x.in.rct, x.in.rwe=x.in.rwe, N=N, method = "BFGS", control = list(maxit=10000))
    (alphas.hat2 <- opt2$par)
    w.tilde <- 1 + 1/(exp(x.in.rct %*% as.matrix(alphas.hat2, length(alphas.hat2), 1)))
  } else if (w.method == 3) {
    ## weight (3)
    left.hand <- apply(x.in.rwe, 2, mean)
    opt <- optim(rep(0, p+1), cost, x.in.rct=x.in.rct, left.hand=left.hand, N=N, method = "BFGS", control = list(maxit=10000))
    (alphas.hat3 <- opt$par)
    w.tilde <- 1 + 1/(exp(x.in.rct %*% as.matrix(alphas.hat3, length(alphas.hat3), 1)))
  } else {
    print("Wrong input for weighting method")
  }
  if (only.weight) {
    ## only need to learn weight and return
    return(w.tilde)
  } else {
    if (method == "linear") {
      x.tilde <- x.in.rct * c(sqrt(w.tilde))
      (beta.hat.star <- solve(t(x.tilde) %*% x.tilde) %*% 
          as.matrix(apply(x.in.rct * y.star.hat.in.rct*c(w.tilde), 2, sum), n, 1))
      action <- (sign(x%*%beta.hat.star)+1)/2
    } else if (method == "tree") {
      R.hat <- ranger(y.in.rct ~ ., data=data.frame(cbind(y.in.rct, x.in.rct[,2:(p+1)])))$predictions
      R.tilde <- y.in.rct-R.hat
      label <- factor(((2*a.in.rct-1)*sign(R.tilde)+1)/2)
      data.classification <- data.frame(label, x.in.rct[,2:(p+1)], y.in.rct = y.in.rct,
                                        A=as.factor(a.in.rct)) 
      ## weighted method(2)
      rule <- rpart(label ~ ., data=data.classification[,1:(p+1)], method = "class",
                    weights=w.tilde*abs(R.tilde/prob.one))
      action <- predict(rule, newdata=data.frame(x[, 2:(p+1)]), type="class")
    } else {
      print("Wrong input for rule-learning method")
    }
    return(list(action=action, w.tilde=w.tilde))
  }
}

value_est <- function(y, x, a, actions, ps, ipw=T) {
  ## actions: actions given by the evaluation policy
  ## a: actions given by the behavior policy
  C <- (actions==a)
  if (ipw) {
    value <- mean(y*C / ps)
  } else {
    ## AIPW
    # lm_fit.out <- lm(y ~ x + a*x)
    (lm_fit.out <- lm(y ~ x + a*x-1-a))
    betas <- lm_fit.out$coefficients
    Q <- cbind(x, actions*x)%*% as.matrix(betas, length(betas),1)
    value <- mean(C*y/ps-(C-ps)/ps*Q)
  }
  return(list(value=value))
}

cv <- function(y.in.rct, x.in.rct, a.in.rct, x.in.rwe, N, nfold, y.star.hat.in.rct, prob.one,
               ipw=T, baseline, contrast, x, calculate.true.val=F) {
  n <- dim(x.in.rct)[1]
  m <- dim(x.in.rwe)[1]
  prop.train <- 0.5 ## proportion of training rct 
  values.multisplit <- values.true.multisplit <- matrix(NA, nrow=nfold, ncol=4)
  for (k in 1:nfold) {
    train.ind <- sample(1:n, size=round(n*prop.train))
    eval.ind <- setdiff(1:n, train.ind)
    train.ind.rwe <- sample(1:m, size=round(m*prop.train))
    eval.ind.rwe <- setdiff(1:m, train.ind.rwe)
    values <- true.values <- c()
    w.tilde <- learn_rule(y.in.rct[eval.ind], x.in.rct[eval.ind, ], a.in.rct[eval.ind], 
                          x.in.rwe[eval.ind.rwe, ], NA, N-round(N*prop.train), NA,
                          w.method=2, method=NA, only.weight=T)
    w.tilde.normalized <- w.tilde/sum(w.tilde) * length(eval.ind)
    for (w.method in 1:2) {
      for (method in c("linear", "tree")) {
        action <- learn_rule(y.in.rct[train.ind], x.in.rct[train.ind, ], a.in.rct[train.ind], 
                   x.in.rwe[train.ind.rwe, ], x.in.rct[eval.ind, ], round(N*prop.train), 
                   y.star.hat.in.rct[train.ind], w.method, method, prob.one=prob.one)$action
        value <- value_est(y.in.rct[eval.ind]*w.tilde.normalized, sweep(x.in.rct[eval.ind, ], MARGIN=1, w.tilde.normalized, `*`), 
                           a.in.rct[eval.ind], as.numeric(as.character(action)), ps=rep(0.5, n-round(n*prop.train)), ipw=ipw)$value
        if (is.na(value)) {
          print(w.tilde)
          print(paste(method, w.method))
          stop("invalid#####")
        }
          
        values <- c(values, value)
        ## true values
        if (calculate.true.val) {
          action <- learn_rule(y.in.rct[train.ind], x.in.rct[train.ind, ], a.in.rct[train.ind], 
                               x.in.rwe[train.ind.rwe, ], x, round(N*prop.train), y.star.hat.in.rct[train.ind], 
                               w.method, method, prob.one=prob.one)$action
          true.values <- c(true.values, mean(baseline + as.numeric(as.character(action))*contrast))
        } else {
          true.values <- c(true.values, 0)
        }
      }
    }
    values.multisplit[k, ] <- values
    values.true.multisplit[k, ] <- true.values
  }
  values <- colMeans(values.multisplit)
  return(list(num=which.max(values), values=values, true.values= colMeans(values.true.multisplit)))
}

kern_cov_balance <- function(x.in.rct, x.in.rwe) {
  n <- dim(x.in.rct)[1]
  treat <- c(rep(1, n), rep(0, dim(x.in.rwe)[1]))
  X <- rbind(x.in.rct, x.in.rwe)[,-1]
  #### \tilde{\delta}=1 ####
  # Sobolev kernel
  Xstd <- transform.sob(X)$Xstd # standardize X to [0,1]^p
  K <- getGram(Xstd) # get Gram matrix using Sobolev kernel
  # design a grid for the tuning parameter
  nlam <- 50
  lams <- exp(seq(log(1e-8), log(1), len=nlam))
  # compute weights for \tilde{\delta}=1
  fit1 <- ATE.ncb.SN(treat, K, lam1s=lams, traceit = FALSE )
  if (sum(fit1$warns)) cat("lambda bound warning!\n")
  return(fit1$w[1:n])
}

entr_balance <- function(tol, x.in.rct, x.in.rwe) {
  n <- dim(x.in.rct)[1] 
  w <- Variable(n)
  # entr: is the element-wise entropy atom
  objective <- Maximize(sum(entr(w))) 
  constraints <- list(w >= 0, sum(w) == 1)
  for (k in 1:p) {
    constraints <- c(constraints, abs(sum(w*x.in.rct[,k+1]) - mean(x.in.rwe[,k+1])) <= tol*sd(x.in.rwe[,k+1]))
  }
  prob <- Problem(objective, constraints)
  w.tilde <- rep(0, n)
  tryCatch(
    expr = {
      result <- solve(prob)
      w.tilde <- result$getValue(w)
      if (is.na(w.tilde)[1]) {
        w.tilde <- rep(0,n)
      }
    },
    error = function(e) {
      # print("error")
    }
  )
  return(w.tilde)
}

weight_tuned_entr <- function(tols, x.in.rct, x.in.rwe, prop=0.5, smps=10) {
  n <- dim(x.in.rct)[1]
  sds <- apply(x.in.rwe, 2, sd)
  means <- apply(x.in.rwe, 2, mean)
  cov.diff.bars <- c()
  w.tildes <- matrix(NA, n, length(tols))
  for (i in 1:length(tols)) {
    tol <- tols[i]
    w.tilde <- entr_balance(tol, x.in.rct, x.in.rwe)
    if (sum(w.tilde)==0) {
      cov.diff.bar <- 1e+10
    } else {
      w.tildes[,i] <- w.tilde
      cov.diffs <- c()
      for (s in 1:smps) {
        boot.ind <- sample.int(n, round(prop*n), replace=TRUE)
        x.in.rct.boot <- x.in.rct[boot.ind,]
        w.tilde.boot <- w.tilde[boot.ind]
        cov.diff <- (apply(x.in.rct.boot*w.tilde.boot, 2, sum)/sum(w.tilde.boot) -  means)[-1] / sds[-1]
        cov.diffs <- c(cov.diffs, sum(cov.diff**2))
      }
      cov.diff.bar <- mean(cov.diffs)
    }
    cov.diff.bars <- c(cov.diff.bars, cov.diff.bar)
  }
  tuned.tol <- tols[which.min(cov.diff.bars)]
  if (min(cov.diff.bars) < 1e+10) {
    return(w.tildes[, which.min(cov.diff.bars)])
  } else {
    stop('************ optimization not feasible *****************')
  }
}

est.contrast <- function(a, x, y, delta.tilde.is.one.ind, prob.one, meth=1) {
  ## estimate contrast function
  if (meth==1) {
    y.star.hat <- a*y/prob.one - (1-a)*y/(1-prob.one) # outcome adjusted
  } else if (meth==2) {
    N<-length(y)
    x<-as.matrix(x[delta.tilde.is.one.ind,]); x <- x[,-1]
    y <- y[delta.tilde.is.one.ind]
    a <- a[delta.tilde.is.one.ind]
    
    # loc.a1<-which(a==1)
    # loc.a0<-which(a==0)
    # dat <- data.frame(cbind(y, x))
    # ## use random forest to est Q
    # fit1 <- ranger(y~., data=dat[loc.a1,])
    # reg.mu1.rf <- predict(fit1, data = dat)$predictions
    # fit0 <- ranger(y~., data=dat[loc.a0,])
    # reg.mu0.rf <- predict(fit0, data = dat)$predictions
    # y.star.hat <- rep(0, N)
    # y.star.hat[delta.tilde.is.one.ind] <- reg.mu1.rf - reg.mu0.rf
    
    n <- length(delta.tilde.is.one.ind)
    dat <- data.frame(cbind(y,x,a))
    fit <- ranger(y~., data=dat)
    a <- rep(1,n); dat1 <- data.frame(cbind(y,x,a))
    reg.mu1.rf <- predict(fit, data = dat1)$predictions
    a <- rep(0,n); dat0 <- data.frame(cbind(y,x,a))
    reg.mu0.rf <- predict(fit, data = dat0)$predictions
    y.star.hat <- rep(0, N)
    y.star.hat[delta.tilde.is.one.ind] <- reg.mu1.rf - reg.mu0.rf
  } else if (meth==3) {
    N<-length(y)
    x<-as.matrix(x[delta.tilde.is.one.ind,]); x <- x[,-1]
    y <- y[delta.tilde.is.one.ind]
    A <- a[delta.tilde.is.one.ind]
    n <- length(delta.tilde.is.one.ind)
    aipw <- v.aipw(x,y,A)$aipw
    tau_jk <- rep(0, n)
    for (i in 1:n) {
      out <- v.aipw(x[-i,],y[-i],A[-i])
      tau_jk[i] <- n*aipw-(n-1)* out$aipw
    }
    y.star.hat <- rep(0, N)
    y.star.hat[delta.tilde.is.one.ind] <- tau_jk
  }
  return(y.star.hat)
}

## Value search methods involve maximization of nonsmooth objective functions, 
## require special techniques; e.g., a genetic algorithm (as in R rgenoud)
fn <- function(Betas, x.in.rct, y.star.hat, delta.tilde.is.one.ind, weight, obj.value=TRUE) {
  ## obj.value=TRUE: ojective function is for value; o.w., for optimal action rate
  Betas <- as.matrix(Betas, length(Betas), 1)
  obj <- ((x.in.rct %*% Betas > 0) - ( y.star.hat[delta.tilde.is.one.ind] > 0))**2 * c(weight)
  if (obj.value) obj = obj * abs(y.star.hat[delta.tilde.is.one.ind])
  obj <- sum(obj)
  return(obj)
}

v.aipw <- function(x,y,A,psm=1,Qm=2) {
  A.temp <- A
  ## aipw
  if (psm==1) {
    phat <- mean(A)
  } else if (psm==2) {
    ## use logistic regression to est PS
    alpha.est<-glm(A~x,family="binomial")$coeff
    lpshat<-cbind(1,x)%*%alpha.est
    phat<-1/(1+exp(-lpshat ))
  } 
  # else if (psm==3) {
  #   ## use gam to est PS
  #   Gam.object <- gam(A ~ s(x1) + s(x2) + s(x3) + s(x4), family=binomial, data=dat2)
  #   phat = predict(Gam.object, type="response", newdata=dat2)
  # }
  dat2 <- data.frame(cbind(y, x, A,  A*x))
  n <- length(y)
  if (Qm==1) {
    ## use RF to est Q in aipw
    dat <- data.frame(cbind(y,x,A))
    fit <- ranger(y~., data=dat)
    A <- A*0+1; dat1 <- data.frame(cbind(y,x,A))
    reg.mu1 <- predict(fit, data = dat1)$predictions
    A <- A*0; dat0 <- data.frame(cbind(y,x,A))
    reg.mu0 <- predict(fit, data = dat0)$predictions
  } else if (Qm==2) {
    ## use linear model to est Q in aipw
    lin <- lm(y ~ ., data=dat2)
    newdf <- data.frame(x, rep(0, n), x*0)
    names(newdf) <- names(dat2)[-1]
    reg.mu0 <- predict(lin, newdata = newdf)
    newdf <-  data.frame(x, rep(1, n), x)
    names(newdf) <- names(dat2)[-1]
    reg.mu1 <- predict(lin, newdata = newdf)
  } 
  # else if (Qm==3) {
  #   ## use GAM to est Q in aipw
  #   reg.mu1 <- reg.mu1.gam
  #   reg.mu0 <- reg.mu0.gam
  # }
  A <- A.temp
  psi.aug <- y*A/phat-y*(1-A)/(1-phat)+(reg.mu1*(1-A/phat)-reg.mu0*(1-(1-A)/(1-phat)) )
  aipw<-mean(psi.aug)
  return(list(aipw=aipw))
}
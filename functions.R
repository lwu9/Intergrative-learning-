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

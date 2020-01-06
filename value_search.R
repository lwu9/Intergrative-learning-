## plots to check what (estimated) contrast VS each dim of X looks like in RCT and RWE respectively
par(mfrow=c(p,3))
for(i in 1:p) {
  plot(x.in.rct[,i+1], y.star.hat[delta.tilde.is.one.ind])
  plot(x.in.rct[,i+1], contrast[delta.tilde.is.one.ind])
  plot(x.in.rwe[,i+1], contrast[delta.is.one.ind] )
}


for (i in 1:(p+1))
print(c(mean(x[,i]*x[,i]), 
            sum(x[delta.tilde.is.one.ind,i]*x[delta.tilde.is.one.ind,i]*c(w.tilde2))/N,
            sum(x[delta.tilde.is.one.ind,i]*x[delta.tilde.is.one.ind,i]*c(w.tilde3))/N,
            sum(x[delta.tilde.is.one.ind,i]*x[delta.tilde.is.one.ind,i]*c(w.tilde4)), 
            mean(x.in.rct[,i]*x.in.rct[,i]),
            sum(x[delta.tilde.is.one.ind,i]*x[delta.tilde.is.one.ind,i]*c(w.tilde.true))/N))
(aa <- t(x) %*% x/N);(bb <- apply(x * c(y.star), 2, mean));cc <- solve(aa) %*% as.matrix(bb, p+1, 1)
w.tilde <- 1 + 1/(exp(x.in.rct %*% as.matrix(alphas.hat2, length(alphas.hat2), 1))); w.tilde2 <- w.tilde
x.tilde <- x.in.rct * c(sqrt(w.tilde))
(aa2 <- t(x.tilde) %*% x.tilde / N);(bb2 <- apply(x.in.rct * y.star.hat[delta.tilde.is.one.ind]*c(w.tilde2), 2, sum)/N)
cc2 <- solve(aa2) %*% as.matrix(bb2, p+1, 1)
w.tilde <- 1 + 1/(exp(x.in.rct %*% as.matrix(alphas.hat3, length(alphas.hat3), 1))); w.tilde3 <- w.tilde
x.tilde <- x.in.rct * c(sqrt(w.tilde))
(aa3 <- t(x.tilde) %*% x.tilde / N);(bb3 <- apply(x.in.rct * y.star.hat[delta.tilde.is.one.ind]*c(w.tilde3), 2, sum)/N)
cc3 <- solve(aa3) %*% as.matrix(bb3, p+1, 1)
x.tilde <- x.in.rct * c(sqrt(w.tilde4))
(aa4 <- t(x.tilde) %*% x.tilde);(bb4 <- apply(x.in.rct * y.star.hat[delta.tilde.is.one.ind]*c(w.tilde4), 2, sum))
cc4 <- solve(aa4) %*% as.matrix(bb4, p+1, 1)
w.tilde = rep(1, n);x.tilde <- x.in.rct * c(sqrt(w.tilde))
(aa1 <- t(x.tilde) %*% x.tilde / n);(bb1 <- apply(x.in.rct * y.star.hat[delta.tilde.is.one.ind], 2, sum)/n)
cc1 <- solve(aa1) %*% as.matrix(bb1, p+1, 1)
x.tilde <- x.in.rct * c(sqrt(w.tilde.true))
(aa0 <- t(x.tilde) %*% x.tilde / N);(bb0 <- apply(x.in.rct * y.star.hat[delta.tilde.is.one.ind]*c(w.tilde.true), 2, sum)/N)
cc0 <- solve(aa0) %*% as.matrix(bb0, p+1, 1)
print(c(sum(abs(aa-aa2)),sum(abs(aa-aa3)),sum(abs(aa-aa4)),sum(abs(aa-aa1)),sum(abs(aa-aa0))) )
print(c(det(aa),det(aa2),det(aa3),det(aa4),det(aa1),det(aa0)))
print(c(sum(abs(solve(aa)-solve(aa2))),sum(abs(solve(aa)-solve(aa3))),sum(abs(solve(aa)-solve(aa4))),
        sum(abs(solve(aa)-solve(aa1))),sum(abs(solve(aa)-solve(aa0)))))
print(c(sum(abs(bb-bb2)),sum(abs(bb-bb3)),sum(abs(bb-bb4)),sum(abs(bb-bb1)),sum(abs(bb-bb0))) )
print(c(sum(abs(cc-cc2)),sum(abs(cc-cc3)),sum(abs(cc-cc4)),sum(abs(cc-cc1)),sum(abs(cc-cc0))) )


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
gen <- genoud(fn, nvars=p+1, max=FALSE,print.level=0,pop.size=5000,x.in.rct=x.in.rct, y.star.hat= y.star.hat,
        delta.tilde.is.one.ind=delta.tilde.is.one.ind, weight=w.tilde4, obj.value=T)
action <-  (sign(x %*% gen$par)+1)/2
(quality.linear.gen <- quality(action, opt.actions, opt.value, baseline, contrast))
gen <- genoud(fn, nvars=p+1, max=FALSE,print.level=0,pop.size=1000,x.in.rct=x, y.star.hat= y.star,
              delta.tilde.is.one.ind=1:N, weight=rep(1,N), obj.value=T)
action <-  (sign(x %*% gen$par)+1)/2
(quality.linear.gen <- quality(action, opt.actions, opt.value, baseline, contrast))


######
A <- a[delta.tilde.is.one.ind]
loc.a1<-which(A==1)
loc.a0<-which(A==0)
dat <- data.frame(cbind(y.in.rct, x.in.rct))
dat2 <- data.frame(cbind(y.in.rct, x.in.rct, A,  A*x.in.rct))

## use random forest to est Q
lm.y<- y.in.rct[loc.a1]
lm.x<- x.in.rct[loc.a1,]
lm.x<- as.matrix(lm.x)
# Gam.object <- gam(y ~ s(x1) + s(x2) + s(x3) + s(x4), data=dat[loc.a1, ])
# reg.mu1.gam = predict(Gam.object, type="response", newdata=dat)
fit1 <- ranger(y.in.rct~., data=dat[loc.a1,])
reg.mu1.rf <- predict(fit1, data = dat)$predictions

lm.y<- y.in.rct[loc.a0]
lm.x<- x.in.rct[loc.a0,]
lm.x<- as.matrix(lm.x)
# Gam.o bject <- gam(y ~ s(x1) + s(x2) + s(x3) + s(x4), data=dat[loc.a0, ])
# reg.mu0.gam = predict(Gam.object, type="response", newdata=dat)
fit0 <- ranger(y.in.rct~., data=dat[loc.a0,])
reg.mu0.rf <- predict(fit0, data = dat)$predictions
y.star.hat.rct <- reg.mu1.rf-reg.mu0.rf
y.star.hat2 <- rep(0, length(y.star))
y.star.hat2[delta.tilde.is.one.ind] <- y.star.hat.rct

## proposed method
method <- fit(x.in.rct[,-1],y.in.rct, A, psm=1, Qm = 2)
aipw <- method$aipw
tau_jk <- rep(0, n)
for (i in 1:n) {
  out <- fit(x.in.rct[-i,-1],y.in.rct[-i],A[-i], psm=1, Qm = 2)
  tau_jk[i] <- n*aipw-(n-1)* out$aipw
}
y.star.hat3 <- rep(0, length(y.star))
y.star.hat3[delta.tilde.is.one.ind] <- tau_jk

## plots to compare contrast estimation
par(mfrow=c(p,5))
for(i in 1:p) {
  plot(x.in.rct[,i+1], y.star.hat[delta.tilde.is.one.ind])
  plot(x.in.rct[,i+1], y.star.hat2[delta.tilde.is.one.ind])
  plot(x.in.rct[,i+1], y.star.hat3[delta.tilde.is.one.ind])
  plot(x.in.rct[,i+1], y.star[delta.tilde.is.one.ind])
  plot(x.in.rwe[,i+1], y.star[delta.is.one.ind] )
}

y.star.hat <- y.star.hat2

y.star.hat <- y.star.hat3
## plot estimated y.star VS y.star
plot(y.star[delta.tilde.is.one.ind]-(a*y/prob.one - (1-a)*y/(1-prob.one) )[delta.tilde.is.one.ind])
plot(y.star[delta.tilde.is.one.ind]- y.star.hat)

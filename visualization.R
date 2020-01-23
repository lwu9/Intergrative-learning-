## visualization for the boxplots
library(ggplot2)
ggplot(Results, aes(x=weight, y=Value.MSE, fill=rule)) + 
  facet_wrap( ~ n, ncol = 1) +# theme(legend.position = "none") +
  geom_boxplot()
ggsave("S:\\Documents\\lab2018\\transfer learning of dtr\\vMSE.eps")
ggplot(Results, aes(x=weight, y=opt.action.rate, fill=rule)) + 
  facet_wrap( ~ n, ncol = 1) +# theme(legend.position = "none") +
  geom_boxplot()

ggsave("S:\\Documents\\lab2018\\transfer learning of dtr\\accuracy.eps")


## visualization for CV
load("S:/Documents/lab2018/transfer learning of dtr/results/results_cv.Rdata")
library(arm)
par(mfrow=c(1,5))
for (i in 1:length(Results)) {
  results <- Results[[i]]
  learned.rule.mse <- matrix(unlist(results[,1]), byrow = T, ncol=4)
  cv.rule.value.est <- matrix(unlist(results[,2]), byrow = T, ncol=4)
  orders <- apply(learned.rule.mse, 1, order)
  cv.rule <- apply(cv.rule.value.est, 1, which.max)
  num <- 1:dim(orders)[1]
  p <- c(mean(cv.rule==orders[1,]), mean(cv.rule==orders[2,]), 
                                        mean(cv.rule==orders[3,]), mean(cv.rule==orders[4,]))
  discrete.histogram (num,p, main=paste0("n=", round(mean(unlist(results[,3])))), ylim=c(0,1))
}
## only compare weighted and unweighted methods instead of across all methods
par(mfrow=c(2,5))
for (i in 1:length(Results)) {
  results <- Results[[i]]
  learned.rule.mse <- matrix(unlist(results[,1]), byrow = T, ncol=4)
  cv.rule.value.est <- matrix(unlist(results[,2]), byrow = T, ncol=4)
  learned.rule.mse <- learned.rule.mse[,c(1,3)]
  cv.rule.value.est <- cv.rule.value.est[,c(1,3)]
  orders <- apply(learned.rule.mse, 1, order)
  cv.rule <- apply(cv.rule.value.est, 1, which.max)
  num <- 1:dim(orders)[1]
  p <- c(mean(cv.rule==orders[1,]), mean(cv.rule==orders[2,]))
  # p <- c(mean(cv.rule==orders[1,]), mean(cv.rule==orders[2,]), 
  #        mean(cv.rule==orders[3,]), mean(cv.rule==orders[4,]))
  discrete.histogram (num,p, main=paste0("n=", round(mean(unlist(results[,3])))), ylim=c(0,1))
}
for (i in 1:length(Results)) {
  results <- Results[[i]]
  learned.rule.mse <- matrix(unlist(results[,1]), byrow = T, ncol=4)
  cv.rule.value.est <- matrix(unlist(results[,2]), byrow = T, ncol=4)
  learned.rule.mse <- learned.rule.mse[,c(1,3)+1]
  cv.rule.value.est <- cv.rule.value.est[,c(1,3)+1]
  orders <- apply(learned.rule.mse, 1, order)
  cv.rule <- apply(cv.rule.value.est, 1, which.max)
  num <- 1:dim(orders)[1]
  p <- c(mean(cv.rule==orders[1,]), mean(cv.rule==orders[2,]))
  # p <- c(mean(cv.rule==orders[1,]), mean(cv.rule==orders[2,]), 
  #        mean(cv.rule==orders[3,]), mean(cv.rule==orders[4,]))
  discrete.histogram (num,p, main=paste0("n=", round(mean(unlist(results[,3])))), ylim=c(0,1))
}
par(mfrow=c(5,1))
Value.MSE <- rule <- n <- weight <- c()
for (i in 1:length(Results)) {
  results <- Results[[i]]
  learned.rule.mse <- matrix(unlist(results[,1]), byrow = T, ncol=4)
  cv.rule.value.est <- matrix(unlist(results[,2]), byrow = T, ncol=4)
  cv.rule <- apply(cv.rule.value.est, 1, which.max)
  cv.rule.mse <- c()
  nn <- length(cv.rule)
  for (k in 1:nn) {
    cv.rule.mse <- c(cv.rule.mse, learned.rule.mse[k,cv.rule[k]])    
  }
  Value.MSE <- c(Value.MSE, learned.rule.mse[,1], learned.rule.mse[,2],
                    learned.rule.mse[,3], learned.rule.mse[,4], cv.rule.mse)
  rule <- c(rule, rep("linear", nn), rep("tree", nn), rep("linear", nn), rep("tree", nn), rep("cv",nn))
  n <- c(n, rep(round(mean(unlist(results[,3]))), nn*5))
  weight <- c(weight, rep("unweight",2*nn),rep("weight",2*nn), rep("cv",nn))
  # boxplot(cv.rule.mse, main=paste0("n=", round(mean(unlist(results[,3])))), ylim=c(0,11))
}
df <- data.frame(weight = weight, Value.MSE=Value.MSE, rule=rule, n=n)
library(ggplot2)
ggplot(df, aes(x=weight, y=Value.MSE, fill=rule)) + 
  facet_wrap( ~ n, ncol = 1) +# theme(legend.position = "none") +
  geom_boxplot()


## visualization for CV1
load("S:/Documents/lab2018/transfer learning of dtr/results_all/results_cv_linear.Rdata")
Value.MSE <- rule <- n <- weight <- c()
for (i in 1:length(Results)) {
  results <- Results[[i]]
  learned.rule.mse <- matrix(unlist(results[,1]), byrow = T, ncol=2)
  cv.rule.value.est <- matrix(unlist(results[,2]), byrow = T, ncol=2)
  cv.rule <- apply(cv.rule.value.est, 1, which.max)
  cv.rule.mse <- c()
  nn <- length(cv.rule)
  for (k in 1:nn) {
    cv.rule.mse <- c(cv.rule.mse, learned.rule.mse[k,cv.rule[k]])    
  }
  Value.MSE <- c(Value.MSE, learned.rule.mse[,1], learned.rule.mse[,2],cv.rule.mse)
                 # learned.rule.mse[,3], learned.rule.mse[,4], cv.rule.mse)
  rule <- c(rule, rep("linear", nn), rep("linear", nn),  rep("cv",nn))
  n <- c(n, rep(round(mean(unlist(results[,3]))), nn*3))
  weight <- c(weight, rep("unweight", 1*nn),rep("weight",1*nn), rep("cv",nn))
  # boxplot(cv.rule.mse, main=paste0("n=", round(mean(unlist(results[,3])))), ylim=c(0,11))
}
df <- data.frame(weight = weight, Value.MSE=Value.MSE, rule=rule, n=n)
library(ggplot2)
ggplot(df, aes(x=weight, y=Value.MSE, fill=rule)) + 
  facet_wrap( ~ n, ncol = 1) +# theme(legend.position = "none") +
  geom_boxplot()


############## visualization for CV only compare weighted and unweighted
library(ggplot2)
load("S:/Documents/lab2018/transfer learning of dtr/results_all/cv_Results_linear.RData")
Value.MSE <- rule <- n <- weight <- c()
for (i in 1:length(Results)) {
  # for (method in c("linear","tree")) {
    results <- Results[[i]]
    learned.rule.mse <- matrix(unlist(results[,1]), byrow = T, ncol=4)
    cv.rule.value.est <- matrix(unlist(results[,2]), byrow = T, ncol=4)
    cv.rule <- apply(cv.rule.value.est, 1, which.max)
    cv.rule.mse <- c()
    nn <- length(cv.rule)
    for (k in 1:nn) {
      cv.rule.mse <- c(cv.rule.mse, learned.rule.mse[k,cv.rule[k]])    
    }
    Value.MSE <- c(Value.MSE, learned.rule.mse[,1], learned.rule.mse[,2],
                   learned.rule.mse[,3], learned.rule.mse[,4], cv.rule.mse)
    n <- c(n, rep(round(mean(unlist(results[,3]))), nn*(dim(learned.rule.mse)[2]+1)))
    weight <- c(weight, rep("unweight",nn),rep("w1",nn),rep("w2",nn), rep("nonpara",nn),rep("cv",nn))
    # boxplot(cv.rule.mse, main=paste0("n=", round(mean(unlist(results[,3])))), ylim=c(0,11))
  # }
}
weight <- factor(weight, levels = c("unweight","w1","w2","nonpara","cv"))
df <- data.frame(weight = weight, Value.MSE=Value.MSE, n=n)
ggplot(df, aes(x=weight, y=Value.MSE, fill=rule)) + 
  facet_wrap( ~ n, ncol = 1) +# theme(legend.position = "none") +
  geom_boxplot()


###### Note: when method is "tree", there is no need to use different y.hat.meth, since it doesn't need estimated contrasts
library(ggplot2)
path <- "S:/Documents/lab2018/transfer learning of dtr/Intergrative-learning--master/Intergrative-learning--master/"
setwd(path)
for (dbn in c("similar", "diff")) {
  for (method in c("linear.LS", "linear.VS")) {
    load(paste0(path,"results_cv_nfold30_method",method,"_",dbn,"_dist.Rdata"))
    Value.MSE <- rule <- n <- weight <- contrast.est.meth <- alpha1s <- c()
    for (i in 1:length(Results)) {
      # for (method in c("linear","tree")) {
      results <- Results[[i]]
      learned.rule.mse <- matrix(unlist(results[,1]), byrow = T, ncol=4)
      cv.rule.mse <- unlist(results[,3])
      nn <- length(cv.rule.mse)
      Value.MSE <- c(Value.MSE, learned.rule.mse[,1], learned.rule.mse[,2],
                     learned.rule.mse[,3], learned.rule.mse[,4], cv.rule.mse)
      n <- c(n, rep(round(mean(unlist(results[,4]))), nn*(dim(learned.rule.mse)[2]+1)))
      weight <- c(weight, rep("unweight",nn),rep("w1",nn),rep("w2",nn), rep("nonpara",nn),rep("cv",nn))
      Name <- names(Results)[i]
      alpha1 <- paste0("alpha1=-", stringr::str_extract(Name, "\\d+"))
      alpha1s <- c(alpha1s, rep(alpha1,  nn*(dim(learned.rule.mse)[2]+1)))
      yhatmeth <- paste0("contrast.est",tail(strsplit(Name,split='')[[1]],1))
      contrast.est.meth <- c(contrast.est.meth, rep(yhatmeth, nn*(dim(learned.rule.mse)[2]+1)))
      # boxplot(cv.rule.mse, main=paste0("n=", round(mean(unlist(results[,3])))), ylim=c(0,11))
      # }
    }
    weight <- factor(weight, levels = c("unweight","w1","w2","nonpara","cv"))
    contrast.est.meth <- factor(contrast.est.meth, levels=c("contrast.est1","contrast.est2","contrast.est3"))
    df <- data.frame(weight = weight, Value.MSE=Value.MSE, n=n, contrast.est.meth=contrast.est.meth,alpha1=alpha1s)
    gp <- ggplot(df, aes(x=weight, y=Value.MSE, fill=contrast.est.meth)) + scale_fill_manual(values=c("#999999", "#E69F00", "#56B4E9"))+
      facet_wrap( ~ alpha1, nrow = 1) + geom_boxplot() +ggtitle(paste0(method,"; ", dbn, " dbn"))+
      ## check how to add mean for each box
      stat_summary(fun.y = mean, color = "red", position = position_dodge(0.75),
                   geom = "point", shape = 20, size = 2.5,
                   show.legend = FALSE) ## add mean value in boxplot
    print(gp)
    ggsave(paste0(method,"_", dbn, ".png"))
  }
}

  
i=5
results <- Results[[i]]
learned.rule.mse <- matrix(unlist(results[,1]), byrow = T, ncol=4)
cv.rule.mse <- unlist(results[,3])
Value.MSE <- c(learned.rule.mse[,1], learned.rule.mse[,2],
               learned.rule.mse[,3], learned.rule.mse[,4], cv.rule.mse)
nn <- length(cv.rule.mse)
weight <- c(rep("unweight",nn),rep("w1",nn),rep("w2",nn), rep("nonpara",nn),rep("cv",nn))
df <- data.frame(weight = weight, Value.MSE=Value.MSE)
ggplot(df, aes(x=weight, y=Value.MSE))+geom_boxplot()+ggtitle(names(Results)[i])+
  stat_summary(fun.y=mean, geom="point", shape=20, size=4, color="red", fill="red") ## add mean value in boxplot

############### use RF in all Q #######################
dbn <- "diff"; method="linear.VS"
## fit RF of Y on all X and all A, get one Q(.,.) 
load(paste0(path,"results_cv_nfold30_method",method,"_",dbn,"_dist_RF.Rdata"))
## fit RF of Y on X for different A sperately, get Q_1(.) and Q_0(.)
# load(paste0(path,"results_cv_nfold1_method",method,"_",dbn,"_dist.Rdata"))
Value.MSE <- rule <- n <- weight <- contrast.est.meth <- alpha1s <- c()
for (i in 1:length(Results)) {
  # for (method in c("linear","tree")) {
  results <- Results[[i]]
  learned.rule.mse <- matrix(unlist(results[,1]), byrow = T, ncol=4)
  cv.rule.mse <- unlist(results[,3])
  nn <- length(cv.rule.mse)
  Value.MSE <- c(Value.MSE, learned.rule.mse[,1], learned.rule.mse[,2],
                 learned.rule.mse[,3], learned.rule.mse[,4], cv.rule.mse)
  n <- c(n, rep(round(mean(unlist(results[,4]))), nn*(dim(learned.rule.mse)[2]+1)))
  weight <- c(weight, rep("unweight",nn),rep("w1",nn),rep("w2",nn), rep("nonpara",nn),rep("cv",nn))
  Name <- names(Results)[i]
  alpha1 <- paste0("alpha1=-", stringr::str_extract(Name, "\\d+"))
  alpha1s <- c(alpha1s, rep(alpha1,  nn*(dim(learned.rule.mse)[2]+1)))
  yhatmeth <- paste0("contrast.est",tail(strsplit(Name,split='')[[1]],1))
  contrast.est.meth <- c(contrast.est.meth, rep(yhatmeth, nn*(dim(learned.rule.mse)[2]+1)))
  # boxplot(cv.rule.mse, main=paste0("n=", round(mean(unlist(results[,3])))), ylim=c(0,11))
  # }
}
weight <- factor(weight, levels = c("unweight","w1","w2","nonpara","cv"))
contrast.est.meth <- factor(contrast.est.meth, levels=c("contrast.est1","contrast.est2","contrast.est3"))
df <- data.frame(weight = weight, Value.MSE=Value.MSE, n=n, contrast.est.meth=contrast.est.meth,alpha1=alpha1s)
gp <- ggplot(df, aes(x=weight, y=Value.MSE, fill=contrast.est.meth)) + scale_fill_manual(values=c( "#E69F00", "#56B4E9"))+
  facet_wrap( ~ alpha1, nrow = 1) + geom_boxplot() +ggtitle(paste0(method,"; ", dbn, " dbn; RF"))+
  ## check how to add mean for each box
  stat_summary(fun.y = mean, color = "red", position = position_dodge(0.75),
               geom = "point", shape = 20, size = 2.5,
               show.legend = FALSE) ## add mean value in boxplot
print(gp)
ggsave(paste0(method,"_", dbn, "_RF.png"))
## ctr.mse
Value.MSE <- rule <- n <- ctr.mses <- contrast.est.meth <- alpha1s <- c()
for (i in 1:length(Results)) {
  # for (method in c("linear","tree")) {
  results <- Results[[i]]
  ctr.mse <- unlist(results[,5])
  nn <- length(ctr.mse)
  ctr.mses <- c(ctr.mses, ctr.mse)
  n <- c(n, rep(round(mean(unlist(results[,4]))), nn))
  Name <- names(Results)[i]
  alpha1 <- paste0("alpha1=-", stringr::str_extract(Name, "\\d+"))
  alpha1s <- c(alpha1s, rep(alpha1,  nn))
  yhatmeth <- paste0("contrast.est",tail(strsplit(Name,split='')[[1]],1))
  contrast.est.meth <- c(contrast.est.meth, rep(yhatmeth, nn))
  # boxplot(cv.rule.mse, main=paste0("n=", round(mean(unlist(results[,3])))), ylim=c(0,11))
  # }
}
contrast.est.meth <- factor(contrast.est.meth, levels=c("contrast.est2","contrast.est3"))
df <- data.frame(ctr.mses=ctr.mses, n=n, contrast.est.meth=contrast.est.meth,alpha1=alpha1s)
gp <- ggplot(df, aes(x=contrast.est.meth, y=ctr.mses)) + #scale_fill_manual(values=c( "#E69F00", "#56B4E9"))+
  facet_wrap( ~ alpha1, nrow = 1) + geom_boxplot() +ggtitle(paste0(method,"; ", dbn, " dbn; RF"))+
  ## check how to add mean for each box
  stat_summary(fun.y = mean, color = "red", position = position_dodge(0.75),
               geom = "point", shape = 20, size = 2.5,
               show.legend = FALSE) ## add mean value in boxplot
print(gp)
ggsave(paste0(method,"_", dbn, "_RF_ctrMSE.png"))

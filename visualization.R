## visualization for the boxplots
library(ggplot2)
path <- "/Users/lwu9/Documents/transfer_learner_dtr/"
setwd(path)
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


#path <- "S:/Documents/lab2018/transfer learning of dtr/Intergrative-learning--master/Intergrative-learning--master/"
path <- "/Users/lwu9/Documents/transfer_learner_dtr/"
for (method in c("linear.LS","linear.VS","tree")) {
  load(paste0(path,"results_cv_nfold30_method",method,"_diff_dist.Rdata"))
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
  # means <- aggregate(Value.MSE ~ weight + contrast.est.meth, df, mean)
  gp <- ggplot(df, aes(x=weight, y=Value.MSE, fill=contrast.est.meth)) + 
    facet_wrap( ~ alpha1, nrow = 1) + geom_boxplot() +ggtitle(method)+
    ## check how to add mean for each box
    # stat_summary(fun.y = mean, geom = "errorbar", 
                 # aes(ymax = ..y.., ymin = ..y.., group = contrast.est.meth),
                 # width = 0.75, linetype = "dashed", position = position_dodge())
    # stat_summary(fun.y=mean, geom="point",
    #     aes(ymax = ..y.., ymin = ..y.., group = contrast.est.meth),
    # #     shape=20, size=4, color="red", fill="red") ## add mean value in boxplot
    # geom_point(data = means, aes(y = Value.MSE, x =  weight), 
    #            position=position_dodge(width=.75), color = "white")
    # 
  print(gp)
}


## contrast estimates MSE
Value.MSE <- rule <- n <- weight <- contrast.est.meth <- alpha1s <- ctr.mse <- c()
for (i in 1:length(Results)) {
  # for (method in c("linear","tree")) {
  results <- Results[[i]]
  cv.rule.mse <- unlist(results[,3])
  nn <- length(cv.rule.mse)
  ctr.mse <- c(ctr.mse, unlist(results[,5]))
  Name <- names(Results)[i]
  alpha1 <- paste0("alpha1=-", stringr::str_extract(Name, "\\d+"))
  alpha1s <- c(alpha1s, rep(alpha1,  nn))
  yhatmeth <- paste0("contrast.est",tail(strsplit(Name,split='')[[1]],1))
  contrast.est.meth <- c(contrast.est.meth, rep(yhatmeth, nn))
  # boxplot(cv.rule.mse, main=paste0("n=", round(mean(unlist(results[,3])))), ylim=c(0,11))
  # }
}
contrast.est.meth <- factor(contrast.est.meth, levels=c("contrast.est1","contrast.est2","contrast.est3"))
df <- data.frame(ctr.mes=ctr.mse, contrast.est.meth=contrast.est.meth,alpha1=alpha1s)
alpha1s <- alpha1s[-which(df$contrast.est.meth=="contrast.est1")]
gp <- ggplot(df, aes(x=contrast.est.meth, y=ctr.mse)) + 
  facet_wrap( ~ alpha1, nrow = 1) + geom_boxplot() +ggtitle(method)+
  ## check how to add mean for each box
  stat_summary(fun.y=mean, geom="point", shape=20, size=4, color="red", fill="red") ## add mean value in boxplot
print(gp)


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

############ contain NN result ############
path <- "/Users/lwu9/Documents/transfer_learner_dtr/"
setwd(path)
nfold <- 5; method <- "linear.VS"; dbn <- "diff"
load(paste0(path,"results_cv_nfold", nfold, "_method",method,"_", dbn, "_dist_all_except_NN.Rdata"))
allResults <- Results
for (alpha1 in c(-8, -10)) {
  load(paste0(path, "results_cv_nfold", nfold, "_method", method, "_alpha", abs(alpha1), "_", dbn, "_dist_NN.Rdata"))
  allResults <- c(allResults, Results)
}
for (alpha1 in c(-8, -10)) {
  load(paste0(path, "results_cv_nfoldNA", "_method", method, "_alpha", abs(alpha1), "_", dbn, "_benchmark.Rdata"))
  allResults <- c(allResults, Results)
}
Results <- allResults; rm(allResults); print(names(Results))

## Qm=1: fit RF of Y on X for different A sperately, get Q_1(.) and Q_0(.)
## Qm=2: linear; Qm=3: NN
Value.MSE <- rule <- n <- weight <- contrast.est.meth <- alpha1s <- Q.model <- ctr.Q<- c()
for (i in 1:(length(Results)-2)) {
  results <- Results[[i]]
  learned.rule.mse <- matrix(unlist(results[,1]), byrow = T, ncol=4)
  cv.rule.mse <- unlist(results[,3])
  nn <- length(cv.rule.mse)
  Value.MSE <- c(Value.MSE, learned.rule.mse[,4], learned.rule.mse[,3],
                 learned.rule.mse[,2], learned.rule.mse[,1], cv.rule.mse)
  n <- c(n, rep(round(mean(unlist(results[,4]))), nn*(dim(learned.rule.mse)[2]+1)))
  weight <- c(weight, rep("unweight",nn),rep("w1",nn),rep("w2",nn), rep("nonpara",nn),rep("cv",nn))
  Name <- names(Results)[i]
  alpha1 <- paste0("alpha1=-", stringr::str_extract(Name, "\\d+"))
  alpha1s <- c(alpha1s, rep(alpha1,  nn*(dim(learned.rule.mse)[2]+1)))
  contrst.est_Q.est <- paste0("contrast", strsplit(Name,split='yhatmeth')[[1]][2])
  ctr.Q <- c(ctr.Q , rep( contrst.est_Q.est, nn*(dim(learned.rule.mse)[2]+1) ))
  # yhatmeth_Qm <- strsplit(strsplit(Name,split='yhatmeth')[[1]][2], split="Qm")[[1]]
  # yhatmeth <- paste0("contrast.est", yhatmeth_Qm[1])
  # contrast.est.meth <- c(contrast.est.meth, rep(yhatmeth, nn*(dim(learned.rule.mse)[2]+1)))
  # Qm <- paste0("Q.est", yhatmeth_Qm[2])
  # Q.model <- c(Q.model, rep(Qm, nn*(dim(learned.rule.mse)[2]+1) ))
}
for (i in (length(Results)-1):length(Results)) {
  ## the last two is for benchmark
  results <- Results[[i]]
  bm.mse <- unlist(results[,1])
  Value.MSE <- c(Value.MSE, bm.mse)
  nn <- length(bm.mse)
  n <- c(n, rep(round(mean(unlist(results[,2]))), nn))
  weight <-  c(weight, rep("benchmark", nn))
  Name <- names(Results)[i]
  alpha1 <- paste0("alpha1=-", stringr::str_extract(Name, "\\d+"))
  alpha1s <- c(alpha1s, rep(alpha1, nn))
  ctr.Q <- c(ctr.Q, rep("true.ctr_w", nn))
}
weight <- factor(weight, levels = c("unweight","w1","w2","nonpara","cv","benchmark"))
ctr.Q <- factor(ctr.Q, levels=c("contrast1Qm0","contrast2Qm1","contrast2Qm2","contrast2Qm3","contrast3Qm1",
                                            "contrast3Qm2","contrast3Qm3","contrast4Qm1","contrast4Qm2","contrast4Qm3",
                                            "true.ctr_w"))
df <- data.frame(weight = weight, Value.MSE=Value.MSE, n=n, ctr.Q=ctr.Q,alpha1=alpha1s)
gp <- ggplot(df, aes(x=weight, y=Value.MSE, fill= ctr.Q)) + #scale_fill_manual(values=c("#999999", "#E69F00", "#56B4E9"))+
  facet_wrap( ~ alpha1, nrow = 2) + geom_boxplot() +ggtitle(paste0(method,"; ", dbn, " dbn"))+
  ## check how to add mean for each box
  stat_summary(fun.y = mean, color = "red", position = position_dodge(0.75),
               geom = "point", shape = 20, size = 2.5,
               show.legend = FALSE) ## add mean value in boxplot
print(gp)
method <- gsub("\\.","_",method)
ggsave(paste0(method,"_", dbn, "_all.png"))


## ctr.mse
Value.MSE <- rule <- n <- ctr.mses <- contrast.est.meth <- alpha1s <- ctr.Q <- c()
for (i in 1:length(Results)) {
  # for (method in c("linear","tree")) {
  results <- Results[[i]]
  ctr.mse <- unlist(results[,5])
  nn <- length(ctr.mse)
  ctr.mses <- c(ctr.mses, ctr.mse)
  n <- c(n, rep(round(mean(unlist(results[,4]))), nn))
  Name <- names(Results)[i]
  
  print(paste("contrast est MSE:", mean(ctr.mse), "with sd", sd(ctr.mse), Name))
  alpha1 <- paste0("alpha1=-", stringr::str_extract(Name, "\\d+"))
  alpha1s <- c(alpha1s, rep(alpha1,  nn))
  yhatmeth <- paste0("contrast.est",tail(strsplit(Name,split='')[[1]],1))
  contrast.est.meth <- c(contrast.est.meth, rep(yhatmeth, nn))
  contrst.est_Q.est <- paste0("contrast", strsplit(Name,split='yhatmeth')[[1]][2])
  ctr.Q <- c(ctr.Q, rep( contrst.est_Q.est, nn) )
  # boxplot(cv.rule.mse, main=paste0("n=", round(mean(unlist(results[,3])))), ylim=c(0,11))
  # }
}
ctr.Q <- factor(ctr.Q, levels=c("contrast1Qm0","contrast2Qm1","contrast2Qm2","contrast2Qm3","contrast3Qm1",
                                            "contrast3Qm2","contrast3Qm3","contrast4Qm1","contrast4Qm2","contrast4Qm3"))
df <- data.frame(ctr.mses=ctr.mses, n=n, ctr.Q=ctr.Q,alpha1=alpha1s)
gp <- ggplot(df, aes(x=ctr.Q, y=ctr.mses)) + # scale_fill_manual(values=c( "#999999", "#E69F00", "#56B4E9"))+
  facet_wrap( ~ alpha1, nrow = 2) + geom_boxplot() +ggtitle(paste0(method,"; ", dbn, " dbn"))+
  ## check how to add mean for each box
  stat_summary(fun.y = mean, color = "red", position = position_dodge(0.75),
               geom = "point", shape = 20, size = 2.5,
               show.legend = FALSE) ## add mean value in boxplot
print(gp)
ggsave(paste0(method,"_", dbn, "_RF_ctrMSE.png"))



############# Mis-specify ##########
## what is true w.tilde in benchmark? Oh, I think it should use true delta.tilde model to get true w.tilde
library(ggplot2)
path <- "/Users/lwu9/Documents/transfer_learner_dtr/"
setwd(path)
nfold <- 5; method <- "linear.VS"; dbn <- "diff"
load(paste0(path,"results_cv_nfold", nfold, "_method",method,"_alpha7_", dbn, "_dist_all_except_NN_mis.Rdata"))
allResults <- Results
# for (alpha1 in c(-7)) {
#   load(paste0(path, "results_cv_nfold", nfold, "_method", method, "_alpha", abs(alpha1), "_", dbn, "_dist_NN.Rdata"))
#   allResults <- c(allResults, Results)
# }
for (alpha1 in c(-7)) {
  load(paste0(path, "results_cv_nfoldNA", "_method", method, "_alpha", abs(alpha1), "_", dbn, "_benchmark_mis.Rdata"))
  allResults <- c(allResults, Results)
}
Results <- allResults; rm(allResults); print(names(Results))

Value.MSE <- rule <- n <- weight <- contrast.est.meth <- alpha1s <- Q.model <- ctr.Q<- c()
for (i in 1:(length(Results)-1)) {
  results <- Results[[i]]
  learned.rule.mse <- matrix(unlist(results[,1]), byrow = T, ncol=4)
  cv.rule.mse <- unlist(results[,3])
  nn <- length(cv.rule.mse)
  Value.MSE <- c(Value.MSE, learned.rule.mse[,4], learned.rule.mse[,3],
                 learned.rule.mse[,2], learned.rule.mse[,1], cv.rule.mse)
  n <- c(n, rep(round(mean(unlist(results[,4]))), nn*(dim(learned.rule.mse)[2]+1)))
  weight <- c(weight, rep("unweight",nn),rep("w1",nn),rep("w2",nn), rep("nonpara",nn),rep("cv",nn))
  Name <- names(Results)[i]
  alpha1 <- paste0("alpha1=-", stringr::str_extract(Name, "\\d+"))
  alpha1s <- c(alpha1s, rep(alpha1,  nn*(dim(learned.rule.mse)[2]+1)))
  contrst.est_Q.est <- paste0("contrast", strsplit(Name,split='yhatmeth')[[1]][2])
  ctr.Q <- c(ctr.Q , rep( contrst.est_Q.est, nn*(dim(learned.rule.mse)[2]+1) ))
  # yhatmeth_Qm <- strsplit(strsplit(Name,split='yhatmeth')[[1]][2], split="Qm")[[1]]
  # yhatmeth <- paste0("contrast.est", yhatmeth_Qm[1])
  # contrast.est.meth <- c(contrast.est.meth, rep(yhatmeth, nn*(dim(learned.rule.mse)[2]+1)))
  # Qm <- paste0("Q.est", yhatmeth_Qm[2])
  # Q.model <- c(Q.model, rep(Qm, nn*(dim(learned.rule.mse)[2]+1) ))
}
for (i in (length(Results)):length(Results)) {
  ## the last two is for benchmark
  results <- Results[[i]]
  bm.mse <- unlist(results[,1])
  Value.MSE <- c(Value.MSE, bm.mse)
  nn <- length(bm.mse)
  n <- c(n, rep(round(mean(unlist(results[,2]))), nn))
  weight <-  c(weight, rep("benchmark", nn))
  Name <- names(Results)[i]
  alpha1 <- paste0("alpha1=-", stringr::str_extract(Name, "\\d+"))
  alpha1s <- c(alpha1s, rep(alpha1, nn))
  ctr.Q <- c(ctr.Q, rep("true.ctr_w", nn))
}

weight <- factor(weight, levels = c("unweight","w1","w2","nonpara","cv","benchmark"))
# ctr.Q <- factor(ctr.Q, levels=c("contrast1Qm0","contrast2Qm1","contrast2Qm2","contrast2Qm3","contrast3Qm1",
#                                 "contrast3Qm2","contrast3Qm3","contrast4Qm1","contrast4Qm2","contrast4Qm3",
#                                 "true.ctr_w"))
ctr.Q <- factor(ctr.Q, levels=c("contrast1Qm0","contrast2Qm1","contrast2Qm2","contrast3Qm1",
                                "contrast3Qm2","contrast4Qm1","contrast4Qm2","true.ctr_w"))
df <- data.frame(weight = weight, Value.MSE=Value.MSE, n=n, ctr.Q=ctr.Q,alpha1=alpha1s)
gp <- ggplot(df, aes(x=weight, y=Value.MSE, fill= ctr.Q)) + #scale_fill_manual(values=c("#999999", "#E69F00", "#56B4E9"))+
  facet_wrap( ~ alpha1, nrow = 2) + geom_boxplot() +ggtitle(paste0(method,"; ", dbn, " dbn"))+
  ## check how to add mean for each box
  stat_summary(fun.y = mean, color = "red", position = position_dodge(0.75),
               geom = "point", shape = 20, size = 2.5,
               show.legend = FALSE) ## add mean value in boxplot
print(gp)
method <- gsub("\\.","_",method)
ggsave(paste0(method,"_", dbn, "_except_NN_mis.png"))


########## when using "genoud" to optimize, use it three times and choose the most one #########
path <- "/Users/lwu9/Documents/transfer_learner_dtr/results/"
setwd(path)
nfold <- 5; method <- "linear.VS"; dbn <- "diff";  misspecify = T
allResults <- list()
for (alpha1 in c(-8, -10)) {
  load(paste0(path, "results_cv_nfold", nfold, "_method",method,"_alpha",abs(alpha1),
                     "_", dbn, "_dist_all_except_NN_mis", misspecify, "2.Rdata"))
  allResults <- c(allResults, Results)
}
Results <- allResults; rm(allResults); print(names(Results))

## Qm=1: fit RF of Y on X for different A sperately, get Q_1(.) and Q_0(.)
## Qm=2: linear; Qm=3: NN
Value.MSE <- rule <- n <- weight <- contrast.est.meth <- alpha1s <- Q.model <- ctr.Q<- c()
for (i in 1:length(Results)) {
  results <- Results[[i]]
  learned.rule.mse <- matrix(unlist(results[,1]), byrow = T, ncol=4)
  cv.rule.mse <- unlist(results[,3])
  nn <- length(cv.rule.mse)
  Value.MSE <- c(Value.MSE, learned.rule.mse[,4], learned.rule.mse[,3],
                 learned.rule.mse[,2], learned.rule.mse[,1], cv.rule.mse)
  n <- c(n, rep(round(mean(unlist(results[,4]))), nn*(dim(learned.rule.mse)[2]+1)))
  weight <- c(weight, rep("unweight",nn),rep("w1",nn),rep("w2",nn), rep("nonpara",nn),rep("cv",nn))
  Name <- names(Results)[i]
  alpha1 <- paste0("alpha1=-", stringr::str_extract(Name, "\\d+"))
  alpha1s <- c(alpha1s, rep(alpha1,  nn*(dim(learned.rule.mse)[2]+1)))
  contrst.est_Q.est <- paste0("contrast", strsplit(Name,split='yhatmeth')[[1]][2])
  ctr.Q <- c(ctr.Q , rep( contrst.est_Q.est, nn*(dim(learned.rule.mse)[2]+1) ))
  # yhatmeth_Qm <- strsplit(strsplit(Name,split='yhatmeth')[[1]][2], split="Qm")[[1]]
  # yhatmeth <- paste0("contrast.est", yhatmeth_Qm[1])
  # contrast.est.meth <- c(contrast.est.meth, rep(yhatmeth, nn*(dim(learned.rule.mse)[2]+1)))
  # Qm <- paste0("Q.est", yhatmeth_Qm[2])
  # Q.model <- c(Q.model, rep(Qm, nn*(dim(learned.rule.mse)[2]+1) ))
}

weight <- factor(weight, levels = c("unweight","w1","w2","nonpara","cv"))
ctr.Q <- factor(ctr.Q, levels=c("contrast1Qm0","contrast2Qm1","contrast2Qm2","contrast3Qm1",
                                "contrast3Qm2","contrast4Qm1","contrast4Qm2"))
df <- data.frame(weight = weight, Value.MSE=Value.MSE, n=n, ctr.Q=ctr.Q,alpha1=alpha1s)
gp <- ggplot(df, aes(x=weight, y=Value.MSE, fill= ctr.Q)) + #scale_fill_manual(values=c("#999999", "#E69F00", "#56B4E9"))+
  facet_wrap( ~ alpha1, nrow = 2) + geom_boxplot() +ggtitle(paste0(method,"; ", dbn, " dbn",  "; Specify=", !misspecify))+
  ## check how to add mean for each box
  stat_summary(fun.y = mean, color = "red", position = position_dodge(0.75),
               geom = "point", shape = 20, size = 2.5,
               show.legend = FALSE) ## add mean value in boxplot
print(gp)
method <- gsub("\\.","_",method)
ggsave(paste0(method,"_", dbn, "_specify_", !misspecify,  "_all_except_NN.png"))



########## when using "genoud" to optimize, use it three times and choose the most one, weighted.para=TRUE, the evaluation part in CV has changed to the same weighting methods--nonpara #########
path <- "/Users/lwu9/Documents/transfer_learner_dtr/results/"
setwd(path)
nfold <- 5; method <- "linear.VS"; dbn <- "diff";  misspecify = F
allResults <- list()
for (alpha1 in c(-8, -10)) {
  load(paste0(path, "results_cv_nfold", nfold, "_method",method,"_alpha",abs(alpha1),
              "_", dbn, "dist_all_except_NN_mis", misspecify, "_CVTRUE.Rdata"))
  allResults <- c(allResults, Results)
}
for (alpha1 in c(-8, -10)) {
  load(paste0(path, "results_cv_nfoldNA", "_method", method, "_alpha", abs(alpha1), "_", dbn, "_benchmark.Rdata"))
  allResults <- c(allResults, Results) 
}
Results <- allResults; rm(allResults); print(names(Results))

## Qm=1: fit RF of Y on X for different A sperately, get Q_1(.) and Q_0(.)
## Qm=2: linear; Qm=3: NN
Value.MSE <- rule <- n <- weight <- contrast.est.meth <- alpha1s <- Q.model <- ctr.Q<- c()
for (i in 1:(length(Results)-2)) {
  results <- Results[[i]]
  learned.rule.mse <- matrix(unlist(results[,1]), byrow = T, ncol=4)
  cv.rule.mse <- unlist(results[,4])
  nn <- length(cv.rule.mse)
  Value.MSE <- c(Value.MSE, learned.rule.mse[,4], learned.rule.mse[,3],
                 learned.rule.mse[,2], learned.rule.mse[,1], cv.rule.mse)
  n <- c(n, rep(round(mean(unlist(results[,2]))), nn*(dim(learned.rule.mse)[2]+1)))
  weight <- c(weight, rep("unweight",nn),rep("w1",nn),rep("w2",nn), rep("nonpara",nn),rep("cv",nn))
  Name <- names(Results)[i]
  alpha1 <- paste0("alpha1=-", stringr::str_extract(Name, "\\d+"))
  alpha1s <- c(alpha1s, rep(alpha1,  nn*(dim(learned.rule.mse)[2]+1)))
  contrst.est_Q.est <- paste0("contrast", strsplit(Name,split='yhatmeth')[[1]][2])
  ctr.Q <- c(ctr.Q , rep( contrst.est_Q.est, nn*(dim(learned.rule.mse)[2]+1) ))
  # yhatmeth_Qm <- strsplit(strsplit(Name,split='yhatmeth')[[1]][2], split="Qm")[[1]]
  # yhatmeth <- paste0("contrast.est", yhatmeth_Qm[1])
  # contrast.est.meth <- c(contrast.est.meth, rep(yhatmeth, nn*(dim(learned.rule.mse)[2]+1)))
  # Qm <- paste0("Q.est", yhatmeth_Qm[2])
  # Q.model <- c(Q.model, rep(Qm, nn*(dim(learned.rule.mse)[2]+1) ))
}
for (i in (length(Results)-1):length(Results)) {
  ## the last two is for benchmark
  results <- Results[[i]][1:200, ] ## to have the same num of reps with the other results
  bm.mse <- unlist(results[,1])
  Value.MSE <- c(Value.MSE, bm.mse)
  nn <- length(bm.mse)
  n <- c(n, rep(round(mean(unlist(results[,2]))), nn))
  weight <-  c(weight, rep("benchmark", nn))
  Name <- names(Results)[i]
  alpha1 <- paste0("alpha1=-", stringr::str_extract(Name, "\\d+"))
  alpha1s <- c(alpha1s, rep(alpha1, nn))
  ctr.Q <- c(ctr.Q, rep("true.ctr_w", nn))
}

weight <- factor(weight, levels = c("unweight", "w1", "w2", "nonpara", "cv", "benchmark"))
ctr.Q <- factor(ctr.Q, levels=c("contrast2Qm1weighted.paraTRUE","contrast2Qm2weighted.paraTRUE","contrast3Qm1weighted.paraTRUE",
                                "contrast3Qm2weighted.paraTRUE","contrast4Qm1weighted.paraTRUE","contrast4Qm2weighted.paraTRUE",
                                "true.ctr_w"))
df <- data.frame(weight = weight, Value.MSE=Value.MSE, n=n, ctr.Q=ctr.Q, alpha1=alpha1s)
gp <- ggplot(df, aes(x=weight, y=Value.MSE, fill= ctr.Q)) + #scale_fill_manual(values=c("#999999", "#E69F00", "#56B4E9"))+
  facet_wrap( ~ alpha1, nrow = 2, scales = "free_y") + geom_boxplot() +ggtitle(paste0(method,"; ", dbn, " dbn",  "; Specify=", !misspecify))+
  ## check how to add mean for each box
  stat_summary(fun.y = mean, color = "red", position = position_dodge(0.75),
               geom = "point", shape = 20, size = 2.5,
               show.legend = FALSE) ## add mean value in boxplot
print(gp)
method <- gsub("\\.","_",method)
ggsave(paste0(method,"_", dbn, "_specify_", !misspecify,  "_all_except_NN_weighted_paraTRUE.png"))




########## when using "genoud" to optimize, use it 30 times and choose the most one, weighted.para=TRUE, the evaluation part in CV has changed to the same weighting methods--nonpara #########
path <- "/Users/lwu9/Documents/transfer_learner_dtr/results/"
setwd(path)
nfold <- 5; method <- "linear.VS"; dbn <- "diff";  misspecify = T
allResults <- list()
for (alpha1 in c(-8, -10)) {
  load(paste0(path, "results_cv_nfold", nfold, "_method",method,"_alpha",abs(alpha1),
              "_", dbn, "_iter30_", misspecify, "_CVTRUE.Rdata"))
  allResults <- c(allResults, Results)
}
Results <- allResults; rm(allResults); print(names(Results))

## Qm=1: fit RF of Y on X for different A sperately, get Q_1(.) and Q_0(.)
## Qm=2: linear; Qm=3: NN
Value.MSE <- rule <- n <- weight <- contrast.est.meth <- alpha1s <- Q.model <- ctr.Q<- c()
for (i in 1:(length(Results))) {
  results <- Results[[i]]
  learned.rule.mse <- matrix(unlist(results[,1]), byrow = T, ncol=4)
  cv.rule.mse <- unlist(results[,4])
  nn <- length(cv.rule.mse)
  Value.MSE <- c(Value.MSE, learned.rule.mse[,4], learned.rule.mse[,3],
                 learned.rule.mse[,2], learned.rule.mse[,1], cv.rule.mse)
  n <- c(n, rep(round(mean(unlist(results[,2]))), nn*(dim(learned.rule.mse)[2]+1)))
  weight <- c(weight, rep("unweight",nn),rep("w1",nn),rep("w2",nn), rep("nonpara",nn),rep("cv",nn))
  Name <- names(Results)[i]
  alpha1 <- paste0("alpha1=-", stringr::str_extract(Name, "\\d+"))
  alpha1s <- c(alpha1s, rep(alpha1,  nn*(dim(learned.rule.mse)[2]+1)))
  contrst.est_Q.est <- paste0("contrast", strsplit(Name,split='yhatmeth')[[1]][2])
  ctr.Q <- c(ctr.Q , rep( contrst.est_Q.est, nn*(dim(learned.rule.mse)[2]+1) ))
  # yhatmeth_Qm <- strsplit(strsplit(Name,split='yhatmeth')[[1]][2], split="Qm")[[1]]
  # yhatmeth <- paste0("contrast.est", yhatmeth_Qm[1])
  # contrast.est.meth <- c(contrast.est.meth, rep(yhatmeth, nn*(dim(learned.rule.mse)[2]+1)))
  # Qm <- paste0("Q.est", yhatmeth_Qm[2])
  # Q.model <- c(Q.model, rep(Qm, nn*(dim(learned.rule.mse)[2]+1) ))
}

weight <- factor(weight, levels = c("unweight", "w1", "w2", "nonpara", "cv"))
ctr.Q <- factor(ctr.Q, levels=c("contrast1Qm0weighted.paraTRUE","contrast2Qm1weighted.paraTRUE","contrast2Qm2weighted.paraTRUE",
                                "contrast3Qm1weighted.paraTRUE","contrast3Qm2weighted.paraTRUE"))
df <- data.frame(weight = weight, Value.MSE=Value.MSE, n=n, ctr.Q=ctr.Q, alpha1=alpha1s)
gp <- ggplot(df, aes(x=weight, y=Value.MSE, fill= ctr.Q)) + #scale_fill_manual(values=c("#999999", "#E69F00", "#56B4E9"))+
  facet_wrap( ~ alpha1, nrow = 2, scales = "free_y") + geom_boxplot() +ggtitle(paste0(method,"; ", dbn, " dbn",  "; Specify=", !misspecify))+
  ## check how to add mean for each box
  stat_summary(fun.y = mean, color = "red", position = position_dodge(0.75),
               geom = "point", shape = 20, size = 2.5,
               show.legend = FALSE) ## add mean value in boxplot
print(gp)
method <- gsub("\\.","_",method)
ggsave(paste0(method,"_", dbn, "_specify_", !misspecify, "_iter30_weighted_paraTRUE.png"))


########## Smooth ramp loss to optimize, weighted.para=TRUE, the evaluation part in CV has changed to the same weighting methods--nonpara #########
library(ggplot2)
path <- "/Users/lwu9/Documents/transfer_learner_dtr/results/"
setwd(path)
nfold <- 10; method <- "linear.VS"; dbn <- "diff";  misspecify = T
allResults <- list()
for (alpha1 in c(-8, -9,-10)) {
  load(paste0(path, "sramp_results_cv_nfold", nfold, "_method",method,"_alpha",abs(alpha1),
              "_", dbn, "dist_mis", misspecify, "_CVTRUE.Rdata"))
  allResults <- c(allResults, Results)
}
load("results_benchmark.Rdata")
Results <- Results[grepl(dbn, names(Results))]
allResults <- c(allResults, Results) 
for (alpha1 in c(-8, -9,-10)) {
  load(paste0("results_alpha", abs(alpha1), "_dist_", dbn, "_Rep200_imai.Rdata"))
  allResults <- c(allResults, Results)
}
Results <- allResults; rm(allResults); print(names(Results))

## Qm=1: fit RF of Y on X for different A sperately, get Q_1(.) and Q_0(.)
## Qm=2: linear
Value.MSE <- rule <- n <- weight <- contrast.est.meth <- alpha1s <- Q.model <- ctr.Q<- c()
for (i in 1:15) {
  results <- Results[[i]]
  learned.rule.mse <- matrix(unlist(results[,1]), byrow = T, ncol=4)
  cv.rule.mse <- unlist(results[,4])
  nn <- length(cv.rule.mse)
  Value.MSE <- c(Value.MSE, learned.rule.mse[,4], learned.rule.mse[,3],
                 learned.rule.mse[,2], learned.rule.mse[,1], cv.rule.mse)
  n <- c(n, rep(round(mean(unlist(results[,2]))), nn*(dim(learned.rule.mse)[2]+1)))
  weight <- c(weight, rep("unweight",nn),rep("w1",nn),rep("w2",nn), rep("nonpara",nn),rep("cv",nn))
  Name <- names(Results)[i]
  alpha1 <- paste0("alpha1 = -", stringr::str_extract(Name, "\\d+"))
  alpha1s <- c(alpha1s, rep(alpha1,  nn*(dim(learned.rule.mse)[2]+1)))
  contrst.est_Q.est <- paste0("ctr", strsplit(strsplit(Name,split='yhatmeth')[[1]][2], "weighted")[[1]][1])
  ctr.Q <- c(ctr.Q , rep( contrst.est_Q.est, nn*(dim(learned.rule.mse)[2]+1) ))
}
for (i in 16:18) {
  ## benchmark
  results <- Results[[i]] ## to have the same num of reps with the other results
  bm.mse <- unlist(results[,1])
  Value.MSE <- c(Value.MSE, bm.mse)
  nn <- length(bm.mse)
  n <- c(n, rep(round(mean(unlist(results[,2]))), nn))
  weight <-  c(weight, rep("benchmark", nn))
  Name <- names(Results)[i]
  alpha1 <- paste0("alpha1 = -", stringr::str_extract(Name, "\\d+"))
  alpha1s <- c(alpha1s, rep(alpha1, nn))
  ctr.Q <- c(ctr.Q, rep("true.ctr_w", nn))
}
for (i in 19:21) {
  ## Imai as the comparison method
  results <- Results[[i]] ## to have the same num of reps with the other results
  imai.mse <- unlist(results[,1])
  Value.MSE <- c(Value.MSE, imai.mse)
  nn <- length(imai.mse)
  n <- c(n, rep(round(mean(unlist(results[,2]))), nn))
  weight <-  c(weight, rep("Imai", nn))
  Name <- names(Results)[i]
  alpha1 <- paste0("alpha1 = -", stringr::str_extract(Name, "\\d+"))
  alpha1s <- c(alpha1s, rep(alpha1, nn))
  ctr.Q <- c(ctr.Q, rep("Imai", nn))
}

weight <- factor(weight, levels = c("unweight", "w1", "w2", "nonpara", "cv", "benchmark", "Imai"))
ctr.Q <- factor(ctr.Q, levels=c( "ctr1Qm0", "ctr2Qm1", "ctr2Qm2",
                                 "ctr3Qm1", "ctr3Qm2", "true.ctr_w", "Imai"))
alpha1s <- factor(alpha1s, levels=c("alpha1 = -10","alpha1 = -9","alpha1 = -8" ))
df <- data.frame(weight = weight, Value.MSE=Value.MSE, n=n, ctr.Q=ctr.Q, alpha1=alpha1s)
gp <- ggplot(df, aes(x=weight, y=Value.MSE, fill= ctr.Q)) + #scale_fill_manual(values=c("#999999", "#E69F00", "#56B4E9"))+
  theme(legend.position="bottom", legend.margin=margin(t=-20), axis.text.x = element_text(angle = 30), legend.box = "horizontal", axis.title.x=element_blank())+
  scale_fill_discrete(NULL)+ guides(fill = guide_legend(nrow = 1))+
  facet_wrap( ~ alpha1, nrow = 1, scales = "free_y") + geom_boxplot() +ggtitle(paste0(method,"; ", dbn, " dbn",  "; Specify=", !misspecify))+
  ## check how to add mean for each box
  stat_summary(fun.y = mean, color = "red", position = position_dodge(0.75),
               geom = "point", shape = 20, size = 2.5,
               show.legend = FALSE) ## add mean value in boxplot
print(gp)
method <- gsub("\\.","_",method)
ggsave(paste0(method,"_", dbn, "_specify_", !misspecify,  "_all_except_NN_weighted_paraTRUE.png"))






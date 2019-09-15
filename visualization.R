## visualization
library(ggplot2)
ggplot(Results, aes(x=weight, y=Value.MSE, fill=rule)) + 
  facet_wrap( ~ n, ncol = 1) +# theme(legend.position = "none") +
  geom_boxplot()
ggsave("S:\\Documents\\lab2018\\transfer learning of dtr\\vMSE.eps")
ggplot(Results, aes(x=weight, y=opt.action.rate, fill=rule)) + 
  facet_wrap( ~ n, ncol = 1) +# theme(legend.position = "none") +
  geom_boxplot()

ggsave("S:\\Documents\\lab2018\\transfer learning of dtr\\accuracy.eps")

setwd("~/Google\ Drive/Studium/Math_VO/Abschlussaufgabe/modelling_VO/")
rm(list=ls())

library(ggplot2)

# load data 
tab <- read.table("micro.dat",header=TRUE)
ntrial <- max(tab$trial)

# time window
t1 = -250
t2 = 600
alpha = 1/30

# generate plot
p <- ggplot(tab,aes(x=onset,y=trial,color=factor(subj))) + 
  geom_point(size=1.5,shape=1) + xlim(c(-400,600)) + ylim(c(0,ntrial)) + 
  xlab("Microsaccade onset time t [ms]") + ylab("Trial number") 
print(p)

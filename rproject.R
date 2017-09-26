#setwd("~/Documents/papers/MathModBook/Exercises-Chap2/R/SAWmodel")
#setwd("~/Google\ Drive/Studium/Math_VO/Abschlussaufgabe/modelling_VO/")
rm(list=ls())

library(ggplot2)
library(grid)
library (stats)

sacRate <- dget("sacRate.R")

T = 600 + 200 #post run time for plotting
L <- 51
pre_run_trans <- 400 #pre run time for plotting
trans <- 10000 + pre_run_trans 
gamma <- 0.940 #0.954 worked quite good #decrease to lower microsacade rate
slope <- 1*L

tp <- seq(-pre_run_trans+200, T) 

#saccades
sac <- matrix(, nrow = 0, ncol = 3)
trial <- 1
subj <- 1 

max_trial <- 1000

#time-dependent parameters
tau_p <- 150
tau_A <- 30
K <- 1
p_1 <- 0.0002
p_2 <- 0.00008 #originally 0.2 determines the width of the parabola of spike in a_A
beta <- 0.3 #orginally 0.3
lambda1 <- 1.2 #originally 0.2 --> peak size of a_p
lambda2 <- 40.5 #orginally 0.7 --> adjusted so that a_A has peak size to approx. 0.8

#time-varying functions and parameters

t <- seq(1, T) 
trans_vec <- c(rep(1,trans))
a_p1 <- 1/(1 + lambda1 * exp( - p_1 * (t - tau_p)^2)) # a=b
a_p <- c(trans_vec, a_p1)

f <- 1/(1 + lambda1 * exp( - p_1 * (t - 2 * tau_p)^2))

D <- lambda2 * ( (p_2^K) / factorial(K) ) * (t - (tau_p + tau_A))^K * exp(- p_2 * (t - (tau_p + tau_A))^2)

a_A1 <- 1/(D + 1)
a_A2 <- c(c(rep(1,179)), a_A1[180:T])
a_A <- c(trans_vec, a_A2)

a <- a_p

thres1 = 1/(1 + beta * ((1 - f) + ( 1 - a_A2)))
thres <- c(trans_vec, thres1)

plotf <- function(x, m, pot){
  # create data frames
  traj.df <- data.frame(x=x[,1],y=x[,2])
  im <- m + pot
  df <- expand.grid(x=15:35, y=15:35)
  df$z <- as.vector(im[15:35,15:35])
  
  # create plot
  p1 <- ggplot(data=traj.df,aes(x=x,y=y)) + 
    geom_raster(data=df,aes(x, y, fill=z)) +
    scale_fill_gradientn(colours=c("blue","red","orange","yellow")) +
    geom_path(color="white") + 
    coord_fixed(ratio=1) + scale_x_continuous(limits=c(14.5,35.5)) + 
    scale_y_continuous(limits=c(14.5,35.5))
  
  # print plot
  pushViewport(viewport(layout=grid.layout(1,1)))
  print(p1,vp=viewport(layout.pos.row=1,layout.pos.col=1))
}

for (trial in 1:max_trial)
{
  # prepare matrices
  pot <- matrix(0, nrow=L, ncol=L)
  c <- (L-1)/2
  for (i in 1:L){
    for (j in 1:L){
      pot[i,j] = slope*( ((i-c)/c)^2 + ((j-c)/c)^2 ) #potential without time-varying factor
    }
  }
  vec <- matrix( c(-1,0,1,0,0,-1,0,1), nrow=4, ncol=2, byrow = TRUE)
  
  # prepare variables
  m <- matrix(runif(L*L),L)
  m1 <- m
  x <-  matrix(0, nrow=trans+1, ncol=2)
  x[1,] <- c(1,1)*((L-1)/2)
  i <- 0
  
  # iteration loop
  for ( j in 1:(trans+T) )
  {
    if ( j==trans+1 )
    {    
      vi <- x[i+1,]
      x <-  matrix(0, nrow=T+1, ncol=2)
      x[1,] <- vi
      i <- 0
    }
    debug_var <- m1[x[i+1,1],x[i+1,2]] #XXX delete
    
    #if(j >= trans+1)
    #{
    #  plotf(x, m, pot)
    #  print("yeah")
    #}
    
    pot1 <- pot * a[j] #apply time-varying factor a=b to potential
    i <- i+1
     
    mx_temp <- m1[x[i,1],x[i,2]] #save value where no no relexation is applied
    m1 <- gamma*m #relaxation
    
    m1[x[i,1],x[i,2]] <- mx_temp + 1 #increase activation at h_ij
    m2 <- m1 + pot1
    
    ri <- sample(4)
    idx <- matrix(x[i,],  nrow=4, ncol=2, byrow = TRUE)+ vec[ri,]
    vi <- which(m2[idx] == min(m2[idx]))
      
    next_x <-  c(x[i,1] + vec[ri[vi],1], x[i,2] + vec[ri[vi],2] ) #get next possible step
    
    if(m1[next_x[1],next_x[2]] > thres[j])
    {
      #select and go to the global minimum
      vi <- which(m2 == min(m2), arr.ind = TRUE) #returns vector with row & column indices of global minimum
      x[i+1,] <-vi #save position
      #save microsaccade for plotting
      if(j >= (trans-pre_run_trans+1))
      {
        sac <- rbind(sac, c(j-trans, subj, trial))
      }
    }  
    else
    { 
      x[i+1,] <-  c(x[i,1] + vec[ri[vi[1]],1], x[i,2] + vec[ri[vi[1]],2] )
    }
    
    m <- m1
  }
}
# jitter trajectory (for better visualization)
#x <- x + matrix(data=0.5*runif(2*length(x[,1])),ncol=2)

plotf(x, m, pot)
sacRate(sac, -pre_run_trans, T)



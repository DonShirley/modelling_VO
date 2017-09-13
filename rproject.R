setwd("~/Documents/papers/MathModBook/Exercises-Chap2/R/SAWmodel")
rm(list=ls())

library(ggplot2)
library(grid)
library (stats)

T = 30
L <- 51
trans <- 0
gamma <- 0.999
slope <- 1.0*L

#time-dependent parameters
tau_p <- 150
tau_A <- 30
K <- 1
p_1 <- 0.0002
p_2 <- 0.02
beta <- 0.3
lambda1 <- 0.2
lambda2 <- 0.7

#time-varying functions and parameters

t <- seq(1, T + trans)


a <- 1/(1 + lambda1 * exp( - p_1 * (t - tau_p)^2)) # a=b

f <- 1/(1 + lambda1 * exp( - p_1 * (t - 2 * tau_p)^2))

D <- lambda2 * ( (p_2^K) / factorial(K) ) * (t - (tau_p + tau_A))^K * exp(- p_2 * (t - (tau_p + tau_A)))

a_A <- 1/(1 + D)

thres = 1/(1 + beta * ((1 - f + (  1 - a_A))))

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
for ( j in 1:(trans+T) ){
  if ( j==trans+1 ){    
    vi <- x[i+1,]
    x <-  matrix(0, nrow=T+1, ncol=2)
    x[1,] <- vi
    i <- 0
  }
  pot1 <- pot * a[j] #apply time-varying factor a=b to potential
  i <- i+1
   
  m1 <- gamma*m #relaxation
  m1[x[i,1],x[i,2]] <- m[x[i,1],x[i,2]]/gamma + 1 #increase activation at h_ij
  m2 <- m1 + pot1
  if( m1[x[i,1],x[i,2]] > thres[j]){ #check threshold
    #select and go to the global minimum
    vi <- which(m2 == min(m2), arr.ind = TRUE) #returns vector with row & column indices of global minimum
    x[i+1,] <-vi #save position
    
  } 
  else{ #go to local minimum
    ri <- sample(4)
    idx <- matrix(x[i,],  nrow=4, ncol=2, byrow = TRUE)+ vec[ri,]
    vi <- which(m2[idx] == min(m2[idx]))
    x[i+1,] <-  c(x[i,1] + vec[ri[vi],1], x[i,2] + vec[ri[vi],2] )
  }
  m <- m1
}

# jitter trajectory (for better visualization)
x <- x + matrix(data=0.5*runif(2*length(x[,1])),ncol=2)

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


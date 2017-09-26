function(sac, t1, t2)
{
  library(ggplot2)
  library(grid)
  
  alpha = 1/30
  
  # definitions
  t <- seq(from=t1,to=t2,by=1)
  L = length(t)
  crate <- rep(0,L)
  
  # compute causal window per time step
  for ( w in 1:L ) {
    tau <- t[w] - sac[,1]
    cw <- alpha^2*tau*exp(-alpha*tau)
    idx <- which(cw>0)
    crate[w] <- sum(cw[idx])
  }
  # divide by number of trials
  Ntrials <- max(sac[,3]) - min(sac[,3])+1
  crate <- crate/Ntrials*1000
  
  # define data frame
  causalrate.df <- data.frame(time=t,causalrate=crate)
  
  # generate plot
  b <- ggplot(causalrate.df) + geom_line(aes(x=time,y=causalrate)) + 
    xlab("Microsaccade onset time  t[ms]") + ylab("Microsaccade rate [1/s]") 
  print(b)
}
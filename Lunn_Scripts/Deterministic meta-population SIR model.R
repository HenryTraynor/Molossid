## R script for deterministic meta-poulation SIR model
## Author: Unknown (2012, edited by Tamika Lunn, SI -> SIR, 2018)
## Sourced from: Welcome Trust "Mathematical Models for Infectious Disease Dynamics" course, Cambridge University (2018)

#Load packages
library(deSolve)

#Define the function to run the model
SIRmat <- function(time, state, pars) { 
  with(as.list(c(state, pars)), 
       { n=pars$n
       beta=pars$beta
       gamma=pars$gamma
       S=state[1:n]
       I=state[(n+1):(2*n)]
       R=state[(2*n+1):(3*n)]
       
       newinfections = colSums(beta*I*S) #sum across the columns of the matrix beta*I*S:
       dS <- -newinfections
       dI <- newinfections - gamma * I
       dR <- gamma * I
       return(list(c(dS, dI, dR)))
       })}

#Define conditions for model
n=5 #number of age groups or populations

#set up the transmission matrix:
d1=5 #within-group B
d2=0.1 #between-group B
B=d1*diag(n) + rbind(cbind(rep(0,n-1),d2*diag(n-1)),rep(0,n)) + rbind(rep(0,n),cbind(d2*diag(n-1),rep(0,n-1))) #matrix with contact between direct neighbours only

#initial conditions:
S=rep(1,n)
I=rep(0,n)
R=rep(0,n)

#start with a small amount of infection in the first group/population:
S[1]=S[1]-1e-6
I[1]=I[1]+1e-6
init <- c(S=S, I=I, R=R)
params <- list(n=n, beta=B, gamma=rep(3,n))
times <- seq(0, 35, by = 0.01)

#run the model:
out <- as.data.frame(ode(y = init, times = times, func= SIRmat, parms = params))
mycolours=rainbow(n)
time=out[,1]
S=out[,2:(n+1)]
I=out[,(n+2):(2*n+1)]
R=out[,(2*n+2):(3*n+1)]
matplot(time,I,type="l",lty=1,col=mycolours,lwd=3)
legend('topright',paste("group",1:n), col=mycolours, pch=19)


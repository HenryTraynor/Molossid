## R script for simple deterministic SIR model
## Author: Hamish McCallum (date unknown, unedited)
## Sourced from: Hamish McCallum (2018)

#Load packages
library(deSolve)

#Define the function to run the model
SIR<-function(t,y,parameters){
  #t is the times for which I want the solution. It is a vector
  #y is a vector with as many values as there are variables in the model
  #parameters is a vector of the parameter values
  S<-y[1]
  I<-y[2]
  R<-y[3]
  
  N<-S+I+R
  a<-parameters[1]
  b<-parameters[2]
  beta<-parameters[3]
  nu<-parameters[4]
  alpha<-parameters[5]
  gamma<-parameters[6]
  
  dy<-numeric(3) #vector called dy, which will contain the three equations
  dy[1]<-a*N-b*S-beta*S*I+nu*R
  dy[2]<-beta*S*I-I*(b+alpha+gamma)
  dy[3]<-gamma*I-R*(b+nu)
  
  return(list(dy))
}

#Define the parameters
pars  <- c(a= 0.06,    # birth rate
           b=0.04,    # death rate
           beta= 0.4,   # transmission rate
           nu= 0.02,    #loss of immunity
           alpha=0.5,  #disease induced death
           gamma= 0.5)     # recovery rate

#Set up the initial values
yini  <- c(S =10, I = 0.1,R=0)
times <- seq(0, 100, by = .1) #output from t=0 to t=100, in steps of 0.1

#Solve
output   <- ode(yini, times, SIR, pars) #needs to be given (in order) initial values of the variables, times for output, the name of the function to solve, the parameters
summary(output) #returns a matrix with columns for time and each of the variables
head(output)

#Visualise
plot(output)
plot(output[,"S"],output[,"I"],type="l")

#Plot all three variables on the one set of axes
output.data<-as.data.frame(output)
par(mfrow=c(1,1),cex=2,lwd=2,las=1)

plot(S~time,data=output.data,type="l",ylim=c(0,10),ylab="N")
lines(I~time,data=output.data,col="red")
lines(R~time,data=output.data,col="blue")

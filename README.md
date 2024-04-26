# Fast-and-Efficient-Bayesian-Inference-for-Epidemic-Models
This repository contains the R code used in:

Fast and Efficient Bayesian Inference for Epidemic Models (2024). Written by Elliot Lloyd. Supervised by Andrew Golightly. Produced as a report for the module 'Project III' within BSc Mathematics (G100) at Durham University.

# R code
The full R scripts with all code including elementary plotting can be found in the labelled files above. The main algorithms used in the report will be described here.

## The deterministic SIR model

The system requires the parameter vector theta=(beta, gamma)', alongside initial conditions such as: total time (T), time step (dt) and the initial population (x0). The reaction network formulas can be expressed in a function:
```
alphaSIR<-function(x,theta)
{
  c(-theta[1]*x[1]*x[2],theta[1]*x[1]*x[2]-theta[2]*x[2])
}
```
which can then be used in the following function to forward simulate the SIR ODE model, for which plotting is straightforward:
```
odesim<-function(T=15,dt=0.0025,x0=c(762,1),theta=c(exp(-6),0.5), afun=alphaSIR)
{
  n<-T/dt
  d<-length(x0)
  x<-matrix(0,nrow=n+1,ncol=d)
  x[1,]<-x0 
  for (i in 2:(n+1)) 
  {
    x[i,]<-x[i-1,]+afun(x[i-1,],theta)*dt
  }
  x
}
```

##The stochastic SIR model


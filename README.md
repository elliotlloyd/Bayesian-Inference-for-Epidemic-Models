# Fast-and-Efficient-Bayesian-Inference-for-Epidemic-Models
This repository contains the R code used in:

Fast and Efficient Bayesian Inference for Epidemic Models (2024). Written by Elliot Lloyd. Supervised by Andrew Golightly. Produced as a report for the module 'Project III' within BSc Mathematics (G100) at Durham University.

# R code
The full R scripts with all code including elementary plotting can be found in the labelled files above. The main algorithms used in the report will be described here.

## The deterministic SIR model

Forward simulation of the SIR ODE model requires the parameter vector theta=(beta, gamma)', alongside initial conditions such as: total time (T), time step (dt) and the initial population (x0). The reaction network formulas can be expressed in a function:
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

## The stochastic SIR model

In order to forward simulate the stochastic SIR model using the Gillespie algorithm and the Poisson leap model, I sampled code from (https://github.com/csgillespie/In-silico-Systems-Biology). This made forward simulation very simple thanks to the *issb* package containing a useful model generator and simulators such as:
```
g = gillespie(model, maxtime=15)
pl = pleap(model, maxtime=15, ddt=0.01)
```
Plots can easily be produced using the matrices returned.

The Euler-Mayurama method is used to forward simulate from the SDE SIR model. It required creating 2 functions:
```

EM_det_ez <- function(x, theta, smat){
  c(-theta[1]*x[1]*x[2],theta[1]*x[1]*x[2]-theta[2]*x[2])
}

EM_stoch_ez <- function(x, theta, smat){
  stoch_m <- matrix(c(theta[1]*x[1]*x[2], -theta[1]*x[1]*x[2],
                        -theta[1]*x[1]*x[2], theta[1]*x[1]*x[2]+theta[2]*x[2]),
                      nrow = 2)
  stoch_sqrt <- t(chol(stoch_m))
  stoch_sqrt
}
```
which calculate the deterministic and stochastic matrices inside the *Euler-Maruyama discretisation formula* which are used in the final function to generate trajectories:
```
sdesim <- function(T, dt, x0, theta, smat, afun1, afun2){
  n <- T/dt
  d <- length(x0)
  x <- matrix(0,nrow=n+1,ncol=d)
  x[1,] <- x0
  for (i in 2:(n+1)) 
  {
    dwt <- c(rnorm(1,0,sqrt(dt)),rnorm(1,0,sqrt(dt)))
    x[i,] <- x[i-1,]+afun1(x[i-1,],theta, smat)*dt + 
      afun2(x[i-1,],theta, smat)%*%dwt
    if(x[i,1] < 0){
      x[i,1] <- 0.0001
    }
    if(x[i,2] < 0){
      x[i,2] <- 0.0001
    }
  }
  x
}
```
Inside `sdesim` the `dwt` variable stores the *Brownian Motion* vector and we ensure that, given the randomness of the model, we do not produce trajectories where compartment counts become negative.

## Bayesian Inference for the ODE model

To generate the *synthetic* data from the specified observation model, normal noise (with standard deviation of 2) was added to the 100 samples from the infective data:
```
Idata <- outODE[1+(0:100)*60,2]
Idata <- rnorm(101,Idata,2)
```
Then given the observation model and the "observed" data points, the parameter MLEs were calculated using the `optim` function. Since `optim` minimises a function passed into it, the `loglikeNeg` function was created to feed the reciprocal of the log-likelihood function into `optim`.
```
loglikeNeg<-function(param, data, tstep, dt, init)
{
  n<-length(data)
  T<-(n-1)*dt
  incr<-dt/tstep
  out<-odesim(T, dt=tstep, x0=init, theta=exp(param), afun=alphaSIR)
  out<-out[1+(0:(n-1))*incr,2] 
  -1*sum(dnorm(data, out, 2, log=TRUE)) 
}

mle <- optim(log(c(exp(-6),0.5)),loglikeNeg,data=Idata,tstep=0.0025,dt=0.15,init=c(762,1))
exp(mle$par)
```
### RWM MCMC sampler 

To perform the Metropolis-Hastings algorithm with a symmetric random walk proposal mechanism 

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

## Bayesian inference for the ODE model

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

To perform the Metropolis-Hastings algorithm with a symmetric random walk proposal mechanism, we first needed functions to evaluate the prior and posterior densities:
```
loglike <- function(param, data, tstep, deltat, init){
  n <- length(data)
  T <- (n-1)*deltat 
  incr <- deltat/tstep
  out <- odesim(T=15, dt=tstep, x0=init, theta=exp(param),afun=alphaSIR)
  out <- out[1+(0:(n-1))*incr,2] 
  sum(dnorm(data, out, exp(param[3]), log=TRUE))
}

logprior <- function(psi){
  sum(dnorm(psi, log = TRUE))
}

lpost <- function(psi){
  lprior <- logprior(psi)
  llike <- loglike(c(psi[1], psi[2], 2), data=Idata, tstep=0.025, deltat=0.15, init=c(762,1)) 
  return(lprior+llike)
}
```
These are necessary to accurately calculate the acceptance probability within the sampler.
```
rwm <- function(N, init, Sigma){
  theta_mat <- matrix(0, nrow = N, ncol = 2)
  theta_mat[1,] <- init
  psi <- init
  count <- 0
  for(i in 2:N){
    can <- mvrnorm(1, psi, Sigma)
    laprob <- lpost(can) - lpost(psi)
    if (log(runif(1)) < laprob){
      psi <- can
      count <- count + 1
    }
    theta_mat[i,] <- psi
  }
  print(count/(N-1))
  return(theta_mat)
}
```
The sampler outputs the logged posterior matrix for theta and also prints the empirical acceptance rate to quickly check efficiency. The sampler was tuned using the optimal proposal variance:
```
SigmaOpt <- (2.38^2)/2*var(out)
```
ESS values were obtained using the `effectiveSize` function found in the `coda` package and the time taken could also be determined using the `system.time` function.

## ABC techniques for the ODE model

Summary statistics are required in order to compare *pseudo-data* to the true data within the ABC samplers. These were chosen to be the mode and the time at which the total number of infectives falls below a certain threshold, calculated by `stop_time`.
```
stop_time <- function(SIR_vec, thresh){
  time <- which(SIR_vec < thresh)
  n <- length(time)
  time_after <- rep(0,n)
  for(i in 1:n){
    if(time[i] > which.max(SIR_vec)){
      time_after[i] <- time[i]
    }
  }
  time_after[time_after == 0] <- NA
  time_after_2 <- time_after[complete.cases(time_after)]
  if(length(time_after_2) == 0){
    end_time <- length(SIR_vec)
  } else {
    end_time <- time_after_2[1]
  }
  end_time
}

summary_stats_2 <- c(which.max(Idata_ode), stop_time(Idata_ode, 10))
```
The ABC samplers here use the Euclidean distance metric to compare pseudo-data to the observed data.
```
euclidean <- function(a,b){
  distance <- sqrt(sum((a - b)^2))
  distance
}
```
The ABC rejection sampler makes use of a `while` loop to ensure we obtained a desired number of approximate posterior draws. It also proposes values from the new prior distributions established.
```
abc_norm_2 <- function(samples, epsilon, v1, v2){
  posterior <- matrix(0, nrow = samples, ncol = 2)
  accepted_samples <- 1
  while(accepted_samples < samples+1){
    theta_prop <- exp(mvrnorm(1, c(-6, log(0.5)), diag(c(v1,v2))))
    pseudo_data <- odesim(theta=theta_prop)
    pseudo_Idata <- pseudo_data[1+(0:100)*60,2]
    pseudo_Idata_obs <- pseudo_Idata + rnorm(101,0,2)
    pseudo_sum <- c(which.max(pseudo_Idata_obs), stop_time(pseudo_Idata_obs, 10))
    distance <- euclidean(summary_stats_2, pseudo_sum)
    if (distance < epsilon){
      posterior[accepted_samples,] <- theta_prop
      accepted_samples <- accepted_samples + 1
    }
  }
  posterior
}
```
When creating an ABC-MCMC function, it was easier to create an ABC acceptance function which could be used inside the main sampler to make it less bulky.
```
abc_acceptance <- function(can){
  pseudo_data <- odesim(theta = exp(can))
  pseudo_Idata <- pseudo_data[1+(0:100)*60,2]
  pseudo_Idata_obs <- pseudo_Idata + rnorm(101,0,2)
  pseudo_sum <- c(which.max(pseudo_Idata_obs), stop_time(pseudo_Idata_obs, 10))
  distance <- euclidean(summary_stats, pseudo_sum)
  distance
}
```
This simplified the ABC-MCMC sampler which performs the ABC step after the Metropolis-Hastings step.
```
abc_mcmc <- function(N, init, Sigma, epsilon){
  theta_mat <- matrix(0, nrow = N, ncol = 2)
  theta_mat[1,] <- init
  psi <- init
  accepted <- 0
  for(i in 2:N){
    can <- mvrnorm(1, psi, Sigma)
    laprob <- logprior(can)-logprior(psi)
    if (abc_acceptance(can) < epsilon){
      if (log(runif(1)) < laprob){
        psi <- can
        theta_mat[i,] <- psi
        accepted <- accepted + 1
      } else{
        theta_mat[i,] <- psi
      }
    } else{
      theta_mat[i,] <- psi
    }
  }
  print(accepted/(N-1))
  return(theta_mat)
}
```
## ABC techniques for stochastic models

Both the ABC rejection and MCMC samplers in a stochastic framwork were largely similar to those used in the deterministic setting. The only changes that were needed was the way in which the pseudo-data was generated. For the Poisson leap model, we used
```
model = create_model(smat, h, initial, theta_prop, get_f)
pseudo_data <- pleap(model, maxtime=15, ddt=0.0025)
```
within the ABC rejection step. Similarly for the SDE model, the pseudo-data was generated by:
```
pseudo_data <- sdesim(T = 15, dt=0.0025, x0=c(762,1),theta=theta_prop,
                          smat=matrix(c(-1,0,1,-1),nrow=2),afun1=EM_det_ez,afun2=EM_stoch_ez)
```
Point estimates and 95% credible intervals could be easily calculated after generating approximate posteriors.
```
c(mean(post_mat_opt_sde_mcmc2[,1]), mean(post_mat_opt_sde_mcmc2[,2]))
sde_conf_intervals <- apply(post_mat_opt_sde_mcmc2, 2, function(column) {
  c(lower = quantile(column, probs = 0.025), 
    upper = quantile(column, probs = 0.975))
})
```
Slight modifications of this code were used on each of the other samplers to obtain the same statistics.

  

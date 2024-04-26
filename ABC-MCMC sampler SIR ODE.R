##ABC-MCMC sampler for ODE model
library(coda)
library(MASS)

#generate data
alphaSIR<-function(x,theta)
{
  c(-theta[1]*x[1]*x[2],theta[1]*x[1]*x[2]-theta[2]*x[2])
}

odesim<-function(T = 15, dt=0.0025, x0=c(762,1),theta, afun=alphaSIR)
{
  n <- T/dt
  d <- length(x0)
  x <- matrix(0,nrow=n+1,ncol=d)
  x[1,] <- x0 
  for (i in 2:(n+1)) 
  {
    x[i,] <- x[i-1,]+afun(x[i-1,],theta)*dt 
  }
  x
}

outODE <- odesim(T = 15, dt=0.0025, x0=c(762,1),theta=c(exp(-6),0.5),afun=alphaSIR)

Idata_ode2 <- outODE[1+(0:100)*60,2]
Idata_ode2 <- rnorm(101,Idata_ode2,2)
plot(ts(outODE[,2],start=0,deltat=0.0025))
lines(ts(Idata_ode2,start=0,deltat=0.15),type="p")

#set up prior pdf evaluation function for acceptance probability
logprior <- function(psi){
  sum(dnorm(psi, log = TRUE))
}

#the sampler itself
euclidean <- function(a,b){
  distance <- sqrt(sum((a - b)^2))
  distance
}

stop_time <- function(SIR_vec, thresh){
  peak_time <- which.max(SIR_vec)
  post_peak_vec <- SIR_vec[(peak_time+1):length(SIR_vec)]
  below_thresh_time <- which(post_peak_vec < thresh)
  if (length(below_thresh_time) > 0) {
    return((peak_time + below_thresh_time[1]))
  } else {
    return(length(SIR_vec))
  }
}

summary_stats <- c(which.max(Idata_ode2), stop_time(Idata_ode2, 10))

#ABC acceptance function
abc_acceptance <- function(can){
  pseudo_data <- odesim(theta = exp(can))
  pseudo_Idata <- pseudo_data[1+(0:100)*60,2]
  pseudo_Idata_obs <- pseudo_Idata + rnorm(101,0,2)
  pseudo_sum <- c(which.max(pseudo_Idata_obs), stop_time(pseudo_Idata_obs, 10))
  distance <- euclidean(summary_stats, pseudo_sum)
  distance
}

#ABC-MCMC sampler
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

out <- abc_mcmc(1000, init=c(-6, log(0.5)), Sigma = diag(rep(0.000033,2)), 2)
par(mfrow = c(2,1))
plot(ts(exp(out[100:1000,1])), xlab = "Iteration", ylab = "beta")
plot(ts(exp(out[100:1000,2])), xlab = "Iteration", ylab = "gamma")
par(mfrow = c(1,1))

SigmaOpt <- (2.38^2)/2*var(out)

outopt <- abc_mcmc(10000, init=c(-6, log(0.5)), SigmaOpt, 2)
effectiveSize(outopt)
cpu_time_3 <- system.time(abc_mcmc(10000, init=c(-6, log(0.5)), SigmaOpt, 2)) #139.81
plot(ts(outopt))

post_mat <- exp(outopt)
effectiveSize(post_mat) #535.1404 and 243.6046
par(mfrow = c(2,3))
plot(ts(post_mat[100:10000,1]), xlab = "iteration", ylab = "value", main = "beta")
acf(post_mat[100:10000,1], main = "")
hist(post_mat[100:10000,1], freq = F, xlab = "value", main = "beta", ylim = c(0,70000))
lines(density(post_mat[100:10000,1]), lwd = 2, col = "red")
plot(ts(post_mat[100:10000,2]), xlab = "iteration", ylab = "value", main = "gamma")
acf(post_mat[100:10000,2], main = "")
hist(post_mat[100:10000,2], freq = F, xlab = "value", main = "gamma", ylim = c(0,150))
lines(density(post_mat[100:10000,2]), lwd = 2, col = "red")
par(mfrow = c(1,1))

#posterior predictive
N <- 10000
T <- 15
dt <- 0.0025
I.mat_ode_mcmc <- matrix(0,nrow=T/dt+1,ncol=N)
for(i in 1:N)
{
  I.mat_ode_mcmc[,i] <- odesim(T=15, dt=0.0025, x0=c(762,1), theta=post_mat[i,], afun=alphaSIR)[,2]
  I.mat_ode_mcmc[,i] <- rnorm((T/dt)+1, I.mat_ode_mcmc[,i],2)
}

I.mean_ode_mcmc <- apply(I.mat_ode_mcmc,1,mean)
I.lq_ode_mcmc <- apply(I.mat_ode_mcmc,1,quantile,0.025)
I.uq_ode_mcmc <- apply(I.mat_ode_mcmc,1,quantile,0.975)

par(mfrow = c(1,1))
plot(ts(I.mean_ode_mcmc,start=0,deltat=0.0025),ylim=c(0,400), xlab = "Time(Days)", ylab = "It")
lines(ts(I.lq_ode_mcmc,start=0,deltat=0.0025),col=2)
lines(ts(I.uq_ode_mcmc,start=0,deltat=0.0025),col=2)
lines(ts(Idata_ode2,start=0,deltat=0.15),type="p")

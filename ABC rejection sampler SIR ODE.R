library(MASS)
library(dplyr)

##get synthetic data
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

Idata_ode <- outODE[1+(0:100)*60,2]
#Add noise!
Idata_ode <- rnorm(101,Idata_ode,2) #sigma=2
#Plot
plot(ts(outODE[,2],start=0,deltat=0.0025))
lines(ts(Idata_ode,start=0,deltat=0.15),type="p")

##define distance metric
euclidean <- function(a,b){
  distance <- sqrt(sum((a - b)^2))
  distance
}

##calculate summary vector
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

##rejection sampler for normal priors
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

post_normal_2 <- abc_norm_2(100, 2, 0.02, 0.05)
cpu_time <- system.time(abc_norm_2(100,2,0.02,0.05)) #119.15
par(mfrow = c(1,2))
hist(post_normal_2[,1], freq = F, xlab = "Value", main = "beta", 
     ylim = c(0,60000))
lines(density(post_normal_2[,1]), lwd = 2, col = 2)
lines(seq(0.002,0.003,0.00001), dlnorm(seq(0.002,0.003,0.00001), -6, 0.02), 
      lwd = 2, col = 4)
hist(post_normal_2[,2], freq = F, xlab = "Value", main = "gamma", 
     ylim = c(0,125))
lines(density(post_normal_2[,2]), lwd = 2, col = 2)
lines(seq(0.4,0.6,0.001), dlnorm(seq(0.4,0.6,0.001), log(0.5), 0.05), 
      lwd = 2, col = 4)
par(mfrow = c(1,1))

#posterior predictive
N <- 100
T <- 15
dt <- 0.0025
I.mat_ode_rej <- matrix(0,nrow=T/dt+1,ncol=N)
for(i in 1:N)
{
  I.mat_ode_rej[,i] <- odesim(T=15, dt=0.0025, x0=c(762,1), theta=post_normal_2[i,], afun=alphaSIR)[,2]
  I.mat_ode_rej[,i] <- rnorm((T/dt)+1, I.mat_ode_rej[,i],2)
}

I.mean_ode_rej <- apply(I.mat_ode_rej,1,mean)
I.lq_ode_rej <- apply(I.mat_ode_rej,1,quantile,0.025)
I.uq_ode_rej <- apply(I.mat_ode_rej,1,quantile,0.975)

par(mfrow = c(1,2))
plot(ts(I.mean_ode_rej,start=0,deltat=0.0025),ylim=c(0,400), xlab = "Time(Days)", ylab = "It")
lines(ts(I.lq_ode_rej,start=0,deltat=0.0025),col=2)
lines(ts(I.uq_ode_rej,start=0,deltat=0.0025),col=2)
lines(ts(Idata_ode,start=0,deltat=0.15),type="p")

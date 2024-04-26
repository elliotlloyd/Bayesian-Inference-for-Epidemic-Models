#ABC on SDE output
library(MASS)
library(coda)

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

sdesim <- function(T, dt, x0, theta, smat, afun1, afun2){
  n <- T/dt
  d <- length(x0)
  x <- matrix(0,nrow=n+1,ncol=d)
  x[1,] <- x0
  dwt <- c(rnorm(1,0,dt),rnorm(1,0,dt))
  for (i in 2:(n+1)) 
  {
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

set.seed(76)

outSDE <- sdesim(T = 15, dt=0.0025, x0=c(762,1),theta=c(exp(-6),0.5),
                 smat=matrix(c(-1,0,1,-1),nrow=2),afun1=EM_det_ez,afun2=EM_stoch_ez)

Idata_SDE <- outSDE[1+(0:100)*60,2]
Idata_SDE_obs <- rnorm(101,Idata_SDE,2)
plot(ts(outSDE[,2],start=0,deltat=0.01),xlab="time",ylab="It",main="")
lines(ts(Idata_SDE_obs,start=0,deltat=0.6),type="p")

#sampler setup
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

summary_stats_sde <- c(which.max(Idata_SDE_obs), stop_time(Idata_SDE_obs, 10))

abc_acceptance_sde <- function(can){
  pseudo_data <- sdesim(T = 15, dt=0.0025, x0=c(762,1),theta=exp(can),
                        smat=matrix(c(-1,0,1,-1),nrow=2),afun1=EM_det_ez,afun2=EM_stoch_ez)
  pseudo_Idata <- pseudo_data[1+(0:100)*60,2]
  pseudo_Idata_obs <- pseudo_Idata + rnorm(101,0,2)
  pseudo_sum <- c(which.max(pseudo_Idata_obs), stop_time(pseudo_Idata_obs, 10))
  distance <- euclidean(summary_stats_sde, pseudo_sum)
  distance
}

#rejection sampler
abc_rej_sde <- function(samples, epsilon, v1, v2){
  posterior <- matrix(0, nrow = samples, ncol = 2)
  accepted_samples <- 1
  while(accepted_samples < samples+1){
    theta_prop <- exp(mvrnorm(1, c(-6, log(0.5)), diag(c(v1,v2))))
    pseudo_data <- sdesim(T = 15, dt=0.0025, x0=c(762,1),theta=theta_prop,
                          smat=matrix(c(-1,0,1,-1),nrow=2),afun1=EM_det_ez,afun2=EM_stoch_ez)
    pseudo_Idata <- pseudo_data[1+(0:100)*60,2]
    pseudo_Idata_obs <- pseudo_Idata + rnorm(101,0,2)
    pseudo_sum <- c(which.max(pseudo_Idata_obs), stop_time(pseudo_Idata_obs, 10))
    distance <- euclidean(summary_stats_sde, pseudo_sum)
    if (distance < epsilon){
      posterior[accepted_samples,] <- theta_prop
      accepted_samples <- accepted_samples + 1
    }
  }
  posterior
}

post_normal_sde_rej <- abc_rej_sde(100, 2, 0.02, 0.05)
c(mean(post_normal_sde_rej[,1]), mean(post_normal_sde_rej[,2]))
sde_conf_intervals <- apply(post_normal_sde_rej, 2, function(column) {
  c(lower = quantile(column, probs = 0.025), 
    upper = quantile(column, probs = 0.975))
})
cpu_time_sde_rej <- system.time(abc_rej_sde(100, 2, 0.02, 0.05))
par(mfrow = c(1,2))
hist(post_normal_sde_rej[,1], freq = F, xlab = "Value", main = "Beta",
     xlim = c(0.00175,0.00325), ylim = c(0,1400))
lines(density(post_normal_sde_rej[,1]), lwd = 2, col = 2)
hist(post_normal_sde_rej[,2], freq = F, xlab = "Value", main = "Gamma", 
     xlim = c(0.3, 0.8))
lines(density(post_normal_sde_rej[,2]), lwd = 2, col = 2)
par(mfrow = c(1,1))

##updated prior distribution
logpriornew <- function(psi) {
  log_density_beta <- dnorm(psi[1], exp(-6), 0.02, log = TRUE)
  log_density_gamma <- dnorm(psi[2], 0.5, 0.05, log = TRUE)
  return(log_density_beta + log_density_gamma)
}

abc_mcmc_sde_new <- function(N, init, Sigma, epsilon){
  theta_mat <- matrix(0, nrow = N, ncol = 2)
  theta_mat[1,] <- init
  psi <- init
  accepted <- 0
  for(i in 2:N){
    can <- mvrnorm(1, psi, Sigma)
    laprob <- logpriornew(exp(can))-logpriornew(exp(psi))
    if (abc_acceptance_sde(can) < epsilon){
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

post_normal_sde_mcmc2 <- abc_mcmc_sde_new(10000, init=c(-6, log(0.5)), 
                                     Sigma = diag(c(0.0002, 0.0005)), 3)

par(mfrow = c(2,1))
plot(ts(exp(post_normal_sde_mcmc2[1000:10000,1])), xlab = "Iteration", ylab = "beta")
plot(ts(exp(post_normal_sde_mcmc2[1000:10000,2])), xlab = "Iteration", ylab = "gamma")
par(mfrow = c(1,1))

SigmaOpt_sde2 <- (2.38^2)/2*var(post_normal_sde_mcmc2)

post_normal_opt_sde_mcmc2 <- abc_mcmc_sde_new(100000, init=c(-6, log(0.5)), SigmaOpt_sde2, 3)
effectiveSize(post_normal_opt_sde_mcmc2)
cpu_time_sde_mcmc <- system.time(abc_mcmc_sde_new(100000, init=c(-6, log(0.5)), SigmaOpt_sde2, 3))

post_mat_opt_sde_mcmc2 <- exp(post_normal_opt_sde_mcmc2)
effectiveSize(post_mat_opt_sde_mcmc2)
c(mean(post_mat_opt_sde_mcmc2[,1]), mean(post_mat_opt_sde_mcmc2[,2]))
sde_conf_intervals <- apply(post_mat_opt_sde_mcmc2, 2, function(column) {
  c(lower = quantile(column, probs = 0.025), 
    upper = quantile(column, probs = 0.975))
})

par(mfrow = c(2,3))
plot(ts(post_mat_opt_sde_mcmc2[10000:100000,1]), xlab = "iteration", ylab = "value", 
     main = "beta")
acf(post_mat_opt_sde_mcmc2[10000:100000,1], main = "")
hist(post_mat_opt_sde_mcmc2[10000:100000,1], freq = F, xlab = "value", main = "beta", ylim = c(0,1200))
lines(density(post_mat_opt_sde_mcmc2[10000:100000,1]), lwd = 2, col = "red")
plot(ts(post_mat_opt_sde_mcmc2[10000:100000,2]), xlab = "iteration", ylab = "value", 
     main = "gamma")
acf(post_mat_opt_sde_mcmc2[10000:100000,2], main = "")
hist(post_mat_opt_sde_mcmc2[10000:100000,2], freq = F, xlab = "value", main = "gamma")
lines(density(post_mat_opt_sde_mcmc2[10000:100000,2]), lwd = 2, col = "red")
par(mfrow = c(1,1))

#predictive distribution for rejection sampler
N <- 100
T <- 15
dt <- 0.0025
I.mat_sde_rej <- matrix(0,nrow=T/dt+1,ncol=N)
for(i in 1:N)
{
  I.mat_sde_rej[,i] <- sdesim(T = 15, dt=0.0025, x0=c(762,1),theta=post_normal_sde_rej[i,],
                              smat=matrix(c(-1,0,1,-1),nrow=2),afun1=EM_det_ez,afun2=EM_stoch_ez)[,2]
  I.mat_sde_rej[,i] <- rnorm((T/dt)+1, I.mat_sde_rej[,i],2)
}

I.mean_sde_rej <- apply(I.mat_sde_rej,1,mean)
I.median_sde_rej <- apply(I.mat_sde_rej,1,quantile,0.5)
I.lq_sde_rej <- apply(I.mat_sde_rej,1,quantile,0.025)
I.uq_sde_rej <- apply(I.mat_sde_rej,1,quantile,0.975)

par(mfrow = c(1,1))
plot(ts(I.median_sde_rej,start=0,deltat=0.0025),ylim=c(0,400), xlab = "Time(Days)", ylab = "It")
lines(ts(I.lq_sde_rej,start=0,deltat=0.0025),col=2)
lines(ts(I.uq_sde_rej,start=0,deltat=0.0025),col=2)
lines(ts(Idata_SDE_obs,start=0,deltat=0.15),type="p")

#predictive distribution for MCMC sampler
N <- 100000
T <- 15
dt <- 0.0025
I.mat_sde_mcmc <- matrix(0,nrow=T/dt+1,ncol=N)
for(i in 1:N)
{
  I.mat_sde_mcmc[,i] <- sdesim(T = 15, dt=0.0025, x0=c(762,1),theta=post_mat_opt_sde_mcmc2[i,],
                      smat=matrix(c(-1,0,1,-1),nrow=2),afun1=EM_det_ez,afun2=EM_stoch_ez)[,2]
  I.mat_sde_mcmc[,i] <- rnorm((T/dt)+1, I.mat_sde_mcmc[,i],2)
}

I.mean_sde_mcmc <- apply(I.mat_sde_mcmc,1,mean)
I.median_sde_mcmc <- apply(I.mat_sde_mcmc,1,quantile,0.5)
I.lq_sde_mcmc <- apply(I.mat_sde_mcmc,1,quantile,0.025)
I.uq_sde_mcmc <- apply(I.mat_sde_mcmc,1,quantile,0.975)

par(mfrow = c(1,1))
plot(ts(I.median_sde_mcmc,start=0,deltat=0.0025),ylim=c(0,400), xlab = "Time(Days)", ylab = "It")
lines(ts(I.lq_sde_mcmc,start=0,deltat=0.0025),col=2)
lines(ts(I.uq_sde_mcmc,start=0,deltat=0.0025),col=2)
lines(ts(Idata_SDE_obs,start=0,deltat=0.15),type="p")


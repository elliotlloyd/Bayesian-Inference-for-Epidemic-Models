##ABC on pleap output
library(devtools)
library(MASS)
library(issb)
library(coda)

h = function(x, pars) {
  hazs = numeric(length(pars))
  hazs[1] = pars[1]*x[1]*x[2]
  hazs[2] = pars[2]*x[2]
  return(hazs)
}

smat = matrix(0,nrow=2,ncol=2)
smat[1,1] = -1; smat[1,2] = 0
smat[2,1] = 1; smat[2,2] = -1
rownames(smat) = c("S", "I")

initial = c(762, 1)

pars = c(exp(-6), 0.5)

model = create_model(smat, h, initial, pars)

get_f = function(x, pars)
{
  fmat = matrix(0, nrow=2, ncol=2)
  fmat[1,1] = pars[1]
  fmat[1,2] = pars[2]*x[2]
  fmat[2,2] = pars[2]*x[1]
  return(fmat)
}

f=get_f

set.seed(70)

model = create_model(smat, h, initial, pars, get_f)
pl = pleap(model, maxtime=15, ddt=0.0025)

plot(pl[,1], pl[,2], type="l", xlab = "Time", ylab = "Number of individuals")
lines(pl[,1],pl[,3], type="l", col = "red")
par(mfrow = c(1,1))

Idatap <- pl[1+(0:100)*60,3]
Idatap_obs <- rnorm(101,Idatap,2)
plot(ts(pl[,3],start=0,deltat=0.01),xlab="time",ylab="It",main="")
lines(ts(Idatap_obs,start=0,deltat=0.6),type="p")

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

summary_stats_p <- c(which.max(Idatap_obs), stop_time(Idatap_obs, 10))

abc_acceptance_p <- function(can){
  model = create_model(smat, h, initial, exp(can), get_f)
  pseudo_data <- pleap(model, maxtime=15, ddt=0.0025)
  pseudo_Idata <- pseudo_data[1+(0:100)*60,3]
  pseudo_Idata_obs <- pseudo_Idata + rnorm(101,0,2)
  pseudo_sum <- c(which.max(pseudo_Idata_obs), stop_time(pseudo_Idata_obs, 10))
  distance <- euclidean(summary_stats_p, pseudo_sum)
  distance
}

#rejection sampler
abc_rej_p <- function(samples, epsilon, v1, v2){
  posterior <- matrix(0, nrow = samples, ncol = 2)
  accepted_samples <- 1
  while(accepted_samples < samples+1){
    theta_prop <- exp(mvrnorm(1, c(-6, log(0.5)), diag(c(v1,v2))))
    model = create_model(smat, h, initial, theta_prop, get_f)
    pseudo_data <- pleap(model, maxtime=15, ddt=0.0025)
    pseudo_Idata <- pseudo_data[1+(0:100)*60,3]
    pseudo_Idata_obs <- pseudo_Idata + rnorm(101,0,2)
    pseudo_sum <- c(which.max(pseudo_Idata_obs), stop_time(pseudo_Idata_obs, 10))
    distance <- euclidean(summary_stats_p, pseudo_sum)
    if (distance < epsilon){
      posterior[accepted_samples,] <- theta_prop
      accepted_samples <- accepted_samples + 1
    }
  }
  posterior
}

post_normal_p_rej <- abc_rej_p(100, 2, 0.02, 0.05)
c(mean(post_normal_p_rej[,1]), mean(post_normal_p_rej[,2]))
pleap_conf_intervals <- apply(post_normal_p_rej, 2, function(column) {
  c(lower = quantile(column, probs = 0.025), 
    upper = quantile(column, probs = 0.975))
})
cpu_time_p_rej <- system.time(abc_rej_p(100, 2, 0.02, 0.05))
par(mfrow = c(2,2))
hist(post_normal_p_rej[,1], freq = F, xlab = "Value", main = "Beta", ylim=c(0,1500))
lines(density(post_normal_p_rej[,1]), lwd = 2, col = 2)
hist(post_normal_p_rej[,2], freq = F, xlab = "Value", main = "Gamma", 
     xlim=c(0.4,0.65), ylim=c(0,8))
lines(density(post_normal_p_rej[,2]), lwd = 2, col = 2)
par(mfrow = c(1,1))

#updated prior distribution
logpriornew <- function(psi) {
  log_density_beta <- dnorm(psi[1], exp(-6), 0.02, log = TRUE)
  log_density_gamma <- dnorm(psi[2], 0.5, 0.05, log = TRUE)
  return(log_density_beta + log_density_gamma)
}

abc_mcmc_p <- function(N, init, Sigma, epsilon){
  theta_mat <- matrix(0, nrow = N, ncol = 2)
  theta_mat[1,] <- init
  psi <- init
  accepted <- 0
  for(i in 2:N){
    can <- mvrnorm(1, psi, Sigma)
    laprob <- logpriornew(exp(can))-logpriornew(exp(psi))
    if (abc_acceptance_p(can) < epsilon){
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

post_normal_p_mcmc <- abc_mcmc_p(10000, init=c(-6, log(0.5)), 
                                 Sigma = diag(c(0.002, 0.005)), 3)
par(mfrow = c(2,1))
plot(ts(exp(post_normal_p_mcmc[1000:10000,1])), xlab = "Iteration", ylab = "beta")
plot(ts(exp(post_normal_p_mcmc[1000:10000,2])), xlab = "Iteration", ylab = "gamma")
par(mfrow = c(1,1))

SigmaOpt_p <- (2.38^2)/2*var(post_normal_p_mcmc)

post_normal_opt_p_mcmc <- abc_mcmc_p(100000, init=c(-6, log(0.5)), SigmaOpt_p, 3)
effectiveSize(post_normal_opt_p_mcmc)
cpu_time_p_mcmc <- system.time(abc_mcmc_p(100000, init=c(-6, log(0.5)), SigmaOpt_p, 3))

post_mat_opt_p_mcmc <- exp(post_normal_opt_p_mcmc)
effectiveSize(post_mat_opt_p_mcmc)
par(mfrow = c(2,3))
plot(ts(post_mat_opt_p_mcmc[10000:100000,1]), xlab = "iteration", ylab = "value", 
     main = "beta")
acf(post_mat_opt_p_mcmc[10000:100000,1], main = "")
hist(post_mat_opt_p_mcmc[10000:100000,1], freq = F, xlab = "value", main = "beta", ylim=c(0,1400))
lines(density(post_mat_opt_p_mcmc[10000:100000,1]), lwd = 2, col = "red")
plot(ts(post_mat_opt_p_mcmc[10000:100000,2]), xlab = "iteration", ylab = "value", 
     main = "gamma")
acf(post_mat_opt_p_mcmc[10000:100000,2], main = "")
hist(post_mat_opt_p_mcmc[10000:100000,2], freq = F, xlab = "value", main = "gamma",
     ylim = c(0,11))
lines(density(post_mat_opt_p_mcmc[10000:100000,2]), lwd = 2, col = "red")
par(mfrow = c(1,1))

#predictive distribution for rejection sampler
N <- 100
T <- 15
dt <- 0.0025
I.mat_p_rej <- matrix(0,nrow=T/dt+1,ncol=N)
for(i in 1:N)
{
  model = create_model(smat, h, initial, post_normal_p_rej[i,], get_f)
  I.mat_p_rej[,i] <- pleap(model, maxtime=15, ddt=0.0025)[,3]
  I.mat_p_rej[,i] <- rnorm((T/dt)+1, I.mat_p_rej[,i],2)
}

I.mean_p_rej <- apply(I.mat_p_rej,1,mean)
I.median_p_rej <- apply(I.mat_p_rej,1,quantile,0.5)
I.lq_p_rej <- apply(I.mat_p_rej,1,quantile,0.025)
I.uq_p_rej <- apply(I.mat_p_rej,1,quantile,0.975)

par(mfrow = c(1,2))
plot(ts(I.median_p_rej,start=0,deltat=0.0025),ylim=c(0,400), xlab = "Time(Days)", ylab = "It")
lines(ts(I.lq_p_rej,start=0,deltat=0.0025),col=2)
lines(ts(I.uq_p_rej,start=0,deltat=0.0025),col=2)
lines(ts(Idatap_obs,start=0,deltat=0.15),type="p")

#predictive distribution for MCMC sampler
N <- 100000
T <- 15
dt <- 0.0025
I.mat_p_mcmc <- matrix(0,nrow=T/dt+1,ncol=N)
for(i in 1:N)
{
  model = create_model(smat, h, initial, post_mat_opt_p_mcmc[i,], get_f)
  I.mat_p_mcmc[,i] <- pleap(model, maxtime=15, ddt=0.0025)[,3]
  I.mat_p_mcmc[,i] <- rnorm((T/dt)+1, I.mat_p_mcmc[,i],2)
}

I.mean_p_mcmc <- apply(I.mat_p_mcmc,1,mean)
I.median_p_mcmc <- apply(I.mat_p_mcmc,1,quantile,0.5)
I.lq_p_mcmc <- apply(I.mat_p_mcmc,1,quantile,0.025)
I.uq_p_mcmc <- apply(I.mat_p_mcmc,1,quantile,0.975)

par(mfrow = c(1,2))
plot(ts(I.median_p_mcmc,start=0,deltat=0.0025),ylim=c(0,400), xlab = "Time(Days)", ylab = "It")
lines(ts(I.lq_p_mcmc,start=0,deltat=0.0025),col=2)
lines(ts(I.uq_p_mcmc,start=0,deltat=0.0025),col=2)
lines(ts(Idatap_obs,start=0,deltat=0.15),type="p")


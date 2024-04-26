##RWM for SIR ODE model
library(MASS)
library(coda)

alphaSIR<-function(x,theta)
{
  c(-theta[1]*x[1]*x[2],theta[1]*x[1]*x[2]-theta[2]*x[2])
}

odesim<-function(T,dt,x0,theta, afun)
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

Idata <- outODE[1+(0:100)*60,2]
#Add noise!
Idata <- rnorm(101,Idata,2) #sigma=2
#Plot
plot(ts(outODE[,2],start=0,deltat=0.0025), ylab = "Value")
lines(ts(Idata,start=0,deltat=0.15),type="p")

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

out <- rwm(10000, init=c(-6, log(0.5)), Sigma = diag(rep(0.000033,2)))

#trace plots for beta and gamma after 1st run with 100 iteration burn in
par(mfrow = c(2,2))
plot(ts(exp(out[100:10000,1])), xlab = "Iteration", ylab = "value", main = "beta")
plot(ts(exp(out[100:10000,2])), xlab = "Iteration", ylab = "value", main = "gamma")
hist(exp(out[100:10000,1]), freq = F, xlab = "value", main = "")
hist(exp(out[100:10000,2]), freq = F, xlab = "value", main = "")
par(mfrow = c(1,1))

#calculate optimal V
SigmaOpt <- (2.38^2)/2*var(out)

#re-run with optimal V
outopt <- rwm(10000, init=c(-6, log(0.5)), SigmaOpt)
effectiveSize(outopt)
cpu_time_2 <- system.time(rwm(10000, init=c(-6, log(0.5)), SigmaOpt)) #58.52
plot(ts(outopt))

post_mat <- exp(outopt)

par(mfrow = c(2,3))
#trace plots for beta and gamma
plot(ts(post_mat[100:10000,1]), xlab = "iteration", ylab = "value", main = "beta")
plot(ts(post_mat[100:10000,2]), xlab = "iteration", ylab = "value", main = "gamma")

par(mfrow = c(1,3))
#plot marginal posterior for beta
hist(post_mat[100:10000,1], freq = F, xlab = "value", main = "beta")
lines(density(post_mat[100:10000,1]), lwd = 2, col = "red")
abline(v= exp(-6), lwd = 3, col = "blue")
#plot marginal posterior for gamma
hist(post_mat[100:10000,2], freq = F, xlab = "value", main = "gamma")
lines(density(post_mat[100:10000,2]), lwd = 2, col = "red")
abline(v=0.5, lwd = 3, col = "blue")
par(mfrow = c(1,1))

#summaries of beta
mean(post_mat[,1])
sd(post_mat[,1])
quantile(post_mat[,1],c(0.025,0.975))

#adding R0 in
R0 <- post_mat[,1]/post_mat[,2]*763
plot(ts(R0[100:10000]), xlab = "iteration", ylab = "value", main = "R0")
hist(R0[100:10000], freq = F, xlab = "value", main = "R0")
lines(density(R0[100:10000]), lwd = 2, col = "red")
abline(v=((exp(-6))/0.5)*763, lwd = 3, col = "blue")
quantile(R0[100:10000],c(0.025,0.975))
mean(R0)

#predictive for It
N <- 10000
T <- 15
dt <- 0.0025
I.mat <- matrix(0,nrow=T/dt+1,ncol=N)
for(i in 1:N)
{
  I.mat[,i] <- odesim(T=15, dt=0.0025, x0=c(762,1), theta=exp(outopt[i,]), afun=alphaSIR)[,2]
  I.mat[,i] <- rnorm((T/dt)+1, I.mat[,i],2)
}

I.mean <- apply(I.mat,1,mean)
I.lq <- apply(I.mat,1,quantile,0.025)
I.uq <- apply(I.mat,1,quantile,0.975)

par(mfrow = c(1,1))
plot(ts(I.mean,start=0,deltat=0.0025),ylim=c(0,400), xlab = "Time(Days)", ylab = "It")
lines(ts(I.lq,start=0,deltat=0.0025),col=2)
lines(ts(I.uq,start=0,deltat=0.0025),col=2)
lines(ts(Idata,start=0,deltat=0.15),type="p")



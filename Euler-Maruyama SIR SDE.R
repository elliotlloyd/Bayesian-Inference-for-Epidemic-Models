library(matrixStats)

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

outSDE <- sdesim(T = 15, dt=0.0025, x0=c(762,1),theta=c(exp(-6),0.5),
                 smat=matrix(c(-1,0,1,-1),nrow=2),afun1=EM_det_ez, afun2=EM_stoch_ez)
par(mfrow = c(1,1))
plot(ts(outSDE[,1],start=0,deltat=0.0025), ylab = "No. of individuals", main = "")
lines(ts(outSDE[,2],start=0,deltat=0.0025), col = "red")

set.seed(32)
real_plot <- function(n){
  st_mat <- matrix(0, nrow = 6001, ncol = n)
  it_mat <- matrix(0, nrow = 6001, ncol = n)
  for(i in 1:n){
    sde <- sdesim(T = 15, dt=0.0025, x0=c(762,1),theta=c(exp(-6),0.5),
                 smat=matrix(c(-1,0,1,-1),nrow=2),afun1=EM_det_ez, afun2=EM_stoch_ez)
    st_mat[,i] <- sde[,1]
    it_mat[,i] <- sde[,2]
    i <- i+1
  }
  st_mat_sum <- matrix(0, nrow = 6001, ncol = 3)
  it_mat_sum <- matrix(0, nrow = 6001, ncol = 3)
  st_mat_sum[,1] <- apply(st_mat, 1, quantile, probs = c(0.25))
  st_mat_sum[,2] <- apply(st_mat, 1, quantile, probs = c(0.5))
  st_mat_sum[,3] <- apply(st_mat, 1, quantile, probs = c(0.75))
  it_mat_sum[,1] <- apply(it_mat, 1, quantile, probs = c(0.25))
  it_mat_sum[,2] <- apply(it_mat, 1, quantile, probs = c(0.5))
  it_mat_sum[,3] <- apply(it_mat, 1, quantile, probs = c(0.75))
  
  par(mfrow=c(1,2))
  plot(seq(0,15, by=0.0025), st_mat[,1], type="l", xlab = "Time", ylab = "",main = "St")
  for(i in 2:n){
    lines(seq(0,15, by=0.0025), st_mat[,i], type="l", col = i)
    i <- i+1
  }
  plot(seq(0,15, by=0.0025), it_mat[,1], type="l", xlab = "Time", ylab = "",main = "It", ylim = c(0, max(it_mat)))
  for(i in 2:n){
    lines(seq(0,15, by=0.0025), it_mat[,i], type="l", col = i)
    i <- i+1
  }
  par(mfrow = c(1,1))
  
  par(mfrow = c(1,2))
  plot(seq(0,15, by=0.0025), st_mat_sum[,1], t="l", 
       xlab = "Time", ylab = "", main = "St", col = 2)
  lines(seq(0,15, by=0.0025), st_mat_sum[,2], lwd = 2)
  lines(seq(0,15, by=0.0025), st_mat_sum[,3], col = 2)
  plot(seq(0,15, by=0.0025), it_mat_sum[,1], t="l", ylim = range(it_mat_sum[,3]),
       xlab = "Time", ylab = "", main = "It", col = 2)
  lines(seq(0,15, by=0.0025), it_mat_sum[,2], lwd = 2)
  lines(seq(0,15, by=0.0025), it_mat_sum[,3], col = 2)
  par(mfrow = c(1,1))
  
  st_mat_sum2 <- matrix(0, nrow = 6001, ncol = 3)
  it_mat_sum2 <- matrix(0, nrow = 6001, ncol = 3)
  st_mat_sum2[,1] <- apply(st_mat, 1, quantile, probs = c(0.025))
  st_mat_sum2[,2] <- apply(st_mat, 1, mean)
  st_mat_sum2[,3] <- apply(st_mat, 1, quantile, probs = c(0.975))
  it_mat_sum2[,1] <- apply(it_mat, 1, quantile, probs = c(0.025))
  it_mat_sum2[,2] <- apply(it_mat, 1, mean)
  it_mat_sum2[,3] <- apply(it_mat, 1, quantile, probs = c(0.975))
  
  par(mfrow = c(1,2))
  plot(seq(0,15, by=0.0025), st_mat_sum2[,1], t="l", 
       xlab = "Time", ylab = "", main = "St", col = 2)
  lines(seq(0,15, by=0.0025), st_mat_sum2[,2], lwd = 2)
  lines(seq(0,15, by=0.0025), st_mat_sum2[,3], col = 2)
  plot(seq(0,15, by=0.0025), it_mat_sum2[,1], t="l", ylim = range(it_mat_sum2[,3]),
       xlab = "Time", ylab = "", main = "It", col = 2)
  lines(seq(0,15, by=0.0025), it_mat_sum2[,2], lwd = 2)
  lines(seq(0,15, by=0.0025), it_mat_sum2[,3], col = 2)
  par(mfrow = c(1,1))
}

real_plot(5000)

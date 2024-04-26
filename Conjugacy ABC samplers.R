##Easy conjugacy code to show ABC works

euclidean <- function(z,x){
  sqrt(sum((z-x)^2))
}

abc_rejection <- function(runs, xbar, n, alpha, beta, epsilon){
  theta_post <- c()
  for (i in 1:runs){
    theta_prop <- rbeta(1, alpha, beta)
    pseudo_data <- rbinom(n, 1, theta_prop)
    pseudo_sum <- mean(pseudo_data)
    distance <- euclidean(pseudo_sum, xbar)
    if (distance <= epsilon){
      theta_post <- append(theta_post, theta_prop)
    }
  }
  print(length(theta_post)/runs)
  theta_post
}

out <- abc_rejection(10000, 0.25, 50, 2, 5, epsilon = 1)
out2 <- abc_rejection(10000, 0.25, 50, 2, 5, epsilon = 0.5)
out3 <- abc_rejection(10000, 0.25, 50, 2, 5, epsilon = 0.1)
out4 <- abc_rejection(10000, 0.25, 50, 2, 5, epsilon = 0.025)
par(mfrow = c(2,2))
hist(out, freq = F, xlab = "value", main = "epsilon = 1", ylim = c(0,7))
lines(seq(0,1,0.005), dbeta(seq(0,1,0.005),2,5), col = 2, lwd = 2)
lines(seq(0,1,0.005), dbeta(seq(0,1,0.005), 2+50*0.25, 5+50-50*0.25), col = 3, lwd = 2)
hist(out2, freq = F, xlab = "value", main = "epsilon = 0.5", ylim = c(0,7))
lines(seq(0,1,0.005), dbeta(seq(0,1,0.005),2,5), col = 2, lwd = 2)
lines(seq(0,1,0.005), dbeta(seq(0,1,0.005), 2+50*0.25, 5+50-50*0.25), col = 3, lwd = 2)
hist(out3, freq = F, xlab = "value", main = "epsilon = 0.1", ylim = c(0,7))
lines(seq(0,1,0.005), dbeta(seq(0,1,0.005),2,5), col = 2, lwd = 2)
lines(seq(0,1,0.005), dbeta(seq(0,1,0.005), 2+50*0.25, 5+50-50*0.25), col = 3, lwd = 2)
hist(out4, freq = F, xlab = "value", main = "epsilon = 0.025", ylim = c(0,7))
lines(seq(0,1,0.005), dbeta(seq(0,1,0.005),2,5), col = 2, lwd = 2)
lines(seq(0,1,0.005), dbeta(seq(0,1,0.005), 2+50*0.25, 5+50-50*0.25), col = 3, lwd = 2)
par(mfrow = c(1,1))

abc_mcmc <- function(N, v, xbar, n, alpha, beta, epsilon){
  theta_post <- c()
  theta <- alpha/(alpha+beta)
  theta_post[1] <- theta
  accepted <- 0
  for(i in 2:N){
    can <- rnorm(1, theta, v)
    if (can < 0.0001){
      can <- 0.0001
    } else if (can > 0.9999){
      can <- 0.9999
    }
    pseudo_data <- rbinom(n, 1, can)
    pseudo_sum <- mean(pseudo_data)
    distance <- euclidean(pseudo_sum, xbar)
    laprob <- dbeta(can,alpha,beta,log=T)-dbeta(theta,alpha,beta,log=T)
    if (distance <= epsilon){
      if (log(runif(1)) < laprob){
        theta <- can
        theta_post[i] <- theta
        accepted <- accepted + 1
      } else{
        theta_post[i] <- theta
      }
    } else{
      theta_post[i] <- theta
    }
  }
  print(accepted/(N-1))
  return(theta_post)
}

out5 <- abc_mcmc(10000, 0.125, 0.25, 50, 2, 5, epsilon = 0.5)
out6 <- abc_mcmc(10000, 0.125, 0.25, 50, 2, 5, epsilon = 0.1)
optv1 <- (2.38^2)*var(out6)
out_opt6 <- out6 <- abc_mcmc(10000, optv1, 0.25, 50, 2, 5, epsilon = 0.1)
out7 <- abc_mcmc(10000, 0.04, 0.25, 50, 2, 5, epsilon = 0.05)
optv2 <- (2.38^2)*var(out7)
out_opt7 <- abc_mcmc(10000, optv2, 0.25, 50, 2, 5, epsilon = 0.025)

par(mfrow = c(2,3))
plot(ts(out6), main = "", xlab = "iteration", ylab = "value")
acf(out6, main = "")
hist(out6, freq = F, main = "", xlab = "value", ylim = c(0,7))
lines(seq(0,1,0.005), dbeta(seq(0,1,0.005), 2+50*0.25, 5+50-50*0.25), col = 3, lwd = 2)
plot(ts(out7), main = "", xlab = "iteration", ylab = "value")
acf(out7, main = "")
hist(out7, freq = F, main = "", xlab = "value", ylim = c(0,7))
lines(seq(0,1,0.005), dbeta(seq(0,1,0.005), 2+50*0.25, 5+50-50*0.25), col = 3, lwd = 2)
par(mfrow = c(1,1))

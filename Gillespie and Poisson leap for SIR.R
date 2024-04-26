library(devtools)
install_github("csgillespie/In-silico-Systems-Biology", subdir = "issb", force = TRUE)
library(MASS)
library(issb)


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

set.seed(23) #for reproducability

model = create_model(smat, h, initial, pars, get_f)
g = gillespie(model, maxtime=15)
pl = pleap(model, maxtime=15, ddt=0.01)

plot(l1[,1], l1[,2], type="l", xlab = "Time", ylab = "Number of individuals")
lines(l1[,1],l1[,3], type="l", col = "red")

par(mfrow=c(1,2))
plot(g[,1], g[,2], type="l", xlab = "Time", ylab = "Number of individuals")
lines(g[,1],g[,3], type="l", col = "red")
plot(pl[,1], pl[,2], type="l", xlab = "Time", ylab = "Number of individuals")
lines(pl[,1],pl[,3], type="l", col = "red")
par(mfrow = c(1,1))

real_plot <- function(n){
  st_mat <- matrix(0, nrow = 1501, ncol = n)
  it_mat <- matrix(0, nrow = 1501, ncol = n)
  for(i in 1:n){
    pl = pleap(model, maxtime=15, ddt=0.01)
    st_mat[,i] <- pl[,2]
    it_mat[,i] <- pl[,3]
    i <- i + 1
  }
  st_rid <- c()
  it_rid <- c()
  for(i in 1:n){
    if(min(st_mat[,i]) == 762){
      st_rid <- append(st_rid, i)
    }
    if(max(it_mat[,i]) == 1){
      it_rid <- append(it_rid, i)
    }
    i <- i + 1
  }
  st_mat_2 <- st_mat[,-st_rid]
  it_mat_2 <- it_mat[,-it_rid]
  
  st_mat_sum <- matrix(0, nrow = 1501, ncol = 3)
  it_mat_sum <- matrix(0, nrow = 1501, ncol = 3)
  st_mat_sum[,1] <- apply(st_mat_2, 1, quantile, probs = 0.25)
  st_mat_sum[,2] <- apply(st_mat_2, 1, quantile, probs = 0.5)
  st_mat_sum[,3] <- apply(st_mat_2, 1, quantile, probs = 0.75)
  it_mat_sum[,1] <- apply(it_mat_2, 1, quantile, probs = 0.25)
  it_mat_sum[,2] <- apply(it_mat_2, 1, quantile, probs = 0.5)
  it_mat_sum[,3] <- apply(it_mat_2, 1, quantile, probs = 0.75)
  
  par(mfrow=c(1,2))
  plot(seq(0,15, by=0.01), st_mat[,1], type="l", xlab = "Time", ylab = "",main = "St")
  for(i in 2:n){
    lines(seq(0,15, by=0.01), st_mat[,i], type="l", col = i)
    i <- i+1
  }
  plot(seq(0,15, by=0.01), it_mat[,1], type="l", xlab = "Time", ylab = "",main = "It", ylim = c(0,max(it_mat)))
  for(i in 2:n){
    lines(seq(0,15, by=0.01), it_mat[,i], type="l", col = i)
    i <- i+1
  }
  plot(seq(0,15, by=0.01), st_mat_sum[,1], t="l", col = 2, xlab = "Time", 
       ylab = "", main = "St")
  lines(seq(0,15, by=0.01), st_mat_sum[,2], lwd=2)
  lines(seq(0,15, by=0.01), st_mat_sum[,3], col=2)
  plot(seq(0,15, by=0.01), it_mat_sum[,1], t="l", col = 2, ylim = range(it_mat_sum[,3]), 
       xlab = "Time", ylab = "", main = "It")
  lines(seq(0,15, by=0.01), it_mat_sum[,2], lwd=2)
  lines(seq(0,15, by=0.01), it_mat_sum[,3], col=2)
  par(mfrow = c(1,1))
}

real_plot(5000)

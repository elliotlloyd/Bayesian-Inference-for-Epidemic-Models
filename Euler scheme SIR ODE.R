alphaSIR<-function(x,theta)
{
  c(-theta[1]*x[1]*x[2],theta[1]*x[1]*x[2]-theta[2]*x[2])
}

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

out<-odesim()
N<-sum(out[1,])
par(mfrow = c(1,1))
plot(ts(out[,1],start=0,deltat=0.0025), ylab = "Number of individuals")
lines(ts(out[,2], start=0, deltat=0.0025), col = "red")
lines(ts(N-(out[,1]+out[,2]),start=0,deltat=0.0025), col = "green")



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








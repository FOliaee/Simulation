library(GoFKernel)
mu.x0=mu.t0=10 
s0=c(1,2,5)
deltax=seq(1.1,2,.1)
deltat=seq(.5,.95,.05)
ATSx=c()
ATSt=c()
ATSxt=matrix(0, nrow=10, ncol=10)

#The value of parameter " a " for weibull distribution

a0=c()
for(i in 1:3){
  f=function(a){
    (gamma((2/a)+1)/((gamma((1/a)+1))^2))-(s0[i]/(mu.t0))^2-1
  }
  a0[i]=uniroot(f,lower = 0.9999 , upper = 1e10 , extendInt = "yes" )$root
}

a1=matrix(0, nrow=3, ncol=10)
for(i in 1:3){
  for(j in 1:10){
    f=function(a){
      (gamma((2/a)+1)/((gamma((1/a)+1))^2))-(s0[i]/(deltat[j]*mu.t0))^2-1
    }
    
    a1[i,j]=uniroot(f,lower =.999, upper =1e10, extendInt = "yes" )$root
  }
}


a1t=matrix(0, nrow=3, ncol=10)
for(i in 1:3){
  for(j in 1:10){
    f=function(a){
      (gamma((2/a)+1)/((gamma((1/a)+1))^2))-(s0[i]/(deltat[j]*mu.t0))^2-1
    }
    
    a1t[i,j]=uniroot(f,lower =.999, upper =1e10, extendInt = "yes" )$root
  }
}

a1x=matrix(0, nrow=3, ncol=10)
for(i in 1:3){
  for(l in 1:10){
    f=function(a){
      (gamma((2/a)+1)/((gamma((1/a)+1))^2))-(s0[i]/(deltax[l]*mu.x0))^2-1
    }
    
    a1x[i,l]=uniroot(f,lower =.999, upper =1e10, extendInt = "yes" )$root
  }
}
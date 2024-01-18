#*****UCL values for Statistic Z1
#gamma-gamma
ucl.gg1=matrix(0, nrow=3, ncol=3, dimnames = list(c("1","2","5"),c("1","2","5")))

for (i in 1:3){
  for (j in 1:3) {  
    F.Z1=function(z){
      func=function(x,mu.t0,mu.x0){
        pgamma((x-z)*mu.t0,shape=(mu.t0/s0[i])^2,scale=(s0[i]^2/mu.t0))*dgamma(x*mu.x0,shape=(mu.x0/s0[j])^2,scale=(s0[j]^2/mu.x0))
      }
      1-mu.x0*integrate(func,0,Inf,mu.t0=mu.t0,mu.x0=mu.x0)$value
    }
    F.inv=inverse(F.Z1,lower=0,upper=100)
    ucl.gg1[i,j]=F.inv(1-0.027)
  }   
} 

#gamma-lognormal
ucl.gl1=matrix(0, nrow=3, ncol=3, dimnames = list(c("1","2","5"),c("1","2","5")))

for (i in 1:3){
  for (j in 1:3) {  
    F.Z1=function(z){
      func=function(x,mu.t0,mu.x0){
        pgamma((x-z)*mu.t0,shape=(mu.t0/s0[i])^2,scale=(s0[i]^2/mu.t0))*dlnorm(x*mu.x0,meanlog=log(mu.x0^2/sqrt(s0[j]^2+mu.x0^2)),sdlog=sqrt(log(1+(s0[j]^2/mu.x0^2))))
      }
      1-mu.x0*integrate(func,0,Inf,mu.t0=mu.t0,mu.x0=mu.x0)$value
    }
    F.inv=inverse(F.Z1,lower=0,upper=100)
    ucl.gl1[i,j]=F.inv(1-0.027)
  }   
}  

#gamma-normal
ucl.gn1=matrix(0, nrow=3, ncol=2, dimnames = list(c("1","2","5"),c("1","2")))

for (i in 1:3){
  for (j in 1:2) {  
    F.Z1=function(z){
      func=function(x,mu.t0,mu.x0){
        pgamma((x-z)*mu.t0,shape=(mu.t0/s0[i])^2,scale=(s0[i]^2/mu.t0))*dnorm(x*mu.x0,mean=mu.x0,sd=s0[j])
      }
      1-mu.x0*integrate(func,0,Inf,mu.t0=mu.t0,mu.x0=mu.x0)$value
    }
    F.inv=inverse(F.Z1,lower=0,upper=100)
    ucl.gn1[i,j]=F.inv(1-0.027)
  }   
} 

#gamma-weibull
ucl.gw1=matrix(0, nrow=3, ncol=3, dimnames = list(c("1","2","5"),c("1","2","5")))

for (i in 1:3){
  for (j in 1:3) {  
    F.Z1=function(z){
      func=function(x,mu.t0,mu.x0){
        pgamma((x-z)*mu.t0,shape=(mu.t0/s0[i])^2,scale=(s0[i]^2/mu.t0))*dweibull(x*mu.x0,shape=a0[j],scale=mu.x0/gamma((1/a0[j])+1))
      }
      1-mu.x0*integrate(func,0,Inf,mu.t0=mu.t0,mu.x0=mu.x0)$value
    }
    F.inv=inverse(F.Z1,lower=0,upper=100)
    ucl.gw1[i,j]=F.inv(1-0.027)
  }   
} 

#lognormal-gamma
ucl.lg1=matrix(0, nrow=3, ncol=3, dimnames = list(c("1","2","5"),c("1","2","5")))

for (i in 1:3){
  for (j in 1:3) {  
    F.Z1=function(z){
      func=function(x,mu.t0,mu.x0){
        plnorm((x-z)*mu.t0,meanlog=log(mu.t0^2/sqrt(s0[i]^2+mu.t0^2)),sdlog=sqrt(log(1+(s0[i]^2/mu.t0^2))))*dgamma(x*mu.x0,shape=(mu.x0/s0[j])^2,scale=(s0[j]^2/mu.x0))
      }
      1-mu.x0*integrate(func,0,Inf,mu.t0=mu.t0,mu.x0=mu.x0)$value
    }
    F.inv=inverse(F.Z1,lower=0,upper=100)
    ucl.lg1[i,j]=F.inv(1-0.027)
  }   
} 

#lognormal-lognormal
ucl.ll1=matrix(0, nrow=3, ncol=3, dimnames = list(c("1","2","5"),c("1","2","5")))
    
for (i in 1:3){
  for (j in 1:3) {  
    F.Z1=function(z){
      func=function(x,mu.t0,mu.x0){
          plnorm((x-z)*mu.t0,meanlog=log(mu.t0^2/sqrt(s0[i]^2+mu.t0^2)),sdlog=sqrt(log(1+(s0[i]^2/mu.t0^2))))*dlnorm(x*mu.x0,meanlog=log(mu.x0^2/sqrt(s0[j]^2+mu.x0^2)),sdlog=sqrt(log(1+(s0[j]^2/mu.x0^2))))
        }
        1-mu.x0*integrate(func,0,Inf,mu.t0=mu.t0,mu.x0=mu.x0)$value
      }
      F.inv=inverse(F.Z1,lower=0,upper=100)
      ucl.ll1[i,j]=F.inv(1-0.027)
  }   
}  

#lognormal-normal
ucl.ln1=matrix(0, nrow=3, ncol=2, dimnames = list(c("1","2","5"),c("1","2")))

for (i in 1:3){
  for (j in 1:2) {  
    F.Z1=function(z){
      func=function(x,mu.t0,mu.x0){
        plnorm((x-z)*mu.t0,meanlog=log(mu.t0^2/sqrt(s0[i]^2+mu.t0^2)),sdlog=sqrt(log(1+(s0[i]^2/mu.t0^2))))*dnorm(x*mu.x0,mean=mu.x0,sd=s0[j])
      }
      1-mu.x0*integrate(func,0,Inf,mu.t0=mu.t0,mu.x0=mu.x0)$value
    }
    F.inv=inverse(F.Z1,lower=0,upper=100)
    ucl.ln1[i,j]=F.inv(1-0.027)
  }   
}  

#lognormal-weibull
ucl.lw1=matrix(0, nrow=3, ncol=3, dimnames = list(c("1","2","5"),c("1","2","5")))

for (i in 1:3){
  for (j in 1:3) {  
    F.Z1=function(z){
      func=function(x,mu.t0,mu.x0){
        plnorm((x-z)*mu.t0,meanlog=log(mu.t0^2/sqrt(s0[i]^2+mu.t0^2)),sdlog=sqrt(log(1+(s0[i]^2/mu.t0^2))))*dweibull(x*mu.x0,shape=a0[j],scale=mu.x0/gamma((1/a0[j])+1))
      }
      1-mu.x0*integrate(func,0,Inf,mu.t0=mu.t0,mu.x0=mu.x0)$value
    }
    F.inv=inverse(F.Z1,lower=0,upper=100)
    ucl.lw1[i,j]=F.inv(1-0.027)
  }   
}  

#weibull-gamma
ucl.wg1=matrix(0, nrow=3, ncol=3, dimnames = list(c("1","2","5"),c("1","2","5")))

for (i in 1:3){
  for (j in 1:3) {  
    F.Z1=function(z){
      func=function(x,mu.t0,mu.x0){
        pweibull((x-z)*mu.t0,shape=a0[i],scale=mu.t0/gamma((1/a0[i])+1))*dgamma(x*mu.x0,shape=(mu.x0/s0[j])^2,scale=(s0[j]^2/mu.x0))
      }
      1-mu.x0*integrate(func,0,Inf,mu.t0=mu.t0,mu.x0=mu.x0)$value
    }
    F.inv=inverse(F.Z1,lower=0,upper=100)
    ucl.wg1[i,j]=F.inv(1-0.027)
  }   
}      

#weibull-lognormal
ucl.wl1=matrix(0, nrow=3, ncol=3, dimnames = list(c("1","2","5"),c("1","2","5")))

for (i in 1:3){
  for (j in 1:3) {  
    F.Z1=function(z){
      func=function(x,mu.t0,mu.x0){
        pweibull((x-z)*mu.t0,shape=a0[i],scale=mu.t0/gamma((1/a0[i])+1))*dlnorm(x*mu.x0,meanlog=log(mu.x0^2/sqrt(s0[j]^2+mu.x0^2)),sdlog=sqrt(log(1+(s0[j]^2/mu.x0^2))))
      }
      1-mu.x0*integrate(func,0,Inf,mu.t0=mu.t0,mu.x0=mu.x0)$value
    }
    F.inv=inverse(F.Z1,lower=0,upper=100)
    ucl.wl1[i,j]=F.inv(1-0.027)
  }   
}  

#weibull-normal
ucl.wn1=matrix(0, nrow=3, ncol=2, dimnames = list(c("1","2","5"),c("1","2")))

for (i in 1:3){
  for (j in 1:2) {  
    F.Z1=function(z){
      func=function(x,mu.t0,mu.x0){
        pweibull((x-z)*mu.t0,shape=a0[i],scale=mu.t0/gamma((1/a0[i])+1))*dnorm(x*mu.x0,mean=mu.x0,sd=s0[j])
      }
      1-mu.x0*integrate(func,0,Inf,mu.t0=mu.t0,mu.x0=mu.x0)$value
    }
    F.inv=inverse(F.Z1,lower=0,upper=100)
    ucl.wn1[i,j]=F.inv(1-0.027)
  }   
}  

#weibull-weibull
ucl.ww1=matrix(0, nrow=3, ncol=3, dimnames = list(c("1","2","5"),c("1","2","5")))

for (i in 1:3){
  for (j in 1:3) {  
    F.Z1=function(z){
      func=function(x,mu.t0,mu.x0){
        pweibull((x-z)*mu.t0,shape=a0[i],scale=mu.t0/gamma((1/a0[i])+1))*dweibull(x*mu.x0,shape=a0[j],scale=mu.x0/gamma((1/a0[j])+1))
      }
      1-mu.x0*integrate(func,0,Inf,mu.t0=mu.t0,mu.x0=mu.x0)$value
    }
    F.inv=inverse(F.Z1,lower=0,upper=100)
    ucl.ww1[i,j]=F.inv(1-0.027)
  }
}
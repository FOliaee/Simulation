#*****UCL values for Statistic Z3
#gamma-gamma
ucl.gg3=matrix(0, nrow=3, ncol=3, dimnames = list(c("1","2","5"),c("1","2","5")))

for (i in 1:3){
  for (j in 1:3) {  
    F.Z3=function(z){
      func=function(x,mu.t0,mu.x0){
        pgamma((1/(z-x))*mu.t0,shape=(mu.t0/s0[i])^2,scale=(s0[i]^2/mu.t0))*dgamma(x*mu.x0,shape=(mu.x0/s0[j])^2,scale=(s0[j]^2/mu.x0))
      }
      Fx=pgamma(z*mu.x0,shape=(mu.x0/s0[j])^2,scale=(s0[j]^2/mu.x0))
      Fx-mu.x0*integrate(func,0,z,mu.t0=mu.t0,mu.x0=mu.x0)$value
    }
    F.inv=inverse(F.Z3,lower=0,upper=100)
    ucl.gg3[i,j]=F.inv(1-0.027)
  }
}

#gamma-lognormal
ucl.gl3=matrix(0, nrow=3, ncol=3, dimnames = list(c("1","2","5"),c("1","2","5")))

for (i in 1:3){
  for (j in 1:3) {  
    F.Z3=function(z){
      func=function(x,mu.t0,mu.x0){
        pgamma((1/(z-x))*mu.t0,shape=(mu.t0/s0[i])^2,scale=(s0[i]^2/mu.t0))*dlnorm(x*mu.x0,meanlog=log(mu.x0^2/sqrt(s0[j]^2+mu.x0^2)),sdlog=sqrt(log(1+(s0[j]^2/mu.x0^2))))
      }
      Fx=plnorm(z*mu.x0,meanlog=log(mu.x0^2/sqrt(s0[j]^2 + mu.x0^2)),sdlog=sqrt(log(1+(s0[j]^2/mu.x0^2))))
      Fx-mu.x0*integrate(func,0,z,mu.t0=mu.t0,mu.x0=mu.x0)$value
    }
    F.inv=inverse(F.Z3,lower=0,upper=100)
    ucl.gl3[i,j]=F.inv(1-0.027)
  }
}

#gamma-normal
ucl.gn3=matrix(0, nrow=3, ncol=2, dimnames = list(c("1","2","5"),c("1","2")))

for (i in 1:3){
  for (j in 1:2) {  
    F.Z3=function(z){
      func=function(x,mu.t0,mu.x0){
        pgamma((1/(z-x))*mu.t0,shape=(mu.t0/s0[i])^2,scale=(s0[i]^2/mu.t0))*dnorm(x*mu.x0,mean=mu.x0,sd=s0[j])
      }
      Fx=pnorm(z*mu.x0,mean=mu.x0,sd=s0[j])
      Fx-mu.x0*integrate(func,0,z,mu.t0=mu.t0,mu.x0=mu.x0)$value
    }
    F.inv=inverse(F.Z3,lower=0,upper=100)
    ucl.gn3[i,j]=F.inv(1-0.027)
  }
}

#gamma-weibull
ucl.gw3=matrix(0, nrow=3, ncol=3, dimnames = list(c("1","2","5"),c("1","2","5")))

for (i in 1:3){
  for (j in 1:3) {  
    F.Z3=function(z){
      func=function(x,mu.t0,mu.x0){
        pgamma((1/(z-x))*mu.t0,shape=(mu.t0/s0[i])^2,scale=(s0[i]^2/mu.t0))*dweibull(x*mu.x0,shape=a0[j],scale=mu.x0/gamma((1/a0[j])+1))
      }
      Fx=pweibull(z*mu.x0,shape=a0[j],scale=mu.x0/gamma((1/a0[j])+1))
      Fx-mu.x0*integrate(func,0,z,mu.t0=mu.t0,mu.x0=mu.x0)$value
    }
    F.inv=inverse(F.Z3,lower=0,upper=100)
    ucl.gw3[i,j]=F.inv(1-0.027)
  }
}

#lognormal-gamma
ucl.lg3=matrix(0, nrow=3, ncol=3, dimnames = list(c("1","2","5"),c("1","2","5")))

for (i in 1:3){
  for (j in 1:3) {  
    F.Z3=function(z){
      func=function(x,mu.t0,mu.x0){
        plnorm((1/(z-x))*mu.t0,meanlog=log(mu.t0^2/sqrt(s0[i]^2+mu.t0^2)),sdlog=sqrt(log(1+(s0[i]^2/mu.t0^2))))*dgamma(x*mu.x0,shape=(mu.x0/s0[j])^2,scale=(s0[j]^2/mu.x0))
      }
      Fx=pgamma(z*mu.x0,shape=(mu.x0/s0[j])^2,scale=(s0[j]^2/mu.x0))
      Fx-mu.x0*integrate(func,0,z,mu.t0=mu.t0,mu.x0=mu.x0)$value
    }
    F.inv=inverse(F.Z3,lower=0,upper=100)
    ucl.lg3[i,j]=F.inv(1-0.027)
  }
}

#lognormal-lognormal
ucl.ll3=matrix(0, nrow=3, ncol=3, dimnames = list(c("1","2","5"),c("1","2","5")))

for (i in 1:3){
  for (j in 1:3) {  
    F.Z3=function(z){
      func=function(x,mu.t0,mu.x0){
        plnorm((1/(z-x))*mu.t0,meanlog=log(mu.t0^2/sqrt(s0[i]^2+mu.t0^2)),sdlog=sqrt(log(1+(s0[i]^2/mu.t0^2))))*dlnorm(x*mu.x0,meanlog=log(mu.x0^2/sqrt(s0[j]^2+mu.x0^2)),sdlog=sqrt(log(1+(s0[j]^2/mu.x0^2))))
      }
      Fx=plnorm(z*mu.x0,meanlog=log(mu.x0^2/sqrt(s0[j]^2 + mu.x0^2)),sdlog=sqrt(log(1+(s0[j]^2/mu.x0^2))))
      Fx-mu.x0*integrate(func,0,z,mu.t0=mu.t0,mu.x0=mu.x0)$value
    }
    F.inv=inverse(F.Z3,lower=0,upper=100)
    ucl.ll3[i,j]=F.inv(1-0.027)
  }
}

#lognormal-normal
ucl.ln3=matrix(0, nrow=3, ncol=2, dimnames = list(c("1","2","5"),c("1","2")))

for (i in 1:3){
  for (j in 1:2) {  
    F.Z3=function(z){
      func=function(x,mu.t0,mu.x0){
        plnorm((1/(z-x))*mu.t0,meanlog=log(mu.t0^2/sqrt(s0[i]^2+mu.t0^2)),sdlog=sqrt(log(1+(s0[i]^2/mu.t0^2))))*dnorm(x*mu.x0,mean=mu.x0,sd=s0[j])
      }
      Fx=pnorm(z*mu.x0,mean=mu.x0,sd=s0[j])
      Fx-mu.x0*integrate(func,0,z,mu.t0=mu.t0,mu.x0=mu.x0)$value
    }
    F.inv=inverse(F.Z3,lower=0,upper=100)
    ucl.ln3[i,j]=F.inv(1-0.027)
  }
}

#lognormal-weibull
ucl.lw3=matrix(0, nrow=3, ncol=3, dimnames = list(c("1","2","5"),c("1","2","5")))

for (i in 1:3){
  for (j in 1:3) {  
    F.Z3=function(z){
      func=function(x,mu.t0,mu.x0){
        plnorm((1/(z-x))*mu.t0,meanlog=log(mu.t0^2/sqrt(s0[i]^2+mu.t0^2)),sdlog=sqrt(log(1+(s0[i]^2/mu.t0^2))))*dweibull(x*mu.x0,shape=a0[j],scale=mu.x0/gamma((1/a0[j])+1))
      }
      Fx=pweibull(z*mu.x0,shape=a0[j],scale=mu.x0/gamma((1/a0[j])+1))
      Fx-mu.x0*integrate(func,0,z,mu.t0=mu.t0,mu.x0=mu.x0)$value
    }
    F.inv=inverse(F.Z3,lower=0,upper=100)
    ucl.lw3[i,j]=F.inv(1-0.027)
  }
}

#weibull-gamma
ucl.wg3=matrix(0, nrow=3, ncol=3, dimnames = list(c("1","2","5"),c("1","2","5")))

for (i in 1:3){
  for (j in 1:3) {  
    F.Z3=function(z){
      func=function(x,mu.t0,mu.x0){
        pweibull((1/(z-x))*mu.t0,shape=a0[i],scale=mu.t0/gamma((1/a0[i])+1))*dgamma(x*mu.x0,shape=(mu.x0/s0[j])^2,scale=(s0[j]^2/mu.x0))
      }
      Fx=pgamma(z*mu.x0,shape=(mu.x0/s0[j])^2,scale=(s0[j]^2/mu.x0))
      Fx-mu.x0*integrate(func,0,z,mu.t0=mu.t0,mu.x0=mu.x0)$value
    }
    F.inv=inverse(F.Z3,lower=0,upper=100)
    ucl.wg3[i,j]=F.inv(1-0.027)
  }
}

#weibull-lognormal
ucl.wl3=matrix(0, nrow=3, ncol=3, dimnames = list(c("1","2","5"),c("1","2","5")))

for (i in 1:3){
  for (j in 1:3) {  
    F.Z3=function(z){
      func=function(x,mu.t0,mu.x0){
        pweibull((1/(z-x))*mu.t0,shape=a0[i],scale=mu.t0/gamma((1/a0[i])+1))*dlnorm(x*mu.x0,meanlog=log(mu.x0^2/sqrt(s0[j]^2+mu.x0^2)),sdlog=sqrt(log(1+(s0[j]^2/mu.x0^2))))
      }
      Fx=plnorm(z*mu.x0,meanlog=log(mu.x0^2/sqrt(s0[j]^2 + mu.x0^2)),sdlog=sqrt(log(1+(s0[j]^2/mu.x0^2))))
      Fx-mu.x0*integrate(func,0,z,mu.t0=mu.t0,mu.x0=mu.x0)$value
    }
    F.inv=inverse(F.Z3,lower=0,upper=100)
    ucl.wl3[i,j]=F.inv(1-0.027)
  }
}

#weibull-normal
ucl.wn3=matrix(0, nrow=3, ncol=2, dimnames = list(c("1","2","5"),c("1","2")))

for (i in 1:3){
  for (j in 1:2) {  
    F.Z3=function(z){
      func=function(x,mu.t0,mu.x0){
        pweibull((1/(z-x))*mu.t0,shape=a0[i],scale=mu.t0/gamma((1/a0[i])+1))*dnorm(x*mu.x0,mean=mu.x0,sd=s0[j])
      }
      Fx=pnorm(z*mu.x0,mean=mu.x0,sd=s0[j])
      Fx-mu.x0*integrate(func,0,z,mu.t0=mu.t0,mu.x0=mu.x0)$value
    }
    F.inv=inverse(F.Z3,lower=0,upper=100)
    ucl.wn3[i,j]=F.inv(1-0.027)
  }
}

#weibull-weibull
ucl.ww3=matrix(0, nrow=3, ncol=3, dimnames = list(c("1","2","5"),c("1","2","5")))

for (i in 1:3){
  for (j in 1:3) {  
    F.Z3=function(z){
      func=function(x,mu.t0,mu.x0){
        pweibull((1/(z-x))*mu.t0,shape=a0[i],scale=mu.t0/gamma((1/a0[i])+1))*dweibull(x*mu.x0,shape=a0[j],scale=mu.x0/gamma((1/a0[j])+1))
      }
      Fx=pweibull(z*mu.x0,shape=a0[j],scale=mu.x0/gamma((1/a0[j])+1))
      Fx-mu.x0*integrate(func,0,z,mu.t0=mu.t0,mu.x0=mu.x0)$value
    }
    F.inv=inverse(F.Z3,lower=0,upper=100)
    ucl.ww3[i,j]=F.inv(1-0.027)
  }
}    
#*****EATSt values for Statictic Z3
#gamma-gamma
EATSt.gg3=matrix(0, nrow=3, ncol=3, dimnames = list(c("1","2","5"),c("1","2","5")))

for (i in 1:3){
  for (j in 1:3) {  
    for (k in 1:10){
      F.Z3=function(z){
        func=function(x,mu.t0,mu.x0){
          pgamma((1/(z-x))*mu.t0,shape=(deltat[k]*mu.t0/s0[i])^2,scale=(s0[i]^2/(deltat[k]*mu.t0)))*dgamma(x*mu.x0,shape=(mu.x0/s0[j])^2,scale=(s0[j]^2/mu.x0))
        }
        Fx=pgamma(z*mu.x0,shape=(mu.x0/s0[j])^2,scale=(s0[j]^2/mu.x0))
        Fx-mu.x0*integrate(func,0,z,mu.t0=mu.t0,mu.x0=mu.x0)$value
      }
      beta=F.Z3(ucl.gg3[i,j])
      ATSt[k]=(deltat[k]*mu.t0)/(1-beta)
    }
    EATSt.gg3[i,j]= mean(ATSt)
  }
}

#gamma-lognormal
EATSt.gl3=matrix(0, nrow=3, ncol=3, dimnames = list(c("1","2","5"),c("1","2","5")))

for (i in 1:3){
  for (j in 1:3) {  
    for (k in 1:10){
      F.Z3=function(z){
        func=function(x,mu.t0,mu.x0){
          pgamma((1/(z-x))*mu.t0,shape=(deltat[k]*mu.t0/s0[i])^2,scale=(s0[i]^2/(deltat[k]*mu.t0)))*dlnorm(x*mu.x0,meanlog=log(mu.x0^2/sqrt(s0[j]^2+mu.x0^2)),sdlog=sqrt(log(1+(s0[j]^2/mu.x0^2))))
        }
        Fx=plnorm(z*mu.x0,meanlog=log(mu.x0^2/sqrt(s0[j]^2 + mu.x0^2)),sdlog=sqrt(log(1+(s0[j]^2/mu.x0^2))))
        Fx-mu.x0*integrate(func,0,z,mu.t0=mu.t0,mu.x0=mu.x0)$value
      }
      beta=F.Z3(ucl.gl3[i,j])
      ATSt[k]=(deltat[k]*mu.t0)/(1-beta)
    }
    EATSt.gl3[i,j]= mean(ATSt)
  }
}

#gamma-normal
EATSt.gn3=matrix(0, nrow=3, ncol=2, dimnames = list(c("1","2","5"),c("1","2")))

for (i in 1:3){
  for (j in 1:2) {  
    for (k in 1:10){
      F.Z3=function(z){
        func=function(x,mu.t0,mu.x0){
          pgamma((1/(z-x))*mu.t0,shape=(deltat[k]*mu.t0/s0[i])^2,scale=(s0[i]^2/(deltat[k]*mu.t0)))*dnorm(x*mu.x0,mean=mu.x0,sd=s0[j])
        }
        Fx=pnorm(z*mu.x0,mean=mu.x0,sd=s0[j])
        Fx-mu.x0*integrate(func,0,z,mu.t0=mu.t0,mu.x0=mu.x0)$value
      }
      beta=F.Z3(ucl.gn3[i,j])
      ATSt[k]=(deltat[k]*mu.t0)/(1-beta)
    }
    EATSt.gn3[i,j]= mean(ATSt)
  }
}

#gamma-weibull
EATSt.gw3=matrix(0, nrow=3, ncol=3, dimnames = list(c("1","2","5"),c("1","2","5")))

for (i in 1:3){
  for (j in 1:3) {  
    for (k in 1:10){
      F.Z3=function(z){
        func=function(x,mu.t0,mu.x0){
          pgamma((1/(z-x))*mu.t0,shape=(deltat[k]*mu.t0/s0[i])^2,scale=(s0[i]^2/(deltat[k]*mu.t0)))*dweibull(x*mu.x0,shape=a0[j],scale=mu.x0/gamma((1/a0[j])+1))
        }
        Fx=pweibull(z*mu.x0,shape=a0[j],scale=mu.x0/gamma((1/a0[j])+1))
        Fx-mu.x0*integrate(func,0,z,mu.t0=mu.t0,mu.x0=mu.x0)$value
      }
      beta=F.Z3(ucl.gw3[i,j])
      ATSt[k]=(deltat[k]*mu.t0)/(1-beta)
    }
    EATSt.gw3[i,j]= mean(ATSt)
  }
}

#lognormal-gamma
EATSt.lg3=matrix(0, nrow=3, ncol=3, dimnames = list(c("1","2","5"),c("1","2","5")))

for (i in 1:3){
  for (j in 1:3) {  
    for (k in 1:10){
      F.Z3=function(z){
        func=function(x,mu.t0,mu.x0){
          plnorm((1/(z-x))*mu.t0,meanlog=log((deltat[k]*mu.t0)^2/sqrt(s0[i]^2+(deltat[k]*mu.t0)^2)),sdlog=sqrt(log(1+(s0[i]^2/(deltat[k]*mu.t0)^2))))*dgamma(x*mu.x0,shape=(mu.x0/s0[j])^2,scale=(s0[j]^2/mu.x0))
        }
        Fx=pgamma(z*mu.x0,shape=(mu.x0/s0[j])^2,scale=(s0[j]^2/mu.x0))
        Fx-mu.x0*integrate(func,0,z,mu.t0=mu.t0,mu.x0=mu.x0)$value
      }
      beta=F.Z3(ucl.lg3[i,j])
      ATSt[k]=(deltat[k]*mu.t0)/(1-beta)
    }
    EATSt.lg3[i,j]= mean(ATSt)
  }
}

#lognormal-lognormal
EATSt.ll3=matrix(0, nrow=3, ncol=3, dimnames = list(c("1","2","5"),c("1","2","5")))

for (i in 1:3){
  for (j in 1:3) {  
    for (k in 1:10){
      F.Z3=function(z){
        func=function(x,mu.t0,mu.x0){
          plnorm((1/(z-x))*mu.t0,meanlog=log((deltat[k]*mu.t0)^2/sqrt(s0[i]^2+(deltat[k]*mu.t0)^2)),sdlog=sqrt(log(1+(s0[i]^2/(deltat[k]*mu.t0)^2))))*dlnorm(x*mu.x0,meanlog=log(mu.x0^2/sqrt(s0[j]^2+mu.x0^2)),sdlog=sqrt(log(1+(s0[j]^2/mu.x0^2))))
        }
        Fx=plnorm(z*mu.x0,meanlog=log(mu.x0^2/sqrt(s0[j]^2 + mu.x0^2)),sdlog=sqrt(log(1+(s0[j]^2/mu.x0^2))))
        Fx-mu.x0*integrate(func,0,z,mu.t0=mu.t0,mu.x0=mu.x0)$value
      }
      beta=F.Z3(ucl.ll3[i,j])
      ATSt[k]=(deltat[k]*mu.t0)/(1-beta)
    }
    EATSt.ll3[i,j]= mean(ATSt)
  }
}

#lognormal-normal
EATSt.ln3=matrix(0, nrow=3, ncol=2, dimnames = list(c("1","2","5"),c("1","2")))

for (i in 1:3){
  for (j in 1:2) {  
    for (k in 1:10){
      F.Z3=function(z){
        func=function(x,mu.t0,mu.x0){
          plnorm((1/(z-x))*mu.t0,meanlog=log((deltat[k]*mu.t0)^2/sqrt(s0[i]^2+(deltat[k]*mu.t0)^2)),sdlog=sqrt(log(1+(s0[i]^2/(deltat[k]*mu.t0)^2))))*dnorm(x*mu.x0,mean=mu.x0,sd=s0[j])
        }
        Fx=pnorm(z*mu.x0,mean=mu.x0,sd=s0[j])
        Fx-mu.x0*integrate(func,0,z,mu.t0=mu.t0,mu.x0=mu.x0)$value
      }
      beta=F.Z3(ucl.ln3[i,j])
      ATSt[k]=(deltat[k]*mu.t0)/(1-beta)
    }
    EATSt.ln3[i,j]= mean(ATSt)
  }
}

#lognormal-weibull
EATSt.lw3=matrix(0, nrow=3, ncol=3, dimnames = list(c("1","2","5"),c("1","2","5")))

for (i in 1:3){
  for (j in 1:3) {  
    for (k in 1:10){
      F.Z3=function(z){
        func=function(x,mu.t0,mu.x0){
          plnorm((1/(z-x))*mu.t0,meanlog=log((deltat[k]*mu.t0)^2/sqrt(s0[i]^2+(deltat[k]*mu.t0)^2)),sdlog=sqrt(log(1+(s0[i]^2/(deltat[k]*mu.t0)^2))))*dweibull(x*mu.x0,shape=a0[j],scale=mu.x0/gamma((1/a0[j])+1))
        }
        Fx=pweibull(z*mu.x0,shape=a0[j],scale=mu.x0/gamma((1/a0[j])+1))
        Fx-mu.x0*integrate(func,0,z,mu.t0=mu.t0,mu.x0=mu.x0)$value
      }
      beta=F.Z3(ucl.lw3[i,j])
      ATSt[k]=(deltat[k]*mu.t0)/(1-beta)
    }
    EATSt.lw3[i,j]= mean(ATSt)
  }
}

#weibull-gamma
EATSt.wg3=matrix(0, nrow=3, ncol=3, dimnames = list(c("1","2","5"),c("1","2","5")))

for (i in 1:3){
  for (j in 1:3) {  
    for (k in 1:10){
      F.Z3=function(z){
        func=function(x,mu.t0,mu.x0){
          pweibull((1/(z-x))*mu.t0,shape=a1[i,k],scale=deltat[k]*mu.t0/gamma((1/a1[i,k])+1))*dgamma(x*mu.x0,shape=(mu.x0/s0[j])^2,scale=(s0[j]^2/mu.x0))
        }
        Fx=pgamma(z*mu.x0,shape=(mu.x0/s0[j])^2,scale=(s0[j]^2/mu.x0))
        Fx-mu.x0*integrate(func,0,z,mu.t0=mu.t0,mu.x0=mu.x0)$value
      }
      beta=F.Z3(ucl.wg3[i,j])
      ATSt[k]=(deltat[k]*mu.t0)/(1-beta)
    }
    EATSt.wg3[i,j]= mean(ATSt)
  }
}

#weibull-lognormal
EATSt.wl3=matrix(0, nrow=3, ncol=3, dimnames = list(c("1","2","5"),c("1","2","5")))

for (i in 1:3){
  for (j in 1:3) {  
    for (k in 1:10){
      F.Z3=function(z){
        func=function(x,mu.t0,mu.x0){
          pweibull((1/(z-x))*mu.t0,shape=a1[i,k],scale=deltat[k]*mu.t0/gamma((1/a1[i,k])+1))*dlnorm(x*mu.x0,meanlog=log(mu.x0^2/sqrt(s0[j]^2+mu.x0^2)),sdlog=sqrt(log(1+(s0[j]^2/mu.x0^2))))
        }
        Fx=plnorm(z*mu.x0,meanlog=log(mu.x0^2/sqrt(s0[j]^2 + mu.x0^2)),sdlog=sqrt(log(1+(s0[j]^2/mu.x0^2))))
        Fx-mu.x0*integrate(func,0,z,mu.t0=mu.t0,mu.x0=mu.x0)$value
      }
      beta=F.Z3(ucl.wl3[i,j])
      ATSt[k]=(deltat[k]*mu.t0)/(1-beta)
    }
    EATSt.wl3[i,j]= mean(ATSt)
  }
}

#weibull-normal
EATSt.wn3=matrix(0, nrow=3, ncol=2, dimnames = list(c("1","2","5"),c("1","2")))

for (i in 1:3){
  for (j in 1:2) {  
    for (k in 1:10){
      F.Z3=function(z){
        func=function(x,mu.t0,mu.x0){
          pweibull((1/(z-x))*mu.t0,shape=a1[i,k],scale=deltat[k]*mu.t0/gamma((1/a1[i,k])+1))*dnorm(x*mu.x0,mean=mu.x0,sd=s0[j])
        }
        Fx=pnorm(z*mu.x0,mean=mu.x0,sd=s0[j])
        Fx-mu.x0*integrate(func,0,z,mu.t0=mu.t0,mu.x0=mu.x0)$value
      }
      beta=F.Z3(ucl.wn3[i,j])
      ATSt[k]=(deltat[k]*mu.t0)/(1-beta)
    }
    EATSt.wn3[i,j]= mean(ATSt)
  }
}

#weibull-weibull
EATSt.ww3=matrix(0, nrow=3, ncol=3, dimnames = list(c("1","2","5"),c("1","2","5")))

for (i in 1:3){
  for (j in 1:3) {  
    for (k in 1:10){
      F.Z3=function(z){
        func=function(x,mu.t0,mu.x0){
          pweibull((1/(z-x))*mu.t0,shape=a1[i,k],scale=deltat[k]*mu.t0/gamma((1/a1[i,k])+1))*dweibull(x*mu.x0,shape=a0[j],scale=mu.x0/gamma((1/a0[j])+1))
        }
        Fx=pweibull(z*mu.x0,shape=a0[j],scale=mu.x0/gamma((1/a0[j])+1))
        Fx-mu.x0*integrate(func,0,z,mu.t0=mu.t0,mu.x0=mu.x0)$value
      }
      beta=F.Z3(ucl.ww3[i,j])
      ATSt[k]=(deltat[k]*mu.t0)/(1-beta)
    }
    EATSt.ww3[i,j]= mean(ATSt)
  }
}
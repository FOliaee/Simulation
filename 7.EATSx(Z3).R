#*****EATSx values for Statictic Z3
#gamma-gamma
EATSx.gg3=matrix(0, nrow=3, ncol=3, dimnames = list(c("1","2","5"),c("1","2","5")))

for (i in 1:3){
  for (j in 1:3) {  
    for (k in 1:10){
      F.Z1=function(z){
        func=function(x,mu.t0,mu.x0){
          pgamma((1/(z-x))*mu.t0,shape=(mu.t0/s0[i])^2,scale=(s0[i]^2/mu.t0))*dgamma(x*mu.x0,shape=(deltax[k]*mu.x0/s0[j])^2,scale=(s0[j]^2/(deltax[k]*mu.x0)))
        }
        Fx=pgamma(z*mu.x0,shape=(deltax[k]*mu.x0/s0[j])^2,scale=(s0[j]^2/(deltax[k]*mu.x0)))
        Fx-mu.x0*integrate(func,0,z,mu.t0=mu.t0,mu.x0=mu.x0)$value
      }
      beta=F.Z1(ucl.gg3[i,j])
      ATSx[k]=mu.t0/(1-beta)
    }
    EATSx.gg3[i,j]= mean(ATSx)
  }
}

#gamma-lognormal
EATSx.gl3=matrix(0, nrow=3, ncol=3, dimnames = list(c("1","2","5"),c("1","2","5")))

for (i in 1:3){
  for (j in 1:3) {  
    for (k in 1:10){
      F.Z3=function(z){
        func=function(x,mu.t0,mu.x0){
          pgamma((1/(z-x))*mu.t0,shape=(mu.t0/s0[i])^2,scale=(s0[i]^2/mu.t0))*dlnorm(x*mu.x0,meanlog=log((deltax[k]*mu.x0)^2/sqrt(s0[j]^2+(deltax[k]*mu.x0)^2)),sdlog=sqrt(log(1+(s0[j]^2/(deltax[k]*mu.x0)^2))))
        }
        Fx=plnorm(z*mu.x0,meanlog=log((deltax[k]*mu.x0)^2/sqrt(s0[j]^2+(deltax[k]*mu.x0)^2)),sdlog=sqrt(log(1+(s0[j]^2/(deltax[k]*mu.x0)^2))))
        Fx-mu.x0*integrate(func,0,z,mu.t0=mu.t0,mu.x0=mu.x0)$value
      }
      beta=F.Z3(ucl.gl3[i,j])
      ATSx[k]=mu.t0/(1-beta)
    }
    EATSx.gl3[i,j]= mean(ATSx)
  }
}

#gamma-normal
EATSx.gn3=matrix(0, nrow=3, ncol=2, dimnames = list(c("1","2","5"),c("1","2")))

for (i in 1:3){
  for (j in 1:2) {  
    for (k in 1:10){
      F.Z3=function(z){
        func=function(x,mu.t0,mu.x0){
          pgamma((1/(z-x))*mu.t0,shape=(mu.t0/s0[i])^2,scale=(s0[i]^2/mu.t0))*dnorm(x*mu.x0,mean=deltax[k]*mu.x0,sd=s0[j])
        }
        Fx=pnorm(z*mu.x0,mean=deltax[k]*mu.x0,sd=s0[j])
        Fx-mu.x0*integrate(func,0,z,mu.t0=mu.t0,mu.x0=mu.x0)$value
      }
      beta=F.Z3(ucl.gn3[i,j])
      ATSx[k]=mu.t0/(1-beta)
    }
    EATSx.gn3[i,j]= mean(ATSx)
  }
}

#gamma-weibull
EATSx.gw3=matrix(0, nrow=3, ncol=3, dimnames = list(c("1","2","5"),c("1","2","5")))

for (i in 1:3){
  for (j in 1:3) {  
    for (k in 1:10){
      F.Z3=function(z){
        func=function(x,mu.t0,mu.x0){
          pgamma((1/(z-x))*mu.t0,shape=(mu.t0/s0[i])^2,scale=(s0[i]^2/mu.t0))*dweibull(x*mu.x0,shape=a1[j,k],scale=(deltax[k]*mu.x0)/gamma((1/a1[j,k])+1))
        }
        Fx=pweibull(z*mu.x0,shape=a1[j,k],scale=deltax[k]*mu.x0/gamma((1/a1[j,k])+1))
        Fx-mu.x0*integrate(func,0,z,mu.t0=mu.t0,mu.x0=mu.x0)$value
      }      
      beta=F.Z3(ucl.gw3[i,j])
      ATSx[k]=mu.t0/(1-beta)
    }
    EATSx.gw3[i,j]= mean(ATSx)
  }
}

#lognormal-gamma
EATSx.lg3=matrix(0, nrow=3, ncol=3, dimnames = list(c("1","2","5"),c("1","2","5")))

for (i in 1:3){
  for (j in 1:3) {  
    for (k in 1:10){
      F.Z3=function(z){
        func=function(x,mu.t0,mu.x0){
          plnorm((1/(z-x))*mu.t0,meanlog=log(mu.t0^2/sqrt(s0[i]^2+mu.t0^2)),sdlog=sqrt(log(1+(s0[i]^2/mu.t0^2))))*dgamma(x*mu.x0,shape=(deltax[k]*mu.x0/s0[j])^2,scale=(s0[j]^2/(deltax[k]*mu.x0)))
        }
        Fx=pgamma(z*mu.x0,shape=(deltax[k]*mu.x0/s0[j])^2,scale=(s0[j]^2/(deltax[k]*mu.x0)))
        Fx-mu.x0*integrate(func,0,z,mu.t0=mu.t0,mu.x0=mu.x0)$value
      }
      beta=F.Z3(ucl.lg3[i,j])
      ATSx[k]=mu.t0/(1-beta)
    }
    EATSx.lg3[i,j]= mean(ATSx)
  }
}

#lognormal-lognormal
EATSx.ll3=matrix(0, nrow=3, ncol=3, dimnames = list(c("1","2","5"),c("1","2","5")))

for (i in 1:3){
  for (j in 1:3) {  
    for (k in 1:10){
      F.Z3=function(z){
        func=function(x,mu.t0,mu.x0){
          plnorm((1/(z-x))*mu.t0,meanlog=log(mu.t0^2/sqrt(s0[i]^2+mu.t0^2)),sdlog=sqrt(log(1+(s0[i]^2/mu.t0^2))))*dlnorm(x*mu.x0,meanlog=log((deltax[k]*mu.x0)^2/sqrt(s0[j]^2+(deltax[k]*mu.x0)^2)),sdlog=sqrt(log(1+(s0[j]^2/(deltax[k]*mu.x0)^2))))
        }
        Fx=plnorm(z*mu.x0,meanlog=log((deltax[k]*mu.x0)^2/sqrt(s0[j]^2+(deltax[k]*mu.x0)^2)),sdlog=sqrt(log(1+(s0[j]^2/(deltax[k]*mu.x0)^2))))
        Fx-mu.x0*integrate(func,0,z,mu.t0=mu.t0,mu.x0=mu.x0)$value
      }
      beta=F.Z3(ucl.ll3[i,j])
      ATSx[k]=mu.t0/(1-beta)
    }
    EATSx.ll3[i,j]= mean(ATSx)
  }
}

#lognormal-normal
EATSx.ln3=matrix(0, nrow=3, ncol=2, dimnames = list(c("1","2","5"),c("1","2")))

for (i in 1:3){
  for (j in 1:2) {  
    for (k in 1:10){
      F.Z3=function(z){
        func=function(x,mu.t0,mu.x0){
          plnorm((1/(z-x))*mu.t0,meanlog=log(mu.t0^2/sqrt(s0[i]^2+mu.t0^2)),sdlog=sqrt(log(1+(s0[i]^2/mu.t0^2))))*dnorm(x*mu.x0,mean=deltax[k]*mu.x0,sd=s0[j])
        }
        Fx=pnorm(z*mu.x0,mean=deltax[k]*mu.x0,sd=s0[j])
        Fx-mu.x0*integrate(func,0,z,mu.t0=mu.t0,mu.x0=mu.x0)$value
      }
      beta=F.Z3(ucl.ln3[i,j])
      ATSx[k]=mu.t0/(1-beta)
    }
    EATSx.ln3[i,j]= mean(ATSx)
  }
}

#lognormal-weibull
EATSx.lw3=matrix(0, nrow=3, ncol=3, dimnames = list(c("1","2","5"),c("1","2","5")))

for (i in 1:3){
  for (j in 1:3) {  
    for (k in 1:10){
      F.Z3=function(z){
        func=function(x,mu.t0,mu.x0){
          plnorm((1/(z-x))*mu.t0,meanlog=log(mu.t0^2/sqrt(s0[i]^2+mu.t0^2)),sdlog=sqrt(log(1+(s0[i]^2/mu.t0^2))))*dweibull(x*mu.x0,shape=a1[j,k],scale=(deltax[k]*mu.x0)/gamma((1/a1[j,k])+1))
        }
        Fx=pweibull(z*mu.x0,shape=a1[j,k],scale=deltax[k]*mu.x0/gamma((1/a1[j,k])+1))
        Fx-mu.x0*integrate(func,0,z,mu.t0=mu.t0,mu.x0=mu.x0)$value
      } 
      beta=F.Z3(ucl.lw3[i,j])
      ATSx[k]=mu.t0/(1-beta)
    }
    EATSx.lw3[i,j]= mean(ATSx)
  }
}

#weibull-gamma
EATSx.wg3=matrix(0, nrow=3, ncol=3, dimnames = list(c("1","2","5"),c("1","2","5")))

for (i in 1:3){
  for (j in 1:3) {  
    for (k in 1:10){
      F.Z3=function(z){
        func=function(x,mu.t0,mu.x0){
          pweibull((1/(z-x))*mu.t0,shape=a0[i],scale=mu.t0/gamma((1/a0[i])+1))*dgamma(x*mu.x0,shape=(deltax[k]*mu.x0/s0[j])^2,scale=(s0[j]^2/(deltax[k]*mu.x0)))
        }
        Fx=pgamma(z*mu.x0,shape=(deltax[k]*mu.x0/s0[j])^2,scale=(s0[j]^2/(deltax[k]*mu.x0)))
        Fx-mu.x0*integrate(func,0,z,mu.t0=mu.t0,mu.x0=mu.x0)$value
      }
      beta=F.Z3(ucl.wg3[i,j])
      ATSx[k]=mu.t0/(1-beta)
    }
    EATSx.wg3[i,j]= mean(ATSx)
  }
}

#weibull-lognormal
EATSx.wl3=matrix(0, nrow=3, ncol=3, dimnames = list(c("1","2","5"),c("1","2","5")))
ucl.wl3=matrix(0, nrow=3, ncol=3, dimnames = list(c("1","2","5"),c("1","2","5")))

for (i in 1:3){
  for (j in 1:3) {  
    for (k in 1:10){
      F.Z3=function(z){
        func=function(x,mu.t0,mu.x0){
          pweibull((1/(z-x))*mu.t0,shape=a0[i],scale=mu.t0/gamma((1/a0[i])+1))*dlnorm(x*mu.x0,meanlog=log((deltax[k]*mu.x0)^2/sqrt(s0[j]^2+(deltax[k]*mu.x0)^2)),sdlog=sqrt(log(1+(s0[j]^2/(deltax[k]*mu.x0)^2))))
        }
        Fx=plnorm(z*mu.x0,meanlog=log((deltax[k]*mu.x0)^2/sqrt(s0[j]^2+(deltax[k]*mu.x0)^2)),sdlog=sqrt(log(1+(s0[j]^2/(deltax[k]*mu.x0)^2))))
        Fx-mu.x0*integrate(func,0,z,mu.t0=mu.t0,mu.x0=mu.x0)$value
      }
      beta=F.Z3(ucl.wl3[i,j])
      ATSx[k]=mu.t0/(1-beta)
    }
    EATSx.wl3[i,j]= mean(ATSx)
  }
}

#weibull-normal
EATSx.wn3=matrix(0, nrow=3, ncol=2, dimnames = list(c("1","2","5"),c("1","2")))

for (i in 1:3){
  for (j in 1:2) {  
    for (k in 1:10){
      F.Z3=function(z){
        func=function(x,mu.t0,mu.x0){
          pweibull((1/(z-x))*mu.t0,shape=a0[i],scale=mu.t0/gamma((1/a0[i])+1))*dnorm(x*mu.x0,mean=deltax[k]*mu.x0,sd=s0[j])
        }
        Fx=pnorm(z*mu.x0,mean=deltax[k]*mu.x0,sd=s0[j])
        Fx-mu.x0*integrate(func,0,z,mu.t0=mu.t0,mu.x0=mu.x0)$value
      }
      beta=F.Z3(ucl.wn3[i,j])
      ATSx[k]=mu.t0/(1-beta)
    }
    EATSx.wn3[i,j]= mean(ATSx)
  }
}

#weibull-weibull
EATSx.ww3=matrix(0, nrow=3, ncol=3, dimnames = list(c("1","2","5"),c("1","2","5")))

for (i in 1:3){
  for (j in 1:3) {  
    for (k in 1:10){
      F.Z3=function(z){
        func=function(x,mu.t0,mu.x0){
          pweibull((1/(z-x))*mu.t0,shape=a0[i],scale=mu.t0/gamma((1/a0[i])+1))*dweibull(x*mu.x0,shape=a1[j,k],scale=(deltax[k]*mu.x0)/gamma((1/a1[j,k])+1))
        }
        Fx=pweibull(z*mu.x0,shape=a1[j,k],scale=deltax[k]*mu.x0/gamma((1/a1[j,k])+1))
        Fx-mu.x0*integrate(func,0,z,mu.t0=mu.t0,mu.x0=mu.x0)$value
      } 
      beta=F.Z3(ucl.ww3[i,j])
      ATSx[k]=mu.t0/(1-beta)
    }
    EATSx.ww3[i,j]= mean(ATSx)
  }
}
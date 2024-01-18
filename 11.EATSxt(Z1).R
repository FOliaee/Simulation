#*****EATSxt values for Statictic Z1
#gamma-gamma
EATSxt.gg1=matrix(0, nrow=3, ncol=3, dimnames = list(c("1","2","5"),c("1","2","5")))

for (i in 1:3){
  for (j in 1:3) {  
    for (k in 1:10){
      for (l in 1:10){
        F.Z1=function(z){
          func=function(x,mu.t0,mu.x0){
            pgamma((x-z)*mu.t0,shape=(deltat[k]*mu.t0/s0[i])^2,scale=(s0[i]^2/(deltat[k]*mu.t0)))*dgamma(x*mu.x0,shape=(deltax[l]*mu.x0/s0[j])^2,scale=(s0[j]^2/(deltax[l]*mu.x0)))
          }
          1-mu.x0*integrate(func,0,Inf,mu.t0=mu.t0,mu.x0=mu.x0)$value
        }
        beta=F.Z1(ucl.gg1[i,j])
        ATSxt[k,l]=(deltat[k]*mu.t0)/(1-beta)
      }
    }
    EATSxt.gg1[i,j]=mean(mean(ATSxt)) 
  }
}

#gamma-lognormal
EATSxt.gl1=matrix(0, nrow=3, ncol=3, dimnames = list(c("1","2","5"),c("1","2","5")))

for (i in 1:3){
  for (j in 1:3) {  
    for (k in 1:10){
      for (l in 1:10){
        F.Z1=function(z){
          func=function(x,mu.t0,mu.x0){
            pgamma((x-z)*mu.t0,shape=(deltat[k]*mu.t0/s0[i])^2,scale=(s0[i]^2/(deltat[k]*mu.t0)))*dlnorm(x*mu.x0,meanlog=log((deltax[l]*mu.x0)^2/sqrt(s0[j]^2+(deltax[l]*mu.x0)^2)),sdlog=sqrt(log(1+(s0[j]^2/(deltax[l]*mu.x0)^2))))
          }
          1-mu.x0*integrate(func,0,Inf,mu.t0=mu.t0,mu.x0=mu.x0)$value
        }
        beta=F.Z1(ucl.gl1[i,j])
        ATSxt[k,l]=(deltat[k]*mu.t0)/(1-beta)
      }
    }
    EATSxt.gl1[i,j]=mean(mean(ATSxt)) 
  }
}

#gamma-normal
EATSxt.gn1=matrix(0, nrow=3, ncol=2, dimnames = list(c("1","2","5"),c("1","2")))

for (i in 1:3){
  for (j in 1:2) {  
    for (k in 1:10){
      for (l in 1:10){
        F.Z1=function(z){
          func=function(x,mu.t0,mu.x0){
            pgamma((x-z)*mu.t0,shape=(deltat[k]*mu.t0/s0[i])^2,scale=(s0[i]^2/(deltat[k]*mu.t0)))*dnorm(x*mu.x0,mean=deltax[l]*mu.x0,sd=s0[j])
          }
          1-mu.x0*integrate(func,0,Inf,mu.t0=mu.t0,mu.x0=mu.x0)$value
        }
        beta=F.Z1(ucl.gn1[i,j])
        ATSxt[k,l]=(deltat[k]*mu.t0)/(1-beta)
      }
    }
    EATSxt.gn1[i,j]=mean(mean(ATSxt)) 
  }
}

#gamma-weibull
EATSxt.gw1=matrix(0, nrow=3, ncol=3, dimnames = list(c("1","2","5"),c("1","2","5")))

for (i in 1:3){
  for (j in 1:3) {  
    for (k in 1:10){
      for (l in 1:10){
        F.Z1=function(z){
          func=function(x,mu.t0,mu.x0){
            pgamma((x-z)*mu.t0,shape=(deltat[k]*mu.t0/s0[i])^2,scale=(s0[i]^2/(deltat[k]*mu.t0)))*dweibull(x*mu.x0,shape=a1[j,l],scale=(deltax[l]*mu.x0)/gamma((1/a1[j,l])+1))
          }
          1-mu.x0*integrate(func,0,Inf,mu.t0=mu.t0,mu.x0=mu.x0)$value
        }
        beta=F.Z1(ucl.gw1[i,j])
        ATSxt[k,l]=(deltat[k]*mu.t0)/(1-beta)
      }
    }
    EATSxt.gw1[i,j]=mean(mean(ATSxt)) 
  }
}

#lognormal-gamma
EATSxt.lg1=matrix(0, nrow=3, ncol=3, dimnames = list(c("1","2","5"),c("1","2","5")))

for (i in 1:3){
  for (j in 1:3) {  
    for (k in 1:10){
      for (l in 1:10){
        F.Z1=function(z){
          func=function(x,mu.t0,mu.x0){
            plnorm((x-z)*mu.t0,meanlog=log((deltat[k]*mu.t0)^2/sqrt(s0[i]^2+(deltat[k]*mu.t0)^2)),sdlog=sqrt(log(1+(s0[i]^2/(deltat[k]*mu.t0)^2))))*dgamma(x*mu.x0,shape=(deltax[l]*mu.x0/s0[j])^2,scale=(s0[j]^2/(deltax[l]*mu.x0)))
          }
          1-mu.x0*integrate(func,0,Inf,mu.t0=mu.t0,mu.x0=mu.x0)$value
        }
        beta=F.Z1(ucl.lg1[i,j])
        ATSxt[k,l]=(deltat[k]*mu.t0)/(1-beta)
      }
    }
    EATSxt.lg1[i,j]=mean(mean(ATSxt)) 
  }
}

#lognormal-lognormal
EATSxt.ll1=matrix(0, nrow=3, ncol=3, dimnames = list(c("1","2","5"),c("1","2","5")))

for (i in 1:3){
  for (j in 1:3) {  
    for (k in 1:10){
      for (l in 1:10){
        F.Z1=function(z){
          func=function(x,mu.t0,mu.x0){
            plnorm((x-z)*mu.t0,meanlog=log((deltat[k]*mu.t0)^2/sqrt(s0[i]^2+(deltat[k]*mu.t0)^2)),sdlog=sqrt(log(1+(s0[i]^2/(deltat[k]*mu.t0)^2))))*dlnorm(x*mu.x0,meanlog=log((deltax[l]*mu.x0)^2/sqrt(s0[j]^2+(deltax[l]*mu.x0)^2)),sdlog=sqrt(log(1+(s0[j]^2/(deltax[l]*mu.x0)^2))))
          }
          1-mu.x0*integrate(func,0,Inf,mu.t0=mu.t0,mu.x0=mu.x0)$value
        }
        beta=F.Z1(ucl.ll1[i,j])
        ATSxt[k,l]=(deltat[k]*mu.t0)/(1-beta)
      }
    }
    EATSxt.ll1[i,j]=mean(mean(ATSxt)) 
  }
}

#lognormal-normal
EATSxt.ln1=matrix(0, nrow=3, ncol=2, dimnames = list(c("1","2","5"),c("1","2")))

for (i in 1:3){
  for (j in 1:2) {  
    for (k in 1:10){
      for (l in 1:10){
        F.Z1=function(z){
          func=function(x,mu.t0,mu.x0){
            plnorm((x-z)*mu.t0,meanlog=log((deltat[k]*mu.t0)^2/sqrt(s0[i]^2+(deltat[k]*mu.t0)^2)),sdlog=sqrt(log(1+(s0[i]^2/(deltat[k]*mu.t0)^2))))*dnorm(x*mu.x0,mean=deltax[l]*mu.x0,sd=s0[j])
          }
          1-mu.x0*integrate(func,0,Inf,mu.t0=mu.t0,mu.x0=mu.x0)$value
        }
        beta=F.Z1(ucl.ln1[i,j])
        ATSxt[k,l]=(deltat[k]*mu.t0)/(1-beta)
      }
    }
    EATSxt.ln1[i,j]=mean(mean(ATSxt)) 
  }
}

#lognormal-weibull
EATSxt.lw1=matrix(0, nrow=3, ncol=3, dimnames = list(c("1","2","5"),c("1","2","5")))

for (i in 1:3){
  for (j in 1:3) {  
    for (k in 1:10){
      for (l in 1:10){
        F.Z1=function(z){
          func=function(x,mu.t0,mu.x0){
            plnorm((x-z)*mu.t0,meanlog=log((deltat[k]*mu.t0)^2/sqrt(s0[i]^2+(deltat[k]*mu.t0)^2)),sdlog=sqrt(log(1+(s0[i]^2/(deltat[k]*mu.t0)^2))))*dweibull(x*mu.x0,shape=a1[j,l],scale=(deltax[l]*mu.x0)/gamma((1/a1[j,l])+1))
          }
          1-mu.x0*integrate(func,0,Inf,mu.t0=mu.t0,mu.x0=mu.x0)$value
        }
        beta=F.Z1(ucl.lw1[i,j])
        ATSxt[k,l]=(deltat[k]*mu.t0)/(1-beta)
      }
    }
    EATSxt.lw1[i,j]=mean(mean(ATSxt)) 
  }
}

#weibull-gamma
EATSxt.wg1=matrix(0, nrow=3, ncol=3, dimnames = list(c("1","2","5"),c("1","2","5")))

for (i in 1:3){
  for (j in 1:3) {  
    for (k in 1:10){
      for (l in 1:10){
        F.Z1=function(z){
          func=function(x,mu.t0,mu.x0){
            pweibull((x-z)*mu.t0,shape=a1[i,k],scale=deltat[k]*mu.t0/gamma((1/a1[i,k])+1))*dgamma(x*mu.x0,shape=(deltax[l]*mu.x0/s0[j])^2,scale=(s0[j]^2/(deltax[l]*mu.x0)))
          }
          1-mu.x0*integrate(func,0,Inf,mu.t0=mu.t0,mu.x0=mu.x0)$value
        }
        beta=F.Z1(ucl.wg1[i,j])
        ATSxt[k,l]=(deltat[k]*mu.t0)/(1-beta)
      }
    }
    EATSxt.wg1[i,j]=mean(mean(ATSxt)) 
  }
}

#weibull-lognormal
EATSxt.wl1=matrix(0, nrow=3, ncol=3, dimnames = list(c("1","2","5"),c("1","2","5")))

for (i in 1:3){
  for (j in 1:3) {  
    for (k in 1:10){
      for (l in 1:10){
        F.Z1=function(z){
          func=function(x,mu.t0,mu.x0){
            pweibull((x-z)*mu.t0,shape=a1[i,k],scale=deltat[k]*mu.t0/gamma((1/a1[i,k])+1))*dlnorm(x*mu.x0,meanlog=log((deltax[l]*mu.x0)^2/sqrt(s0[j]^2+(deltax[l]*mu.x0)^2)),sdlog=sqrt(log(1+(s0[j]^2/(deltax[l]*mu.x0)^2))))
          }
          1-mu.x0*integrate(func,0,Inf,mu.t0=mu.t0,mu.x0=mu.x0)$value
        }
        beta=F.Z1(ucl.wl1[i,j])
        ATSxt[k,l]=(deltat[k]*mu.t0)/(1-beta)
      }
    }
    EATSxt.wl1[i,j]=mean(mean(ATSxt)) 
  }
}

#weibull-normal
EATSxt.wn1=matrix(0, nrow=3, ncol=2, dimnames = list(c("1","2","5"),c("1","2")))

for (i in 1:3){
  for (j in 1:2) {  
    for (k in 1:10){
      for (l in 1:10){
        F.Z1=function(z){
          func=function(x,mu.t0,mu.x0){
            pweibull((x-z)*mu.t0,shape=a1[i,k],scale=deltat[k]*mu.t0/gamma((1/a1[i,k])+1))*dnorm(x*mu.x0,mean=deltax[l]*mu.x0,sd=s0[j])
          }
          1-mu.x0*integrate(func,0,Inf,mu.t0=mu.t0,mu.x0=mu.x0)$value
        }
        beta=F.Z1(ucl.wn1[i,j])
        ATSxt[k,l]=(deltat[k]*mu.t0)/(1-beta)
      }
    }
    EATSxt.wn1[i,j]=mean(mean(ATSxt)) 
  }
}

#weibull-weibull
EATSxt.ww1=matrix(0, nrow=3, ncol=3, dimnames = list(c("1","2","5"),c("1","2","5")))

for (i in 1:3){
  for (j in 1:3) {  
    for (k in 1:10){
      for (l in 1:10){
        F.Z1=function(z){
          func=function(x,mu.t0,mu.x0){
            pweibull((x-z)*mu.t0,shape=a1t[i,k],scale=deltat[k]*mu.t0/gamma((1/a1t[i,k])+1))*dweibull(x*mu.x0,shape=a1x[j,l],scale=(deltax[l]*mu.x0)/gamma((1/a1x[j,l])+1))
          }
          1-mu.x0*integrate(func,0,Inf,mu.t0=mu.t0,mu.x0=mu.x0)$value
        }
        beta=F.Z1(ucl.ww1[i,j])
        ATSxt[k,l]=(deltat[k]*mu.t0)/(1-beta)
      }
    }
    EATSxt.ww1[i,j]=mean(mean(ATSxt)) 
  }
}
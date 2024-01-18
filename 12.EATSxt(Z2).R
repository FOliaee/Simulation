#*****EATSxt values for Statictic Z2
#gamma-gamma
EATSxt.gg2=matrix(0, nrow=3, ncol=3, dimnames = list(c("1","2","5"),c("1","2","5")))

for (i in 1:3){
  for (j in 1:3) {  
    for (k in 1:10){
      for (l in 1:10){
        F.Z2=function(z){
          func=function(x,mu.t0,mu.x0){
            pgamma((x/z)*mu.t0,shape=(deltat[k]*mu.t0/s0[i])^2,scale=(s0[i]^2/(deltat[k]*mu.t0)))*dgamma(x*mu.x0,shape=(deltax[l]*mu.x0/s0[j])^2,scale=(s0[j]^2/(deltax[l]*mu.x0)))
          }
          1-mu.x0*integrate(func,0,Inf,mu.t0=mu.t0,mu.x0=mu.x0)$value
        }
        beta=F.Z2(ucl.gg2[i,j])
        ATSxt[k,l]=(deltat[k]*mu.t0)/(1-beta)
      }
    }
    EATSxt.gg2[i,j]=mean(mean(ATSxt)) 
  }
}

#gamma-lognormal
EATSxt.gl2=matrix(0, nrow=3, ncol=3, dimnames = list(c("1","2","5"),c("1","2","5")))

for (i in 1:3){
  for (j in 1:3) {  
    for (k in 1:10){
      for (l in 1:10){
        F.Z2=function(z){
          func=function(x,mu.t0,mu.x0){
            pgamma((x/z)*mu.t0,shape=(deltat[k]*mu.t0/s0[i])^2,scale=(s0[i]^2/(deltat[k]*mu.t0)))*dlnorm(x*mu.x0,meanlog=log((deltax[l]*mu.x0)^2/sqrt(s0[j]^2+(deltax[l]*mu.x0)^2)),sdlog=sqrt(log(1+(s0[j]^2/(deltax[l]*mu.x0)^2))))
          }
          1-mu.x0*integrate(func,0,Inf,mu.t0=mu.t0,mu.x0=mu.x0)$value
        }
        beta=F.Z2(ucl.gl2[i,j])
        ATSxt[k,l]=(deltat[k]*mu.t0)/(1-beta)
      }
    }
    EATSxt.gl2[i,j]=mean(mean(ATSxt)) 
  }
}

#gamma-normal
EATSxt.gn2=matrix(0, nrow=3, ncol=2, dimnames = list(c("1","2","5"),c("1","2")))

for (i in 1:3){
  for (j in 1:2) {  
    for (k in 1:10){
      for (l in 1:10){
        F.Z2=function(z){
          func=function(x,mu.t0,mu.x0){
            pgamma((x/z)*mu.t0,shape=(deltat[k]*mu.t0/s0[i])^2,scale=(s0[i]^2/(deltat[k]*mu.t0)))*dnorm(x*mu.x0,mean=deltax[l]*mu.x0,sd=s0[j])
          }
          1-mu.x0*integrate(func,0,Inf,mu.t0=mu.t0,mu.x0=mu.x0)$value
        }
        beta=F.Z2(ucl.gn2[i,j])
        ATSxt[k,l]=(deltat[k]*mu.t0)/(1-beta)
      }
    }
    EATSxt.gn2[i,j]=mean(mean(ATSxt)) 
  }
}

#gamma-weibull
EATSxt.gw2=matrix(0, nrow=3, ncol=3, dimnames = list(c("1","2","5"),c("1","2","5")))

for (i in 1:3){
  for (j in 1:3) {  
    for (k in 1:10){
      for (l in 1:10){
        F.Z2=function(z){
          func=function(x,mu.t0,mu.x0){
            pgamma((x/z)*mu.t0,shape=(deltat[k]*mu.t0/s0[i])^2,scale=(s0[i]^2/(deltat[k]*mu.t0)))*dweibull(x*mu.x0,shape=a1[j,l],scale=(deltax[l]*mu.x0)/gamma((1/a1[j,l])+1))
          }
          1-mu.x0*integrate(func,0,Inf,mu.t0=mu.t0,mu.x0=mu.x0)$value
        }
        beta=F.Z2(ucl.gw2[i,j])
        ATSxt[k,l]=(deltat[k]*mu.t0)/(1-beta)
      }
    }
    EATSxt.gw2[i,j]=mean(mean(ATSxt)) 
  }
}

#lognormal-gamma
EATSxt.lg2=matrix(0, nrow=3, ncol=3, dimnames = list(c("1","2","5"),c("1","2","5")))

for (i in 1:3){
  for (j in 1:3) {  
    for (k in 1:10){
      for (l in 1:10){
        F.Z2=function(z){
          func=function(x,mu.t0,mu.x0){
            plnorm((x/z)*mu.t0,meanlog=log((deltat[k]*mu.t0)^2/sqrt(s0[i]^2+(deltat[k]*mu.t0)^2)),sdlog=sqrt(log(1+(s0[i]^2/(deltat[k]*mu.t0)^2))))*dgamma(x*mu.x0,shape=(deltax[l]*mu.x0/s0[j])^2,scale=(s0[j]^2/(deltax[l]*mu.x0)))
          }
          1-mu.x0*integrate(func,0,Inf,mu.t0=mu.t0,mu.x0=mu.x0)$value
        }
        beta=F.Z2(ucl.lg2[i,j])
        ATSxt[k,l]=(deltat[k]*mu.t0)/(1-beta)
      }
    }
    EATSxt.lg2[i,j]=mean(mean(ATSxt)) 
  }
}

#lognormal-lognormal
EATSxt.ll2=matrix(0, nrow=3, ncol=3, dimnames = list(c("1","2","5"),c("1","2","5")))

for (i in 1:3){
  for (j in 1:3) {  
    for (k in 1:10){
      for (l in 1:10){
        F.Z2=function(z){
          func=function(x,mu.t0,mu.x0){
            plnorm((x/z)*mu.t0,meanlog=log((deltat[k]*mu.t0)^2/sqrt(s0[i]^2+(deltat[k]*mu.t0)^2)),sdlog=sqrt(log(1+(s0[i]^2/(deltat[k]*mu.t0)^2))))*dlnorm(x*mu.x0,meanlog=log((deltax[l]*mu.x0)^2/sqrt(s0[j]^2+(deltax[l]*mu.x0)^2)),sdlog=sqrt(log(1+(s0[j]^2/(deltax[l]*mu.x0)^2))))
          }
          1-mu.x0*integrate(func,0,Inf,mu.t0=mu.t0,mu.x0=mu.x0)$value
        }
        beta=F.Z2(ucl.ll2[i,j])
        ATSxt[k,l]=(deltat[k]*mu.t0)/(1-beta)
      }
    }
    EATSxt.ll2[i,j]=mean(mean(ATSxt)) 
  }
}

#lognormal-normal
EATSxt.ln2=matrix(0, nrow=3, ncol=2, dimnames = list(c("1","2","5"),c("1","2")))

for (i in 1:3){
  for (j in 1:2) {  
    for (k in 1:10){
      for (l in 1:10){
        F.Z2=function(z){
          func=function(x,mu.t0,mu.x0){
            plnorm((x/z)*mu.t0,meanlog=log((deltat[k]*mu.t0)^2/sqrt(s0[i]^2+(deltat[k]*mu.t0)^2)),sdlog=sqrt(log(1+(s0[i]^2/(deltat[k]*mu.t0)^2))))*dnorm(x*mu.x0,mean=deltax[l]*mu.x0,sd=s0[j])
          }
          1-mu.x0*integrate(func,0,Inf,mu.t0=mu.t0,mu.x0=mu.x0)$value
        }
        beta=F.Z2(ucl.ln2[i,j])
        ATSxt[k,l]=(deltat[k]*mu.t0)/(1-beta)
      }
    }
    EATSxt.ln2[i,j]=mean(mean(ATSxt)) 
  }
}

#lognormal-weibull
EATSxt.lw2=matrix(0, nrow=3, ncol=3, dimnames = list(c("1","2","5"),c("1","2","5")))

for (i in 1:3){
  for (j in 1:3) {  
    for (k in 1:10){
      for (l in 1:10){
        F.Z2=function(z){
          func=function(x,mu.t0,mu.x0){
            plnorm((x/z)*mu.t0,meanlog=log((deltat[k]*mu.t0)^2/sqrt(s0[i]^2+(deltat[k]*mu.t0)^2)),sdlog=sqrt(log(1+(s0[i]^2/(deltat[k]*mu.t0)^2))))*dweibull(x*mu.x0,shape=a1[j,l],scale=(deltax[l]*mu.x0)/gamma((1/a1[j,l])+1))
          }
          1-mu.x0*integrate(func,0,Inf,mu.t0=mu.t0,mu.x0=mu.x0)$value
        }
        beta=F.Z2(ucl.lw2[i,j])
        ATSxt[k,l]=(deltat[k]*mu.t0)/(1-beta)
      }
    }
    EATSxt.lw2[i,j]=mean(mean(ATSxt)) 
  }
}

#weibull-gamma
EATSxt.wg2=matrix(0, nrow=3, ncol=3, dimnames = list(c("1","2","5"),c("1","2","5")))

for (i in 1:3){
  for (j in 1:3) {  
    for (k in 1:10){
      for (l in 1:10){
        F.Z2=function(z){
          func=function(x,mu.t0,mu.x0){
            pweibull((x/z)*mu.t0,shape=a1[i,k],scale=deltat[k]*mu.t0/gamma((1/a1[i,k])+1))*dgamma(x*mu.x0,shape=(deltax[l]*mu.x0/s0[j])^2,scale=(s0[j]^2/(deltax[l]*mu.x0)))
          }
          1-mu.x0*integrate(func,0,Inf,mu.t0=mu.t0,mu.x0=mu.x0)$value
        }
        beta=F.Z2(ucl.wg2[i,j])
        ATSxt[k,l]=(deltat[k]*mu.t0)/(1-beta)
      }
    }
    EATSxt.wg2[i,j]=mean(mean(ATSxt)) 
  }
}

#weibull-lognormal
EATSxt.wl2=matrix(0, nrow=3, ncol=3, dimnames = list(c("1","2","5"),c("1","2","5")))

for (i in 1:3){
  for (j in 1:3) {  
    for (k in 1:10){
      for (l in 1:10){
        F.Z2=function(z){
          func=function(x,mu.t0,mu.x0){
            pweibull((x/z)*mu.t0,shape=a1[i,k],scale=deltat[k]*mu.t0/gamma((1/a1[i,k])+1))*dlnorm(x*mu.x0,meanlog=log((deltax[l]*mu.x0)^2/sqrt(s0[j]^2+(deltax[l]*mu.x0)^2)),sdlog=sqrt(log(1+(s0[j]^2/(deltax[l]*mu.x0)^2))))
          }
          1-mu.x0*integrate(func,0,Inf,mu.t0=mu.t0,mu.x0=mu.x0)$value
        }
        beta=F.Z2(ucl.wl2[i,j])
        ATSxt[k,l]=(deltat[k]*mu.t0)/(1-beta)
      }
    }
    EATSxt.wl2[i,j]=mean(mean(ATSxt)) 
  }
}

#weibull-normal
EATSxt.wn2=matrix(0, nrow=3, ncol=2, dimnames = list(c("1","2","5"),c("1","2")))

for (i in 1:3){
  for (j in 1:2) {  
    for (k in 1:10){
      for (l in 1:10){
        F.Z2=function(z){
          func=function(x,mu.t0,mu.x0){
            pweibull((x/z)*mu.t0,shape=a1[i,k],scale=deltat[k]*mu.t0/gamma((1/a1[i,k])+1))*dnorm(x*mu.x0,mean=deltax[l]*mu.x0,sd=s0[j])
          }
          1-mu.x0*integrate(func,0,Inf,mu.t0=mu.t0,mu.x0=mu.x0)$value
        }
        beta=F.Z2(ucl.wn2[i,j])
        ATSxt[k,l]=(deltat[k]*mu.t0)/(1-beta)
      }
    }
    EATSxt.wn2[i,j]=mean(mean(ATSxt)) 
  }
}

#weibull-weibull
EATSxt.ww2=matrix(0, nrow=3, ncol=3, dimnames = list(c("1","2","5"),c("1","2","5")))

for (i in 1:3){
  for (j in 1:3) {  
    for (k in 1:10){
      for (l in 1:10){
        F.Z2=function(z){
          func=function(x,mu.t0,mu.x0){
            pweibull((x/z)*mu.t0,shape=a1t[i,k],scale=deltat[k]*mu.t0/gamma((1/a1t[i,k])+1))*dweibull(x*mu.x0,shape=a1x[j,l],scale=(deltax[l]*mu.x0)/gamma((1/a1x[j,l])+1))
          }
          1-mu.x0*integrate(func,0,Inf,mu.t0=mu.t0,mu.x0=mu.x0)$value
        }
        beta=F.Z2(ucl.ww2[i,j])
        ATSxt[k,l]=(deltat[k]*mu.t0)/(1-beta)
      }
    }
    EATSxt.ww2[i,j]=mean(mean(ATSxt)) 
  }
}
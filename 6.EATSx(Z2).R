#*****EATSx values for Statictic Z2
#gamma-gamma
EATSx.gg2=matrix(0, nrow=3, ncol=3, dimnames = list(c("1","2","5"),c("1","2","5")))

for(i in 1:3){
  for(j in 1:3){
    for (k in 1:10){
      F.Z2=function(z){
        func=function(x,mu.t0,mu.x0){
          pgamma((x/z)*mu.t0,shape=(mu.t0/s0[i])^2,scale=(s0[i]^2/mu.t0))*dgamma(x*mu.x0,shape=(deltax[k]*mu.x0/s0[j])^2,scale=(s0[j]^2/(deltax[k]*mu.x0)))
        }
        1-mu.x0*integrate(func,0,Inf,mu.t0=mu.t0,mu.x0=mu.x0)$value
      }
      beta=F.Z2(ucl.gg2[i,j])
      ATSx[k]=mu.t0/(1-beta)
    }
    EATSx.gg2[i,j]= mean(ATSx)
  }
}
    
#gamma-lognormal
EATSx.gl2=matrix(0, nrow=3, ncol=3, dimnames = list(c("1","2","5"),c("1","2","5")))

for(i in 1:3){
  for(j in 1:3){
    for (k in 1:10){
      F.Z2=function(z){
        func=function(x,mu.t0,mu.x0){
          pgamma((x/z)*mu.t0,shape=(mu.t0/s0[i])^2,scale=(s0[i]^2/mu.t0))*dlnorm(x*mu.x0,meanlog=log((deltax[k]*mu.x0)^2/sqrt(s0[j]^2+(deltax[k]*mu.x0)^2)),sdlog=sqrt(log(1+(s0[j]^2/(deltax[k]*mu.x0)^2))))
        }
        1-mu.x0*integrate(func,0,Inf,mu.t0=mu.t0,mu.x0=mu.x0)$value
      }
      beta=F.Z2(ucl.gl2[i,j])
      ATSx[k]=mu.t0/(1-beta)
    }
    EATSx.gl2[i,j]= mean(ATSx)
  }
}

#gamma-normal
EATSx.gn2=matrix(0, nrow=3, ncol=2, dimnames = list(c("1","2","5"),c("1","2")))

for(i in 1:3){
  for(j in 1:2){
    for (k in 1:10){
      F.Z2=function(z){
        func=function(x,mu.t0,mu.x0){
          pgamma((x/z)*mu.t0,shape=(mu.t0/s0[i])^2,scale=(s0[i]^2/mu.t0))*dnorm(x*mu.x0,mean=deltax[k]*mu.x0,sd=s0[j])
        }
        1-mu.x0*integrate(func,0,Inf,mu.t0=mu.t0,mu.x0=mu.x0)$value
      }
      beta=F.Z2(ucl.gn2[i,j])
      ATSx[k]=mu.t0/(1-beta)
    }
    EATSx.gn2[i,j]= mean(ATSx)
  }
}

#gamma-weibull
EATSx.gw2=matrix(0, nrow=3, ncol=3, dimnames = list(c("1","2","5"),c("1","2","5")))

for(i in 1:3){
  for(j in 1:3){
    for (k in 1:10){
      F.Z2=function(z){
        func=function(x,mu.t0,mu.x0){
          pgamma((x/z)*mu.t0,shape=(mu.t0/s0[i])^2,scale=(s0[i]^2/mu.t0))*dweibull(x*mu.x0,shape=a1[j,k],scale=(deltax[k]*mu.x0)/gamma((1/a1[j,k])+1))
        }
        1-mu.x0*integrate(func,0,Inf,mu.t0=mu.t0,mu.x0=mu.x0)$value
      }
      beta=F.Z2(ucl.gw2[i,j])
      ATSx[k]=mu.t0/(1-beta)
    }
    EATSx.gw2[i,j]= mean(ATSx)
  }
}

#lognormal-gamma
EATSx.lg2=matrix(0, nrow=3, ncol=3, dimnames = list(c("1","2","5"),c("1","2","5")))

for(i in 1:3){
  for(j in 1:3){
    for (k in 1:10){
      F.Z2=function(z){
        func=function(x,mu.t0,mu.x0){
          plnorm((x/z)*mu.t0,meanlog=log(mu.t0^2/sqrt(s0[i]^2+mu.t0^2)),sdlog=sqrt(log(1+(s0[i]^2/mu.t0^2))))*dgamma(x*mu.x0,shape=(deltax[k]*mu.x0/s0[j])^2,scale=(s0[j]^2/(deltax[k]*mu.x0)))
        }
        1-mu.x0*integrate(func,0,Inf,mu.t0=mu.t0,mu.x0=mu.x0)$value
      }
      beta=F.Z2(ucl.lg2[i,j])
      ATSx[k]=mu.t0/(1-beta)
    }
    EATSx.lg2[i,j]= mean(ATSx)
  }
}

#lognormal-lognormal
EATSx.ll2=matrix(0, nrow=3, ncol=3, dimnames = list(c("1","2","5"),c("1","2","5")))

for(i in 1:3){
  for(j in 1:3){
    for (k in 1:10){
      F.Z1=function(z){
        func=function(x,mu.t0,mu.x0){
          plnorm((x/z)*mu.t0,meanlog=log(mu.t0^2/sqrt(s0[i]^2+mu.t0^2)),sdlog=sqrt(log(1+(s0[i]^2/mu.t0^2))))*dlnorm(x*mu.x0,meanlog=log((deltax[k]*mu.x0)^2/sqrt(s0[j]^2+(deltax[k]*mu.x0)^2)),sdlog=sqrt(log(1+(s0[j]^2/(deltax[k]*mu.x0)^2))))
        }
        1-mu.x0*integrate(func,0,Inf,mu.t0=mu.t0,mu.x0=mu.x0)$value
      }
      beta=F.Z1(ucl.ll2[i,j])
      ATSx[k]=mu.t0/(1-beta)
    }
    EATSx.ll2[i,j]= mean(ATSx)
  }
}

#lognormal-normal
EATSx.ln2=matrix(0, nrow=3, ncol=2, dimnames = list(c("1","2","5"),c("1","2")))

for(i in 1:3){
  for(j in 1:2){
    for (k in 1:10){
      F.Z2=function(z){
        func=function(x,mu.t0,mu.x0){
          plnorm((x/z)*mu.t0,meanlog=log(mu.t0^2/sqrt(s0[i]^2+mu.t0^2)),sdlog=sqrt(log(1+(s0[i]^2/mu.t0^2))))*dnorm(x*mu.x0,mean=deltax[k]*mu.x0,sd=s0[j])
        }
        1-mu.x0*integrate(func,0,Inf,mu.t0=mu.t0,mu.x0=mu.x0)$value
      }
      beta=F.Z2(ucl.ln2[i,j])
      ATSx[k]=mu.t0/(1-beta)
    }
    EATSx.ln2[i,j]= mean(ATSx)
  }
}

#lognormal-weibull
EATSx.lw2=matrix(0, nrow=3, ncol=3, dimnames = list(c("1","2","5"),c("1","2","5")))

for(i in 1:3){
  for(j in 1:3){
    for (k in 1:10){
      F.Z2=function(z){
        func=function(x,mu.t0,mu.x0){
          plnorm((x/z)*mu.t0,meanlog=log(mu.t0^2/sqrt(s0[i]^2+mu.t0^2)),sdlog=sqrt(log(1+(s0[i]^2/mu.t0^2))))*dweibull(x*mu.x0,shape=a1[j,k],scale=(deltax[k]*mu.x0)/gamma((1/a1[j,k])+1))
        }
        1-mu.x0*integrate(func,0,Inf,mu.t0=mu.t0,mu.x0=mu.x0)$value
      }
      beta=F.Z2(ucl.lw2[i,j])
      ATSx[k]=mu.t0/(1-beta)
    }
    EATSx.lw2[i,j]= mean(ATSx)
  }
}  

#weibull-gamma
EATSx.wg2=matrix(0, nrow=3, ncol=3, dimnames = list(c("1","2","5"),c("1","2","5")))

for(i in 1:3){
  for(j in 1:3){
    for (k in 1:10){
      F.Z2=function(z){
        func=function(x,mu.t0,mu.x0){
          pweibull((x/z)*mu.t0,shape=a0[i],scale=mu.t0/gamma((1/a0[i])+1))*dgamma(x*mu.x0,shape=(deltax[k]*mu.x0/s0[j])^2,scale=(s0[j]^2/(deltax[k]*mu.x0)))
        }
        1-mu.x0*integrate(func,0,Inf,mu.t0=mu.t0,mu.x0=mu.x0)$value
      }
      beta=F.Z2(ucl.wg2[i,j])
      ATSx[k]=mu.t0/(1-beta)
    }
    EATSx.wg2[i,j]= mean(ATSx)
  }
}

#weibull-lognormal
EATSx.wl2=matrix(0, nrow=3, ncol=3, dimnames = list(c("1","2","5"),c("1","2","5")))

for(i in 1:3){
  for(j in 1:3){
    for (k in 1:10){
      F.Z2=function(z){
        func=function(x,mu.t0,mu.x0){
          pweibull((x/z)*mu.t0,shape=a0[i],scale=mu.t0/gamma((1/a0[i])+1))*dlnorm(x*mu.x0,meanlog=log((deltax[k]*mu.x0)^2/sqrt(s0[j]^2+(deltax[k]*mu.x0)^2)),sdlog=sqrt(log(1+(s0[j]^2/(deltax[k]*mu.x0)^2))))
        }
        1-mu.x0*integrate(func,0,Inf,mu.t0=mu.t0,mu.x0=mu.x0)$value
      }
      beta=F.Z2(ucl.wl2[i,j])
      ATSx[k]=mu.t0/(1-beta)
    }
    EATSx.wl2[i,j]= mean(ATSx)
  }
}

#weibull-normal
EATSx.wn2=matrix(0, nrow=3, ncol=2, dimnames = list(c("1","2","5"),c("1","2")))

for(i in 1:3){
  for(j in 1:2){
    for (k in 1:10){
      F.Z2=function(z){
        func=function(x,mu.t0,mu.x0){
          pweibull((x/z)*mu.t0,shape=a0[i],scale=mu.t0/gamma((1/a0[i])+1))*dnorm(x*mu.x0,mean=deltax[k]*mu.x0,sd=s0[j])
        }
        1-mu.x0*integrate(func,0,Inf,mu.t0=mu.t0,mu.x0=mu.x0)$value
      }
      beta=F.Z2(ucl.wn2[i,j])
      ATSx[k]=mu.t0/(1-beta)
    }
    EATSx.wn2[i,j]= mean(ATSx)
  }
}   

#weibull-weibull
EATSx.ww2=matrix(0, nrow=3, ncol=3, dimnames = list(c("1","2","5"),c("1","2","5")))

for(i in 1:3){
  for(j in 1:3){
    for (k in 1:10){
      F.Z2=function(z){
        func=function(x,mu.t0,mu.x0){
          pweibull((x/z)*mu.t0,shape=a0[i],scale=mu.t0/gamma((1/a0[i])+1))*dweibull(x*mu.x0,shape=a1[j,k],scale=(deltax[k]*mu.x0)/gamma((1/a1[j,k])+1))
        }
        1-mu.x0*integrate(func,0,Inf,mu.t0=mu.t0,mu.x0=mu.x0)$value
      }
      beta=F.Z2(ucl.ww2[i,j])
      ATSx[k]=mu.t0/(1-beta)
    }
    EATSx.ww2[i,j]= mean(ATSx)
  }
}
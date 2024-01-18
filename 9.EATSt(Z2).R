#*****EATSt values for Statictic Z2
#gamma-gamma
EATSt.gg2=matrix(0, nrow=3, ncol=3, dimnames = list(c("1","2","5"),c("1","2","5")))

for (i in 1:3){
  for (j in 1:3) {  
    for (k in 1:10){
      F.Z2=function(z){
        func=function(x,mu.t0,mu.x0){
          pgamma((x/z)*mu.t0,shape=(deltat[k]*mu.t0/s0[i])^2,scale=(s0[i]^2/(deltat[k]*mu.t0)))*dgamma(x*mu.x0,shape=(mu.x0/s0[j])^2,scale=(s0[j]^2/mu.x0))
        }
        1-mu.x0*integrate(func,0,Inf,mu.t0=mu.t0,mu.x0=mu.x0)$value
      }
      beta=F.Z2(ucl.gg2[i,j])
      ATSt[k]=(deltat[k]*mu.t0)/(1-beta)
    }
    EATSt.gg2[i,j]= mean(ATSt)
  }
}

#gamma-lognormal
EATSt.gl2=matrix(0, nrow=3, ncol=3, dimnames = list(c("1","2","5"),c("1","2","5")))

for (i in 1:3){
  for (j in 1:3) {  
    for (k in 1:10){
      F.Z2=function(z){
        func=function(x,mu.t0,mu.x0){
          pgamma((x/z)*mu.t0,shape=(deltat[k]*mu.t0/s0[i])^2,scale=(s0[i]^2/(deltat[k]*mu.t0)))*dlnorm(x*mu.x0,meanlog=log(mu.x0^2/sqrt(s0[j]^2+mu.x0^2)),sdlog=sqrt(log(1+(s0[j]^2/mu.x0^2))))
        }
        1-mu.x0*integrate(func,0,Inf,mu.t0=mu.t0,mu.x0=mu.x0)$value
      }
      beta=F.Z2(ucl.gl2[i,j])
      ATSt[k]=(deltat[k]*mu.t0)/(1-beta)
    }
    EATSt.gl2[i,j]= mean(ATSt)
  }
}

#gamma-normal
EATSt.gn2=matrix(0, nrow=3, ncol=2, dimnames = list(c("1","2","5"),c("1","2")))

for (i in 1:3){
  for (j in 1:2) {  
    for (k in 1:10){
      F.Z2=function(z){
        func=function(x,mu.t0,mu.x0){
          pgamma((x/z)*mu.t0,shape=(deltat[k]*mu.t0/s0[i])^2,scale=(s0[i]^2/(deltat[k]*mu.t0)))*dnorm(x*mu.x0,mean=mu.x0,sd=s0[j])
        }
        1-mu.x0*integrate(func,0,Inf,mu.t0=mu.t0,mu.x0=mu.x0)$value
      }
      beta=F.Z2(ucl.gn2[i,j])
      ATSt[k]=(deltat[k]*mu.t0)/(1-beta)
    }
    EATSt.gn2[i,j]= mean(ATSt)
  }
}

#gamma-weibull
EATSt.gw2=matrix(0, nrow=3, ncol=3, dimnames = list(c("1","2","5"),c("1","2","5")))

for (i in 1:3){
  for (j in 1:3) {  
    for (k in 1:10){
      F.Z2=function(z){
        func=function(x,mu.t0,mu.x0){
          pgamma((x/z)*mu.t0,shape=(deltat[k]*mu.t0/s0[i])^2,scale=(s0[i]^2/(deltat[k]*mu.t0)))*dweibull(x*mu.x0,shape=a0[j],scale=mu.x0/gamma((1/a0[j])+1))
        }
        1-mu.x0*integrate(func,0,Inf,mu.t0=mu.t0,mu.x0=mu.x0)$value
      }
      beta=F.Z2(ucl.gw2[i,j])
      ATSt[k]=(deltat[k]*mu.t0)/(1-beta)
    }
    EATSt.gw2[i,j]= mean(ATSt)
  }
}

#lognormal-gamma
EATSt.lg2=matrix(0, nrow=3, ncol=3, dimnames = list(c("1","2","5"),c("1","2","5")))

for (i in 1:3){
  for (j in 1:3) {  
    for (k in 1:10){
      F.Z2=function(z){
        func=function(x,mu.t0,mu.x0){
          plnorm((x/z)*mu.t0,meanlog=log((deltat[k]*mu.t0)^2/sqrt(s0[i]^2+(deltat[k]*mu.t0)^2)),sdlog=sqrt(log(1+(s0[i]^2/(deltat[k]*mu.t0)^2))))*dgamma(x*mu.x0,shape=(mu.x0/s0[j])^2,scale=(s0[j]^2/mu.x0))
        }
        1-mu.x0*integrate(func,0,Inf,mu.t0=mu.t0,mu.x0=mu.x0)$value
      }
      beta=F.Z2(ucl.lg2[i,j])
      ATSt[k]=(deltat[k]*mu.t0)/(1-beta)
    }
    EATSt.lg2[i,j]= mean(ATSt)
  }
}

#lognormal-lognormal
EATSt.ll2=matrix(0, nrow=3, ncol=3, dimnames = list(c("1","2","5"),c("1","2","5")))

for (i in 1:3){
  for (j in 1:3) {  
    for (k in 1:10){
      F.Z2=function(z){
        func=function(x,mu.t0,mu.x0){
          plnorm((x/z)*mu.t0,meanlog=log((deltat[k]*mu.t0)^2/sqrt(s0[i]^2+(deltat[k]*mu.t0)^2)),sdlog=sqrt(log(1+(s0[i]^2/(deltat[k]*mu.t0)^2))))*dlnorm(x*mu.x0,meanlog=log(mu.x0^2/sqrt(s0[j]^2+mu.x0^2)),sdlog=sqrt(log(1+(s0[j]^2/mu.x0^2))))
        }
        1-mu.x0*integrate(func,0,Inf,mu.t0=mu.t0,mu.x0=mu.x0)$value
      }
      beta=F.Z2(ucl.ll2[i,j])
      ATSt[k]=(deltat[k]*mu.t0)/(1-beta)
    }
    EATSt.ll2[i,j]= mean(ATSt)
  }
}

#lognormal-normal
EATSt.ln2=matrix(0, nrow=3, ncol=2, dimnames = list(c("1","2","5"),c("1","2")))

for (i in 1:3){
  for (j in 1:2) {  
    for (k in 1:10){
      F.Z2=function(z){
        func=function(x,mu.t0,mu.x0){
          plnorm((x/z)*mu.t0,meanlog=log((deltat[k]*mu.t0)^2/sqrt(s0[i]^2+(deltat[k]*mu.t0)^2)),sdlog=sqrt(log(1+(s0[i]^2/(deltat[k]*mu.t0)^2))))*dnorm(x*mu.x0,mean=mu.x0,sd=s0[j])
        }
        1-mu.x0*integrate(func,0,Inf,mu.t0=mu.t0,mu.x0=mu.x0)$value
      }
      beta=F.Z2(ucl.ln2[i,j])
      ATSt[k]=(deltat[k]*mu.t0)/(1-beta)
    }
    EATSt.ln2[i,j]= mean(ATSt)
  }
}

#lognormal-weibull
EATSt.lw2=matrix(0, nrow=3, ncol=3, dimnames = list(c("1","2","5"),c("1","2","5")))

for (i in 1:3){
  for (j in 1:3) {  
    for (k in 1:10){
      F.Z2=function(z){
        func=function(x,mu.t0,mu.x0){
          plnorm((x/z)*mu.t0,meanlog=log((deltat[k]*mu.t0)^2/sqrt(s0[i]^2+(deltat[k]*mu.t0)^2)),sdlog=sqrt(log(1+(s0[i]^2/(deltat[k]*mu.t0)^2))))*dweibull(x*mu.x0,shape=a0[j],scale=mu.x0/gamma((1/a0[j])+1))
        }
        1-mu.x0*integrate(func,0,Inf,mu.t0=mu.t0,mu.x0=mu.x0)$value
      }
      beta=F.Z2(ucl.lw2[i,j])
      ATSt[k]=(deltat[k]*mu.t0)/(1-beta)
    }
    EATSt.lw2[i,j]= mean(ATSt)
  }
}

#weibull-gamma
EATSt.wg2=matrix(0, nrow=3, ncol=3, dimnames = list(c("1","2","5"),c("1","2","5")))

for (i in 1:3){
  for (j in 1:3) {  
    for (k in 1:10){
      F.Z2=function(z){
        func=function(x,mu.t0,mu.x0){
          pweibull((x/z)*mu.t0,shape=a1[i,k],scale=deltat[k]*mu.t0/gamma((1/a1[i,k])+1))*dgamma(x*mu.x0,shape=(mu.x0/s0[j])^2,scale=(s0[j]^2/mu.x0))
        }
        1-mu.x0*integrate(func,0,Inf,mu.t0=mu.t0,mu.x0=mu.x0)$value
      }
      beta=F.Z2(ucl.wg2[i,j])
      ATSt[k]=(deltat[k]*mu.t0)/(1-beta)
    }
    EATSt.wg2[i,j]= mean(ATSt)
  }
}

#weibull-lognormal
EATSt.wl2=matrix(0, nrow=3, ncol=3, dimnames = list(c("1","2","5"),c("1","2","5")))

for (i in 1:3){
  for (j in 1:3) {  
    for (k in 1:10){
      F.Z2=function(z){
        func=function(x,mu.t0,mu.x0){
          pweibull((x/z)*mu.t0,shape=a1[i,k],scale=deltat[k]*mu.t0/gamma((1/a1[i,k])+1))*dlnorm(x*mu.x0,meanlog=log(mu.x0^2/sqrt(s0[j]^2+mu.x0^2)),sdlog=sqrt(log(1+(s0[j]^2/mu.x0^2))))
        }
        1-mu.x0*integrate(func,0,Inf,mu.t0=mu.t0,mu.x0=mu.x0)$value
      }
      beta=F.Z2(ucl.wl2[i,j])
      ATSt[k]=(deltat[k]*mu.t0)/(1-beta)
    }
    EATSt.wl2[i,j]= mean(ATSt)
  }
}

#weibull-normal
EATSt.wn2=matrix(0, nrow=3, ncol=2, dimnames = list(c("1","2","5"),c("1","2")))

for (i in 1:3){
  for (j in 1:2) {  
    for (k in 1:10){
      F.Z2=function(z){
        func=function(x,mu.t0,mu.x0){
          pweibull((x/z)*mu.t0,shape=a1[i,k],scale=deltat[k]*mu.t0/gamma((1/a1[i,k])+1))*dnorm(x*mu.x0,mean=mu.x0,sd=s0[j])
        }
        1-mu.x0*integrate(func,0,Inf,mu.t0=mu.t0,mu.x0=mu.x0)$value
      }
      beta=F.Z2(ucl.wn2[i,j])
      ATSt[k]=(deltat[k]*mu.t0)/(1-beta)
    }
    EATSt.wn2[i,j]= mean(ATSt)
  }
}

#weibull-weibull
EATSt.ww2=matrix(0, nrow=3, ncol=3, dimnames = list(c("1","2","5"),c("1","2","5")))

for (i in 1:3){
  for (j in 1:3) {  
    for (k in 1:10){
      F.Z2=function(z){
        func=function(x,mu.t0,mu.x0){
          pweibull((x/z)*mu.t0,shape=a1[i,k],scale=deltat[k]*mu.t0/gamma((1/a1[i,k])+1))*dweibull(x*mu.x0,shape=a0[j],scale=mu.x0/gamma((1/a0[j])+1))
        }
        1-mu.x0*integrate(func,0,Inf,mu.t0=mu.t0,mu.x0=mu.x0)$value
      }
      beta=F.Z2(ucl.ww2[i,j])
      ATSt[k]=(deltat[k]*mu.t0)/(1-beta)
    }
    EATSt.ww2[i,j]= mean(ATSt)
  }
}
#*****EATSt values for Statictic Z1
#gamma-gamma
EATSt.gg1=matrix(0, nrow=3, ncol=3, dimnames = list(c("1","2","5"),c("1","2","5")))

for (i in 1:3){
  for (j in 1:3) {  
    for (k in 1:10){
      F.Z1=function(z){
        func=function(x,mu.t0,mu.x0){
          pgamma((x-z)*mu.t0,shape=(deltat[k]*mu.t0/s0[i])^2,scale=(s0[i]^2/(deltat[k]*mu.t0)))*dgamma(x*mu.x0,shape=(mu.x0/s0[j])^2,scale=(s0[j]^2/mu.x0))
        }
        1-mu.x0*integrate(func,0,Inf,mu.t0=mu.t0,mu.x0=mu.x0)$value
      }
      beta=F.Z1(ucl.gg1[i,j])
      ATSt[k]=(deltat[k]*mu.t0)/(1-beta)
    }
    EATSt.gg1[i,j]= mean(ATSt)
  }
}

#gamma-lognormal
EATSt.gl1=matrix(0, nrow=3, ncol=3, dimnames = list(c("1","2","5"),c("1","2","5")))

for (i in 1:3){
  for (j in 1:3) {  
    for (k in 1:10){
      F.Z1=function(z){
        func=function(x,mu.t0,mu.x0){
          pgamma((x-z)*mu.t0,shape=(deltat[k]*mu.t0/s0[i])^2,scale=(s0[i]^2/(deltat[k]*mu.t0)))*dlnorm(x*mu.x0,meanlog=log(mu.x0^2/sqrt(s0[j]^2+mu.x0^2)),sdlog=sqrt(log(1+(s0[j]^2/mu.x0^2))))
        }
        1-mu.x0*integrate(func,0,Inf,mu.t0=mu.t0,mu.x0=mu.x0)$value
      }
      beta=F.Z1(ucl.gl1[i,j])
      ATSt[k]=(deltat[k]*mu.t0)/(1-beta)
    }
    EATSt.gl1[i,j]= mean(ATSt)
  }
}

#gamma-normal
EATSt.gn1=matrix(0, nrow=3, ncol=2, dimnames = list(c("1","2","5"),c("1","2")))

for (i in 1:3){
  for (j in 1:2) {  
    for (k in 1:10){
      F.Z1=function(z){
        func=function(x,mu.t0,mu.x0){
          pgamma((x-z)*mu.t0,shape=(deltat[k]*mu.t0/s0[i])^2,scale=(s0[i]^2/(deltat[k]*mu.t0)))*dnorm(x*mu.x0,mean=mu.x0,sd=s0[j])
        }
        1-mu.x0*integrate(func,0,Inf,mu.t0=mu.t0,mu.x0=mu.x0)$value
      }
      beta=F.Z1(ucl.gn1[i,j])
      ATSt[k]=(deltat[k]*mu.t0)/(1-beta)
    }
    EATSt.gn1[i,j]= mean(ATSt)
  }
}

#gamma-weibull
EATSt.gw1=matrix(0, nrow=3, ncol=3, dimnames = list(c("1","2","5"),c("1","2","5")))

for (i in 1:3){
  for (j in 1:3) {  
    for (k in 1:10){
      F.Z1=function(z){
        func=function(x,mu.t0,mu.x0){
          pgamma((x-z)*mu.t0,shape=(deltat[k]*mu.t0/s0[i])^2,scale=(s0[i]^2/(deltat[k]*mu.t0)))*dweibull(x*mu.x0,shape=a0[j],scale=mu.x0/gamma((1/a0[j])+1))
        }
        1-mu.x0*integrate(func,0,Inf,mu.t0=mu.t0,mu.x0=mu.x0)$value
      }
      beta=F.Z1(ucl.gw1[i,j])
      ATSt[k]=(deltat[k]*mu.t0)/(1-beta)
    }
    EATSt.gw1[i,j]= mean(ATSt)
  }
}

#lognormal-gamma
EATSt.lg1=matrix(0, nrow=3, ncol=3, dimnames = list(c("1","2","5"),c("1","2","5")))

for (i in 1:3){
  for (j in 1:3) {  
    for (k in 1:10){
      F.Z1=function(z){
        func=function(x,mu.t0,mu.x0){
          plnorm((x-z)*mu.t0,meanlog=log((deltat[k]*mu.t0)^2/sqrt(s0[i]^2+(deltat[k]*mu.t0)^2)),sdlog=sqrt(log(1+(s0[i]^2/(deltat[k]*mu.t0)^2))))*dgamma(x*mu.x0,shape=(mu.x0/s0[j])^2,scale=(s0[j]^2/mu.x0))
        }
        1-mu.x0*integrate(func,0,Inf,mu.t0=mu.t0,mu.x0=mu.x0)$value
      }
      beta=F.Z1(ucl.lg1[i,j])
      ATSt[k]=(deltat[k]*mu.t0)/(1-beta)
    }
    EATSt.lg1[i,j]= mean(ATSt)
  }
}

#lognormal-lognormal
EATSt.ll1=matrix(0, nrow=3, ncol=3, dimnames = list(c("1","2","5"),c("1","2","5")))

for (i in 1:3){
  for (j in 1:3) {  
    for (k in 1:10){
      F.Z1=function(z){
        func=function(x,mu.t0,mu.x0){
          plnorm((x-z)*mu.t0,meanlog=log((deltat[k]*mu.t0)^2/sqrt(s0[i]^2+(deltat[k]*mu.t0)^2)),sdlog=sqrt(log(1+(s0[i]^2/(deltat[k]*mu.t0)^2))))*dlnorm(x*mu.x0,meanlog=log(mu.x0^2/sqrt(s0[j]^2+mu.x0^2)),sdlog=sqrt(log(1+(s0[j]^2/mu.x0^2))))
        }
        1-mu.x0*integrate(func,0,Inf,mu.t0=mu.t0,mu.x0=mu.x0)$value
      }
      beta=F.Z1(ucl.ll1[i,j])
      ATSt[k]=(deltat[k]*mu.t0)/(1-beta)
    }
    EATSt.ll1[i,j]= mean(ATSt)
  }
}

#lognormal-normal
EATSt.ln1=matrix(0, nrow=3, ncol=2, dimnames = list(c("1","2","5"),c("1","2")))

for (i in 1:3){
  for (j in 1:2) {  
    for (k in 1:10){
      F.Z1=function(z){
        func=function(x,mu.t0,mu.x0){
          plnorm((x-z)*mu.t0,meanlog=log((deltat[k]*mu.t0)^2/sqrt(s0[i]^2+(deltat[k]*mu.t0)^2)),sdlog=sqrt(log(1+(s0[i]^2/(deltat[k]*mu.t0)^2))))*dnorm(x*mu.x0,mean=mu.x0,sd=s0[j])
        }
        1-mu.x0*integrate(func,0,Inf,mu.t0=mu.t0,mu.x0=mu.x0)$value
      }
      beta=F.Z1(ucl.ln1[i,j])
      ATSt[k]=(deltat[k]*mu.t0)/(1-beta)
    }
    EATSt.ln1[i,j]= mean(ATSt)
  }
}

#lognormal-weibull
EATSt.lw1=matrix(0, nrow=3, ncol=3, dimnames = list(c("1","2","5"),c("1","2","5")))

for (i in 1:3){
  for (j in 1:3) {  
    for (k in 1:10){
      F.Z1=function(z){
        func=function(x,mu.t0,mu.x0){
          plnorm((x-z)*mu.t0,meanlog=log((deltat[k]*mu.t0)^2/sqrt(s0[i]^2+(deltat[k]*mu.t0)^2)),sdlog=sqrt(log(1+(s0[i]^2/(deltat[k]*mu.t0)^2))))*dweibull(x*mu.x0,shape=a0[j],scale=mu.x0/gamma((1/a0[j])+1))
        }
        1-mu.x0*integrate(func,0,Inf,mu.t0=mu.t0,mu.x0=mu.x0)$value
      }
      beta=F.Z1(ucl.lw1[i,j])
      ATSt[k]=(deltat[k]*mu.t0)/(1-beta)
    }
    EATSt.lw1[i,j]= mean(ATSt)
  }
}

#weibull-gamma
EATSt.wg1=matrix(0, nrow=3, ncol=3, dimnames = list(c("1","2","5"),c("1","2","5")))

for (i in 1:3){
  for (j in 1:3) {  
    for (k in 1:10){
      F.Z1=function(z){
        func=function(x,mu.t0,mu.x0){
          pweibull((x-z)*mu.t0,shape=a1[i,k],scale=deltat[k]*mu.t0/gamma((1/a1[i,k])+1))*dgamma(x*mu.x0,shape=(mu.x0/s0[j])^2,scale=(s0[j]^2/mu.x0))
        }
        1-mu.x0*integrate(func,0,Inf,mu.t0=mu.t0,mu.x0=mu.x0)$value
      }
      beta=F.Z1(ucl.wg1[i,j])
      ATSt[k]=(deltat[k]*mu.t0)/(1-beta)
    }
    EATSt.wg1[i,j]= mean(ATSt)
  }
}

#weibull-lognormal
EATSt.wl1=matrix(0, nrow=3, ncol=3, dimnames = list(c("1","2","5"),c("1","2","5")))

for (i in 1:3){
  for (j in 1:3) {  
    for (k in 1:10){
      F.Z1=function(z){
        func=function(x,mu.t0,mu.x0){
          pweibull((x-z)*mu.t0,shape=a1[i,k],scale=deltat[k]*mu.t0/gamma((1/a1[i,k])+1))*dlnorm(x*mu.x0,meanlog=log(mu.x0^2/sqrt(s0[j]^2+mu.x0^2)),sdlog=sqrt(log(1+(s0[j]^2/mu.x0^2))))
        }
        1-mu.x0*integrate(func,0,Inf,mu.t0=mu.t0,mu.x0=mu.x0)$value
      }
      beta=F.Z1(ucl.wl1[i,j])
      ATSt[k]=(deltat[k]*mu.t0)/(1-beta)
    }
    EATSt.wl1[i,j]= mean(ATSt)
  }
}

#weibull-normal
EATSt.wn1=matrix(0, nrow=3, ncol=2, dimnames = list(c("1","2","5"),c("1","2")))

for (i in 1:3){
  for (j in 1:2) {  
    for (k in 1:10){
      F.Z1=function(z){
        func=function(x,mu.t0,mu.x0){
          pweibull((x-z)*mu.t0,shape=a1[i,k],scale=deltat[k]*mu.t0/gamma((1/a1[i,k])+1))*dnorm(x*mu.x0,mean=mu.x0,sd=s0[j])
        }
        1-mu.x0*integrate(func,0,Inf,mu.t0=mu.t0,mu.x0=mu.x0)$value
      }
      beta=F.Z1(ucl.wn1[i,j])
      ATSt[k]=(deltat[k]*mu.t0)/(1-beta)
    }
    EATSt.wn1[i,j]= mean(ATSt)
  }
}

#weibull-weibull
EATSt.ww1=matrix(0, nrow=3, ncol=3, dimnames = list(c("1","2","5"),c("1","2","5")))

for(i in 1:3){
  for(j in 1:3){
    for (k in 1:10){
      F.Z1=function(z){
        func=function(x,mu.t0,mu.x0){
          pweibull((x-z)*mu.t0,shape=a1[i,k],scale=deltat[k]*mu.t0/gamma((1/a1[i,k])+1))*dweibull(x*mu.x0,shape=a0[j],scale=mu.x0/gamma((1/a0[j])+1))
        }
        1-mu.x0*integrate(func,0,Inf,mu.t0=mu.t0,mu.x0=mu.x0)$value
      }
      beta=F.Z1(ucl.ww1[i,j])
      ATSt[k]=(deltat[k]*mu.t0)/(1-beta)
    }
    EATSt.ww1[i,j]= mean(ATSt)
  }
}
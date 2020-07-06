



######################################################################################################################################################################
##### TSNP.OC() generates the OC of TSNP design through simulation studies ###########################################################################################
##### TSNP.OC() will report the utility values, OBD selection percentages, patients allocation; overall toxicity rate and efficacy rate ############################## 
##### For the OBD selection percentages, the selection of 0 means the trial is early terminated and no OBD is selected ###############################################
######################################################################################################################################################################



TSNP.OC=function(pi,cohortsize,ncohort,tau,m,phi.t,c.t,phi.e,c.e,omega,start,ntrial){
##### pi: the joint toxicity-efficacy probability matrix; row=dose level; column=(tox=0,eff=0; tox=0,eff=1; tox=1,eff=0; tox=1,eff=1) 
##### cohortsize: the size (number of patients) for each cohort 
##### ncohort: the number of cohort for a phase I/II trial 
##### tau: the margin of meaningful difference of toxicity rate used in the CCD design
##### m: if m patients have been treated at the current dose, switch the trial to the second stage 
##### phi.t: the maximum tolerable toxicity rate 
##### c.t: nominal level of toxicity for admissible set 
##### phi.e: the minimum tolerable efficacy rate 
##### c.e: nominal level of futility for admissible set 
##### omega: the vector of desirability scores for T=0 and E=0; T=0 and E=1; T=1 and E=0; T=1 and E=1 
##### start: the starting dose level 
##### ntrial: the number of simulated trials  
  
  
  library(Iso)
  
  adm.tox=function(n,res,phi.t,c.t){
    p=NULL
    L=length(n)
    for (i in 1:L){
      p[i]=binom.test(res[i],n[i],phi.t,alternative="greater")[[3]]
    }
    if( any(p<c.t)==1  ){
      re=min( which(p<c.t)-1  )
    } else {re=L}
    return(re)
  }
  
  adm.eff=function(n,res,phi.e,c.e){
    p=NULL
    L=length(n)
    for (i in 1:L){
      p[i]=binom.test(res[i],n[i],phi.e, alternative="less")[[3]]
    }
    re=which(p>=c.e)
    if( length(re)==0 ){
      return(0)
    } else {    return(re)   }
  }
  
  
  
  stageI=function(pi,cohortsize,ncohort,tau,m,phi.t,c.t,start){
    ndose=dim(pi)[2]
    y00=rep(0,ndose)
    y01=rep(0,ndose)
    y10=rep(0,ndose)
    y11=rep(0,ndose)
    yt=rep(0,ndose)
    n=rep(0,ndose)
    d=start
    for (i in 1: ncohort){
      if (n[d]>=m){break}
      n[d]=n[d]+cohortsize
      res=rmultinom(1,cohortsize,pi[,d])
      y00[d]=y00[d]+res[1]
      y01[d]=y01[d]+res[2]
      y10[d]=y10[d]+res[3]
      y11[d]=y11[d]+res[4]
      yt[d]=y10[d]+y11[d]
      try=which(n>0) ## all the tired dose
      try.l=min(try) ## low bound for try
      try.u=max(try) ## up bound for try
      at=adm.tox(n[try],yt[try],phi.t,c.t)
      adm.high=at+try.l-1
      if (adm.high==0){
        break
      }else if (adm.high<try.l){
        d=adm.high
      } else{
        ptox=(yt[try])/n[try]
        piso=pava(ptox,w=n[try])
        if (adm.high==try.u){
          if (piso[d-try.l+1]<(phi.t-tau)){
            d=min((d+1), ndose   )
          } else if (piso[d-try.l+1]>(phi.t+tau)){
            d=max(  (d-1),1 )
          } else{d=d}
        }else {
          if (piso[d-try.l+1]<(phi.t-tau)){
            d=min((d+1), adm.high  )
          } else if (piso[d-try.l+1]>(phi.t+tau)){
            d=min(max(  (d-1),1 ), adm.high)
          } else {d=min(d,adm.high)}
        }
      }
      
      
    }
    
    return( list( "y00"=y00,"y01"=y01,"y10"=y10,"y11"=y11   )  )
  }
  
  
  
  
  findobd=function(y00,y01,y10,y11,omega,admset){
    n=y00+y01+y10+y11
    q=(y10+y11)/n
    q=pava(q,w=n)
    p00=NULL
    p01=NULL
    p10=NULL
    p11=NULL
    u=NULL
    for(i in 1: length(n)){
      if( (y00[i]+y01[i])==0   ){
        p00[i]=(1-q[i])/2
        p01[i]=(1-q[i])/2
      }else{
        p00[i]=y00[i]*(1-q[i])/(y00[i]+y01[i])
        p01[i]=y01[i]*(1-q[i])/(y00[i]+y01[i])              
      } 
      if ((y10[i]+y11[i])==0){
        p10[i]=q[i]/2
        p11[i]=q[i]/2
      }else{
        p10[i]=q[i]*y10[i]/(y10[i]+y11[i])
        p11[i]=q[i]*y11[i]/(y10[i]+y11[i])
      }
      u[i]=omega[1]*p00[i]+omega[2]*p01[i]+omega[3]*p10[i]+omega[4]*p11[i]
    }
    obd=which(u[admset]==max(u[admset]))[1]
    return(obd)
  }
  
  dosefind=function(y00,y01,y10,y11,omega,phi.t,c.t,phi.e,c.e){
    yt=y10+y11
    ye=y01+y11
    n=y00+y01+y10+y11
    ndose=length(n)
    try=which(n>0) ## all the tired dose
    try.l=min(try) ## low bound for try
    try.u=max(try)
    at=adm.tox(n[try],yt[try],phi.t,c.t)
    adm.high=at+try.l-1
    if (adm.high<try.l){
      dselect=adm.high
    } else if ((adm.high==try.u)&(try.u<ndose)) {
      dselect=adm.high+1
    } else {
      ae=adm.eff(n[try.l: adm.high], ye[try.l: adm.high], phi.e, c.e)
      if (sum(ae)==0){
        dselect=try.l-1
      } else {
        temp=findobd(y00[try],y01[try],y10[try],y11[try],omega,ae)
        dselect=ae[temp]+try.l-1
      }
    }
    return(dselect)
    
  }
  
  stageII=function(y00,y01,y10,y11,pi,cohortsize,ncohort,omega,phi.t,c.t,phi.e,c.e){
    n=y00+y01+y10+y11
    if(sum(n)==cohortsize*ncohort){
      yt=y10+y11
      ye=y01+y11
      ndose=length(n)
      try=which(n>0) ## all the tired dose
      try.l=min(try) ## low bound for try
      at=adm.tox(n[try],yt[try],phi.t,c.t)
      adm.high=at+try.l-1
      if (adm.high<try.l){
        dselect=0
      } else {
        ae=adm.eff(n[try.l: adm.high], ye[try.l: adm.high], phi.e, c.e)
        if (sum(ae)==0){
          dselect=0
        } else {
          temp=findobd(y00[try],y01[try],y10[try],y11[try],omega,ae)
          dselect=ae[temp]+try.l-1
        }
      }
      return(  list( "OBD"=dselect,"n"=n,"toxrate"= (sum(y10)+sum(y11)    )/sum(n),"effrate"=( sum(y01)+sum(y11)  )/sum(n)   ) )
    } else{
      estop=0
      ncohort2=ncohort-sum(n)/cohortsize
      for (i in 1: ncohort2){
        d=dosefind(y00,y01,y10,y11,omega,phi.t,c.t,phi.e,c.e)
        if(d==0){
          estop=1
          break
        } else {
          n[d]=n[d]+cohortsize
          res=rmultinom(1,cohortsize,pi[,d])
          y00[d]=y00[d]+res[1]
          y01[d]=y01[d]+res[2]
          y10[d]=y10[d]+res[3]
          y11[d]=y11[d]+res[4]
        }
        
      }
      if (estop==1){
        dselect=0
      } else {
        yt=y10+y11
        ye=y01+y11
        ndose=length(n)
        try=which(n>0) ## all the tired dose
        try.l=min(try) ## low bound for try
        at=adm.tox(n[try],yt[try],phi.t,c.t)
        adm.high=at+try.l-1
        if (adm.high<try.l){
          dselect=0
        } else {
          ae=adm.eff(n[try.l: adm.high], ye[try.l: adm.high], phi.e, c.e)
          if (sum(ae)==0){
            dselect=0
          } else {
            temp=findobd(y00[try],y01[try],y10[try],y11[try],omega,ae)
            dselect=ae[temp]+try.l-1
          }
        }
        
        
        
      }
      return(  list( "OBD"=dselect,"n"=n,"toxrate"= (sum(y10)+sum(y11)    )/sum(n),"effrate"=( sum(y01)+sum(y11)  )/sum(n)   ) )    
    }
  }
  
  ndose=dim(pi)[2]
  tu=NULL
  for(i in 1: ndose){
    tu[i]=sum(pi[,i]*omega)
  }
  N=matrix( rep(0,ndose*ntrial),ncol=ndose    )
  OBD=rep(0,ntrial)
  Toxrate=rep(0,ntrial)
  Effrate=rep(0,ntrial)
  for(i in 1:ntrial){
    data1=stageI(pi,cohortsize,ncohort,tau,m,phi.t,c.t,start)
    y00=data1[[1]]
    y01=data1[[2]]
    y10=data1[[3]]
    y11=data1[[4]]
    data2=stageII(y00,y01,y10,y11,pi,cohortsize,ncohort,omega,phi.t,c.t,phi.e,c.e)
    OBD[i]=data2[[1]]
    N[i,]=data2[[2]]
    Toxrate[i]=data2[[3]]
    Effrate[i]=data2[[4]]
  }
  return(list("True Utility Value"=round(tu,digits=1), "OBD Selection Percentages"=round(table(OBD)*100/ntrial,digits=1) , "Average Number of Patients"=round(apply(N,2,mean),digits=1), "Average Efficacy Rate"=round(mean(Effrate)*100,digits=1), "Average Toxicity Rate"=round(mean(Toxrate)*100,digits=1)   ))   
  
}

##### Example 1: To generate the simulation result of scenario 1 in Table 1

crossratio=function(ttox,teff,gamma){
  ndose=length(ttox) ## dose level
  out=matrix(rep(0,4*ndose),nrow=4) ## joint tox-eff probability matrix; row is dose level, column=(tox=0,eff=0; tox=0,eff=1; tox=1,eff=0; tox=1,eff=1)
  for(j in 1:ndose){
    if (gamma[j]==1){
      out[1,j]=(1-ttox[j])*(1-teff[j])
      out[2,j]=(1-ttox[j])*(teff[j])
      out[3,j]=(ttox[j])*(1-teff[j])
      out[4,j]=(ttox[j])*(teff[j])
    }else{
      a=1+(ttox[j]+teff[j])*(gamma[j]-1)
      b=-4*gamma[j]*(gamma[j]-1)*ttox[j]*teff[j]
      out[4,j]=(a-sqrt(a^2+b))/(2*(gamma[j]-1))
      out[3,j]=ttox[j]-out[4,j]
      out[2,j]=teff[j]-out[4,j]
      out[1,j]=1-ttox[j]-teff[j]+out[4,j]
    }
  }
  return(out)
}

pi=crossratio(ttox=c(0.05,0.1,0.15,0.25,0.4), teff=c(0.15,0.2,0.4,0.3,0.2), gamma=c(1,2,2,3,3))

set.seed(1)

TSNP.OC(pi=pi,cohortsize=3,ncohort=20,tau=0.1,m=12,phi.t=0.3, c.t=0.1, phi.e=0.4, c.e=0, omega=c(25,100,0,50), start=1,ntrial=10000)



##### Example 2: To gerate the simulation result of scenario 1 in Table 3


Gumbel=function(ttox,teff,c){
  ## ttox: marginal toxicity probability
  ## teff: marginal efficacy probability
  ## c: association parameter between tox and eff, c>0 indicates a positive correlation
  ndose=length(ttox) ## dose level
  out=matrix(rep(0,4*ndose),nrow=4) ## joint tox-eff probability matrix; row is dose level, column=(tox=0,eff=0; tox=0,eff=1; tox=1,eff=0; tox=1,eff=1)
  for (j in 1: ndose){
    out[1,j]=(1-ttox[j])*(1-teff[j])+ttox[j]*teff[j]*(1-ttox[j])*(1-teff[j])*( exp(c)-1   )/(exp(c)+1)
    out[2,j]=(1-ttox[j])*(teff[j])-ttox[j]*teff[j]*(1-ttox[j])*(1-teff[j])*( exp(c)-1   )/(exp(c)+1)
    out[3,j]=(ttox[j])*(1-teff[j])-ttox[j]*teff[j]*(1-ttox[j])*(1-teff[j])*( exp(c)-1   )/(exp(c)+1)
    out[4,j]=(ttox[j])*(teff[j])+ttox[j]*teff[j]*(1-ttox[j])*(1-teff[j])*( exp(c)-1   )/(exp(c)+1)
  }
  return(out)
}


pi=Gumbel(ttox=c(0.05,0.1,0.15,0.4,0.5), teff=c(0.1,0.3,0.5,0.6,0.7), c=0.2)

set.seed(1)

TSNP.OC(pi=pi,cohortsize=3,ncohort=20,tau=0.1,m=12,phi.t=0.3, c.t=0.1, phi.e=0.2, c.e=0.2, omega=c(25,100,0,50), start=1,ntrial=10000)






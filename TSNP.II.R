



#########################################################################################################
##### TSNP.II() implements stage II of the TSNP design #####################################################
#########################################################################################################


TSNP.II=function(y00,y01,y10,y11,omega,phi.t,c.t,phi.e,c.e,ss){ 
##### y00: number of patients with tox=0, eff=0 at different doses
##### y01: number of patients with tox=0, eff=1 at different doses
##### y10: number of patients with tox=1, eff=0 at different doses
##### y11: number of patients with tox=1, eff=1 at different doses
##### omega: the vector of desirability scores for T=0 and E=0; T=0 and E=1; T=1 and E=0; T=1 and E=1  
##### phi.t: the maximum tolerable toxicity rate   
##### c.t: nominal level of toxicity for admissible set 
##### phi.e: the minimum tolerable efficacy rate  
##### c.e: nominal level of futility for admissible set  
##### ss: maximum sam[le size  
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
  
  ute=function(y00,y01,y10,y11,omega){
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
    return(u)
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
  
  
  n=y00+y01+y10+y11
  ndose=length(n)
  if(sum(n)>=ss){
    yt=y10+y11
    ye=y01+y11
    ndose=length(n)
    try=which(n>0) ## all the tired dose
    pt=NULL
    pe=NULL
    L=length(try)
    for (i in 1:L){
      pt[i]=binom.test(yt[try[i]],n[try[i]],phi.t,alternative="greater")[[3]]
      pe[i]=binom.test(ye[try[i]],n[try[i]],phi.e, alternative="less")[[3]]
    }
    
    
    try.l=min(try) ## low bound for try
    try.u=max(try)
    pt2=c( rep(NA, try.l-1), pt, rep(NA, ndose-try.u)    )
    pe2=c( rep(NA, try.l-1), pe, rep(NA, ndose-try.u)    )
    at=adm.tox(n[try],yt[try],phi.t,c.t)
    adm.high=at+try.l-1
    if (adm.high<try.l){
      print("P-values testing for overly toxic ") 
      cat(formatC(pt2, digits=2, format="f"), sep="  ", "\n")
      print("P-values testing for less efficacious ") 
      cat(formatC(pe2, digits=2, format="f"), sep="  ", "\n")
      print("No dose is selected.")
      
    } else {
      ae=adm.eff(n[try.l: adm.high], ye[try.l: adm.high], phi.e, c.e)
      if (sum(ae)==0){
        print("P-values testing for overly toxic ") 
        cat(formatC(pt2, digits=2, format="f"), sep="  ", "\n")
        print("P-values testing for less efficacious ") 
        cat(formatC(pe2, digits=2, format="f"), sep="  ", "\n")
        print("No dose is selected.")
        
      } else {
        temp=findobd(y00[try],y01[try],y10[try],y11[try],omega,ae)
        dselect=ae[temp]+try.l-1
        utility=ute(y00[try],y01[try],y10[try],y11[try],omega)
        utility2=c( rep(NA, try.l-1), utility, rep(NA, ndose-try.u)    )
        
        print("P-values testing for overly toxic ") 
        cat(formatC(pt2, digits=2, format="f"), sep="  ", "\n")
        print("P-values testing for less efficacious ") 
        cat(formatC(pe2, digits=2, format="f"), sep="  ", "\n")
        print("Utility value estimates ")
        cat(formatC(utility2, digits=2, format="f"), sep="  ", "\n")
        print("The final recommended OBD")
        cat(formatC(dselect, digits=0, format="f"), sep="  ", "\n")
      }
    }
  } else{
    yt=y10+y11
    ye=y01+y11
    
    try=which(n>0) ## all the tired dose
    pt=NULL
    pe=NULL
    try.l=min(try) ## low bound for try
    try.u=max(try)
   
    
    
    L=length(try)
    for (i in 1:L){
      pt[i]=binom.test(yt[try[i]],n[try[i]],phi.t,alternative="greater")[[3]]
      pe[i]=binom.test(ye[try[i]],n[try[i]],phi.e, alternative="less")[[3]]
    }
    pt2=c( rep(NA, try.l-1), pt, rep(NA, ndose-try.u)    )
    pe2=c( rep(NA, try.l-1), pe, rep(NA, ndose-try.u)    )
    
    utility=ute(y00[try],y01[try],y10[try],y11[try],omega)
    utility2=c( rep(NA, try.l-1), utility, rep(NA, ndose-try.u)    )
    d=dosefind(y00,y01,y10,y11,omega,phi.t,c.t,phi.e,c.e)
    if(d==0){
      print("P-values testing for overly toxic ") 
      cat(formatC(pt2, digits=2, format="f"), sep="  ", "\n")
      print("P-values testing for less efficacious ") 
      cat(formatC(pe2, digits=2, format="f"), sep="  ", "\n")
      print("No dose is selected.")
    }else{
      print("P-values testing for overly toxic ") 
      cat(formatC(pt2, digits=2, format="f"), sep="  ", "\n")
      print("P-values testing for less efficacious ") 
      cat(formatC(pe2, digits=2, format="f"), sep="  ", "\n")
      print("Utility value estimates ")
      cat(formatC(utility2, digits=2, format="f"), sep="  ", "\n")
      print("The dose level for next cohort of patients")
      cat(formatC(d, digits=0, format="f"), sep="  ", "\n")
    }
  }
}


##### Example: Practical trial implementation as shown in Figure 4

## cohort 6
data=matrix(  c(3,0,0,0,3,0,0,0,2,1,0,0,2,2,0,2,1,0,2,0)  ,nrow=4  )
TSNP.II(y00=data[1,],y01=data[2,],y10=data[3,],y11=data[4,],omega=c(25,100,0,50),phi.t=0.3,c.t=0.1,phi.e=0.2,c.e=0,ss=30)


## cohort 7
data=matrix(  c(3,0,0,0,3,0,0,0,2,1,0,0,3,2,1,3,1,0,2,0)  ,nrow=4  )
TSNP.II(y00=data[1,],y01=data[2,],y10=data[3,],y11=data[4,],omega=c(25,100,0,50),phi.t=0.3,c.t=0.1,phi.e=0.2,c.e=0,ss=30)


## cohort 8
data=matrix(  c(3,0,0,0,3,0,0,0,3,1,0,2,3,2,1,3,1,0,2,0)  ,nrow=4  )
TSNP.II(y00=data[1,],y01=data[2,],y10=data[3,],y11=data[4,],omega=c(25,100,0,50),phi.t=0.3,c.t=0.1,phi.e=0.2,c.e=0,ss=30)




## cohort 9
data=matrix(  c(3,0,0,0,3,0,0,0,3,1,0,2,5,3,1,3,1,0,2,0)  ,nrow=4  )
TSNP.II(y00=data[1,],y01=data[2,],y10=data[3,],y11=data[4,],omega=c(25,100,0,50),phi.t=0.3,c.t=0.1,phi.e=0.2,c.e=0,ss=30)


## cohort 10
data=matrix(  c(3,0,0,0,3,0,0,0,3,1,0,2,6,4,2,3,1,0,2,0)  ,nrow=4  )
TSNP.II(y00=data[1,],y01=data[2,],y10=data[3,],y11=data[4,],omega=c(25,100,0,50),phi.t=0.3,c.t=0.1,phi.e=0.2,c.e=0,ss=30)















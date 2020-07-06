





#########################################################################################################
##### TSNP.I() implements stage I of the TSNP design #####################################################
#########################################################################################################


TSNP.I=function(y00,y01,y10,y11,tau,m,phi.t,c.t,current){
##### y00: number of patients with tox=0, eff=0 at different doses
##### y01: number of patients with tox=0, eff=1 at different doses
##### y10: number of patients with tox=1, eff=0 at different doses
##### y11: number of patients with tox=1, eff=1 at different doses
##### tau: margin of meaningful difference of toxicity rate used in the CCD design
##### m: if m patients have been treated at the current dose, switch the trial to the second stage
##### phi.t: the maximum tolerable toxicity rate 
##### c.t: nominal level of toxicity for admissible set   
##### current: current dose level  
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
    return(list(p, re))
  }
  n=y00+y01+y10+y11
  yt=y10+y11
  ndose=length(n)
  try=which(n>0) ## all the tired dose
  try.l=min(try) ## low bound for try
  try.u=max(try) ## up bound for try
  temp=adm.tox(n[try],yt[try],phi.t,c.t)
  pv=temp[[1]]
  pv2=c( rep(NA, try.l-1), pv, rep(NA, ndose-try.u)    )
  at=temp[[2]]
  adm.high=at+try.l-1
  if (adm.high==0){
    print("P-values testing for overly toxic ") 
    cat(formatC(pv2, digits=2, format="f"), sep="  ", "\n")
    print("The trial is early terminated.")
  } else if (adm.high<try.l){
    d=adm.high
    print("P-values testing for overly toxic") 
    cat(formatC(pv2, digits=2, format="f"), sep="  ", "\n")
    if (n[current]>=m){
      print ("This is the end of stage I. Please use TSNP.II() for second stage of dose-finding.")
    } else{
      print("The dose level for next cohort of patients")
      cat(formatC(d, digits=0, format="f"), sep="  ", "\n")
    }    
  } else{
    ptox=(yt[try])/n[try]
    piso=pava(ptox,w=n[try])
    if (adm.high==try.u){
      if (piso[current-try.l+1]<(phi.t-tau)){
        d=min((current+1), ndose   )
      } else if (piso[current-try.l+1]>(phi.t+tau)){
        d=max(  (current-1),1 )
      } else{d=current}
    }else {
      if (piso[current-try.l+1]<(phi.t-tau)){
        d=min((current+1), adm.high  )
      } else if (piso[current-try.l+1]>(phi.t+tau)){
        d=min(max(  (current-1),1 ), adm.high)
      } else {d=min(current,adm.high)}
    }
    
    print("P-values testing for overly toxic ") 
    cat(formatC(pv2, digits=2, format="f"), sep="  ", "\n")
    if (n[current]>=m){
      print ("This is the end of stage I. Please use TSNP.II() for second stage of dose-finding.")
    } else{
      print ("The  isotonic estimate of the toxicity rate of current dose level")
      cat(formatC(piso[current-try.l+1], digits=2, format="f"), sep="  ", "\n")
      print("The dose level for next cohort of patients")
      cat(formatC(d, digits=0, format="f"), sep="  ", "\n")
    }   
  }
  
  
  
}

##### Example: Practical trial implementation as shown in Figure 4

## cohort 1
data=matrix(  c(3,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0)  ,nrow=4  )
TSNP.I(y00=data[1,],y01=data[2,],y10=data[3,],y11=data[4,],tau=0.1,m=6,phi.t=0.3,c.t=0.1,current=1)
  
## cohort 2  
data=matrix(  c(3,0,0,0,3,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0)  ,nrow=4  )
TSNP.I(y00=data[1,],y01=data[2,],y10=data[3,],y11=data[4,],tau=0.1,m=6,phi.t=0.3,c.t=0.1,current=2)

## cohort 3
data=matrix(  c(3,0,0,0,3,0,0,0,2,1,0,0,0,0,0,0,0,0,0,0)  ,nrow=4  )
TSNP.I(y00=data[1,],y01=data[2,],y10=data[3,],y11=data[4,],tau=0.1,m=6,phi.t=0.3,c.t=0.1,current=3)


## cohort 4
data=matrix(  c(3,0,0,0,3,0,0,0,2,1,0,0,1,2,0,0,0,0,0,0)  ,nrow=4  )
TSNP.I(y00=data[1,],y01=data[2,],y10=data[3,],y11=data[4,],tau=0.1,m=6,phi.t=0.3,c.t=0.1,current=4)


## cohort 5
data=matrix(  c(3,0,0,0,3,0,0,0,2,1,0,0,1,2,0,0,1,0,2,0)  ,nrow=4  )
TSNP.I(y00=data[1,],y01=data[2,],y10=data[3,],y11=data[4,],tau=0.1,m=6,phi.t=0.3,c.t=0.1,current=5)










  
  
  
  

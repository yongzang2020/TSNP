 

##############################################################################################################################################################
##### crossratio() generates the joint toxicity-efficacy probabilities based on the cross ratio model ########################################################
##### The output is a joint toxicity-efficacy probability matrix #############################################################################################
#### The rows of the matrxs represent dose level #############################################################################################################
#### The columns of the matrix represnt the probability for different outcomes; column=(tox=0,eff=0; tox=0,eff=1; tox=1,eff=0; tox=1,eff=1) ##################
##############################################################################################################################################################




crossratio=function(ttox,teff,gamma){
  ## ttox: marginal toxicity probability
  ## teff: marginal efficacy probability
  ## gamma: cross ratio between tox and eff, gamma>1 indicates a positive correlation
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


##### Example: to generate the joint probabilities for scenario 1 in Table 1.

pi=crossratio(ttox=c(0.05,0.1,0.15,0.25,0.4), teff=c(0.15,0.2,0.4,0.3,0.2), gamma=c(1,2,2,3,3))




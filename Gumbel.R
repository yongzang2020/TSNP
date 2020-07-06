


##############################################################################################################################################################
##### Gumbel() generates the joint toxicity-efficacy probabilities based on the Gumbel copula model ########################################################
##### The output is a joint toxicity-efficacy probability matrix #############################################################################################
#### The rows of the matrix represent dose level #############################################################################################################
#### The columns of the matrix represnt the probability for different outcomes; column=(tox=0,eff=0; tox=0,eff=1; tox=1,eff=0; tox=1,eff=1) ##################
##############################################################################################################################################################



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


##### Example: to generate the joint probabilities for scenario 1 in Table 3.

pi=Gumbel(ttox=c(0.05,0.1,0.15,0.4,0.5), teff=c(0.1,0.3,0.5,0.6,0.7), c=0.2)

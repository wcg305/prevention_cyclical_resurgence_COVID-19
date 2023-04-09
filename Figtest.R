
#####################################################
#
#    Code for paper "Prevention of cyclical resurgences of COVID-19-like pandemics in the long term:What are the trade-offs?"
#    Author: Ichiro Nakamoto
     Date: 04/09/2023
#####################################################
library(deSolve)
library(ggplot2)
library(reshape2)
library(pracma)
library(grid)
library(gridExtra)
library(plotly)
library(gplots)
library(foreach)
library(doParallel)
library(abind)
library(RColorBrewer)
library(ggpattern)
library(plotly)
library(tikzDevice)
library(cowplot)
library(patchwork)
library(ggpubr)
library(grid)
library(gridtext)
library(svglite)
library(latex2exp)
library(RConics)
library(metR)
library(cwhmisc)




##########################################################################################
######################################ODE MODELS##########################################
##########################################################################################

############################################################
##Seasonality and Social Distancing and Vaccination
############################################################

ode_function<-function (t,x,parameters) {

  SP=x[1]
  IP=x[2]
  R=x[3]
  SS=x[4]
  IS=x[5]
  V1=x[6]
  V2=x[7]
  V3=x[8]
  IV=x[9]
  #V=x[7]
  SS1=x[10]
  SS2=x[11]
  SS3=x[12]
  IS1=x[13]
  IS2=x[14]
  IS3=x[15]

  N=parameters[["N"]]
  mu=parameters[["mu"]]
  gamma=parameters[["gamma"]]
  alpha=parameters[["alpha"]]
  R0.list=parameters[["R0.list"]]
  R0red=parameters[["R0red"]]
  #############################

  alphaV=parameters[["alphaV"]]
  tred1=parameters[["tred1"]]
  tred2=parameters[["tred2"]]
  tred3=parameters[["tred3"]]
  tred4=parameters[["tred4"]]
  delta=parameters[["delta"]]
  rho2=parameters[["rho2"]]
  tvax=parameters[["tvax"]]
  tvax2=parameters[["tvax2"]]
  tvax3=parameters[["tvax3"]]
#######################################

  nu=parameters[["nu"]]
  primary.burden=parameters[["primary.burden"]]
  secondary.burden=parameters[["secondary.burden"]]
  epsilon=parameters[["epsilon"]]

###########################################

  epsilonV1=parameters[["epsilonV1"]]
  epsilonV2=parameters[["epsilonV2"]]
  epsilon1=parameters[["epsilon1"]]
  epsilon2=parameters[["epsilon2"]]
  omega=parameters[["omega"]]
  rho1=parameters[["rho1"]]

###########################################

  one.dose.burden=parameters[["one.dose.burden"]]
  two.doses.burden=parameters[["two.doses.burden"]]
  three.doses.burden=parameters[["three.doses.burden"]]
  vaccine.burden=parameters[["vaccine.burden"]]
  d1=parameters[["d1"]]
  d2=parameters[["d2"]]
  indicator=parameters[["indicator"]]

  ###########################################

  start.value=parameters[["start.value"]]
  max.y.cases=parameters[["max.y.cases"]]
  rho3=parameters[["rho3"]]
  epsilon3=parameters[["epsilon3"]]
  epsilonV3=parameters[["epsilonV3"]]
  omega1=parameters[["omega1"]]
  ###########################################


  alpha1=parameters[["alpha1"]]
  alpha2=parameters[["alpha2"]]
  alpha3=parameters[["alpha3"]]

  R0red.increase=parameters[["R0red.increase"]]


  beta=R0.list[t]*gamma



  #############################
  #epsilonV3=epsilonV3[t]
  #epsilon3=epsilon3[t]

  if((t > tred1) & (t < tred2)){
    beta = beta * R0red
  }
  if((t >= tred2) & (t < tred3)){
    beta = beta * R0red
  }
  if((t >= tred3) & (t < tred4)){
    beta = beta * R0red * R0red.increase

  }
  svax=0
  if(t>=tvax){
    svax=1
  }
  svax2=0
  if(t>=tvax2){
    svax2=1
  }
  svax3=0
  if(t>=tvax3){
    svax3=1
  }

  #ODE Equations test
  #dSPdt=mu-(beta*SP*(IP+alpha*IS+alphaV*IV))-mu*SP-svax*nu*SP
  #dIPdt=(beta*SP*(IP+alpha*IS+alphaV*IV))-(gamma+mu)*IP
  #dRdt=gamma*(IP+IS+IV)-(delta+mu)*R
  #dSSdt=delta*R-epsilon*(beta*SS*(IP+alpha*IS+alphaV*IV))-mu*SS-svax*nu*SS
  #dISdt=epsilon*(beta*SS*(IP+alpha*IS+alphaV*IV))-(gamma+mu)*IS
  #dIVdt=(epsilonV*beta*V)*(IP+alpha*IS+alphaV*IV)-(gamma+mu)*IV
  #dVdt=svax*nu*(SP+d1*SS)-epsilonV*beta*V*(IP+alpha*IS+alphaV*IV)-(rho+mu)*V


  #ODE Equations
  dSPdt=mu-(beta*SP*(IP+alpha*IS+alphaV*IV+alpha1*IS1+alpha2*IS2+alpha3*IS3))-mu*SP-svax*nu*SP
  dIPdt=(beta*SP*(IP+alpha*IS+alphaV*IV+alpha1*IS1+alpha2*IS2+alpha3*IS3))-(gamma+mu)*IP
  dRdt=gamma*(IP+IS+IV+IS1+IS2+IS3)-(delta+mu)*R
  dSSdt=delta*R-epsilon*(beta*SS*(IP+alpha*IS+alphaV*IV+alpha1*IS1+alpha2*IS2+alpha3*IS3))-mu*SS-svax*nu*SS
  dISdt=epsilon*(beta*SS*(IP+alpha*IS+alphaV*IV+alpha1*IS1+alpha2*IS2+alpha3*IS3))-(gamma+mu)*IS
  dV1dt=svax*nu*(SP+d1*SS)-epsilonV1*beta*V1*(IP+alpha*IS+alphaV*IV+alpha1*IS1+alpha2*IS2+alpha3*IS3)-(omega*svax2+rho1+mu)*V1
  dV2dt=d2*svax*nu*SS+omega*svax2*V1-epsilonV2*beta*V2*(IP+alpha*IS+alphaV*IV+alpha1*IS1+alpha2*IS2+alpha3*IS3)-(omega1*svax3+rho2+mu)*V2
  dV3dt=(1-d1-d2)*svax*nu*SS+omega1*svax3*V2-epsilonV3*beta*V3*(IP+alpha*IS+alphaV*IV+alpha1*IS1+alpha2*IS2+alpha3*IS3)-(rho3+mu)*V3
  dIVdt=(epsilonV1*beta*V1+epsilonV2*beta*V2+epsilonV3*beta*V3)*(IP+alpha*IS+alphaV*IV+alpha1*IS1+alpha2*IS2+alpha3*IS3)-(gamma+mu)*IV
  dSS1dt=rho1*V1-epsilon1*beta*SS1*(IP+alpha*IS+alphaV*IV+alpha1*IS1+alpha2*IS2+alpha3*IS3)-mu*SS1
  dSS2dt=rho2*V2-epsilon2*beta*SS2*(IP+alpha*IS+alphaV*IV+alpha1*IS1+alpha2*IS2+alpha3*IS3)-mu*SS2
  dSS3dt=rho3*V3-epsilon3*beta*SS3*(IP+alpha*IS+alphaV*IV+alpha1*IS1+alpha2*IS2+alpha3*IS3)-mu*SS3
  dIS1dt=epsilon1*beta*SS1*(IP+alpha*IS+alphaV*IV+alpha1*IS1+alpha2*IS2+alpha3*IS3)-(gamma+mu)*IS1
  dIS2dt=epsilon2*beta*SS2*(IP+alpha*IS+alphaV*IV+alpha1*IS1+alpha2*IS2+alpha3*IS3)-(gamma+mu)*IS2
  dIS3dt=epsilon3*beta*SS3*(IP+alpha*IS+alphaV*IV+alpha1*IS1+alpha2*IS2+alpha3*IS3)-(gamma+mu)*IS3

   dxdt <- c(dSPdt,dIPdt,dRdt,dSSdt,dISdt,dV1dt,dV2dt,dV3dt,dIVdt,dSS1dt,dSS2dt,dSS3dt,dIS1dt,dIS2dt,dIS3dt)

   #dxdt <- c(dSPdt,dIPdt,dRdt,dSSdt,dISdt,dIVdt,dVdt)


  list(dxdt)

}


##########################################################################################
######################################Data input##########################################
##########################################################################################
NYdata=read.table("data.csv",header=T,sep=",")
hku1=NYdata$hku1

##########################################################################################
######################################PLOTS###############################################
##########################################################################################

######################################
###########Run the ODEs###############
######################################

noweeks = 52
gridsize = 50

N=1
mu=1/(20*52)
alpha=alpha1=alpha2=alpha3=1

gamma=7/5
R0red=0.45
R0red.increase=1.325
#Social distancing
##For Soc Dis on want to use R0red=0.6
#R0red=0.6
#tred1=15
#tred2=56
#tred3=0
#tred4=0
#Seasonality
times<-seq(from=1,to=4*noweeks,by=1)

R0.listOrig=read.table("NY.csv",header=T,sep=",")
R0.list.a=2.3*R0.listOrig$hku1/(mean(R0.listOrig$hku1))
R0.list=rep(R0.list.a,length=length(times))
#R0.list <- rep(hku1,length=length(times))
##If you want to shift seasonal cycle do it here
#seasonlag <- 0
#R0.list <- circshift(R0.list,-seasonlag)


#epsilon=0.7
ImmuneTime=3
#delta=1/ImmuneTime

#Run the primary scenario
#epsilon3=0.1
rho3=1/(0.5*52)

omega1 =0

parameters_primary=list(N=N,mu=mu,gamma=gamma,alpha=alpha,R0.list = R0.list, R0red=R0red,R0red.increase=R0red.increase,
                        alphaV=1,tred1=7,tred2=35,tred3=48,tred4=80,delta=1/4,rho2=1/(1.5*52),tvax=44,tvax2=52,tvax3=76,
                        nu=0.01,primary.burden=0.14,secondary.burden=0.07,epsilon=0.5,
                        epsilonV1=0.2,epsilonV2=0.05,epsilon1=0.5,epsilon2=0.3,omega=1/20,rho1=1/(0.125*52),
                        alpha1=alpha1,alpha2=alpha2,alpha3=alpha3,vaccine.burden=0.1,d1=0.01,d2=0.01,indicator="no",
                        start.value=0,max.y.cases=.04,rho3=rho3,epsilon3=0.1,epsilonV3=0.05,omega1=omega1)




times<-seq(from=1,to=noweeks*4,by=1)
I0=1e-9

xstart=c(SP=1-I0,IP=I0,R=0,SS=0,IS=0,V1=0,V2=0,V3=0,IV=0,SS1=0,SS2=0,SS3=0,IS1=0,IS2=0,IS3=0)

#xstart=c(SP=1-I0,IP=I0,R=0,SS=0,IS=0,IV=0,V=0)

out_primary<-as.data.frame(rk(xstart,times, ode_function,parameters_primary))



#epsilon3=seq(from=0,to=1,length.out=gridsize)
#epsilonV3=seq(from=0,to=1,length.out=gridsize)

allIvax=function(omega1,rho3){



  #epsilon3=0.1
  #omega1=4

  delta <- 1/4

  parameters_booster_faster<-list(N=N,mu=mu,gamma=gamma,alpha=alpha,R0.list = R0.list, R0red=R0red,R0red.increase=R0red.increase,
                                      alphaV=1,tred1=7,tred2=35,tred3=48,tred4=80,delta=delta,rho2=1/(1.5*52),tvax=44,tvax2=52,tvax3=76,
                                      nu=0.01,primary.burden=0.14,secondary.burden=0.07,epsilon=0.5,
                                      epsilonV1=0.2,epsilonV2=0.05,epsilon1=0.5,epsilon2=0.3,omega=1/20,rho1=1/(0.125*52),
                                      alpha1=alpha1,alpha2=alpha2,alpha3=alpha3,vaccine.burden=0.1,d1=0.01,d2=0.01,indicator="no",
                                      start.value=0,max.y.cases=.04,rho3=rho3,epsilon3=0.1,epsilonV3=0.05,omega1=omega1)


  times<-seq(from=1,to=noweeks*4,by=1)
  I0=1e-9
  xstart=c(SP=1-I0,IP=I0,R=0,SS=0,IS=0,V1=0,V2=0,V3=0,IV=0,SS1=0,SS2=0,SS3=0,IS1=0,IS2=0,IS3=0)

  out_booster_faster<-as.data.frame(rk(xstart,times, ode_function, parameters_booster_faster))


  #Total
  Itotseries_booster_faster <- out_booster_faster$IP + out_booster_faster$IS1+ out_booster_faster$IS2+ out_booster_faster$IS3+
    out_booster_faster$IS+out_booster_faster$IV

  Itot_booster_faster <- gamma * sum(Itotseries_booster_faster[(44):(4*noweeks)])


  Itotseries_primary <- out_primary$IP + out_primary$IS1+ out_primary$IS2+ out_primary$IS3+
    out_primary$IS+ out_primary$IV

  Itot_primary <- gamma * sum(Itotseries_primary[(44):(4*noweeks)])

  #Primary
  Ipseries_booster_faster <- out_booster_faster$IP
  Ip_booster_faster <- gamma * sum(Ipseries_booster_faster[(44):(4*noweeks)])
  Ipseries_primary <- out_primary$IP
  Ip_primary <- gamma * sum(Ipseries_primary[(44):(4*noweeks)])

  #Secondary
  Isseries_booster_faster <- out_booster_faster$IS
  Is_booster_faster <- gamma * sum(Isseries_booster_faster[(44):(4*noweeks)])
  Isseries_primary <- out_primary$IS
  Is_primary <- gamma * sum(Isseries_primary[(44):(4*noweeks)])


  #Vaccine
  Ivseries_booster_faster <- out_booster_faster$IV+out_booster_faster$IS1+out_booster_faster$IS2+out_booster_faster$IS3
  Iv_booster_faster <- gamma * sum(Ivseries_booster_faster[(44):(4*noweeks)])
  Ivseries_primary <- out_primary$IV+out_primary$IS1+out_primary$IS2+out_primary$IS3
  Iv_primary <- gamma * sum(Ivseries_primary[(44):(4*noweeks)])




  totalgap <- (Itot_primary-Itot_booster_faster)*1e9
  primarygap <- (Ip_primary-Ip_booster_faster)*1e9
  secondarygap  <-  (Is_primary-Is_booster_faster)*1e9
  vaccinegap  <-  (Iv_primary-Iv_booster_faster)*1e9

  symp.cases=primarygap
  asymp.cases= vaccinegap + secondarygap *0.82203
  hospitalization.cases=primarygap*0.09789 + secondarygap*0.04237
  icu.cases=(primarygap*0.09789 + secondarygap*0.04237)*0.2
  death.cases= primarygap*0.01941+secondarygap*0.05932




  value <- list(totalgap,symp.cases,asymp.cases,hospitalization.cases,icu.cases,death.cases)

  return(value)
}

#nu=seq(from=0,to=1,length.out=gridsize)
#ImmuneTime_vax=seq(from=0.1*noweeks,to=2*noweeks,length.out=gridsize)
omega1=seq(from=0.0,to=0.1,length.out=gridsize)
#timeInterval=seq(from=0*noweeks,to=1*noweeks,length.out=gridsize)
rho3=seq(from=0.0,to=0.1,length.out=gridsize)


#nu=5e-3
#ImmuneTime_vax=0.5*noweeks


out_totalgap <- matrix(0,length(delta),length(rho3))
out_symp.cases <- matrix(0,length(delta),length(rho3))
out_asymp.cases <- matrix(0,length(delta),length(rho3))
out_hospitalization.cases <- matrix(0,length(delta),length(rho3))
out_icu.cases <- matrix(0,length(delta),length(rho3))
out_death.cases <- matrix(0,length(delta),length(rho3))

for(i in 1:length(delta)){
  for(j in 1:length(rho3)){
    temp <- allIvax(delta[i],rho3[j])
    out_totalgap[i,j] <- as.numeric(temp[[1]])
    out_symp.cases[i,j] <- as.numeric(temp[[2]])
    out_asymp.cases[i,j] <- as.numeric(temp[[3]])
    out_hospitalization.cases[i,j] <- as.numeric(temp[[4]])
    out_icu.cases[i,j] <- as.numeric(temp[[5]])
    out_death.cases[i,j] <- as.numeric(temp[[6]])
  }
}

##################################################
###########Plot with filled contour###############
##################################################



surface <- reshape2::melt(out_totalgap)
plot.total <- ggplot(surface, aes(Var1, Var2, z = value)) +
  geom_contour_fill(aes(fill = after_stat(level_d))) +
  geom_contour(color = "white", size = 0.1)+theme_classic()+
  #xlab(expression(paste(epsilon[3])))+ylab(expression(paste(epsilon[V][3])))+
  xlab("Rollout rate of booster doses")+ylab("VS")+
  scale_fill_brewer(palette = "BrBG", direction = - 1)+
  theme(text=element_text(size=24), #change font size of all text
        axis.text=element_text(size=24), #change font size of axis text
        axis.title=element_text(size=24), #change font size of axis titles
        plot.title=element_text(size=24), #change font size of plot title
        legend.text=element_text(size=24), #change font size of legend text
        legend.title=element_text(size=24))



surface <- reshape2::melt(out_symp.cases)
plot.symp <- ggplot(surface, aes(Var1, Var2, z = value)) +
  geom_contour_fill(aes(fill = after_stat(level_d))) +
  geom_contour(color = "white", size = 0.1)+theme_classic()+
  #xlab(expression(paste(epsilon[3])))+ylab(expression(paste(epsilon[V][3])))+
  xlab("Rollout rate of booster doses")+ylab("VS")+
  scale_fill_brewer(palette = "PiYG", direction = - 1)+
  theme(text=element_text(size=24), #change font size of all text
        axis.text=element_text(size=24), #change font size of axis text
        axis.title=element_text(size=24), #change font size of axis titles
        plot.title=element_text(size=24), #change font size of plot title
        legend.text=element_text(size=24), #change font size of legend text
        legend.title=element_text(size=24))+
  labs(fill='Symptomatic cases\n averted:')



surface <- reshape2::melt(out_asymp.cases)
plot.asymp <- ggplot(surface, aes(Var1, Var2, z = value)) +
  geom_contour_fill(aes(fill = after_stat(level_d))) +
  geom_contour(color = "white", size = 0.1)+theme_classic()+
  #xlab(expression(paste(epsilon[3])))+ylab(expression(paste(epsilon[V][3])))+
  xlab("Rollout rate of booster doses")+ylab("VS")+
  scale_fill_brewer(palette = "PuOr", direction = - 1)+
  theme(text=element_text(size=24), #change font size of all text
        axis.text=element_text(size=24), #change font size of axis text
        axis.title=element_text(size=24), #change font size of axis titles
        plot.title=element_text(size=24), #change font size of plot title
        legend.text=element_text(size=24), #change font size of legend text
        legend.title=element_text(size=24))+
  labs(fill='Asymptomatic cases\n averted:')



surface <- reshape2::melt(out_hospitalization.cases)
plot.hosp <- ggplot(surface, aes(Var1, Var2, z = value)) +
  geom_contour_fill(aes(fill = after_stat(level_d))) +
  geom_contour(color = "white", size = 0.1)+theme_classic()+
  #xlab(expression(paste(epsilon[3])))+ylab(expression(paste(epsilon[V][3])))+
  xlab("Rollout rate of booster doses")+ylab("VS")+
  scale_fill_brewer(palette ="PRGn", direction = - 1)+
  theme(text=element_text(size=24), #change font size of all text
        axis.text=element_text(size=24), #change font size of axis text
        axis.title=element_text(size=24), #change font size of axis titles
        plot.title=element_text(size=24), #change font size of plot title
        legend.text=element_text(size=24), #change font size of legend text
        legend.title=element_text(size=24))+
  labs(fill='Hospitalized cases\n averted:')





surface <- reshape2::melt(out_icu.cases)
plot.icu <- ggplot(surface, aes(Var1, Var2, z = value)) +
  geom_contour_fill(aes(fill = after_stat(level_d))) +
  geom_contour(color = "white", size = 0.1)+theme_classic()+
  #xlab(expression(paste(epsilon[3])))+ylab(expression(paste(epsilon[V][3])))+
  xlab("Rollout rate of booster doses")+ylab("VS")+
  scale_fill_brewer(palette ="PRGn", direction = - 1)+
  theme(text=element_text(size=24), #change font size of all text
        axis.text=element_text(size=24), #change font size of axis text
        axis.title=element_text(size=24), #change font size of axis titles
        plot.title=element_text(size=24), #change font size of plot title
        legend.text=element_text(size=24), #change font size of legend text
        legend.title=element_text(size=24))




surface <- reshape2::melt(out_death.cases)
plot.death <- ggplot(surface, aes(Var1, Var2, z = value)) +
  geom_contour_fill(aes(fill = after_stat(level_d))) +
  geom_contour(color = "white", size = 0.1)+theme_classic()+
  #xlab(expression(paste(epsilon[3])))+ylab(expression(paste(epsilon[V][3])))+
  xlab("Rollout rate of booster doses")+ylab("VS")+
  scale_fill_brewer(palette ="PRGn", direction = - 1)+
  theme(text=element_text(size=24), #change font size of all text
        axis.text=element_text(size=24), #change font size of axis text
        axis.title=element_text(size=24), #change font size of axis titles
        plot.title=element_text(size=24), #change font size of plot title
        legend.text=element_text(size=24), #change font size of legend text
        legend.title=element_text(size=24))



ggarrange(plot.total,plot.symp,plot.asymp,
          plot.hosp,plot.icu,plot.death,
          labels = c("A","B","C","D","E","F"),
          font.label = list(size = 24),ncol = 2, nrow = 3)
ggsave("Figtest.pdf",width = 20,height = 18)













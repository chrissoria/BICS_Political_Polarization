# adapted from Harris et al.
# https://github.com/mjharris95/divided-disease/blob/main/aware-eqs.R
# https://www.cambridge.org/core/journals/evolutionary-human-sciences/article/social-divisions-and-risk-perception-drive-divergent-epidemics-and-large-later-waves/71613CAF7D03C9F20A8864A59B9BF311#article

library(deSolve)
library(tidyverse)
library(reshape2)

full_aware<-function(h=NA, epsilon=NA, mu=NA, kappa=NA, beta=NA, rho=NA, sigma = NA, omega = NA,
                     theta=5, phi=.1, kappa_a=.1,  kappa_b=.1, rho_a=.1,  
                     mu_a=.05, mu_b=.05, beta_a=.5, beta_b=.5, rho_b=.1,
                     ha=.99, hb=.99, epsilon_a=.5, epsilon_b=.5, sigma_a = .5, sigma_b = .5, 
                     ell=1, I0_a=1, I0_b=1, N0_a = 0.8*10000000, N0_b = 0.2*10000000,
                     v_val=0, v_start=0, time=200, theta_a=.2, theta_b=.2, omega_a = .5, omega_b = .5,
                     get_params=FALSE){
  if(!is.na(beta)){
    beta_a<-beta
    beta_b<-beta
  }
  
  if(!is.na(kappa)){
    kappa_a<-kappa
    kappa_b<-kappa
  }
  
  if(!is.na(mu)){
    mu_a<-mu
    mu_b<-mu
  }
  
  if(!is.na(rho)){
    rho_a<-rho
    rho_b<-rho
  }
  
  if(!is.na(h)){
    ha<-h
    hb<-h
  }
  
  if(!is.na(epsilon)){
    epsilon_a<-epsilon
    epsilon_b<-epsilon
  }
  
  if(!is.na(sigma)){
    sigma_a<-sigma
    sigma_b<-sigma
  }
  
  if(!is.na(theta)){
    theta_a<-theta
    theta_b<-theta
  }
  
  if(!is.na(omega)){
    omega_a<-omega
    omega_b<-omega
  }
  
  N_a = N0_a
  N_b = N0_b
  
  params<-c("kappa_a"=kappa_a, "kappa_b"=kappa_b,  "theta"=theta, 
            "beta_a"=beta_a, "beta_b"=beta_b, "rho"=rho, "phi"=phi, 
            "rho_a"=rho_a, "rho_b"=rho_b,
            "mu_a"=mu_a, "mu_b"=mu_b, "epsilon_a"=epsilon_a, "epsilon_b"=epsilon_b,  "ha"=ha, "hb"=hb, "ell"=ell, "v_start"=v_start, "v_val"=v_val, "time"=time, 
            "I0_a"=I0_a, "I0_b"=I0_b, 
            "N0_a" = N0_a, "N0_b"=N0_b,
            "theta_a"=theta_a, "theta_b"=theta_b,
            "sigma_a" = sigma_a, "sigma_b" = sigma_b,
            "omega_a" = omega_a, "omega_b" = omega_b) 
  state<-c( 
    
    #group a
    SUa<-(N0_a - I0_a)*1,
    SPa<-(N0_a - I0_a)*0,
    IUa<-I0_a,
    IPa<-0,
    RUa<-0,
    RPa<-0,
    DUa<-0,
    DPa<-0,

    #group b
    SUb<-(N0_b - I0_b)*1,
    SPb<-(N0_b - I0_b)*0,
    IUb<-I0_b,
    IPb<-0,
    RUb<-0,
    RPb<-0,
    DUb<-0,
    DPb<-0
  )
  
  names(state)<- c("SUa", "SPa", "IUa", "IPa", "RUa", "RPa", "DUa", "DPa",
                   "SUb", "SPb", "IUb", "IPb", "RUb", "RPb", "DUb", "DPb")
  
  sir_upu<-function(t, state, parameter){
    
    with(as.list(c(parameter, state)), {
      
      
      if(t<ell){
        lag<-rep(0, 16)
      }
      else{
        lag<-lagvalue(t-ell)
      }
      
      #transitions for group a
      if(t>v_start){
        v <- v_val
      }
      else{
        v<-0
      }
      

      dSUa <- -sqrt(beta_a)*SUa*(sqrt(beta_a)*ha*(IUa/N_a+kappa_a*IPa/N_a)+sqrt(beta_b)*(1-ha)*(IUb/N_b+kappa_b*IPb/N_b)) - 
        theta_a * SUa * (epsilon_a * (DUa+DPa-lag[7]-lag[8]) + (1-epsilon_a) * (DUb+DPb-lag[15]-lag[16]) ) - 
        omega_a * SUa * (sigma_a *(SPa + IPa + RPa)/N_a + (1 - sigma_a)* (SPb + IPb + RPb)/N_b) + 
        phi * SPa  
      
      dSPa <- -sqrt(beta_a)*kappa_a*SPa*(sqrt(beta_a)*ha*(IUa/N_a+kappa_a*IPa/N_a)+sqrt(beta_b)*(1-ha)*(IUb/N_b+kappa_b*IPb/N_b)) + 
        theta_a * SUa * (epsilon_a * (DUa+DPa-lag[7]-lag[8]) + (1-epsilon_a) * (DUb+DPb-lag[15]-lag[16]) ) +  
        omega_a * SUa * (sigma_a *(SPa + IPa + RPa)/N_a + (1 - sigma_a)* (SPb + IPb + RPb)/N_b) - 
        phi * SPa  - v*SPa
      
      dIUa <- sqrt(beta_a)*SUa*(sqrt(beta_a)*ha*(IUa/N_a+kappa_a*IPa/N_a)+sqrt(beta_b)*(1-ha)*(IUb/N_b+kappa_b*IPb/N_b)) - 
        theta_a * IUa * (epsilon_a * (DUa+DPa-lag[7]-lag[8]) + (1-epsilon_a) * (DUb+DPb-lag[15]-lag[16]) ) -
        omega_a * IUa * (sigma_a *(SPa + IPa + RPa)/N_a + (1 - sigma_a)* (SPb + IPb + RPb)/N_b) + 
        phi * IPa -
        (rho_a) * IUa
      
      dIPa <- sqrt(beta_a)*kappa_a*SPa*(sqrt(beta_a)*ha*(IUa/N_a+kappa_a*IPa/N_a)+sqrt(beta_b)*(1-ha)*(IUb/N_b+kappa_b*IPb/N_b)) + 
        theta_a * IUa * (epsilon_a * (DUa+DPa-lag[7]-lag[8]) + (1-epsilon_a) * (DUb+DPb-lag[15]-lag[16])  ) +
        omega_a * IUa * (sigma_a *(SPa + IPa + RPa)/N_a + (1 - sigma_a)* (SPb + IPb + RPb)/N_b)- 
        phi * IPa -
        (rho_a) * IPa
      
      
      #dSUa <- -beta*SUa*(C_aa*(IUa+kappa_a*IPa)/N_a + C_ba*(IUb+kappa_b*IPb)/N_b) - 
      #  theta_a * SUa * (epsilon_a * (DUa+DPa-lag[7]-lag[8]) + (1-epsilon_a) * (DUb+DPb-lag[16]-lag[17]) ) + phi * SPa  
      
      #dSPa <- -beta*kappa_a*SPa*(C_aa*(IUa+kappa_a*IPa)/N_a + C_ba*(IUb +kappa_b*IPb)/N_b) + 
      #  theta_a * SUa * (epsilon_a * (DUa+DPa-lag[7]-lag[8]) + (1-epsilon_a) * (DUb+DPb-lag[16]-lag[17]) ) - phi * SPa  -
      #  v*SPa
      
      #dIUa <- beta*SUa*(C_aa*(IUa + kappa_a*IPa)/N_a + C_ba*(IUb + kappa_b*IPb)/N_b) - 
      #  theta_a * IUa * (epsilon_a * (DUa+DPa-lag[7]-lag[8]) + (1-epsilon_a) * (DUb+DPb-lag[16]-lag[17]) ) + phi * IPa -
      #  (rho_a) * IUa
      
      #dIPa <- beta*kappa_a*SPa*(C_aa*(IUa + kappa_a*IPa)/N_a + C_ba*(IUb + kappa_b*IPb)/N_b) + 
      #  theta_a * IUa * (epsilon_a * (DUa+DPa-lag[7]-lag[8]) + (1-epsilon_a) * (DUb+DPb-lag[16]-lag[17])  ) - phi * IPa -
      #  (rho_a) * IPa
      
      dRUa <- rho_a * (1-mu_a) * IUa - 
        theta_a * RUa * (epsilon_a * (DUa+DPa-lag[7]-lag[8]) + (1-epsilon_a) * (DUb+DPb-lag[15]-lag[16]) ) -
        omega_a * RUa * (sigma_a *(SPa + IPa + RPa)/N_a + (1 - sigma_a)* (SPb + IPb + RPb)/N_b) + phi * RPa  
      
      dRPa <- rho_a * (1-mu_a) * IPa + v * SPa +
        theta_a * RUa * (epsilon_a * (DUa+DPa-lag[7]-lag[8]) + (1-epsilon_a) * (DUb+DPb-lag[15]-lag[16])  ) +
        omega_a * RUa * (sigma_a *(SPa + IPa + RPa)/N_a + (1 - sigma_a)* (SPb + IPb + RPb)/N_b) - phi * RPa  
      
      dDUa <- rho_a*mu_a*IUa
      
      dDPa <- rho_a*mu_a*IPa

      
      #transitions for group b
      
      dSUb <- -sqrt(beta_b)*SUb*(sqrt(beta_a)*(1-hb)*(IUa/N_a+kappa_a*IPa/N_a)+sqrt(beta_b)*(hb)*(IUb/N_b+kappa_b*IPb/N_b)) - 
        theta_b * SUb * ((1-epsilon_b) * (DUa+DPa-lag[7]-lag[8]) + (epsilon_b) * (DUb+DPb-lag[15]-lag[16]) ) -
        omega_b * SUb * ((1- sigma_b) *(SPa + IPa + RPa)/N_a + (sigma_b)* (SPb + IPb + RPb)/N_b) + phi * SPb
      
      dSPb <- -sqrt(beta_b)*kappa_b*SPb*(sqrt(beta_a)*(1-hb)*(IUa/N_a+kappa_a*IPa/N_a)+sqrt(beta_b)*(hb)*(IUb/N_b+kappa_b*IPb/N_b)) + 
        theta_b * SUb * ((1-epsilon_b) * (DUa+DPa-lag[7]-lag[8]) + (epsilon_b) * (DUb+DPb-lag[15]-lag[16])  ) +
        omega_b * SUb * ((1- sigma_b) *(SPa + IPa + RPa)/N_a + (sigma_b)* (SPb + IPb + RPb)/N_b) - phi * SPb -
        v * SPb
      
      dIUb <- sqrt(beta_b)*SUb*(sqrt(beta_a)*(1-hb)*(IUa/N_a+kappa_a*IPa/N_a)+sqrt(beta_b)*(hb)*(IUb/N_b+kappa_b*IPb/N_b)) - 
        theta_b * IUb * ((1-epsilon_b) * (DUa+DPa-lag[7]-lag[8]) + (epsilon_b) * (DUb+DPb-lag[15]-lag[16]) ) -
        omega_b * IUb * ((1- sigma_b) *(SPa + IPa + RPa)/N_a + (sigma_b)* (SPb + IPb + RPb)/N_b) + phi * IPb -
        (rho_b) * IUb
      
      dIPb <- sqrt(beta_b)*kappa_b*SPb*(sqrt(beta_a)*(1-hb)*(IUa/N_a+kappa_a*IPa/N_a)+sqrt(beta_b)*(hb)*(IUb/N_b+kappa_b*IPb/N_b)) +
        theta_b * IUb * ((1-epsilon_b) * (DUa+DPa-lag[7]-lag[8]) + (epsilon_b) * (DUb+DPb-lag[15]-lag[16]) ) + 
        omega_b * IUb * ((1- sigma_b) *(SPa + IPa + RPa)/N_a + (sigma_b)* (SPb + IPb + RPb)/N_b)- phi * IPb -
        (rho_b) * IPb
      
      dRUb <- (1-mu_b) * rho_b * IUb - 
        theta_b * RUb * ((1-epsilon_b) * (DUa+DPa-lag[7]-lag[8]) + (epsilon_b) * (DUb+DPb-lag[15]-lag[16]) ) -
        omega_b * RUb * ((1- sigma_b) *(SPa + IPa + RPa)/N_a + (sigma_b)* (SPb + IPb + RPb)/N_b) + phi * RPb
      
      dRPb <- (1-mu_b) * rho_b * IPb + v * SPb +
        theta_b * RUb * ((1-epsilon_b) * (DUa+DPa-lag[7]-lag[8]) + (epsilon_b) * (DUb+DPb-lag[15]-lag[16])  ) +
        omega_b * RUb * ((1- sigma_b) *(SPa + IPa + RPa)/N_a + (sigma_b)* (SPb + IPb + RPb)/N_b) - phi * RPb
      
      dDUb <- rho_b*mu_b*IUb
      
      dDPb <- rho_b*mu_b*IPb
      

      return(list(c(dSUa, dSPa, dIUa, dIPa, dRUa, dRPa, dDUa, dDPa,
                    dSUb, dSPb, dIUb, dIPb, dRUb, dRPb, dDUb, dDPb
      )))
    })
  }
  
  times<-seq(from=0, to=time, by=1)
  
  as.data.frame(dede(state, times, sir_upu, params))->sim
  
  if(get_params){
    return(list(sim=sim, params=params))
  }
  
  else{
    return(sim)
  }
  
}


f = 0.5

N_a = f*10000000; N_b = (1-f)*10000000

C_a = 5; C_b = 5
ls<-full_aware(h=.5,  mu=.01, kappa=0.03, beta_a = 0.2*C_a, beta_b = 0.2*C_b, rho=1/10, 
               theta_a=0.3, theta_b = 0.001, epsilon=.5,
               omega_a = 0.3, omega_b = 0.001, sigma = .5,
               N0_a = N_a, N0_b = N_b, time = 500,
               ell=1, phi=0, get_params=TRUE)
sim<-ls$sim




plot(SUa/N_a~time, data=sim, type="l", col="#3D85BD", ylab="Susceptible", main="a")
lines(SPa/N_a~time, data=sim, col="#3D85BD", lty=2)
plot(SUb/N_b~time, data=sim, col="#3D85BD", type="l", main="b", ylab="Susceptible")
lines(SPb/N_b~time, data=sim, col="#3D85BD", lty=2)


plot(IUa/N_a~time, data=sim, type="l", col="blue", ylim = c(0, 0.0001), ylab="Infected")
lines(IPa/N_a~time, data=sim, col="blue", lty=2)
lines(IUb/N_b~time, data=sim, col="red", type="l", ylab="Infected")
lines(IPb/N_b~time, data=sim, col="red", lty=2)


plot((IUa + IPa)~time, data=sim, type="l", col="blue", ylab="Infected")
plot((IUb+IPb)~time, data=sim, type="l", col="blue", ylab="Infected")

plot((IUa+IPa)/N_a~time, data=sim, type="l", col="blue", ylim = c(0, 0.0001), ylab="Infected")
lines((IUb+IPb)/N_b~time, data=sim, col="red", type="l", ylab="Infected")

plot((DUa+DPa)~time, data=sim, type="l", col="blue", ylab="Deaths")
plot((DUb+DPb)~time, data=sim, type="l", col="blue", ylab="Deaths")


plot((SPa+IPa+RPa)/N_a~time, data=sim, type="l", col="blue", ylab="Protected")
lines((SPb+IPb+RPb)/N_b~time, data = sim, col = "red")

#can change the values of h and epsilon. Both vectors can be any length.
h<-c(.5, .99)
epsilon<-c(.5, .99)

expand.grid(h=h, epsilon=epsilon)  %>%
  data.frame() %>%
  apply(1, function(par) full_aware(theta=100, beta=.2, rho=.1, mu=.01, h=par[["h"]],time=200, kappa=.3, I0_a=.001, I0_b=0, ell=1, phi=0, epsilon=par[["epsilon"]]) %>% 
          cbind(h=par[["h"]], epsilon=par[["epsilon"]])) %>% 
  do.call(rbind, .) %>%
  mutate(eps_descript = ifelse(epsilon==.5, "Uniform Awareness", 
                               ifelse(epsilon==.99,"Separated Awareness",
                                      " "))) %>%
  mutate(h_descript = ifelse(h==.5, "Uniform \nMixing", 
                             ifelse(h==.99,"Separated \nMixing",
                                    " "))) -> df

line_sz<-.8

ggplot(df)+geom_line(aes(x=time, y=IUa+IPa), size=line_sz, color="#f4a896")+geom_line(aes(x=time, y=IUb+IPb, size=as.factor(h)), color="#358597")+
  annotate("segment", x=-Inf, xend=Inf, y=-Inf, yend=-Inf, size=1)+
  annotate("segment", x=-Inf, xend=-Inf, y=-Inf, yend=Inf, size=1)+
  geom_hline(yintercept=c(0,.04,.08, .12), linetype="dotted", color="gray")+
  facet_grid(h+h_descript~epsilon+eps_descript, 
             labeller=label_bquote(
               cols=atop(atop(phantom(), .(eps_descript)), atop("("~epsilon==.(epsilon)~")", phantom())),
               rows=atop(atop(.(h_descript), ""), atop("("~h==.(h)~")", phantom()))))+
  ylab("Infections")+
  theme(strip.text.y = element_text(angle = 0))+scale_y_continuous(breaks=c(0, .04, .08, .12))+
  scale_size_manual(values=c("0.5"=.3*line_sz, "0.99"=line_sz), guide="none", na.value = line_sz)+
  scale_colour_manual(name="", 
                      labels = c(expression("group"~italic(a)), expression("group"~italic(b))),
                      values = c("group *a*"="#f4a896", "group *b*"="#358597"))+
  theme(strip.text=element_text(size=20))




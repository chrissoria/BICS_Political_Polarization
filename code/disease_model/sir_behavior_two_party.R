# adapted from Harris et al.
# https://github.com/mjharris95/divided-disease/blob/main/aware-eqs.R
# https://www.cambridge.org/core/journals/evolutionary-human-sciences/article/social-divisions-and-risk-perception-drive-divergent-epidemics-and-large-later-waves/71613CAF7D03C9F20A8864A59B9BF311#article

library(deSolve)
library(tidyverse)
library(reshape2)

# two group SIR model (A and B); each compartment is split into "Protected (P)" and "Unprotected (U)"
# both groups can adopt protective behavior (i.e. mask usage) in two ways: 1) in response to deaths over a specified time window; 2) through contact with other "protected" individuals
# protective behavior is protective of both transmission and infection
# groups can differ in their rate of contact and their propensity to adopt protective behaviors
# fraction of the population in each group can vary
# homophily (degree to which members of a group mix within group versus outside their group) can vary


sir_two_group_pu <- function(C = NA, C_a = 3, C_b = 3, 
                          trans_p = 0.05, rho=1/10, mu = 0.01, 
                          h_a=.5, kappa=0.03, phi = 0,
                          I0_a=1, I0_b=1, N0 = 10000000, frac_a = 0.5, time = 500,
                          ell = 1, theta = NA, theta_a = 100, theta_b  = 100, 
                          epsilon = NA, epsilon_a = 0.5, epsilon_b = 0.5,
                          omega = NA, omega_a = 0.01, omega_b = 0.01,
                          get_params=FALSE) {
  
  
  if(!is.na(epsilon)){
    epsilon_a<-epsilon
    epsilon_b<-epsilon
  }
  
  if(!is.na(theta)){
    theta_a<-theta
    theta_b<-theta
  }
  
  if(!is.na(omega)){
    omega_a<-omega
    omega_b<-omega
  }
  
  
  if(!is.na(C)){
    C_a<-C
    C_b<-C
  }
  
  N0_a = N0 * frac_a #letting population be constant for now
  N0_b = N0 * (1 - frac_a)
  
  # create contact matrix entries : C_ba = contacts members of a have with members of b
  C_aa = C_a * h_a
  C_ba = C_a * (1 - h_a)
  # by reciprocity
  C_ab = (frac_a * C_ba ) / (1 - frac_a)
  C_bb = max(0, C_b - C_ab)
  
  N_a = N0_a
  N_b = N0_b
  
  #print(c(C_aa, C_ba, C_ab, C_bb))
  
  params<-c("trans_p" = trans_p, 
            "h_a"= h_a, 
            "C_a" = C_a, "C_b" = C_b,
            "C_aa" = C_aa, "C_ab" = C_ab,
            "C_ba" = C_ba, "C_bb" = C_bb,
            "kappa" = kappa,
            "phi"=phi, 
            "rho"=rho, 
            "mu"=mu, "time"=time, 
            "I0_a"=I0_a, "I0_b"=I0_b, 
            "N0_a" = N0_a, "N0_b"= N0_b,
            "ell"=ell,
            "epsilon_a"=epsilon_a, "epsilon_b"=epsilon_b,
            "theta_a"=theta_a, "theta_b"=theta_b) 
  
  state<-c( 
    
    #group a
    SUa<-(N0_a - I0_a),
    SPa<-0,
    IUa<-I0_a,
    IPa<-0,
    RUa<-0,
    RPa<-0,
    DUa<-0,
    DPa<-0,
    
    #group b
    SUb<-N0_b - I0_b,
    SPb<-0,
    IUb<-I0_b,
    IPb<-0,
    RUb<-0,
    RPb<-0,
    DUb<-0,
    DPb<-0
  )
  
  names(state)<- c("SUa", "SPa", "IUa", "IPa", "RUa", "RPa", "DUa", "DPa",
                   "SUb", "SPb", "IUb", "IPb", "RUb", "RPb", "DUb", "DPb")
  
  sir_up<-function(t, state, parameter){
    
    with(as.list(c(parameter, state)), {
      
      
      if(t<ell){
        lag<-rep(0, length(state))
      }
      else{
        lag<-lagvalue(t-ell)
      }
      
      dSUa <- -SUa*trans_p*(C_aa*(IUa/N_a + kappa*IPa/N_a) + C_ba*(IUb/N_b + kappa*IPb/N_b)) + 
        phi * SPa  - 
        theta_a * SUa * (epsilon_a * ((DUa+DPa-lag[which(names(state) == "DUa")]-lag[which(names(state) == "DPa")])/N_a) + 
                           (1-epsilon_a) * ((DUb+DPb-lag[which(names(state) == "DUb")]-lag[which(names(state) == "DPb")])/N_b) ) -
        omega_a * SUa * (C_aa *((SPa + IPa + RPa)/N_a) + (C_ba)*((SPb + IPb + RPb)/N_b))
      
 
      dSPa <- -SPa*kappa*trans_p*(C_aa*(IUa/N_a + kappa*IPa/N_a) + C_ba*(IUb/N_b + kappa*IPb/N_b)) - 
        phi * SPa +
        theta_a * SUa * (epsilon_a * ((DUa+DPa-lag[which(names(state) == "DUa")]-lag[which(names(state) == "DPa")])/N_a) + 
                           (1-epsilon_a) * ((DUb+DPb-lag[which(names(state) == "DUb")]-lag[which(names(state) == "DPb")])/N_b) ) +
        omega_a * SUa * (C_aa *((SPa + IPa + RPa)/N_a) + (C_ba)*((SPb + IPb + RPb)/N_b))
      
      dIUa <- SUa*trans_p*(C_aa*(IUa/N_a + kappa*IPa/N_a) + C_ba*(IUb/N_b + kappa*IPb/N_b))  + 
        phi * IPa -
        (rho) * IUa -
        theta_a * IUa * (epsilon_a * ((DUa+DPa-lag[which(names(state) == "DUa")]-lag[which(names(state) == "DPa")])/N_a) + 
                           (1-epsilon_a) * ((DUb+DPb-lag[which(names(state) == "DUb")]-lag[which(names(state) == "DPb")])/N_b) ) -
        omega_a * IUa * (C_aa *((SPa + IPa + RPa)/N_a) + (C_ba)*((SPb + IPb + RPb)/N_b))
      
      dIPa <- SPa*trans_p*kappa*(C_aa*(IUa/N_a + kappa*IPa/N_a) + C_ba*(IUb/N_b + kappa*IPb/N_b)) - 
        phi * IPa -
        (rho) * IPa +
        theta_a * IUa * (epsilon_a * ((DUa+DPa-lag[which(names(state) == "DUa")]-lag[which(names(state) == "DPa")])/N_a) + 
                           (1-epsilon_a) * ((DUb+DPb-lag[which(names(state) == "DUb")]-lag[which(names(state) == "DPb")])/N_b) ) +
        omega_a * IUa * (C_aa *((SPa + IPa + RPa)/N_a) + (C_ba)*((SPb + IPb + RPb)/N_b))
      
      dRUa <- rho * (1-mu) * IUa  + phi * RPa  -
        theta_a * RUa * (epsilon_a * ((DUa+DPa-lag[which(names(state) == "DUa")]-lag[which(names(state) == "DPa")])/N_a) + 
                           (1-epsilon_a) * ((DUb+DPb-lag[which(names(state) == "DUb")]-lag[which(names(state) == "DPb")])/N_b) ) -
        omega_a * RUa * (C_aa *((SPa + IPa + RPa)/N_a) + (C_ba)*((SPb + IPb + RPb)/N_b))
      
      dRPa <- rho * (1-mu) * IPa  - phi * RPa  +
        theta_a * RUa * (epsilon_a * ((DUa+DPa-lag[which(names(state) == "DUa")]-lag[which(names(state) == "DPa")])/N_a) + 
                           (1-epsilon_a) * ((DUb+DPb-lag[which(names(state) == "DUb")]-lag[which(names(state) == "DPb")])/N_b) ) +
        omega_a * RUa * (C_aa *((SPa + IPa + RPa)/N_a) + (C_ba)*((SPb + IPb + RPb)/N_b))
      
      dDUa <- rho*mu*IUa
      
      dDPa <- rho*mu*IPa
      
      
      #transitions for group b
      
      dSUb <- -SUb*trans_p*(C_ab*(IUa/N_a + kappa*IPa/N_a) + C_bb*(IUb/N_b + kappa*IPb/N_b))  + phi * SPb - 
        theta_b * SUb * ((1-epsilon_b) * ((DUa+DPa-lag[which(names(state) == "DUa")]-lag[which(names(state) == "DPa")])/N_a) + 
                           (epsilon_b) * ((DUb+DPb-lag[which(names(state) == "DUb")]-lag[which(names(state) == "DPb")]) / N_b)) -
        omega_b * SUb * (C_ab *((SPa + IPa + RPa)/N_a) + (C_bb)*((SPb + IPb + RPb)/N_b))
      
      dSPb <- -trans_p*kappa*SPb*(C_ab*(IUa/N_a + kappa*IPa/N_a) + C_bb*(IUb/N_b + kappa*IPb/N_b))  - phi * SPb +
        theta_b * SUb * ((1-epsilon_b) * ((DUa+DPa-lag[which(names(state) == "DUa")]-lag[which(names(state) == "DPa")])/N_a) + 
                           (epsilon_b) * ((DUb+DPb-lag[which(names(state) == "DUb")]-lag[which(names(state) == "DPb")]) / N_b)) +
        omega_b * SUb * (C_ab *((SPa + IPa + RPa)/N_a) + (C_bb)*((SPb + IPb + RPb)/N_b))
      
      dIUb <- trans_p*SUb*(C_ab*(IUa/N_a + kappa*IPa/N_a) + C_bb*(IUb/N_b + kappa*IPb/N_b)) + phi * IPb -
        (rho) * IUb -
        theta_b * IUb * ((1-epsilon_b) * ((DUa+DPa-lag[which(names(state) == "DUa")]-lag[which(names(state) == "DPa")])/N_a) + 
                           (epsilon_b) * ((DUb+DPb-lag[which(names(state) == "DUb")]-lag[which(names(state) == "DPb")]) / N_b)) -
        omega_b * IUb * (C_ab *((SPa + IPa + RPa)/N_a) + (C_bb)*((SPb + IPb + RPb)/N_b))
      
      dIPb <- trans_p*kappa*SPb*(C_ab*(IUa/N_a + kappa*IPa/N_a) + C_bb*(IUb/N_b + kappa*IPb/N_b)) - phi * IPb -
        (rho) * IPb +
        theta_b * IUb * ((1-epsilon_b) * ((DUa+DPa-lag[which(names(state) == "DUa")]-lag[which(names(state) == "DPa")])/N_a) + 
                           (epsilon_b) * ((DUb+DPb-lag[which(names(state) == "DUb")]-lag[which(names(state) == "DPb")]) / N_b)) +
        omega_b * IUb * (C_ab *((SPa + IPa + RPa)/N_a) + (C_bb)*((SPb + IPb + RPb)/N_b))
      
      dRUb <- (1-mu) * rho * IUb + phi * RPb -
        theta_b * RUb * ((1-epsilon_b) * ((DUa+DPa-lag[which(names(state) == "DUa")]-lag[which(names(state) == "DPa")])/N_a) + 
                           (epsilon_b) * ((DUb+DPb-lag[which(names(state) == "DUb")]-lag[which(names(state) == "DPb")]) / N_b)) -
        omega_b * RUb * (C_ab *((SPa + IPa + RPa)/N_a) + (C_bb)*((SPb + IPb + RPb)/N_b))
      
      dRPb <- (1-mu) * rho * IPb  - phi * RPb +
        theta_b * RUb * ((1-epsilon_b) * ((DUa+DPa-lag[which(names(state) == "DUa")]-lag[which(names(state) == "DPa")])/N_a) + 
                           (epsilon_b) * ((DUb+DPb-lag[which(names(state) == "DUb")]-lag[which(names(state) == "DPb")]) / N_b)) +
        omega_b * RUb * (C_ab *((SPa + IPa + RPa)/N_a) + (C_bb)*((SPb + IPb + RPb)/N_b))
      
      dDUb <- rho*mu*IUb
      
      dDPb <- rho*mu*IPb
      
      
      return(list(c(dSUa, dSPa, dIUa, dIPa, dRUa, dRPa, dDUa, dDPa,
                    dSUb, dSPb, dIUb, dIPb, dRUb, dRPb, dDUb, dDPb
      )))
    })
  }
  
  times<-seq(from=0, to=time, by=1)
  
  as.data.frame(dede(state, times, sir_up, params))->sim
  
  if(get_params){
    return(list(sim=sim, params=params))
  }
  
  else{
    return(sim)
  }
  
  
  
}

#### test run ####

N0 = 10000000; frac_a = 0.5
N_a = frac_a*N0; N_b = N0 - N_a

ls<-sir_two_group_pu(C = NA, C_a = 5, C_b = 3, 
                     trans_p = 0.05, rho=1/10, mu = 0.01, 
                     h_a=0.5, kappa=0.3, phi = 0,
                     I0_a=1, I0_b=1, N0 = N0, frac_a = frac_a, time = 500,
                     theta_a = 20, theta_b = 100, epsilon = 0.5,
                     omega_a = 0.02, omega_b = 0.08,
                     get_params=TRUE)
sim<-ls$sim


plot(SUa/N_a~time, data=sim, type="l", col="#3D85BD", ylab="Susceptible", main="a", ylim = c(0,1))
lines(SPa/N_a~time, data=sim, col="#3D85BD", lty=2)
plot(SUb/N_b~time, data=sim, col="#3D85BD", type="l", main="b", ylab="Susceptible", ylim = c(0,1))
lines(SPb/N_b~time, data=sim, col="#3D85BD", lty=2)


plot(IUa/N_a~time, data=sim, type="l", col="red", ylab="Infected")
lines(IPa/N_a~time, data=sim, col="red", lty=2)
lines(IUb/N_b~time, data=sim, col="blue", type="l", ylab="Infected")
lines(IPb/N_b~time, data=sim, col="blue", lty=2)


plot((IUa + IPa)~time, data=sim, type="l", col="red", ylab="Infected")
plot((IUb+IPb)~time, data=sim, type="l", col="blue", ylab="Infected")

plot((IUa+IPa)/N_a~time, data=sim, type="l", col="red", ylab="Infected")
lines((IUb+IPb)/N_b~time, data=sim, col="blue", type="l", ylab="Infected")

#### simulate across parameters ####

# intial conditions
N0 = 10000000 # population size
frac_a = 0.5 # fraction in group A
N_a = frac_a*N0; N_b = N0 - N_a
I0_a = 1 #intial infected in A
I0_b = 1 # intial infected in B

# parameters
C_a = 5; C_b = 3 # average number of contacts per day 
trans_p = 0.05 # probability of transmission given contact
rho = 1/10 # 1 / infectious period
mu = 0.01 # probability of dying following infection
kappa = 0.3 # reducting in transmission and infection due to protective behavior
phi = 0 # waning rate of protective behavior

time = 365 # time steps for simulation
theta_a = 100 # responsiveness to deaths for adopting protective behavior in group A
theta_b = 100 # responsiveness to deaths for adopting protective behavior in group B 
omega_a = 0.02 # group A's rate of adopting protective behavior based on contact with protected individuals 
omega_b = 0.02 # group B's rate of adopting protective behavior based on contact with protected individuals 

h_a = c(.5, .99) # measure of homiphily in group A
epsilon = c(.5, .99) #


expand.grid(h_a=h_a, epsilon=epsilon)  %>%
  data.frame() %>%
  apply(1, function(par) sir_two_group_pu(C = NA, C_a = C_a, C_b = C_b, 
                                          trans_p = trans_p, rho=rho, mu = mu, 
                                          h_a=par[["h_a"]],
                                          kappa=kappa, phi = phi,
                                          I0_a=I0_a, I0_b=I0_b, N0 = N0, frac_a = frac_a, 
                                          time = time,
                                          theta_a = theta_a, 
                                          theta_b = theta_b, 
                                          epsilon=par[["epsilon"]],
                                          omega_a = omega_a, omega_b = omega_b) %>% 
          cbind(h_a=par[["h_a"]], epsilon=par[["epsilon"]])) %>% 
  do.call(rbind, .) %>%
  mutate(eps_descript = ifelse(epsilon==.5, "Uniform Awareness", 
                               ifelse(epsilon==.99,"Separated Awareness",
                                      " "))) %>%
  mutate(h_descript = ifelse(h_a==.5, "Uniform \nMixing", 
                             ifelse(h_a==.99,"Separated \nMixing",
                                    " "))) -> df

line_sz<-.8

ggplot(df) +
  geom_line(aes(x=time, y=(IUa+IPa)/N_a), size=line_sz, color="red") + 
  geom_line(aes(x=time, y=(IUb+IPb)/N_b, size=as.factor(h_a)), color="blue")+
  annotate("segment", x=-Inf, xend=Inf, y=-Inf, yend=-Inf, size=1)+
  annotate("segment", x=-Inf, xend=-Inf, y=-Inf, yend=Inf, size=1)+
  geom_hline(yintercept=c(0,.04,.08, .12), linetype="dotted", color="gray")+
  facet_grid(h_a+h_descript~epsilon+eps_descript, 
             labeller=label_bquote(
               cols=atop(atop(phantom(), .(eps_descript)), atop("("~epsilon==.(epsilon)~")", phantom())),
               rows=atop(atop(.(h_descript), ""), atop("("~h_a==.(h_a)~")", phantom()))))+
  ylab("Infections")+
  theme(strip.text.y = element_text(angle = 0))+scale_y_continuous(breaks=c(0, .04, .08, .12))+
  scale_size_manual(values=c("0.5"=.3*line_sz, "0.99"=line_sz), guide="none", na.value = line_sz)+
  scale_colour_manual(name="", 
                      labels = c(expression("group"~italic(a)), expression("group"~italic(b))),
                      values = c("group *a*"="red", "group *b*"="blue"))+
  theme(strip.text=element_text(size=20))


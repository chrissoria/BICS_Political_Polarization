# adapted from Harris et al.
# https://github.com/mjharris95/divided-disease/blob/main/aware-eqs.R
# https://www.cambridge.org/core/journals/evolutionary-human-sciences/article/social-divisions-and-risk-perception-drive-divergent-epidemics-and-large-later-waves/71613CAF7D03C9F20A8864A59B9BF311#article

library(deSolve)
library(tidyverse)
library(reshape2)

# two group SIR model (A and B); each compartment is split into "Protected (P)" and "Unprotected (U)"
# both groups can adopt protective behavior (i.e. mask usage) in two ways: 
#     1) in response to deaths over a specified time window; 2) in response to the proportion of "protected" individuals in the population
#     at the moment these are both linear relationships, i.e. individuals are more likely to adopt protective behavior if there are more deaths and there are a greater proportion of protected individuals in the population
# protective behavior is protective of both transmission and infection
# groups can differ in their rate of contact and their propensity to adopt protective behaviors
# fraction of the population in each group can vary
# preference for homophily (degree to which members of a group mix within group versus outside their group) can vary


sir_two_group_pu <- function(c = NA, c_a = 3, c_b = 3, 
                          trans_p = 0.05, rho=1/10, mu = 0.01, 
                          h_a=.5, kappa=0.03, phi = 0,
                          I0_a=1, I0_b=1, N0 = 10000000, frac_a = 0.5, time = 500,
                          ell = 1, theta = NA, theta_a = 100, theta_b  = 100, 
                          epsilon = NA, epsilon_a = 0.5, epsilon_b = 0.5,
                          omega = NA, omega_a = 0.1, omega_b = 0.1,
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
  
  
  if(!is.na(c)){
    c_a<-c
    c_b<-c
  }
  
  N0_a = N0 * frac_a #letting population be constant for now
  N0_b = N0 * (1 - frac_a)
  
  # create contact matrix entries : C_ba = total contacts members of a have with members of b
  C_a = c_a * N0_a
  C_b = c_b *N0_b
  
  if (C_a * (1 - h_a) < C_b) {
    C_ba = C_ab = C_a * (1 - h_a)
  } else {
    C_ba = C_ab = C_b
  }
  
  C_aa = C_a - C_ba
  C_bb = max(0, C_b - C_ab)
  
  # get averages
  
  c_aa = C_aa / N0_a
  c_ba = C_ba / N0_a
  c_bb = C_bb/ N0_b
  c_ab = C_ab / N0_b
  
  
  N_a = N0_a
  N_b = N0_b
  

 
  params<-c("trans_p" = trans_p, # probability of transmission given contact
            "h_a"= h_a, # # proportion of group A's total contact with members of their own group (bounded by population size and total contacts in group B)
            "c_a" = c_a, "c_b" = c_b, # average number of contacts per day 
            "c_aa" = c_aa, "c_ab" = c_ab, # average number of contacts per day 
            "c_ba" = c_ba, "c_bb" = c_bb, # average number of contacts per day 
            "kappa" = kappa, # reduction in probability of transmission given contact resulting from protective behavior (masking, distancing)
            "phi"=phi, # waning rate of protective behavior
            "rho"=rho, # 1 / infectious period
            "mu"=mu, # probability of dying following infection
            "time"=time, # time steps in simulation
            "I0_a"=I0_a, "I0_b"=I0_b, # starting number infected in each group
            "N0_a" = N0_a, "N0_b"= N0_b, # population size in each group
            "ell"=ell, # time window for considering deaths that influence adoption of protective behavior
            "epsilon_a"=epsilon_a, "epsilon_b"=epsilon_b, #measure of assortativeness in influence (i.e. are people equally aware and influenced by deaths/number of protected individuals in their own group versus in the out group)
            "theta_a"=theta_a, "theta_b"=theta_b, # responsiveness to deaths for adopting protective behavior 
            "omega_a"=omega_a, "omega_b"=omega_b) # responsiveness to proportion of protected individuals for adopting protective behavior
  
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
      
      dSUa <- -SUa*trans_p*(c_aa*(IUa/N_a + kappa*IPa/N_a) + c_ba*(IUb/N_b + kappa*IPb/N_b)) + 
        phi * SPa  - 
        theta_a * SUa * (epsilon_a * ((DUa+DPa-lag[which(names(state) == "DUa")]-lag[which(names(state) == "DPa")])/N_a) + 
                           (1-epsilon_a) * ((DUb+DPb-lag[which(names(state) == "DUb")]-lag[which(names(state) == "DPb")])/N_b) ) -
        omega_a * SUa * (epsilon_a *((SPa + IPa + RPa)/N_a) + (1-epsilon_a)*((SPb + IPb + RPb)/N_b))
      
 
      dSPa <- -SPa*kappa*trans_p*(c_aa*(IUa/N_a + kappa*IPa/N_a) + c_ba*(IUb/N_b + kappa*IPb/N_b)) - 
        phi * SPa +
        theta_a * SUa * (epsilon_a * ((DUa+DPa-lag[which(names(state) == "DUa")]-lag[which(names(state) == "DPa")])/N_a) + 
                           (1-epsilon_a) * ((DUb+DPb-lag[which(names(state) == "DUb")]-lag[which(names(state) == "DPb")])/N_b) ) +
        omega_a * SUa * (epsilon_a *((SPa + IPa + RPa)/N_a) + (1 - epsilon_a)*((SPb + IPb + RPb)/N_b))
      
      dIUa <- SUa*trans_p*(c_aa*(IUa/N_a + kappa*IPa/N_a) + c_ba*(IUb/N_b + kappa*IPb/N_b))  + 
        phi * IPa -
        (rho) * IUa -
        theta_a * IUa * (epsilon_a * ((DUa+DPa-lag[which(names(state) == "DUa")]-lag[which(names(state) == "DPa")])/N_a) + 
                           (1-epsilon_a) * ((DUb+DPb-lag[which(names(state) == "DUb")]-lag[which(names(state) == "DPb")])/N_b) ) -
        omega_a * IUa * (epsilon_a *((SPa + IPa + RPa)/N_a) + (1-epsilon_a)*((SPb + IPb + RPb)/N_b))
      
      dIPa <- SPa*trans_p*kappa*(c_aa*(IUa/N_a + kappa*IPa/N_a) + c_ba*(IUb/N_b + kappa*IPb/N_b)) - 
        phi * IPa -
        (rho) * IPa +
        theta_a * IUa * (epsilon_a * ((DUa+DPa-lag[which(names(state) == "DUa")]-lag[which(names(state) == "DPa")])/N_a) + 
                           (1-epsilon_a) * ((DUb+DPb-lag[which(names(state) == "DUb")]-lag[which(names(state) == "DPb")])/N_b) ) +
        omega_a * IUa * (epsilon_a *((SPa + IPa + RPa)/N_a) + (1 - epsilon_a)*((SPb + IPb + RPb)/N_b))
      
      dRUa <- rho * (1-mu) * IUa  + phi * RPa  -
        theta_a * RUa * (epsilon_a * ((DUa+DPa-lag[which(names(state) == "DUa")]-lag[which(names(state) == "DPa")])/N_a) + 
                           (1-epsilon_a) * ((DUb+DPb-lag[which(names(state) == "DUb")]-lag[which(names(state) == "DPb")])/N_b) ) -
        omega_a * RUa * (epsilon_a *((SPa + IPa + RPa)/N_a) + (1-epsilon_a)*((SPb + IPb + RPb)/N_b))
      
      dRPa <- rho * (1-mu) * IPa  - phi * RPa  +
        theta_a * RUa * (epsilon_a * ((DUa+DPa-lag[which(names(state) == "DUa")]-lag[which(names(state) == "DPa")])/N_a) + 
                           (1-epsilon_a) * ((DUb+DPb-lag[which(names(state) == "DUb")]-lag[which(names(state) == "DPb")])/N_b) ) +
        omega_a * RUa * (epsilon_a *((SPa + IPa + RPa)/N_a) + (1-epsilon_a)*((SPb + IPb + RPb)/N_b))
      
      dDUa <- rho*mu*IUa
      
      dDPa <- rho*mu*IPa
      
      
      #transitions for group b
      
      dSUb <- -SUb*trans_p*(c_ab*(IUa/N_a + kappa*IPa/N_a) + c_bb*(IUb/N_b + kappa*IPb/N_b))  + phi * SPb - 
        theta_b * SUb * ((1-epsilon_b) * ((DUa+DPa-lag[which(names(state) == "DUa")]-lag[which(names(state) == "DPa")])/N_a) + 
                           (epsilon_b) * ((DUb+DPb-lag[which(names(state) == "DUb")]-lag[which(names(state) == "DPb")]) / N_b)) -
        omega_b * SUb * ((1-epsilon_b) *((SPa + IPa + RPa)/N_a) + (epsilon_b)*((SPb + IPb + RPb)/N_b))
      
      dSPb <- -trans_p*kappa*SPb*(c_ab*(IUa/N_a + kappa*IPa/N_a) + c_bb*(IUb/N_b + kappa*IPb/N_b))  - phi * SPb +
        theta_b * SUb * ((1-epsilon_b) * ((DUa+DPa-lag[which(names(state) == "DUa")]-lag[which(names(state) == "DPa")])/N_a) + 
                           (epsilon_b) * ((DUb+DPb-lag[which(names(state) == "DUb")]-lag[which(names(state) == "DPb")]) / N_b)) +
        omega_b * SUb * ((1-epsilon_b) *((SPa + IPa + RPa)/N_a) + (epsilon_b)*((SPb + IPb + RPb)/N_b))
      
      dIUb <- trans_p*SUb*(c_ab*(IUa/N_a + kappa*IPa/N_a) + c_bb*(IUb/N_b + kappa*IPb/N_b)) + phi * IPb -
        (rho) * IUb -
        theta_b * IUb * ((1-epsilon_b) * ((DUa+DPa-lag[which(names(state) == "DUa")]-lag[which(names(state) == "DPa")])/N_a) + 
                           (epsilon_b) * ((DUb+DPb-lag[which(names(state) == "DUb")]-lag[which(names(state) == "DPb")]) / N_b)) -
        omega_b * IUb * ((1-epsilon_b) *((SPa + IPa + RPa)/N_a) + (epsilon_b)*((SPb + IPb + RPb)/N_b))
      
      dIPb <- trans_p*kappa*SPb*(c_ab*(IUa/N_a + kappa*IPa/N_a) + c_bb*(IUb/N_b + kappa*IPb/N_b)) - phi * IPb -
        (rho) * IPb +
        theta_b * IUb * ((1-epsilon_b) * ((DUa+DPa-lag[which(names(state) == "DUa")]-lag[which(names(state) == "DPa")])/N_a) + 
                           (epsilon_b) * ((DUb+DPb-lag[which(names(state) == "DUb")]-lag[which(names(state) == "DPb")]) / N_b)) +
        omega_b * IUb * ((1-epsilon_b) *((SPa + IPa + RPa)/N_a) + (epsilon_b)*((SPb + IPb + RPb)/N_b))
      
      dRUb <- (1-mu) * rho * IUb + phi * RPb -
        theta_b * RUb * ((1-epsilon_b) * ((DUa+DPa-lag[which(names(state) == "DUa")]-lag[which(names(state) == "DPa")])/N_a) + 
                           (epsilon_b) * ((DUb+DPb-lag[which(names(state) == "DUb")]-lag[which(names(state) == "DPb")]) / N_b)) -
        omega_b * RUb * ((1-epsilon_b) *((SPa + IPa + RPa)/N_a) + (epsilon_b)*((SPb + IPb + RPb)/N_b))
      
      dRPb <- (1-mu) * rho * IPb  - phi * RPb +
        theta_b * RUb * ((1-epsilon_b) * ((DUa+DPa-lag[which(names(state) == "DUa")]-lag[which(names(state) == "DPa")])/N_a) + 
                           (epsilon_b) * ((DUb+DPb-lag[which(names(state) == "DUb")]-lag[which(names(state) == "DPb")]) / N_b)) +
        omega_b * RUb * ((1-epsilon_b) *((SPa + IPa + RPa)/N_a) + (epsilon_b)*((SPb + IPb + RPb)/N_b))
      
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
#Average daily total contacts in BICS: Democrats = 5.618922; Republicans = 7.188501

N0 = 10000000; frac_a = 0.5; h_a = 0.5
N_a = frac_a*N0; N_b = N0 - N_a

ls<-sir_two_group_pu(c = NA, c_a = 7.1, c_b = 5.6, 
                     trans_p = 0.05, rho=1/10, mu = 0.01, 
                     h_a=h_a, kappa=0.3, phi = 0.02,
                     I0_a=1, I0_b=1, N0 = N0, frac_a = frac_a, time = 200,
                     theta_a = 100, theta_b = 200, epsilon = 0.5,
                     omega_a = 0.1, omega_b = 0.2,
                     get_params=TRUE)
sim<-ls$sim


plot(SUa/N_a~time, data=sim, type="l", col="red", ylab="Susceptible", main="a", ylim = c(0,1))
lines(SPa/N_a~time, data=sim, col="red", lty=2)
lines(SUb/N_b~time, data=sim, col="blue", type="l", main="b", ylab="Susceptible", ylim = c(0,1))
lines(SPb/N_b~time, data=sim, col="blue", lty=2)


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
c_a = 7.1; c_b = 5.6 # average number of contacts per day 
trans_p = 0.05 # probability of transmission given contact
rho = 1/10 # 1 / infectious period
mu = 0.01 # probability of dying following infection
kappa = 0.3 # reduction in probability of transmission given contact resulting from protective behavior
phi = 0 # waning rate of protective behavior


time = 365 # time steps for simulation
theta_a = 100 # responsiveness to deaths for adopting protective behavior in group A
theta_b = 200 # responsiveness to deaths for adopting protective behavior in group B 
omega_a = 0.1 # responsiveness to proportion of protected individuals for adopting protective behavior in group A
omega_b = 0.2 # responsiveness to proportion of protected individuals for adopting protective behavior in group A

h_a = c(.5, .8) # proportion of group A's total contact with members of their own group (bounded by population size and total contacts in group B)
epsilon = c(.5, .99) #measure of assortativeness in influence (i.e. are people equally aware and influenced by deaths/number of protected individuals in their own group versus in the out group)


expand.grid(h_a=h_a, epsilon=epsilon)  %>%
  data.frame() %>%
  apply(1, function(par) sir_two_group_pu(c = NA, c_a = c_a, c_b = c_b, 
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


ggplot(df %>% filter(epsilon == 0.99)) +
  geom_line(aes(x=time, y=(IUa+IPa)/N_a), size=line_sz, color="red") + 
  geom_line(aes(x=time, y=(IUb+IPb)/N_b), size=line_sz, color="blue") + 
  facet_grid(~h_a) +
  ylab("Infections")+
  theme_minimal()
  

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


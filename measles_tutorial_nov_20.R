#==================================================================================================
# This is a measles tutorial that uses the library for partially observed Markov processes. All of
# this code was written by Professor Aaron King at the University of Michigan. I copied all of this
# code from his tutorial found at this website: https://kingaa.github.io/sbied/stochsim/stochsim.html
# You can find additional information about what is going on in the code on his website. I have
# tried to annotate this code to the best of my abilities but please reach out if there are
# problems at justinksheen@gmail.com. Thank you. -Justin Sheen November 2019
#==================================================================================================

#==================================================================================================
# - Read in libraries
# - Read in the raw Measles data
# - Take a look at what the raw Measles data looks like
#==================================================================================================

library(dplyr)
library(tidyverse)
library(pomp)
library(foreach)
library(doParallel)
library(doRNG)
read_csv("https://kingaa.github.io/sbied/stochsim/Measles_Consett_1948.csv") %>%
  select(week,reports=cases) -> meas
as.data.frame(meas)

#==================================================================================================
# PROCESS MODEL
# - Create a function sir_step that tells how the states of the model change at each iteration of
#   the process.
# - Create initial conditions of the process.
# - Load everything into a "pomp" object (variable) named "measSIR"
#==================================================================================================
sir_step <- function (S, I, R, H, N, Beta, mu_IR, delta.t, ...) {
  dN_SI <- rbinom(n=1,size=S,prob=1-exp(-Beta*I/N*delta.t))
  dN_IR <- rbinom(n=1,size=I,prob=1-exp(-mu_IR*delta.t))
  S <- S - dN_SI
  I <- I + dN_SI - dN_IR
  R <- R + dN_IR
  H <- H + dN_IR;
  c(S = S, I = I, R = R, H = H)
}

sir_init <- function (N, eta, ...) {
  c(S = round(N*eta), I = 1, R = round(N*(1-eta)), H = 0)
}

meas %>% 
  pomp(
    times="week",t0=0,
    rprocess=euler(sir_step,delta.t=1/7),
    rinit=sir_init,accumvars="H"
  ) -> measSIR

#==================================================================================================
# OBSERVATION MODEL
# - Create a function that will evaluate how likely the observed number of reports is given that
#   this number of reports is assumed to have originated from a binomial distribution with size=H
#   and probability=rho.
# - Create a function that will model the number of reports received of Measles as a random draw
#   from a binomial distribution with size=H and probability=rho
# - Load everything into the "measSIR" pomp object
#==================================================================================================
dmeas <- function (reports, H, rho, log, ...) {
  dbinom(x=reports, size=H, prob=rho, log=log)
}

rmeas <- function (H, rho, ...) {
  c(reports=rbinom(n=1, size=H, prob=rho))
}

measSIR %>% pomp(rmeasure=rmeas,dmeasure=dmeas) -> measSIR

#==================================================================================================
# - Rewrite everything we just did above in C so that everything runs faster
#==================================================================================================
sir_step <- Csnippet("
  double dN_SI = rbinom(S,1-exp(-Beta*I/N*dt));
                     double dN_IR = rbinom(I,1-exp(-mu_IR*dt));
                     S -= dN_SI;
                     I += dN_SI - dN_IR;
                     R += dN_IR;
                     H += dN_IR;
                     ")

sir_init <- Csnippet("
                     S = nearbyint(eta*N);
                     I = 1;
                     R = nearbyint((1-eta)*N);
                     H = 0;
                     ")

dmeas <- Csnippet("
                  lik = dbinom(reports,H,rho,give_log);
                  ")

rmeas <- Csnippet("
                  reports = rbinom(H,rho);
                  ")

measSIR %>%
  pomp(rprocess=euler(sir_step,delta.t=1/7),
       rinit=sir_init,
       rmeasure=rmeas,
       dmeasure=dmeas,
       accumvars="H",
       statenames=c("S","I","R","H"),
       paramnames=c("Beta","mu_IR","N","eta","rho"),
       params=c(Beta=15,mu_IR=0.5,rho=0.5,eta=0.06,N=38000) # You should not worry about this line.
                                                            # This is an initial guess, but pomp does
                                                            # not use this explicitly in the below code
                                                            # for simulation nor particle filtering.
  ) -> measSIR

#==================================================================================================
# SIMULATION
# - Here, we will simulate the model with some parameters. They are intentionally bad as we will see
#   when we plot the simulations, and the idea is that a student should play around with these 
#   parameter values until they arrive at simulations that look "reasonably" like the data.
#==================================================================================================
measSIR %>%
  simulate(params=c(Beta=7.5, # Value of beta
                    mu_IR=0.5, # Value of transition from I to R box
                    rho=0.5, # Recovery value
                    eta=0.03, # Initial number of infected
                    N=38000), # Number of persons
           nsim=20, # Number of simulations to run
           format="data.frame",include.data=TRUE) -> sims

# Everything will be saved in the "sims" data frame variable.

#==================================================================================================
# PLOT SIMULATIONS
# - Here, we will plot the simulations in red, and the real measles data in blue.
#==================================================================================================
sims %>%
  ggplot(aes(x=week,y=reports,group=.id,color=.id=="data"))+
  geom_line()+
  guides(color=FALSE)


#==================================================================================================
# PARTICLE FILTER
# - Create a "sliceDesign" that will vary two parameters: Beta and mu_IR. You can look at what
#   sliceDesign looks like by simply typing in "p" in the console.
# - Prepare for parallel computing, then run everything in a big "foreach" loop. Save the likelihood#
#   of each parameter value in "p"
# - Plot results. With Np=5000 particles, it is incredibly good. With Np=500 particles less so. 
#==================================================================================================
sliceDesign(
  center=coef(measSIR),
  Beta=rep(seq(from=5,to=20,length=40),each=3), # What this line is saying, is that it wants to 
                                                # create a sequence of length 40 of real numbers 
                                                # from 5 to 20. It wants to repeat ("rep") this 
                                                # sequence 3 times, hench the "each=3" argument.
  
  mu_IR=rep(seq(from=0.2,to=2,length=40),each=3) # What this line is saying, is that it wants to 
                                                 # create a sequence of length 40 of real numbers 
                                                 # from 0.2 to 2. It wants to repeat ("rep") this 
                                                 # sequence 3 times, hench the "each=3" argument.
) -> p


registerDoParallel()
registerDoRNG(108028909)

foreach (theta=iter(p,"row"),
         .combine=rbind,.inorder=FALSE) %dopar% {
           library(pomp)
           
           measSIR %>% pfilter(params=theta,Np=5000) -> pf
           
           theta$loglik <- logLik(pf)
           theta
         } -> p

p %>% 
  gather(variable,value,Beta,mu_IR) %>%
  filter(variable==slice) %>%
  ggplot(aes(x=value,y=loglik,color=variable))+
  geom_point()+
  facet_grid(~variable,scales="free_x")+
  guides(color=FALSE)+
  labs(x="parameter value",color="")+
  theme_bw()


# Started: September 25, 2019
# Justin K. Sheen and Michael Z. Levy
#	Chagas house infestation removal POMP

#===================================================================================================
#	Load in Libraries and Data
#    -load in Libraries
#    -bring in fake data for pomp object
#===================================================================================================

library(dplyr)
library(readr)
library(pomp)
library(ggplot2)
library(tidyverse)
library(foreach)
library(doParallel)
library(doRNG)
library(boot)


#chagas <- read.csv("~/PETM-shiny/pomp/real_infestation_data.csv")
chagas <- as.data.frame(1:600)
chagas <- cbind(chagas, chagas, chagas)
colnames(chagas) <- c("month", "TRUE", "reports")

#===================================================================================================
#	Process model
#    -sir_step is a stochastic SIHR model
#     -S suceptible houses
#     -I infested and infectious houses
#     -H a counter of how many houses were detected and treated at each time sir_step
#     -R Detected and Treated
# Parameters
#      -Beta transmission
#      -mu rate of detection AND treatment
#      -gamma rate of wearing off of insecticide
#      -delta.t is for rescaling time steps. Not used (set to 1) presently
# Initial Conditions
#      -eta the unobserved Prevalence of Infested houses at time=0
#===================================================================================================

sir_step <- function (S, I, R, H, N, Beta, mu, gamma, delta.t, ...) {
  dN_SI <- rbinom(n=1, size=S, prob=1-exp(-inv.logit(Beta)*(I)*delta.t))
  dN_IR <- rbinom(n=1, size=I, prob=inv.logit(mu) * delta.t)
  dN_RS <- rbinom(n=1, size=R, prob=inv.logit(gamma) * delta.t)
  S <- S + dN_RS - dN_SI
  I <- I + dN_SI - dN_IR
  H <- dN_IR
  R <- R + H - dN_RS
  c(S = S, I = I, R = R, H = H)
}

sir_init <- function (N, eta, ...) {
  c(S = round(N*(1-inv.logit(eta))), I = round(N*inv.logit(eta)), R = 0, H = 0)
}

#===================================================================================================
#	Observation Model
#    -rmeas assigns H as the number of observed reports at the time step
#    -dmeas probability that the observed number of reports (H) originated from a binomial distribution with expected value I*mu
#       -uses a binomial distribution
#===================================================================================================
rmeas <- function (H, ...) {
  c(reports=H)
}

dmeas <- function (reports, I, mu, log, ...) {
  dbinom(x=reports, size=I, prob=inv.logit(mu), log=log)
}

#===================================================================================================
#	Create a pomp object Chagas SIR
#     -Define Process model
#     -Set initial conditions
#     -Define Observation model
#===================================================================================================
chagas %>%
  pomp(times="month",t0=0,
       rprocess=euler(sir_step, delta.t=1),
       rinit=sir_init,
       dmeasure=dmeas,
       rmeasure=rmeas,
       statenames=c("S", "I", "R", "H"),
       paramnames=c("Beta", "mu", "gamma", "eta", "N"),
       params=c(Beta=1,mu=1,gamma=1,eta=1,N=1) # Need to input this argument for some reason
  ) -> chagasSIR

#===================================================================================================
#	Csnippets of the above
#     -sir_step Process Model  (note dt is the time scale delta.t above)
#     -sir_init Initial conditions
#     -rmeas sets observed values to H
#     -dmeas gets probability of observed number, H, from binomial with mean I*mu
#===================================================================================================

sir_step <- Csnippet("
                     double dN_SI = rbinom(S, 1 - exp(- (pow(2.71828, Beta) / (pow(2.71828, Beta) + 1)) * I * dt));
                     double dN_IR = rbinom(I, (pow(2.71828, mu) / (pow(2.71828, mu) + 1)) * dt);
                     double dN_RS = rbinom(R, gamma * dt);
                     S += dN_RS - dN_SI;
                     I += dN_SI - dN_IR;
                     H = dN_IR;
                     R += H - dN_RS;
                     ")
sir_init <- Csnippet("
                     S = nearbyint(N * (1 - (pow(2.71828, eta) / (pow(2.71828, eta) + 1)) ));
                     I = nearbyint(N * (pow(2.71828, eta) / (pow(2.71828, eta) + 1)) );
                     H = 0;
                     R = 0;
                     ")

rmeas <- Csnippet("
                 reports = H;
                 ")

dmeas <- Csnippet("
                  lik = dbinom(reports, I, pow(2.71828, mu) / (pow(2.71828, mu) + 1), give_log);
                  ")

chagas %>%
  pomp(times="month",t0=0,
       rprocess=euler(sir_step, delta.t=1),
       rinit=sir_init,
       dmeasure=dmeas,
       rmeasure=rmeas,
       statenames=c("S", "I", "R", "H"),
       paramnames=c("Beta", "mu", "gamma", "eta", "N"),
       params=c(Beta=1,mu=1,gamma=1,eta=1,N=1) # Need to input this argument for some reason
  ) -> chagasSIR

#===================================================================================================
# Run a for loop of parameters in the chagas SIR model and save to a data frame named "sims." Write
# each .csv of the simulations, as well as the corresponding plot for easy understanding.
#   - Beta = transmission rate
#   - mu = rate of detection and treatment
#   - gamma = rate of wearing off of insecticide
#   - eta = starting prevalence of unobserved infestations
#   - N = number of houses in the system
# 
# Plot the results of the parameters
#   - Loop through each unique ".id" of the data frame "sims"
#   - Print the iteration number of the loop
#   - Check that the unique ".id" does not equal "data." (i.e. if it is a simulation, the ".id" should
#      be an integer value)
#   - Within the loop, subset to one set of simulations (e.g. the simulation ".id" == 2)
#   - Within the loop, if this is the first iteration, use the command "plot" to plot the results. 
#     Else, if this is not the first iteration, use the command "lines" to plot the results.
#       - Plot the "I" box at each time step in red
#       - Plot the "D" box at each time step in black
#===================================================================================================

betas <- c(1.95e-06, 3.89e-06, 5.835e-06)
mus <- c(1/6, 1/12, 1/24)
etas <- c(50/192529, 150/192529, 200/192529)

for (beta in betas) {
  for (mu in mus) {
    for (eta in etas) {
      chagasSIR %>%
        simulate(params=c(Beta=logit(beta),
                          mu=logit(mu),
                          gamma=1/6,
                          eta=logit(eta),
                          N=192529),
                 nsim=3, format="data.frame", include.data = TRUE) -> sims
      
      jpeg(paste0("~/PETM-shiny/pomp/recovery/simulated_real/sim_eta_", round(eta, digits = 4), "_beta_", round(beta, digits = 10), "_mu_", round(mu, digits = 4), ".jpg"))
      # For loop for each unique ".id" of the data frame "sims" 
      for (i in 1:length(unique(sims$.id))) {
        # Check that the unique ".id" does not equal "data." (i.e. if it is a simulation, the ".id" should
        # be an integer value)
        if (as.character(unique(sims$.id)[i] != "data")) {
          # Within the loop, subset to one set of simulations (e.g. the simulation ".id" == 2)
          temp <- sims[which(sims$.id == unique(sims$.id)[i]),]
          # If this is the first iteration, "plot" the results. Else, use the command "lines"
          # Plot the "I" box at each time step in red
          # Plot the "D" box at each time step in black
          if (i == 2) {
            plot(temp$month, temp$I,pch=".", ylim=c(0,max(temp$I, temp$R)),
                 main="Stochastic model", ylab="# of houses", xlab="month")
            lines(temp$month, temp$R, pch=".")
            lines(temp$month, temp$reports, col="blue")
          } else {
            lines(temp$month, temp$I, pch=".", col="red")
            lines(temp$month, temp$R, pch=".")
            lines(temp$month, temp$reports, col="blue")
          }
        }
      }
      dev.off()
      write.csv(sims, paste0("~/PETM-shiny/pomp/recovery/simulated_real/sim_eta_", round(eta, digits = 4), "_beta_", round(beta, digits = 10), "_mu_", round(mu, digits = 4), ".csv"), row.names=F)
    }
  }
}

#===================================================================================================
# "Particle filtering"
#  - Create a list of parameter combinations you would like to get the likelihoods for, and assign
#    it to a data frame "p"
#  - We plan to use the "foreach" function in order to use parallel computing, and to do so, we will
#    use the registerDoParallel() and registerDoRNG() commands to register with the backend of the
#    computer
#  - Using the "foreach" function, we assign iterate over each row of "p," assigning the parameter
#    values to the variable "theta." We then use a particle filter for this parameter set of size 
#    Np = 500, and save this likelihood value in the variable pf. We then save the "logLik(pf)" into
#    the "loglik" column of theta. The end result is that the data frame "p" now has a new column of
#    the log likelihood values of each parameter set.
#===================================================================================================
betas <- c(1.95e-07, 3.89e-07)
mus <- c(1/6, 1/12, 1/24)
etas <- c(50/192529, 150/192529)

#betas <- c(1.95e-07)
#mus <- c(1/6)
#etas <- c(50/192529)
for (beta in betas) {
  for (mu in mus) {
    for (eta in etas) {
      print(paste0("Beta: ", beta))
      print(paste0("mu: ", mu))
      print(paste0("eta: ", eta))
      
      #===========================================================================================
      # Read in the .csv that of the simulated infestation. Take the first simulation (there are three for
      # each dataset) and use this as the new infestation data for pomp. Also, we only need the "month" and
      # "reports" columns. Prepare the same pomp object as before.
      #===========================================================================================
      simulated_infestation <- read.csv(paste0("~/PETM-shiny/pomp/recovery/simulated_real/sim_eta_", round(eta, digits = 4), "_beta_", round(beta, digits = 10), "_mu_", round(mu, digits = 4), ".csv"))
      simulated_infestation <- simulated_infestation[which(simulated_infestation$.id == 1),]
      simulated_infestation <- as.data.frame(cbind(simulated_infestation$month, simulated_infestation$reports))
      colnames(simulated_infestation) <- c("month", "reports")
      
      simulated_infestation %>%
        pomp(times="month",t0=0,
             rprocess=euler(sir_step, delta.t=1),
             rinit=sir_init,
             dmeasure=dmeas,
             rmeasure=rmeas,
             statenames=c("S", "I", "R", "H"),
             paramnames=c("Beta", "mu", "gamma", "eta", "N")
        ) -> chagasSIR
      
      #===========================================================================================
      # - Create a list of parameter combinations you would like to get the likelihoods for, and assign
      #   it to a data frame "p"
      # - Ready the computer backend to do parallel computing using registerDoParallel() and registerDoRNG
      # - Create grid of points to scan over
      # - For each parameter combination, get the log likelihood using the particle filter, and assign the
      #   log likelihood value to the column "loglik" of the data frame "p"
      #===========================================================================================
      registerDoParallel()
      registerDoRNG(421776444)
      
      expand.grid(
        center=1,
        Beta=seq(from=0, to=5.835e-06,length=50),
        mu=seq(from=0,to=0.5,length=50),
        eta=seq(from=0,to=0.005,length=50),
        gamma=1/6,
        N=192529
      ) -> p
      
      foreach (theta=iter(p,"row"),
               .combine=rbind,.inorder=FALSE) %dopar% {
                 library(pomp)
                 
                 chagasSIR %>% pfilter(params=theta,Np=5000) -> pf
                 
                 theta$loglik <- logLik(pf)
                 theta
               } -> p
      
      write.csv(p, paste0("~/PETM-shiny/pomp/recovery/log_likelihoods/log_lik_eta_", round(eta, digits = 4),"_beta_", round(beta, digits = 10), "_mu_", round(mu, digits = 4), ".csv"))
    }
  }
}

#===================================================================================================
# Create nice pdfs of the results:
# - Each row will be a different value of beta
# - Each column will be a different value of mu
#===================================================================================================
betas <- c(1.95e-07, 3.89e-07)
mus <- c(1/6, 1/12, 1/24)
etas <- c(50/192529, 150/192529)
for (eta in etas) {
  for (beta in betas) {
    for (mu in mus) {
      pdf(file=paste0("~/PETM-shiny/pomp/recovery/log_lik_plots/eta_", round(eta, digits=4), "_beta_", round(beta, digits=10), "_mu_", round(mu, digits=4), ".pdf"),width=20, height=15)
      p <- read.csv(paste0("~/PETM-shiny/pomp/recovery/log_likelihoods/log_lik_eta_", round(eta, digits=4), "_beta_", round(beta, digits = 10), "_mu_", round(mu, digits = 4), ".csv"))
      plots <- list()
      i <- 1
      for (recov_eta in unique(p$eta)) {
        temp_p <- p[which(p$eta == recov_eta),]
        temp_p %>%
          mutate(loglik=ifelse(loglik>max(loglik)-50,loglik,NA)) -> temp_p
        plotted <-
          ggplot(aes(x=Beta,y=mu,z=loglik,fill=loglik),data=temp_p)+
          geom_tile(color=NA)+
          scale_fill_gradient()+
          annotate("point", x = beta, y = mu, colour = "red", size=3)+
          labs(x=expression(beta),y=expression(mu[ID]), title = paste0("eta_", round(recov_eta, digits=4)))
        plots[[i]] <- plotted
        i <- i + 1
      }
      
      grid.arrange(plots[[1]],plots[[2]],plots[[3]],
                   plots[[4]],plots[[5]],plots[[6]],
                   plots[[7]],plots[[8]],plots[[9]],
                   plots[[10]],plots[[11]],plots[[12]],
                   plots[[13]],plots[[14]],plots[[15]],
                   plots[[16]],plots[[17]],plots[[18]],
                   plots[[19]],plots[[20]],plots[[21]],
                   plots[[22]],plots[[23]],plots[[24]],
                   plots[[25]],plots[[26]],plots[[27]],
                   plots[[28]],plots[[29]],plots[[30]],
                   plots[[31]],plots[[32]],
                   plots[[33]],plots[[34]],plots[[35]],
                   plots[[36]],plots[[37]],plots[[38]],
                   plots[[39]],plots[[40]],plots[[41]],
                   plots[[42]],plots[[43]],plots[[44]],
                   plots[[45]],plots[[46]],plots[[47]],
                   plots[[48]],plots[[49]],plots[[50]],
                   nrow=10, ncol=5)
      dev.off()
    }
  }
}

#===================================================================================================
# Find 95% density of parameters
# - For each "real" infestation, do the following:
# - First, calculate the total density of the grid (exponentiate all log likelihoods)
# - Then, calculate the percent of the density each parameter combination has
# - Then, sort them by the percent of the density, and iterate until 95% of the density has
#   been recaptured
#===================================================================================================
betas <- c(1.95e-07, 3.89e-07)
mus <- c(1/6, 1/12, 1/24)
etas <- c(50/192529, 150/192529)
for (eta in etas) {
  for (beta in betas) {
    for (mu in mus) {
      p <- read.csv(paste0("~/PETM-shiny/pomp/recovery/log_likelihoods/log_lik_eta_", round(eta, digits=4), "_beta_", round(beta, digits = 10), "_mu_", round(mu, digits = 4), ".csv"))
      p <- p[-which(p$mu == 0 | p$eta == 0),]
      p$lik <- exp(p$loglik)
      p$lik_dens <- p$lik / sum(p$lik)
      for (i in 1:nrow(p)) {
        if (p$lik_dens[i] > 1) {
          print("Error in normalized density.")
        }
      }
      p <- p[order(-p$lik_dens),]
      tot <- 0
      for (row_num in 1:nrow(p)) {
        tot <- tot + p$lik_dens[row_num]
        if (tot > 0.95) {
          break
        }
      }
      p <- p[1:row_num,]
      # Plot distribution
      pdf(paste0("~/PETM-shiny/pomp/recovery/ninety_five_dens/eta_", round(eta, digits=4), "_beta_", round(beta, digits = 10), "_mu_", round(mu, digits = 4), '.pdf'), width=9, height=3)
      par(mfrow=c(1,3))
      plot(p$Beta, p$lik_dens, main="Beta")
      plot(p$mu, p$lik_dens, main="mu")
      plot(p$eta, p$lik_dens, main="eta")
      dev.off()
      
      min_beta <- min(p$Beta)
      min_eta <- min(p$eta)
      min_mu <- min(p$mu)
      max_beta <- max(p$Beta)
      max_eta <- max(p$eta)
      max_mu <- max(p$mu)
      # For each of the min and max parameters, go through the sequence and see which block it is
      # a part of
      # Beta
      nf_min_beta <- NA
      nf_max_beta <- NA
      bs <- seq(from=0, to=5.835e-06,length=50)
      for (b_i in 1:length(bs)) {
        if (all.equal(bs[b_i], min_beta) == TRUE) {
          nf_min_beta <- bs[b_i]
        }
        if (all.equal(bs[b_i], max_beta) == TRUE) {
          nf_max_beta <- bs[b_i]
        }
      }
      # Eta
      nf_min_eta <- NA
      nf_max_eta <- NA
      es <- seq(from=0,to=0.005,length=50)
      es <- es[2:50]
      for (e_i in 1:length(es)) {
        if (all.equal(es[e_i], min_eta) == TRUE) {
          nf_min_eta <- es[e_i]
        }
        if (all.equal(es[e_i], max_eta) == TRUE) {
          nf_max_eta <- es[e_i]
        }
      }
      # Mu
      nf_min_mu <- NA
      nf_max_mu <- NA
      ms <- seq(from=0,to=0.5,length=50)
      ms <- ms[2:50]
      for (m_i in 1:length(ms)) {
        if (all.equal(ms[m_i], min_mu) == TRUE) {
          nf_min_mu <- ms[m_i]
        }
        if (all.equal(ms[m_i], max_mu) == TRUE) {
          nf_max_mu <- ms[m_i]
        }
      }
      
      o0 <- paste0("Num. combinations: ", row_num)
      o1 <- paste0("True Beta: ", round(beta, digits=10))
      o2 <- paste0("Estimated Beta: [", round(nf_min_beta, digits = 10), ",", round(nf_max_beta, digits = 10), "]")
      o3 <- paste0("Captured Beta: ", beta >= nf_min_beta & beta <= nf_max_beta)
      o4 <- paste0("True eta: ", eta * 192529)
      o5 <- paste0("Estimated eta: [", round(nf_min_eta * 192529), ",", round(nf_max_eta * 192529), "]")
      o6 <- paste0("Captured eta: ", eta >= nf_min_eta & eta <= nf_max_eta)
      o7 <- paste0("True mu: ", round(1 / mu))
      o8 <- paste0("Estimated mu: [", round(1 / nf_max_mu), ",", round(1 / nf_min_mu), "]")
      o9 <- paste0("Captured mu: ", round(1/mu) <= round(1/nf_min_mu) & round(1/mu) >= round(1/nf_max_mu))
      res<-list(o0, o1, o2, o3, o4, o5, o6, o7, o8, o9)
      filenameRes <- file.path(paste0("~/PETM-shiny/pomp/recovery/ninety_five_dens/ninety_five_dens_eta_", round(eta, digits=4), "_beta_", round(beta, digits = 10), "_mu_", round(mu, digits = 4), ".txt"))
      capture.output(res, file = filenameRes)
    }
  }
}

#===================================================================================================
# Create pdf that will tell us how well this methodology will do
# - Plot the simulation traces (reports(t), I(t), D(t))
# - Report the parameter combinations that have 95% of density
#===================================================================================================
betas <- c(1.95e-07, 3.89e-07)
mus <- c(1/6, 1/12, 1/24)
etas <- c(50/192529, 150/192529)
pdf(file=paste0("~/PETM-shiny/pomp/recovery/success_plots/ninety_five_success.pdf"),width=30, height=40)
par(mfrow=c(4,3))
for (eta in etas) {
  for (beta in betas) {
    for (mu in mus) {
      simulated_infestation <- read.csv(paste0("~/PETM-shiny/pomp/recovery/simulated_real/sim_eta_", round(eta, digits = 4), "_beta_", round(beta, digits = 10), "_mu_", round(mu, digits = 4), ".csv"))
      simulated_infestation <- simulated_infestation[which(simulated_infestation$.id == 1),]
      simulated_infestation <- as.data.frame(cbind(simulated_infestation$month, simulated_infestation$reports, simulated_infestation$I, simulated_infestation$R))
      colnames(simulated_infestation) <- c("month", "reports", "I", "R")
      plot(simulated_infestation$month, simulated_infestation$I, col="red",pch=".", ylim=c(0,100),
           main="# of infested and detected houses over time", ylab="# of houses", xlab="month")
      lines(simulated_infestation$month, simulated_infestation$I, col="red")
      lines(simulated_infestation$month, simulated_infestation$R, col="blue")
      lines(simulated_infestation$month, simulated_infestation$reports, col="black")
      
      success <- read.table(paste0("~/PETM-shiny/pomp/recovery/ninety_five_dens/ninety_five_dens_eta_", round(eta, digits=4), "_beta_", round(beta, digits = 10), "_mu_", round(mu, digits = 4), ".txt"), sep="\n")
      # Beta
      if (success$V1[8] == "[1] Captured Beta: TRUE") {
        text(cex = 1.5, x=50, y=100, labels=paste0(success$V1[6], "\n",success$V1[4]), col="green", pos=1)
      } else {
        text(cex = 1.5, x=50, y=100, labels=paste0(success$V1[6], "\n",success$V1[4]), col="grey", pos=1)
      }
      # Mu
      if (success$V1[20] == "[1] Captured mu: TRUE") {
        text(cex = 1.5, x=50, y=80, labels=paste0(success$V1[18], "\n",success$V1[16]), col="green", pos=1)
      } else {
        text(cex = 1.5, x=50, y=80, labels=paste0(success$V1[18], "\n",success$V1[16]), col="grey", pos=1)
      }
      # Eta
      # Mu
      if (success$V1[14] == "[1] Captured eta: TRUE") {
        text(cex = 1.5, x=50, y=60, labels=paste0(success$V1[12], "\n",success$V1[10]), col="green", pos=1)
      } else {
        text(cex = 1.5, x=50, y=60, labels=paste0(success$V1[12], "\n",success$V1[10]), col="grey", pos=1)
      }
    }
  }
}
dev.off()

#===================================================================================================
# Create pdf that show the best guess (NOT the 95% density surface)
#===================================================================================================
betas <- c(1.95e-07, 3.89e-07)
mus <- c(1/6, 1/12, 1/24)
etas <- c(50/192529, 150/192529)
pdf(file=paste0("~/PETM-shiny/pomp/recovery/success_plots/best_success.pdf"),width=30, height=40)
par(mfrow=c(4,3))
for (eta in etas) {
  for (beta in betas) {
    for (mu in mus) {
      simulated_infestation <- read.csv(paste0("~/PETM-shiny/pomp/recovery/simulated_real/sim_eta_", round(eta, digits = 4), "_beta_", round(beta, digits = 10), "_mu_", round(mu, digits = 4), ".csv"))
      simulated_infestation <- simulated_infestation[which(simulated_infestation$.id == 1),]
      simulated_infestation <- as.data.frame(cbind(simulated_infestation$month, simulated_infestation$reports, simulated_infestation$I, simulated_infestation$R))
      colnames(simulated_infestation) <- c("month", "reports", "I", "R")
      plot(simulated_infestation$month, simulated_infestation$I, col="red",pch=".", ylim=c(0,100),
           main="# of infested and detected houses over time", ylab="# of houses", xlab="month")
      lines(simulated_infestation$month, simulated_infestation$I, col="red")
      lines(simulated_infestation$month, simulated_infestation$R, col="blue")
      lines(simulated_infestation$month, simulated_infestation$reports, col="black")
      
      p <- read.csv(paste0("~/PETM-shiny/pomp/recovery/log_likelihoods/log_lik_eta_", round(eta, digits=4), "_beta_", round(beta, digits = 10), "_mu_", round(mu, digits = 4), ".csv"))
      p <- p[-which(p$mu == 0 | p$eta == 0),]
      max_beta <- p$Beta[which(p$loglik == max(p$loglik))]
      max_eta <- p$eta[which(p$loglik == max(p$loglik))]
      max_mu <- p$mu[which(p$loglik == max(p$loglik))]

      # Beta
      text(cex = 1.5, x=50, y=100, labels=paste0("True beta: ", beta, "\n Best est. beta: ",max_beta, "\n diff: ", abs(max_beta - beta)), col="green", pos=1)

      # Eta
      text(cex = 1.5, x=50, y=80, labels=paste0("True eta: ", round(192529*eta), "\n Best est. eta: ",round(max_eta * 192529), "\n diff: ", round(192529*abs(max_eta - eta))), col="red", pos=1)

      # Mu
      text(cex = 1.5, x=50, y=60, labels=paste0("True mu: ", round(1/mu), "\n Best est. mu: ", round(1/(max_mu)), "\n diff: ", round(abs(1/max_mu - 1/mu))), col="blue", pos=1)
    }
  }
}
dev.off()

#===================================================================================================
# PMCMC Global search
#===================================================================================================
# Need to set up prior distributions for the parameters
hyperparams <- list(min = coef(chagasSIR), max = coef(chagasSIR))
hyperparams$min <- c(logit(1.95e-06),#logit(1e-8), 
                     logit(1/72), 
                     1/6, 
                     logit(1/192529), 
                     192529)
hyperparams$max <- c(logit(1.95e-06),#logit(0.00005), 
                     logit(1/1.1),
                     1/6, 
                     logit(10000/192529), 
                     192529)
names(hyperparams$min) <- c("Beta", "mu", "gamma", "eta", "N")
names(hyperparams$max) <- c("Beta", "mu", "gamma", "eta", "N")
chagasSIR.dprior <- function (Beta, mu, eta, log,...) {
  f_beta <- log(1) 
    #dunif(Beta, min = hyperparams$min[1], max = hyperparams$max[1],log = TRUE)
  f_mu <- dunif(mu, min = hyperparams$min[2], max = hyperparams$max[2],log = TRUE)
  f_gamma <- log(1)
  f_eta <- dpois(eta, logit(50 / 192529), log = TRUE)
    #dunif(eta, min = hyperparams$min[4], max = hyperparams$max[4],log = TRUE)
  f_N <- log(1)
  f <- sum(f_beta + f_mu + f_gamma + f_eta + f_N)
  if (log) f else exp(f)
}

# Starting values
expand.grid(
  Beta=logit(1.95e-06),#rep((logit(1e-8) + logit(0.00005)) / 2, 1),
  mu=rep((logit(1/72) + logit(1/1.1)) / 2, 2),
  eta=rep(logit(50/192529), 1),
  #eta=rep((logit(1/192529) + logit(10000/192529)) / 2, 1),
  gamma=1/6,
  N=192529
) -> p


# Real simulated infestation
betas <- c(1.95e-06)#, 3.89e-06)
mus <- c(1/6, 1/12, 1/24)
etas <- c(50/192529, 150/192529)
for (t_beta in betas) {
  for (t_mu in mus) {
    for (t_eta in etas) {
      simulated_infestation <- read.csv(paste0("~/PETM-shiny/pomp/recovery/simulated_real/sim_eta_", round(t_eta, digits = 4), "_beta_", round(t_beta, digits = 10), "_mu_", round(t_mu, digits = 4), ".csv"))
      simulated_infestation <- simulated_infestation[which(simulated_infestation$.id == 1),]
      simulated_infestation <- as.data.frame(cbind(simulated_infestation$month, simulated_infestation$reports))
      colnames(simulated_infestation) <- c("month", "reports")
      
      # PMCMC algorithm
      foreach (theta=iter(p,"row"),
               .combine=rbind,
               .inorder=FALSE) %dopar% {
                 library(pomp)
                 
                 simulated_infestation %>%
                   pomp(times="month",t0=0,
                        rprocess=euler(sir_step, delta.t=1),
                        rinit=sir_init,
                        dmeasure=dmeas,
                        rmeasure=rmeas,
                        statenames=c("S", "I", "R", "H"),
                        paramnames=c("Beta", "mu", "gamma", "eta", "N"),
                        dprior=chagasSIR.dprior,
                        params=theta
                   ) -> chagasSIR
                 
                 rw.sd <- c(Beta=0,#-logit((5.835e-06 - 5.835e-08)) / 10, # Divide by 10 for reasonable steps of the standard deviation
                            mu=logit((1/1.5 - 1/48)) / 10,
                            eta=-logit(0.0002597011) / 2,
                            gamma=0,
                            N=0)
                 
                 chagasSIR %>% pmcmc(
                   start=theta,
                   Nmcmc = 40000,
                   Np = 500,
                   proposal = mvn.diag.rw(rw.sd)
                 ) -> pmcmc
                 
                 as.data.frame(pmcmc@traces)
                 
               } -> pmcmc_res
      
      # The idea is to understand what the posteriors look like for each variable
      write.csv(pmcmc_res, paste0("~/PETM-shiny/pomp/traces_dprior_eta_", round(t_eta, digits = 4), "_beta_", round(t_beta, digits = 10), "_mu_", round(t_mu, digits = 4), ".csv"))
    }
  }
}

pdf(paste0("~/PETM-shiny/pomp/traces_full.pdf"), width=8, height=48)
par(mfrow=c(18,3))
for (t_beta in betas) {
  for (t_eta in etas) {
    for (t_mu in mus) {
      simulated_infestation <- read.csv(paste0("~/PETM-shiny/pomp/recovery/simulated_real/sim_eta_", round(t_eta, digits = 4), "_beta_", round(t_beta, digits = 10), "_mu_", round(t_mu, digits = 4), ".csv"))
      simulated_infestation <- simulated_infestation[which(simulated_infestation$.id == 1),]
      simulated_infestation <- as.data.frame(cbind(simulated_infestation$month, simulated_infestation$reports))
      colnames(simulated_infestation) <- c("month", "reports")

      # Load in results of pmcmc
      pmcmc_res <- read.csv(paste0("~/PETM-shiny/pomp/traces_dprior_eta_", round(t_eta, digits = 4), "_beta_", round(t_beta, digits = 10), "_mu_", round(t_mu, digits = 4), ".csv"))
      
      # Plot 2 posterior simulations
      best_log_lik_1 <- which(pmcmc_res$loglik[20000:40001] == max(pmcmc_res$loglik[20000:40001]))[1]
      best_log_lik_2 <- which(pmcmc_res$loglik[60000:80002] == max(pmcmc_res$loglik[60000:80002]))[1]
      best_log_lik <- NA
      add <- NA
      if (pmcmc_res$loglik[best_log_lik_1 + 20000] > pmcmc_res$loglik[best_log_lik_2 + 60000]) {
        best_log_lik <- best_log_lik_1
        add <- 20000
      } else {
        best_log_lik <- best_log_lik_2
        add <- 60000
      }
      posterior_beta <- inv.logit(pmcmc_res$Beta[best_log_lik + add])
      posterior_eta <- inv.logit(pmcmc_res$eta[best_log_lik + add])
      posterior_mu <- inv.logit(pmcmc_res$mu[best_log_lik + add])
      
      simulated_infestation %>%
        pomp(times="month",t0=0,
             rprocess=euler(sir_step, delta.t=1),
             rinit=sir_init,
             dmeasure=dmeas,
             rmeasure=rmeas,
             statenames=c("S", "I", "R", "H"),
             paramnames=c("Beta", "mu", "gamma", "eta", "N"),
             dprior=chagasSIR.dprior,
             params=theta
        ) -> chagasSIR
      
      chagasSIR %>%
        simulate(params=c(Beta=logit(posterior_beta),
                          mu=logit(posterior_mu),
                          gamma=1/6,
                          eta=logit(posterior_eta),
                          N=192529),
                 nsim=2, format="data.frame", include.data = TRUE) -> sims
      
      plot(simulated_infestation$month, simulated_infestation$reports, pch=".", main="Simulated infestation", xlab="month", ylab="# reports", ylim=c(0, max(simulated_infestation$reports, sims$reports[61:180])))
      lines(simulated_infestation$month, simulated_infestation$reports)
      
      # For loop for each unique ".id" of the data frame "sims" 
      for (i in 1:length(unique(sims$.id))) {
        # Check that the unique ".id" does not equal "data." (i.e. if it is a simulation, the ".id" should
        # be an integer value)
        if (as.character(unique(sims$.id)[i] != "data")) {
          # Within the loop, subset to one set of simulations (e.g. the simulation ".id" == 2)
          temp <- sims[which(sims$.id == unique(sims$.id)[i]),]
          # If this is the first iteration, "plot" the results. Else, use the command "lines"
          # Plot the "I" box at each time step in red
          # Plot the "D" box at each time step in black
          plot(temp$month, temp$reports, pch=".",
               main=paste0("Posterior sim ", i - 1, " reports"), ylab="# reports", xlab="month", col="grey", ylim=c(0, max(simulated_infestation$reports, sims$reports[61:180])))
          lines(temp$month, temp$reports, col="grey")
        }
      }
      
      plot(1:40001, inv.logit(pmcmc_res$Beta[1:40001]), pch=".", xlab="i", ylab="Beta", main="Beta PMCMC", ylim=c(0, max(inv.logit(pmcmc_res$Beta))))
      lines(1:40001, inv.logit(pmcmc_res$Beta[1:40001]))
      lines(1:40001, inv.logit(pmcmc_res$Beta[40002:80002]))
      abline(h = t_beta, col="red")
      
      plot(1:40001, inv.logit(pmcmc_res$eta[1:40001]), pch=".", xlab="i", ylab="eta", main="eta PMCMC", ylim=c(0, max(inv.logit(pmcmc_res$eta))))
      lines(1:40001, inv.logit(pmcmc_res$eta[1:40001]))
      lines(1:40001, inv.logit(pmcmc_res$eta[40002:80002]))
      abline(h = t_eta, col="red")
      
      plot(1:40001, inv.logit(pmcmc_res$mu[1:40001]), pch=".", xlab="i", ylab="mu", main="mu PMCMC", ylim=c(0, max(inv.logit(pmcmc_res$mu))))
      lines(1:40001, inv.logit(pmcmc_res$mu[1:40001]))
      lines(1:40001, inv.logit(pmcmc_res$mu[40002:80002]))
      abline(h = t_mu, col="red")
      
      hist(c(inv.logit(pmcmc_res$Beta[20000:40001]), inv.logit(pmcmc_res$Beta[60000:80002])), main="Beta posterior (i=20000:40000)", xlab="Est. Beta")
      abline(v = t_beta, col="red")
      
      hist(c(inv.logit(pmcmc_res$eta[20000:40001]), inv.logit(pmcmc_res$eta[60000:80002])), main="eta posterior (i=20000:40000)", xlab="Est. eta")
      abline(v = t_eta, col="red")
      
      hist(c(inv.logit(pmcmc_res$mu[20000:40001]), inv.logit(pmcmc_res$mu[60000:80002])), main="mu posterior (i=20000:40000)", xlab="Est. mu")
      abline(v = t_mu, col="red")
    }
  }
}
dev.off()

#===================================================================================================
# (INCOMPLETE) Testing out random data and seeing how well the particle filtering does
#===================================================================================================
data_random <- round(runif(34, 0, 8))
data_random <- as.data.frame(data_random)
data_random <- cbind(c(1:34), data_random, data_random)
colnames(data_random) <- c("day", "true", "reports")
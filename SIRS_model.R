library(deSolve)

sir <- function(time, state, parameters) 
	{
	with(
		 as.list(c(state, parameters)), 
			{
			dS <- -beta * S * I
			dI <- beta * S * I - gamma * I
			dR <- gamma * I
			return(list(c(dS, dI, dR)))
		    }
		 )
}


sis <- function(time, state, parameters) 
{
  with(
    as.list(c(state, parameters)), 
    {
      dS <- -beta * S * I + gamma * I
      dI <- beta * S * I - gamma * I
      return(list(c(dS, dI)))
    }
  )
}

sirs <- function(time, state, parameters) 
{
  with(
    as.list(c(state, parameters)), 
    {
      dS <- -beta * S * I + psi * R
      dI <- beta * S * I - gamma * I
      dR <- gamma * I - psi * R
      return(list(c(dS, dI, dR)))
    }
  )
}

init <- c(S = round(192529*(1-50/192529)), I =round(192529*50/192529), R = 0)
parameters <- c(beta = 5.835e-7, gamma = 1/12, psi=1/(100/30))
times <- seq(0, 60, by = 1)

out <- as.data.frame(ode(y = init, times = times, func = sirs, parms = parameters))
S<-out$S
I<-out$I
R<-out$R
RESULTS<-data.frame(out$S,out$I,out$R)

plot(times,I, ylim=c(0,max(I, R)),xlim=c(0,max(times)),  xlab = "Time (months)", ylab = "# houses", main = "Compartmental Model", typ="l",col="red")
lines(times,R,col="black")


# SIDR
sidr <- function(time, state, parameters) 
{
  with(
    as.list(c(state, parameters)), 
    {
      dS <- psi * R - beta * S * (I + D) 
      dI <- beta * S * (I + D) - mu * I
      dD <- mu * I - phi * D
      dR <- phi * D - psi * R
      return(list(c(dS, dI, dD, dR)))
    }
  )
}

init <- c(S = 192335, I = 193, D = 10, R = 100)
parameters <- c(beta = .0000000035, mu = 0.0005, phi = 0.00695, psi=0.00695)
times <- seq(0, 1000000, by = 1)

out <- as.data.frame(ode(y = init, times = times, func = sidr, parms = parameters))
S<-out$S
I<-out$I
D<-out$D
R<-out$R
RESULTS<-data.frame(out$S,out$I,out$R,out$D)

plot(times,I, ylim=c(0,200000),xlim=c(0,max(times)),  xlab = "Time", ylab = "I", main = "SIDR Model", typ="l",col="green")
lines(times, D)
D[length(D)]

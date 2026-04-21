# Deterministic SIR
#still has to be checked. 

SIR.ode <- function(t, state, par) {
  
  with(as.list(c(state, par)), {
    
    dS <- -beta * S * I
    
    dI <-  beta * S * I - gamma * I
    
    list(c(dS, dI))
    
  })
  
}



SIR.1LD.ode <- function(t, state, par) {
  
  with(as.list(c(state, par)), {
    
    u <- as.numeric(start < t & t <= (start + dur)) * lev
    
    dS <- -beta * (1 - u) * S * I
    
    dI <-  beta * (1 - u) * S * I - gamma * I
    
    list(c(dS, dI))
    
  })
  
}



SIR.WMR.ode <- function(t, state, par) {
  
  with(as.list(c(state, par)), {
    
    u <- 0
    
    if (t > t1 && t < t2) {
      
      u <- 1 - 1 / (1 + b * beta * (t2 - t))
      
    }
    
    dS <- -(1 - u) * beta * S * I
    
    dI <-  (1 - u) * beta * S * I - gamma * I
    
    list(c(dS, dI))
    
  })
  
}



run.deterministic <- function(ode_func, epar, init_prop, tf = 100, dt = 0.1) {
  
  times <- seq(0, tf, by = dt)
  
  path  <- ode(c(S = init_prop$S, I = init_prop$I), times, ode_func, epar)
  
  as.data.frame(path)
  
}


# Rate and intervention functions
# This file contains the rate function that has to be forwarded to the ssa.adaptivetau function 
#and the intervention functions needed to use the rate function

u_func_WMR <- function(t, params){
  # the function implements the wait maintain relax function of Miclo et al.  
  # params must contain: t1 (start of intervention), t2 (stop of intervention), beta, gamma, 
  u <-(t>params$t1)*(t<params$t2)* (1 - (1 / ( (params$beta/params$gamma) * params$S)))
  return(u)
  
}

u_func_constant <- function(t, params){
  # the function implements the constant function of Britton et al. 2023  
  # params must contain:t1 (start of intervention), t2 (stop of intervention), lev (intensity of intervention)
  
    u <- params$lev*(t>params$t1)*(t<params$t2)
  return(u)
}


ratefunc.SIR <- function(state, params, t){
  
  # rate function for the ssa.adaptivetau function
  # parameters are: state (S, I , R), params (parameters for the u function and additional parameters needed), t (time), u_func (intervention function, parameters of intervention function are t and params)
  # params must include: t1 (start of intervention), t2 (stop of intervention), N, 
  #beta, gamma, lambda (influx rate), lev (intensity of intervention),  u_func= function(t, params){return(0)}, influx = FALSE
 
  I <- state["I"]
  S <- params$S <- state["S"]
  u <- params$u_func(t, params)
  
  infection_rate <- params$beta*(1-u)*S*I/params$N
  recovery_rate <- params$gamma*I

  if(!params$influx){
  return(c(infection_rate, recovery_rate))
  }else{
    return(c(infection_rate, recovery_rate, influx_rate = params$lambda))
  }
  
}

ratefunc.SIR.age <- function(state, params, t) {
  # rate function for the ssa.adaptivetau function
  # parameters are: state (S, I , R), params (parameters for the u function and additional parameters needed), t (time)
  # params must include: t1 (start of intervention), t2 (stop of intervention), N, 
  # beta, gamma, lambda (influx rate), lev (intensity of intervention),  u_func= function(t, params){return(0)}, influx = FALSE,
  # u_func (intervention function, parameters of intervention function are t and params)
  
  n  <- params$n_age
  NN <- params$N
  
  S <- state[1:n]
  I <- state[(n + 1):(2 * n)]
  # R not needed for rates
  
  # Group population sizes: N_a = f_a * N
  N_a <- params$f * NN
  
  # Update total S in params (so u_func can access it)
  params$S_total <- sum(S)
  
  # Get intervention level (scalar)
  u <- params$u_func(t, params)
  
  # Force of infection per group: 
  # lambda_a = beta * (1-u) * sum_b C_ab * I_b / N_b
  # This is the age-structured generalization of beta*(1-u)*S*I/N
  lambda <- params$beta_scalar * (1 - u) * (params$C %*% (I / N_a))
  
  # Transition rates
  infection_rates <- as.numeric(lambda * S)  # length n_age
  recovery_rates  <- params$gamma * I        # length n_age
  
  # Ensure non-negative (can happen if S or I slightly negative in ODE)
  infection_rates <- pmax(infection_rates, 0)
  recovery_rates  <- pmax(recovery_rates, 0)
  
  if (!params$influx) {
    return(c(infection_rates, recovery_rates))
  } else {
    return(c(infection_rates, recovery_rates, params$lambda_vec))
  }
}

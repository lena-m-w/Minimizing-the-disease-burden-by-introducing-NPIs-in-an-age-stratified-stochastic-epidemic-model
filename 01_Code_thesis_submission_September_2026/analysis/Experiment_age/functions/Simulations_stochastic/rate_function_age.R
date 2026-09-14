# Age-structured rate and intervention functions
# This file contains the rate function that has to be forwarded to the 
# ssa.adaptivetau function and the intervention functions needed to use 
# the rate function.
# Extends rate_and_intervention_functions.R to n age groups.
#
# Key difference from the non-structured version:
#   - beta is a MATRIX (beta_ab = const * C_per_capita), not a scalar
#   - state is a named vector: (S1, S2, ..., I1, I2, ..., R1, R2, ...)
#   - the rate function returns n_age infection rates + n_age recovery rates
#
# In proportions (deterministic ODE): lambda_foi = beta %*% I
# In counts (stochastic SSA):         lambda_foi = beta %*% (I / N)
# because I_proportion = I_count / N


# Rate function for ssa.adaptivetau

ratefunc.SIR<- function(state, params, t) {
  
  # Rate function for the ssa.adaptivetau function (age-structured SIR)
  #
  # parameters are: 
  #   state  : named vector (S1,..,Sn, I1,..,In, R1,..,Rn) in COUNTS
  #   params : list with model parameters (see below)
  #   t      : time
  #
  # params must include:
  #   n_age      : number of age groups
  #   beta       : transmission MATRIX (n_age x n_age), e.g. const * C_per_capita
  #   gamma      : recovery rate (scalar, same for all groups)
  #   N          : total population size
  #   t1, t2     : intervention window
  #   lev        : intervention level (for constant strategy)
  #   R0_target  : basic reproduction number (for WMR strategy)
  #   u_func     : intervention function, signature u_func(t, params)
  #   W          : named n_age x n_age scaling matrix (default: NULL = all ones)
  #   influx     : logical
  #   lambda_vec : influx rates per group (length n_age, if influx = TRUE) 
  
  n <- params$n_age
  # Unpack state (counts)
  # Update S in params so the intervention function can use it
  
  S <- params$S <- state[1:n]
  I <- state[(n + 1):(2 * n)]
  
  # Get intervention level (scalar)
  U_results <- params$u_func(state, params, t)
  one_minus_U <- U_results$one_minus_U 
 
  # Force of infection per group:
  # In the deterministic (proportions): lambda_foi = beta %*% I
  # In stochastic (counts):            lambda_foi = beta %*% (I / N)
  # because I_proportion = I_count / N
  lambda_foi <- (one_minus_U*params$beta) %*% (I)#/ params$N)# Überlegen ob das muss
  
  # Transition rates (must be non-negative)
  infection_rates <- pmax(as.numeric(lambda_foi) * S/params$N, 0)  # length n_age
  recovery_rates  <- pmax(params$gamma * I, 0)             # length n_age
  influx_rates <- params$lambda_influx*S
  
    result <- c(infection_rates, recovery_rates)
  if (params$influx) { # influx proportional zur gemeinschaft machen
    result <-c(result, influx_rates)
  }
    
  if(params$budget_count){
    cost_rate <- U_results$cost_rate#*10000
    result <-c(result, cost_rate)
  }
  
  if(params$intervention_activation){
    intervention_activation <- U_results$intervention_activation*10000
    result <-c(result,intervention_activation)
  }
  if(params$last_time_tracking){
    duration <- U_results$duration#*10000
    if(duration<0){
      duration <- 0
    }
    
    result <-c(result,duration)
  }
    # cat("\n rate_functions:",result)
  return(result)
}


# ══════════════════════════════════════════════════════════════════════════════
# EXAMPLE USAGE
# ══════════════════════════════════════════════════════════════════════════════
#
# library(adaptivetau)
#
# # --- beta matrix (from your existing code) ---
# # beta is already defined as: beta <- const * C_per_capita
# # (a 3x3 matrix, NOT a scalar)
#
# # --- Parameters ---
# N <- 10e6
# n_age <- 3
# 
# params <- list(
#   n_age      = n_age,
#   beta       = beta,          # 3x3 matrix
#   gamma      = gamma,         # 0.2
#   N          = N,
#   t1         = 25,
#   t2         = 45,
#   lev        = 0.75,
#   R0_target  = R0_target,     # 3.0 (needed for WMR)
#   u_func     = u_func_constant_age,
#   influx     = FALSE
# )
#
# # --- Initial state (absolute counts) ---
# I0_total <- 1e-4 * N         # = 1000 for N=10e6 (like Britton)
# I0 <- round(I0_total * f)
# I0[I0 == 0] <- 1
# S0 <- round(f * N) - I0
# R0_init <- rep(0, n_age)
#
# state0 <- c(S0, I0, R0_init)
# names(state0) <- c(paste0("S", 1:n_age),
#                     paste0("I", 1:n_age),
#                     paste0("R", 1:n_age))
#
# # --- Build transition matrix ---
# trans <- build_transitions_age(n_age = 3, influx = TRUE)
#
# # --- Run ---
# result_ssa <- ssa.adaptivetau(state0, trans, ratefunc.SIR.age, params, tf = 200)
# Transition matrix:
# For ssa.adaptivetau: each row is a transition, each column is a state variable.
# State order: (S1, ..., Sn, I1, ..., In, R1, ..., Rn)
#
# Transitions:
#   Row a        : infection in group a  -> S_a - 1, I_a + 1
#   Row n_age + a: recovery  in group a -> I_a - 1, R_a + 1
#   (Row 2*n_age + a: influx into S_a -> S_a + 1,  if influx = TRUE)

# This script was adapted from a code skeleton generated with AI Claude Opus 4.7.
# I have either written or carefully checked and tested the code and take responsibility for its correctness.


build_transitions_age <- function(n_age, influx = FALSE) {
  
  n_trans <- 2 * n_age
  if (influx) n_trans <- n_trans + n_age
  n_state <- 3 * n_age
  
  trans <- matrix(0, nrow = n_trans, ncol = n_state)
  colnames(trans) <- c(paste0("S", 1:n_age),
                       paste0("I", 1:n_age),
                       paste0("R", 1:n_age))
  
  for (a in 1:n_age) {
    # Infection in group a
    trans[a, a]           <- -1   # S_a decreases
    trans[a, n_age + a]   <- +1   # I_a increases
    
    # Recovery in group a
    trans[n_age + a, n_age + a]     <- -1   # I_a decreases
    trans[n_age + a, 2 * n_age + a] <- +1   # R_a increases
  }
  
  if (influx) {
    for (a in 1:n_age) {
      trans[2 * n_age + a, a] <- -1   # S_a decreases
      trans[2 * n_age + a, n_age +a] <- 1   # I_a increases
    }
  }
  
  return(trans)
}


# trans <- build_transitions_age(n_age = 3, influx = TRUE)


#################################################################
# not relevant anymore
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
###############################################################
create_intervention_matrix <- function(u, W = NULL) {
  # Creates the intervention matrix U from age-group-specific levels u
  # and an optional scaling matrix W.
  #
  # Without W (default):  U_ab = sqrt(1 - u_a) * sqrt(1 - u_b)
  # With W:               U_ab = sqrt(1 - W_ab * u_a) * sqrt(1 - W_ab * u_b)
  #
  # u must be a NAMED vector. W (if provided) must be a NAMED symmetric 
  # matrix with row/colnames matching the names of u.
  #
  # u : named vector of length n_age, each entry in [0, 1]
  # W : named n_age x n_age scaling matrix (default: NULL = all ones)
  # tol: tolerance for symmetry check
  # ?????Vielleicht später W nicht nur Matrix sondern Funktion die Matrix macht.
  #
  # Returns: named n_age x n_age matrix U
  
  age_names <- names(u)
  
  s <- sqrt(1 - u)
  U <- outer(s, s)
  rownames(U) <- age_names
  colnames(U) <- age_names
  if (is.null(W)) {
    return(U)
  }
  
  return(U*W)
  }

# ── Rate function for ssa.adaptivetau ─────────────────────────────────────────

ratefunc.SIR <- function(state, params, t) {
  
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
  u <- params$u_func(t, params)
  U <- create_intervention_matrix(u = u, W = params$W) 
  # Force of infection per group:
  # In the deterministic (proportions): lambda_foi = beta %*% I
  # In stochastic (counts):            lambda_foi = beta %*% (I / N)
  # because I_proportion = I_count / N
  lambda_foi <- (U*params$beta) %*% (I)# / params$N)) Überlegen ob das muss
  
  # Transition rates (must be non-negative)
  infection_rates <- pmax(as.numeric(lambda_foi) * S, 0)  # length n_age
  recovery_rates  <- pmax(params$gamma * I, 0)             # length n_age
  
  if (!params$influx) { # influx proportional zur gemeinschaft machen
    return(c(infection_rates, recovery_rates))
  } else {
    return(c(infection_rates, recovery_rates, params$lambda_vec))
  }
}


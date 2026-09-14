# check again

init_age_state <- function(I0_total, f, age_names, budget_count = TRUE,
                           intervention_activation = TRUE, last_time_tracking = TRUE ) {
  # Creates a named initial state vector for the age-structured SIR model
  #
  # I0_total : total initial infected as proportion (e.g. 1e-4, like Britton)
  # f        : population fractions per age group (must sum to 1)
  # age_names: character vector of age group names
  #
  # Returns: named integer vector (S1,..,Sn, I1,..,In, R1,..,Rn) in COUNTS
  
  n <- length(f)
  
  # Distribute infected proportionally across groups
  I0 <- (I0_total * f )
  I0[I0 == 0] <- 1   # at least 1 infected per group to seed epidemic
  
  S0 <- (f) - I0
  R0 <- rep(0, n)
  
  state0 <- c(S0, I0, R0)
  names(state0) <- c(paste0("S", 1:n),
                     paste0("I", 1:n),
                     paste0("R", 1:n))
  if(budget_count){
    state0 <- c(state0,budget_count = 0)
  }
  if(intervention_activation){
    state0 <- c(state0, intervention_activation = 0)
  }
  if(last_time_tracking){
    state0 <- c(state0, last_time_tracking = 0)
  }
  # Sanity checks
  if (any(state0 < 0)) {
    stop("Negative S0 — I0_total is too large for the population size")
  }
  if (abs(sum(f) - 1) > 1e-6) {
    warning("Population fractions f do not sum to 1")
  }
  
  cat("Initial state:\n")
 # for (a in 1:n) {
   # cat(sprintf("  %s: S=%d, I=%d, R=%d (N_a=%d)\n",
   #             age_names[a], S0[a], I0[a], R0[a], S0[a] + I0[a]))
#  }
 # cat(sprintf("  Total: S=%d, I=%d, R=%d, N=%d\n",
  #            sum(S0), sum(I0), sum(R0), sum(state0)))
  
  return(state0)
}
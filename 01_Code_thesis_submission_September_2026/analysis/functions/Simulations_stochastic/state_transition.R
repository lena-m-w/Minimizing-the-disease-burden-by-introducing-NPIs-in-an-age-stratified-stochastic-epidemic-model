# Transition function

# For ssa.adaptivetau: each row is a transition, each column is a state variable.
# State order: (S1, ..., Sn, I1, ..., In, R1, ..., Rn)
#
# Transitions:
#   Row a        : infection in group a  → S_a - 1, I_a + 1
#   Row n_age + a: recovery  in group a  → I_a - 1, R_a + 1
#   (Row 2*n_age + a: influx into S_a    → S_a + 1,  if influx = TRUE)

build_transitions_age <- function(n_age, influx = FALSE, budget_count = TRUE,
                                  intervention_activation = TRUE, last_time_tracking = TRUE ){

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
  if(budget_count){
    trans <- cbind(trans, budget_count = rep(0,n_trans),
                   intervention_activation = rep(0,n_trans), 
                   last_time_tracking = rep(0,n_trans))
    matrix_temp <- cbind(matrix(rep(0,n_state*3), nrow = 3), diag(c(1,1,1)))
    trans <- rbind(trans, matrix_temp)
    
  }

  return(t(trans))
}

# ## check again
# build_transitions_age <- function(n_age, influx = FALSE) {
# 
#   state_names <- c(paste0("S", 1:n_age),
#                    paste0("I", 1:n_age),
#                    paste0("R", 1:n_age))
# 
#   transitions <- list()
# 
#   for (a in 1:n_age) {
#     # Infection in group a: S_a -1, I_a +1
#     tr <- setNames(c(-1, +1), c(paste0("S", a), paste0("I", a)))
#     transitions[[length(transitions) + 1]] <- tr
#   }
# 
#   for (a in 1:n_age) {
#     # Recovery in group a: I_a -1, R_a +1
#     tr <- setNames(c(-1, +1), c(paste0("I", a), paste0("R", a)))
#     transitions[[length(transitions) + 1]] <- tr
#   }
# 
#   if (influx) {
#     for (a in 1:n_age) {
#       # Influx into S_a: S_a +1
#       tr <- setNames(c(+1), paste0("S", a))
#       transitions[[length(transitions) + 1]] <- tr
#     }
#   }
# 
#   return(transitions)
# }



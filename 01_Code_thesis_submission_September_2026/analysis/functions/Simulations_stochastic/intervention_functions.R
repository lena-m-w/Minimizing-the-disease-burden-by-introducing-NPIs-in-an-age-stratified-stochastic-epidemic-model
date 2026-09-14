# Intervention functions
#bei u bei denen buget und so nicht in return error einbauen

# Intervention functions
#bei u bei denen buget und so nicht in return error einbauen

u_func_constant_analytical_rr_age2 <- function(state, params, t) {
  # Stateless analytical intervention for up to 2 age groups, with a
  # symmetric matrix-valued `lev`, per-group budgets, and per-cell start times.
  #
  # Special case: when both diagonals of lev are zero (only the cross term
  # exists), alloc_12 is overridden -- both groups must spend everything on
  # the cross, and the duration is capped by whichever group has less
  # per-capita budget.
  # params also includes the boolen parameter double_non_diagonal_costs which indicates 
  # if the costs are calculated as  \sum_i\sum_j||u||_1^{(ij)} (TRUE) or as \sum_{i\leq j}||u||_1^{(ij)} (FALSE)
  
  n <- params$n_age
  if (n > 2) stop("This intervention function only supports n_age <= 2.")
  
  lev      <- params$lev
  pi_v     <- params$f
  if (lev[1, 1] == 0 && lev[2, 2] == 0) {
    B <- B_total  <- sum(params$budget)
    
  }else{
    B        <- params$budget
    if(lev[1, 1] == 0 & lev[1, 2] == 0 ){
      B[1] <- 0
      B[2] <- sum(params$budget)
    }
  }
  
  t1       <- params$t1
  alloc_12 <- if (n == 2) params$alloc_12 else 0
  
  # ---- promote scalars to matrices --------------------------------------
  if (length(lev) == 1) lev <- matrix(lev, n, n)
  if (length(t1)  == 1) t1  <- matrix(t1,  n, n)
  
  # ---- sanity checks ----------------------------------------------------
  if (!is.matrix(lev) || any(dim(lev) != n))
    stop("lev must be n_age x n_age or a scalar.")
  if (!is.matrix(t1) || any(dim(t1) != n))
    stop("t1 must be n_age x n_age or a scalar.")
  if (n == 2 && !isSymmetric(lev))
    stop("lev must be symmetric.")
  if (n == 2 && !isSymmetric(t1))
    stop("t1 must be symmetric.")
  if (length(params$budget) != n)
    stop("budget must have length n_age.")
  if (n == 2 && (alloc_12 < 0 || alloc_12 > 1))
    stop("alloc_12 must be in [0, 1].")
  if (any(pi_v <= 0))
    stop("All entries of pi (params$f) must be > 0.")
  
  # ---- compute durations -----------------------------------------------
  dur <- matrix(0, n, n)
  
  if (n == 1) {
    if (lev[1, 1] <= 0)
      stop("lev[1,1] must be > 0 to spend budget on u_11.")
    dur[1, 1] <- B[1] / (lev[1, 1] * pi_v[1])
    
  } else {  # n == 2
    
    # -------- special case: cross-only (both diagonals zero) -------------
    if (lev[1, 1] == 0 && lev[2, 2] == 0) {
      if (lev[1, 2] <= 0){
        U_mat <- lev * 0
        
        dimnames(U_mat) <- list(params$age_names, params$age_names)
        return(list(one_minus_U = 1 - U_mat))
      }
        
      # both groups must spend everything on the cross; duration is capped
      # by whichever group runs out first.
      # when only the interaction of the two groups is intervened
      # then dur_12 = budget_1/(lev_12*pi_1)= budget_2/(lev_12*pi_2)
      # and budget_1 has to be budget_total* (pi_1/(pi_2+pi_1))
      x = pi_v[1]/sum(pi_v)
      B <- c(B_total*x, B_total*(1-x)) 
      if((B[1] / (lev[1, 2] * pi_v[1])-B[2] / (lev[1, 2] * pi_v[2]))>0.0001){
        stop("Error in budget calculations, change function.")
      }
      #if(params$double_non_diagonal_costs)
      t_12 <- B[1] / (lev[1, 2] * pi_v[1])
      
        
    
      dur[1, 2] <- t_12
      dur[2, 1] <- t_12
      
    } else {
      # general case: alloc_12 governs group 1's split
      t_12 <- 0
      if (alloc_12 > 0) {
        if (lev[1, 2] <= 0) {
          
          if (B[1] == 0) {
            t_12 <- 0
          } else {
            stop("alloc_12 > 0 but lev[1,2] = 0: cannot spend budget on a zero-level intervention.")
          }
        } else {
          t_12 <- (B[1] * alloc_12) / (lev[1, 2] * pi_v[1])
        }
        dur[1, 2] <- t_12
        dur[2, 1] <- t_12
      }  
      # group 1 diagonal
      if (alloc_12 < 1) {
        if (lev[1, 1] <= 0)
          stop("alloc_12 < 1 but lev[1,1] = 0: leftover group-1 budget cannot fund u_11.")
        dur[1, 1] <- (B[1] * (1 - alloc_12)) / (lev[1, 1] * pi_v[1])
      }
      
      # group 2 diagonal
      B2_remaining <- B[2] - t_12 * lev[1, 2] * pi_v[2]
      if (B2_remaining < 0)
        stop("Cross term (1,2) exhausts more than group 2's budget. Reduce alloc_12, lev[1,2], or increase B_2.")
      if (B2_remaining > 0) {
        if (lev[2, 2] <= 0)
          stop("Group 2 has remaining budget but lev[2,2] = 0: cannot fund u_22.")
        dur[2, 2] <- B2_remaining / (lev[2, 2] * pi_v[2])
      }
    }
  }
  
  # ---- build U(t) -------------------------------------------------------
  active_mat <- (t > t1) & (t < t1 + dur)
  U_mat <- lev * active_mat
  
  dimnames(U_mat) <- list(params$age_names, params$age_names)
  return(list(one_minus_U = 1 - U_mat, dur = dur))
}


u_func_constant_rr <- function(state, params, t) {
  lev <- params$lev
  if(!((is.matrix(lev)&length(lev)>1)|length(lev)==1)|!all(lev ==t(lev))){
    stop("lev must be a symetric matrix!")
  }
  if(is.matrix(lev)&length(lev)>1){
    cost_rate <- rowSums((lev*params$f)%*%diag(c(1,rep(0.5,(params$n_age-1)))))
  }else{
    cost_rate <- lev*params$f
  }       # shared budget
  t2        <- params$t1 + params$budget / cost_rate
  active    <- (t > params$t1) & (t < t2)
  u         <- params$lev * active
  names(u)  <- params$age_names
  n <- params$n_age

  return(list(one_minus_U  = 1 - u))
}

# # Intervention functions
# #bei u bei denen buget und so nicht in return error einbauen
# u_func_constant_rr_S <- function(state, params, t) {
#   lev <- params$lev
#   if(is.matrix(lev)&length(lev)>1){
#     cost_rate <- rowSums(lev*params$f%*%diag(c(1,rep(0.5,(params$n_age-1)))))
#   }else{
#     cost_rate <- lev*params$f
#   }       # shared budget
#   intervention_activation <- state[2*n+3]
#   t2        <- params$t1 + params$budget / cost_rate
#   active    <- (t > params$t1) & (t < t2)
#   u         <- params$lev * active
#   names(u)  <- params$age_names
#   n <- params$n_age
#   
#   return(list(one_minus_U  = 1 - u))
# }
# 
# u_func_constant_rr <- function(state, params,t) {
#   # Age-dependent constant intervention (Britton et al. 2023)
#   # relative reduction
#   #
#   # The intervention runs at level lev from t1 until the budget is
#   # exhausted. t2 set when no budget is left.
#   #
#   # params must contain:
#   #   t1        : start time(s) — scalar or named vector
#   #   lev       : level(s) — scalar or named vector
#   #   budget    : cumulative cost budget — see budget_type
#   #   budget_type : "per_group" or "shared"
#   #   age_names : character vector of age group names
#   #
#   # params$state_tb must be an ENVIRONMENT (mutable) containing:
#   #   t_old       : time of the previous call
#   #   budget_used : per-group cumulative cost (named numeric vector)
#   #
#   # Budget types:
#   #   "per_group" — each group has its own budget
#   #       budget = c("0-19" = 10, "20-64" = 5, "65+" = 15)
#   #       duration_a = budget_a / lev_a, groups stop independently
#   #
#   #   "shared" — one total budget, cost = sum of all u_a per unit time
#   #       budget = 15 (scalar)
#   #       cost rate = sum(lev) when active, all groups stop together
# 
#   # wie ich die Funktion umschreibe:
#   # 1. schauen ob die startzeit erreicht ist.
#   # suche nach zeitspeicherung
#   # 2. berechnen wie hoch die kosten bis jetzt sind.
#   # 3. schauen, ob die kosten das budget überschreiten
#   # 4. Wenn kosten das budget überschreiten würden stoppen.
#   # 5. environment kreieren für params$state_tb
#   #Reset between runs: because the state_tb lives in an environment attached to params, recreate params$state_tb before each simulation run, otherwise budget from the previous run carries over.
# 
#   n <- params$n_age
# 
#   if(params$influx){
#   budget_used <- state[2*n+1]#/10000
# 
#   t_old <- state[2*n+3]#/10000
#   #cat("\n state:", budget_used, t_old)
#   }else{
#     budget_used <- state[2*n+2]#/10000
# 
#     t_old <- state[2*n+4]#/10000
#    # cat("\n state_infl:", budget_used, t_old, params$t1)
#   }
#   intervention_activation <- (t > params$t1)
#   cat("\n state_infl:", budget_used, t_old, t, params$t1, intervention_activation)
#   active <- intervention_activation&(budget_used<params$budget)
#   u <- params$lev * active
#   if(any(u>1|u<0)){
#     stop("Error: intervention is not between 0 and 1" )
#   }
#   names(u) <- params$age_names
# 
#   dur <- t-t_old
#   if(dur >= 0){
#     if(is.matrix(u)&length(u)>1){
#       cost_rate <- rowSums(u*params$f%*%diag(c(1,rep(0.5,(params$n_age-1)))))
#     }else{
#       cost_rate <- u*params$f
#     }
#     budget_count <- budget_used + cost_rate*dur
# 
# 
#   }else{
#     cost_rate <- 0
#   }
#   active <- intervention_activation&(budget_used<=params$budget)
#   u <- params$lev * active
#   #cat("state:", state)
#   return(list(one_minus_U  = 1 - u, cost_rate = cost_rate, intervention_activation = intervention_activation/10000, duration = dur))
# 
#   }

u_func_none <- function(state, params, t) {
  # No intervention (baseline / do-nothing scenario)
  n_age <- length(params$age_names)
  u <- matrix(rep(0, n_age*n_age), ncol = n_age)
  names(u) <- params$age_names

  return(list(one_minus_U  = 1-u))

}

# ## old intervention function
# 
# u_func_constant_age_calculate_budget_once <- function(t, params) {
#   # Age-dependent constant intervention (Britton et al. 2023)
#   #
#   # The intervention runs at level lev from t1 until the budget is
#   # exhausted. t2 is computed from the budget, not specified directly.
#   #
#   # params must contain:
#   #   t1        : start time(s) — scalar or named vector
#   #   lev       : level(s) — scalar or named vector
#   #   budget    : cumulative cost budget — see budget_type
#   #   budget_type : "per_group" or "shared"
#   #   age_names : character vector of age group names
#   #
#   # Budget types:
#   #   "per_group" — each group has its own budget
#   #       budget = c("0-19" = 10, "20-64" = 5, "65+" = 15)
#   #       duration_a = budget_a / lev_a, groups stop independently
#   #
#   #   "shared" — one total budget, cost = sum of all u_a per unit time
#   #       budget = 15 (scalar)
#   #       cost rate = sum(lev) when active, all groups stop together
# 
#   # wie ich die Funktion umschreibe: 1. schauen ob die startzeit erreicht ist.
#   # 2. berechnen wie hoch die kosten bis jetzt sind.
#   # 3. schauen, ob die kosten das budget überschreiten
#   # 4. Wenn kosten das budget überschreiten würden stoppen.
# 
#   t2 <- params$t1 + params$budget / params$lev
# 
#   active <- (t > params$t1) * (t < t2)
#   u <- params$lev * active
#   u <- pmax(0, pmin(u, 1))
#   names(u) <- params$age_names
#   return(u)
# }
# 
# u_func_none <- function(t, params) {
#   # No intervention (baseline / do-nothing scenario)
# 
#   u <- rep(0, length(params$age_names))
#   names(u) <- params$age_names
# 
#   return(1-u)
# 
# }
# 
# create_intervention_matrix <- function(u, W = NULL) {
#   # Creates the intervention matrix U from age-group-specific levels u
#   # and an optional scaling matrix W.
#   #
#   # Without W (default):  U_ab = sqrt(1 - u_a) * sqrt(1 - u_b)
#   # With W:               U_ab = sqrt(1 - W_ab * u_a) * sqrt(1 - W_ab * u_b)
#   #
#   # u must be a NAMED vector. W (if provided) must be a NAMED symmetric 
#   # matrix with row/colnames matching the names of u.
#   #
#   # u : named vector of length n_age, each entry in [0, 1]
#   # W : named n_age x n_age scaling matrix (default: NULL = all ones)
#   # tol: tolerance for symmetry check
#   # ?????Vielleicht später W nicht nur Matrix sondern Funktion die Matrix macht.
#   #
#   # Returns: named n_age x n_age matrix U
#   
#   age_names <- names(u)
#   
#   s <- (1 - u)
#   U <- outer(s, s)
#   rownames(U) <- age_names
#   colnames(U) <- age_names
#   if (is.null(W)) {
#     return(U)
#   }
#   
#   return(U*W)
# }
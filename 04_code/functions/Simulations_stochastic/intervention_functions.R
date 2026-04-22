# Intervention functions 

u_func_constant_age <- function(t, params) {
  # Age-dependent constant intervention (Britton et al. 2023)
  #
  # The intervention runs at level lev from t1 until the budget is
  # exhausted. t2 is computed from the budget, not specified directly.
  #
  # params must contain:
  #   t1        : start time(s) — scalar or named vector
  #   lev       : level(s) — scalar or named vector
  #   budget    : cumulative cost budget — see budget_type
  #   budget_type : "per_group" or "shared"
  #   age_names : character vector of age group names
  #
  # Budget types:
  #   "per_group" — each group has its own budget
  #       budget = c("0-19" = 10, "20-64" = 5, "65+" = 15)
  #       duration_a = budget_a / lev_a, groups stop independently
  #
  #   "shared" — one total budget, cost = sum of all u_a per unit time
  #       budget = 15 (scalar)
  #       cost rate = sum(lev) when active, all groups stop together
  

  t2 <- params$t1 + params$budget / params$lev

  
  active <- (t > params$t1) * (t < t2)
  u <- params$lev * active
  u <- pmax(0, pmin(u, 1))
  names(u) <- params$age_names
  return(u)
}

u_func_none <- function(t, params) {
  # No intervention (baseline / do-nothing scenario)
  
  u <- rep(0, length(params$age_names))
  names(u) <- params$age_names
  
  return(u)
  
}
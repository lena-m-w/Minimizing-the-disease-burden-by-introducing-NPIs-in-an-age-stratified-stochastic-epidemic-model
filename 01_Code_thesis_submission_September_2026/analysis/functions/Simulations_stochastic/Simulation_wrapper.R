# Simulation Wrapper

simulate.stochastic <- function( params, params_adaptive_tau) {
  # params must include: t1 (start of intervention), t2 (stop of intervention), N, 
  #beta, gamma, lambda (influx rate), lev (intensity of intervention),  u_func= function(t, params){return(0)}, influx = FALSE
  
  #influx (TRUE /FALSE)
  # define a maxtau, wie finde ich raus, was ein sinnvolles maxtau ist. Was sind die Werte von tau? 
  # if influx = TRUE (default: FALSE) for one infectious import one susceptible less in population 
  # under default parameters, see transitions_influx. 
  
  if (is.null(params$u_func)) {
    params$u_func <- u_func_none
  }
  if (is.null(params$influx)) {
    params$influx <- FALSE
  }
  
  n <- params$n_age

  # ── Run ────────────────────────────────────────────────────────────

  result <- ssa.adaptivetau(
    init.values       = params_adaptive_tau$init,
    transitions       = params_adaptive_tau$transitions,
    rateFunc          = params_adaptive_tau$ratefunc,
    params            = params,
    tf                = params_adaptive_tau$tf,
    jacobianFunc      = params_adaptive_tau$jacobianFunc,
    maxTauFunc        = params_adaptive_tau$maxTauFunc,
    deterministic     = params_adaptive_tau$deterministic,
    halting           = params_adaptive_tau$halting,
    relratechange     = params_adaptive_tau$relratechange,
    tl.params         = params_adaptive_tau$tl.params,
    reportTransitions = params_adaptive_tau$reportTransitions
  )
  
  # ── Convert to data.frame with proportions ─────────────────────────
  df <- as.data.frame(result)
  
  age_names <- params$age_names
  
  # Per-group proportions
  for (a in 1:n) {
    S_col <- paste0("S", a)
    I_col <- paste0("I", a)
    R_col <- paste0("R", a)
    N_a   <- params$f[a] * params$N
    
    df[[paste0("S_", age_names[a], "_prop")]] <- df[[S_col]] / N
    df[[paste0("I_", age_names[a], "_prop")]] <- df[[I_col]] / N
    df[[paste0("R_", age_names[a], "_prop")]] <- df[[R_col]] / N
  }
  
  # Totals
  S_cols <- paste0("S", 1:n)
  I_cols <- paste0("I", 1:n)
  R_cols <- paste0("R", 1:n)
  
  if(is.null(ncol(df[, S_cols]))){
    df$S_total     <- df[, S_cols]
    df$I_total     <- df[, I_cols]
    df$R_total     <- df[, R_cols]
  } else{
    df$S_total     <- rowSums(df[, S_cols])
    df$I_total     <- rowSums(df[, I_cols])
    df$R_total     <- rowSums(df[, R_cols])
  }
  
  df$S_total_prop <- df$S_total / params$N
  df$I_total_prop <- df$I_total / params$N
  df$R_total_prop <- df$R_total / params$N
  
  return(df)
}
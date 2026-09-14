# Vielleciht extra funktion summary statistics 
#
run.replicates <- function(n_rep = 200, init,transitions_influx = rbind(S = c(-1,  0,  -1), I = c(+1, -1, +1), R = c( 0, +1,  0)),
                           ratefunc, params, tf = 100, dt = 0.5, tl.params = list(epsilon = 0.05)) {
 # params, params_adaptive_tau)
  # Function to run a single replicate
  
  run_one <- function(i) {
    
    sim <- simulate.stochastic(params = params, params_adaptive_tau = params_adaptive_tau)
    
    simS = sim$S
    simI_prop = sim$I_prop
    
    ti <- 1 - tail(simS, 1) / simS[1]  # total number of infected
    
    pk <- max(simI_prop)
    
    pt <- sim$time[which.max(simI_prop)]
    
    fs <- tail(simS, 1)
    
    head(sim)
    
    list(simulation_results = sim,
         
         total_number_infected = ti, peak_I = pk, peak_time = pt, final_S = fs)
    
  }
  
  
  # Summary statistics per replicate
  
  total_number_infected <- numeric(n_rep)
  peak_I          <- numeric(n_rep)
  peak_time       <- numeric(n_rep)
  final_S         <- numeric(n_rep)
  
  
  
  # Run (parallel or serial)
  
  cat(sprintf("Running %d replicates (N=%d)...\n", n_rep, params$N))
  
  t_start <- Sys.time()
  
  
  
#  if (n_cores > 1 && .Platform$OS.type == "unix") {results_list <- mclapply(1:n_rep, run_one, mc.cores = n_cores)
    
    pb <- txtProgressBar(min = 0, max = n_rep, style = 3)
    
    results_list <- lapply(1:n_rep, function(i) {
      
      res <- run_one(i)
      
      setTxtProgressBar(pb, i)
      
      res
      
    })
    
    close(pb)
    
  
  
  
  elapsed <- as.numeric(difftime(Sys.time(), t_start, units = "secs"))
  
  cat(sprintf("Completed in %.1f seconds (%.2f s/replicate)\n",
              
              elapsed, elapsed / n_rep))
  
  
  
  # Unpack
  
  for (i in 1:n_rep) {
    
    total_number_infected[i] <- results_list[[i]]$total_number_infected
    
    peak_I[i]          <- results_list[[i]]$peak_I
    
    peak_time[i]       <- results_list[[i]]$peak_time
    
    final_S[i]         <- results_list[[i]]$final_S
    
  }
  
  all_simulation_results <- lapply(results_list, function(x) x$simulation_results)
  
  return(list(
    
    total_number_infected = total_number_infected,
    
    peak_I          = peak_I,
    
    peak_time       = peak_time,
    
    final_S         = final_S,
    
    params          = params,
    
    n_rep           = n_rep,
    
    tf              = tf,
    
    all_simulation_results = all_simulation_results
    
  ))
  
}



# Run replicates (age-structured)
run.replicates.age <- function(n_rep = 200, params, params_adaptive_tau) {
  
  n_age <- params$n_age
  age_names <- params$age_names
  
  run_one <- function(i) {
    
    sim <- simulate.stochastic(params = params, 
                               params_adaptive_tau = params_adaptive_tau)
    
    # ── Total incidence: 1 - S_total(end) / S_total(0) ──────────
    S_total_start <- sim$S_total[1]
    S_total_end   <- tail(sim$S_total, 1)
    ti <- 1 - S_total_end / S_total_start # total number of infected
    
    # ── Per-group incidence ──────────────────────────────────────
    ti_group <- numeric(n_age)
    names(ti_group) <- age_names
    for (a in 1:n_age) {
      S_col <- paste0("S", a)
      ti_group[a] <- 1 - tail(sim[[S_col]], 1) / sim[[S_col]][1]
    }
    
    # Peak total prevalence
    pk       <- max(sim$I_total_prop) # peak hight 
    pt       <- sim$time[which.max(sim$I_total_prop)] # time of peak
    
    # Per-group peak prevalence 
    pk_group <- numeric(n_age)
    pt_group <- numeric(n_age)
    names(pk_group) <- names(pt_group) <- age_names
    for (a in 1:n_age) {
      I_prop_col <- paste0("I_", age_names[a], "_prop")
      pk_group[a] <- max(sim[[I_prop_col]])
      pt_group[a] <- sim$time[which.max(sim[[I_prop_col]])]
    }
    
    # Final susceptibles
    fs <- S_total_end / params$N
    fs_group <- numeric(n_age)
    
    for (a in 1:n_age) {
      S_col <- paste0("S", a)
      fs_group[a] <- tail(sim[[S_col]], 1) / params$N
    }
    
    list(simulation_results    = sim,
         total_incidence       = ti,
         incidence_by_group    = ti_group,
         peak_I                = pk,
         peak_time             = pt,
         peak_I_by_group       = pk_group,
         peak_time_by_group    = pt_group,
         final_S_total_prop    = fs,
         final_S_group_prop    = fs_group)
  }
  
  # ── Run ──────────────────────────────────────────────────────────
  cat(sprintf("Running %d replicates (N=%d, %d age groups)...\n", 
              n_rep, params$N, n_age))
  t_start <- Sys.time()
  
  pb <- txtProgressBar(min = 0, max = n_rep, style = 3)
  results_list <- lapply(1:n_rep, function(i) {
    res <- run_one(i)
    setTxtProgressBar(pb, i)
    res
  })
  close(pb)
  
  elapsed <- as.numeric(difftime(Sys.time(), t_start, units = "secs"))
  cat(sprintf("Completed in %.1f seconds (%.2f s/replicate)\n",
              elapsed, elapsed / n_rep))
  
  # ── Unpack scalars ───────────────────────────────────────────────
  total_incidence <- sapply(results_list, `[[`, "total_incidence")
  peak_I          <- sapply(results_list, `[[`, "peak_I")
  peak_time       <- sapply(results_list, `[[`, "peak_time")
  final_S_total_prop  <- sapply(results_list, `[[`, "final_S_total_prop")
  
  # ── Unpack per-group (each becomes a matrix: n_rep x n_age) ─────
  incidence_by_group   <- t(sapply(results_list, `[[`, "incidence_by_group"))
  peak_I_by_group      <- t(sapply(results_list, `[[`, "peak_I_by_group"))
  peak_time_by_group   <- t(sapply(results_list, `[[`, "peak_time_by_group"))
  final_S_group_prop   <- sapply(results_list, `[[`, "final_S_group_prop")
  
  all_simulation_results <- lapply(results_list, `[[`, "simulation_results")
  
  return(list(
    total_incidence        = total_incidence,
    peak_I                 = peak_I,
    peak_time              = peak_time,
    final_S_total_prop     = final_S_total_prop,
    incidence_by_group     = incidence_by_group,     # matrix: n_rep x n_age
    peak_I_by_group        = peak_I_by_group,        # matrix: n_rep x n_age
    peak_time_by_group     = peak_time_by_group,     # matrix: n_rep x n_age
    params                 = params,
    params_adaptive_tau    = params_adaptive_tau,
    n_rep                  = n_rep,
    all_simulation_results = all_simulation_results
  ))
}


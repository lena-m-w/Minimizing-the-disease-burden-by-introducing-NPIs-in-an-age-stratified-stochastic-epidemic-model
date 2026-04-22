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



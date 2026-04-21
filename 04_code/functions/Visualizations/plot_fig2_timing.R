# 9. Stochastic optimization: grid search for optimal lockdown start time

lev <- 0.75

dur <- 20

start_grid <- seq(0, 35, by = 0.5)

n_rep_opt  <- 100  # replicates per grid point (increase for more precision)



opt_results <- data.frame(
  
  start = start_grid,
  
  mean_incidence   = NA,
  
  sd_incidence     = NA,
  
  median_incidence = NA,
  
  q25_incidence   = NA,
  
  q75_incidence   = NA
  
)



cat(sprintf("Evaluating %d start times x %d replicates each...\n",
            
            length(start_grid), n_rep_opt))



for (j in seq_along(start_grid)) {
  
  s <- start_grid[j]
  
  params_j <- list(beta = beta, gamma = gamma, N = N,
                   
                   start = s, dur = dur, lev = lev)
  
  
  
  res_j <- run.replicates(n_rep_opt, init, transitions, ratefunc.SIR.1LD,
                          
                          params_j, tf = tf, dt = dt, n_cores = n_cores)
  
  
  
  opt_results$mean_incidence[j]   <- mean(res_j$total_incidence)
  
  opt_results$sd_incidence[j]     <- sd(res_j$total_incidence)
  
  opt_results$median_incidence[j] <- median(res_j$total_incidence)
  
  opt_results$q25_incidence[j]   <- quantile(res_j$total_incidence, 0.25)
  
  opt_results$q75_incidence[j]   <- quantile(res_j$total_incidence, 0.75)
  
  opt_results$q025_incidence[j]   <- quantile(res_j$total_incidence, 0.025)
  
  opt_results$q975_incidence[j]   <- quantile(res_j$total_incidence, 0.975)
  
  
  
  cat(sprintf("  start=%.1f: mean_incidence=%.4f (SD=%.4f)\n",
              
              s, opt_results$mean_incidence[j], opt_results$sd_incidence[j]))
  
}



# Deterministic comparison: total incidence for each start time

det_incidence_by_start <- sapply(start_grid, function(s) {
  
  det_j <- run.deterministic(
    
    SIR.1LD.ode,
    
    epar = list(beta = beta, gamma = gamma, start = s, dur = dur, lev = lev),
    
    init_prop, tf = tf, dt = 0.1
    
  )
  
  1 - tail(det_j$S, 1) / det_j$S[1]
  
})



# Plot

par(mfrow = c(1, 1), mar = c(5, 5, 3, 1))



plot(opt_results$start, opt_results$median_incidence, type = "b",
     
     pch = 16, col = "blue", lwd = 2,
     
     xlab = "Lockdown start time (days)",
     
     ylab = "Total incidence",
     
     ylim = range(c(opt_results$q25_incidence, opt_results$q75_incidence,
                    
                    det_incidence_by_start)),
     
     main = sprintf("Optimal start time: stochastic (N=%s) vs deterministic",
                    format(N, big.mark = ",")
     ))



# IQR envelope

polygon(c(opt_results$start, rev(opt_results$start)),
        
        c(opt_results$q25_incidence, rev(opt_results$q75_incidence)),
        
        col = rgb(0, 0, 1, 0.15), border = NA)



# Deterministic curve

lines(start_grid, det_incidence_by_start, col = "red", lwd = 2, lty = 2)



# Mark deterministic optimum

det_opt_idx <- which.min(det_incidence_by_start)

abline(v = start_grid[det_opt_idx], col = "red", lty = 3)



# Mark stochastic optimum

stoch_opt_idx <- which.min(opt_results$median_incidence)

abline(v = start_grid[stoch_opt_idx], col = "blue", lty = 3)







legend("bottomright",
       
       legend = c("median, IQR",
                  
                  "Deterministic",
                  
                  sprintf("Det. optimum: t*=%.1f", start_grid[det_opt_idx]),
                  
                  sprintf("Stoch. optimum: t*=%.1f", start_grid[stoch_opt_idx])),
       
       col = c("blue", "red", "red", "blue"),
       
       lwd = c(2, 2, 1, 1), lty = c(1, 2, 3, 3), pch = c(16, NA, NA, NA),
       
       bty = "n", cex = 0.8)



cat(sprintf("\nDeterministic optimal start: %.1f (incidence=%.4f)\n",
            
            start_grid[det_opt_idx], min(det_incidence_by_start)))

cat(sprintf("Stochastic optimal start:   %.1f (mean incidence=%.4f)\n",
            
            start_grid[stoch_opt_idx], min(opt_results$mean_incidence)))


# one run, lockdown in the very beginning 

params_early <- list(beta = beta, gamma = gamma, N = N,
                     start = 0, dur = 20, lev = 0.75)

res_early <- run.replicates(100, init, transitions, ratefunc.SIR.1LD,
                            params_early, tf = tf, dt = dt, n_cores = n_cores)

sort(res_early$total_incidence)

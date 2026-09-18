# Combined experiment 
#   (3) N = 10,000      no influx
#   (4) N = 10,000      with influx
#   (5) N = 100,000     no influx
#   (6) N = 100,000     with influx
#
# AI use: parts of the script were written with the help of AI assistant (Anthropic Claude)
# based on detailed description from the author (code structure, variable naming, plotting,
# some individual lines, debugging support, and the skeleton of the parameter
# file.). All code was written or reviewed and tested by the author.


final_plot_file_name <- "/Exp3-6_combined_2x2_0409.png"
final_data_file_name <- "/Exp3-6_combined_2x2_0409.rds"

# Parameters
# source(paste0(path_wd, path_experiments, path_influx,
#              "/setup_parameters_exp3-6_combined.R"))
# scenarios

scenarios <- list(
  list(N = 10000,  influx = FALSE, label = "N = 10,000 - no influx"),
  list(N = 10000,  influx = TRUE,  label = "N = 10,000 - with influx"),
  list(N = 100000, influx = FALSE, label = "N = 100,000 - no influx"),
  list(N = 100000, influx = TRUE,  label = "N = 100,000 - with influx")
)

all_results <- vector("list", length(scenarios))

# functions
stochastic_run <- function(start_time){
  s <- start_grid[start_time]
  
  params_j <- list(
    n_age      = n_age,
    age_names  = age_names,
    beta       = beta,
    gamma      = gamma,
    N          = N,
    f          = f,
    t1         = s,
    lev        = lev,
    budget     = budget_vec,
    u_func     = u_func_constant_rr,
    influx     = influx,
    lambda_influx = lambda_influx_vec,
    W          = NULL,
    budget_count            = budget_count,
    intervention_activation = intervention_activation,
    last_time_tracking      = last_time_tracking
  )
  
  res_j <- run.replicates.age(
    n_rep               = n_rep_opt,
    params              = params_j,
    params_adaptive_tau = params_adaptive_tau
  )
  
  ti <- res_j$total_incidence
  return(ti)
}

deterministic_run <- function(start_time) {
  
  epar_j <- list(
    n_age     = n_age,
    age_names = age_names,
    beta      = beta,
    gamma     = gamma,
    f         = f,
    t1        = start_time,
    lev       = lev,
    budget    = budget_vec,
    u_func    = u_func_constant_analytical,
    W         = NULL,
    state_tb  = NULL,
    influx        = influx,
    lambda_influx = lambda_influx_vec
  )

  epar_j$params <- list(influx = influx, lambda_influx = lambda_influx_vec)
  
  epar_j$state_tb             <- new.env()
  epar_j$state_tb$t_old       <- 0
  epar_j$state_tb$budget_used <- setNames(rep(0, length(age_names)),
                                          age_names)
  
  det_j <- run.age.deterministic(
    SIR.age.intervention.ode,
    epar      = epar_j,
    init_prop = init,
    tf        = tf,
    dt        = 0.1
  )
  
  S_cols  <- 2:(n_age + 1)
  S_start <- sum(det_j[1, S_cols])
  S_end   <- sum(det_j[nrow(det_j), S_cols])
  1 - S_end / S_start
}

# run simulations
for (sc_index in seq_along(scenarios)) {
  
  sc     <- scenarios[[sc_index]]
  N      <- sc$N
  influx <- sc$influx
  
  source(paste0(path_wd, path_experiments, path_influx,
                "/setup_parameters_exp3-6_combined.R"))
  
  pb <- txtProgressBar(min = 0, max = length(start_grid), style = 3)
  
  for (j in seq_along(start_grid)) {
  ti=  stochastic_run(j)
  
  opt_results$mean_incidence[j]   <- mean(ti)
  opt_results$sd_incidence[j]     <- sd(ti)
  opt_results$median_incidence[j] <- median(ti)
  opt_results$q25_incidence[j]    <- quantile(ti, 0.25)
  opt_results$q75_incidence[j]    <- quantile(ti, 0.75)
  opt_results$q025_incidence[j]   <- quantile(ti, 0.025)
  opt_results$q975_incidence[j]   <- quantile(ti, 0.975)
  opt_results$q95_incidence[j]    <- quantile(ti, 0.95)
  
   setTxtProgressBar(pb, j)
  }
  close(pb)
  
  
# deterministic 
  det_incidence_by_start <- sapply(start_grid, deterministic_run)
  

# Store everything for this scenario
  all_results[[sc_index]] <- list(
    N          = N,
    influx     = influx,
    label      = sc$label,
    lev        = lev,
    budget_vec = budget_vec,
    start_grid = start_grid,
    opt_results            = opt_results,
    det_incidence_by_start = det_incidence_by_start
  )
}

saveRDS(all_results,
        file = paste0(path_wd, path_results_data, final_data_file_name))


# Plot 

png(paste0(path_wd, path_results_plot, final_plot_file_name),
    width = 12, height = 9, units = "in", res = 300)

par(mfrow = c(2, 2),
    mar   = c(4.5, 4.5, 4.5, 1),   # was c(4.5, 4.5, 3, 1) -- more top room
    oma   = c(0, 0, 0, 0))         # was c(0,0,2,0) -- no outer top margin

for (sc_index in seq_along(all_results)) {
  
  r          <- all_results[[sc_index]]
  opt_res    <- r$opt_results
  det_inc    <- r$det_incidence_by_start
  start_grid <- r$start_grid
  
  # Build the panel title with N already formatted (no 1e+05)
  panel_title <- sprintf("N = %s - %s",
                         format(r$N, big.mark = ",", scientific = FALSE),
                         ifelse(r$influx, "with influx", "no influx"))
  
  plot(opt_res$start, opt_res$median_incidence, type = "b",
       pch = 16, col = "blue", lwd = 2,
       xlab = "Intervention start time (days)",
       ylab = "Final size (at t = 100)",
       ylim = c(0, 1),
       main = "")                              # title drawn manually below
  
  # Panel title (drawn higher to leave room for subheader below it)
  mtext(panel_title, side = 3, line = 2.2, cex = 1, font = 2)
  
  # Subheader: drop N, keep c_infinity and c_1 only
  mtext(bquote(paste(c[infinity] * "=" * .(r$lev) *
                       ",  " * c[1] * "=" * .(r$budget_vec))),
        side = 3, line = 0.6, cex = 0.75, col = "grey40")
  
  # IQR envelope
  polygon(c(opt_res$start, rev(opt_res$start)),
          c(opt_res$q25_incidence, rev(opt_res$q75_incidence)),
          col = rgb(0, 0, 1, 0.15), border = NA)
  
  # 95% quantile
  lines(opt_res$start, opt_res$q95_incidence,
        col = "blue", lwd = 1, lty = 4)
  
  # Deterministic curve
  lines(start_grid, det_inc, col = "red", lwd = 2, lty = 2)
  
  # Optima
  det_opt_index   <- which.min(det_inc)
  stoch_opt_index <- which.min(opt_res$median_incidence)
  q95_opt_index   <- which.min(opt_res$q95_incidence)
  
  abline(v = start_grid[det_opt_index],   col = "red",  lty = 3)
  abline(v = start_grid[stoch_opt_index], col = "blue", lty = 3)
  abline(v = opt_res$start[q95_opt_index], col = "blue", lty = 4)
  
  legend("bottomright",
         legend = c("Stoch. median, IQR",
                    "Stoch. 95% quantile",
                    "Deterministic",
                    sprintf("Det. opt: t*=%.1f",    start_grid[det_opt_index]),
                    sprintf("Median opt: t*=%.1f",  start_grid[stoch_opt_index]),
                    sprintf("Q95 opt: t*=%.1f",     opt_res$start[q95_opt_index])),
         col = c("blue", "blue", "red", "red", "blue", "blue"),
         lwd = c(2, 1, 2, 1, 1, 1),
         lty = c(1, 4, 2, 3, 3, 4),
         pch = c(16, NA, NA, NA, NA, NA),
         bty = "n", cex = 0.7)
  
  cat(sprintf("\n[%s]  det opt = %.1f (%.4f)   stoch median opt = %.1f (%.4f)\n",
              panel_title,
              start_grid[det_opt_index],   min(det_inc),
              start_grid[stoch_opt_index], min(opt_res$median_incidence)))
}

dev.off()

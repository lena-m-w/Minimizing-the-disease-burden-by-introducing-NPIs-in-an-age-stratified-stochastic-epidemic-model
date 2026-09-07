################################################################################
# Three-cell budget optimisation across a level grid (Figure 8 optimum)
# 
# Only one level is used in this run.
#
# Values in Table 3 and 4 in the Thesis were run with a exactness of 0.001 
# After that code was slightly changed, and Nelder-Mead was added and due to little time it was 
# run again with a coarser value of 0.005 as now stated in the thesis. Figures in the thesis were updated. 
# Table values were not updated.
# Corrected values for table: 
# Table 3: Cumulative fatalities for 70+: 92.7 => 92.8
#          Cumulative fatalities for Total: 159.4 => 159.5                      
#
#                         Budget    |   Start     |  Duration   |    End
# Table 4: within < 70: 5.88 => 5.93| 25.1 == 25.1| 10.0 => 10.1| 35.0 => 35.1
#          within  70+: 0.81 => 0.75|  7.8 =>  8.6| 83.3 => 77.1| 91.0 => 85.7
#          between    : 8.31 => 8.33| 10.2 => 10.1| 54.9 => 55.0| 65.1 == 65.1
#
# AI statement:
# This script was adapted from a code skeleton generated with AI Claude Opus 4.8.
# I contributed with the modelling decisions and the search and optimisation logic.
# I have carefully checked and tested all code and take responsibility for its correctness.
################################################################################

# run : run_code.R to get basic settings

# 1. Setup

# 1.1 libaries
library(parallel)

# 1.1 paths 
path_results_data <- "/03_results/results_data"
path_results_plot <- "/03_results/results_plot"
source(paste0(path_wd, path_experiments,
              "/Experiment_age/setup_parameters_6.R"))
out_dir_data <- file.path(path_wd, path_results_data)
out_dir_plot <- file.path(path_wd, path_results_plot)
stopifnot(exists("params"), exists("init_prop"),
          exists("u_func_constant_analytical_rr_age2"),
          exists("compute_burden"),
          exists("run.age.deterministic"))

# 1.2 fixed parameters

n_age    <- params$n_age
B_TOTAL  <- 15
DT     <- dt  <- 0.1 # 0.1 more exact, man köntne das vielleciht bei einem feineren machen, hier könnte ein fehler herkommen# potential error source, 0.1 more exact but takes much longer
TF_LOCAL <- params$tf
T1_LOWER <- 0
T1_UPPER <- TF_LOCAL - 5
LEVEL_GRID <- 0.75 #c(1, 0.75, 0.5)

# 1.3 run code parallel settings

run_parallel = FALSE
# Workers that can be used by this machine and files needed to run the code 
n_workers <- 1
helper_files <- c( paste0(path_wd,path_preparation,"/read_packages.R"),
                   paste0(path_wd, path_functions,"/compute_burden.R"),
                   paste0(path_wd,path_deterministic_functions,"/deterministic_SIR.R"),
                   paste0(path_wd,path_stochastic_functions,"/state_transition.R"),
                   paste0(path_wd,path_stochastic_functions,"/initiaization_initial_state.R"),
                   paste0(path_wd,path_stochastic_functions,"/intervention_functions.R"),
                   paste0(path_wd,path_stochastic_functions,"/rate_function_age.R"),
                   paste0(path_wd,path_data,"/contact_matrix.R"),
                   paste0(path_wd, path_experiments,
                          "/Experiment_age/setup_parameters_6.R")
)

################################
# Help functions
###############################

# Core evaluator: three independent budgets per cell.
evaluate_3cell <- function(c_11, c_22, c_12,
                           lev_11, lev_22, lev_12,
                           t1_11, t1_22, t1_12,
                           params, init_prop,
                           dt) {
  pi_v <- params$f
  dur_11 <- if (c_11 > 0 && lev_11 > 0) c_11 / (lev_11 *pi_v[1]* pi_v[1]) else 0
  dur_22 <- if (c_22 > 0 && lev_22 > 0) c_22 / (lev_22 * pi_v[2]* pi_v[2]) else 0
  dur_12 <- if (c_12 > 0 && lev_12 > 0) c_12 / (lev_12*2*pi_v[2]* pi_v[1]) else 0
  
  u_func_3cell <- function(state, params, t) {
    active_11 <- (t > t1_11) && (t < t1_11 + dur_11)
    active_22 <- (t > t1_22) && (t < t1_22 + dur_22)
    active_12 <- (t > t1_12) && (t < t1_12 + dur_12)
    U_mat <- matrix(0, 2, 2)
    if (active_11) U_mat[1, 1] <- lev_11
    if (active_22) U_mat[2, 2] <- lev_22
    if (active_12) { U_mat[1, 2] <- lev_12; U_mat[2, 1] <- lev_12 }
    list(one_minus_U = 1 - U_mat)
  }
  
  p <- params
  p$u_func <- u_func_3cell
  p$start_and_end_times <- c(t1_11, t1_22, t1_12, t1_11+dur_11, t1_22+dur_22, t1_12+dur_12)
  out <- tryCatch(
    run.age.deterministic(SIR.age.intervention.ode, p, init_prop,
                          tf = TF_LOCAL, dt = dt),
    error = function(e) {
      warning("run.age.deterministic failed for: ", conditionMessage(e), "\n", 
              "params: ", paste(c_11, c_22, c_12,
                                lev_11, lev_22, lev_12,
                                t1_11, t1_22, t1_12))
      NULL
    }
  )
  if (is.null(out)) return(list(F_total = 100000, res = NULL))
  res <- compute_burden(out, p)
  list(F_total = res$summary$F_total, res = res)
}
  

optimise_timings <- function(c_11, c_22, c_12, lev_11, lev_22, lev_12,
                             t1_init, t1_box) {
  obj <- function(tv){
    ev <- evaluate_3cell(c_11, c_22, c_12, lev_11, lev_22, lev_12,
                   tv[1], tv[2], tv[3], params, init_prop, dt = dt)
    if (is.na(ev$F_total)){
      return(1e6)}
    return(ev$F_total)
  }
  
  # L-BFGS-B
  fit_l <- tryCatch(
    optim(par = t1_init, fn = obj, method = "L-BFGS-B",
          lower = t1_box$lower, upper = t1_box$upper,
          control = list(factr = 1e7, maxit = 200)),
    error = function(e) {
      warning("L-BFGS-B failed: ", conditionMessage(e), "\n",
              "params: ", paste(c_11, c_22, c_12, lev_11, lev_22, lev_12))
      NULL
    }
  )
  
  # Nelder-Mead (always run it too), clamped to the box afterward
  fit_n <- tryCatch(
    optim(par = t1_init, fn = obj, method = "Nelder-Mead",
          control = list(reltol = 1e-8, maxit = 500)),
    error = function(e) {
      warning("Nelder-Mead failed: ", conditionMessage(e), "\n",
              "params: ", paste(c_11, c_22, c_12, lev_11, lev_22, lev_12))
      NULL
    }
  )
  if (!is.null(fit_n)) {
    fit_n$par   <- pmin(pmax(fit_n$par, t1_box$lower), t1_box$upper)
    fit_n$value <- obj(fit_n$par)   # re-evaluate after clamping
  }
  
  # keep the better of the two
  cands <- Filter(Negate(is.null), list(fit_l, fit_n))
  if (length(cands) == 0) {
    warning("Both optimizers failed for params: ",
            paste(c_11, c_22, c_12, lev_11, lev_22, lev_12))
    return(list(t1 = t1_init, F_total = 100000))
  }
  
  best <- cands[[which.min(sapply(cands, `[[`, "value"))]]
  return(list(t1 = best$par, F_total = best$value))
}

search_share12 <- function(share_22, share_12_grid,
                           lev_11, lev_22, lev_12,
                           t1_init, t1_box,
                           pb = NULL, cnt = NULL, pb_max = NULL) {
  
  log_df <- data.frame()
  best <- list(share_12 = NA, share_11 =NA, F_total = Inf, t1 = t1_init)
  for(s12 in share_12_grid){
    s11 <- 1-share_22-s12
   s11 <- max(s11, 0)
   c_11 <- s11 * B_TOTAL
   c_22 <- share_22 * B_TOTAL
   c_12 <- s12 * B_TOTAL
   
   # hier vielleicht anstatt t1_init best$t1 einfügen dann ist hier auch warm start
   ti <- optimise_timings(c_11, c_22, c_12, lev_11, lev_22, lev_12,
                          #t1_init, t1_box)
                          best$t1, t1_box)
   
   if (!is.null(pb)) {                                  
     cnt$n <- cnt$n + 1L
     setTxtProgressBar(pb, min(cnt$n, pb_max))
   }
   log_df <- rbind(log_df,
                   data.frame(share_22 = share_22, share_12 = s12,
                              share_11 = s11,
                              F_total = ti$F_total,
                              t1_11 = ti$t1[1], t1_22 = ti$t1[2],
                              t1_12 = ti$t1[3]))
   if (ti$F_total < best$F_total)
     best <- list(share_12 = s12, share_11 = s11,
                  F_total = ti$F_total, t1 = ti$t1)
  }
  
  list(log = log_df, best = best)
}

###################################
# run single pass
###################################
# noch ziemlich viel an dieser funktion zu basteln. 
# progress bar hinzufügen
# 
run_pass <- function(share_22_grid, share_12_step,
                     lev_11, lev_22, lev_12,
                     t1_init_first,
                     pass_name, complete_search, best1_s12 = NA,  warm_start_range = 0.1) {
 
  # initial values  (for start time and )
  log_share22 <- data.frame()
  prev_t1 <- t1_init_first
  # prev_s12  <- NULL
  
  # progressbar muss noch geöffnet werden. 
  n_s22     <- length(share_22_grid)
  first_pts <- length(seq(0, max(0, 1 - share_22_grid[1]), by = share_12_step))
  warm_pts  <- floor(0.20 / share_12_step) + 1
  est_total <- max(1L, first_pts + (n_s22 - 1) * warm_pts)
  cat(sprintf("\n[%s] ~%d evaluations\n", pass_name, est_total))
  pb  <- txtProgressBar(min = 0, max = est_total, style = 3)
  cnt <- new.env(); cnt$n <- 0L
  # search along the share between ealderly and the rest
  for (k in seq_along(share_22_grid)){
    s22 <- share_22_grid[k]
    
    # share between group 11 and 12
    max_s12 <- 1- s22
    
    if(complete_search){
    s12_grid <- seq(0, max_s12, by = share_12_step)
    }else{
      lo <- max(0, best1_s12 - warm_start_range)
      hi <- min(max_s12, best1_s12 + warm_start_range)
      s12_grid <- if (lo > hi) seq(0, max_s12, by = share_12_step) else seq(lo, hi, by = share_12_step)
    }
    # s12_grid <- if(is.null(prev_s12)){
    #   seq(0, max_s12, by = share_12_step)
    # }else{
    #   ##### noch mal überlegen, ob ich wirklich nur 0.1 möchte, ich glaube das ist nicht ganz wie ich will
    #   lo <- max(0, prev_s12 - 0.10)
    #   hi <- min(max_s12, prev_s12 + 0.10)
    #   if (lo > hi) seq(0, max_s12, by = share_12_step) else seq(lo, hi, by = share_12_step)
    #   # vielleicht noch einbauen, falls der rand das optimum ist, dass noch mal über den rand getestet wird
  
    
    # noch mal die box durchdenken
    t1_box <- list(
      lower = T1_LOWER,
      upper = T1_UPPER
    )
    
    inner <- search_share12(s22, s12_grid, lev_11, lev_22, lev_12,
                            prev_t1, t1_box,
                            pb = pb, cnt = cnt, pb_max = est_total) 
  
  log_share22 <- rbind(log_share22,
                     data.frame(share_22 = s22,
                                share_12 = inner$best$share_12,
                                share_11 = inner$best$share_11,
                                F_total  = inner$best$F_total,
                                t1_11    = inner$best$t1[1],
                                t1_22    = inner$best$t1[2],
                                t1_12    = inner$best$t1[3]))
  
  prev_t1  <- inner$best$t1
  F_history <- log_share22$F_total
  n_stop <- max(2, ceiling(0.15 * length(share_22_grid)))
  n_check <- (which.min(F_history)+ n_stop)
  if (length(F_history) >= n_check &&
      all(tail(F_history, (n_stop)) > min(F_history, na.rm = TRUE))) { # bisschen unnötig
     
     break
  }
}
close(pb)                                              # <- close bar
return(log_share22)
}


###################################
# Grid search for one level combination
###################################
search_for_levels <- function(lev_11, lev_22, lev_12,
                              checkpoint_prefix){
  t1_init <- c(23.5, 23.5, 23.5)
  
  # pass 1
  share_22_grid_p1 <- seq(0, 1, by = 0.10)
  pass1 <- run_pass(share_22_grid_p1, share_12_step = 0.10,
                    lev_11, lev_22, lev_12,
                    t1_init_first = t1_init,
                    pass_name = "P1", complete_search = TRUE)
  
  bad1 <- pass1[is.na(pass1$F_total)|pass1$F_total >1, ]
  if(nrow(bad1)>0){
    warning(sprintf("[%s] Pass 1: %d of %d rows had F_total NA or >1:\n%s",
                    checkpoint_prefix, nrow(bad1), nrow(pass1),
                    paste(capture.output(print(bad1)), collapse = "\n")))
  }
  
  pass1 <- pass1[!is.na(pass1$F_total) & pass1$F_total <= 1, ]
  
  
  if (nrow(pass1) == 0) {
    warning(sprintf("[%s] Pass 1 empty after filter (no F_total <= 1); skipping combo.",
                    checkpoint_prefix))
    return(list(pass1 = pass1, pass2 = NULL, 
                best_overall = data.frame(share_22 = NA, share_12 = NA, share_11 = NA,
                                          F_total = 100000,
                                          t1_11 = NA, t1_22 = NA, t1_12 = NA)))
  }
  
  best1 <- pass1[which.min(pass1$F_total),]
  saveRDS(list(pass1 = pass1, best1 = best1),
          paste0(checkpoint_prefix, "_pass1_fr_0409.rds"))
  
# Pass 2
  share_22_grid_p2 <- seq(max(0, best1$share_22 - 0.10),
                          min(1, best1$share_22 + 0.10), by = 0.01)
  pass2 <- run_pass(share_22_grid_p2, share_12_step = 0.005,
                    lev_11, lev_22, lev_12,
                    t1_init_first = c(best1$t1_11, best1$t1_22, best1$t1_12),
                    pass_name = "P2", complete_search = FALSE, best1_s12 = best1$share_12 , 
                    warm_start_range = 0.1)
  
  bad2 <- pass2[is.na(pass2$F_total) | pass2$F_total > 1, ]
  if (nrow(bad2) > 0) {
    warning(sprintf("[%s] Pass 2: %d of %d rows had F_total NA or >1:\n%s",
                    checkpoint_prefix, nrow(bad2), nrow(pass2),
                    paste(capture.output(print(bad2)), collapse = "\n")))
  }
  pass2 <- pass2[!is.na(pass2$F_total) & pass2$F_total <= 1,]
  
  if (nrow(pass2) == 0) {
    warning(sprintf("[%s] Pass 2 empty after filter; using Pass 1 result.",
                    checkpoint_prefix))
    return(list(pass1 = pass1, pass2 = NULL, 
                best_overall = best1, pass2_ran = NULL))
  }
  
  best2 <- pass2[which.min(pass2$F_total), ]
  saveRDS(list(pass1 = pass1, best1 = best1, pass2 = pass2, best2 = best2),
          paste0(checkpoint_prefix, "_pass2_0409.rds"))
  

  
  # summary
  
  all_logs <- rbind(pass1, pass2)
  best_overall <- all_logs[which.min(all_logs$F_total), ]
  
  list(pass1 = pass1, pass2 = pass2,
       best_overall = best_overall,
       pass2_ran = !is.null(pass2))
}

###################################
# Code to run per level combination
####################################
evaluate_level_combos <- function(i){
  # function to evaluate and save one level combination 
  # i: an integer indicating which level combination (row of level combos) should be evaluated
  
  lv <- level_combos[i,]
  cat(sprintf("[%s] combo %d/%d\n", format(Sys.time(), "%H:%M:%S"), i, n_combos)); flush.console() # output also current on Windows and Mac
  
  combo_name <- sprintf("u11_%.2f_u22_%.2f_u12_%.2f", lv$lev_11, lv$lev_22, lv$lev_12)
  path_prefix <- file.path(out_dir_data, paste0("\\Exp8_", combo_name))

  res <- tryCatch(
    search_for_levels(lv$lev_11, lv$lev_22, lv$lev_12, path_prefix),
    error = function(e) {
      warning(sprintf("[%s] combo failed entirely: %s", combo_name, conditionMessage(e)))
      list(pass1 = NULL, pass2 = NULL,
           best_overall = data.frame(share_22 = NA, share_12 = NA, share_11 = NA,
                                     F_total = 100000,
                                     t1_11 = NA, t1_22 = NA, t1_12 = NA, pass2_ran = NULL))
    }
  )
  out <- list(combo_name = combo_name,
              lev_11 = lv$lev_11, lev_22 = lv$lev_22, lev_12 = lv$lev_12,
              result = res)
  return(out)
  }


#============================
#Create level combinations
#==========================
level_combos <- expand.grid(lev_11 = LEVEL_GRID,
                            lev_22 = LEVEL_GRID,
                            lev_12 = LEVEL_GRID)
n_combos <- nrow(level_combos)
cat(sprintf("\n[GLOBAL] %d level combinations to evaluate.\n", n_combos))


#=============================
# Parallelize and run code
#=============================

t_start <- Sys.time()

# Source all helper files
for (f in helper_files) source(f)

if(run_parallel){
  # Put session into workers
  cl <- makeCluster(n_workers, type = "FORK")
  results_list <- parLapply(cl, seq_len(n_combos), evaluate_level_combos)
  stopCluster(cl)
}else{
  results_list <- lapply(seq_len(n_combos), evaluate_level_combos)
}
t_total <- as.numeric(Sys.time() - t_start, units = "hours")
cat(sprintf("\nAll %d combinations done in %.2f hours\n", n_combos, t_total))

####################
# Combine results
#######################

all_results <- setNames(
  lapply(results_list, function(r) {
    list(lev_11 = r$lev_11, lev_22 = r$lev_22, lev_12 = r$lev_12,
         result = r$result)
  }),
  sapply(results_list, `[[`, "combo_name")
)

combo_summary <- do.call(rbind, lapply(names(all_results), function(nm) {
  r <- all_results[[nm]]
  data.frame(
    combo    = nm,
    lev_11   = r$lev_11, lev_22 = r$lev_22, lev_12 = r$lev_12,
    share_22 = r$result$best_overall$share_22,
    share_12 = r$result$best_overall$share_12,
    share_11 = r$result$best_overall$share_11,
    F_total  = r$result$best_overall$F_total,
    t1_11    = r$result$best_overall$t1_11,
    t1_22    = r$result$best_overall$t1_22,
    t1_12    = r$result$best_overall$t1_12,
    pass2    = r$result$pass2_ran
  )
}))

best_combo <- combo_summary[which.min(combo_summary$F_total), ]

cat("\n==================== GLOBAL OPTIMUM ====================\n")
cat(sprintf("  Level combination: lev_11 = %.2f, lev_22 = %.2f, lev_12 = %.2f\n",
            best_combo$lev_11, best_combo$lev_22, best_combo$lev_12))
cat(sprintf("  Budget shares:     s_11 = %.4f, s_22 = %.4f, s_12 = %.4f\n",
            best_combo$share_11, best_combo$share_22, best_combo$share_12))
cat(sprintf("  Start times:       t1_11 = %.2f, t1_22 = %.2f, t1_12 = %.2f\n",
            best_combo$t1_11, best_combo$t1_22, best_combo$t1_12))
cat(sprintf("  F_total            = %.6f\n", best_combo$F_total))
cat("========================================================\n\n")

saveRDS(list(all_results = all_results,
             combo_summary = combo_summary,
             best_combo = best_combo,
             meta = list(B_TOTAL = B_TOTAL, LEVEL_GRID = LEVEL_GRID,
                         tf = TF_LOCAL, dt = DT, n_workers = n_workers,
                         total_hours = t_total)),
        file.path(out_dir_data, "\\Exp8_final_fr_0409.rds"))
cat("Saved final consolidated -> Exp8_final_fr_0409.rds\n")




# Computes group specific fatalities of disease
#
# AI used to check correctness of code, write some lines and comments and to debug

compute_burden <- function(out, params) {
  n  <- params$n_age
  Nv <- params$f *params$N
  fr <- params$ifr            # infection-fatality ratio per group, length n_age
  
  S_mat <- as.matrix(out[, paste0("S", 1:n), drop = FALSE])   # T x n
  R_mat <- as.matrix(out[, paste0("R", 1:n), drop = FALSE])
  
  # Cumulative infected per group: X_a(t) = N_a - S_a(t)
  X_mat <- sweep(-S_mat, 2, Nv, FUN = "+")
  colnames(X_mat) <- paste0("X", 1:n)
  
  # Cumulative fatalities per group: ifr_a * R_a(t)
  F_mat <- sweep(R_mat, 2, fr, FUN = "*")
  colnames(F_mat) <- paste0("F", 1:n)
  
  # Incident fatalities per group (new deaths between consecutive output rows).
  # diff() drops one row; pad with 0 at t = 0 so dimensions match the trajectory.
  dF_mat <- rbind(0, apply(F_mat, 2, diff))
  colnames(dF_mat) <- paste0("dF", 1:n)
  
  traj <- data.frame(
    time       = out$time,
    X_mat,
    F_mat,
    dF_mat,
    X_total    = rowSums(X_mat),
    final_size = rowSums(X_mat) / sum(Nv),
    F_total    = rowSums(F_mat),       # cumulative deaths
    dF_total   = rowSums(dF_mat)       # incident deaths per output step
  )
  
  final <- tail(traj, 1)
  summary <- list(
    X_by_group = setNames(as.numeric(X_mat[nrow(X_mat), ]), paste0("X", 1:n)),
    F_by_group = setNames(as.numeric(F_mat[nrow(F_mat), ]), paste0("F", 1:n)),
    final_size = final$final_size,
    F_total    = final$F_total
  )
  
  list(trajectory = traj, summary = summary)
}

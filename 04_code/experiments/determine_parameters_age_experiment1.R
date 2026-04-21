age_labels <- c("0-19", "20-64", "65+")

N   <- 100000
dur <- 20
t1  <- 20
n_age = length(age_labels) 
params_age <- list(
  n_age     = n_age,
  beta      = beta,             # 3x3 matrix (const * C_per_capita)
  gamma     = gamma,            # 0.2
  N         = N,
  R0_target = R0_target,        # 3.0 (needed for WMR)
  t1        = t1,
  t2        = t1 + dur,
  lev       = c("0-19" = 0.75, "20-64" = 0.75, "65+" = 0.75),  # per-group levels
  u_func    = u_func_constant_age,
  W         = NULL,             # NULL = full symmetric effect everywhere
  influx    = FALSE,
  lambda_vec = NULL             # not needed when influx = FALSE
)

# ── Initial state (counts, named) ─────────────────────────────────
I0 <- round(1e-4 * f * N)
I0[I0 == 0] <- 1
S0 <- round(f * N) - I0
R0_init <- rep(0, n_age)

state0 <- c(S0, I0, R0_init)
names(state0) <- c(paste0("S", 1:n_age),
                   paste0("I", 1:n_age),
                   paste0("R", 1:n_age))

# ── Run ───────────────────────────────────────────────────────────
trans <- build_transitions_age(n_age)
result_ssa <- ssa.adaptivetau(state0, trans, ratefunc.SIR.age, params_age, tf = 200)
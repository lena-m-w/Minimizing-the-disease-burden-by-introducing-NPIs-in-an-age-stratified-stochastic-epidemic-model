# simulation wrapper ausprobieren
# simulation wrapper mit intervention ausprobieren
# multi replicate simulation schreiben
# plotten

# 1. simulation wrapper ausprobieren

# params und params_adaptivetau initialisieren
# unteranderem init initialisieren
age_names <- colnames(beta)
N        <- 100000
n_age    <- length(f)
transitions <- build_transitions_age(n_age)

state0 <- init_age_state(I0_total = 0.0001, f = f, age_names = age_names)
state0 <- round(state0*N)
params_adaptive_tau <- list(
  init             = state0,
  transitions      = t(transitions),
  ratefunc         = ratefunc.SIR,
  tf               = 200,
  jacobianFunc     = NULL,
  maxTauFunc       = NULL,
  deterministic    = FALSE,
  halting          = NULL,
  relratechange    = rep(1, length(state0)),
  tl.params        = list(epsilon = 0.05),
  reportTransitions = FALSE
)

params <- list(
  n_age       = n_age,
  age_names   = age_names,
  beta        = beta,
  gamma       = gamma,
  N           = N,
  f           = f,
  R0_target   = R0_target,
  t1          = c("0-19" = 20, "20-64" = 20, "65+" = 20),
  lev         = c("0-19" = 0.75, "20-64" = 0.75, "65+" = 0.75),
  budget      = c("0-19" = 15*f[1], "20-64" = 15*f[2], "65+" = 15*f[3]),
  u_func      = u_func_constant_age,
  W           = NULL,
  influx      = FALSE,
  lambda_vec  = NULL
)


# dann simulation wrapper laufen lassen
result <- simulate.stochastic(params= params, params_adaptive_tau = params_adaptive_tau)


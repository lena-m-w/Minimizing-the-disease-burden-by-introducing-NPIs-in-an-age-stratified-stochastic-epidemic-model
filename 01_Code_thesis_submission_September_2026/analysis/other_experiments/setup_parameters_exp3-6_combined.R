# Setup parameters that are SHARED across all 4 scenarios.
# The scenario-specific parts (N, influx) are set inside the loop in the
# main experiment file before this script is sourced.
#
# Expects in the workspace before sourcing:
#   N         : population size  (e.g. 10000 or 100000)
#   influx    : logical, TRUE/FALSE

# Epidemic parameters
gamma     <- 0.2
R0_target <- 3

# Intervention parameters
lev        <- 0.75
dur        <- 20
budget_vec <- lev * dur
t1         <- 23.5

# Initial parameters
I0_total <- 0.0001


# Population parameters
age_limits <- c(0)

data(polymod)
survey    <- polymod
countries <- "Finland"

# Simulation parameters
tf         <- 365
dt         <- 0.5
start_grid <- seq(0, 35, by = 0.5)
n_rep_opt  <- 500

budget_count            <- FALSE
intervention_activation <- FALSE
last_time_tracking      <- FALSE
budget_activ_time <- budget_count + intervention_activation + last_time_tracking

# Build epidemiology parameters
ep        <- create_beta(survey = survey, countries = countries,
                         age_limits = age_limits,
                         gamma = gamma, R0_target = R0_target)
beta      <- ep$beta
f         <- ep$pi
n_age     <- ep$n_age
age_names <-ep$age_names

# Deterministic flag for adaptivetau (account for influx column if used)
if (budget_activ_time > 0) {
  deterministic <- c(rep(FALSE, n_age * 2 + as.integer(influx)),
                     rep(TRUE,  budget_activ_time))
} else {
  deterministic <- NULL
}

# Influx rate (scales with 1/N so the expected influx rate is the same)
lambda_influx_vec <- rep(1 / N, n_age)

# Initial state
init_state <- init_age_state(I0_total, f, age_names,
                             budget_count            = budget_count,
                             intervention_activation = intervention_activation,
                             last_time_tracking      = last_time_tracking) * N
init <- list(
  S = f - I0_total * f,
  I = I0_total * f,
  R = rep(0, n_age)
)

# adaptivetau setup -- transitions depend on influx flag
transitions <- build_transitions_age(n_age, influx = influx,
                                     budget_count            = budget_count,
                                     intervention_activation = intervention_activation,
                                     last_time_tracking      = last_time_tracking)
params_adaptive_tau <- list(
  init          = init_state,
  transitions   = transitions,
  ratefunc      = ratefunc.SIR,
  tf            = tf,
  jacobianFunc  = NULL, maxTauFunc = NULL,
  deterministic = deterministic, halting = NULL,
  relratechange = rep(1, length(init_state)),
  tl.params     = list(epsilon = 0.05),
  reportTransitions = FALSE
)

# Empty results frame
opt_results <- data.frame(
  start            = start_grid,
  mean_incidence   = NA,
  sd_incidence     = NA,
  median_incidence = NA,
  q25_incidence    = NA,
  q75_incidence    = NA,
  q025_incidence   = NA,
  q975_incidence   = NA,
  q95_incidence    = NA
)

cat(sprintf("\n--- Scenario: N=%s, influx=%s ---\n",
            format(N, big.mark = ","), influx))
cat(sprintf("Evaluating %d start times x %d replicates each...\n",
            length(start_grid), n_rep_opt))

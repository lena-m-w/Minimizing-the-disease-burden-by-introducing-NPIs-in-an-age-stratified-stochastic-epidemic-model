library(deSolve)
# ODE function 
# par must contain: beta, gamma, C (contact matrix), f (population fractions), n_age
SIR.age.ode <- function(t, state, par) {
  
  with(par, {
    # Unpack state into matrices: S, I, R each of length n_age
    S <- state[1:n_age]
    I <- state[(n_age + 1):(2 * n_age)]
    R <- state[(2 * n_age + 1):(3 * n_age)]
    
    # N_a = f_a (population fraction per group, constant since no demography)
    # N <- f
    
    # Force of infection: lambda_foi_a = beta * sum_b C_ab * I_b / N_b
    
    # (beta_CC * I_C + beta_CA * I_A) terms
    lambda_foi <- (beta %*% (I ))#/ N))
    
    # ODEs (screenshot model with nu = mu = l_C = 0)
    dS <- -lambda_foi * S
    dI <-  lambda_foi * S - gamma * I
    dR <-  gamma * I
    
    list(c(dS, dI, dR))
  })
}

# Run function (extends run.deterministic to age-structured case)
run.age.deterministic <- function(ode_func, epar, init_prop, tf = 200, dt = 0.1) {
  
  times <- seq(0, tf, by = dt)
  
  # Build initial state vector: (S1, S2, S3, I1, I2, I3, R1, R2, R3)
  state0 <- c(init_prop$S, init_prop$I, init_prop$R)
  
  path <- ode(state0, times, ode_func, epar, method = "lsoda")
  
  as.data.frame(path)
}

# Set initial conditions and run
# Small seed of infection in each group, proportional to group size
I0_total <- 1e-4
I0 <- I0_total * f       # distribute seed across groups by population share
S0 <- f - I0             # remaining susceptibles (proportions)
R0_init <- rep(0, n_age)  # no recovered initially

init <- list(S = S0, I = I0, R = R0_init)

epar <- list(beta = beta, gamma = gamma, C = C, f = f, n_age = n_age)

result <- run.age.deterministic(SIR.age.ode, epar, init, tf = 200, dt = 0.1)

init_simple <- list(S = sum(S0), I = sum(I0))
epar_simple <- list(beta = 0.6, gamma = 0.2)
result_simple <- run.deterministic(SIR.ode, epar_simple, init_simple, tf = 200, dt = 0.1)

# Post-process and label columns 
age_labels <- c("0-19", "20-64", "65+")
colnames(result) <- c("time",
                      paste0("S_", age_labels),
                      paste0("I_", age_labels),
                      paste0("R_", age_labels))

# Total prevalence (sum of I across groups)
result$I_total <- rowSums(result[, paste0("I_", age_labels)])

# Plot
par(mfrow = c(1, 2), mar = c(4, 4, 2, 1))

# Panel 1: Infected by age group
cols <- c("steelblue", "darkorange", "firebrick")
matplot(result$time, result[, paste0("I_", age_labels)], type = "l", lty = 1, lwd = 2,
        col = cols, xlab = "Time (days)", ylab = "Proportion infected",
        main = "Infected by age group")
legend("topright", age_labels, col = cols, lty = 1, lwd = 2, bty = "n")

# Panel 2: Total S, I, R
plot(result$time, rowSums(result[, paste0("S_", age_labels)]), type = "l", 
     lwd = 2, col = "steelblue", xlab = "Time (days)", ylab = "Proportion",
     main = "Total S, I, R", ylim = c(0, 1))
lines(result$time, result$I_total, lwd = 2, col = "firebrick")
lines(result$time, rowSums(result[, paste0("R_", age_labels)]), lwd = 2, col = "forestgreen")
legend("right", c("S", "I", "R"), col = c("steelblue", "firebrick", "forestgreen"),
       lty = 1, lwd = 2, bty = "n")


#########comparison with results simple############

par(mfrow = c(1, 3), mar = c(4, 4, 2, 1))

# Panel 1: Infected by age group + simple total
cols <- c("steelblue", "darkorange", "firebrick")
matplot(result$time, result[, paste0("I_", age_labels)], type = "l", lty = 1, lwd = 2,
        col = cols, xlab = "Time (days)", ylab = "Proportion infected",
        main = "Infected by age group", ylim = range(result_simple$I))
lines(result_simple$time, result_simple$I, lwd = 2, col = "black", lty = 2)
legend("topright", c(age_labels, "Simple (no age)"), 
       col = c(cols, "black"), lty = c(1, 1, 1, 2), lwd = 2, bty = "n")

# Panel 2: Total S, I, R (age-structured)
plot(result$time, rowSums(result[, paste0("S_", age_labels)]), type = "l", 
     lwd = 2, col = "steelblue", xlab = "Time (days)", ylab = "Proportion",
     main = "Age-structured S, I, R", ylim = c(0, 1))
lines(result$time, result$I_total, lwd = 2, col = "firebrick")
lines(result$time, rowSums(result[, paste0("R_", age_labels)]), lwd = 2, col = "forestgreen")
legend("right", c("S", "I", "R"), col = c("steelblue", "firebrick", "forestgreen"),
       lty = 1, lwd = 2, bty = "n")

# Panel 3: S, I, R (simple model)
plot(result_simple$time, result_simple$S, type = "l",
     lwd = 2, col = "steelblue", xlab = "Time (days)", ylab = "Proportion",
     main = "Simple S, I, R", ylim = c(0, 1))
lines(result_simple$time, result_simple$I, lwd = 2, col = "firebrick")
lines(result_simple$time, 1 - result_simple$S - result_simple$I, lwd = 2, col = "forestgreen")
legend("right", c("S", "I", "R"), col = c("steelblue", "firebrick", "forestgreen"),
       lty = 1, lwd = 2, bty = "n")

par(mfrow = c(1, 1), mar = c(4, 4, 2, 1))

# Plot age and simle same time
plot(result$time, rowSums(result[, paste0("S_", age_labels)]), type = "l", 
     lwd = 2, col = "steelblue", xlab = "Time (days)", ylab = "Proportion",
     main = "Age-structured S, I, R", ylim = c(0, 1))
lines(result$time, result$I_total, lwd = 2, col = "firebrick")
lines(result$time, rowSums(result[, paste0("R_", age_labels)]), lwd = 2, col = "forestgreen")
lines(result_simple$time, result_simple$I, lwd = 2, col = "firebrick")
lines(result_simple$time, 1 - result_simple$S - result_simple$I, lwd = 2, col = "forestgreen")
lines(result_simple$time, result_simple$S, type = "l",
      lwd = 2, col = "steelblue")
legend("right", c("S", "I", "R"), col = c("steelblue", "firebrick", "forestgreen"),
       lty = 1, lwd = 2, bty = "n")
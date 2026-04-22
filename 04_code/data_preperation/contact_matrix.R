data(polymod)

cm <- contact_matrix(polymod, countries = "Finland",  
                      age_limits = c(0, 20, 65),
                      symmetric = TRUE, per_capita= TRUE)

cm2 <- contact_matrix(polymod, countries = "Finland",  
                      age_limits = c(0, 20, 65),
                      symmetric = TRUE, split = TRUE) # split causes other matrix

# Calibrate beta to achieve R0 = 3
# unsing eigenvalue approach Keeling & Rohani chapter 3 page 60 Box 3.1
# same recovery rate in all groups

C_per_capita <- cm$matrix.per.capita
f <- cm$demography$proportion
R0_target <- 3
gamma     <- 0.2


# Next-generation matrix: R = -log(1-p)*1/gamma*C*N^bar/N
# Wenn C schon durch N geteilt wird beim per capita, bin mir nicht sicher, 
# ob es hier noch mal durch N geteilt werden muss
# biggest eigenvalue of R is R_0
# beta = -log(1-p)*C
# Definitions see report

R_unscaled <- ((1/gamma)*t(C_per_capita) *(f))
C_per_capita *(f[1])

rho <- max(Re(eigen(R_unscaled)$values))  # spectral radius

# R0 = (beta/gamma) * rho(C_scaled)
# => beta = R0_target * gamma / rho(C_scaled)
const <- R0_target / rho
p = 1-exp(-const)
# Verify
R <- const* R_unscaled
R0_check <- max(Re(eigen(R)$values))

beta <- const * C_per_capita
# Formula to calculate R_0 according to socialmixr 

d <- cm2$contacts 
A <- cm2$matrix
n <- cm2$demography$proportion
M2 <- (3*A*outer(d,n))

# check if biggest value is the targeted R_0
max(eigen(M2)$value)



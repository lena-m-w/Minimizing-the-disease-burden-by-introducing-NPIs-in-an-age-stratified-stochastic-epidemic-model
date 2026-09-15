# Contact matrix
#
# Calculate transmission rate and calibrate contact matrix and check for validity of calculations. 
#
# AI used to check correctness of code and to debugg

#d1 = data(polymod)
data(polymod)
create_beta <- function(survey = polymod, countries = "Finland", age_limits = c(0, 20, 65), 
                        gamma = 0.2, R0_target = 3) {


cm <- contact_matrix(survey = survey, countries = countries,  
                      age_limits = age_limits,
                      symmetric = TRUE, per_capita= TRUE)


# Calibrate beta to achieve R0 = 3
# unsing eigenvalue approach Keeling & Rohani chapter 3 page 60 Box 3.1
# same recovery rate in all groups

C_per_capita <- cm$matrix.per.capita
pi <- cm$demography$proportion
R0_target <- R0_target
gamma     <- gamma 


# Next-generation matrix: R = -log(1-p)*1/gamma*C*N^bar/N

# biggest eigenvalue of R is R_0
# beta = -log(1-p)*C
# Definitions see report

R_unscaled <- ((1/gamma)*t(C_per_capita) *(pi))
C_per_capita *(pi[1])

rho <- max(Re(eigen(R_unscaled)$values))  # spectral radius

# R0 = (const/gamma) * rho(C_per_capita*n)
# => const = R0_target * gamma / rho(C_per_capita*n)
const <- R0_target / rho
p = 1-exp(-const)
# Verify
R <- const* R_unscaled
R0_check <- max(Re(eigen(R)$values))

beta <- const * C_per_capita#/sum(cm$demography$population)

return(list(
  beta         = beta,
  gamma        = gamma,
  R0_target    = R0_target,
  C_per_capita = C_per_capita,
  pi            = pi,
  n_age        = length(pi),
  age_names    = as.character(cm$demography$age.group),
  const        = const,
  cm           = cm
))

}

#cm2 <- contact_matrix(polymod, countries = "Finland",  
#age_limits = c(0, 20, 65),
#symmetric = TRUE, split = TRUE) # split causes other matrix

# Formula to calculate R_0 according to socialmixr 
# check for validity of other functions
#d <- cm2$contacts 
#A <- cm2$matrix
#n <- cm2$demography$proportion
#M2 <- (3*A*outer(d,n))

# check if biggest value is the targeted R_0
#max(eigen(M2)$value)

#cm2$matrix*cm2$contacts*cm2$normalisation


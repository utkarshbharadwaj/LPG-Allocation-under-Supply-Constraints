# Set your working directory
setwd("/home/diazonium1109/Downloads/Assignments/Semester - II/Optimization and Numerical Methods/Project")

# Load necessary packages
library(dplyr)
library(quadprog)

# ---------------------------------------------------------
# 1. Load Baseline Demand Data
# ---------------------------------------------------------
baseline_demand <- read.csv("baseline_demands_calculated.csv")

# ---------------------------------------------------------
# 2. Define Global Constraints and Penalties
# ---------------------------------------------------------
S <- 2045.5       # Constrained national supply in TMT 
c_h <- 532        # Health penalty weight
c_e <- 604        # Economic penalty weight

# ---------------------------------------------------------
# 3. Parameter Synthesis
# ---------------------------------------------------------
# NFHS-5 "Households using clean fuel for cooking (%)"
clean_fuel_pct <- c(97.4, 68.3, 68.2, 40.1, 71.9, 58.6, 41.5, 84.8, 32.7, 40.6)

baseline_demand <- baseline_demand %>%
  mutate(
    C_i = clean_fuel_pct,
    w_i = (100 - C_i) / 100
  )

min_survival_kg <- 10      
conversion_to_TMT <- 1e-6  

# Highly vulnerable household demographics
vulnerable_households <- c(150000, 2200000, 800000, 1800000, 1450000, 15000000, 1600000, 1950000, 750000, 1700000)

# Chance-Constrained Programming (95% guarantee)
cv <- 0.15
z_score <- qnorm(1 - 0.05)

baseline_demand <- baseline_demand %>%
  mutate(
    vulnerable_households = vulnerable_households,
    mu_i = (vulnerable_households * min_survival_kg) * conversion_to_TMT,
    sigma_i = mu_i * cv,
    theta_i = mu_i + (z_score * sigma_i)
  )

# ---------------------------------------------------------
# 4. Pure Matrix Formulation for quadprog (Normalized)
# ---------------------------------------------------------
D_id <- baseline_demand$D_id
D_ic <- baseline_demand$D_ic
w_i  <- baseline_demand$w_i
theta_i <- baseline_demand$theta_i
n <- nrow(baseline_demand)

# NORMALIZED OBJECTIVE: min sum( w_i * c * (D - x)^2 / D )
# Dmat: The quadratic penalty matrix is now divided by baseline demands (D_id and D_ic)
diag_elements <- c( (2 * w_i * c_h) / D_id, (2 * w_i * c_e) / D_ic )
Dmat <- diag(diag_elements)
# Add a microscopic constant to ensure strict positive-definiteness (solver requirement)
Dmat <- Dmat + diag(1e-8, 2 * n) 

# dvec: The linear terms for the expanded squared equation
dvec <- c( rep(2 * c_h, n) * w_i, rep(2 * c_e, n) * w_i )

# Amat and bvec: The Constraints
# Constraint 1: Sum of all allocations == S
A_sum <- rep(1, 2 * n)
b_sum <- S

# Constraint 2: x_id >= theta_i
A_theta <- rbind(diag(n), matrix(0, nrow = n, ncol = n))
b_theta <- theta_i

# Constraint 3: x_id <= D_id  =>  -x_id >= -D_id
A_did <- rbind(-diag(n), matrix(0, nrow = n, ncol = n))
b_did <- -D_id

# Constraint 4: x_ic >= 0
A_xic0 <- rbind(matrix(0, nrow = n, ncol = n), diag(n))
b_xic0 <- rep(0, n)

# Constraint 5: x_ic <= 0.70 * D_ic  =>  -x_ic >= -0.70 * D_ic
A_dic <- rbind(matrix(0, nrow = n, ncol = n), -diag(n))
b_dic <- -0.70 * D_ic

# Combine all constraints into one matrix for the solver
Amat <- cbind(A_sum, A_theta, A_did, A_xic0, A_dic)
bvec <- c(b_sum, b_theta, b_did, b_xic0, b_dic)

# ---------------------------------------------------------
# 5. Execute Solver & Formatting
# ---------------------------------------------------------
# meq = 1 tells the solver that the FIRST constraint in Amat (the sum) is a strict Equality.
result <- solve.QP(Dmat, dvec, Amat, bvec, meq = 1)

# Extract optimal allocations
optimal_x <- result$solution
baseline_demand <- baseline_demand %>%
  mutate(
    Optimal_x_id = round(optimal_x[1:n], 2),
    Optimal_x_ic = round(optimal_x[(n+1):(2*n)], 2),
    Domestic_Pct_Met = round((Optimal_x_id / D_id) * 100, 2),
    Commercial_Pct_Met = round((Optimal_x_ic / D_ic) * 100, 2)
  )

# Display and Export
write.csv(baseline_demand, "final_optimized_allocations.csv", row.names = FALSE)
print("Optimization successful! Equity Normalization Applied.")
print(baseline_demand %>% select(State, Optimal_x_id, Domestic_Pct_Met, Optimal_x_ic, Commercial_Pct_Met))

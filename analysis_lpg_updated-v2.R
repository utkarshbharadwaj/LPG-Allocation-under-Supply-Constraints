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
# 2. Define Global Constraints, Reserves, and Penalties
# ---------------------------------------------------------
S <- 2045.5                 # Total constrained national supply in TMT
rho <- 0.05                 # 5% Reserve policy
R_min <- rho * S            # Reserve parameter (102.275 TMT)
S_usable <- S - R_min       # Effectively allocable supply (1943.225 TMT)

c_h <- 532                  # Health penalty weight
c_e <- 604                  # Economic penalty weight

# ---------------------------------------------------------
# 3. Parameter Synthesis
# ---------------------------------------------------------
# NFHS-5 "Households using clean fuel for cooking (%)"
clean_fuel_pct <- c(97.4, 68.3, 68.2, 40.1, 71.9, 58.6, 41.5, 84.8, 32.7, 40.6)

baseline_demand <- baseline_demand %>%
  mutate(
    C_i = clean_fuel_pct,
    # NEW: Vulnerability Floor applied using pmax()
    w_i = pmax(0.20, (100 - C_i) / 100)
  )

min_survival_kg <- 10      
conversion_to_TMT <- 1e-6  

# Highly vulnerable household demographics
vulnerable_households <- c(150000, 2200000, 800000, 1800000, 1450000, 15000000, 1600000, 1950000, 750000, 1700000)

# Chance-Constrained Programming (95% guarantee - Log-Normal Distribution)
cv <- 0.15

baseline_demand <- baseline_demand %>%
  mutate(
    vulnerable_households = vulnerable_households,
    mu_i = (vulnerable_households * min_survival_kg) * conversion_to_TMT,
    # Convert standard mean/variance to Log-Normal parameters
    sigma_sq_log = log(1 + cv^2),
    mu_log = log(mu_i) - 0.5 * sigma_sq_log,
    # Extract the 95th percentile from the Log-Normal curve
    theta_i = qlnorm(0.95, meanlog = mu_log, sdlog = sqrt(sigma_sq_log))
  ) %>%
  select(-sigma_sq_log, -mu_log) # Clean up temporary columns

# ---------------------------------------------------------
# 4. Pure Matrix Formulation for quadprog 
# ---------------------------------------------------------
D_id <- baseline_demand$D_id
D_ic <- baseline_demand$D_ic
w_i  <- baseline_demand$w_i
theta_i <- baseline_demand$theta_i
n <- nrow(baseline_demand)

# Dmat and dvec (Proportional Normalization applied)
diag_elements <- c( (2 * w_i * c_h) / D_id, (2 * w_i * c_e) / D_ic )
Dmat <- diag(diag_elements)
Dmat <- Dmat + diag(1e-8, 2 * n) 

dvec <- c( rep(2 * c_h, n) * w_i, rep(2 * c_e, n) * w_i )

# Amat and bvec: The Constraints
# Constraint 1: Sum of all allocations == S_usable 
A_sum <- rep(1, 2 * n)
b_sum <- S_usable

# Constraint 2: x_id >= theta_i (Stochastic Survival)
A_theta <- rbind(diag(n), matrix(0, nrow = n, ncol = n))
b_theta <- theta_i

# Constraint 3: x_id <= D_id => -x_id >= -D_id (No Domestic Over-allocation)
A_did <- rbind(-diag(n), matrix(0, nrow = n, ncol = n))
b_did <- -D_id

# Constraint 4: x_ic >= 0.25 * D_ic (25% Commercial Survival Floor)
A_xic_floor <- rbind(matrix(0, nrow = n, ncol = n), diag(n))
b_xic_floor <- 0.25 * D_ic

# Constraint 5: x_id >= 0.40 * D_id (40% Domestic Equity Guarantee)
A_xid_floor <- rbind(diag(n), matrix(0, nrow = n, ncol = n))
b_xid_floor <- 0.40 * D_id

# Constraint 6: x_ic <= 0.70 * D_ic => -x_ic >= -0.70 * D_ic (Commercial Cap)
A_dic <- rbind(matrix(0, nrow = n, ncol = n), -diag(n))
b_dic <- -0.70 * D_ic

# Combine all constraints
Amat <- cbind(A_sum, A_theta, A_did, A_xic_floor, A_xid_floor, A_dic)
bvec <- c(b_sum, b_theta, b_did, b_xic_floor, b_xid_floor, b_dic)

# ---------------------------------------------------------
# 5. Execute Solver & Formatting
# ---------------------------------------------------------
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
write.csv(baseline_demand, "final_allocations_v4_fully_constrained.csv", row.names = FALSE)
print("Optimization successful! Equity Guarantee, Commercial Floors, and Vulnerability Floors Applied.")
print(baseline_demand %>% select(State, Optimal_x_id, Domestic_Pct_Met, Optimal_x_ic, Commercial_Pct_Met))


# =========================================================================
# 6. VISUALIZATION AND PLOTTING FOR REPORT
# =========================================================================
# Load visualization libraries
library(ggplot2)
library(tidyr)
library(ggrepel)

# Set a consistent, clean academic theme for all plots
report_theme <- theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(face = "bold", size = 16, margin = margin(b = 10)),
    plot.subtitle = element_text(size = 12, color = "gray30", margin = margin(b = 15)),
    plot.caption = element_text(size = 10, color = "gray50", hjust = 0),
    legend.position = "bottom",
    legend.title = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank()
  )

# ---------------------------------------------------------
# Plot 1: Percentage of Demand Met by Sector (The Policy Impact)
# ---------------------------------------------------------
# This plot visualizes the 40% floor, 70% cap, and how pain is distributed.

# Reshape data for grouped bar chart
plot_data_pct <- baseline_demand %>%
  select(State, w_i, Domestic_Pct_Met, Commercial_Pct_Met) %>%
  pivot_longer(cols = c(Domestic_Pct_Met, Commercial_Pct_Met), 
               names_to = "Sector", 
               values_to = "Pct_Met") %>%
  mutate(
    Sector = ifelse(Sector == "Domestic_Pct_Met", "Domestic Sector", "Commercial Sector"),
    # Reorder states by vulnerability weight (Descending)
    State = reorder(State, w_i) 
  )

p1 <- ggplot(plot_data_pct, aes(x = State, y = Pct_Met, fill = Sector)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8), width = 0.7) +
  # Add constraint reference lines
  geom_hline(yintercept = 40, linetype = "dashed", color = "#2c3e50", size = 0.8) +
  geom_hline(yintercept = 70, linetype = "dashed", color = "#c0392b", size = 0.8) +
  geom_hline(yintercept = 100, linetype = "solid", color = "black", size = 0.5) +
  
  # Annotations for the constraint lines
  annotate("text", x = 1, y = 43, label = "40% Equity Floor", color = "#2c3e50", size = 4, hjust = 0) +
  annotate("text", x = 1, y = 73, label = "70% Commercial Cap", color = "#c0392b", size = 4, hjust = 0) +
  
  coord_flip() + # Flip coordinates for readable state names
  scale_fill_manual(values = c("Domestic Sector" = "#2980b9", "Commercial Sector" = "#e67e22")) +
  scale_y_continuous(labels = function(x) paste0(x, "%"), breaks = seq(0, 100, 20)) +
  labs(
    title = "Optimization Output: Percentage of Baseline Demand Met",
    subtitle = "Highlighting the activation of the 40% domestic equity floor and 70% commercial policy cap.",
    x = "State (Ordered by Vulnerability \u2192)",
    y = "Percentage of Normal Demand Allocated",
    caption = "Source: Optimization Solver Output | TMT Data: 2026 Projections"
  ) +
  report_theme

# Save Plot 1
ggsave("plot1_allocation_percentages.png", plot = p1, width = 10, height = 7, dpi = 300)

# ---------------------------------------------------------
# Plot 2: The Vulnerability Trade-off (Scatter Plot) - LEGEND VERSION
# ---------------------------------------------------------

p2 <- ggplot(baseline_demand, aes(x = w_i, y = Domestic_Pct_Met)) +
  # Map color to State and size to Demand
  geom_point(aes(size = D_id, color = State), alpha = 0.85) +
  
  # Trend line
  geom_smooth(method = "lm", se = FALSE, color = "gray50", linetype = "dashed") +
  
  # 40% Floor Line
  geom_hline(yintercept = 40, linetype = "dashed", color = "#2c3e50") +
  annotate("text", x = max(baseline_demand$w_i), y = 38.5, label = "40% Binding Floor", color = "#2c3e50", hjust = 1) +
  
  scale_x_continuous(breaks = seq(0, 1, 0.1)) +
  scale_y_continuous(labels = function(x) paste0(x, "%"), limits = c(35, 100)) +
  scale_size_continuous(name = "Baseline Demand (TMT)", range = c(4, 14)) +
  
  # Use a clear, distinct color palette for the 10 states
  scale_color_brewer(palette = "Paired") + 
  
  labs(
    title = "The Health Penalty Trade-off: Protection vs. Vulnerability",
    subtitle = "States with higher biomass-reversion risk (w_i) receive geometrically higher supply protection.",
    x = "Vulnerability Weight (w_i)",
    y = "Domestic Allocation Met (%)",
    caption = "Note: Bubble size represents baseline domestic demand."
  ) +
  report_theme +
  theme(
    legend.position = "right",
    legend.box = "vertical",
    # Bring the legend title back just for this plot so it's clear
    legend.title = element_text(face = "bold", size = 12) 
  ) +
  # Force the legend dots to be large enough to see the colors clearly
  guides(color = guide_legend(override.aes = list(size = 6)))

# Save Plot 2
ggsave("plot2_vulnerability_scatter_legend.png", plot = p2, width = 10, height = 7, dpi = 300)

# ---------------------------------------------------------
# Plot 3: Absolute Shortage (Macro Visualization)
# ---------------------------------------------------------
# Shows the raw physical numbers (TMT) to give the reader a sense of scale.

plot_data_abs <- baseline_demand %>%
  select(State, w_i, D_id, Optimal_x_id) %>%
  mutate(Shortfall = D_id - Optimal_x_id) %>%
  select(State, w_i, Optimal_x_id, Shortfall) %>%
  pivot_longer(cols = c(Optimal_x_id, Shortfall), 
               names_to = "Status", 
               values_to = "TMT") %>%
  mutate(
    Status = factor(Status, levels = c("Shortfall", "Optimal_x_id"), 
                    labels = c("Unmet Demand (Shortfall)", "Allocated Supply")),
    State = reorder(State, w_i)
  )

p3 <- ggplot(plot_data_abs, aes(x = State, y = TMT, fill = Status)) +
  geom_bar(stat = "identity", width = 0.7) +
  coord_flip() +
  scale_fill_manual(values = c("Allocated Supply" = "#27ae60", "Unmet Demand (Shortfall)" = "#c0392b")) +
  labs(
    title = "Physical Supply Allocation by State (Domestic)",
    subtitle = "Total stacked bars represent baseline demand; green indicates fulfilled allocations.",
    x = "",
    y = "Thousand Metric Tonnes (TMT)",
    caption = "Data: Optimization Solver Output"
  ) +
  report_theme +
  theme(legend.position = "top", legend.justification = "left")

# Save Plot 3
ggsave("plot3_absolute_shortage.png", plot = p3, width = 10, height = 7, dpi = 300)

print("Plots successfully generated and saved to your working directory!")



# ---------------------------------------------------------
# Plot 4: Residual Societal Penalty by State (Objective Function Vis)
# ---------------------------------------------------------
# Calculates the actual numerical penalty value left over after optimization
plot_data_penalty <- baseline_demand %>%
  mutate(
    # Applying the exact objective function formula
    Health_Penalty = w_i * c_h * ((D_id - Optimal_x_id)^2) / D_id,
    Economic_Penalty = w_i * c_e * ((D_ic - Optimal_x_ic)^2) / D_ic
  ) %>%
  select(State, w_i, Health_Penalty, Economic_Penalty) %>%
  pivot_longer(cols = c(Health_Penalty, Economic_Penalty), 
               names_to = "Penalty_Type", 
               values_to = "Penalty_Value") %>%
  mutate(
    Penalty_Type = ifelse(Penalty_Type == "Health_Penalty", "Health Damage (DALYs Monetized)", "Economic Damage (Commercial Loss)"),
    State = reorder(State, w_i) # Keep ordered by vulnerability
  )

p4 <- ggplot(plot_data_penalty, aes(x = State, y = Penalty_Value, fill = Penalty_Type)) +
  geom_bar(stat = "identity", width = 0.7) +
  coord_flip() +
  scale_fill_manual(values = c("Economic Damage (Commercial Loss)" = "#e74c3c", "Health Damage (DALYs Monetized)" = "#8e44ad")) +
  scale_y_continuous(labels = scales::comma) +
  labs(
    title = "Residual Societal Penalty Under Optimal Allocation",
    subtitle = "Visualizing the minimized objective function components (Millions of INR).",
    x = "",
    y = "Residual Penalty (Millions INR)",
    caption = "Lower bars indicate successful penalty minimization by the solver."
  ) +
  report_theme +
  theme(legend.position = "top", legend.title = element_blank())

# Save Plot 4
ggsave("plot4_residual_penalty.png", plot = p4, width = 10, height = 7, dpi = 300)

# ---------------------------------------------------------
# Plot 5: Macro Resource Distribution (Donut Chart)
# ---------------------------------------------------------
# Shows exactly where the 2045.5 TMT went (Domestic, Commercial, Reserve)

macro_data <- data.frame(
  Category = c("Allocated: Domestic", "Allocated: Commercial", "Strategic Reserve (Unallocated)"),
  Amount = c(sum(baseline_demand$Optimal_x_id), 
             sum(baseline_demand$Optimal_x_ic), 
             102.275) # The 5% reserve you hardcoded
)

# Calculate percentages for the labels
macro_data <- macro_data %>%
  mutate(
    Percentage = Amount / sum(Amount) * 100,
    Label = sprintf("%s\n%.1f TMT (%.1f%%)", Category, Amount, Percentage)
  )

p5 <- ggplot(macro_data, aes(x = 2, y = Amount, fill = Category)) +
  geom_bar(stat = "identity", color = "white", width = 1) +
  coord_polar(theta = "y", start = 0) +
  xlim(0.5, 2.5) + # Creates the donut hole
  scale_fill_manual(values = c("Allocated: Domestic" = "#2980b9", 
                               "Allocated: Commercial" = "#e67e22", 
                               "Strategic Reserve (Unallocated)" = "#7f8c8d")) +
  theme_void() + # Removes gridlines and axes for a clean pie chart
  labs(
    title = "National LPG Crisis Allocation Summary",
    subtitle = paste0("Total Constrained Supply: ", S, " TMT")
  ) +
  theme(
    plot.title = element_text(face = "bold", size = 16, hjust = 0.5, margin = margin(b = 10)),
    plot.subtitle = element_text(size = 12, color = "gray30", hjust = 0.5, margin = margin(b = 20)),
    legend.position = "right",
    legend.title = element_blank(),
    legend.text = element_text(size = 12)
  )

# Save Plot 5
ggsave("plot5_macro_distribution.png", plot = p5, width = 9, height = 6, dpi = 300)

print("Plots 4 and 5 generated successfully!")

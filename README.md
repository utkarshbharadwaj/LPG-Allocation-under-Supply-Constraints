# A Convex Optimization Approach to Constrained LPG Allocation

![R](https://img.shields.io/badge/Language-R-276DC3?style=flat-square&logo=r)
![Optimization](https://img.shields.io/badge/Methodology-Convex_Optimization-8E44AD?style=flat-square)
![Status](https://img.shields.io/badge/Status-Completed-27AE60?style=flat-square)

## 📌 Project Overview
This repository contains the data, R scripts, and LaTeX documentation for a data-driven convex optimization model designed to manage national liquefied petroleum gas (LPG) rationing during severe supply chain disruptions. 

Simulating a hypothetical 2026 geopolitical shock that constrains India's LPG supply to 75% of normal availability (2,045.5 TMT), this framework optimally allocates the scarce fuel across domestic and commercial sectors in 10 major geographical entities. The algorithm's primary objective is to minimize the combined societal penalties of **commercial economic loss** (measured in Gross Value Added) and **household biomass pollution** (measured in monetized DALYs).

## 🧠 Key Features & Innovations
* **Dual-Penalty Objective Function:** Standardizes health damages (c_h = 532 Million INR/TMT) and economic bankruptcies (c_e = 604 Million INR/TMT) into a single minimization space.
* **Proportional Normalization:** Eliminates the "Big State Bias" by penalizing percentage shortfalls rather than absolute tonnage, ensuring equitable distribution across states with vastly different populations.
* **Log-Normal Chance Constraints:** Utilizes a strictly positive, right-skewed Log-Normal distribution to guarantee a 95% stochastic minimum survival baseline, accounting for the "fat-tailed" risk of crisis hoarding.
* **Binding Real-World Policies:** Enforces strict political constraints, including a 40% domestic equity floor, a 70% commercial policy cap, and a protected 5% national strategic reserve.

## 🧮 Mathematical Formulation
The core allocation problem is modeled as a strictly convex Quadratic Programming (QP) problem:

min f(x) = sum [ w_i * ( c_h * D_i^d * ((D_i^d - x_i^d)/D_i^d)^2 + c_e * D_i^c * ((D_i^c - x_i^c)/D_i^c)^2 ) ]

Subject to:
1. Supply Constraint: sum (x_i^d + x_i^c) = S_usable (accounting for a 5% reserve)
2. Equity Floor: x_i^d >= 0.40 * D_i^d
3. Policy Cap: x_i^c <= 0.70 * D_i^c
4. Stochastic Survival: x_i^d >= theta_i

## 📁 Repository Structure
* analysis.R: The core R script containing the quadprog solver, data synthesis, parameter definitions, and ggplot2 visualization code.
* Data/: Directory containing the datasets used in our analysis:
  * baseline_demands_calculated.csv: Contains average domestic and commercial demands for 9 states and a 10th National average "Others"
  * lpg_consumption_monthly&statewise_2023-2025.csv: Contains monthly state wise domestic and commercial demand data from January 2023 to December 2025. This was used to calculated values in baseline_demands_calculated.csv
  * LPG_supply_2025.csv: Contains domestically produced as well as imported LPG supply values from January 2025 to December 2025. It average was used in the calculation of the value of S
  * NFHS-5 Percentage of household using solid fuel for cooking.pdf: Used to calculate vulnerability weights
  * 2019-21 India National Family Health Survey.pdf: NFHS-5 detailed report
  * final_allocations_fully_constrained.csv: Final allocation values as obtained from the model
 * Report/: Directory contains the final report as submitted:
   * project-lpg-allocation.pdf: Final report file listing our research methodology, and findings along with references
 * README.md

## 🚀 How to Run the Code
1. Clone the repository:
   git clone https://github.com/utkarshbharadwaj/LPG-Allocation-under-Supply-Constraints.git
   cd LPG-Allocation-under-Supply-Constraints

2. Install Required R Packages:
   Ensure you have an active R environment. Open your R console and run:
   install.packages(c("quadprog", "dplyr", "ggplot2", "tidyr", "ggrepel"))

3. Execute the Solver:
   Run the analysis.R script in RStudio or via the command line. The script will automatically load the CSV, compute the matrices, execute the Quadratic Programming solver, and output the final allocation CSV alongside the 5 analytical plots.

## 📊 Analytical Findings

### 1. The Policy Impact (Sector-wise Allocation)
![Allocation Percentages](plots/plot1_allocation_percentages.png)
* **Protective Floors:** The 40% equity floor actively prevented the convex objective from zeroing out highly developed states (e.g., Delhi, Tamil Nadu), proving the model's political viability.
* **The Squeezed Middle:** Mid-tier states settled in the ~55-60% allocation range, representing the pure gradient equilibrium where the marginal economic penalty perfectly balances the marginal health penalty.

### 2. Physical Supply Allocation vs. Shortfall
![Absolute Shortage](plots/plot3_absolute_shortage.png)
* **Prioritization of Vulnerability:** Highly vulnerable states (e.g., Uttar Pradesh, West Bengal) successfully retained ~80% of their domestic supply, minimizing respiratory distress.

### 3. Residual Societal Penalty
![Residual Penalty](plots/plot4_residual_penalty.png)
* Visualizing the minimized objective function components, demonstrating the exact financial and epidemiological damages remaining under the optimal allocation.

## 🎓 Authors
* **Utkarsh Bharadwaj**
* **Neerav Bhuyan**
* **Saptarshi Datta**

*Indian Statistical Institute (April 2026)*

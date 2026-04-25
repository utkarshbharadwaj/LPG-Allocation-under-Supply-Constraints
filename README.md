# A Convex Optimization Approach to Constrained LPG Allocation

![R](https://img.shields.io/badge/R-276DC3?style=for-the-badge&logo=r&logoColor=white)
![Optimization](https://img.shields.io/badge/Convex_Optimization-Quadprog-success?style=for-the-badge)
![Institution](https://img.shields.io/badge/Indian_Statistical_Institute-Delhi-red?style=for-the-badge)

## Overview
This repository contains the data, methodology, and R scripts for a data-driven convex optimization model designed to manage national liquefied petroleum gas (LPG) rationing during severe supply chain disruptions. 

Simulating a hypothetical 2026 geopolitical conflict that reduces India's LPG supply to 75% of its normal capacity (2045.5 TMT), this model optimally allocates the constrained fuel across domestic and commercial sectors across 10 geographical entities. The objective is to minimize the combined societal penalties of commercial economic loss and household biomass pollution.

## Key Mathematical Innovations
Unlike standard linear rationing policies, this model employs several advanced numerical methods to ensure equitable and realistic distribution:
* **Strictly Convex Quadratic Programming (QP):** Formulated and solved using the `quadprog` package to balance competing sector penalties.
* **Proportional Normalization:** The objective function minimizes the *percentage* shortfall `(D - x)^2 / D` rather than absolute volume. This critical adjustment eliminates "Big State Bias," ensuring highly populated states do not disproportionately consume the supply at the expense of smaller states.
* **Chance-Constrained Programming:** Human survival requirements are modeled stochastically. We implemented a 95% probabilistic guarantee of survival, mathematically transforming a normally distributed random variable into a deterministic linear constraint using Z-scores.

## Empirical Data Integration
The parameters for the optimization solver were rigorously synthesized from real-world datasets:
* **Supply & Demand constraints:** Averaged historical consumption (2023-2025) and physical supply data from the Petroleum Planning & Analysis Cell (PPAC).
* **Vulnerability Weights (w_i):** Derived from the 2019-2021 National Family Health Survey (NFHS-5) clean fuel usage indicators, mathematically prioritizing states heavily reliant on solid biomass.
* **Penalty Weights (c_h & c_e):** * *Health Penalty:* Monetized using the Global Burden of Disease (GBD) Disability-Adjusted Life Years (DALYs) and the Value of a Statistical Life (VSL) for India.
  * *Economic Penalty:* Synthesized using the Gross Value Added (GVA) for the Accommodation and Food Services sector from MoSPI.

## Repository Structure
* `analysis_lpg_updated-v2.R` : The core optimization script containing the complete matrix formulation for the `quadprog` solver.
* `baseline_demands_calculated.csv` : The aggregated empirical input dataset containing baseline demands and state-level parameters.
* `final_optimized_allocations.csv` : The solver output detailing the optimal allocation in Thousand Metric Tonnes (TMT) and percentage of normal demand met.
* `ONM_Project_Report.pdf` : The comprehensive academic report detailing the problem definition, parameter synthesis methodologies, mathematical proofs, and analytical findings.

## Prerequisites & Installation
To run the optimization model locally, you will need **R** installed on your machine. The script utilizes the `quadprog` package for matrix solving, which bypasses the system compiler requirements of heavier wrappers like `CVXR`.

    # Install required packages
    install.packages("dplyr")
    install.packages("quadprog")

## Usage
1. Clone this repository to your local machine.
2. Open `analysis_lpg_updated-v2.R` in RStudio or your preferred R environment.
3. Update the `setwd()` path at the top of the script to match your local directory.
4. Run the script. The solver will ingest the baseline CSV, synthesize the chance-constrained parameters, execute the QP optimization, and output the `final_optimized_allocations.csv` file.

## Results Summary
The model successfully dynamicized prioritization based on demographic vulnerability. In the simulated 75% supply crisis:
* Highly vulnerable, biomass-reliant states (e.g., Uttar Pradesh, West Bengal) retained ~82% to 84% of their normal domestic supply.
* Wealthy states with near-universal clean fuel access (e.g., Delhi, Tamil Nadu) were mathematically restricted to their absolute survival baselines.
* The proportional normalization ensured equitable commercial sector funding up to a strict 70% governmental policy cap for mid-tier and highly vulnerable states.

## Authors
* **Utkarsh Bharadwaj**
* **Saptarshi Datta**
* **Neerav Bhuyan**

*Indian Statistical Institute, Delhi Center*
*Optimization and Numerical Methods (Semester II)*

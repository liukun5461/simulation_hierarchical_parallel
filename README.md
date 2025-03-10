# README: Parallel Simulation of Hierarchical Data & ATE Estimation

## Overview
This R script performs a **parallel simulation** of hierarchical (multilevel) data and estimates the **Average Treatment Effect (ATE)** along with its **standard deviation (SD)** using the **Monte Carlo method**. The script is designed for efficient computation using **parallel processing** to speed up large-scale simulations.

## Methodology
The simulation follows a **multilevel model**, where data is generated hierarchically with varying group effects. The Monte Carlo method is used to:

1. Generate hierarchical data.
2. Compute treatment effects across multiple simulations.
3. Estimate the **ATE (Average Treatment Effect)**.
4. Compute the **standard deviation (SD) of ATE estimates**.

## Parallel Processing Implementation
1. Uses **foreach** with **doParallel** to distribute simulations across multiple CPU cores.
2. Efficiently aggregates results from parallel simulations.

## Author
**Kun Liu**

# Area-Level Dirichlet Mixed Models (ADMM) for SAE

![R](https://img.shields.io/badge/Language-R-276DC3) ![License](https://img.shields.io/badge/License-MIT-green)

This repository contains the **R** implementation routines developed for the research paper on **Small area estimation of proportions and rates under area-level Dirichlet mixed models**.
The code provides the necessary tools for Small Area Estimation (SAE) of proportions and rates under a frequentist approach, using the Dirichlet distribution for the target variables.

## 📂 Repository Structure

* `R/functions_MLE_Laplace.R`: Functions for model fitting via Laplace approximation (ML-Laplace).
* `R/functions_MLE_GQ.R`: Functions for fitting via Gaussian Quadrature (GHQ).
* `R/functions_predictors.R`: Functions for obtaining EBPs and plug-in predictors.
* `R/boot_Dirichlet.R`: Function for obtaining the MSEs of the plug-in predictors.
* `R/example_of_use.R`: Script with an example of use with real-world data.
* `Simulation/`: Scripts for reproducing the simulation studies presented in the manuscript.

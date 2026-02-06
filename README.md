# Area-Level Dirichlet Mixed Models (ADMM) for SAE

![R](https://img.shields.io/badge/Language-R-276DC3) ![License](https://img.shields.io/badge/License-MIT-green)

This repository contains the **R** implementation routines developed for the research paper on **Small area estimation of proportions and rates under area-level Dirichlet mixed models**.
The code provides the necessary tools for Small Area Estimation (SAE) of proportions and rates under a frequentist approach, using the Dirichlet distribution for the target variables.

## 📂 Repository Structure

* `R/functions_MLE_Laplace.R`: Functions for model fitting via Laplace approximation (ML-Laplace).
* `R/functions_MLE_GQ.R`: Functions for fitting via Gaussian Quadrature (GHQ).
* `R/functions_predictors.R`: Functions for obtaining EBPs and plug-in predictors.
* `R/boot_Dirichlet.R`: Function for obtaining the MSEs of the plug-in predictors.
* `R/RealData.R`: Script with an example of use with real-world data.
* `R/real.data.xlsx`: Excel file with the real data.
* `Simulation/`: Scripts for reproducing the simulation studies presented in the manuscript.

  
## 📊 Data Description

The repository includes the file `real.data.xlsx`, which contains the dataset used for the practical application of the model. Each row represents a domain or small area.

The variables are structured as follows:

* **Domain Identifiers:**
    * `Prov`, `Age`, `Sex`: Categorical variables defining the domains $d$ (cross-classification of Province, Age Group, and Sex).

* **Response Variable (Proportions):**
    * `Y1`, `Y2`, `Y3`: Direct estimates for the $K=3$ categories of the target variable in each domain (proportions of employed, unemployed and inactive).
    * `cont.Y1`, `cont.Y2`, `cont.Y3`: Sample counts or frequencies associated with each category.

* **Covariates (Matrix X):**
    * Variables ending in `.MEAN` (e.g., `nac1.MEAN`, `edu2.MEAN`, `work1.MEAN`): Represent the means or proportions of auxiliary variables (nationality, education, labor status) within the domain. These are used to construct the design matrix $\mathbf{X}_d$ for the fixed part of the model.

* **Design and Precision Information:**
    * `nd`: Sample size in the domain $d$
    * `Nd`: Estimated population size of the domain $d$
    * `varY*` and `covY*`: Variance and covariance estimates of the direct estimators.

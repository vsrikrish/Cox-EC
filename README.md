# Cox-EC
Analysis of the Cox et al (2018) emergent constraint for equilibrium climate sensitivity

This repository contains R code for an analysis of the sensitivity of the predicted equilibrium climate sensitivity (ECS) distribution produced by the Cox et al (2018) emergent constraint to several alternative model specifications and time window lengths.

## Required Data and Packages
For data and processing scripts, please contact the corresponding author of Cox et al (2018), Peter Cox. To run these scripts, the data for a given time window should be processed and saved into a csv (by default, in data/ and named "ECS-<windowlength>.csv") with columns <index>, name, ECS, Psi, and dPsi. The first row should contain the Psi and dPsi values for the observed data, while the remaining rows should include the ECS and Psi values for each processed CMIP5 model. 
  
To produce MCMC chains, mcmc_driver.R requires the following R packages: DEOptim, adaptMCMC, parallel. To sample from the posterior predictive distribution, sample_ECS.R requires coda and reshape2, and the plotting functions require reshape2, RColorBrewer, and ggplot2.

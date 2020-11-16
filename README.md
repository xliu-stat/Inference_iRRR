# Inference_iRRR
This repository contains code for the paper 'Multivariate Log-Contrast Regression with Sub-Compositional Predictors: Testing the Association Between Preterm Infants’ Gut Microbiome and Neurobehavioral Outcomes.' All the functions used in the paper and the code for simulation (part 1 and part 2) are included. As the dataset used in the application section is not available for public distribution, this part of code (i.e., application and simulation part 3) is not included. 

R file descriptions:

* `functions_inference.R`: contains all the functions used in the paper
* `Simu_part1_notADMM_multivariate.R`: simulation with normally distributed predictors. Inference results are obtained from multivariate analysis and the scaled iRRR is fitted with blockwise coordinate descent algorithm. (Setting 1 of simulation part 1 is fitted with blockwise coordinate descent algorithm)
* `Simu_part1_ADMM_multivariate.R`: simulation with normally distributed predictors. Inference results are obtained from multivariate analysis and the scaled iRRR is fitted with ADMM algorithm. (Setting 2 of simulation part 1 is fitted with ADMM)
* `Simu_part1_ADMM_univariate.R`:  simulation with normally distributed predictors. Inference results are obtained from univariate analysis and the scaled iRRR is fitted with ADMM algorithm.
* `Simu_part2_ADMM_multivariate.R`: simulation with generated compositional predictors. Inference results are obtained from multivariate analysis and the scaled iRRR is fitted with ADMM algorithm.
* `Simu_part2_ADMM_univariate.R`: simulation with generated compositional predictors. Inference results are obtained from univariate analysis and the scaled iRRR is fitted with ADMM algorithm.




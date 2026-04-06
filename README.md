R package for "Adjusting for selection bias in effect estimation from nonprobability samples: A latent heterogeneity‑adjusted approach with application to health monitoring survey data"

Implement functions to estimate association parameters of interest for nonprobability sample.

**Installation**
```
source(main_function_NHEI.R)
source(simu_data_NHEI.R)
```

**Example for application based on simuated data**  

Estimate the effect size and 95%CI of with simulated data using the NHEI method.   

1.Generate the simulated data using `simu_data_NHEI()`
```
library(NHDR)
data <- simu_data_NHEI(N=100000,n_A=500,n_B=2000,n_LC=2,icc=0=0.5,n_Z=1,n_P=1,dist=1,path_scenario=9,R2_ZtoR=0.5,R2_ZtoX=0.5,R2_ZtoY=0.5,coef_ZXtoY=0.1,coef_XtoR=0.1,coef_YtoR=0.1)
data_SA <- data$nonprobability_sample
data_SB <- data$reference_sample
```
2.Estimate the effect size and the 95%CI, with the upper limit of latent classes set as 5
```
NHEI <- NPasso_ALC(data_NP=data_SA,data_P=data_SB,
                    CIV=c('Z1','Z2','Z3','Z4','Z5'),
                    type_CIV=c('binomial','gaussian','gaussian','gaussian','gaussian'),
                    sample_weight='sample_weight',maxclassnum=5)

weight <- NHEI$kwweight_dorm


```
**Example for simulation study**  

1.Load the `simucode`
```
source(simulation.R)
```

2.Set the simulation parameters and run the simulation. For example, take the following as the illustrative parameter setup and the code shown below:
```
simucode(
  B = 3000,           # number of replications
  N = 100000,         # finite population size
  n_A = 500,          # size of the nonprobability sample
  n_B = 2000,         # size of the reference sample
  R2 = 0.5,           # correlation parameter between key variables
  n_LC = 2,           # number of latent classes
  icc = 0.5,          # intraclass correlation coefficient
  n_Z = 1,            # number of covariates with sign differences 
  n_P = 1,            # number of subpopulations with opposing covariate effects
  Ytype = 'gaussian', # outcome distribution
  dist = 1,           # how Z_i distributions vary across subpopulations
  maxclassnum = 4     # upper limit for latent classes
)

R package for "Adjusting for selection bias in effect size estimation from nonprobability samples: A latent heterogeneity‑adjusted approach with application to health monitoring survey data"

Implement functions to estimate effect size parameters of interest for nonprobability sample.

**Installation**
```
source(main_function_NHEI.R)
source(simu_data_NHEI.R)
```

**Example for application based on simuated data**  

Estimate the effect size and 95%CI of with simulated data using the NHEI method.   

1.Generate the simulated data using `simu_data_NHEI()`
```
data <- simu_data_NHEI(N=100000,n_A=500,n_B=2000,n_LC=2,icc=0.5,n_Z=1,n_P=1,Ytype='gaussian',dist=1,path_scenario=9,R2_ZtoR=0.5,R2_ZtoX=0.5,R2_ZtoY=0.5,coef_ZXtoY=0.1,coef_XtoR=0.1,coef_YtoR=0.1)
data_SA <- data$nonprobability_sample
data_SB <- data$reference_sample
```
2.Estimate the weights for adjusting selection bias, with the upper limit of latent classes set as 5 in the proposed NHEI method
```
weight <- NPasso_ALC(data_NP=data_SA,data_P=data_SB,
                   CIV=c('Z1','Z2','Z3','Z4','Z5'),
                   type_CIV=c('binomial','gaussian','gaussian','gaussian','gaussian'),
                   sample_weight='sample_weight',maxclassnum=5)
```
3.Fit the weighted effect size model
```
fit_model <- glm(Y~X+Z1+Z2+Z3+Z4+Z5, data=data_SA, family = gaussian, weights = weight)
```
4.Extract point estimate of the adjusted effect size
```
point_estimate <- round(coef(fit_model)[2],3)
```
5.Estimate variance of the point estimate via Huber Sandwich Estimator
```
vcov <- vcovHC(fit_model)
se <- (diag(vcov)[2])^0.5
```
6.Compute the 95% confidence interval based the above estimate of variance
```
lower <- round(coef(fit_model)[2]-1.96*se,3)
upper <- round(coef(fit_model)[2]+1.96*se,3)
output <- list('point_estimate'=point_estimate,'95% CI'=paste('(',lower,", ",upper,')',sep=''))
```

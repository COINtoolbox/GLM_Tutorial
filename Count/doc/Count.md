================================================================================
GLM III - Bayesian Negative Binomial Regression and Globular Cluster Populations
================================================================================

```
GLM.pois<-model{
#Priors for regression coefficients

beta.0~dnorm(0,0.000001)
beta.1~dnorm(0,0.000001)

#Poisson GLM Likelihood

for (i in 1:N){
eta[i]<-beta.0+beta.1*x[i]
log(mu[i])<-eta[i]
y[i]~dpois(mu[i])
              }
}
```

![Example figure] 
(https://github.com/COINtoolbox/GLM_Tutorial/blob/master/Count/figures/logit3D.png)

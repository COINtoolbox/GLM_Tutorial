
# GLM III - Bayesian Negative Binomial Regression and Globular Cluster Populations
================================================================================


## Poisson 

The basic JAGS syntax  for a Poisson GLM model:

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

## Negative Binomial

The basic JAGS syntax  for a NB GLM model:

```
GLM.NB<-model{
#Priors for regression coefficients

beta.0~dnorm(0,0.000001)
beta.1~dnorm(0,0.000001)
k~dunif(0.001,10)

#NB GLM Likelihood

for (i in 1:N){
eta[i]<-beta.0+beta.1*x[i]
log(mu[i])<-eta[i]
p[i]<-k/(k+mu[i])
y[i]~dnegbin(p[i],k)
              }
}
```

Another approach to fit a NB model in JAGS is via a combination of a Gamma distribution with a Poisson distribution in the form:

```
GLM.NB<-model{
#Priors for regression coefficients

beta.0~dnorm(0,0.000001)
beta.1~dnorm(0,0.000001)
k~dunif(0.001,10)

#NB GLM Likelihood

for (i in 1:N){
eta[i]<-beta.0+beta.1*x[i]
log(mu[i])<-eta[i]
rateParm[i]<-k/mu[i]
g[i]~dgamma(k,rateParm[i])
y[i]~dpois(g[i],k)
              }
}
```

![Example figure] 
(https://github.com/COINtoolbox/GLM_Tutorial/blob/master/Count/figures/logit3D.png)

## Bibtex entry

If you use this tutorial in your research, we kindly as you to cite the original paper:

[de Souza, R. S.,  *et al.*,  The overlooked potential of generalized linear models in astronomy - III. Bayesian negative binomial regression and globular cluster populations, MNRAS, vol. 453, p.1928-1940](http://adsabs.harvard.edu/abs/2015MNRAS.453.1928D)

The corresponding bibitex entry is:

```

@ARTICLE{2015MNRAS.453.1928D,
   author = {{de Souza}, R.~S. and {Hilbe}, J.~M. and {Buelens}, B. and {Riggs}, J.~D. and 
	{Cameron}, E. and {Ishida}, E.~E.~O. and {Chies-Santos}, A.~L. and 
	{Killedar}, M.},
    title = "{The overlooked potential of generalized linear models in astronomy - III. Bayesian negative binomial regression and globular cluster populations}",
  journal = {\mnras},
archivePrefix = "arXiv",
   eprint = {1506.04792},
 primaryClass = "astro-ph.IM",
 keywords = {methods: data analysis, methods: statistical, globular clusters: general},
     year = 2015,
    month = oct,
   volume = 453,
    pages = {1928-1940},
      doi = {10.1093/mnras/stv1825},
   adsurl = {http://adsabs.harvard.edu/abs/2015MNRAS.453.1928D},
  adsnote = {Provided by the SAO/NASA Astrophysics Data System}
}
```



# Generalized Linear Models in Astronomy-Binomial Regression
[![arxiv](http://img.shields.io/badge/arXiv-1503.07736-lightgrey.svg?style=plastic)](http://arxiv.org/abs/1409.7696)


We elucidate the potential of a particular class of GLMs for handling binary/binomial data, the so-called logit and probit regression techniques, from both a maximum likelihood and a Bayesian perspective. As a case in point, we present the use of these GLMs to explore the conditions of star formation activity and metal enrichment in primordial minihaloes from cosmological hydro-simulations including detailed chemistry, gas physics, and stellar feedback.

## Install R and Rstudio from 

* http://www.r-project.org
* http://www.rstudio.com

## Bayesian GLM with probit link 
```{r,results='hide',message=FALSE, cache=FALSE}
  library(arm) 
  #Output identical to ML logit
   blr1 <- bayesglm(y ~ x1+x2+ ..., 
   family=binomial(link="logit"),
   prior.scale=Inf, prior.df=Inf,
   data=<datafile>)
   display(blr1)          
                
  #Bayes GLM  with default binomial
  #logit link and  Cauchy prior
  #with scale=2.5
   
   blr2 <- bayesglm(y~x1+x2+...,          
   family=binomial,
   data=<datafile>)
   display(blr2)
   
  #Bayes logit with normal prior 
  #with scale=2.5                    
   blr3 <- bayesglm(y~x1+x2+..., 
   family=binomial,
   prior.scale=2.5, prior.df=Inf,
   data=<datafile>)
   display(blr3)
```
## Bayesian GLM with logit link 
```{r,results='hide',message=FALSE, cache=FALSE}
 library(arm)
    bpr <- bayesglm(y~x1+x2+...,
    family=binomial(link="probit"), 
    prior.scale=2.5, prior.df=Inf,
    data=<datafile>)
    display(bpr)
```

![Logit 3D] 
(https://github.com/COINtoolbox/GLM_Tutorial/blob/master/Logit/figures/logit3D.png)



## Bibtex entry

If you use this tutorial in your research, we kindly as you to cite the original paper:

[de Souza, R. S.,  *et al.*,  The overlooked potential of Generalized Linear Models in astronomy, I: Binomial regression, A&C, vol. 12, p.21-32, 2015](http://adsabs.harvard.edu/abs/2015A%26C....12...21D)

The corresponding bibitex entry is:

```

@ARTICLE{2015A&C....12...21D,
   author = {{de Souza}, R.~S. and {Cameron}, E. and {Killedar}, M. and {Hilbe}, J. and 
	{Vilalta}, R. and {Maio}, U. and {Biffi}, V. and {Ciardi}, B. and 
	{Riggs}, J.~D.},
    title = "{The overlooked potential of Generalized Linear Models in astronomy, I: Binomial regression}",
  journal = {Astronomy and Computing},
archivePrefix = "arXiv",
   eprint = {1409.7696},
 primaryClass = "astro-ph.IM",
 keywords = {Cosmology, First stars, Methods, Statistical, Stars, Population III},
     year = 2015,
    month = sep,
   volume = 12,
    pages = {21-32},
      doi = {10.1016/j.ascom.2015.04.002},
   adsurl = {http://adsabs.harvard.edu/abs/2015A%26C....12...21D},
  adsnote = {Provided by the SAO/NASA Astrophysics Data System}
}

```





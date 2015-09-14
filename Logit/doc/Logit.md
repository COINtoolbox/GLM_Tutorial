# Generalized Linear Models in Astronomy
[![arxiv](http://img.shields.io/badge/arXiv-1503.07736-lightgrey.svg?style=plastic)](http://arxiv.org/abs/1503.07736)
[![ascl](http://img.shields.io/badge/ascl-1503.006-blue.svg?style=plastic)](http://ascl.net/1503.006)

Welcome to the AMADA - Analysis of Muldimensional Astronomical DAtasets 

AMADA allows an iterative exploration and information retrieval of high-dimensional data sets.
This is done by performing a hierarchical clustering analysis for different choices of correlation matrices and by doing a principal components analysis
in the original data. Additionally, AMADA provides a set of modern  visualization data-mining diagnostics.  The user can switch between them using the different tabs. 

## Install R and Rstudio from 

http://www.r-project.org
http://www.rstudio.com


## Install Required libraries
```{r,results='hide',message=FALSE, cache=FALSE}
install.packages('ape',dependencies=TRUE)
```
![Example figure] 
(https://github.com/COINtoolbox/GLM_Tutorial/blob/master/Logit/figures/logit3D.png)



## Install AMADA R package from github
```{r,results='hide',message=FALSE, cache=FALSE}
require(devtools)

install_github("RafaelSdeSouza/AMADA")
```

## Run Shiny App
```{r,results='hide',message=FALSE, cache=FALSE}

require(shiny)
runUrl('https://github.com/RafaelSdeSouza/AMADA_shiny/archive/master.zip')
### If the above does not work, try this
### options("download.file.extra" = "--no-check-certificate") 
```



###  Data Input
  
AMADA allows the users to either use available datasets or upload their own.   Check the bottom of the 'Import Dataset' panel to see if the data have been properly imported. The data can be seen on the screen by clicking in the tab "Dataset" on the main page. 

#### Available datasets

The available  datasets  follow the same nomenclature of their respective source articles. I recommend  the user to check the original articles or catalogs for a better understanding of their meaning.

* Supernova host properties ([Table 1, Sako, M. et al. 2014](http://adsabs.harvard.edu/abs/2014arXiv1401.3317S)): Type Ia and II  Supernova host-galaxy  properties  from  Sloan Digital Sky Survey  multi-band photometry.

* Mock galaxy catalog ([Guo et al. 2011](http://adsabs.harvard.edu/abs/2011MNRAS.413..101G)): Galaxy semi-analytic  formation models build on top of  the Millennium  ([Springel et al. 2005](http://adsabs.harvard.edu/abs/2003MNRAS.339..312S)) and Millenium II simulations ([Boylan-Kolchin et al. 2009](http://adsabs.harvard.edu/abs/2009MNRAS.398.1150B)). 

* [ZENS catalog](http://www.astro.ethz.ch/carollo/research/ZENS) ([Carollo et al. 2012](http://arxiv.org/abs/1206.5807), [Cibinel et al. 2012](http://arxiv.org/abs/1206.6108), [Cibinel et al. 2013](http://adsabs.harvard.edu/cgi-bin/bib_query?arXiv:1206.6496)): The Zurich ENvironmental Study (ZENS) is a survey of galaxy groups in the local universe.  The  sample consists of 141 groups in the  redshift range 0.05 < z < 0.0585.
    


#### Import dataset

 Data must be imported as a CSV/TXT format, columns are named and  separated by  spaces.
It may contain an arbitrary number of columns and rows. If missing data is present, it should be marked as NA. An example of how a dataset should be formatted can be found by clicking the tab "Dataset" on the main page.

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





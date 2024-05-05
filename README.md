
## **reReg**

[![Project Status: Active – The project has reached a stable, usable
state and is being actively
developed.](https://www.repostatus.org/badges/latest/active.svg)](https://www.repostatus.org/#active)
[![minimal R
version](https://img.shields.io/badge/R%3E%3D-3.4.0-6666ff.svg)](https://cran.r-project.org/)
[![CRAN_Status_Badge](https://www.r-pkg.org/badges/version/reReg)](https://cran.r-project.org/package=reReg)
[![packageversion](https://img.shields.io/badge/Package%20version-1.4.4-orange.svg?style=flat-square)](commits/master)
[![AppVeyor Build
Status](https://ci.appveyor.com/api/projects/status/github/stc04003/reReg?branch=master&svg=true)](https://ci.appveyor.com/project/stc04003/reReg)
[![Travis-CI Build
Status](https://travis-ci.org/stc04003/reReg.svg?branch=master)](https://travis-ci.org/stc04003/reReg)
[![Last-changedate](https://img.shields.io/badge/last%20change-2024--05--04-yellowgreen.svg)](/commits/master)
<!-- [![Build Status](https://travis-ci.org/user/pkg.svg?branch=master)](https://travis-ci.org/user/pkg) -->
<!-- README.md is generated from README.Rmd. Please edit that file -->

### Regression models for recurrent event data

***reReg*** implements a collection of regression models for recurrent
event process and failure time. The package is still under active
development.

### Installation

You can install and load **reReg** from CRAN using

``` r
install.packages("reReg")
library(reReg)
```

You can install reReg from github with:

``` r
## install.packages("devtools")
devtools::install_github("stc04003/reReg", ref = "main")
```

### Citation

Cite ***reReg*** with `citation("reReg")`.

``` r
citation("reReg")
#> To cite reReg in publications use:
#> 
#>   Chiou SH, Xu G, Yan J, Huang C (2023). "Regression Modeling for
#>   Recurrent Events Possibly with an Informative Terminal Event Using R
#>   Package reReg." _Journal of Statistical Software_, *105*(5), 1-34.
#>   doi:10.18637/jss.v105.i05 <https://doi.org/10.18637/jss.v105.i05>.
#> 
#> A BibTeX entry for LaTeX users is
#> 
#>   @Article{,
#>     title = {Regression Modeling for Recurrent Events Possibly with an Informative Terminal Event Using {R} Package {reReg}},
#>     author = {Sy Han Chiou and Gongjun Xu and Jun Yan and Chiung-Yu Huang},
#>     journal = {Journal of Statistical Software},
#>     year = {2023},
#>     volume = {105},
#>     number = {5},
#>     pages = {1--34},
#>     doi = {10.18637/jss.v105.i05},
#>   }
```

### Online documentation

[Online document](https://www.sychiou.com/reReg/index.html) includes:

- Package vignette on [Creating `Recur`
  objects](https://www.sychiou.com/reReg/articles/reReg-Recur.html).
- Package vignette on [visualization of recurrent event
  data](https://www.sychiou.com/reReg/articles/reReg-plot.html).
- Package vignette on [simulating recurrent event
  data](https://www.sychiou.com/reReg/articles/reReg-sims.html).
- Package vignette on [regression analysis for recurrent event
  data](https://www.sychiou.com/reReg/articles/reReg-reg.html).
- Package vignette on [improving semiparametric estimation of the
  Cox-type proportional rate
  model](https://www.sychiou.com/reReg/articles/reReg-cppl.html).

### References:

Lin, D., Wei, L., Yang, I. and Ying, Z. (2000). Semiparametric
Regression for the Mean and Rate Functions of Recurrent Events. *Journal
of the Royal Statistical Society: Series B (Methodological)*, **62**:
711-730.

Wang, M.-C., Qin, J., and Chiang, C.-T. (2001). Analyzing Recurrent
Event Data with Informative Censoring. *Journal of the American
Statistical Association* **96**(455): 1057-1065.

Ghosh, D. and D.Y. Lin (2002). Marginal Regression Models for Recurrent
and Terminal Events. *Statistica Sinica*, 663-688.

Ghosh, D. and D.Y. Lin (2003). Semiparametric Analysis of Recurrent
Events Data in the Presence of Dependent Censoring. *Biometrics*,
**59**: 877-885.

Huang, C.-Y. and Wang, M.-C. (2004). Joint Modeling and Estimation for
Recurrent Event Processes and Failure Time Data. *Journal of the
American Statistical Association* **99**(468), 1153-1165.

Xu, G., Chiou, S.H., Huang, C.-Y., Wang, M.-C. and Yan, J. (2017). Joint
Scale-change Models for Recurrent Events and Failure Time. *Journal of
the American Statistical Association* **112**(518): 796-805.

Xu, G., Chiou, S.H., Yan, J., Marr, K., and Huang, C.-Y. (2020)
Generalized Scale-Change Models for Recurrent Event Processes under
Informative Censoring. *Statistica Sinica* **30**, 1773–1795.

Huang, M.-Y. and Huang, C.Y. (2022) Improved semiparametric estimation
of the proportional rate model with recurrent event data. **Under
revision**.

Chiou, S.H., Xu, G., Yan, J. and Huang, C.-Y. (2022) Regression modeling
for recurrent events possibly with an informative terminal event using R
package reReg **Provisionally accepted**

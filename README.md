[![Travis Build Status](https://travis-ci.org/Pomona.svg?branch=master)](https://travis-ci.org/Pomona)
[![AppVeyor Build Status](https://ci.appveyor.com/api/projects/status/github/Pomona?branch=master&svg=true)](https://ci.appveyor.com/project/Ponoma)
[![Coverage Status](https://coveralls.io/repos/github/Pomona/badge.svg?branch=master)](https://coveralls.io/github/Pomona?branch=master)
## Pomona
Silke Szymczak and Cesaire J.K. Fouodo

### Introduction
This package provides different methods for identifying relevant
    variables in omics data sets using Random Forests. It implements the following
    approaches: empirical and parametric permutation (Altmann), Boruta, Vita,
    r2VIM (recurrent relative veriable importance), RFE (recursive feature
    elimination) and Hybrid, combining Vita and Boruta. All approaches use unscaled permutation variable importance and
    the R package ranger to generate the forests. The
    package also includes a function to simulate correlated gene expression data.

### Installation
Installation from Github:
```R
devtools::install_github("silkeszy/Pomona")
```

CRAN release coming soon.

### Usage
For usage in R, see ?Pomona in R. Most importantly, see the Examples section. As a first example you could try 

```R  
data <- simulation.data.cor(no.samples = 100, group.size = rep(10, 6), no.var.total = 200)
res <- var.sel.hybrid(x = data[, -1], y = data[, 1])
```

### References
* Nembrini, S., Koenig, I. R. & Wright, M. N. (2018). The revival of the Gini Importance? Bioinformatics. https://doi.org/10.1093/bioinformatics/bty373.
* Janitza, S, Celik, E, Boulesteix, AL. (2018). A computationally fast variable importance test for random forests for high-dimensional data. Adv Data Anal Classif.; doi.org: 10.1007/s11634-016-0276-4
* Kursa, M. B. and Rudnicki, W. R. (2010). Feature Selection with the Boruta Package. Journal of Statistical Software. \emph{Journal of Statistical Software, 36(11)}, p. 1-13. URL: \url{http://www.jstatsoft.org/v36/i11/}.
* Wright, M. N. and Ziegler, A. (2017). ranger: A fast implementation of random forests for high dimensional data in C++ and R. Journal of Statistical Software, 77(1), 1–17.
* Szymczak, S., Holzinger, E., Dasgupta, A., Malley, J. D., Molloy, A. M., Mills, J. L., Brody, L. C., Stambolian, D., and Bailey-Wilson, J. E. (2016). r2VIM: A new variable selection method for random forests in genome-wide association studies. BioData Mining, 9(1), 7.

## Introduction about nCV package
Robust oscillation of clock genes is a core feature of the circadian system. Relative amplitude (rAMP) measures the robustness of clock gene oscillations, but only works for longitudinal samples. We lack a method for estimating robust oscillations from human samples without labeled time. Therefore, we developed the nCV package for addressing this challenge. 

The nCV package has two functions: nCVnet and nCVgene. nCVnet can test whether there is a functional clock network in population scale data. nCVgene can evaluate the robustness of clock genes in population scale data. 

## Installation
Use **devtools** to install this version from Github:

  ```r
# install 'devtools' in R(>3.0.2)
# install.packages("devtools")
# install nCV
devtools::install_github('gangwug/nCV')
```

## Usage
```r
library(nCV)
# see detail introduction and associated application examples
?nCVnet
?nCVgene
```

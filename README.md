[![License](https://img.shields.io/badge/license-GPL%20%28%3E=%203%29-lightgrey.svg?style=flat)](http://www.gnu.org/licenses/gpl-3.0.html)
<img align="right" src="cenith.png" alt="drawing" width="350">


## Installation
Install via devtools from github:

``` r
devtools::install_github("SchoenbergA/CENITH")
```

## Workflow TreeCrown Segmentation

- load a canopy height model (chm) for AOI
- clip a representativ subarea
- use QGis or equal to generate a Pointlayer and set points on Treepositions
- load the PointsLayer 

``` r
require(CENITH
chm <- #load chm
tp <- #load Points with Treepositions 

define sequences to iterate over to find best fitting parameters to compute Crownareas
seqa <- seq(0.1,0.8,0.05)
seqb <- seq(0.4,0.8,0.1)
seqh <- c(10,13) # maxheight

# run iterations
x <- BestSegVal(chm,seqa,seqb,seqh,vp)
```


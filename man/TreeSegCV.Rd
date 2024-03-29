% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/TreeSegCV.R
\name{TreeSegCV}
\alias{TreeSegCV}
\title{N-Fold Cross-validation for TreeCrown Segments}
\usage{
TreeSegCV(sites, a, b, h, MIN, MAX, CHMfilter = 1, vps)
}
\arguments{
\item{sites}{list - a list of chm RasterLayers (see details)}

\item{a}{numeric - single value for MovingWindow.}

\item{b}{numeric - single value for MovingWindow.}

\item{h}{numeric - maximum height of trees (in meter) to detect trees.}

\item{MIN}{numeric - the minimum area for crowns. Smaller poylgons are cropped.}

\item{MAX}{numeric - the maximum area for crowns. Larger polygons are cropped.}

\item{CHMfilter}{numeric - uses a sum filter on the chm with a MovingWindow of (x*x), which must be odd. Default=1 no filter.}

\item{vps}{list - a list of PointLayers with estimated positions of trees (see details)}
}
\value{
returns a table with quality values for each fold (site) and calculated overall performance (mean values) for all sites.
}
\description{
performs an n-fold cross validation to estimate the performance of a segmentation 'model' for an AOI
}
\details{
\itemize{
\item 'sites' and 'vps' must be lists with the same order of chm and the respective validation points.
\item parameters of the model: a, b, h, MIN, MAX should be calculated by \code{\link{BestSegVal}}.
\item The overall performance is the mean of all sites.
\item For the details of the table values see \code{\link{BestSegVal}}.
}
}
\note{
The overall performance helps to estimate the precision for an AOI but does NOT give the "exact" precision. More folds will increase the
expressiveness but will need more time to set validation points. Further the supervised setting of validation points is highly subjectiv and does not have to correlate with the real amount and or position of trees.
}
\examples{
require(CENITH)
require(mapview)
require(raster)
require(rgdal)
# load data
chm  <- raster::raster(system.file("extdata","lau_chm.tif",package = "CENITH"))
chm2 <- raster::raster(system.file("extdata","lau_chm_side2.tif",package = "CENITH"))
chm3 <- raster::raster(system.file("extdata","lau_chm_side3.tif",package = "CENITH"))
vp <- rgdal::readOGR(system.file("extdata","lau_vp.shp",package = "CENITH"))
vp2 <- rgdal::readOGR(system.file("extdata","lau_vp_side2.shp",package = "CENITH"))
vp3 <- rgdal::readOGR(system.file("extdata","lau_vp_side3.shp",package = "CENITH"))
# handle CRS string
crs(vp) <-crs(chm)
crs(vp2)<-crs(chm)
crs(vp3)<-crs(chm)

# list all sites and validation points
chmlist <- list(chm,chm2,chm3)
vplist <- list(vp,vp2,vp3)

# run 3 fold cross validation with parameters computed by 'BestSegVal' (from example)
cv <- CENITH::TreeSegCV(sites=chmlist,a=0.3,b=0.5,h=0.5,MIN=5,MAX=1000,CHMfilter=3,vps=vplist)
cv
### the model trained with BestSegVal on site 1 reaches a overall performance of 0.77 @ 0.12 for all tree sites.
### Note that the performance on sites 2 and 3 is even better than on site 1 where it was trained. This effect is caused probably because sites 2 and 3 are more homogenious.
### For sure this is just an example and it can happen that the performance on one site is worst. The overall mean performance is used to estimate the quality of an segmentation for an AOI by testing some subareas.
}
\seealso{
\code{\link{BestSegVal}}
}
\author{
Andreas Schönberg
}

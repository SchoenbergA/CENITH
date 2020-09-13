#' CENITH Tree Segmentation Tool
#'
#' @name LEGION - package
#' @docType package
#' @title CENITH Tree Segmentation Tool
#' @description In general the \code{CENITH} package is a wrapper for the \code{ForestTools} package (Andrew Plowright and Jean-Romain Roussel,2020).
#' \code{CENITH} is designed to test for best fitting parameters to performe a Tree Segmentation. Furher it provides a cross-validation for tree segmentation.
#'
#' @note For a more independent and stable package the source code of some fucntions provided by \code{ForestTools} is used instead of the respective function.
#' Possible changes of \code{ForestTools} in the future will NOT change the workflow of \code{CENITH}
#'
#' @author Andreas Schönberg
#' @import raster
#' @import ForestTools
#' @keywords package
NULL
#' @docType data
#' @name lau_AOI_chm - data
#' @title A canopy height model (chm) of the AOI for the tutorial
#' @description  chm of the AOI in the Lautaret valley in the French Alps. Resolution 0.15 m.
#' @format \code{"raster::raster"}
NULL
#' @docType data
#' @name lau_AOI_rgb - data
#' @title An RGB Raster image of the AOI for the tutorial
#' @description RGB Raster of the AOI in the Lautaret valley in the French Alps. Resolution 0.15 m.
#' @format \code{"raster::raster"}
NULL
#' @docType data
#' @name lau_chm - data
#' @title A canopy height model (chm)
#' @description CHM of trees in site 1 of 3 in the Lautaret valley in the French Alps. Resolution 0.15 m.
#' @format \code{"raster::raster"}
NULL
#' @docType data
#' @name lau_chm2 - data
#' @title A canopy height model (chm)
#' @description CHM of trees in site 2 of 3 in the Lautaret valley in the French Alps. Resolution 0.15 m.
#' @format \code{"raster::raster"}
NULL
#' @docType data
#' @name lau_chm3 - data
#' @title A canopy height model (chm)
#' @description CHM of trees in site 3 of 3 in the Lautaret valley in the French Alps. Resolution 0.15 meter.
#' @format \code{"raster::raster"}
NULL
#' @docType data
#' @name lau_vp - data
#' @title Supervised computed Treepositions
#' @description Estimated position of trees for the first example chm in the Lautaret Valley. It is used to validate the segmentation.
#' @format \code{"rgdal::readOGR"}
NULL
#' @docType data
#' @name lau_vp_site2 - data
#' @title Supervised computed Treepositions
#' @description Estimated position of trees for the second example chm in the Lautaret Valley. It is used to validate the segmentation.
#' @format \code{"rgdal::readOGR"}
NULL
#' @docType data
#' @name lau_vp_site3 - data
#' @title Supervised computed Treepositions
#' @description Estimated position of trees for the third example chm in the Lautaret Valley. It is used to validate the segmentation.
#' @format \code{"rgdal::readOGR"}
NULL


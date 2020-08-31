#' CENITH Tree Segmentation Tool
#'
#' @name LEGION - package
#' @docType package
#' @title CENITH Tree Segmentation Tool
#' @description In general the \code{CENITH} package is a wrapper for the \code{ForestTools} package (Andrew Plowright and Jean-Romain Roussel,2020).
#' \code{CENITH} is designed to test for best fitting parameters to performe a Tree Segmentation. Furher it provides a CrossValidation for Tree Segmentation.
#'
#' @note For a more independent and stabel package the source code of some fucntions provided by \code{ForestTools} is used instead of the respectiv function.
#' Possible changes of \code{ForestTools} in the future will NOT change the workflow of \code{CENITH}
#'
#' @author Andreas Schönberg
#' @import raster
#' @import ForestTools
#' @keywords package
NULL
#' @docType data
#' @name lau_AOI_chm - data
#' @title A canopy height model (chm) for the tutorial AOI
#' @description AOI chm of Tutorial in the Lautaret vally in the france alps. Resolution 0.15 meter.
#' @format \code{"raster::raster"}
NULL
#' @docType data
#' @name lau_AOI_rgb - data
#' @title An rgb image for the tutorial AOI
#' @description AOI rgb of Tutorial in the Lautaret vally in the france alps. Resolution 0.15 meter.
#' @format \code{"raster::raster"}
NULL
#' @docType data
#' @name lau_chm - data
#' @title A canopy height model (chm)
#' @description CHM of some trees in side 1/3 in the Lautaret vally in the france alps. Resolution 0.15 meter.
#' @format \code{"raster::raster"}
NULL
#' @docType data
#' @name lau_chm2 - data
#' @title A canopy height model (chm)
#' @description CHM of some trees in side 2/3 in the Lautaret vally in the france alps. Resolution 0.15 meter.
#' @format \code{"raster::raster"}
NULL
#' @docType data
#' @name lau_chm3 - data
#' @title A canopy height model (chm)
#' @description CHM of some trees in side 3/3 in the Lautaret vally in the france alps. Resolution 0.15 meter.
#' @format \code{"raster::raster"}
NULL
#' @docType data
#' @name lau_vp - data
#' @title Supervised computed Treepositions
#' @description Estimated Position of trees for the example chm area in the Lautaret Vally. Used to validate Segmentation.
#' @format \code{"rgdal::readOGR"}
NULL
#' @docType data
#' @name lau_vp_side2 - data
#' @title Supervised computed Treepositions
#' @description Estimated Position of trees for the second example chm area in the Lautaret Vally. Used to validate Segmentation.
#' @format \code{"rgdal::readOGR"}
NULL
#' @docType data
#' @name lau_vp_side3 - data
#' @title Supervised computed Treepositions
#' @description Estimated Position of trees for the thrid example chm area in the Lautaret Vally. Used to validate Segmentation.
#' @format \code{"rgdal::readOGR"}
NULL



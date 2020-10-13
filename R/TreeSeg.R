#' Compute TreeCrown Segments
#' @description computes polygon segments for TreeCrowns based on watershed algorithm
#' @param chm raster -  Canopy Height Model RasterLayer.
#' @param a numeric - single value for MovingWindow.
#' @param b numeric - single value for MovingWindow.
#' @param h numeric - maximum height of trees (in meter) to detect trees.
#' @param MIN numeric - the minimum area for crowns. Smaller polygons are cropped.
#' @param MAX numeric - the maximum area for crowns. Larger polygons are cropped.
#' @param CHMfilter numeric - uses a sum filter on the chm with a MovingWindow of (x*x), which must be odd. Default=1 no filter.
#' @param silent bolean - if TRUE the function will not print any progress messages (default=FALSE).
#' @return returns a PolygonLayer with segments
#' @details uses a MovingWindow of x*a+b to detect local maxima in a chm to compute TreeCrown Segments
#' * parameter selection - use \code{\link{BestSegVal}} to automated detect best fitting parameters for a, b, h, MIN and filter.
#' @note A 'brute force' segmentation with random parameters is not recommended. TreeSeg is mainly used to compute segments AFTER the validation of best fitting parameters with \code{\link{BestSegVal}}.
#' Further to estimate the quality of the computed polygons it is recommended to use \code{\link{TreeSegCV}} for a x-fold CrossValdiation over x different subareas. For full workflow see the 'Tutorial'.
#' @seealso \code{\link{BestSegVal}}
#' @seealso \code{\link{TreeSegCV}}
#' @author Andreas Schönberg
#' @examples
#' # load data
#' chmpath <-system.file("extdata","lau_chm.tif",package = "CENITH")
#' chm <- raster::raster(chmpath)
#' # take a look on the data
#' plot(chm)
#' # NOTE: the exmple should NOT show to get optimal results (for this see 'BestSegVal')
#' # start segmentation
#' x <- TreeSeg(chm,a=0.3,b=0.7,h=13)
#' length(x)# amount of trees
#' # compare result with chm
#' mapview(chm)+x
#' # clip min and or max polygons
#' y <-TreeSeg(chm,a=0.3,b=0.7,h=13,MIN=10)
#' length(y)# amount of trees
#' mapview(chm)+y
#' @export TreeSeg
#' @aliases TreeSeg

TreeSeg <- function(chm=NULL,a,b,h,MIN=0,MAX=1000,CHMfilter=1,silent=FALSE){
  # function ForestTools vwf cleaned from cat code #############################
  vwf_clean <-function (CHM, winFun, minHeight = NULL, maxWinDiameter = NULL,
                              minWinNeib = "queen", verbose = FALSE)
  {
    if (verbose)
      if (!minWinNeib %in% c("queen", "rook"))
        stop("Invalid input for 'minWinNeib'. Set to 'queen' or 'rook'")
    CHM.crs <- as.character(raster::crs(CHM))
    CHM.prj <- regmatches(CHM.crs, regexpr("(?<=proj=).*?(?=\\s)",
                                           CHM.crs, perl = TRUE))
    if (length(CHM.prj) > 0 && CHM.prj %in% c(c("latlong", "latlon",
                                                "longlat", "lonlat"))) {
      warning("'CHM' map units are in degrees. Ensure that 'winFun' outputs values in this unit.")
    }
    roundRes <- round(raster::res(CHM), 5)
    if (roundRes[1] != roundRes[2])
      stop("Input 'CHM' does not have square cells")
    if (roundRes[1] == 0)
      stop("The map units of the 'CHM' are too small")
    if (!is.null(minHeight) && minHeight <= 0)
      stop("Minimum canopy height must be set to a positive value.")
    CHM.rng <- suppressWarnings(raster::cellStats(CHM, range))
    names(CHM.rng) <- c("min", "max")
    if (is.infinite(CHM.rng["max"]) | is.infinite(CHM.rng["min"])) {
      stop("Input 'CHM' does not contain any usable values. Check input data.")
    }
    if (!is.null(minHeight)) {
      if (minHeight >= CHM.rng["max"])
        stop("'minHeight' is set to a value higher than the highest cell value in 'CHM'")
      if (minHeight > CHM.rng["min"]) {
        CHM[CHM < minHeight] <- NA
        CHM.rng["min"] <- minHeight
      }
    }
    if (verbose)
      radii <- seq(floor(winFun(CHM.rng["min"])), ceiling(winFun(CHM.rng["max"])),
                   by = roundRes[1])
    radii <- radii[radii >= roundRes[1]]
    if (length(radii) == 0) {
      warning("The maximum window radius computed with 'winFun' is smaller than the CHM's resolution",
              "\nA 3x3 cell search window will be uniformly applied",
              "\nUse a higher resolution 'CHM' or adjust 'winFun' to produce wider dynamic windows")
      radii <- roundRes[1]
    }
    maxDimension <- (max(radii)/roundRes[1]) * 2 + 1
    if (!is.null(maxWinDiameter) && maxDimension > maxWinDiameter) {
      stop("Input function for 'winFun' yields a window diameter of ",
           maxDimension, " cells, which is wider than the maximum allowable value in 'maxWinDiameter'.",
           "\nChange the 'winFun' function or set 'maxWinDiameter' to a higher value (or to NULL).")
    }
    windows <- lapply(radii, function(radius) {
      win.mat <- raster::focalWeight(raster::raster(resolution = roundRes),
                                     radius, type = "circle")
      if (nrow(win.mat) == 3 && minWinNeib == "queen")
        win.mat[] <- 1
      win.pad <- raster::extend(raster::raster(win.mat), (maxDimension -
                                                            ncol(win.mat))/2, value = 0)
      win.vec <- as.vector(win.pad != 0)
      return(win.vec)
    })
    names(windows) <- radii
    .vwMax <- function(x, ...) {
      centralValue <- x[length(x)/2 + 0.5]
      if (is.na(centralValue)) {
        return(NA)
      }
      else {
        radius <- winFun(centralValue)
        window <- windows[[which.min(abs(as.numeric(names(windows)) -
                                           radius))]]
        return(if (max(x[window], na.rm = TRUE) == centralValue) 1 else 0)
      }
    }
    lm.pts <- raster::rasterToPoints(raster::focal(CHM, w = matrix(1,
                                                                   maxDimension, maxDimension), fun = .vwMax, pad = TRUE,
                                                   padValue = NA), fun = function(x) x == 1, spatial = TRUE)
    lm.pts[["height"]] <- raster::extract(CHM, lm.pts)
    lm.pts[["winRadius"]] <- winFun(lm.pts[["height"]])
    lm.pts[["treeID"]] <- 1:length(lm.pts)
    return(lm.pts)
  } ############################################################################

  # function ForestTools mcws cleaned from cat code #############################
  mcws_clean <- function (treetops, CHM, minHeight = 0, format = "raster", OSGeoPath = NULL,
                                verbose = FALSE)
  {
    if (verbose)
      if (!toupper(format) %in% c("RASTER", "POLYGONS", "POLYGON",
                                  "POLY")) {
        stop("'format' must be set to either 'raster' or 'polygons'")
      }
    if (verbose)
      CHM.max <- suppressWarnings(raster::cellStats(CHM, "max"))
    if (is.infinite(CHM.max)) {
      stop("Input CHM does not contain any usable values.")
    }
    if (minHeight > CHM.max) {
      stop("'minHeight' is set higher than the highest cell value in 'CHM'")
    }
    raster::crs(CHM) <- raster::crs(treetops)
    treetopsHgts <- raster::extract(CHM, treetops)
    treetops <- treetops[!is.na(treetopsHgts) & treetopsHgts >=
                           minHeight, ]
    if (length(treetops) == 0) {
      stop("No usable treetops. Treetops are either outside of CHM's extent, or are located elow the 'minHeight' value")
    }
    if (!"treeID" %in% names(treetops)) {
      warning("No 'treeID' found for input treetops. New 'treeID' identifiers will be added to segments")
      treetops[["treeID"]] <- 1:length(treetops)
    }
    else {
      if (any(treetops[["treeID"]] == 0))
        stop("'treeID' cannot be equal to 0")
      if (any(duplicated(treetops[["treeID"]])))
        warning("Duplicate 'treeID' identifiers detected")
    }
    if (verbose)
      CHM.mask <- is.na(CHM) | CHM < minHeight
    CHM[CHM.mask] <- 0
    if (verbose)
      ttops.ras <- raster::rasterize(treetops, CHM, "treeID", background = 0)
    CHM.img <- imager::as.cimg(raster::as.matrix(CHM))
    ttops.img <- imager::as.cimg(raster::as.matrix(ttops.ras))
    if (verbose)
      ws.img <- imager::watershed(ttops.img, CHM.img)
    ws.ras <- raster::raster(vals = ws.img[, , 1, 1], nrows = nrow(CHM),
                             ncols = ncol(CHM), ext = raster::extent(CHM), crs = raster::crs(CHM))
    ws.ras[CHM.mask] <- NA
    if (toupper(format) %in% c("POLYGONS", "POLYGON", "POLY")) {
      if (verbose)
        if (is.null(OSGeoPath)) {
          polys <- raster::rasterToPolygons(ws.ras)
          polys <- rgeos::gUnaryUnion(polys, id = polys[["layer"]])
          polys <- sp::disaggregate(polys)
        }
      else {
        polys <- APfun::APpolygonize(ws.ras, OSGeoPath = OSGeoPath)
      }
      if (verbose)
        polys.over <- sp::over(polys, treetops)
      polys.out <- sp::SpatialPolygonsDataFrame(polys, subset(polys.over,
                                                              select = which(names(polys.over) != "treeID")))
      polys.out <- polys.out[match(treetops[["treeID"]], polys.over[,
                                                                    "treeID"]), ]
      if (verbose)
        if ("crownArea" %in% names(polys.out)) {
          i <- 1
          while (paste0("crownArea", i) %in% names(polys.out)) i <- i +
              1
          crownArea <- paste0("crownArea", i)
          warning("Input data already has a 'crownArea' field. Writing new crown area values to the 'crownArea",
                  i, "' field")
        }
      else {
        crownArea <- "crownArea"
      }
      polys.out[[crownArea]] <- rgeos::gArea(polys.out, byid = TRUE)
      if (verbose)
        return(polys.out)
    }
    else {
      if (verbose)
        cat("..Returning crown outlines as a raster\n\nFinished at:",
            format(Sys.time(), "%Y-%m-%d, %X"), "\n")
      return(ws.ras)
    }
  }#############################################################################

  # setup UAVRST chmseg
  chmseg_FT <- function(treepos = NULL,
                        chm = NULL,
                        minTreeAlt = 2,
                        format = "polygons",
                        winRadius = 1.5,
                        verbose = FALSE) {

    if (class(treepos) %in% c("RasterLayer", "RasterStack", "RasterBrick")) {
      # add projection of chm TODO
      pr<-raster::crs(raster::projection(raster::raster(chm)))
      raster::projection(treepos) <- pr
      treepos <- raster::rasterToPoints(treepos,spatial = TRUE)
      # reformat it to the needs of mcws
      treepos@data$layer <- 1
      treepos@data$winRadius <- 1.5
      treepos@data$treeID <- seq(1:nrow(treepos@data))
      names(treepos)<-c("height","treeID","winRadius","layer")
    }
    # Crown segmentation
    crownsFT <- mcws_clean(treetops = treepos,
                                  CHM = chm,
                                  format = format,
                                  minHeight = minTreeAlt,
                                  verbose = verbose)
    return(crownsFT)
  } ############################################################################

  # TreeSeg

  # check inputs
  if(CHMfilter %% 2 == 0){
    stop("selected 'filter' sizes contain even values (use odd values only)")
  }
  # check filter
  if(CHMfilter!=1){
    if(silent==FALSE){
    cat(paste0("### Cenith computes ",CHMfilter,"*",CHMfilter," sum filter for chm ###",sep = "\n"))
    }
    chm <-raster::focal(chm,w=matrix(1/(CHMfilter*CHMfilter),nrow=CHMfilter,ncol=CHMfilter),fun=sum,na.rm=TRUE)
  }

  if(silent==FALSE){
  cat(paste0("### Cenith starts segmentation ###",sep = "\n"))
  }
  # compute Treepositions
  tpos = vwf_clean(chm,
                   winFun = function(x){x * a + b},
                   minHeight = h,
                   verbose = TRUE)

  # compute segments
  seg <- chmseg_FT(chm = chm,
                           treepos = tpos,
                           format = "polygons",
                           minTreeAlt = h,
                           verbose = TRUE)
  if(silent==FALSE){
    cat(paste0("detected ",length(tpos), " trees ###",sep = "\n"))
  }

  # clip min and max
  if(silent==FALSE){
    cat(paste0("### Cenith starts cropping to MIN and MAX ###",sep = "\n"))
  }
  seg_min <- seg[seg$crownArea>MIN,]
  seg_max <- seg_min[seg_min$crownArea<MAX,]
  seg <- seg_max
  if(length(seg)==0){
    stop(" !!! after cliping to MIN and MAX no polygons are returned, setting results to NA",sep = "\n")}

  if(silent==FALSE){
    cat(paste0("clipped: ",length(tpos)-length(seg)," polygons. ", length(seg)," trees remaining",sep = "\n"))
    cat(paste0("### Cenith finished segmentation ###",sep = "\n"))

  }
  return(seg)
}#end of main function

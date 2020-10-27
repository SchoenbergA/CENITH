#' N-Fold Cross-validation for TreeCrown Segments
#' @description performs an n-fold cross validation to estimate the performance of a segmentation 'model' for an AOI
#' @param sites list - a list of chm RasterLayers (see details)
#' @param a numeric - single value for MovingWindow.
#' @param b numeric - single value for MovingWindow.
#' @param h numeric - maximum height of trees (in meter) to detect trees.
#' @param vps list - a list of PointLayers with estimated positions of trees (see details)
#' @param MIN numeric - the minimum area for crowns. Smaller poylgons are cropped.
#' @param MAX numeric - the maximum area for crowns. Larger polygons are cropped.
#' @param CHMfilter numeric - uses a sum filter on the chm with a MovingWindow of (x*x), which must be odd. Default=1 no filter.
#' @return returns a table with quality values for each fold (site) and calculated overall performance (mean values) for all sites.
#' @details
#' * 'sites' and 'vps' must be lists with the same order of chm and the respective validation points.
#' * parameters of the model: a, b, h, MIN, MAX should be calculated by \code{\link{BestSegVal}}.
#' * The overall performance is the mean of all sites.
#' * For the details of the table values see \code{\link{BestSegVal}}.
#' @note The overall performance helps to estimate the precision for an AOI but does NOT give the "exact" precision. More folds will increase the
#' expressiveness but will need more time to set validation points. Further the supervised setting of validation points is highly subjectiv and does not have to correlate with the real amount and or position of trees.
#' @author Andreas Schönberg
#' @seealso \code{\link{BestSegVal}}
#' @examples
#' require(CENITH)
#' require(mapview)
#' require(raster)
#' require(rgdal)
#' # load data
#' chm  <- raster::raster(system.file("extdata","lau_chm.tif",package = "CENITH"))
#' chm2 <- raster::raster(system.file("extdata","lau_chm_side2.tif",package = "CENITH"))
#' chm3 <- raster::raster(system.file("extdata","lau_chm_side3.tif",package = "CENITH"))
#' vp <- rgdal::readOGR(system.file("extdata","lau_vp.shp",package = "CENITH"))
#' vp2 <- rgdal::readOGR(system.file("extdata","lau_vp_side2.shp",package = "CENITH"))
#' vp3 <- rgdal::readOGR(system.file("extdata","lau_vp_side3.shp",package = "CENITH"))
#' # handle CRS string
#' crs(vp) <-crs(chm)
#' crs(vp2)<-crs(chm)
#' crs(vp3)<-crs(chm)
#'
#' # list all sites and validation points
#' chmlist <- list(chm,chm2,chm3)
#' vplist <- list(vp,vp2,vp3)
#'
#' # run 3 fold cross validation with parameters computed by 'BestSegVal' (from example)
#' cv <- CENITH::TreeSegCV(sites=chmlist,a=0.3,b=0.5,h=0.5,MIN=5,MAX=1000,CHMfilter=3,vps=vplist)
#' cv
#' ### the model trained with BestSegVal on site 1 reaches a overall performance of 0.77 @ 0.12 for all tree sites.
#' ### Note that the performance on sites 2 and 3 is even better than on site 1 where it was trained. This effect is caused probably because sites 2 and 3 are more homogenious.
#' ### For sure this is just an example and it can happen that the performance on one site is worst. The overall mean performance is used to estimate the quality of an segmentation for an AOI by testing some subareas.


#' @export TreeSegCV
#' @aliases TreeSegCV

### tests
# test cat codes with lists

TreeSegCV <- function(sites,a,b,h,MIN,MAX,CHMfilter=1,vps){
  # cheking inputs

  # create dataframe to save informations
  cat(paste0("### CENITH starts ",length(sites),"-fold cross-validation ###",sep = "\n"))
  result <- data.frame(matrix(nrow = length(b), ncol = 8)) # ncol = n information stored
  # iteration
  # loop to iterate on varibale
  for (i in seq(1:length(sites))){

    cat(paste0("starting fold ",as.factor(i),"/",as.factor(length(sites))," ",sep = "\n"))
    # TreeSeg
    seg <- try(TreeSeg(chm=sites[[i]],a,b,h,MIN,MAX,CHMfilter,silent=TRUE),silent = TRUE)
    # handle error catch
    if(class(seg)=="try-error"){
      cat(paste0(" !!! iteration a=",a," b=",b," h=",h," leads to an error on site:", names(sites[[i]]), "setting to 'NA' in results ",sep="\n"))
      result[i, 1] <- a
      result[i, 2] <- b
      result[i, 3] <- h
      #absolut results
      result[i, 4] <- NA
      result[i, 5] <- NA
      result[i, 6] <- NA
      result[i, 7] <- NA
      result[i, 8] <- NA
      # rates
      result[i, 9]  <- NA
      result[i, 10] <- NA
      result[i, 11] <- NA
      result[i, 12] <- NA


    } # end if error
    else{


      ### save results in dataframe################################

      # get stats of overlapping
      stat <- ForestTools::sp_summarise(vps[[i]], seg) # compute data points in polygons
      stat[is.na(stat$TreeCount)] <- 0 # na to 0

      # get n trees in poygons
      TCO <- sum(stat$TreeCount<1) # amount polygon without any tree (miss) <- Oversegmented
      TC1 <- sum(stat$TreeCount==1) # amount polygon with exact 1 tree (hit)
      TCU <- sum(stat$TreeCount>1) # amount polygons with more than 1 tree (miss) <- Undersegmented
      TC2 <- sum(stat$TreeCount==2) # amount polygon with 2 Trees (miss)
      TC3 <- sum(stat$TreeCount==3) # amount polygon with 3 Trees (miss)
      TCX <- sum(stat$TreeCount>3)  # amount polygon more tha 3 Trees (miss)

      # calculate validation scores
      # absolute
      hit = paste(TC1,"/",length(vps[[i]]))

      # rates
      hitrate = TC1/length(vps[[i]])        # hits / seg amount of hits in relation to nseg
      over = TCO/length(stat$TreeCount)           # oversegmented, Segnents with no VP
      under = TCU/length(stat$TreeCount) # undersegmented, Segments with more than one tree

      # quality value calculation
      miss <- (over+(2*under))/2
      segQy <- paste0(round(hitrate,2)," @ ",round(miss,2))
      segQy

      # additional informations
      tseg = length(seg)
      area =  sum(seg$crownArea)#

      # write out informations in dataframe
      result[i, 1] <- names(sites[[i]])
      result[i, 2] <- a
      result[i, 3] <- b
      result[i, 4] <- h
      #absolut results
      result[i, 5] <- tseg
      result[i, 6] <- TC1
      result[i, 7] <- length(vps[[i]])
      result[i, 8] <- TCU
      result[i, 9] <- TCO
      result[i, 10] <- area
      # rates
      result[i, 11] <- hitrate
      result[i, 12] <- under
      result[i, 13] <- over
      result[i, 14] <- segQy
    } # end of more than null polygons
            # handle df
            #if(i==1){
            #  res <-result
#
 #           }    else {
  #            res2 <-result
   #           res= rbind(res,res2)
    #        }
    res=result

  }# end of iteration

  # calc means

  # write out information in dataframe
  names(res)<- c("sites","a","b","height","total_seg","hit","vp","under","over","area","hitrate","underrate","overrate","Seg_qualy")
  res[nrow(res)+1, 1] <- "Mean"
  res[nrow(res), 2] <- a
  res[nrow(res), 3] <- b
  res[nrow(res), 4] <- h
  # absolute results
  res[nrow(res), 5] <- mean(res[1:nrow(res)-1,5])
  res[nrow(res), 6] <- mean(res[1:nrow(res)-1,6])
  res[nrow(res), 7] <- mean(res[1:nrow(res)-1,7])
  res[nrow(res), 8] <- mean(res[1:nrow(res)-1,8])
  res[nrow(res), 9] <- mean(res[1:nrow(res)-1,9])
  # rates
  res[nrow(res), 10]  <- mean(res[1:nrow(res)-1,10])
  res[nrow(res), 11] <-  mean(res[1:nrow(res)-1,11])
  res[nrow(res), 12] <-  mean(res[1:nrow(res)-1,12])
  res[nrow(res), 13] <-  mean(res[1:nrow(res)-1,13])
  res[nrow(res), 14] <- paste0(round(mean(res[1:nrow(res)-1,11]),2),
                               " @ ",
                               round(
                                (mean(res[1:nrow(res)-1,13])+ 2* mean(res[1:nrow(res)-1,12]))/2
                                     ,2))

  cat(paste0("### CENITH finished ",length(sites),"fold cross validation ###",sep = "\n"))
  cat(paste0("Overall performance of model: ",res[nrow(res),14] ,sep = "\n"))
  return(res)
}# end of function



#' BETA Detect best input variables for TreeSegmentation by validation points
#' @description uses supervised computed Treepositions to validate best fitting values for a,b and height.
#' @param chm raster - RasterLayer with Canopy height model
#' @param a numeric - function for MovingWindow
#' @param b numeric - function for MovingWindow
#' @param h numeric - maximum height of Trees
#' @param vp Polygon - PointLayer with estimated Positions of Trees.
#' @param MIN numeric - minimum area for Crowns. smaller poylgons are cropped
#' @param MAX numeric - maximum area for Crowns. larger polygons are cropped
#' @param skipCheck development - bolean - if TRUE skips chekcing the inputs
#' @return returns a dataframe with several validation scores
#' @details
#' * if 'skipCheck' = TRUE the input checking is skipped and the itrations starts directly.
#' * has implemneted error catching. if an iteration would cause critical stop error, the loops continue. Returns NA for corrupted iterations
#' * Input for a,b,h supports - numeric, single, combination by c() or sequence by seq() to iterate over.
#' * Validation Point - supervised computed pointlayer. Use a GIS like Qgis to set points where Trees are estimated.
#' @author Andreas Schönberg
#' @examples
#' # load data
#' chmpath <-system.file("extdata","exp_chm.tif",package = "CENITH")
#' chm <- raster::raster(chmpath)
#' plot(chm)
#' vppath <-system.file("extdata","exp_vp.shp",package = "CENITH")
#' vp <- rgdal::readOGR(vppath)
#' # run validation
#' x <-BestSegValBETA(chm,seq(0.7,0.8,0.1),seq(0.5,0.7,0.1),c(13,8),vp)
#' # view best accuracy (diffenrent iterations can have the same hit rate, but differ in other values)
#' maxrow <- x[which.max(x$hit),] # search max vale but rturn only 1 value
#' maxhit <- maxrow$hit
#' x[which(x$hit==maxhit),]
#' # test error causing values
#' test1 <-BestSegValBETA(chm=chm,a=c(0.9,0.91,0.92,0.1),b=0.5,h=1,vp=vp,MIN = 40,MAX = 200,skipCheck = FALSE)
#' # if the function stops, retry.

#' @export BestSegValBETA
#' @aliases BestSegValBETA


BestSegValBETA<- function(chm,a,b,h,vp,MIN=0,MAX=1000,skipCheck=FALSE){
  if(skipCheck==FALSE){
  cat(paste0("### Cenith checking input ###",sep = "\n"))
  #checking projection
  if(compareCRS(chm,vp)==FALSE){
    stop("unequal CRS for 'chm' and 'vp' detected !")
  }
  # check input height is lesser than max height in chm
  if(cellStats(chm,max)<=max(h)){
    stop("input height 'h' values are higher than the highest cell value of teh CHM")
  }
  # estimate time remaining, trying up to tree iterations, if all fail, skipping time estimation
    cat("",sep = "\n")
    cat(paste0("### Cenith calculates estimate time to finish ###",sep = "\n"))
    cat("",sep = "\n")
  start <- Sys.time()
  ETA1 <-try(ForestTools::vwf(chm,
                                 winFun = function(x){x * sample(a,1) + sample(b,1)},
                                 minHeight = sample(h,1),
                                 verbose = TRUE),silent = TRUE)
      stop  <- Sys.time()
              if(class(ETA1)=="try-error"){
                start <- Sys.time()
                ETA2 <-try(ForestTools::vwf(chm,
                                               winFun = function(x){x * sample(a,1) + sample(b,1)},
                                               minHeight = sample(h,1),
                                               verbose = TRUE),silent = TRUE)
                stop  <- Sys.time()
                if(class(ETA2)=="try-error"){
                  start <- Sys.time()
                  ETA3 <-try(ForestTools::vwf(chm,
                                              winFun = function(x){x * sample(a,1) + sample(b,1)},
                                              minHeight = sample(h,1),
                                              verbose = TRUE),silent = TRUE)
                  stop  <- Sys.time()
                  if(class(ETA2)=="try-error"){
                    stop("unable to calculate ETA, tried 3 random iterations, skipping.")}
                }
              }




  # amount of iterations
  na <- length(a)
  nb <- length(b)
  nh <- length(h)
  niter <- nh*nb*na
  # test time for 1 iteration
    timedif <-difftime(stop,start,units = "min")
  # predict estimated time by multiply timedif with n iterations

  eta <- timedif*niter
  eta <-1.2*eta # add 20% because the iteration wrapper needs more time.
  cat(paste("requires",niter,"iterations @",round(timedif,4),units(timedif),"per iteration"),sep="\n")
  cat("",sep = "\n")
  cat(paste("estimated",round(eta,4),units(timedif),"to finish"))
  cat("",sep = "\n")
}# end skip check

  # setup function to iterate over b ############################################################
  loop_b <- function(chm,a,b,h,vp){


    # create dataframe to save informations
    result <- data.frame(matrix(nrow = length(b), ncol = 8)) # ncol = n information stored

    # loop to iterate on varibale
    for (j in seq(1:length(b))){

      cat(paste0("iterating over b ",as.factor(j),"/",as.factor(length(b))," ",sep = "\n"))
      # TreeSeg
      seg <- try(TreeSeg(chm=chm,a,b[j],h,MIN,MAX),silent = TRUE)
                      # handle error catch
                      if(class(seg)=="try-error"){
                        cat(paste0(" !!! iteration a=",a," b=",b[j]," h=",h," leads to an error, setting to 'NA' in results ",sep="\n"))
                        result[j, 1] <- a
                        result[j, 2] <- b[j]
                        result[j, 3] <- h
                        #absolut results
                        result[j, 4] <- NA
                        result[j, 5] <- NA
                        result[j, 6] <- NA
                        result[j, 7] <- NA
                        result[j, 8] <- NA
                        # rates
                        result[j, 9]  <- NA
                        result[j, 10] <- NA
                        result[j, 11] <- NA
                        result[j, 12] <- NA


      } # if error
      else{


      ### save results in dataframe################################

      # get stats of overlapping
      stat <- ForestTools::sp_summarise(vp, seg) # compute data points in polygons
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
      hit = paste(TC1,"/",length(vp))

      # rates
      hitrate = TC1/length(vp)        # hits / seg amount of hits in relation to nseg
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
      result[j, 1] <- a
      result[j, 2] <- b[j]
      result[j, 3] <- h
      #absolut results
      result[j, 4] <- tseg
      result[j, 5] <- hit
      result[j, 6] <- TCU
      result[j, 7] <- TCO
      result[j, 8] <- area
      # rates
      result[j, 9] <- hitrate
      result[j, 10] <- under
      result[j, 11] <- over
      result[j, 12] <- segQy
                } # end of more than null polygons


    }# end loop for b
    # output
    return(result)
  }# end of iteration over b #####################################################

  # setup function to iterate over a ############################################################
  loop_a <- function(chm,a,b,h,vp){

    for (i in seq(1:length(a))){
      cat("",sep = "\n")
      cat(paste0("### Cenith starts iterating over a ",as.factor(i),"/",as.factor(length(a))," ",sep = "\n"))


      if(i==1){
        res <-loop_b(chm,a[i],b,h,vp)

      }    else {
        res2 <-loop_b(chm,a[i],b,h,vp)
        res= rbind(res,res2)
      }




    }# end loop for a

    return(res)
  }# end of iteration a



  # BestSegVal

  for (c in seq(1:length(h))){
    cat("",sep = "\n")
    cat(paste0("### Cenith starts iterating over h ",as.factor(c),"/",as.factor(length(h))," ###",sep = "\n"))
    if(c==1){
      res <-loop_a(chm,a,b,h[c],vp)
    }    else {
      res2 <-loop_a(chm,a,b,h[c],vp)
      res= rbind(res,res2)}

  }# end loop h
  # output
  cat("",sep = "\n")
  cat(paste0("### Cenith finsihed Segmentation ###"),sep = "\n")

  names(res)<- c("a","b","height","total_seg","hit/vp","under","over","area","hitrate","underrate","overrate","Seg_qualy")
  #names(res)<- c("a","b","height","hit","tp/vp_rate","tpos/vp","miss","area","empty")
  return(res)
}# end core fucntion

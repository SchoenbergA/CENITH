#' Test for best performing parameters for TreeSeg TreeSegmentation using Validation Treepositions
#' @description Iterates over a,b,h,MIN and supports iterating over filtered chms. Uses supervised computed Treepositions to validate best fitting values for a,b, height, MIN and CHM filters to detect trees.
#' @param chm raster -  Canopy Height Model RasterLayer derived form LiDAR data
#' @param a numeric - single value, combination of values or sequence for MovingWindow
#' @param b numeric - single value, combination of values or sequence for MovingWindow
#' @param h numeric - single value, combination of values or a sequence for the maximum height of trees (in meter) to detect trees.
#' @param vp shp - PointLayer with estimated Positions of trees (see details).
#' @param MIN numeric - single value, combination of values or a sequence of minimum area for crowns. Smaller polygons are cropped Default= 0
#' @param MAX numeric - the maximum area for crowns. Larger polygons are cropped. Default=1000
#' @param filter numeric - single value, combination of values or a sequence for filtersize, uses a sum filter on the chm with a MovingWindow of (x*x), which must be odd., default= 1 (no filter.)
#' @param skipCheck for development - bolean - if TRUE skips checking input data and estimation of ETA.

#' @return returns a dataframe with several validation scores (see details).
#' @details
#' * 'start and go to sleep' The function has implemented error catching. If an iteration would cause a critical stop error, the loops continue. Returns NA for corrupted iterations.
#' * ETA: before starting the iterations, up to 3 random combinations are tested and the processing time is recorded. The processing time is multiplied with the total amount of iterations to calculate the ETA.
#' * Input for a, b, h, MIN and filter supports - numeric, single and combination: c(), or sequence: seq() to iterate over: (e.g a=0.5 ; a=c(0.3,0.5) ; a=seq(0.1,0.9,0.05).
#' * filter - it could be helpful to filter the raw chm to 'smooth' small peaks and 'holes' in crowns. Helps to suppress the segmentation of many tiny objects (like stones) and leads to polygons with less 'holes'.
#' * vp - validation point: supervised computed pointlayer. Use a GIS like Qgis to set points where Trees are estimated (see Tutorial).
#' * if 'skipCheck' = TRUE the input checking is skipped and the iterations start directly.
#' * result table
#'    + a, b, height, MIN, chm - the used parameters
#'    + total_seg - total computed segments
#'    + hit/vp - amount of segments computed for trees marked with validation points
#'    + under - undersegmnetation in absolut values, amount of segments which contain more than one validation point.
#'    + over - oversegmentation in absolut values, amount of segments without validation point
#'    + area - the total area of segments
#'    + hitrate - percentage of segments which contain exactly one validation point in relation to the total validation points.
#'    + underrate - undersegmentation in relation to the total segments.
#'    + overrate - oversegmentation in relation to the total segments.
#'    + Seg_qualy - quality of segmentation:  hitrate @ combined over- and undersegmentations. A segment with more than one validation points is estimated to be better than segments without validation points. Therefore the 'miss' value is (over+(2*under))/2.
#' @author Andreas Schönberg
#' @note The 'brute force' approach to iterate over many parameters may result in very long time to finish. Preselected smaller samples may be more efficient (see Tutorial).
#' @examples
#' ### NOTE: this example is used to show the usage of 'BestSegVal'. It is NOT used to show a workflow for best Results (see 'Tutorial' vignette for workflow strategies)
#' ### further NOTE: to reduce time usage for the example only small iteartions are used to get an overlook for the functionallities.
#' # load data
#' chmpath <-system.file("extdata","lau_chm.tif",package = "CENITH")
#' chm <- raster::raster(chmpath)
#' vppath <-system.file("extdata","lau_vp.shp",package = "CENITH")
#' vp <- rgdal::readOGR(vppath)
#' # take a look at the data (chm with estimated Treeposition)
#' mapview(chm)+vp
#' # start BestSegVal with parameters for Moving window and minimum height
#' x <- BestSegVal(chm,a=c(0.3,0.5),b=c(0.5,0.7),h=c(0.1,1),vp = vp,filter=c(1,3))
#' # check for best hitrate
#' maxrow <- x[which.max(x$hitrate),] # search max vale but rturn only 1 value
#' maxhit <- maxrow$hitrate
#' x[which(x$hitrate==maxhit),]
#' # visualise best combination (due to Segment Quality)
#' y <- TreeSeg(chm,0.3,0.7,0.1,CHMfilter = 3)
#' # show result segments with validation points and chm
#' mapview(chm)+vp+y
#' # BestSegVal based on this results keep a,b and filter 3 with MIN values to reduce oversegmentation and different heights.
#' z <- BestSegVal(chm,a=0.3,b=0.7,h=c(0.1,0.2,0.5,1),MIN=c(10,100),vp = vp,filter=3)
#' ### Note that row 1 and 2 lead to much lesser oversegmnetation just differ in area.
#' # visualize
#' v <- TreeSeg(chm,0.3,0.7,0.1,MIN=10,CHMfilter = 3)
#' mapview(chm)+vp+v

#' @export BestSegVal
#' @aliases BestSegVal


BestSegVal<- function(chm,a,b,h,vp,MIN=0,MAX=1000,skipCheck=FALSE,filter=1){
  if(skipCheck==FALSE){
    cat(paste0("### Cenith checking input ###",sep = "\n"))
    #check for wrong sizes input
    if(any(filter %% 2 == 0)){
      stop("selected 'filter' sizes contain even values (use odd values only)")
    }


    #checking projection
    if(compareCRS(chm,vp)==FALSE){
      stop("unequal CRS for 'chm' and 'vp' detected !")
    }
    # check input height is lesser than max height in chm
    if(cellStats(chm,max)<=max(h)){
      stop("input height 'h' values are higher than the highest cell value of teh CHM")
    }

    # estimate time remaining, trying up to tree iterations, if all fail, skipping time estimation

    cat(paste0("### Cenith calculates estimate time to finish ###",sep = "\n"))

    start <- Sys.time()

    ETA1 <-try(TreeSeg(chm,a=sample(a,1),b=sample(b,1),h=sample(h,1),MIN=sample(MIN,1),MAX,CHMfilter=1,silent=TRUE),silent = T)
    stop  <- Sys.time()
    if(class(ETA1)=="try-error"){
      start <- Sys.time()
      ETA2 <-try(TreeSeg(chm,a=sample(a,1),b=sample(b,1),h=sample(h,1),MIN=sample(MIN,1),MAX,CHMfilter=1,silent=TRUE),silent = T)
      stop  <- Sys.time()
      if(class(ETA2)=="try-error"){
        start <- Sys.time()
        ETA3 <-try(TreeSeg(chm,a=sample(a,1),b=sample(b,1),h=sample(h,1),MIN=sample(MIN,1),MAX,CHMfilter=1,silent=TRUE),silent = T)
        stop  <- Sys.time()
        if(class(ETA2)=="try-error"){
          stop("unable to calculate ETA, tried 3 random iterations, skipping. Please restart.")}
      }
    }




    # amount of iterations
    na <- length(a)
    nb <- length(b)
    nh <- length(h)
    nm <- length(MIN)
    nc <- length(filter)
    niter <- nh*nb*na*nm*nc

    # test time for 1 iteration
    timedif <-difftime(stop,start,units = "min")
    # predict estimated time by multiply timedif with n iterations

    eta <- timedif*niter
    eta <-1.2*eta # add 20% because the iteration wrapper needs more time.
    cat(paste("requires",niter,"iterations @",round(timedif,4),units(timedif),"per iteration"),sep="\n")
    cat(paste("estimated",round(eta,4),units(timedif),"to finish"),sep = "\n")
    cat("#----------------------------------------#",sep = "\n")


    # warning / wait for input if niter is too big
    #if(niter>50&&niter<100){
    #  cat("",sep = "\n")
    #  cat(paste0("WARNING over 50 iterations detected this may take long time to finsih, maybe abort and start smaller trys (see details)",sep = "\n"))
    #}
    #if(niter>=100){
    #  abort <- readline(prompt="More than 100 iterations detected! Are you sure you want to continue? 1:YES, 2:NO, plz abort ")
    #  if(abort==2){
    #    stop("aborded on users desicion, too much iterations detected")
    #  }
    #}
  }# end skip check

  # setup function to iterate over b ############################################################
  loop_b <- function(chm,a,b,h,vp,MIN,MAX){


    # create dataframe to save informations
    result <- data.frame(matrix(nrow = length(b), ncol = 8)) # ncol = n information stored

    # loop to iterate on varibale
    for (j in seq(1:length(b))){

      cat(paste0("iterating over b ",as.factor(j),"/",as.factor(length(b))," ",sep = "\n"))
      # TreeSeg
      seg <- try(TreeSeg(chm=chm,a,b[j],h,MIN,MAX,silent = TRUE,CHMfilter=1),silent=T)
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
        result[j, 13] <- MIN
        result[j, 14] <- names(chm)


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
        result[j, 13] <- MIN
        result[j, 14] <- names(chm)
      } # end of more than null polygons


    }# end loop for b
    # output
    return(result)
  }# end of iteration over b #####################################################

  # setup function to iterate over a ############################################################
  loop_a <- function(chm,a,b,h,vp,MIN,MAX){
    for (i in seq(1:length(a))){
      cat("",sep = "\n")
      cat(paste0("### Cenith starts iterating over a ",as.factor(i),"/",as.factor(length(a))," ",sep = "\n"))
      if(i==1){
        res <-loop_b(chm,a[i],b,h,vp,MIN,MAX)

      }    else {
        res2 <-loop_b(chm,a[i],b,h,vp,MIN,MAX)
        res= rbind(res,res2)
      }


    }# end loop for a

    return(res)
  }# end of iteration a



  # iteration h
  loop_h <- function(chm,a,b,h,vp,MIN,MAX){
    for (c in seq(1:length(h))){

      cat(paste0("### Cenith starts iterating over h ",as.factor(c),"/",as.factor(length(h))," ###",sep = "\n"))
      if(c==1){
        res <-loop_a(chm,a,b,h[c],vp,MIN,MAX)
      }    else {
        res2 <-loop_a(chm,a,b,h[c],vp,MIN,MAX)
        res= rbind(res,res2)}

    }# end loop h
    return(res)
  }# end iteration h

  # iteration h
  loop_MIN <- function(chm,a,b,h,vp,MIN,MAX){
    for (k in seq(1:length(MIN))){

      cat(paste0("### Cenith starts iterating over MIN ",as.factor(k),"/",as.factor(length(MIN))," ###",sep = "\n"))
      if(k==1){
        res <-loop_h(chm,a,b,h,vp,MIN[k],MAX)
      }    else {
        res2 <-loop_h(chm,a,b,h,vp,MIN[k],MAX)
        res= rbind(res,res2)}

    }# end loop h
    return(res)
  }# end iteration h





  # BestSeg Core


  for(l in seq(1:length(filter))){
    cat("",sep = "\n")
    cat(paste0("### Cenith starts iterating over CHM ",as.factor(l)),"/",(length(filter))," ###")
    cat("",sep = "\n")
    if(l==1){
      #for no filtering with 1
      if(filter[l]==1){
        res <-loop_MIN(chm,a,b,h,vp,MIN,MAX)
      } else {
        chmf <- raster::focal(chm,w=matrix(1/(filter[l]*filter[l]),nrow=filter[l],ncol=filter[l]),fun=sum,na.rm=TRUE)
        names(chmf) <- paste0(names(chm),"_",as.factor(filter[l]))
        res <-loop_MIN(chmf,a,b,h,vp,MIN,MAX)
      }
    }    else {
      chmf <- raster::focal(chm,w=matrix(1/(filter[l]*filter[l]),nrow=filter[l],ncol=filter[l]),fun=sum,na.rm=TRUE)
      names(chmf) <- paste0(names(chm),"_",as.factor(filter[l]))
      res2 <-loop_MIN(chmf,a,b,h,vp,MIN,MAX)
      res= rbind(res,res2)}

  }# end core function


  # output
  cat("",sep = "\n")
  cat(paste0("### Cenith finsihed Best Parameters Validation ###"),sep = "\n")

  names(res)<- c("a","b","height","total_seg","hit/vp","under","over","area","hitrate","underrate","overrate","Seg_qualy","MIN","chm")
  #names(res)<- c("a","b","height","hit","tp/vp_rate","tpos/vp","miss","area","empty")
  return(res)
}# end core fucntion

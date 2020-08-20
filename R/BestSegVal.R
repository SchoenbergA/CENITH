#' NOPAR EXPRMTL Detect best input variables for TreSeg TreeSegmentation using Validation Treepositions
#' @description Iterates over a,b,h,MIN and supports iterating over filtered chms.uses supervised computed Treepositions to validate best fitting values for a,b, height, MIN and CHM filters to detect trees.
#' @param chm raster - RasterLayer with Canopy height model
#' @param a numeric - function for MovingWindow
#' @param b numeric - function for MovingWindow
#' @param h numeric - minimum height of Trees to detect
#' @param vp Polygon - PointLayer with estimated Positions of Trees (see details).
#' @param MIN numeric - minimum area for Crowns. smaller poylgons are cropped
#' @param MAX numeric - maximum area for Crowns. larger polygons are cropped
#' @param filter numeric - uses a sum filter on the chm with Moving Window of (x*x), must be odd (see details)
#' @param skipCheck development - bolean - if TRUE skips checking input data

#' @return returns a dataframe with several validation scores (see details).
#' @details
#' * 'start and go sleep' The function has implemented error catching. If an iteration would cause critical stop error, the loops continue. Returns NA for corrupted iterations.
#' * Input for a,b,h,MIN,filter supports - numeric, single, combination by c() or sequence by seq() to iterate over (e.G a=0.5 ; a=c(0.3,0.5) ; a=seq(0.1,0.9,0.05).
#' * filter - i could be helpful to filter the raw chm to 'smooth' small peaks and suppress the segmentation of many tiny objects.
#' * Validation Point - supervised computed pointlayer. Use a GIS like Qgis to set points where Trees are estimated (see Tutorial).
#' * if 'skipCheck' = TRUE the input checking is skipped and the iterations start directly.
#' * ! even if the function support iteration over all parameters (except MAX due to estimated no need for MAX Crownareas to clip) it should not be used to iterate over sequences which would be useless.
#' @author Andreas Schönberg
#' @note The 'Brute Force' approach to iterate over many parameters may result in very long time to finish. Preselected smaller samples may be more efficient (see Tutorial)
#' @examples
#' # load data
#' chmpath <-system.file("extdata","lau_chm.tif",package = "CENITH")
#' chm <- raster::raster(chmpath)
#' vppath <-system.file("extdata","lau_vp.shp",package = "CENITH")
#' vp <- rgdal::readOGR(vppath)
#' # take a look on the data
#' mapview(chm)+vp
#' ### start iteration runs
#' # simple (core values)
#' x <- BestSegVal(chm,a=seq(0.1,0.7,0.1),b=seq(0.1,0.7,0.2),h=1,vp=vp,filter=1)
#' # check hole table
#' x
#' # check best results for 'hitrate'
#' maxrow <- x[which.max(x$hitrate),] # search max vale but rturn only 1 value
#' maxhit <- maxrow$hitrate
#' x[which(x$hit==maxhit),]
#' # more complex run
#' y <- BestSegVal(a=c(0.1,0.7),b=c(0.1,0.7),h=1,MIN=20,filter=c(1,3),vp=vp)
#' # check hole table
#' y
#' # check best results for 'hitrate'
#' maxrow <- y[which.max(y$hitrate),] # search max value but return only 1 value
#' maxhit <- maxrow$hitrate
#' y[which(y$hitrate==maxhit),]
#' @export BestSegVal
#' @aliases BestSegVal


BestSegVal<- function(chm,a,b,h,vp,MIN=0,MAX=1000,skipCheck=FALSE,filter=NULL){
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

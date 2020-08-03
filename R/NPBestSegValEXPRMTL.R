#' NOPAR EXPRMTL Detect best input variables for TreeSegmentation by validation points
#' @description Experimantal Version of BestSegVal. Iterates over a,b,h,MIN and supports iterating over filtered chms.uses supervised computed Treepositions to validate best fitting values for a,b and height.
#' @param chm raster - RasterLayer with Canopy height model
#' @param a numeric - function for MovingWindow
#' @param b numeric - function for MovingWindow
#' @param h numeric - maximum height of Trees
#' @param vp Polygon - PointLayer with estimated Positions of Trees.
#' @param MIN numeric - minimum area for Crowns. smaller poylgons are cropped
#' @param MAX numeric - maximum area for Crowns. larger polygons are cropped
#' @param skipCheck development - bolean - if TRUE skips chekcing the inputs
#' @param filter numeric values for the moving window, must be odd

#' @return returns a dataframe with several validation scores
#' @details
#' * if 'skipCheck' = TRUE the input checking is skipped and the itrations starts directly.
#' * has implemneted error catching. if an iteration would cause critical stop error, the loops continue. Returns NA for corrupted iterations
#' * Input for a,b,h,MIN,filter supports - numeric, single, combination by c() or sequence by seq() to iterate over.
#' * Validation Point - supervised computed pointlayer. Use a GIS like Qgis to set points where Trees are estimated.
#'
#' * ! even if the function dupport iteration over all parameters (except MAX due to estimated no nneed for MAX Crownareas to clip) it should not be used to iterate over seqences which would be useless.
#' @author Andreas Schönberg
#' @note This function is experimental. The 'Brute Force' abbroach to iterate over many parameters may result in very long time to finsih. Preselected samller sample may be more efficent
#' @examples
#' # load data
#' chmpath <-system.file("extdata","exp_chm.tif",package = "CENITH")
#' chm <- raster::raster(chmpath)
#' plot(chm)
#' vppath <-system.file("extdata","exp_vp.shp",package = "CENITH")
#' vp <- rgdal::readOGR(vppath)
#' # lab test if doParalel works
#' # make sure no open cluster are running
#' stopCluster()
#' # 1st test without doPar
#' start1 <- Sys.time()
#' np <-NPBestSegValEXPRMTL(chm,a=c(0.1,0.3,0.5),b=c(0.1,0.3,0.5),h=0.5,vp,MIN=c(10,30,50),MAX=1000,filter=c(1,3,5))
#' stop1 <- Sys.time()
#' # 2nd test with doPar full cores
#' start2 <- Sys.time()
#' fc <-BestSegValEXPRMTL(chm,a=c(0.1,0.3,0.5),b=c(0.1,0.3,0.5),h=0.5,vp,MIN=c(10,30,50),MAX=1000,filter=c(1,3,5),Cores=1)
#' stop2 <- Sys.time()
#' # 3rd test with doPar but lesser cores (3)
#' start3 <- Sys.time()
#' x1 <-BestSegValEXPRMTL(chm,a=c(0.1,0.3,0.5),b=0.5,h=0.5,vp,MIN=c(10,30,50),MAX=1000,filter=c(1,3,5),Cores=5)
#' stop3 <- Sys.time()
#' # check time differences
#' difftime(start1,stop1) # no Par
#' difftime(start2,stop2) # full cores
#' difftime(start3,stop3) # 3 cores (if 8 availibe)
#' @export NPBestSegValEXPRMTL
#' @aliases NPBestSegValEXPRMTL


NPBestSegValEXPRMTL<- function(chm,a,b,h,vp,MIN=0,MAX=1000,skipCheck=FALSE,filter=NULL){
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
    nm <- length(MIN)
    niter <- nh*nb*na*nm

    # warning / wait for input if niter is too big
    if(niter>50){
      cat("",sep = "\n")
      cat(paste0("WARNING over 50 iterations detected this may take long time to finsih, maybe abort and start smaller trys (see details)",sep = "\n"))
    }
    if(niter>100){
      abort <- readline(prompt="More than 100 iterations detected! Are you sure you want to continue? 1:YES, 2:NO, plz abort ")
      if(abort==2){
        stop("aborded on users desicion, too much iterations detected")
      }
    }

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
  loop_b <- function(chm,a,b,h,vp,MIN,MAX){


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
      cat("",sep = "\n")
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
      cat("",sep = "\n")
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
  cat(paste0("### Cenith finsihed Segmentation ###"),sep = "\n")

  names(res)<- c("a","b","height","total_seg","hit/vp","under","over","area","hitrate","underrate","overrate","Seg_qualy","MIN","chm")
  #names(res)<- c("a","b","height","hit","tp/vp_rate","tpos/vp","miss","area","empty")
  return(res)
}# end core fucntion

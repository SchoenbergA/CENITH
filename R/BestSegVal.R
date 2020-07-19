#' Detect best input variables for TreeSegmentation by validation points
#' @description uses supervised computed Treepositions to validate best fitting values for a,b and height.
#' @param chm raster - RasterLayer with Canopy height model
#' @param a numeric - function for MovingWindow
#' @param b numeric - function for MovingWindow
#' @param h numeric - maximum height of Trees
#' @param vp Polygon - PointLayer with estimated Positions of Trees.
#' @return returns a dataframe with several validation scores
#' @details
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
#' x <-BestSegVal(chm,seq(0.7,0.8,0.1),seq(0.5,0.7,0.1),c(13,8),vp)
#' # view best accuracy (diffenrent iterations can have the same hit rate, but differ in other values)
#' maxrow <- x[which.max(x$hit),] # search max vale but rturn only 1 value
#' maxhit <- maxrow$hit
#' x[which(x$hit==maxhit),]

#' @export BestSegVal
#' @aliases BestSegVal


BestSegVal<- function(chm,a,b,h,vp){
  cat(paste0("### Cenith checking input ###",sep = "\n"))
  # check input height is lesser than max height in chm
  if(cellStats(chm,max)<=max(h)){
    stop("input height 'h' values are higher than the highest cell value of teh CHM")
  }
  # check if vlaues in a and or b are to big
  maxTest <-try(ForestTools::vwf(chm,
                              winFun = function(x){x * max(a) + max(b)},
                              minHeight = max(h),
                              verbose = TRUE),silent = TRUE)
  if(class(maxTest)=="try-error"){
    stop("any input values for 'a' and or 'b' lead to wider window diameter. Reduce Values! ")
  }
  cat("",sep = "\n")
  cat(paste0("### Cenith calculates estimate time to finish ###",sep = "\n"))
  cat("",sep = "\n")
  # estimate time remaining
  # amount of iterations
  na <- length(a)
  nb <- length(b)
  nh <- length(h)
  niter <- nh*nb*na
  # test time for 1 iteration
  start <- Sys.time()
  timecheck <- TreeSeg(chm=chm,a[1],b[1],h[1])
  stop  <- Sys.time()
  timedif <-difftime(stop,start,units = "min")
  # predict estimated time by multiply timedif with n iterations

  eta <- timedif*niter
  eta <-1.2*eta # add 20% because the iteration wrapper needs more time.
  cat(paste("requires",niter,"iterations @",round(timedif,4),units(timedif),"per iteration"),sep="\n")
  cat("",sep = "\n")
  cat(paste("estimated",round(eta,4),units(timedif),"to finish"))
  cat("",sep = "\n")


                                  # setup function to iterate over b ############################################################
                                  loop_b <- function(chm,a,b,h,vp){


                                  # create dataframe to save informations
                                  result <- data.frame(matrix(nrow = length(b), ncol = 8)) # ncol = n information stored

                                  # loop to iterate on varibale
                                  for (j in seq(1:length(b))){

                                    cat(paste0("iterating over b ",as.factor(j),"/",as.factor(length(b))," ",sep = "\n"))
                                  # TreeSeg
                                  seg <- TreeSeg(chm=chm,a,b[j],h)

                                  # get stats of overlapping
                                  stat <- ForestTools::sp_summarise(vp, seg) # compute data points in polygons
                                  stat[is.na(stat$TreeCount)] <- 0 # na to 0


                                  pkb <- sum(stat$TreeCount<1) # amount polygon without tree (miss)
                                  pb <- sum(stat$TreeCount==1) # amount polygon with exact 1 tree (hit)
                                  pmb <-sum(stat$TreeCount>1) # amount polygon with more than 1 tree (miss)

                                  hit = pb/length(stat$TreeCount) # calc hit ration in percent (amount of exact trees
                                  over = pkb/length(stat$TreeCount) #calc empty ration in percent (amount of polygon without trees)
                                  under = pmb/length(stat$TreeCount) # mis.rati (or jan error) miss rate in percent (amount of polygons with more than 1 Tree)
                                  ntree_vp = length(seg)/length(vp)
                                  area =  sum(seg$crownArea)

                                  # write out informations in dataframe
                                  result[j, 1] <- a
                                  result[j, 2] <- b[j]
                                  result[j, 3] <- h
                                  result[j, 4] <- hit
                                  result[j, 5] <- ntree_vp
                                  result[j, 6] <- paste(length(seg),"/",length(vp))
                                  result[j, 7] <- over
                                  result[j, 8] <- area
                                  result[j, 9] <- under

                                  }
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

names(res)<- c("a","b","height","hit","tp/vp_rate","tpos/vp","miss","area","empty")
return(res)
}# end core fucntion

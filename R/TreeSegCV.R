#' TreeSegCV
#' @description cs
#' @param sides list - chm rasterlayer (see details)
#' @param a numeric - function for MovingWindow
#' @param b numeric - function for MovingWindow
#' @param h numeric - maximum height of Trees
#' @param vps list - of PointLayers with estimated Positions of Trees (see details)
#' @param MIN numeric - minimum area for Crowns. smaller poylgons are cropped
#' @param MAX numeric - maximum area for Crowns. larger polygons are cropped
#' @param skipCheck development - bolean - if TRUE skips chekcing the inputs
#' @return returns a dataframe with quality values for the folds and the mean values for a estimted mean quality
#' @details
#' *'sides' and 'vps' must be lists with same order of chm and respectiv validtion points
#' * a,b,h,MIN,MAX values should be calculated by 'TreeSegVal'
#' @author Andreas Schönberg
#' @example
#' # load chm sides and vp layers
#' # list all chm and vp in seperated list
#' sides <-list(chm1,chm2,chm3)
#' vps <-list(vp1,vp2,vp3)
#' # start CV with model (a,b,h,MIN,MAX Values)
#' x <- TreeSegCV(sides,a,b,h,MIN,MAX,vps)


#' @export TreeSegCV
#' @aliases TreeSegCV

### tests
# test cat codes with lists

TreeSegCV <- function(sides,a,b,h,MIN,MAX,vps){
  # cheking inputs

  # create dataframe to save informations
  cat(paste0("### Cenith starts ",length(sides),"fold cross validation ###",sep = "\n"))
  result <- data.frame(matrix(nrow = length(b), ncol = 8)) # ncol = n information stored
  # iteration
  # loop to iterate on varibale
  for (i in seq(1:length(sides))){

    cat(paste0("starting fold ",as.factor(i),"/",as.factor(length(sides))," ",sep = "\n"))
    # TreeSeg
    seg <- try(TreeSeg(chm=sides[[i]],a,b,h,MIN,MAX),silent = TRUE)
    # handle error catch
    if(class(seg)=="try-error"){
      cat(paste0(" !!! iteration a=",a," b=",b," h=",h," leads to an error on side:", side[[i]], "setting to 'NA' in results ",sep="\n"))
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

      # additional informations
      tseg = length(seg)
      area =  sum(seg$crownArea)#

      # write out informations in dataframe
      result[i, 1] <- names(sides[[i]])
      result[i, 2] <- a
      result[i, 3] <- b
      result[i, 4] <- h
      #absolut results
      result[i, 5] <- tseg
      result[i, 6] <- hit
      result[i, 7] <- TCU
      result[i, 8] <- TCO
      result[i, 9] <- area
      # rates
      result[i, 10] <- hitrate
      result[i, 11] <- under
      result[i, 12] <- over
      result[i, 13] <- segQy
    } # end of more than null polygons
            # handle df
            if(i==1){
              res <-result

            }    else {
              res2 <-result
              res= rbind(res,res2)
            }

  }# end of iteration

  # calc means
  res
  nrow(res)+1
  # write out informations in dataframe
  res[nrow(res)+1, 1] <- "Mean"
  res[nrow(res)+1, 2] <- a
  res[nrow(res)+1, 3] <- b
  res[nrow(res)+1, 4] <- h
  #abolut results
  res[nrow(res)+1, 5] <- mean(res$tseg)
  res[nrow(res)+1, 6] <- mean(res$hit)
  res[nrow(res)+1, 7] <- mean(res$TCU)
  res[nrow(res)+1, 8] <- mean(res$TCO)
  res[nrow(res)+1, 9] <- mean(res$area)
  # rates
  res[nrow(res)+1, 10]  <-  mean(res$hitrate)
  res[nrow(res)+1, 11] <-  mean(res$under)
  res[nrow(res)+1, 12] <-  mean(res$over)
  res[nrow(res)+1, 13] <- (mean(res$TCO)+ 2* mean(res$TCU))/2
  names(res)<- c("sides","a","b","height","total_seg","hit/vp","under","over","area","hitrate","underrate","overrate","Seg_qualy")
  return(res)
}# end of function



# Title: "MapSpam Monfreda Code Geoshpere produce CCRI raster 0.5 degree Sweetpotato"
#author: "Y Xing"
#date: "09/07/2020"
#type: functions

#------- 4.  Function : generate data frame with GIS info. and values of harvest cropland fraction -----
#---------- 4.1 Function to generate data frame with GIS info. and values of harvest cropland fraction -----
FuncGIScroplandW <- function(cropRasterW, HCFcutoff){
                    cropRasterWValues <- getValues(cropRasterW)
                    WcellNum <- which(cropRasterWValues > HCFcutoff)
                    cellNumValueW <- cropRasterWValues[WcellNum]
                    Wlength <- length(cellNumValueW)
                    colnumsW <- ncol(cropRasterW)
                    latiW <- c(rep(0, Wlength))
                    logiW <- c(rep(0, Wlength))
                    # looping through file 
                             for (k in c(1:Wlength)) {
                               me<-(WcellNum[k] - 1) %/% colnumsW
                               latiW[k] <- latitoW - me * CellDegree
                               logiW[k] <- (WcellNum[k] - 1 - me * colnumsW) * CellDegree + longifromW
                              }    
                    cropdataW <- cbind(latiW, logiW, cellNumValueW)
                    return(as.data.frame(cropdataW))
}
#---
#---
#-------  5. Inverse Power Law Model -----------
#---
#---
networkbetaW <- function(cropdataW, beta3, cutoffadja3){   # added cropdataW 
                # create adjacency matrix
                rownumber13 <- nrow(cropdataW)
                latilongimatr3 <- cropdataW[1:rownumber13, c(2,1)]# save the latitude and longitude as new matrix
                # Geosphere package 
                # function distVincentyEllipsoid() is used to calculate the distance, defult distance is meter
                dvse <- distVincentyEllipsoid(c(0,0), cbind(1, 0)) # reference of standard distance in meter for one degree
                #dvse is the constant distance of 1 degree fron the Ecuator to the Poles
                # The constant is the (meters) distance for one degree at the Ecuator
                latilongimatr3 <- as.matrix(latilongimatr3)
                TemMat <- matrix(-999, nrow( latilongimatr3), nrow(latilongimatr3))
                # looping thrhoug lat-lon data
                      for (i in 1:nrow(latilongimatr3)) {
                      TemMat[i, ] <- distVincentyEllipsoid(latilongimatr3[i,], latilongimatr3) #/dvse # Divided the Vicenty by the meters at the Ecuator
                      }                                                                       # commenting this would change from degree to meters
               # check the value, if is not 0 [closer to 0] use this value as a denominator TemMat/min(TemMat)
                distancematr3 <- TemMat/as.numeric(names(table(TemMat)[2])) #but not 0 [using some logical]  # This is the Euclidean distance Vicenty ellipsoid = Geodesic
                # Distance in degre lon-lat° 
                distancematrexp3 <- distancematr3^(-beta3)  #use function C=AX^(-beta), here A=1, X=distancematr3
                # Distance matrix in degree, apply inverse- powerlog (beta values)
                # New distance: Proportional to the prob. of spread 
                cropmatr3 <- cropdataW[1:rownumber13, 3] # complete gravity model with crop data
                cropmatr13 <- matrix(cropmatr3, , 1)
                cropmatr23 <- matrix(cropmatr3, 1, )
                cropmatrix3 <- cropmatr13 %*% cropmatr23
                cropmatrix3 <- as.matrix(cropmatrix3)
                cropdistancematr3 <- distancematrexp3 * cropmatrix3 # dij^-beta * cicj
                #  
                logicalmatr3 <- cropdistancematr3 > cutoffadja3
                stan3 <- cropdistancematr3 * logicalmatr3
                cropdistancematrix3 <- graph.adjacency(stan3, mode = c("undirected"), 
                                         diag = FALSE, weighted = TRUE) #change the threshhold to see the difference
                #
                # create network for all the selected nodes
                #
                V(cropdistancematrix3)$label.cex = 0.7
                E(cropdistancematrix3)$color = "red"
                edgeweight3 <- E(cropdistancematrix3)$weight * 10000 
                knnpref0 <- graph.knn(cropdistancematrix3,weights = NA)$knn
                knnpref0[is.na(knnpref0)] <- 0
                degreematr <- degree(cropdistancematrix3)
                knnpref <-knnpref0 * degreematr
                      if(max(knnpref) == 0){knnprefp = 0}else
                        if(max(knnpref) > 0){knnprefp = knnpref / max(knnpref) / 6}
                #
                ####  node degree, node strengh
                #
                nodestrength <- graph.strength(cropdistancematrix3)
                nodestrength[is.na(nodestrength)] <- 0
                      if(max(nodestrength) == 0){nodestr = 0}else
                        if(max(nodestrength) > 0){nodestr = nodestrength / max(nodestrength) / 6}
                #
                ####   betweenness centrality
                #    
                between <- betweenness(cropdistancematrix3)
                between[is.na(between)] <- 0
                      if(max(between) == 0){betweenp = 0}else
                        if(max(between) > 0){betweenp = between / max(between) / 2}
                #
                ####   eigenvector and eigenvalues
                #
                eigenvectorvalues <- evcent(cropdistancematrix3)
                ev <- eigenvectorvalues$vector
                ev[is.na(ev)] <- 0
                    if(max(ev) == 0){evp = 0}else
                        if(max(ev) != 0){evp = ev / max(ev) / 6}
                #
                ####   plot index layer
                #    
                index <- knnprefp + evp + betweenp + nodestr
                indexpre <- cropharvestRasterWaggValues
                indexpre[] <- 0
                indexpre[cellNumW] <- index
                indexv3 <- indexpre
      # notes:                       Eucl       power-law         product between each pair or   
      return(list(indexv3, stan3, distancematr3, distancematrexp3, cropmatrix3))
}
#---
#---
#-------  6. Negative Exponential Model -----------
#---
#---
networkgammaW <- function(cropdataW, gamma3, cutoffadja3){ # Negative Exponential Model 
                  #
                  #### create adjacency matrix
                  #
                  rownumber13 <- nrow(cropdataW)
                  latilongimatr3 <- cropdataW[1:rownumber13, c(2,1)]# save the latitude and longitude as new matrix  #---- use Geosphere package, function distVincentyEllipsoid() is used to calculate the distance, defult distance is meter
                  dvse <- distVincentyEllipsoid(c(0,0), cbind(1, 0)) # reference of standard distance in meter for one degree
                  #
                  latilongimatr3 <- as.matrix(latilongimatr3)
                  TemMat <- matrix(-999, nrow(latilongimatr3), nrow(latilongimatr3))
                  #  
                            for (i in 1:nrow(latilongimatr3)) {
                                TemMat[i, ] <- distVincentyEllipsoid(latilongimatr3[i,], latilongimatr3) #/ dvse
                            }                                                                            # commenting this would change from degree to meters
                  distancematr3 <- TemMat/as.numeric(names(table(TemMat)[2])) #but not 0 [using some logical]  # This is the Euclidean distance Vicenty ellipsoid = Geodesic
                
                  # distancematr3 <- TemMat
                  #
                  #----
                  #
                  eulernumber <- exp(1)
                  distancematrexponential3 <- eulernumber ^ (-gamma3 * distancematr3)# exponential model
                  
                  cropmatr3 <- cropdataW[1:rownumber13, 3] # complete gravity model with crop data
                  cropmatr13 <- matrix(cropmatr3,, 1)
                  cropmatr23 <- matrix(cropmatr3, 1, )
                  cropmatrix3 <- cropmatr13 %*% cropmatr23
                  cropmatrix3 <- as.matrix(cropmatrix3)
                  cropdistancematr3 <- distancematrexponential3 * cropmatrix3
                  logicalmatr3 <- cropdistancematr3 > cutoffadja3
                  stan3 <- cropdistancematr3 * logicalmatr3
                  cropdistancematrix3 <- graph.adjacency(stan3, mode=c("undirected"), 
                                                         diag = FALSE, weighted = TRUE)#change the thresh to see the difference
                  #
                  #### create network for all the selected nodes
                  #
                  V(cropdistancematrix3)$label.cex = 0.7
                  E(cropdistancematrix3)$color = "red"
                  edgeweight3 <- E(cropdistancematrix3)$weight * 10000 
                  knnpref0 <- graph.knn(cropdistancematrix3, weights = NA)$knn
                  knnpref0[is.na(knnpref0)] <- 0
                  degreematr <- degree(cropdistancematrix3)
                  knnpref <- knnpref0 * degreematr
                      if(max(knnpref) == 0) {knnprefp = 0} else
                        if(max(knnpref) > 0) {knnprefp = knnpref / max(knnpref) / 6}
                  #
                  ####  node degree, node strengh
                  #
                  nodestrength <- graph.strength(cropdistancematrix3)
                  nodestrength[is.na(nodestrength)] <- 0
                      if(max(nodestrength) == 0) {nodestr = 0} else
                        if(max(nodestrength) > 0){nodestr = nodestrength / max(nodestrength) / 6}
                  #
                  ####  betweenness centrality
                  #    
                  between<-betweenness(cropdistancematrix3)
                  between[is.na(between)] <- 0
                      if(max(between) == 0) {betweenp = 0} else
                        if(max(between) > 0) {betweenp = between / max(between) / 2}
                  #
                  ####   eigenvector and eigenvalues
                  #
                  eigenvectorvalues <- evcent(cropdistancematrix3)
                  ev <- eigenvectorvalues$vector
                  ev[is.na(ev)] <- 0
                      if(max(ev) == 0) {evp = 0} else
                        if(max(ev) != 0){evp = ev / max(ev) / 6}
                  #
                  ####   plot index layer
                  #    
                  index <- knnprefp + evp + betweenp + nodestr
                  indexpre <- cropharvestRasterWaggValues
                  indexpre[] <- 0
                  indexpre[cellNumW] <- index
                  indexv3 <- indexpre
        # notes:                       Eucl       power-law         product between each pair or  
        return(list(indexv3, stan3, distancematr3, distancematrexponential3, cropmatrix3))
}
#---
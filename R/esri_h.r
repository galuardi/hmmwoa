esri_h <- function(rdf){
  # uses function for bandwidth as suggested by ESRI at:
  # from ESRI: http://pro.arcgis.com/en/pro-app/tool-reference/spatial-analyst/how-kernel-density-works.htm

  sdd <- calc_sdd(calccentre=TRUE, weighted=TRUE, weights=rdf[,1],
                  points=rdf[,c(2,3)], verbose=F)$SDD.radius
  
  mc <- mean_centre(weighted=TRUE, weights=rdf[,1], points=rdf[,c(2,3)])
  Dm <- median(distances(centre.xy = mc[,2:3], destmat = rdf[,c(2,3)]))
  
  term2 <- sqrt(1 / log(2)) * Dm
  
  #n <- cellStats(r, 'sum')
  n <- sum(weights)
  
  if(sdd < term2){
    h <- .9 * sdd * n^(-0.2)
  } else{
    h <- .9 * term2 * n^(-0.2)
  }
  
  h
  
}

calc_sdd <- function (id = 1, centre.xy = NULL, 
                      calccentre = TRUE, weighted = FALSE, weights = NULL, points = activities, 
                      verbose = FALSE) 
{
  errorcode <- 1000
  n <- dim(points)[1]
  if (calccentre) {
    if (length(centre.xy) == 2) {
      errorcode <- 21
      cat("\n\nWARNING: Invalid combination: calccentre=TRUE and centre.xy!=NULL")
      cat("\nERROR CODE: ", errorcode, "\n\n", sep = "")
      return("ERROR")
    }
    else {
      if (weighted) {
        wt.x <- points[, 1] * weights
        wt.y <- points[, 2] * weights
        WMC.x <- c(sum(wt.x)/sum(weights))
        WMC.y <- c(sum(wt.y)/sum(weights))
        centre.xy[1] <- WMC.x
        centre.xy[2] <- WMC.y
      }
      else {
        meanx <- sum(points[, 1])/n
        meany <- sum(points[, 2])/n
        centre.xy[1] <- meanx
        centre.xy[2] <- meany
      }
    }
  }
  dist <- distances(centre.xy, points)
  if (length(dist) >= 3) {
    if (weighted) {
      SDD <- sqrt(sum((weights * dist^2)/((sum(weights)) - 2)))
    }
    else {
      SDD <- sqrt(sum(dist^2/(length(dist) - 2)))
    }
    sddarea <- pi * SDD^2
    B <- min(SDD, SDD)
    A <- max(SDD, SDD)
    d2 <- (A - B) * (A + B)
    phi <- 2 * pi * seq(0, 1, len = 360)
    sp <- sin(phi)
    cp <- cos(phi)
    r <- SDD * SDD/sqrt(B^2 + d2 * sp^2)
    xy <- r * cbind(cp, sp)
    al <- 0 * pi/180
    ca <- cos(al)
    sa <- sin(al)
    coordsSDD <- xy %*% rbind(c(ca, sa), c(-sa, ca)) + cbind(rep(centre.xy[1], 
                                                                 360), rep(centre.xy[2], 360))
    sddloc <- as.data.frame(cbind(id, coordsSDD))
    colnames(sddloc) = c("id", "x", "y")
    #write.table(sddloc, sep = ",", file = filename, col.names = FALSE)
    assign("sddloc", sddloc, pos = 1)
    r.SDD <- list(id = id, points = points, coordsSDD = coordsSDD, 
                  SDD = SDD, calccentre = calccentre, weighted = weighted, 
                  weights = weights, CENTRE.x = centre.xy[1], CENTRE.y = centre.xy[2], 
                  SDD.area = sddarea)
    assign("r.SDD", r.SDD, pos = 1)
    result.sdd <- list(id = id, calccentre = calccentre, 
                       weighted = weighted, CENTRE.x = centre.xy[1], CENTRE.y = centre.xy[2], 
                       SDD.radius = SDD, SDD.area = sddarea)
    print(result.sdd)
    result.sdd <- as.data.frame(result.sdd)
    assign("sddatt", result.sdd, pos = 1)
  }
  else {
    errorcode <- 25
    if (verbose) {
      cat("\n\nWARNING: Not enough values to compute SDD.")
      cat("\nERROR CODE: ", errorcode, "\n\n", sep = "")
    }
    return("ERROR")
  }
}

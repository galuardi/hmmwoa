esri_h <- function(rdf){
  # uses function for bandwidth as suggested by ESRI at:
  # from ESRI: http://pro.arcgis.com/en/pro-app/tool-reference/spatial-analyst/how-kernel-density-works.htm

  # seems like units degrees
  sdd <- calc_sdd(calccentre=TRUE, weighted=TRUE, weights=rdf[,1],
                  points=rdf[,c(2,3)], verbose=F)$SDD.radius
  
  # this is working and yielding weighted center
  mc <- mean_centre(weighted=TRUE, weights=rdf[,1], points=rdf[,c(2,3)])
  
  # also in degrees
  Dm <- median(distances(centre.xy = as.numeric(mc[,2:3]), destmat = rdf[,c(2,3)]))
  
  term2 <- sqrt(1 / log(2)) * Dm
  
  # n, scaled by likelihood values
  n <- sum(rdf[,1])
  
  if(sdd < term2){
    h <- .9 * sdd * n^(-0.2)
  } else{
    h <- .9 * term2 * n^(-0.2)
  }
  
  h
  
}

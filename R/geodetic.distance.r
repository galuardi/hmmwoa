geodetic.distance <- function(point1, point2){
  # function from: http://www.biostat.umn.edu/~sudiptob/Software/distonearth.R
  #dist between 2 points on surface of earth; 
  #points are of the form long,lat; results in kilometers
  
  R <- 6371
  p1rad <- point1 * pi/180
  p2rad <- point2 * pi/180
  d <- sin(p1rad[2])*sin(p2rad[2])+cos(p1rad[2])*cos(p2rad[2])*cos(abs(p1rad[1]-p2rad[1]))
  d <- acos(d)
  R*d
}

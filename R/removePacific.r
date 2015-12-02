removePacific <- function(dat, lat, lon){

  #' script to remove pacific from wod data to prevent matching there
  #' @param dat is output from extract.woa
  #' @param lat is output from extract.woa
  #' @param lon is output from extract.woa
  #' @return dat is WOA data grid with Pacific removed
  #' only tested when area of interest in N Atlantic
  
  
  # run 1
  ylims = c(min(lat), 5)
  xlims = c(min(lon), -60)
  
  ymin = which.min((ylims[1] - lat) ^ 2)
  ymax = which.min((ylims[2] - lat) ^ 2)
  xmin = which.min((xlims[1] - lon) ^ 2)
  xmax = which.min((xlims[2] - lon) ^ 2)
  
  dat[xmin:xmax, ymin:ymax,,] = NA
  
  # run 2
  ylims = c(5, 9)
  xlims = c(min(lon), -70)
  
  ymin = which.min((ylims[1] - lat) ^ 2)
  ymax = which.min((ylims[2] - lat) ^ 2)
  xmin = which.min((xlims[1] - lon) ^ 2)
  xmax = which.min((xlims[2] - lon) ^ 2)
  
  dat[xmin:xmax, ymin:ymax,,] = NA
  
  # run 3
  ylims = c(9, 15)
  xlims = c(min(lon), -84)
  
  ymin = which.min((ylims[1] - lat) ^ 2)
  ymax = which.min((ylims[2] - lat) ^ 2)
  xmin = which.min((xlims[1] - lon) ^ 2)
  xmax = which.min((xlims[2] - lon) ^ 2)
  
  dat[xmin:xmax, ymin:ymax,,] = NA
  
  # run 4
  ylims = c(15, 18)
  xlims = c(min(lon), -90)
  
  ymin = which.min((ylims[1] - lat) ^ 2)
  ymax = which.min((ylims[2] - lat) ^ 2)
  xmin = which.min((xlims[1] - lon) ^ 2)
  xmax = which.min((xlims[2] - lon) ^ 2)
  
  dat[xmin:xmax, ymin:ymax,,] = NA
  
  # run 5
  ylims = c(9, 15)
  xlims = c(min(lon), -84)
  
  ymin = which.min((ylims[1] - lat) ^ 2)
  ymax = which.min((ylims[2] - lat) ^ 2)
  xmin = which.min((xlims[1] - lon) ^ 2)
  xmax = which.min((xlims[2] - lon) ^ 2)
  
  dat[xmin:xmax, ymin:ymax,,] = NA
  
  return(dat)
  
}


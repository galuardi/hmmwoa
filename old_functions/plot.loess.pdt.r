
plot.loess.pdt <- function(days, depth, temperature, year, tempRange = c(0, 32),writePDF = TRUE,
                           filename, DOY.first, zissou = FALSE){

  #' function to plot depth temperature profiles output from grid2dloess.r
  #' depth is plotted through time and data colored by temperature at that depth
  #'
  #' @param DAYS is a row vector of days to plot along x axis (must match length of columns in matrix TEMPERATURE)
  #' @param DEPTH is a row vector of depths to plot along y axis (must match length of ROWS in matrix TEMPERATURE)
  #' @param TEMPERATURE is a DEPTH X DAYS matrix of temperature values interpolated using grid2dloess.r
  #' @param YEAR is row vector of years corresponding to days on x axis (meant only for checking for leap year and labelling accordingly)
  #' @param TEMPRANGE is temperature range of input matrix TEMPERATURE
  #' @param WRITEPDF is argument to determine whether the plot should be output to file as .pdf or just displayed
  #'     this one is not yet functional so keep default of writing out
  #' @param FILENAME is character string to name your output .pdf
  #' @param DOY.FIRST is integer vector of length one indicating first day of year in the data to ensure correct date labeling
  #' @param ZISSOU is logical to determine whether plot is created with zissou color ramp and fonts
  #'     - currently this is done with a zissou color ramp from https://github.com/karthik/wesanderson#wes-anderson-palettes
  #'
  #' @return PDF file of plot as described above
  #'
  #'
  #' Revisions Needed:
  #' 1) add functionality for handling datasets that span > 1 calendar year (e.g. Nov->Feb)
  #' 2) allow use of writePDF = FALSE
  #' 3) add ability to use flag output from grid2dloess() within this function
  #' 4) handle working directory/output file directory
 
  #require(oce);
  
  # setup some plotting variables
  leap.year = seq(1904, 2104, by = 4)
  monthLabels = c('1Jan', '1Feb', '1Mar', '1Apr', '1May', '1Jun', '1Jul', '1Aug', 
                  '1Sep', '1Oct', '1Nov', '1Dec')
  if(any(year == leap.year)){
    axis.days = c(1, 32, 61, 92, 122, 153, 183, 214, 245, 275, 306, 336) #labels for x axis, leap year
  } else{
    axis.days = c(1, 32, 60, 91, 121, 152, 182, 213, 244, 274, 305, 335) #labels for x axis, non-leap year
  }
  zr = tempRange
  
  # write PDF?
  if(writePDF){
    pdf(file = paste(filename, '.pdf', sep = ''), width = 11, height = 8) #write pdf file
  } else{}
  
  # start the plot
  if(zissou){
    Zissou <- c("#3B9AB2", "#78B7C5", "#EBCC2A", "#E1AF00", "#F21A00")
    par(mfrow = c(1, 1), mar = c(3, 5, 3, 5), family = 'Helvetica') #set margins so have space for colorbar legend
    zissouCol = colorRampPalette(Zissou)
    image(days, depth, axes = F, xlab = '', ylab = 'Depth (m)', t(temperature), #draw the image
          ylim = c(max(z), min(z) - 1), zlim = zr, col = zissouCol(200))
    axis(1, at = axis.days - min(DOY.first, na.rm = T), labels = monthLabels) #x axis
    axis(2) #y axis
    contour(days, depth, t(temperature), add = T, nlevels = 20, zlim = zr) #add contours, uses derivation of pretty() to determine intervals
    box(lwd = 1)
    mtext('Water Temp (C)', side = 4, line = 3.5)
    image.plot(legend.only = T, zlim = zr,col = zissouCol(200))
  } else{
    par(mfrow = c(1, 1), mar = c(3, 5, 3, 5)) #set margins so have space for colorbar legend
    image(days, depth, axes = F, xlab = '', ylab = 'Depth (m)', t(temperature), #draw the image
          ylim = c(max(z), min(z) - 1), zlim = zr, col = oceColorsJet(200))
    axis(1, at = axis.days - min(DOY.first), labels = monthLabels) #x axis
    axis(2) #y axis
    contour(days, depth, t(temperature), add = T, nlevels = 20, zlim = zr) #add contours, uses derivation of pretty() to determine intervals
    box(lwd = 1)
    mtext('Water Temp (C)', side = 4, line = 3.5)
    image.plot(legend.only = T, zlim = zr)
  }
  
try(dev.off(), silent = T) #close written pdf

  }

# run it

# load and format tag data
# read in data
setwd('~/Documents/WHOI/RData/WhiteSharks/2013/121325/')
data = read.table('121325-PDTs.csv',sep=',',header=T,blank.lines.skip=F)

pdt = extract.pdt()

# download daily hycom products
get.hycom(lon,lat,time,filename=paste()) # filenames based on dates from above

# calc.ohc
likelihood = calc.ohc()



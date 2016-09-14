####
#
# Basic plot of AXCTD
#
####

library(oce) #oceanographic data/plotting package downloadable through CRAN

myfile = file.choose()  #user menu selection
axctd_columns = c('Frame', 'Data', 'CRC','Depth','Temp','Cond',"Salinity") #fixed column names
axctd.data <- read.table(myfile,skip=6, header=FALSE, na.strings='*****', col.names = axctd_columns)

df <- as.ctd(axctd.data$Salinity,axctd.data$Temp,axctd.data$Depth,axctd.data$Cond)

###setup plot
save_file = gsub('.dta','.png',myfile)
png(filename=save_file, 
    units="in", 
    width=5, 
    height=4, 
    pointsize=12, 
    res=300)

par(oma=c(0,0,2,0))                                           # outer margins: thick at top, or side=3
plot(df, which=1, Slim=c(20,35), Tlim=c(-4,10), plim=c(50,0), keepNA=TRUE)  # set salinity and temperature bounds and make missing values data gaps
title(basename(myfile), outer=TRUE, cex=0.4) # set title based on filename


dev.off()
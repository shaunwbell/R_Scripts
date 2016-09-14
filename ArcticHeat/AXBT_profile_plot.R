####
#
# Basic plot of AXBT
#
####

library(ggplot2) #extended plotting package available from CRAN

#path to axbt file
myfile = file.choose()
axbt.data <- read.table(myfile,skip=4, header=TRUE, na.strings='******') #read data


###setup plot
save_file = gsub('.dta','.png',myfile) #substitutes png for dta file ending
png(filename=save_file, 
    units="in", 
    width=5, 
    height=4, 
    pointsize=8, 
    res=300)

### Using ggplot2
p <- ggplot(axbt.data, aes(x=X.C.,y=Depth))
p <- p + geom_point(color='red', cex=0.25)
p <- p + scale_y_continuous(trans='reverse', limits=c(NaN,0))
p <- p + ggtitle(basename(myfile))
p

dev.off()
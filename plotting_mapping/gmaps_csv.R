library(ggmap)

##load data
df.aeis <- read.csv('/Users/bell/Programs/Python/AtSeaPrograms/PostCruiseRoutines/2013_twolayer_integratedtemp.csv')

# ## using google maps basemap
# myLocation <- c(lon = -160.0, lat = 65)
# myMap <- get_map(location=myLocation, source="google", maptype="satellite", zoom=4, crop=FALSE)
# ggmap(myMap, extent="normal", legend="bottomright")
# 
# ggmap(myMap) + geom_point(data=df.aeis, aes(x=-1*df.aeis$lon..W., y=df.aeis$lat..N., size=df.aeis$IntegratedTemp_top15m, color=df.aeis$IntegratedTemp_top15m)) +
#   scale_color_gradient(low="#fc9272", high="#67000d", name="Upper 15m Integrated Temperature") 

## using stamen maps basemap
# center lat/lon point will prevent map from loading if too far west
myLocation <- c(lon = -150, lat = 65)
myMap <- get_map(location=myLocation, source="stamen", maptype="toner", zoom=4)
ggmap(myMap, extent="normal", legend="bottomright")

ggmap(myMap) + geom_point(data=df.aeis, aes(x=-1*df.aeis$lon..W., y=df.aeis$lat..N., size=df.aeis$IntegratedTemp_top15m, color=df.aeis$IntegratedTemp_top15m)) +
  scale_color_gradient(low="#fc9272", high="#67000d", name="Upper 15m Integrated Temperature") 

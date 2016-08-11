#SAVING WORKSPACE
#save.image("C:/MyWorkSpace.RData")
#load("C:/MyWorkSpace.RData")

#LIBRARIES
?gam
library(mgcv)

library(help=mgcv)

#OBJECTS
ls()

#Vectors: numeric
a<-c(2,3,4,5)
a

b<-seq(2,5)
b

b<-c(2:5)
b

ab<-c(a,b)
ab

length(ab)

#Vectors: character
f<-c('apple','orange','grapes')
f

#Vectors: logical
l<-ab>=3
l
ab[l]

#Class of vectors
class(ab)
class(f)
class(l)

#Sum, length and remove
sum(ab)
length(ab)
rm(a,b)

#Matrices
a<-matrix(1:8, nrow=2,ncol=4)
a
a<-matrix(1:8, nrow=2,ncol=4,byrow=T)
a

dim(a)
sum(a)

#List
l<-list(ab=ab,f=f,a=a)
l

length(l)
unlist(l)
class(unlist(l))

#Data Frames
d<-data.frame(year=1991:2000,temp=rnorm(n=10,mean=6,sd=3),cpue=runif(n=10,min=2,max=3))
d

names(d)

names(d)<-c('col1','col2','col3')
names(d)

names(d)<-c('year','temp','cpue')
names(d)

dim(d)

summary(d)

#FUNCTIONS
mean(ab)

a<-matrix(1:8, nrow=2,ncol=4)
a
a<-matrix(1:8, ncol=4,nrow=2)
a
a<-matrix(1:8,4,2)#Order of nonnamed arguments 
a
a<-matrix(1:8, ncol=4,nrow=2)

#ARITHMENTIC AND MANIPULATION
ab
ab+1

#manipulation and data extraction
#On vectors
ab
ab[4:6]
#On matrices
a[1,2:3]
a
t(a)

#On data frames
round(d[1:5,2],2)
round(d$temp[1:5],2)
attach(d)
round(temp[1:5],2)
detach(d)

#On list
l$f
l[2]

#On data frame
attach(d)
d[year==1995,]
d[year>=1995,]
d[year!=1995,]

round(sort(temp),2)
order(temp)
d[order(temp),]

temp>=4
1*(temp>=4)
detach(d)

#IMPORTING/EXPORTING
age0dist<-read.table('C:/MyDocuments/complex post-doc/Barents Sea Cod/Data Gjert/Cod age0 ab and length lc.csv',header=T,sep=';')
head(age0dist)

age0dist$lnage0<-log(age0dist$catch+1)
head(age0dist)

map<-read.table('C:/MyDocuments/complex post-doc/Barents Sea Cod/spatial analysis/coastline_highrange.txt')
names(map)<-c('x','y')
head(map)

#Map data are from NOAA Coastline extractor:
#http://rimmer.ngdc.noaa.gov/mgg/coast/getcoast.html

#EXPORTING
write.table(d,file='C:/DataExp')

#PLOTTING
#Scatterplot
x<-runif(50)
y<-3*x^2+0.3*rnorm(50)
plot(x,y,main='X Y Plot',xlab='X data',ylab='Y data',cex.lab=1.5,cex=1.5,cex.axis=1.5)

yt<-3*x^2#Underlying true function between x and y
res<-y-yt

par(mfrow=c(2,2))
plot(x,y,main='X Y Plot', xlab='X data', ylab='Ydata')
lines(sort(x),yt[order(x)])#This adds a line to the current plot
points(x[order(abs(res),decreasing=T)[1:10]],y[order(abs(res),decreasing=T)[1:10]],col='red',pch=16)
hist(x,main='Histogram of X',xlab='X data')
hist(res,main='Histogram of Residuals',xlab='Residuals')
qqnorm(res)
qqline(res)

#Pairs plot
newd<-data.frame(ind=x,dep=y,predictions=yt,deviations=res)
pairs(newd)

#coplot
x1<-runif(100)
x2<-runif(100)
y<-x1^2*exp(-2*x2)
coplot(y~x1|x2,pch=16)

#tapply and barplot
annual.average<-tapply(age0dist$lnage0,age0dist$year,mean,na.rm=T)
annual.ss<-tapply(age0dist$lnage0,age0dist$year,length)

par(mfrow=c(2,1))
barplot(annual.average,ylab='ln(cpue)',main='Average cpue')
box()
barplot(annual.ss,ylab='n of tows',main='Sample size')
box()

#Boxplot
boxplot(age0dist$lnage0~age0dist$year,outline=T,range=0)

#SPATIAL DATA
#bubble plots, one year
years<-2005
plot(map,type='l',ylim=range(age0dist$lat),xlim=range(age0dist$lon),
ylab='latitude (degree north)',xlab='longitude (degree east)',
main=paste(years,'Age0'))
symbols(age0dist$lon[age0dist$year==years],age0dist$lat[age0dist$year==years],
circles=age0dist$lnage0[age0dist$year==years],
inches=0.1,add=T,bg='purple')

#Multiple years and 'for' loop for one page only
par(mfrow=c(2,2))
years<-2001:2004
for(i in 1:length(years)){
plot(map,type='l',ylim=range(age0dist$lat),xlim=range(age0dist$lon),
ylab='latitude (degree north)',xlab='longitude (degree east)',
main=paste(years[i],'Age0'))
symbols(age0dist$lon[age0dist$year==years[i]],age0dist$lat[age0dist$year==years[i]],
circles=age0dist$lnage0[age0dist$year==years[i]],
inches=0.06,add=T,bg='purple')}

#Adjust plot area
windows()
par(mfrow=c(2,2),mai=c(0.7,0.6,0.3,0.1))
years<-2001:2004
for(i in 1:length(years)){
plot(map,type='l',ylim=range(age0dist$lat),xlim=range(age0dist$lon),
ylab='latitude (degree north)',xlab='longitude (degree east)',
main=paste(years[i],'Age0'))
symbols(age0dist$lon[age0dist$year==years[i]],age0dist$lat[age0dist$year==years[i]],
circles=age0dist$lnage0[age0dist$year==years[i]],
inches=0.06,add=T,bg='purple')}

#Multiple plotting areas: double for loop
years<-unique(age0dist$year)
tmp1<-1:ceiling(length(years)/4)
for(j in 1:length(tmp1)){
windows()
par(mfcol=c(2,2),mai=c(0.7,0.6,0.3,0.1))
for(i in (4*tmp1[j]-3):min(length(years),(4*tmp1[j]))){
plot(map,type='l',ylim=range(age0dist$lat),xlim=range(age0dist$lon),
ylab='latitude (degree north)',xlab='longitude (degree east)',
main=paste(years[i],'Age0'))
symbols(age0dist$lon[age0dist$year==years[i]],age0dist$lat[age0dist$year==years[i]],
circles=age0dist$lnage0[age0dist$year==years[i]],
inches=0.06,add=T,bg='purple')}}

#Standardize bubble size among years
md<-max(age0dist$lnage0,na.rm=T)
years<-unique(age0dist$year)
tmp1<-1:ceiling(length(years)/4)
for(j in 1:length(tmp1)){
windows()
par(mfcol=c(2,2),mai=c(0.7,0.6,0.3,0.1))
for(i in (4*tmp1[j]-3):min(length(years),(4*tmp1[j]))){
my<-max(age0dist$lnage0[age0dist$year==years[i]])
plot(map,type='l',ylim=range(age0dist$lat),xlim=range(age0dist$lon),
ylab='latitude (degree north)',xlab='longitude (degree east)',
main=paste(years[i],'Age0'))
symbols(age0dist$lon[age0dist$year==years[i]],age0dist$lat[age0dist$year==years[i]],
circles=age0dist$lnage0[age0dist$year==years[i]],
inches=0.06*my/md,add=T,bg='purple')}}

#DEPTH: import data and show contour
#Depth data are from:
#http://www.ngdc.noaa.gov/mgg/gdas/gd_designagrid.html
bathy.dat<-read.table('C:/MyDocuments/COAS/classes/RClassJune2007/BarentsDepth.txt',sep='')
names(bathy.dat)<-c('lon','lat','depth')
bathy.dat<-bathy.dat[bathy.dat$lon<=60,]#Crop the area
bathy.dat$depth[bathy.dat$depth>0]<-NA#Avoid points above water
head(bathy.dat)
bathy.mat<-matrix(bathy.dat$depth,nrow=length(unique(bathy.dat$lon)),ncol=length(unique(bathy.dat$lat)))[,order(unique(bathy.dat$lat))]

plot(map,type='l',ylim=range(age0dist$lat),xlim=range(age0dist$lon),ylab='',xlab='',main='Bottom depth')
contour(unique(bathy.dat$lon),sort(unique(bathy.dat$lat)),bathy.mat,levels=-seq(100,500,by=100),labcex=0.4,add=T,col='gray')
#image.plot(unique(bathy.dat$lon),sort(unique(bathy.dat$lat)),bathy.mat,col=terrain.colors(100))
symbols(age0dist$lon[age0dist$year==2005],age0dist$lat[age0dist$year==2005],
circles=age0dist$lnage0[age0dist$year==2005],
inches=0.05,add=T,bg='purple')

contour(unique(bathy.dat$lon),sort(unique(bathy.dat$lat)),bathy.mat,levels=-seq(100,500,by=100),labcex=0.4,add=T,col='gray')
lines(map,lwd=2)

#TEMPERATURE: import data and show image
years<-2004
temp.surf<-as.matrix(read.table(paste('C:/MyDocuments/complex post-doc/Barents Sea Cod/model temperature/surface/temp_surf_8_',years,'.asc',sep=''),header=F))
lat.t<-as.matrix(read.table('C:/MyDocuments/complex post-doc/Barents Sea Cod/model temperature/lat.asc',header=F))
lon.t<-as.matrix(read.table('C:/MyDocuments/complex post-doc/Barents Sea Cod/model temperature/long.asc',header=F))

library(fields)#needed to show image legend
image.plot(lon.t[1,],lat.t[,1],t(temp.surf),col=tim.colors(100),zlim=range(temp.surf,na.rm=T),ylab='Latitude N',xlab='Longitude E',main=paste('Surface Temperature',years))
symbols(age0dist$lon[age0dist$year==years],age0dist$lat[age0dist$year==years],
circles=age0dist$lnage0[age0dist$year==years],
inches=0.05,add=T,bg='purple')

contour(unique(bathy.dat$lon),sort(unique(bathy.dat$lat)),bathy.mat,levels=-seq(100,500,by=100),labcex=0.4,add=T,col='gray')
lines(map,lwd=2)
box()

#TEMPERATURE: all years
md<-max(age0dist$lnage0,na.rm=T)
years<-unique(age0dist$year)
tmp1<-1:ceiling(length(years)/4)
for(j in 1:length(tmp1)){
windows()
par(mfcol=c(2,2),mai=c(0.7,0.6,0.3,0.1))
for(i in (4*tmp1[j]-3):min(length(years),(4*tmp1[j]))){
temp.surf<-as.matrix(read.table(paste('C:/MyDocuments/complex post-doc/Barents Sea Cod/model temperature/surface/temp_surf_8_',years[i],'.asc',sep=''),header=F))
image.plot(lon.t[1,],lat.t[,1],t(temp.surf),col=tim.colors(100),zlim=range(-2,17),ylab='Latitude N',xlab='Longitude E',main=paste('Surface Temperature',years[i]))
symbols(age0dist$lon[age0dist$year==years[i]],age0dist$lat[age0dist$year==years[i]],
circles=age0dist$lnage0[age0dist$year==years[i]],
inches=0.04*my/md,add=T,bg='purple')
contour(unique(bathy.dat$lon),sort(unique(bathy.dat$lat)),bathy.mat,levels=-seq(100,500,by=100),labcex=0.4,add=T,col='gray')
lines(map,lwd=2)
box()}}


#INTERPOLATION of age0 density: one year
years<-2005
plot(map,type='l',ylim=range(age0dist$lat),xlim=range(age0dist$lon),
ylab='latitude (degree north)',xlab='longitude (degree east)',
main=paste(years,'Age0'))
symbols(age0dist$lon[age0dist$year==years],age0dist$lat[age0dist$year==years],
circles=age0dist$lnage0[age0dist$year==years],
inches=0.1,add=T,bg='purple')

#Build LOESS model: play with SPAN and DEGREE of the LOESS fit 
age0.loess<-loess(lnage0~lon*lat,span=0.08,degree=2,data=age0dist[age0dist$year==years,])
summary(lm(age0.loess$fitted~age0dist$lnage0[age0dist$year==years]))

#Make prediction grid
lond<-unique(bathy.dat$lon)
latd<-sort(unique(bathy.dat$lat))
predict.grid<-expand.grid(lond,latd)
names(predict.grid)<-c('lon','lat')
head(predict.grid)

#Make predictions using LOESS and plot results
age0.pred<-predict(age0.loess,newdata=predict.grid)
image.plot(lond,latd,age0.pred,col=tim.colors(100),zlim=range(age0.pred,na.rm=T),ylab='Latitude N',xlab='Longitude E',main=paste('Age0 density',years))
symbols(age0dist$lon[age0dist$year==years],age0dist$lat[age0dist$year==years],
circles=age0dist$lnage0[age0dist$year==years],
inches=0.07,add=T,bg='purple')
lines(map,lwd=2)
box()

#Remove predictions on land: use bottom depth matrix
bathy.dat<-read.table('C:/MyDocuments/COAS/classes/RClassJune2007/BarentsDepth.txt',sep='')
names(bathy.dat)<-c('lon','lat','depth')
bathy.dat<-bathy.dat[bathy.dat$lon<=60,]#Crop the area
bathy.dat$depth[bathy.dat$depth>0]<-NA#Avoid points above water
head(bathy.dat)
bathy.mat<-matrix(bathy.dat$depth,nrow=length(unique(bathy.dat$lon)),ncol=length(unique(bathy.dat$lat)))[,order(unique(bathy.dat$lat))]
land.mask<-(bathy.mat*0)+1
range(land.mask,na.rm=T)
age0.pred<-age0.pred*land.mask

#Outside of sampling area
library(sgeostat)#for convex hull
age0.hull<-chull(age0dist$lon[age0dist$year==years],age0dist$lat[age0dist$year==years])
age0.hull#gives indices of data hull
age0.poly<-list(x=age0dist$lon[age0dist$year==years][age0.hull],y=age0dist$lat[age0dist$year==years][age0.hull]) 
polygon(age0.poly)
within.hull<-in.chull(predict.grid$lon,predict.grid$lat,age0.poly$x,age0.poly$y)
age0.pred[!within.hull]<-NA;

#Remake Image
windows()
image.plot(lond,latd,age0.pred,col=tim.colors(100),zlim=range(age0dist$lnage0[age0dist$year==years]),ylab='Latitude N',xlab='Longitude E',main=paste('Age0 density',years))
symbols(age0dist$lon[age0dist$year==years],age0dist$lat[age0dist$year==years],
circles=age0dist$lnage0[age0dist$year==years],
inches=0.07,add=T,bg='purple')
lines(map,lwd=2)
box()

#CENTER OF GRAVITY: cpue-weighted longitude and latidue
years<-2005
plot(map,type='l',ylim=range(age0dist$lat),xlim=range(age0dist$lon),
ylab='latitude (degree north)',xlab='longitude (degree east)',
main=paste(years,'Age0'))
symbols(age0dist$lon[age0dist$year==years],age0dist$lat[age0dist$year==years],
circles=age0dist$lnage0[age0dist$year==years],
inches=0.1,add=T,bg='purple')

clon<-sum(age0dist$lon[age0dist$year==years]*age0dist$lnage0[age0dist$year==years])/
sum(age0dist$lnage0[age0dist$year==years])
clat<-sum(age0dist$lat[age0dist$year==years]*age0dist$lnage0[age0dist$year==years])/
sum(age0dist$lnage0[age0dist$year==years])

points(clon,clat,pch=16,col='red',cex=2)


#DISTANCE BETWEEN POINTS
#between center of gravity and all other points sampled
source('C:/MyDocuments/complex post-doc/R/distance.function.R')#load the function
dst<-distance.function(clat,clon,age0dist$lat[age0dist$year==years],age0dist$lon[age0dist$year==years])
#find points within 100km from center of gravity 
index<-dst<=100000
#plot them
points(age0dist$lon[age0dist$year==years][index],age0dist$lat[age0dist$year==years][index],pch='+')


#Import coastline and bottom depth data from NOAA web sites





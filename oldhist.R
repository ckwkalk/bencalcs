nsbenthosdata <-read.csv("~/no@guff.eu/Bucket_Trans/FEdp_20141008122542-csv.csv", header=TRUE)

library(plyr)
nbcount<-join(nsbenthosdata, count(nsbenthosdata, spec))
nbcount<-join(nsbenthosdata, count(nsbenthosdata, nsbenthosdata$spec))
nsbenthosdata$count <- ave(as.numeric(nsbenthosdata$spec), nsbenthosdata$spec, FUN = length)
nsbenthosdata$idx <- as.integer(interaction(nsbenthosdata$x, nsbenthosdata$y))
nb25 <- nsbenthosdata[nsbenthosdata$count > 25, ]
hist(nsbenthosdata$count)
hist(nsbenthosdata$count, breaks=20)
hist(nb25$count, breaks=20)
nb25 <- droplevels(nb25)
length(unique(nb25$idx)
)
length(unique(nsbenthosdata$idx))
nb25$Lon <- nb25$x
nb25$Lat <- nb25$y
bd <- nb25
crs.geo <- CRS("+proj=longlat +ellps=WGS84 +datum=WGS84")  # geographical, datum WGS84
ncrs.geo <- CRS("+proj=utm +zone=31 +ellps=WGS84 +datum=WGS84 +units=m +no_defs")
coordinates(bd) <- ~x+y
proj4string(bd) <- crs.geo  # define projection system of our data
bd <- spTransform(bd, ncrs.geo)
bd <-as.data.frame(bd, xy=TRUE)
library("rgdal")
#define names
bd <- nb25
crs.geo <- CRS("+proj=longlat +ellps=WGS84 +datum=WGS84")  # geographical, datum WGS84
ncrs.geo <- CRS("+proj=utm +zone=31 +ellps=WGS84 +datum=WGS84 +units=m +no_defs")
coordinates(bd) <- ~x+y
proj4string(bd) <- crs.geo  # define projection system of our data
bd <- spTransform(bd, ncrs.geo)
bd <-as.data.frame(bd, xy=TRUE)
out<- lm(bd$depth ~ poly(bd$Lat, 3))
depthbd <- predict(out,bd)
bd$depth2 <- impute(bd$depth, depthbd)
impute <- function (a, a.impute){ifelse (is.na(a), a.impute, a)}
bd$depth2 <- impute(bd$depth, depthbd)



bd <-normalise(bd, "abun", by="spid")
View(bd)
##Exclude entrys with NA for abundance - also xy?
bd<-bd[complete.cases(bd[,c(3)]),]
bd<-droplevels(bd)
bd <-normalise(bd, "abun", by="spid")
bd <- nb25
crs.geo <- CRS("+proj=longlat +ellps=WGS84 +datum=WGS84")  # geographical, datum WGS84
ncrs.geo <- CRS("+proj=utm +zone=31 +ellps=WGS84 +datum=WGS84 +units=m +no_defs")
coordinates(bd) <- ~x+y
proj4string(bd) <- crs.geo  # define projection system of our data
bd <- spTransform(bd, ncrs.geo)
bd <-as.data.frame(bd, xy=TRUE)
###fill depthnas with average values====
##doesnt work as some idx's are empty
#test <- ddply(bd, .(idx), function(bd) {bd$depth[is.na(bd$depth)] <- mean(bd$depth, na.rm=TRUE); return(bd)})
#install.packages("DMwR")
#library("DMwR", lib.loc="~/R/win-library/3.1")
##doesn't work because too many nas
#dn<- knnImputation(depthd)
##3rd order polynomial fit data
out<- lm(bd$depth ~ poly(bd$Lat, 3)) #out<- lm(bd$depth ~ poly(bd$Lat, 3), na.action="na.exclude")
##create output values - include dataframe to get same length of vector
depthbd <- predict(out,bd)
###function
##impute <- function (a, a.impute){ifelse (is.na(a), a.impute, a)}
##substitute na values for imputed polynomial fit values
bd$depth2 <- impute(bd$depth, depthbd)
#convert species name to id
bd$spid <-as.numeric(bd$spec)
bd$spid<-as.factor(bd$spid)
#bdf$spid <-as.factor(bdf$spid)
#bdfs2 <- bdf[complete.cases(bdf),]
#log transform abundance
bd$labun <- log1p(bd$abun)
#bd$llabun <- log(bd$labun+1)
##Exclude entrys with NA for abundance - also xy?
bd<-bd[complete.cases(bd[,c(2)]),]
bd<-droplevels(bd)
bd <-normalise(bd, "abun", by="spid")
range01 <- function(x){(100-1)/(max(x)-min(x))*((x)-min(x))+1}
rankabun<- ave(bd$abun*-1, bd$spec, FUN=rank)
bd$rnglogit<- ave(bd$qlogit, bd$spec, FUN=range01) ##use rnglogit for weighting centroid offsets
bd$rngrnkab<- ave(rankabun, bd$spec, FUN=range01)
bd$rankabun<- ave(bd$abun*-1, bd$spec, FUN=rank)
bd$rngrnkab<- ave(bd$rankabun, bd$spec, FUN=range01)

x_unique <- unique(bd$countsp)
x_ranks <- rank(x_unique)
bd$spfreqrnk<- x_ranks[match(bd$countsp,x_unique)]
bd$qrtsp <- with(bd, cut(countsp, breaks=quantile(countsp, probs=seq(0,1, by=0.2)),include.lowest=TRUE))

quintile_groups <- unique(bd$qrtsp)
bd$qrtsp <- as.numeric(bd$qrtsp)
View(bd)
plot(bd$spfreqrnk~bd$qrtsp)
plot(bd$rngrnkab~bd$qrtsp)
hist(bd$spfreqrnk)
bd$yj <-jitter(bd$y, amount=20000)
bd$xj <-jitter(bd$x, amount=20000)
library("ggplot2", lib.loc="~/R/x86_64-pc-linux-gnu-library/3.1")
ggplot()+ geom_point(data=bd, mapping=aes(x=xj, y=yj, color=rnglogit, size=labun)) +scale_size_continuous(range = c(3,20)) + scale_colour_gradientn(colours = rainbow(7))+ facet_wrap(~spec)
bd$qlogit <- qlogis(bd$norm.abun)
bd$qlogit<-bd$qlogit+min(bd$qlogit*-1)
bd$rnglogit<- ave(bd$qlogit, bd$spec, FUN=range01)
ggplot()+ geom_point(data=bd, mapping=aes(x=xj, y=yj, color=rnglogit, size=labun)) +scale_size_continuous(range = c(3,20)) + scale_colour_gradientn(colours = rainbow(7))+ facet_wrap(~spec)
ggplot()+ geom_point(data=bd, mapping=aes(x=xj, y=yj, color=rnglogit, size=labun)) +scale_size_continuous(range = c(1,10)) + scale_colour_gradientn(colours = rainbow(7))+ facet_wrap(~spec)
ggplot()+ geom_point(data=bd, mapping=aes(x=xj, y=yj, color=spfreqrnk, size=labun)) +scale_size_continuous(range = c(1,10)) + scale_colour_gradientn(colours = rainbow(7))+ facet_wrap(~spec)
install.packages("colorRamps")
library("colorRamps")
ggplot()+ geom_point(data=bd, mapping=aes(x=xj, y=yj, color=spfreqrnk, size=labun)) +scale_size_continuous(range = c(1,10)) + scale_colour_gradientn(colours = colours=matlab.like(10))+ facet_wrap(~spec)


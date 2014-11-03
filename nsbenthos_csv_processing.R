###TAKE SPECISE ANUNDANCE CSV AND IDENTIFY GROUPS OF SPECIEs ACCORDING TO THEIR DISTRIBUTION BEHAVIOUR##
#Load+CLEAN - CHECK FOR DUPLICATES BASED ON SAME COORDINATES,SPECIES AND ABUNDANCE, REMOVE nas
#1 CONVERT TO PROJECTED COORDINATES [csv>spatial>data.frame]
#2 NORMALISE THE ABUNDANCE VALUES
#3 CALCULATE K-NEAREST NEIGHBOUR 3RD, 6TH, 12TH SLOPE OF DECAY
#4 CALCULATE WEIGHTED MEANS
#5 CALCULATE DISTANCE AND HEADING BETWEEN WEIGHTED MEAN AND MEAN
##5a CALCULATE DISTANCE AND HEADING BASED ON UARTILES ?50%
##5b ADD NO OF RECORDS, ABUNDANCE -AVG-SD-SQUEWNESS PER SPECIES


#Load+CLEAN - CHECK FOR DUPLICATES BASED ON SAME COORDINATES,SPECIES AND ABUNDANCE, REMOVE nas
nsbenthosdata <-read.csv("E:/ADMIN_lib/no@guff.eu/Bucket_Trans/FEdp_2014100812254-csv2.csv", header=TRUE)


library(plyr)
##clean data====
#get counts per species
nbcount<-join(nsbenthosdata, count(nsbenthosdata, spec))
nbcount$idx <- as.integer(interaction(nbcount$x, nbcount$y))
nbcount<-ddply(nbcount, .(idx), mutate, count = length(unique(spec)))  
##dt$count <- ave(as.numeric(dt$school), dt$school, FUN = length)

#drop species with fewer than 20 representatives
nb25 <- nbcount[nbcount$countsp > 25, ] 
#nb25 <- nbcount[!(as.numeric(nbcount$spec) %in% which(table(nbcount$countsp)<25)),]
nb25 <- droplevels(nb25)
nb25$Lon <- nb25$x
nb25$Lat <- nb25$y

## counts of combinations of values - perhaps use for species by year?
##ddply(nsbenthosdata, .(spec?date?), mutate, count = length(unique(spec)))

#spec_count<- count(nsbenthosdata, vars = "spec")
#nsbenthosdata <- join(nsbenthosdata, spec_count)
#nsbenthosdata <- droplevels(subset(nsbenthosdata, freq > 20))

#1 CONVERT TO PROJECTED COORDINATES [csv>spatial>data.frame]
library("rgdal", lib.loc="~/R/win-library/3.1")

#define names
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


#log transform abundance
bd$labun <- log1p(bd$abun)
#bd$llabun <- log(bd$labun+1)


#normalise with function(see below)
library("vegan", lib.loc="~/R/win-library/3.1")
bd <-normalise(bd, "abun", by="spid")
##produces a bd$norm.abun column

#convert species name to id
bd$spid <-as.numeric(bd$spec)
bd$spid<-as.factor(bd$spid)
#bdf$spid <-as.factor(bdf$spid)
#bdfs2 <- bdf[complete.cases(bdf),]

##Exclude entrys with NA for abundance - also xy?
bd<-bd[complete.cases(bd[,c(2)]),] 
bd<-droplevels(bd)


##Logit transform
bd$qlogit <- qlogis(bd$norm.abun)
bd$qlogit<-bd$qlogit+min(bd$qlogit*-1)
#bd$norm.labun100 <- bdfs$norm.labun*100

## inverse rank abundance by group
bd$rankabun<- ave(bd$abun*-1, bd$spec, FUN=rank)
bd$rnglogit<- ave(bd$qlogit, bd$spec, FUN=range01) ##use rnglogit for weighting centroid offsets
bd$rngrnkab<- ave(bd$rankabun, bd$spec, FUN=range01)

## rank sum abundance tota
#bd$sumab<-ave(bd$abun, bd$spec, FUN=sum)
bd$rankabun<- ave(bd$abun*-1, bd$spec, FUN=rank)

x_unique <- unique(bd$countsp)
x_ranks <- rank(x_unique)
bd$spfreqrnk<- x_ranks[match(bd$countsp,x_unique)]

bd$qrtsp <- with(bd, cut(countsp, breaks=quantile(countsp, probs=seq(0,1, by=0.25)),include.lowest=TRUE))

#### select factors
#selected<-c("[26,51]", "(51,89]", "(89,159]","(159,534]") #want rows that contain A and B
#out<- data[data$Code %in% selected,]
#temp$quartile <- as.numeric(temp$quartile)
##Subset to plot images in https://docs.google.com/a/guff.eu/document/d/1BSHtIx-OlJPNzjC4lqs1olXXx0dyfqJoDXIvMNi2K0s/edit
ggplot()+ geom_point(data=bd100sp, mapping=aes(x=xj, y=yj, color=rnglogit, size=labun)) +scale_size_continuous(range = c(3,20)) + scale_colour_gradientn(colours = rainbow(7))+ facet_wrap(~spec)

ggplot(kdareas50) +
  geom_polygon(aes(x=long, y=lat, group = group, fill = id, colour = id),
               alpha = .4) +
  theme_bw() +
  coord_equal() + facet_wrap(~group)

o <- over(g, kareas50, returnList=TRUE)
ct <- sapply(o, length)
SGDF <- SpatialGridDataFrame(g, data=data.frame(ct=ct))
spplot(SGDF, "ct", col.regions=bpy.colors(20))


bd25 <- bd[which(bd$spfreqrnk<=25),]
bd50 <- bd[which(bd$spfreqrnk>=25 & bd$spfreqrnk<=50),]
bd75 <- bd[which(bd$spfreqrnk>=50 & bd$spfreqrnk<=75),]
bd100<- bd[which(bd$spfreqrnk>=75),]

ggplot(bd100, aes(rnglogit, fill = spec)) + geom_density(aes(alpha =bd100$countsp))+scale_alpha(range=c(0.1,0.8))+ggtitle("bd100")
ggplot(bd75, aes(rnglogit, fill = spec)) + geom_density(aes(alpha =bd75$countsp))+scale_alpha(range=c(0.1,0.8))+ggtitle("bd75")
ggplot(bd50, aes(rnglogit, fill = spec)) + geom_density(aes(alpha =bd50$countsp))+scale_alpha(range=c(0.1,0.8))+ggtitle("bd50")
ggplot(bd25, aes(rnglogit, fill = spec)) + geom_density(aes(alpha =bd25$countsp))+scale_alpha(range=c(0.1,0.8))+ggtitle("bd25")
ggplot(bd50, aes(k12, fill = spec)) + geom_density(aes(alpha =bd50$countsp))+scale_alpha(range=c(0.1,0.8))+ggtitle("bd50_k12")
ggplot(bd50, aes((k12^.5), fill = spec)) + geom_density(aes(alpha =bd50$countsp))+scale_alpha(range=c(0.1,0.8))+ggtitle("bd50_k12")
ggplot(bd, aes((k12^.3), fill = spec)) + geom_density(aes(alpha =bd$countsp))+scale_alpha(range=c(0.1,0.8))+ggtitle("bd_k12")
ggplot(bd, aes(log(k12), fill = spec)) + geom_density(aes(alpha =bd$countsp))+scale_alpha(range=c(0.1,0.8))+ggtitle("bd_k12")


#####Test on subset but should be run on whole

###i/1/ autofit variogram-----
###create spatial subset for variogram and kernel density only cordinates and spec needed
bdsub <- bd100[which(bd100$spec=="Turbellaria"),]
bdsub <- bdsub[,c("xj","yj","labun")]
bd50sp <-droplevels(bd50)
coordinates(bd50sp) <- ~xj+yj
gpbd50 <- split(bd50sp, bd50sp$spec)

##variogram = autofitVariogram(rnglogit~1,Turb, miscFitOptions = list(min.np.bin = 150)), need to have a large bin size to remove high variability at close distances

## for species with high local patchiness the range is always small at low bin sizes abundand species larger bins are better for less freq smaller are needed

krigeall.list_rnglogit50=lapply(gpbd50,function(i) {
  lapply(names(i@data)[23],function(j) {
    as.data.frame(autofitVariogram(formula(paste("(",j,")~1")),i, , miscFitOptions = list(min.np.bin = 10))$var_model$range[2]) } ) } )

var50r <- melt(krigeall.list_rnglogit50)

krigeall.list_rnglogit50=lapply(gpbd50,function(i) {
  lapply(names(i@data)[23],function(j) {
    as.data.frame(autofitVariogram(formula(paste("(",j,")~1")),i, , miscFitOptions = list(min.np.bin = 10))$var_model$psill[2]) } ) } )

var50s <- melt(krigeall.list_rnglogit50)

######i/2/ KERNEL DENSITY

### nice visualisation for one species at a time (shows spatial variation vvertically) or with 3 d
## representing time so changes over the years can be detected
##library("feature", lib.loc="~/R/win-library/3.1")
##featureSignifGUI(bdsub)bdsub.Anth <- bd100[which(bd100$spec=="Anthozoa"),]
#bdsub.Anth <- bdsub.Anth[,c("xj","yj","labun")]
#bdsub.fsa <- featureSignif(bdsub.Anth, bw=c(50000, 50000, 0.8), signifLevel=0.3)
#plot(bdsub.fsa, addKDE=TRUE, addSignifCurvRegion=TRUE)

###Create duplicate but spatially offset values based on abundance value====
names(bd100)
## better to use range rank abundance/30 instead of labun?
bd100l<-data.frame(bd100[rep(seq_len(dim(bd100)[1]), bd100$labun), c(1,25:26), drop = FALSE], row.names=NULL)
## jitter only duplicates (better kernel density analysis)====
jitter.duplicate <- function(x, only_positive = F)
{
  if(only_positive) {
    ss <- x > 0
  } else {
    ss <- T
  }
  ss_dup <- duplicated(x[ss])
  # ss <- ss & ss_dup
  temp_length <- length(x[ss][ss_dup])
  x[ss][ss_dup] <- x[ss][ss_dup] + seq(from = 1000, to = 2000, length.out = temp_length)
  x
}
#====
bd100l$yj<-jitter.duplicate(bd100l$yj)
bd100l$xj<-jitter.duplicate(bd100l$xj)
bd100l<-droplevels(bd100l)
tagList <- unique(bd100l$spec)
df<-bd100l

coordinates(df) <- ~xj+yj
proj4string(df) <- ncrs.geo
kud <- kernelUD(df[ ,1], grid = 400, same4all=TRUE)
hr <- kernel.area(kud, percent = 75)
plot(hr)
kareas <- getverticeshr(kud, 75)


kud1 <- list()
kud_spdf <- list()
vd <- list()

## running a loop across animals, estimating 95% 
for(i in tagList){
  thing <- SpatialPoints(subset(df[df$ID == i,], select = c(xj, yj)),   
                         proj4string = CRS("+proj=utm +zone=31 +ellps=WGS84"))
  thing <- SpatialPointsDataFrame(coords = thing, data = subset(df[df$ID == i,], select = c(Season)))
  thing$Season <- as.character(thing$Season)
  kud1[[i]] <- kernelUD(thing, grid = 120, extent = 0.2, same4all=TRUE)
  kud_spdf[[i]] <- estUDm2spixdf(kud1[[i]]) 
  fullgrid(kud_spdf[[i]]) <- TRUE
  vd[[i]] <- kud_spdf[[i]]@data 
}


kud <- kernelUD(df[ ,1], grid = 20000, same4all=TRUE)
hr <- kernel.area(kud, percent = 75)

kareas <- getverticeshr(kud, 75)

##i//Get the number of pieces from kernel density polygons
#library(ggplot2)----
kdareas <- fortify(kareas)
ggplot(kareas) + geom_polygon(aes(x=long, y=lat, group = group, fill = id, colour = id), alpha = .4) + theme_bw() + coord_equal()

##add number of pieces to main summry table once created match based on spec
kdareas$max<- ave(as.numeric(kdareas$piece), kdareas$id, FUN=max)


##i//calculate perimiter
data = readShapePoly("shapefile.shp")
>>> data.polydata = as.PolyData(data@data)
>>> calcLength(data.polydata)

##i//calculate area
# total area
sapply(slot(mp, "polygons"), slot, "area")

# get list of individual polys
p <- lapply(mp@polygons , slot , "Polygons")

# areas of individual polygons
lapply(p[[1]], function(x) slot(x, "area"))
install.packages("rgeos")

centroid locations
library("rgeos", lib.loc="~/R/win-library/3.1")
t<-gCentroid(kareas30,byid=TRUE)
plot(t)
outt <- lapply(kareas30@polygons , slot , "area"[[2]] )


###function----
###range01 <- function(x){(100-1)/(max(x)-min(x))*((x)-min(x))+1}
###invrange = function(x){100*(1/((100-1)/(max(x)-min(x))*((x)-min(x))+1))}


####function to apply max value to INF/NA----
##new <- ddply(bd, .(spec), function(bd) {bdfs$qlogit[bdfs$qlogit == Inf] <- max(bdfs$qlogit[is.finite(bdfs$qlogit)], na.rm=TRUE); return(bd)})

#calculate kth nearest neighbour
library("spatstat", lib.loc="~/R/win-library/3.1")
bd$xj <-jitter(bd$x, amount=30000)
bd$yj <-jitter(bd$y, amount=30000)

#calculate kth nearest neighbour
bd$k3 <- nndist(bd[,c("xj","yj")], k=3, by =bd$spid)
#bd$k6 <- nndist(bd[,c("xj","yj")], k=6, by =bd$spid)
#hist(bd$k6)
bd$k12 <- nndist(bd[,c("xj","yj")], k=12, by =bd$spid)
bd$k25 <- nndist(bd[,c("xj","yj")], k=25, by =bd$spid)


bd$k3 <- apply(bd$k3, 1, FUN = function(x) {min(x[x > 0])})
bd$k6 <- apply(bd$k6, 1, FUN = function(x) {min(x[x > 0])})
bd$k12 <- apply(bd$k12, 1, FUN = function(x) {min(x[x > 0])})
bd$k25 <- apply(bd$k25, 1, FUN = function(x) {min(x[x > 0])})

bd100 <-bd[(as.numeric(factor(bd$countsp)))>16,]

###Nice plot showing all abundance disributions faceted (could also colour by quartiles?)
install.packages("colorRamps")
library("colorRamps", lib.loc="~/R/win-library/3.1")
ggplot()+ geom_point(data=bd100sp, mapping=aes(x=xj, y=yj, color=rnglogit, size=labun, alpha=0.1))  +coord_equal() + scale_size_continuous(range = c(3,15)) + scale_colour_gradientn(colours=matlab.like(10))+ facet_wrap(~spec)
ggplot() + geom_point(data=bd100sp, mapping=aes(x=xj, y=yj, color=-rngrnkab, size=labun, alpha=0.1)) +coord_equal() + scale_size_continuous(range = c(3,15)) + scale_colour_gradientn(colours=matlab.like(10))+ facet_wrap(~spec)

##cluster based on values calculated in previous steps - use Tanagra to remove covariables R ClustOfVar package? 
spclust <- Mclust(bd[,c(9,11,16,20,23,27,29)], initialization=list(subset=sample(1:nrow(bd), size=4000)))
plot(spclust)
bd$clust <-spclust$classification


######MAJOR STEP GOING FROM ALL RECORDS TO SUMMARISING BY SPECIES#####################

library("moments", lib.loc="~/R/win-library/3.1")
library("dplyr", lib.loc="~/R/win-library/3.1")

bdfsa <- ddply(bdfs, .(spid,spec), summarise,
               mnx = mean(x),
               mny = mean(y),
               wax = weighted.mean(x, norm.labun),
               way = weighted.mean(y, norm.labun),
               mnlabun=mean(labun),
               sdlabun=sd(labun),
               sklabun=skewness(labun),
               count=length(unique(xj)),
               mnk3=mean(k3),
               mnk6=mean(k6),
               mnk12=mean(k12),
               sdk3=sd(k3, na.rm=TRUE),
               sdk6=sd(k6, na.rm=TRUE),
               sdk12=sd(k12, na.rm=TRUE),
               skk3=skewness(k3, na.rm=TRUE),
               skk6=skewness(k6, na.rm=TRUE),
               skk12=skewness(k12, na.rm=TRUE)
)

bdfsa
library("geosphere", lib.loc="~/R/win-library/3.1")
bdfsa$head <- (270-((atan2((bdfsa[,c(4)]-bdfsa[,c(6)]),(bdfsa[,c(3)]-bdfsa[,c(5)])))*(180/pi)))%%360 
bdfsa$dist <- ((bdfsa$mnx-bdfsa$wax)^2+(bdfsa$mny-bdfsa$way)^2)^.5

library("mclust", lib.loc="~/R/win-library/3.1")
cbdfsan <- na.omit(bdfsa)
cbdfsa <- Mclust(cbdfsan[,c(5:11, 13:17, 20)])
save(cbdfsa)

library("fields", lib.loc="~/R/win-library/3.1")
library("ggplot2", lib.loc="~/R/win-library/3.1")
library("scales", lib.loc="~/R/win-library/3.1")

cbdfsan$arrow <- rescale(cbdfsan$dist, to=c(1,3))
ggplot(data=cbdfsan[which(cbdfsan$dist > 0.05),], aes(x=mnx, y=mny, size=arrow)) + geom_segment(aes(xend=wax, yend=way, colour=head + arrow), arrow = arrow(length = unit(0.5,"cm"))) + scale_colour_gradientn(colours = rainbow(7))

ggplot(data=cbdfsan[which(cbdfsan$count > 25),], aes(x=mnx, y=mny, size=arrow)) + geom_segment(aes(xend=wax, yend=way, colour=head + arrow), arrow = arrow(length = unit(0.5,"cm"))) + scale_colour_gradientn(colours = rainbow(7))



##scale to count weighted xy
geom_point(data=cbdfsan, mapping=aes(x=wax, y=way, color=clust, size=count)) +
  scale_size_continuous(range = c(3,20))+ scale_colour_gradientn(colours = rainbow(7))



##scale to avg abundance weighted xy
geom_point(data=cbdfsan, mapping=aes(x=wax, y=way, color=clust, size=mnlabun)) +
  scale_size_continuous(range = c(3,20))+ scale_colour_gradientn(colours = rainbow(7))

library("xlsx", lib.loc="~/R/win-library/3.1")
write.xlsx(cbdfsan, file="nsbenthoskneibour.xls")


ggplot(data=cnsbavg, aes(x=mnx, y=mny)) + geom_segment(aes(xend=wax, yend=way), arrow = arrow(length = unit(0.3,"cm")))



###extract all species >24 to find which variables contribute # analysis in TANAGRA
selnsb<- subset(bdfs, ave(rep(1,nrow(bdfs)), bdfs$spid, FUN=sum) > 24)
write.csv(selnsb, file="selnsb.csv")


##not run #selnsb<- bdfs[bdfs$spid %in% names(table(bdfs$spid))[table(bdfs$spid) > 24],]

c <- cbind(cnsb$f, cnsb$Tabu)

groups <- split(c, c$f)   
gLengths <-  lapply(groups, FUN=length)   
grange    <- lapply(groups, FUN=range02)            
out.L <- lapply(seq_along(groups), function(i) rep_len(1/gSums[[i]], length.out=gLengths[[i]]) * cumsum(groups[[i]])) 
return(out.L)} 

######ADDITIONAL SUMMARY VALUES TO ADD
##//
##//
##//VARIOGRAMS
##create new dataset for varioram analysis just xy and (?qlogit range) abundance
##???outt <- lapply(kareas30@polygons , slot , "area"[[2]] )


bd50sp <-droplevels(bd50)
coordinates(bd50sp) <- ~xj+yj
gpbd50 <- split(bd50sp, bd50sp$spec)

##AUTOFIT VARIOGRAMS DO TWO onne at min values = 20 and one at 180

##To check individual species settings min.np.bin affects the sensitivity should be ok with low values 10-40 but if there is significant short range variability then larger bins will be needed to get a good value for range 
#variogram = autofitVariogram(rnglogit~1, Turb, model=c("Sph"), fix.values=c(180,NA,NA), miscFitOptions = list(min.np.bin = 140))
#plot(variogram)
#'i@data)[23] refers to the abundance value and will need to be edited
##psil
krigeall.list_rnglogit50=lapply(gpbd50,function(i) {
  lapply(names(i@data)[23],function(j) {
    as.data.frame(autofitVariogram(formula(paste("(",j,")~1")),i, , miscFitOptions = list(min.np.bin = 10))$var_model$psill[2]) } ) } )
##range
krigeall.list_rnglogit50=lapply(gpbd50,function(i) {
  lapply(names(i@data)[23],function(j) {
    as.data.frame(autofitVariogram(formula(paste("(",j,")~1")),i, , miscFitOptions = list(min.np.bin = 10))$var_model$range[2]) } ) } )

##extract the value RANGE and PSILL 
var50s <- melt(krigeall.list_rnglogit50)

##var50<-ddply(var50, .(L1), subset, value == max(value))
View(var50)
plot(outt~var50$value )










##GENERAL NOTE BITS-----
##make wide
##weather3 %>% spread(element, value)

# range function =====
#(new_max - new_min) / (old_max - old_min) * (v - old_min) + new_min

range01 <- function(x){(100-1)/(max(x)-min(x))*((x)-min(x))+1}

invrange = function(x){100*(1/((100-1)/(max(x)-min(x))*((x)-min(x))+1))}


# normalise function====

normalise <- function(dataframe, columns, by, na.rm=TRUE){
  # Risch et al 2007 formula 3
  #
  # y_norm,i = y_i / sqrt(sum(y_i^2))
  #if(length(by)>2) stop("maximum of 2 by variables supported")
  if(length(by)==1){
    by.1 <- list()
    for(i in levels(dataframe[,by[1]])){
      by.1[[i]] <- dataframe[dataframe[,by]==i,]
      by.1[[i]] <- cbind(by.1[[i]], vegan:::decostand(by.1[[i]][,columns], 
                                                      method="normalize", MARGIN=2))
      names(by.1[[i]])[dim(by.1[[i]])[2]] <- paste("norm.",names(by.1[[i]][columns]), sep="")
    }
    for(i in levels(dataframe[,by[1]])){
      if(i == 1) {OUT <- by.1[[1]]} else {
        OUT <- rbind(OUT, by.1[[i]])}}
  }
  
}

###loganderson=====
loganderson <- function(dataframe, columns, by, na.rm=TRUE){
  # Risch et al 2007 formula 3
  #
  # y_norm,i = y_i / sqrt(sum(y_i^2))
  #if(length(by)>2) stop("maximum of 2 by variables supported")
  if(length(by)==1){
    by.1 <- list()
    for(i in levels(dataframe[,by[1]])){
      by.1[[i]] <- dataframe[dataframe[,by]==i,]
      by.1[[i]] <- cbind(by.1[[i]], vegan:::decostand(by.1[[i]][,columns], 
                                                      method="log", MARGIN=2))
      names(by.1[[i]])[dim(by.1[[i]])[2]] <- paste("logand.",names(by.1[[i]][columns]), sep="")
    }
    for(i in levels(dataframe[,by[1]])){
      if(i == 1) {OUT <- by.1[[1]]} else {
        OUT <- rbind(OUT, by.1[[i]])}}
  }
  if(length(by)==2){
    by.1 <- list()
    OUT.1 <- list()
    for(i in levels(dataframe[,by[1]])){
      by.1[[i]] <- dataframe[dataframe[,by[1]]==i,]
      by.2 <- list()
      for(j in levels(by.1[[i]][,by[2]])){
        by.2[[j]] <- by.1[[i]][by.1[[i]][,by[2]]==j,]
        by.2[[j]] <- cbind(by.2[[j]], vegan:::decostand(by.2[[j]][,columns], 
                                                        method="log", MARGIN=2, na.rm=na.rm))
        names(by.2[[j]])[dim(by.2[[j]])[2]] <- paste("logand.",
                                                     names(by.2[[j]][columns]), 
                                                     sep="")
      }
      for(p in 1:length(by.2)){
        if(p==1){OUT.1[[i]] <- by.2[[1]]} else {
          OUT.1[[i]] <- rbind(OUT.1[[i]], by.2[[p]])}
      }
    }
    for(i in 1:length(levels(dataframe[,by[1]])))
      ifelse(i==1,OUT <- OUT.1[[1]],OUT <- rbind(OUT, OUT.1[[i]]))
  } 
  OUT
}


normalise <- function(dataframe, columns, by, na.rm=TRUE){
  # Risch et al 2007 formula 3
  #
  # y_norm,i = y_i / sqrt(sum(y_i^2))
  #if(length(by)>2) stop("maximum of 2 by variables supported")
  if(length(by)==1){
    by.1 <- list()
    for(i in levels(dataframe[,by[1]])){
      by.1[[i]] <- dataframe[dataframe[,by]==i,]
      by.1[[i]] <- cbind(by.1[[i]], vegan:::decostand(by.1[[i]][,columns],
                                                      method="normalize", MARGIN=2))
      names(by.1[[i]])[dim(by.1[[i]])[2]] <- paste("norm.",names(by.1[[i]][columns]), sep="")
    }
    for(i in levels(dataframe[,by[1]])){
      if(i == 1) {OUT <- by.1[[1]]} else {
        OUT <- rbind(OUT, by.1[[i]])}}
  }
}
bd <-normalise(bd, "abun", "spec")
normalise <- function(dataframe, columns, by, na.rm=TRUE){
  # Risch et al 2007 formula 3
  #
  # y_norm,i = y_i / sqrt(sum(y_i^2))
  #if(length(by)>2) stop("maximum of 2 by variables supported")
  if(length(by)==1){
    by.1 <- list()
    for(i in levels(dataframe[,by[1]])){
      by.1[[i]] <- dataframe[dataframe[,by]==i,]
      by.1[[i]] <- cbind(by.1[[i]], vegan:::decostand(by.1[[i]][,columns],
                                                      method="normalize", MARGIN=2))
      names(by.1[[i]])[dim(by.1[[i]])[2]] <- paste("norm.",names(by.1[[i]][columns]), sep="")
    }
    for(i in levels(dataframe[,by[1]])){
      if(i == 1) {OUT <- by.1[[1]]} else {
        OUT <- rbind(OUT, by.1[[i]])}}
  }
  if(length(by)==2){
    by.1 <- list()
    OUT.1 <- list()
    for(i in levels(dataframe[,by[1]])){
      by.1[[i]] <- dataframe[dataframe[,by[1]]==i,]
      by.2 <- list()
      for(j in levels(by.1[[i]][,by[2]])){
        by.2[[j]] <- by.1[[i]][by.1[[i]][,by[2]]==j,]
        by.2[[j]] <- cbind(by.2[[j]], vegan:::decostand(by.2[[j]][,columns],
                                                        method="normalize", MARGIN=2, na.rm=na.rm))
        names(by.2[[j]])[dim(by.2[[j]])[2]] <- paste("norm.",
                                                     names(by.2[[j]][columns]),
                                                     sep="")
      }
      for(p in 1:length(by.2)){
        if(p==1){OUT.1[[i]] <- by.2[[1]]} else {
          OUT.1[[i]] <- rbind(OUT.1[[i]], by.2[[p]])}
      }
    }
    for(i in 1:length(levels(dataframe[,by[1]])))
      ifelse(i==1,OUT <- OUT.1[[1]],OUT <- rbind(OUT, OUT.1[[i]]))
  }
  OUT
}





## geometric mean function========
gm_mean = function(x, na.rm=TRUE){
  exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
}

###Multi-histogramplot=====
multi.hist <- function(x) {nvar <- dim(x)[2]  #number of variables
                           nsize=trunc(sqrt(nvar))+1   #size of graphic
                           old.par <- par(no.readonly = TRUE) # all par settings which can be changed
                           par(mfrow=c(nsize,nsize))       #set new graphic parameters
                           for (i in 1:nvar) {
                             name=names(x)[i]                #get the names for the variables
                             hist(x[,i],main=name,xlab=name) }  #draw the histograms for each variable
                           on.exit(par(old.par))   #set the graphic parameters back to the original
}

##multihistogram ggplot====
chats$from_country <- factor(chats$from_country, 
                             levels = unique(c(chats$from_country, 
                                               chats$to_country)))
bdf$spec <- factor(bdf$spec, 
                   levels = levels(chats$spec))



sp_bdf <- melt(bdfs[,c("spec","abun")])

ggplot(subset(sp_bdfs, ave(rep(1,nrow(sp_bdfs)), sp_bdfs$value, FUN=sum, na.rm=TRUE) > 160), aes(x = value)) + 
  facet_wrap(~spec,scales = "free_x") + 
  geom_histogram()

##replace inf----
bdfs$qlogit[which(apply(bdfs$qlogit, 1, function(x) x == max(x,na.rm=TRUE)))] <- -1

df %>% group_by(Subject) %>% filter(pt == max(pt))

new <- ddply(bdfs, .(spec), function(bdfs) {bdfs$qlogit2[bdfs$qlogit2 == Inf] <- max(bdfs$qlogit2, na.rm=TRUE); return(bdfs)})
mx <- ddply(mx, .(groups), function(df) {df$value[is.na(df$value)] <- mean(df$value, na.rm=TRUE); return(df)})

###Jitter duplicate by single coordiante- needed to do kernel density - package geoR can do both xy=====

jitter.duplicate <- function(x, only_positive = F)
{
  if(only_positive) {
    ss <- x > 0
  } else {
    ss <- T
  }  
  ss_dup <- duplicated(x[ss])
  # ss <- ss & ss_dup
  temp_length <- length(x[ss][ss_dup])	
  x[ss][ss_dup] <- x[ss][ss_dup] + seq(from = 1000, to = 2000, length.out = temp_length)
  x
}

##Function Points to polygons for rebuilding Kernel Density polygons from fortified objects----
points2polygons2 <- function(df,data) {
  get.grpPoly <- function(group,ID,df) {
    Polygon(coordinates(df[df$idn==ID & df$group==group,]))
  }
  get.spPoly  <- function(ID,df) {
    Polygons(lapply(unique(df[df$idn==ID,]$group),get.grpPoly,ID,df),ID)
  }
  spPolygons  <- SpatialPolygons(lapply(unique(df$idn),get.spPoly,df))
  SpatialPolygonsDataFrame(spPolygons,match.ID=T,data=data)
}

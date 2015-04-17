##Take Long list of Species presence by UTM xya (created in QGIS)


SpeciesUTM <- read.csv("~/Documents/ENV_VARS/Null_Model/SpeciesUTM.csv")
##install.packages("spatstat")
uSpec <- unique(SpeciesUTM)
uSpec<-droplevels(uSpec)
library("plyr", lib.loc="/home/monst/R/x86_64-pc-linux-gnu-library/3.1")
uSpec<-join(uSpec,(count(uSpec,"Species")))
library("spatstat")
##average third nearest neighbour
uSpec$kk8<-nndist(uSpec[,c("Long","Lat")], k=2, by=uSpec$Species)
##need to tidy format
uSpec$kk8<-apply(uSpec$kk8,1,FUN=function(x) {min(x[x>0])})
uSpec$kk3<-nndist(uSpec[,c("Long","Lat")], k=3, by=uSpec$Species)
uSpec$kk3<-apply(uSpec$kk3,1,FUN=function(x) {min(x[x>0])})
uSpec$kk2<-nndist(uSpec[,c("Long","Lat")], k=2, by=uSpec$Species)
uSpec$kk2<-apply(uSpec$kk2,1,FUN=function(x) {min(x[x>0])})
## normalised standard deviation
uSpec$sd3z <- ave(uSpec$kk3, uSpec$Species, FUN=scale)
uSpec$sd8z <- ave(uSpec$kk8, uSpec$Species, FUN=scale)


##install.packages("moments")
library("moments")
##create wide summaruy by Species
library("dplyr")
guSpec<-uSpec %>% group_by(Species) %>%
  + summarise(freq=mean(freq), mn2=mean(kk2), mn3=mean(kk3), mn8=mean(kk8), sd2=sd(kk2), sd3=sd(kk3), sd8=sd(kk8), sk2=skewness(kk2), sk3=skewness(kk3), sk8=skewness(kk8), IQR3=IQR(sd3z, na.rm=TRUE, type=8), IQR8=IQR(sd8z, na.rm=TRUE, type=8))


## nice graph showing relationship between frequency groups and nearest neighbour distance
guSpec$lfreq<-log1p(guSpec$freq)
plot((guSpec$sd2)~(sqrt(guSpec$lfreq)))

ggplot(uSpec, aes(x=log(kk8),y=log(kk2),colour=(freq)^.2))+ geom_point(alpha = 0.5, size=3) +
scale_colour_gradientn(colours=rainbow(4))

guSps<-guSpec[which(guSpec$freq>=3),]

#####CLUSTERS TO JOIN TO SPECIES SUMMARY DATA#####
##Group from Gnumeric nnDspec=guSpec based on quartile manual regression fit of nearneighbour$sd3 int 4 groups   
nnDspec <- read.csv("~/Documents/R/nnDspec.csv")
##to plot the maxent rasters according to
nnDspec$filenn<- substr(nnDspec$filen,1,nchar(as.character(nnDspec$filen))-4)

G1list<-nnDspec[nnDspec$Group == 1, c(18)]
G2list<-nnDspec[nnDspec$Group == 2, c(18)]
G3list<-nnDspec[nnDspec$Group == 3, c(18)]
G4list<-nnDspec[nnDspec$Group == 4, c(18)]

G1sub40<-sample(G1list,40)
plot(s, G1sub40)

##Import summaryCSV of MAxent Model of Simplified real variables 
RSimpMxt <- read.csv("~/Documents/ENV_VARS/Maxent_Out_Simp/RSimpMxt.csv")
nnDspec$filenn<-nnDspec$Species
##nnDspec$Species<-nnDspec$filenn
mergeRS<-merge(x = RSimpMxt, y = nnDspec, by = "Species")
plot(mergeRS$Entropy~mergeRS$Group)
boxplot(mergeRS$Entropy~mergeRS$Group)

##Create matrix to calculate dissimilarity
SpeciesUTM$pres <- 1
library("reshape2")
SW <- dcast(SpeciesUTM, Long + Lat ~ Species, value.var="pres")
View(SW)
##Extract coordinates
SWcoo<-SW[,1:2]
SW[,1]<-NULL
SW[,1]<-NULL
SW[is.na(SW)] <- 0
##install.packages("pvclust")
library("pvclust")
resultW <- pvclust(SW, method.dist="binary", method.hclust="ward", nboot=2500)
plot(resultW)
pvrect(resultW, alpha=0.95)
clusterW <- pvpick(resultW)

resultW2 <- pvclust(SW, method.dist="binary", method.hclust="ward.D2", nboot=8000)
plot(resultW2)
pvrect(resultW2, alpha=0.95)
clusterW2 <- pvpick(resultW2)
seplot(resultW2, identify=TRUE)

s2names<-names(s2)
s2names<-gsub("_", " ", s2names)
names(s2)<-s2names

clusternW <- rapply(clusterW, function(x) gsub(" ", ".", x), how = "replace")
clusternW <- rapply(clusternW, function(x) gsub("\\(", ".", x), how = "replace")
clusternW <- rapply(clusternW, function(x) gsub("\\)", ".", x), how = "replace")
##plot(s2, clustern$clusters[[15]])
library("reshape")
clustWm<-melt(clusternW, id=c(1,))
names(clustWm)<-c("filenn","Wclust","delete")
##W1clust <- read.csv("~/Documents/ENV_VARS/Null_Model/Maxent/test/W1clust.csv")
##wclust <- rapply(W1clust, as.character, classes="factor", how="replace")
plot(s2,clusternW[[2]])



##nnDspec <- read.csv("~/Documents/R/nnDspec.csv")
#nnDspec$Species<-nnDspec$filenn
#nnDspec$filen<-nnDspec$Species

#guSpec$filenn<-gsub("\\(", ".", guSpec$filenn)
#guSpec$filenn<-gsub("\\)", ".", guSpec$filenn)
#guSpec$filenn<-gsub("_", ".", guSpec$filenn)
#guSpec$filenn<-gsub(" ", ".", guSpec$filenn)
#clustWmm<-merge(x = clustWm, y = guSpec[,c("filenn","freq")], by="filenn")

nnDspec$filenn<-gsub("\\(", ".", nnDspec$filenn)
nnDspec$filenn<-gsub("\\)", ".", nnDspec$filenn)
nnDspec$filenn<-gsub("_", ".", nnDspec$filenn)
nnDspec$filenn<-gsub(" ", ".", nnDspec$filenn)
##nnDspec$filenn<-nnDspec$Species
nnDspecW<-merge(x = clustWm, y = nnDspec, by="filenn")


guSpsm<-merge(x = nnDspecW, y = guSps, by = "filenn")
plot(guSpsm$nnClust~guSpsm$freq.x)
boxplot(guSpsm$nnClust~guSpsm$Group)

guSpsm$sumF <- ave(guSpsm$freq.x, guSpsm$Wclust, FUN=sum)

ggplot(uSpec, aes(x=log(freq),y=(sd8z),colour=as.numeric(Species)))+ geom_point(alpha = 0.5, size=3) +
scale_colour_gradientn(colours=rainbow(10))

Wclustfreq<-guSpsm %>%group_by(Wclust) %>%
  summarise(mean=mean(freq.x))

c3 <- t(Wclustfreq$mean)
boxplot(guSpsm$nnClust~guSpsm$Wclust, width=c3[1:17], col="gray")


guClust3 <- Mclust(guSps2[c(6:14)])
plot(guClust3)




###Quantile regression -see end of file for example and quantreg library 
###nice idea but didn't work had to do it manually in Gnumeric ~/monst/Documents/R/nnDist_Abundance_Calc.gnumeric


##cell statistics of raster layers
s2St<-cellStats(s , 'sum')
sname<-names(s)
sSum<-data.frame(sname,s2St)
colnames(sSum) <- c("Species", "sSum")
G4sub6<-sample(G4list,6)
sub6G4s<-subset(s,G4sub6)
##install.packages("automap")
library("automap", lib.loc="/home/monst/R/x86_64-pc-linux-gnu-library/3.1")

Vsub6G4s<-autofitVariogram(sub6G4s)
##Fails - perhaps package error
##Fix projections
proj4string(sub6G4s) <- CRS("+init=epsg:4326")
proj4string(sumProb2) <- CRS("+init=epsg:32631")
##Create raster template to harmonise raster layers in terms of projection, resolution and extent

UTM31extent<-projectExtent(sumProb2)
crs31 <- CRS("+init=epsg:32631")
UTM31extent<-projectExtent(sumProb2, crs31)
library("rgdal")
sub6G4sP <- projectRaster(sub6G4s, UTM31extent)
##Fit Variogram to sample - convert to points first
r.spgrd<-as(sub6G4sP,"SpatialPointsDataFrame")
r.spgrd = r.spgrd[!is.na(r.spgrd[[1]]),]
selectedPoints = sample(1:length(r.spgrd[[1]]), 10000)
r.sampled = r.spgrd[selectedPoints,]
Vsub6G4s<-autofitVariogram(r.sampled)
library("gstat", lib.loc="/home/monst/R/x86_64-pc-linux-gnu-library/3.1")
gstat_variogram <- variogram(Echinocardium_flavescens ~ 1, data = rsamp1)
##Still problems with autofit so new approach
##install.packages("usdm")
library("usdm", lib.loc="/home/monst/R/x86_64-pc-linux-gnu-library/3.1")
v1 <- Variogram(sub6G4sP[[1]])
v2 <- Variogram(sub6G4sP[[2]])
v1lag <- Variogram(sub6G4sP[[1]],lag=25000,cutoff=170000)
plot(v1lag)
plot(v1lag,box=TRUE)
v2lag <- Variogram(sub6G4sP[[2]],lag=10000,cutoff=270000)
plot(v2lag,box=TRUE)
plot(v2lag)
sp <- sampleRandom(sub6G4sP[[1]],3500,sp=TRUE)
lisa1<-lisa(x=sub6G4sP[[1]],y=sp,d1=0,d2=30000,statistic="I")
lisa1sp<-sp
lisa1sp@data@lisa1<-lisa1
lisa1sp <- SpatialPointsDataFrame(sp, data.frame(lisa1 = lisa1))
plot(lisa1sp, col=lisa1+5)
sequ<-seq(10000, 200000, by=10000 )
models <- lapply(sequ, function(S) {lisa(x=sub6G4sP[[1]],y=sp,d1=0,d2=S,statistic="I")})
lisa1sp <- SpatialPointsDataFrame(sp, data.frame(sequ=models))
for (i in sequ)
mod <- lisa(x=sub6G4sP[[2]],y=sp,d1=0,d2=i,statistic="I")
names(lisa1sp)<-sequ
xsequ<-make.names(sequ)
names(lisa1sp)<-xsequ
plot(lisa1sp[20], col=lisa1sp$X150000+100)
writeOGR(lisa1sp, ".", "lisa1sp", driver="ESRI Shapefile")


###Clustering Jaccard/Sorenson on presence
##Create matrix to calculate dissimilarity
SpeciesUTM$pres <- 1
library("reshape2")
SW <- dcast(SpeciesUTM, Long + Lat ~ Species, value.var="pres")
View(SW)
##Extract coordinates
SWcoo<-SW[,1:2]
SW[,1]<-NULL
SW[,1]<-NULL
SW[is.na(SW)] <- 0








clusters[[1]][[9]]
lclust<-as.data.frame(clusters)
lclust<-plyr::ldply(clusters[[1]], rbind)
write.csv(lclustn, "llclust.csv")
llclust <- read.csv("~/Documents/ENV_VARS/Null_Model/Maxent/test/llclust.csv")
View(llclust)
llclust$x <- NULL
tlclust<-t(llclust)

m<-melt(tlclust, id=c(1,))
l1<-na.omit(m)
x1<-l1[l1$X2==1, "value"]
x2<-l1[l1$X2==2, "value"]

x3<-l1[l1$X2==3, "value"]
x4<-l1[l1$X2==4, "value"]
plot(s, (x4))
x5<-l1[l1$X2==5, "value"]
plot(s, (x5))
x6<-l1[l1$X2==6, "value"]
x7<-l1[l1$X2==7, "value"]
plot(s, (x7))
View(`l1`)
names(l1)[[4]]<-"filenn"
names(l1)[[3]]<-"filenn"
View(`l1`)
l2<-merge(x = l1, y = nnDspec, by = "filenn")


plot(l2$X2~l2$freq)
plot(l2$X2~l2$Group)
boxplot(l2$X2~l2$Group)
plot(s2, (x4))
plot(stack_simple, (x4))
plot(stack_simple, (x1))
plot(stack_simple, (x2))
plot(stack_simple, (x3))
plot(stack_simple, (x5))
plot(stack_simple, (x6))

###Recluster package not so good
##install.packages("simba")
library("simba")
SWs<-sim(SW,method="jaccard", listout=TRUE)
##Tanspose
tSW<-t(SW)
SWj <- vegdist(tSW, method = "jaccard")
fit <- hclust(SWj, method="ward.D2")
plot(fit)
fit <- hclust(SWj, method="ward.D2", n=12)
plot(fit)


tSWc<-tSW[complete.cases(tSW),]
dat<-rowSums(tSWc)
dat<-rowSums(tSWc)
tSWc<- cbind(tSWc, dat)


##problem with too few samples so remove species with less than 3
cSWc = tSWc[tSWc[,229]>=3,]
sordiss<- recluster.dist(cSWc,dist="sorensen")
points<-metaMDS(sordiss, center=TRUE)$points
col<-recluster.col(points)
recluster.plot.col(col)

library("mclust", lib.loc="/home/monst/R/x86_64-pc-linux-gnu-library/3.1")
sfit <- Mclust(col)
sfit <- Mclust(col[,1:2], G=15)
cSWc<-cbind(cSWc,sfit$classification)

##ADD Reclust classes to original file
write.csv(cSWc[,230], file="Sclust.csv")
Sclust <- read.csv("~/Documents/ENV_VARS/Null_Model/Maxent/test/Sclust.csv")
names(Sclust)[[2]]<-"filen"
names(Sclust)[[1]]<-"filen"
l2<-merge(x = l2, y = Sclust, by = "filen")
plot(l2$Sclust~l2$freq)
plot(l2$Sclust~l2$freq)
boxplot(l2$Sclust~l2$X2)
boxplot(l2$X2~l2$Sclust)
hist(l2$SX2)

###Better cluster method gives nestedness graphs etc..
install.packages("betapart")
library("betapart", lib.loc="/home/monst/R/x86_64-pc-linux-gnu-library/3.1")


#getbetapartobjects
ceram.s.core<-betapart.core(tSW)
#multiplesitemeasures
ceram.s.multi<-beta.multi(ceram.s.core)
#samplingacrossequalsites
ceram.s.samp<-beta.sample(ceram.s.core,
                          sites=10,samples=100)
#plottingthedistributionsofcomponents
dist.s<-ceram.s.samp$sampled.values
plot(density(dist.s$beta.SOR),
     xlim=c(0,0.8),ylim=c(0,19),xlab='Beta
     diversity',main='',lwd=3)
lines(density(dist.s$beta.SNE),lty=1, lwd=2)
lines(density(dist.s$beta.SIM),lty=2,lwd=2)
#pairwiseforsouth
pair.s<-beta.pair(ceram.s.core)
#plottingclusters
dist.s<-ceram.s.samp$sampled.values
plot(hclust(pair.s$beta.sim,method="average"),hang=-1,main='',sub='',xlab='')
title(xlab=expression(beta[sim]),line=0.3)
plot(hclust(pair.s$beta.sne,method="average"),hang=-1,main='',sub='',,xlab='')
title(xlab=expression(beta[sne]),line=03)

sneA<-(hclust(pair.s$beta.sne,method="average"))
sneW<-(hclust(pair.s$beta.sne,method="ward.D2"))
sneA8 <- cutree(sneA, k=8) # cut tree into 5 clusters
# draw dendogram with red borders around the 5 clusters
rect.hclust(sneA, k=8, border="red") 

pvresultA <- pvclust(SW, method.dist="binary", method.hclust="average", nboot=1000)
pvresultA <- pvclust(SW, method.dist="binary", method.hclust="ward.D2", nboot=200)
plot(pvresultA)
pvrect(pvresultA, alpha=0.95)
clustersA <- pvpick(pvresultA)
##Now transfer ward clusters to l2





pvresultA <- pvclust(SW, method.dist="binary", method.hclust="average", nboot=1000)






###Quantile regression - nice idea but didn't work had to do it manually in Gnumeric ~/monst/Documents/R/nnDist_Abundance_Calc.gnumeric


library(quantreg)
#model.rq <- rq(sd3 ~ lfreq+lfreq2+lfreq3, tau=c(0.10, 0.25, 0.75, 0.95), data=guSpec)
model.rq <- rq(sd3 ~ lfreq2, tau=c(0.15, 0.8, 0.995), data=guSpec)
quantile.regressions <- data.frame(t(coef(model.rq)))
#colnames(quantile.regressions) <- c("intercept", "slope", "slope2", "slope3")
colnames(quantile.regressions) <- c("intercept", "slope")
quantile.regressions$quantile <- rownames(quantile.regressions)
quantile.regressions

##intercept     slope  quantile
##tau= 0.25 85.63636 -1.363636 tau= 0.25

library(ggplot2)

scatterplot <- qplot(x=lfreq2, y=sd3, data=guSpec)
scatterplot + geom_abline(aes(intercept=intercept, slope=slope,
                              colour=quantile), data=quantile.regressions)

f <- rq(sd3 ~ rcs(lfreq, 2), tau=c(0.10,  0.80, 0.995), data=guSpec)

##
replace underscore with space
new_str = re.sub(r'[\w_+]', ' ', new_str)






jacdiss<- recluster.dist(cSWc,dist="jaccard")
jacpoints<-metaMDS(jacdiss, center=TRUE)$points
jcol<-recluster.col(jacpoints)
recluster.plot.col(jcol)

stree<-recluster.cons(cSWc, dist="sorensen",tr=20, p=0.5)
expl_div<-recluster.expl.diss(stree$cons,sordiss)
expl_div

#Select cut #4 and group data in RGB space
ncol<-recluster.group.col(col,expl_div$matrix[,4])


sfit <- Mclust(col)
plot(sfit) # plot results
summary(sfit) # display the best model 



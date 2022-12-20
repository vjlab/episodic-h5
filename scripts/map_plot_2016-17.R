library("seraphim")
library(plyr)
library(tidyverse)
library(rgdal)
library(diagram)
library(stringr) 
library(viridis)


# estimating dispersal statistics

localTreesDirectory = "2016-17_extracted_trees"
allTrees = scan(file="clade2344b_v1_HA_host_lat_long.burnin.0.1.trees", what="", sep="\n", quiet=TRUE)
burnIn = 0
randomSampling = FALSE
nberOfTreesToSample = 1000
mostRecentSamplingDatum = 2022.466
coordinateAttributeName = "location"

treeExtractions(localTreesDirectory, allTrees, burnIn, randomSampling, nberOfTreesToSample, mostRecentSamplingDatum, coordinateAttributeName)


# Extracting spatio-temporal information embedded in the MCC tree
source("mccExtractions.r")
mcc_tre = readAnnotatedNexus("clade2344b_v1_HA_host_lat_long.mean.mcc.tree")
mcc_tab = mccExtractions(mcc_tre, mostRecentSamplingDatum)
write.csv(mcc_tab, "2016-17_MCC.csv", row.names=F, quote=F)

# Estimating the HPD region for each time slice
nberOfExtractionFiles = nberOfTreesToSample
#prob = 0.95; 
prob = 0.8; 
prob = 1; 
precision = 0.025
startDatum = min(mcc_tab[,"startYear"])

polygons = suppressWarnings(spreadGraphic2(localTreesDirectory, nberOfExtractionFiles, prob, startDatum, precision))


# spatial boundaries
my_spdf <- readOGR(
  dsn= paste0("../seraphim/shapefile-js-master/thematicmapping/TM_WORLD_BORDERS_SIMPL-0.3.shp") ,
  layer="TM_WORLD_BORDERS_SIMPL-0.3",
  verbose=FALSE
)

# which(my_spdf@data$NAME=="Brazil")
# BRAZILBD <- my_spdf[my_spdf@data$NAME=="Brazil" , ]
plot(my_spdf)

# borders=BRAZILBD
# template_raster=BRAZILBD

borders=my_spdf
template_raster=my_spdf


#Defining the different colour scales
#minYear = 2010; maxYear = 2022.466
minYear = 2014; maxYear = 2017
mcc_tab_this<- mcc_tab %>% data.frame()
mcc_tab_this<- mcc_tab_this %>% dplyr::filter(startYear>=minYear) %>% dplyr::filter(endYear<=maxYear)
endYears_indices = (((mcc_tab_this[,"endYear"]-minYear)/(maxYear-minYear))*100)+1
min(mcc_tab_this[,"startYear"])
max(mcc_tab_this[,"endYear"])

n_number_colours_needed<- max(round(endYears_indices))
# n_repeats_discrete<- 10
# c1<- rev(brewer.pal(4,"PuRd"))
# c2<- (brewer.pal(9,"Blues"))
# colours<- rev(rep(c(c1,c2), each=n_repeats_discrete))
# colour_scale<- colorRampPalette(colours)(n_number_colours_needed)

colour_scale <- magma(n_number_colours_needed)
#colour_scale <- viridis(n_number_colours_needed)


endYears_colours = colour_scale[round(endYears_indices)]
polygons_colours = rep(NA, length(polygons))
for (i in 1:length(polygons))
{
  date = as.numeric(names(polygons[[i]]))
  
  if(date > maxYear) next
  
  polygon_index = round((((date-minYear)/(maxYear-minYear))*100)+1)
  #polygons_colours[i] = paste0(colour_scale[polygon_index],"25")
  polygons_colours[i] = paste0(substr(colour_scale[polygon_index],1,nchar(colour_scale[polygon_index])-2),"15")
}

# 5. Generating the dispersal history plot
pdf(sprintf('map_2016-2017_%s_v2.pdf',maxYear),width=16, height=15,bg="white")

ptsize<- 0.6
pitjit<- 0.08
par(mar=c(0,0,0,0), oma=c(1.2,3.5,1,0), mgp=c(0,0.4,0), lwd=0.2, bty="o")
#plot(template_raster, col="white", box=F, axes=F, colNA="grey90", legend=F)
#plot(borders, add=T, lwd=0.1, border="gray10")
plot(template_raster, col="#bcbebf", box=F, axes=F, colNA="grey90", legend=F)
plot(borders, add=T, lwd=0.1, border="black")
for (i in length(polygons):1)
{
  #if(as.numeric(names(polygons[[i]])) > 2020.8) next
  if(as.numeric(names(polygons[[i]])) > maxYear) next
  plot(polygons[[i]], axes=F, col=polygons_colours[i], add=T, border=NA)
}
for (i in 1:dim(mcc_tab_this)[1])
{
  if(mcc_tab_this[i,"startYear"] > maxYear) next
  curvedarrow(cbind(mcc_tab_this[i,"startLon"],mcc_tab_this[i,"startLat"]), cbind(mcc_tab_this[i,"endLon"],mcc_tab_this[i,"endLat"]), arr.length=0,
              arr.width=0, lwd=2.5, lty=1, lcol="grey22", arr.col=NA, arr.pos=FALSE, curve=0.3, dr=NA, endhead=F)
  curvedarrow(cbind(mcc_tab_this[i,"startLon"],mcc_tab_this[i,"startLat"]), cbind(mcc_tab_this[i,"endLon"],mcc_tab_this[i,"endLat"]), arr.length=0,
              arr.width=0, lwd=2, lty=1, lcol=endYears_colours[i], arr.col=NA, arr.pos=FALSE, curve=0.3, dr=NA, endhead=F)
}
for (i in dim(mcc_tab_this)[1]:1)
{
  if(mcc_tab_this[i,"startYear"] > maxYear) next
  xs<- mcc_tab_this[i,"startLon"]
  ys<- mcc_tab_this[i,"startLat"]
  xe<- jitter(mcc_tab_this[i,"endLon"],pitjit)
  ye<- jitter(mcc_tab_this[i,"endLat"],pitjit)
  if (i == 1)
  {
    points(xs, ys, pch=16, col=colour_scale[1], cex=ptsize)
    points(xs, ys, pch=1, col="gray10", cex=ptsize)
  }
  points(xe, ye, pch=16, col=endYears_colours[i], cex=ptsize)
  points(xe, ye, pch=1, col="gray10", cex=ptsize)
}

xrange<- c(xmin(template_raster), xmax(template_raster))
yrange<- c(ymin(template_raster), ymax(template_raster))
rect(xrange[1], yrange[1], xrange[2], yrange[2], xpd=T, lwd=0.2)
axis(1, c(ceiling(xmin(template_raster)), floor(xmax(template_raster))), pos=ymin(template_raster), mgp=c(0,0.2,0), cex.axis=0.5, lwd=0, lwd.tick=0.2, padj=-0.8, tck=-0.01, col.axis="gray30")
axis(2, c(ceiling(ymin(template_raster)), floor(ymax(template_raster))), pos=xmin(template_raster), mgp=c(0,0.5,0), cex.axis=0.5, lwd=0, lwd.tick=0.2, padj=1, tck=-0.01, col.axis="gray30")
#rast = raster(matrix(nrow=1, ncol=2)); rast[1] = min(mcc_tab_this[,"startYear"]); rast[2] = max(mcc_tab_this[,"endYear"])
rast = raster(matrix(nrow=1, ncol=2)); rast[1] = min(mcc_tab_this[,"startYear"]); rast[2] = 2017
plot(rast, legend.only=T, add=T, col=colour_scale, legend.width=0.5, legend.shrink=0.3, smallplot=c(0.40,0.80,0.14,0.155),
     legend.args=list(text="", cex=0.7, line=0.3, col="gray30"), horizontal=T,
     axis.args=list(cex.axis=0.6, lwd=0, lwd.tick=0.2, tck=-0.5, col.axis="gray30", line=0, mgp=c(0,-0.02,0)))

a<-dev.off()


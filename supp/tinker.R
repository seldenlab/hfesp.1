# Load e.biface image
e.biface <- loadImage("out/TARL-41BW4-6-2-67.png", 
                    CIELab = FALSE,
                    lower = FALSE,
                    upper = FALSE)

# plot in CIE Lab color space
plotPixels(e.biface,
           color.space = "rgb", 
           n = TRUE,
           main = "Pixels (RGB Color Space)")

# plot clusters in 3D space
scatter3dclusters(e.biface.hist, scaling=22, opacity=0.99, 
                  plus=0.03, y.margin.add=1,
                  type="h", lwd=1.5, angle=50, main="Clusters (RGB)")

# get histogram clusters
e.biface.hist <- getImageHist(e.biface,
                            bins = 3,
                            sample.size = TRUE,
                            lower = FALSE,
                            upper = FALSE,
                            hsf = FALSE,
                            title = "Histogram color clusters",
                            ylab="Proportion of image")


# compare colors----
library(here)
library(colordistance)
library(recolorize)
library(tidyverse)
library(vegan)
library(ape)
library(ggConvexHull)
library(wesanderson)

# get image paths
images <- colordistance::getImagePaths('out')

# generate color distance matrix
cdm <- imageClusterPipeline(images,
                            cluster.method = 'hist',
                            distance.method = 'emd',
                            lower = FALSE,
                            upper = FALSE,
                            hist.bins = 3,
                            color.space = 'lab',
                            ref.white = 'D65',
                            plot.bins = TRUE,
                            pausing = FALSE
)

# replace NAs with 0
cdm[is.na(cdm)] = 0
# export color distance matrix
write.csv(cdm, file = 'distance_matrix.csv')

#preliminaries
# import categorical data
hfesp.data <- read.csv('qdata.csv',
                       header = TRUE,
                       as.is = TRUE)

# import color distance matrix
c.data <- read.csv('distance_matrix.csv',
                   header = TRUE,
                   as.is = TRUE)

# cdm rev sans .png
c.data2 <- as.data.frame(c.data %>% 
                           select(c(CNO.Pohler.2002.01.18.png:UArk_BattleMound.55.15.9.png)) %>% 
                           rename('CNO-Pohler-2002-01-18' = CNO.Pohler.2002.01.18.png,
                                  'CNO-Pohler-2002-01-20' = CNO.Pohler.2002.01.20.png,
                                  'CNO-Pohler-2002-01-23' = CNO.Pohler.2002.01.23.png,
                                  'CNO-Pohler-2002-01-27' = CNO.Pohler.2002.01.27.png,
                                  'LSEM-Haley-HFE1' = LSEM.Haley.HFE1.png,
                                  'LSEM-Haley-HFE2' = LSEM.Haley.HFE2.png,
                                  'LSEM-Haley-HFE3' = LSEM.Haley.HFE3.png,
                                  'LSEM-Haley-HFE4' = LSEM.Haley.HFE4.png,
                                  'LSEM-Haley-HFE5' = LSEM.Haley.HFE5.png,
                                  'NSU-AllenPlantation-142' = NSU.AllenPlantation.142.png,
                                  'NSU-Belcher-404' = NSU.Belcher.404.png,
                                  'NSU-Gahagan-955' = NSU.Gahagan.955.png,
                                  'NSU-Gahagan-956' = NSU.Gahagan.956.png,
                                  'NSU-SmithportLanding-152' = NSU.SmithportLanding.152.png,
                                  'NSU-SmithportLanding-157' = NSU.SmithportLanding.157.png,
                                  'NSU-SmithportLanding-95' = NSU.SmithportLanding.95.png,
                                  'NSU-SmithportLanding-96' = NSU.SmithportLanding.96.png,
                                  'TARL-41BW3-6-1-7-FS7' = TARL.41BW3.6.1.7.FS7.png,
                                  'TARL-41BW4-341-427' = TARL.41BW4.341.427.png,
                                  'TARL-41BW4-341-464' = TARL.41BW4.341.464.png,
                                  'TARL-41BW4-6-2-132' = TARL.41BW4.6.2.132.png,
                                  'TARL-41BW4-6-2-67' = TARL.41BW4.6.2.67.png,
                                  'TARL-41BW4-6-2-78' = TARL.41BW4.6.2.78.png,
                                  'TARL-41CE19-2015-1' = TARL.41CE19.2015.1.png,
                                  'TARL-41RR3-2' = TARL.41RR3.2.png,
                                  'UArk-Burrows-55-19-2' = UArk.Burrows.55.19.2.png,
                                  'UArk-Crenshaw-55-1-10' = UArk.Crenshaw.55.1.10.png,
                                  'UArk-Crenshaw-55-1-4' = UArk.Crenshaw.55.1.4.png,
                                  'UArk-Crenshaw-55-1-43' = UArk.Crenshaw.55.1.43.png,
                                  'UArk_BattleMound-55-15-9' = UArk_BattleMound.55.15.9.png)
)

## Neighbor-Joining Tree----
# define add_image function
add_image <- function(obj,
                      x = NULL,
                      y = NULL,
                      width = NULL,
                      interpolate = TRUE,
                      angle = 0) {
  
  # get current plotting window parameters:
  usr <- graphics::par()$usr
  pin <- graphics::par()$pin
  
  # image dimensions and scaling factor:
  imdim <- dim(obj)
  sf <- imdim[1] / imdim[2]
  
  # set the width of the image (relative to x-axis)
  w <- width / (usr[2] - usr[1]) * pin[1]
  h <- w * sf
  hu <- h / pin[2] * (usr[4] - usr[3])
  
  # plot the image
  graphics::rasterImage(image = obj,
                        xleft = x - (width / 2), xright = x + (width / 2),
                        ybottom = y - (hu / 2), ytop = y + (hu/2),
                        interpolate = interpolate,
                        angle = angle)
}

# c.data2 sans 'X'
c.data2 <- c.data2 %>% 
  select('CNO-Pohler-2002-01-18':'UArk_BattleMound-55-15-9')
# column names to ID
rownames(c.data2) <- colnames(c.data2)

# neighbor-joining tree----
tree <- nj(as.dist(c.data2))
# plot tree
plot(tree,
     show.tip.label = FALSE,
     direction = 'upward')

# get parameters from plotting environment
lastPP <- get('last_plot.phylo', 
              envir = .PlotPhyloEnv)

# get xy coordinates of tips
ntip <- lastPP$Ntip

# first n values are tips/remaining values are node coordinates
xy <- data.frame(x = lastPP$xx[1:ntip],
                 y = lastPP$yy[1:ntip])

# get image names
imnames <- tools::file_path_sans_ext(basename(images))

# get tip labels
tipnames <- tree$tip.label

# check that images are in correct order:
image_order <- match(tipnames, imnames)
images <- images[image_order]

## final plot
par(mar = rep(0, 4))
plot(tree, show.tip.label = FALSE, direction = 'upward')

# parameters from plotting environment
lastPP <- get('last_plot.phylo', envir = .PlotPhyloEnv)

# xy coordinates of tips
ntip <- lastPP$Ntip
xy <- data.frame(x = lastPP$xx[1:ntip],
                 y = lastPP$yy[1:ntip]) 

for (i in 1:length(images)) {
  add_image(png::readPNG(images[i]),
            x = xy[i, 1],
            y = xy[i, 2], 
            width = 2)
}

## NMDS Ordination----
nmds_scores <- scores(metaMDS(comm = as.dist(c.data2)))

# set plot parameters
plot(nmds_scores,
     xlim = c(-24, 33),
     ylim = c(-5.5, 8),
     cex = 1
)

# Convex hulls by firing method
ordihull(
  nmds_scores,
  hfesp.data$ColorGroup,
  draw = 'polygon',
  col = c("#798E87","#C27D38"),
  border = c("#798E87","#C27D38")
)

# add images
for (i in 1:length(images)) {
  
  # read image:
  img <- png::readPNG(images[i]) 
  # add image:
  add_image(img,
            x = nmds_scores[i, 1],
            y = nmds_scores[i, 2],
            width = 2)
}

## Permutational MANOVA----
c.data2$X <- colnames(c.data2)
bottle <- full_join(hfesp.data, c.data2, by = 'X') # left join by specimen number
bottle.dist <- bottle[,6:35] # color distance matrix
set.seed(10) # make results reproducible

# model: color as a function of firing category (oxidized v. reduced)
bottle.colour <- adonis2(bottle.dist ~ ColorGroup,
                         data = bottle,
                         permutations = 10000,
                         method = 'bray')

## does color differ by firing category?
bottle.colour

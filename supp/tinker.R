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

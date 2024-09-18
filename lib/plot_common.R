
# color transparency
addAlphaToColor <- Vectorize(function(color, alpha = 0.25) {
    x <- col2rgb(color)
    rgb(x[1], x[2], x[3], max = 255, alpha = alpha * 255)
})
cnColors <- c("black","blue","darkcyan","red3","purple3","orange3") # CN  to 5
cnColorsAlpha <- addAlphaToColor(cnColors)

# set common image properties
width     <- 2.5
height    <- 2.8
units     <- 'in' # w and h in inches
pointsize <- 8
resol     <- 900 #dpi
pch       <- "." #20
cex       <- 1 #0.3

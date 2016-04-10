library(gplots)
library(ggplot2)
library(reshape2)
library(pvclust)
library(RColorBrewer);
library(grDevices)
require(scales)

pngname= "output.png"
E = read.delim("hist2.txt",row.names=1,header=T)

png(pngname,width=1000,height=1000)
par(mfrow=c(1,1))
par(mar=c(5,1,1,1))

X = as.numeric(rownames(E))

plot(X,E[,1],type="l",lwd=5)

dev.off()

stop("bawk")

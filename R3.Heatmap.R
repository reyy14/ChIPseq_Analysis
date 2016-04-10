library(gplots)
#library(ggplot2)
library(reshape2)
library(pvclust)
library(RColorBrewer);
library(grDevices)
require(scales)
colorL=c("red","purple","blue","yellow","green","orange","brown","gray","black","coral","beige","cyan","pink","khaki","magenta")
#clusterNum=12

profilename = "ghist.txt"
pngname = "heatmap.png"
E = read.delim(profilename,header=F)
gene <- rownames(E)
rowsum = apply(E[,2:dim(E)[2]],1,sum)
orderedE <- E[order(rowsum),]
E <- orderedE

png(pngname,width=300,height=800)
par(mfrow=c(1,1))
my <- E
par(mar=c(1,1,5,1))
mylist = c("H3S10ac_A_8_23")

Dat = my[,2:ncol(my)]
pal.1 = colorRampPalette(c("white","red"),space="rgb")
mydat = t(as.matrix(Dat))
mymin = quantile(mydat,0.05)
mymax = quantile(mydat,0.9)
mydat[mydat<mymin]=mymin
mydat[mydat>mymax]=mymax
mydat = rescale(mydat)
image(mydat,col=pal.1(30),axes=F,main="H3S10ac_A_8_23",cex.main=3)
box()

dev.off()

stop("bawk")

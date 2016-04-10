library(gplots)
#library(ggplot2)
library(reshape2)
library(pvclust)
library(RColorBrewer);
library(grDevices)
require(scales)
library(plyr)
colorL=c("red","purple","blue","yellow","green","orange","brown","gray","black","coral","beige","cyan","pink","khaki","magenta")
#clusterNum=12
png("heatmap_peak_revised.png",width=1500,height=1000)
par(mfrow=c(1,5))
par(mar=c(10,2,5,1))
profilelist = c("ghist_peak/H3S10ac_A_8_23_ghist.txt","ghist_peak/H3K9ac_ghist.txt","ghist_peak/H3K9me2_ghist.txt","ghist_peak/H3K9me3_ghist.txt","ghist_peak/H3K27ac_ghist.txt")

sortguide = read.delim("ghist_peak/H3S10ac_A_8_23_ghist.txt",header=F)
rowsum = apply(sortguide[,2:dim(sortguide)[2]],1,sum)
ordered = sortguide[order(rowsum),]
sorted = ordered

for (i in 1:5){
	if (i == 1) {xl="H3S10ac"}
	if (i == 2) {xl="H3K9ac"}
	if (i == 3) {xl="H3K9me2"}
	if (i == 4) {xl="H3K9me3"}
	if (i == 5) {xl="H3K27ac"}

	profilename = profilelist[i]
	E = read.delim(profilename,header=F)
	#rowsum = apply(E[,2:dim(E)[2]],1,sum)
	#orderedE <- E[order(rowsum),]
	#E <- orderedE
	my <- E
	#mylist = c("ahahahahahaha")
	
	my <- my[match(sorted[,1],my[,1]),]
	
	my <- my[apply(my,1,function(x) return(any(!is.na(x)))),]

	Dat = my[,2:ncol(my)]
	pal.1 = colorRampPalette(c("white","red"),space="rgb")
	mydat = t(as.matrix(Dat))
	mymin = quantile(mydat,0.05)
	mymax = quantile(mydat,0.9)
	mydat[mydat<mymin]=mymin
	mydat[mydat>mymax]=mymax
	mydat = rescale(mydat)
	image(mydat,col=pal.1(30),axes=F,xlab=xl,cex.lab=4,cex.main=3,mgp=c(6,3,0))
	if (i == 1) {axis(1,at=c(0,0.5,1),labels=c("-2k","0","2k"),cex.axis=3,line=1.5,lwd=0,lwd.ticks=0)}
	box()

	
}
dev.off()


stop("bawk")

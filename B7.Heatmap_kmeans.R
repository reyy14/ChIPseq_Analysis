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
png("heatmap_compiled_new_trimmed.png",width=3000,height=1000)
par(mfrow=c(1,7))
par(mar=c(10,2,5,1))
profilelist = c("ghist/H3S10ac_A_8_23_ghist.txt","ghist/H3K27me3_ghist1.txt","ghist/H3K9ac_ghist.txt","ghist/SRR1130791_ghist.txt","ghist/H3K9me3_ghist.txt","ghist/H3K27ac_ghist.txt","ghist/H3K4me3_ghist1.txt")

sortguide = read.delim("ghist/H3K27me3_ghist1.txt",header=F)
rowsum = apply(sortguide[,2:dim(sortguide)[2]],1,sum)
ordered = sortguide[order(rowsum),]
sorted = ordered

S10ac = read.delim("ghist/H3S10ac_A_8_23_ghist.txt",header=F)
K27me3 = read.delim("ghist/H3K27me3_ghist1.txt",header=F)
K4me3 = read.delim("ghist/H3K4me3_ghist1.txt",header=F)

K27me3 = K27me3[match(S10ac[,1],K27me3[,1]),]
K4me3 = K4me3[match(S10ac[,1],K4me3[,1]),]
K27me3 <- K27me3[apply(K27me3,1,function(x) return(any(!is.na(x)))),]
K4me3 <- K4me3[apply(K4me3,1,function(x) return(any(!is.na(x)))),]

S10ac = S10ac[match(K27me3[,1],S10ac[,1]),]
K4me3 = K4me3[match(K27me3[,1],K4me3[,1]),]
S10ac <- S10ac[apply(S10ac,1,function(x) return(any(!is.na(x)))),]
K4me3 <- K4me3[apply(K4me3,1,function(x) return(any(!is.na(x)))),]

S10ac = S10ac[match(K4me3[,1],S10ac[,1]),]
K27me3 = K27me3[match(K4me3[,1],K27me3[,1]),]
K27me3 <- K27me3[apply(K27me3,1,function(x) return(any(!is.na(x)))),]
S10ac <- S10ac[apply(S10ac,1,function(x) return(any(!is.na(x)))),]

mat = cbind(apply(S10ac[,2:dim(S10ac)[2]],1,function(x) mean(x)),apply(K4me3[,2:dim(K4me3)[2]],1,function(x) mean(x)), apply(K27me3[,2:dim(K27me3)[2]],1,function(x) mean(x)))

rownames(mat) <- S10ac[,1]
colnames(mat) <- c("S10ac","K4me3","K27me3")

write.table(mat,file="data.txt",quote=F,sep="\t")

#K27me3 <- K27me3[apply(K27me3,1,function(x) return(any(!is.na(x)))),]
#K4me3 <- K4me3[apply(K4me3,1,function(x) return(any(!is.na(x)))),]

H3S10ac = read.delim("ghist/H3S10ac_A_8_23_ghist.txt",header=F)
H3S10ac = H3S10ac[match(sorted[,1],H3S10ac[,1]),]
H3S10ac = H3S10ac[apply(H3S10ac,1,function(x) return(any(!is.na(x)))),]
H3S10ac = H3S10ac[apply(H3S10ac[,2:dim(H3S10ac)[2]],1,function(x) return(all(x!=0))),]

for (i in 1:7){
	if (i == 1) {xl="H3S10ac"}
	if (i == 2) {xl="H3K27me3"}
	if (i == 3) {xl="H3K9ac"}
	if (i == 4) {xl="H3K9me2"}
	if (i == 5) {xl="H3K9me3"}
	if (i == 6) {xl="H3K27ac"}
	if (i == 7) {xl="H3K4me3"}

	profilename = profilelist[i]
	E = read.delim(profilename,header=F)
	#rowsum = apply(E[,2:dim(E)[2]],1,sum)
	#orderedE <- E[order(rowsum),]
	#E <- orderedE
	my <- E
	#mylist = c("ahahahahahaha")
	
	my <- my[match(sorted[,1],my[,1]),]
	my <- my[match(H3S10ac[,1],my[,1]),]
	
	my <- my[apply(my,1,function(x) return(any(!is.na(x)))),]

	Dat = my[,2:ncol(my)]
	pal.1 = colorRampPalette(c("white","red"),space="rgb")
	mydat = t(as.matrix(Dat))
	mymin = quantile(mydat,0.05)
	mymax = quantile(mydat,0.9)
	mydat[mydat<mymin]=mymin
	mydat[mydat>mymax]=mymax
	mydat = rescale(mydat)
	image(mydat,col=pal.1(30),axes=F,xlab=xl,cex.lab=4,cex.main=3,mgp=c(6,3,1))
	if (i ==1) {axis(1,at=c(0,0.5,1),labels=c("-2k","0","2k"),cex.axis=3,line=1.5,lwd=0,lwd.ticks=0)}
	box()

	
}
dev.off()


stop("bawk")

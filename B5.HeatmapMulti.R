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
png("heatmap_ok.png",width=3200,height=1000)
par(mfrow=c(1,7))
par(mar=c(10,2,5,1))
profilelist = c("ghist/H3K9ac_ghist.txt","ghist/H3K27ac_ghist.txt","ghist/SRR1130791_ghist.txt","ghist/H3K9me3_ghist.txt","ghist/H3S10ac_A_8_23_ghist.txt","ghist/H3K27me3_ghist1.txt","ghist/H3K4me3_ghist1.txt")


sortguide = read.delim("ghist/H3S10ac_A_8_23_ghist.txt",header=F)
H3S10ac = read.delim("ghist/H3S10ac_A_8_23_ghist.txt",header=F)
H3K4me3 = read.delim("ghist/H3K4me3_ghist1.txt",header=F)
H3K27me3 = read.delim("ghist/H3K27me3_ghist1.txt",header=F)

for (i in 1:7){
	da = read.delim(profilelist[i],header=F)
	sortguide = sortguide[sortguide[,1] %in% da[,1],]
}

#sortguide = sortguide[sortguide[,1] %in% H3K27me3[,1],]
H3K27me3 = H3K27me3[H3K27me3[,1] %in% sortguide[,1],]
H3S10ac = H3S10ac[H3S10ac[,1] %in% H3K27me3[,1],]
H3K4me3 = H3K4me3[H3K4me3[,1] %in% H3K27me3[,1],]

rownames(sortguide) = c(1:nrow(sortguide))
rownames(H3K27me3) = c(1:nrow(H3K27me3))
rownames(H3S10ac) = c(1:nrow(H3S10ac))
rownames(H3K4me3) = c(1:nrow(H3K4me3))

rowsum = apply(sortguide[,2:dim(sortguide)[2]],1,sum)
ordered = sortguide[order(rowsum),]
sorted = ordered

H3S10ac = H3S10ac[rownames(sorted),]
H3S10ac_1 = H3S10ac[(nrow(H3S10ac)-24999):nrow(H3S10ac),]
H3S10ac_2 = H3S10ac[1:(nrow(H3S10ac)-25000),]

H3K27me3 = H3K27me3[rownames(sorted),]
H3K27me3_1 = H3K27me3[(nrow(H3K27me3)-24999):nrow(H3K27me3),]
H3K27me3_2 = H3K27me3[1:(nrow(H3K27me3)-25000),]

H3K4me3 = H3K4me3[rownames(sorted),]
H3K4me3_1 = H3K4me3[(nrow(H3K4me3)-24999):nrow(H3K4me3),]
H3K4me3_2 = H3K4me3[1:(nrow(H3K4me3)-25000),]

rowsum_k27 = apply(H3K27me3_1[,2:dim(H3K27me3_1)[2]],1,sum)
ordered_k27 = H3K27me3_1[order(rowsum_k27),]
sorted_k27 = ordered_k27

sorted_k4 = H3K4me3_1[rownames(sorted_k27),]
sorted_k4 = sorted_k4[apply(sorted_k4,1,function(x) return(all(!is.na(x)))),]

sorted_s10 = H3S10ac_1[rownames(sorted_k27),]
sorted_s10 = sorted_s10[apply(sorted_s10,1,function(x) return(all(!is.na(x)))),]

s10 = rbind(H3S10ac_2,sorted_s10)
k27 = rbind(H3K27me3_2,sorted_k27)
k4 = rbind(H3K4me3_2,sorted_k4)

Dat = s10[,2:ncol(s10)]
pal.1 = colorRampPalette(c("white","red"),space="rgb")
mydat = t(as.matrix(Dat))
mymin = quantile(mydat,0.05)
mymax = quantile(mydat,0.9)
mydat[mydat<mymin]=mymin
mydat[mydat>mymax]=mymax
mydat = rescale(mydat)
image(mydat,col=pal.1(30),axes=F,xlab="H3S10ac",cex.lab=4,cex.main=3,mgp=c(6,3,1))
axis(1,at=c(0,0.5,1),labels=c("-2k","0","2k"),cex.axis=3,line=1.5,lwd=0,lwd.ticks=0)
box()

Dat = k4[,2:ncol(k4)]
pal.1 = colorRampPalette(c("white","red"),space="rgb")
mydat = t(as.matrix(Dat))
mymin = quantile(mydat,0.05)
mymax = quantile(mydat,0.9)
mydat[mydat<mymin]=mymin
mydat[mydat>mymax]=mymax
mydat = rescale(mydat)
image(mydat,col=pal.1(30),axes=F,xlab="H3K4me3",cex.lab=4,cex.main=3,mgp=c(6,3,1))
box()

Dat = k27[,2:ncol(k27)]
pal.1 = colorRampPalette(c("white","red"),space="rgb")
mydat = t(as.matrix(Dat))
mymin = quantile(mydat,0.05)
mymax = quantile(mydat,0.9)
mydat[mydat<mymin]=mymin
mydat[mydat>mymax]=mymax
mydat = rescale(mydat)
image(mydat,col=pal.1(30),axes=F,xlab="H3K27me3",cex.lab=4,cex.main=3,mgp=c(6,3,1))
box()
#dev.off()
#stop("bawk")
for (i in 1:4){
	if(i == 1){xl = "H3K9ac"}
	if(i == 2){xl = "H3K27ac"}
	if(i == 3){xl = "H3K9me2"}
	if(i == 4){xl = "H3K9me3"}

	data = read.delim(profilelist[i],header=F)
	data = data[data[,1] %in% H3K27me3[,1],]
	rownames(data) = c(1:nrow(data))
	data = data[rownames(sorted),]
	data_1 = data[(nrow(data)-24999):nrow(data),]
	data_2 = data[1:(nrow(data)-25000),]
	sorted_data = data_1[rownames(sorted_k27),]
	sorted_data = sorted_data[apply(sorted_data,1,function(x) return(all(!is.na(x)))),]
	Data = rbind(data_2,sorted_data)

	Dat = Data[,2:ncol(Data)]
	pal.1 = colorRampPalette(c("white","red"),space="rgb")
	mydat = t(as.matrix(Dat))
	mymin = quantile(mydat,0.05)
	mymax = quantile(mydat,0.9)
	mydat[mydat<mymin]=mymin
	mydat[mydat>mymax]=mymax
	mydat = rescale(mydat)
	image(mydat,col=pal.1(30),axes=F,xlab=xl,cex.lab=4,cex.main=3,mgp=c(6,3,1))
	box()

}


dev.off()


stop("bawk")

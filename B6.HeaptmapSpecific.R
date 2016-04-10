library(gplots)
library(reshape2)
library(pvclust)
library(RColorBrewer);
library(grDevices)
require(scales)

colorL=c("red","purple","blue","yellow","green","orange","brown","gray","black","coral","beige","cyan","pink","khaki","magenta")

profilelist = c("ghist/H3S10ac_A_8_23_ghist.txt","ghist/H3K9ac_ghist.txt")


S10ac = read.delim(profilelist[1],header=F,row.names=1)
K9ac = read.delim(profilelist[2],header=F,row.names=1)

png("S10_K9.png",width=600,height=1000)
par(mfrow=c(3,2))
par(mar=c(5,1,5,1))

cutoff = 1

tf = rownames(K9ac) %in% rownames(S10ac)
K9ac = K9ac[tf,]
tf = rownames(S10ac) %in% rownames(K9ac)
S10ac = S10ac[tf,]

##################################

s10_cutoff = S10ac[apply(S10ac,1,function(x) return(any(x>cutoff))),]
k9_temp = K9ac[rownames(K9ac) %in% rownames(S10ac),]
k9_temp = k9_temp[apply(k9_temp,1,function(x) return(any(x>cutoff))),]
s10_cutoff = s10_cutoff[rownames(s10_cutoff) %in% rownames(k9_temp),]

rowsum = apply(s10_cutoff,1,sum)
ordered = s10_cutoff[order(rowsum),]
s10_cutoff = ordered

k9_temp_order = k9_temp[match(rownames(s10_cutoff),rownames(k9_temp)),]

Dat = s10_cutoff
pal.1 = colorRampPalette(c("white","red"),space="rgb")
mydat = t(as.matrix(Dat))
mymin = quantile(mydat,0.05)
mymax = quantile(mydat,0.9)
mydat[mydat<mymin]=mymin
mydat[mydat>mymax]=mymax
mydat = rescale(mydat)
image(mydat,col=pal.1(30),axes=F,xlab="S10",cex.lab=2,cex.main=3)
box()

Dat = k9_temp_order
pal.1 = colorRampPalette(c("white","red"),space="rgb")
mydat = t(as.matrix(Dat))
#mymin = quantile(mydat,0.05)
#mymax = quantile(mydat,0.9)
mydat[mydat<mymin]=mymin
mydat[mydat>mymax]=mymax
mydat = rescale(mydat)
image(mydat,col=pal.1(30),axes=F,xlab="K9",cex.lab=2,cex.main=3)
box()


uniq = apply(s10_cutoff,1,max) > 20*(apply(k9_temp_order,1,max))
s10_u = s10_cutoff[uniq,]
k9_temp_nu = k9_temp_order[uniq,]

Dat = s10_u
pal.1 = colorRampPalette(c("white","red"),space="rgb")
mydat = t(as.matrix(Dat))
#mymin = quantile(mydat,0.05)
#mymax = quantile(mydat,0.9)
mydat[mydat<mymin]=mymin
mydat[mydat>mymax]=mymax
mydat = rescale(mydat)
image(mydat,col=pal.1(30),axes=F,xlab="S10_Unique",cex.lab=2,cex.main=3)
box()

Dat = k9_temp_nu
pal.1 = colorRampPalette(c("white","red"),space="rgb")
mydat = t(as.matrix(Dat))
#mymin = quantile(mydat,0.05)
#mymax = quantile(mydat,0.9)
mydat[mydat<mymin]=mymin
mydat[mydat>mymax]=mymax
mydat = rescale(mydat)
image(mydat,col=pal.1(30),axes=F,xlab="K9 of S10_Unique",cex.lab=2,cex.main=3)
box()

#
k9_cutoff = K9ac[apply(K9ac,1,function(x) return(any(x>cutoff))),]
s10_temp = S10ac[rownames(S10ac) %in% rownames(K9ac),]
s10_temp = s10_temp[apply(s10_temp,1,function(x) return(any(x>cutoff))),]
k9_cutoff = k9_cutoff[rownames(k9_cutoff) %in% rownames(s10_temp),]
s10_temp_order = s10_temp[match(rownames(k9_cutoff),rownames(s10_temp)),]

uniq = apply(k9_cutoff,1,max) > 20*(apply(s10_temp_order,1,max))
k9_u = k9_cutoff[uniq,]
s10_temp_nu = s10_temp_order[uniq,]


Dat = k9_u
pal.1 = colorRampPalette(c("white","red"),space="rgb")
mydat = t(as.matrix(Dat))
#mymin = quantile(mydat,0.05)
#mymax = quantile(mydat,0.9)
mydat[mydat<mymin]=mymin
mydat[mydat>mymax]=mymax
mydat = rescale(mydat)
image(mydat,col=pal.1(30),axes=F,xlab="K9_Unique",cex.lab=2,cex.main=3)
box()

Dat = s10_temp_nu
pal.1 = colorRampPalette(c("white","red"),space="rgb")
mydat = t(as.matrix(Dat))
#mymin = quantile(mydat,0.05)
#mymax = quantile(mydat,0.9)
mydat[mydat<mymin]=mymin
mydat[mydat>mymax]=mymax
mydat = rescale(mydat)
image(mydat,col=pal.1(30),axes=F,xlab="S10 of K9_Unique",cex.lab=2,cex.main=3)
box()

dev.off()

stop("bawk")

# Samantha Ho
# CL only bulkRNAseq denovolyzer result scatter plot
# 2/10


library("readxl")
library(tidyverse)
df2 <- read_csv("UP_DenovolyzerResults_bulksc.csv")

labels.tmp1 <- cumsum(as.vector(table(df2$CHR)))
labels.mid <- as.vector(table(df2$CHR)/2)
tickpos <- c(0,labels.tmp1 - labels.mid)

cols1 <-data.frame(num=as.character(1:23),lab=c(1:23), color=c(rep(c("#8E338A","#d298da"),11),"#8E338A"), stringsAsFactors=FALSE)
cols2 <-data.frame(num=as.character(1:23),lab=c(1:23), color2=c(rep(c("#E92622","#f06c69"),11),"#E92622"), stringsAsFactors=FALSE)

df2 <- merge(df2,cols1,by.x="CHR",by.y="num",all.x=TRUE)
df2 <- merge(df2,cols2,by.x="CHR",by.y="num",all.x=TRUE)
df2 <- df2[order(df2$CHR,df2$TSS_START),]

df2$order <- seq.int(nrow(df2))

df2.LoF <- df2[which(df2$Neg.Log.pvalue.LoF>=1),]

# Making blank plot
plot(NA,NA,xlim=c(1,67),ylim=c(-6,6),axes=FALSE,xlab='',ylab='')
axis(side=2,at=c(0.1,2,4,6),labels=c(0,2,4,6),lwd=1,lwd.ticks=1,cex.axis=.8)
axis(side=2,at=c(-0.1,-2,-4,-6),labels=c('',2,4,6),lwd=1,lwd.ticks=1,cex.axis=.8)
mtext(side=2,text=expression('  -log10(P) \n LOF DNMs'),line=2,at=3.5,cex=0.8)
mtext(side=2,text=expression('            -log10(P) \n Protein-Altering DNMs'),line=2,at=-3,cex=0.8)
axis(1,at=c(tickpos),labels=NA,cex.axis=0.7,las=1,pos=0,lwd=1,lwd.ticks=1,tck=0.01)
axis(1,at=c(tickpos),labels=c("",1,2,5,6,7,9,10,11,12,15,16,17,19,20,21),cex.axis=0.46,las=1,pos=1.2,lwd=0)


# Plot points 
points(df2.LoF$order,df2.LoF$Neg.Log.pvalue.LoF,xlab="Chromosome",ylab="-log10(P)",pch=19, col=df2.LoF$color2) 
points(df2$order,-df2$Neg.Log.pvalue.proteinaltering,xlab="Chromosome",ylab="-log10(P)",pch=19,col=df2$color)
abline(h=-log10(0.0017),col="black",lty=3)
abline(h=-1*(-log10(0.0017)),col="black",lty=3)

# Label significant points 

## LOF
with(subset(df2, Neg.Log.pvalue.LoF>3.12), text(13,3.7,"SUSD1",cex=0.6,font=4))
with(subset(df2, Neg.Log.pvalue.LoF>3.12), text(3,3.7,"ZSWIM5",cex=0.6,font=4))
with(subset(df2, Neg.Log.pvalue.LoF>3.12), text(8,3,"PLEKHA6",cex=0.6,font=4))

## Protein altering
with(subset(df2, Neg.Log.pvalue.proteinaltering>3.12), text(25,-3.6,"UBALD1",cex=0.6,font=4))
with(subset(df2, Neg.Log.pvalue.proteinaltering>3.12), text(21.5,-3,"METRN",cex=0.6,font=4))
with(subset(df2, Neg.Log.pvalue.proteinaltering>3.12), text(29,-3,"CENPV",cex=0.6,font=4))





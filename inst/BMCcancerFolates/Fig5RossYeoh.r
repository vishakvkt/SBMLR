library(SBMLR)
setwd(file.path(.path.package("SBMLR"), "BMCcancerFolates")) #default dump site 
#setwd("C:/cwru/active/Morrison")  # set this to where figs should be dumped, with comment removed
library(rossEset)
pD=subset(pData(ross),subset=(type=="TEL.AML1")|(type=="T.Cell")|(type=="BCR.ABL"))
pD$type=factor(pD$type, levels=c("BCR.ABL","TEL.AML1","T.Cell"))
pD$type
as.numeric(pD$type)
I=order(pD$type)
pD=pD[I,]

cols=1+c(1,3,2)[as.numeric(pD$type)]
pchs=c("b","B","T")[as.numeric(pD$type)]


par(mfrow=c(2,2))

load("rossExprs.Rdata") 
attach(DF)


plot(TYMS,MTHFD1, pch=pchs,col=cols,cex=.6,ylim=c(0,3000),xlim=c(0,30000))
title(main="Ross et al. ",line=1,cex.main = 1,   font.main= 3, col.main= "black")
x=DF
r <- cor.test(x[cols==2,"TYMS"], x[cols==2,"MTHFD1"])
txt1 <- paste("r=",format(c(r$estimate, 0.123456789), digits=2)[1]," P=",format(c(r$p.value, 0.123456789), digits=2)[1], sep="")
text(0, 3000, txt1, col=2,cex = .7,pos=4)
r <- cor.test(x[cols==3,"TYMS"], x[cols==3,"MTHFD1"])
txt2 <- paste("r=",format(c(r$estimate, 0.123456789), digits=2)[1]," P=",format(c(r$p.value, 0.123456789), digits=2)[1], sep="")
text(0, 2700, txt2, col=3,cex = .7,pos=4)
r <- cor.test(x[cols==4,"TYMS"], x[cols==4,"MTHFD1"])
txt3 <- paste("r=",format(c(r$estimate, 0.123456789), digits=2)[1]," P=",format(c(r$p.value, 0.123456789), digits=2)[1], sep="")
text(0, 2850, txt3, col=4,cex = .7,pos=4)
detach(DF)

load("yeohExprs.Rdata") 
attach(DF)
plot(TYMS,MTHFD1, pch=pchs,col=cols,cex=.6,ylim=c(0,2500),xlim=c(0,12000))
title(main="Yeoh et al. ",line=1,cex.main = 1,   font.main= 3, col.main= "black")
x=DF
r <- cor.test(x[cols==2,"TYMS"], x[cols==2,"MTHFD1"])
txt1 <- paste("r=",format(c(r$estimate, 0.123456789), digits=2)[1]," P=",format(c(r$p.value, 0.123456789), digits=2)[1], sep="")
text(0, 2500, txt1, col=2,cex = .7,pos=4)
r <- cor.test(x[cols==3,"TYMS"], x[cols==3,"MTHFD1"])
txt2 <- paste("r=",format(c(r$estimate, 0.123456789), digits=2)[1]," P=",format(c(r$p.value, 0.123456789), digits=2)[1], sep="")
text(0, 2200, txt2, col=3,cex = .7,pos=4)
r <- cor.test(x[cols==4,"TYMS"], x[cols==4,"MTHFD1"])
txt3 <- paste("r=",format(c(r$estimate, 0.123456789), digits=2)[1]," P=",format(c(r$p.value, 0.123456789), digits=2)[1], sep="")
text(0, 2350, txt3, col=4,cex = .7,pos=4)
detach(DF)

load("FmorrRossAll.Rdata") 
DF=flux[1:49,c("MTHFD","TYMS")]
attach(DF)
plot(TYMS,MTHFD/2, pch=pchs,col=cols,cex=.6,xlab="dTMP Flux (uM/hr)",ylab="DNPS Flux (uM/hr)",ylim=c(0,1300),xlim=c(0,14))
title(main="Ross et al. ",line=1,cex.main = 1,   font.main= 3, col.main= "black")
x=DF
r <- cor.test(x[cols==2,1], x[cols==2,2])
txt1 <- paste("r=",format(c(r$estimate, 0.123456789), digits=2)[1]," P=",format(c(r$p.value, 0.123456789), digits=2)[1], sep="")
text(0, 1300, txt1, col=2,cex = .7,pos=4)
r <- cor.test(x[cols==3,1], x[cols==3,2])
txt2 <- paste("r=",format(c(r$estimate, 0.123456789), digits=2)[1]," P=",format(c(r$p.value, 0.123456789), digits=2)[1], sep="")
text(0, 1150, txt2, col=3,cex = .7,pos=4)
r <- cor.test(x[cols==4,1], x[cols==4,2])
txt3 <- paste("r=",format(c(r$estimate, 0.123456789), digits=2)[1]," P=",format(c(r$p.value, 0.123456789), digits=2)[1], sep="")
text(0, 1225, txt3, col=4,cex = .7,pos=4)
detach(DF)


load("FmorrYeohALL.Rdata") 
DF=flux[1:49,c("MTHFD","TYMS")]
attach(DF)
plot(TYMS,MTHFD/2, pch=pchs,col=cols,cex=.6,xlab="dTMP Flux (uM/hr)",ylab="DNPS Flux (uM/hr)",ylim=c(0,1300),xlim=c(0,14))
title(main="Yeoh et al. ",line=1,cex.main = 1,   font.main= 3, col.main= "black")
x=DF
r <- cor.test(x[cols==2,1], x[cols==2,2])
txt1 <- paste("r=",format(c(r$estimate, 0.123456789), digits=2)[1]," P=",format(c(r$p.value, 0.123456789), digits=2)[1], sep="")
text(0, 1300, txt1, col=2,cex = .7,pos=4)
r <- cor.test(x[cols==3,1], x[cols==3,2])
txt2 <- paste("r=",format(c(r$estimate, 0.123456789), digits=2)[1]," P=",format(c(r$p.value, 0.123456789), digits=2)[1], sep="")
text(0, 1150, txt2, col=3,cex = .7,pos=4)
r <- cor.test(x[cols==4,1], x[cols==4,2])
txt3 <- paste("r=",format(c(r$estimate, 0.123456789), digits=2)[1]," P=",format(c(r$p.value, 0.123456789), digits=2)[1], sep="")
text(0, 1225, txt3, col=4,cex = .7,pos=4)
detach(DF)


par(mfrow=c(1,1))
dev.copy(pdf,file="fig5abcd.pdf", width = 7, height = 7)

dev.off()


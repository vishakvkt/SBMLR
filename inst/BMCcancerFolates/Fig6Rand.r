setwd(file.path(.path.package("SBMLR"), "BMCcancerFolates")) #default dump site 
#setwd("C:/cwru/active/Morrison")  # set this to where figs should be dumped, with comment removed
load("FmorrRand.Rdata") 
DF=flux[,c("MTHFD","TYMS")]
attach(DF)
plot(TYMS,MTHFD/2, pch=1,xlab="dTMP Flux (uM/hr)",ylab="DNPS Flux (uM/hr)",ylim=c(0,1300),xlim=c(0,14))
title(main="Random",line=1,cex.main = 1,   font.main= 3, col.main= "black")
x=DF
r <- cor.test(x[,1], x[,2])
txt1 <- paste("r=",format(c(r$estimate, 0.123456789), digits=2)[1]," P=",format(c(r$p.value, 0.123456789), digits=2)[1], sep="")
text(0, 1300, txt1, col=2,cex = .7,pos=4)
detach(DF)
dev.copy(pdf,file="fig6rand.pdf", width = 6, height = 6)
dev.off()



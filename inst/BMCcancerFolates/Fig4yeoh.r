# This is Morrison's folate model driven by Yoeh et al also in [M.E. Ross et al, Blood (2003)]
library(Biobase)
library(odesolve)
library(annotate)
library(hgu95av2)
library(SBMLR)  
setwd(file.path(.path.package("SBMLR"), "BMCcancerFolates")) #default dump site 
#setwd("C:/cwru/active/Morrison")  # set this to where figs should be dumped, with comment removed

morr=readSBMLR(file.path(.path.package("SBMLR"), "models/morrison.r"))  
library(rossEset)
pD=subset(pData(ross),subset=(type=="TEL.AML1")|(type=="T.Cell")|(type=="BCR.ABL"))
pD$type=factor(pD$type, levels=c("BCR.ABL","TEL.AML1","T.Cell"))
pD$type
as.numeric(pD$type)
I=order(pD$type)
pD=pD[I,]
summary(pD$type)

library(yeoh5Eset)
eset=yeoh5[,rownames(pD)]
eset
pD
pData(eset)

morrsym=c('MTHFD1','GART','ATIC','TYMS','DHFR')
morrsym=c('SHMT1','MTHFR','MTR','MTHFD1','GART','ATIC','TYMS','DHFR')
key=c(GARFT="GART",ATIC7="ATIC",MTHFD="MTHFD1",TYMS="TYMS",DHFReductase="DHFR",ATIC12="ATIC")
key=c(MTHFR="MTHFR",MTR="MTR",SHMT="SHMT1",SHMTr="SHMT1",GARFT="GART",ATIC7="ATIC",MTHFD="MTHFD1",TYMS="TYMS",DHFReductase="DHFR",ATIC12="ATIC")

AffyID<- ls(env = hgu95av2SYMBOL)
lsym <- mget(AffyID, env = hgu95av2SYMBOL)
sym <- as.character(lsym)
names(sym)=names(lsym)
mID=names(sort(sym[is.element(sym,morrsym)]))
msym=as.character(mget(mID, env = hgu95av2SYMBOL))
(folate <- eset[mID, ])
mykp <- function(y) kruskal.test(y, pD$type)$p.value
om <- esApply(folate[mID,], 1, mean) 
tam <- esApply(folate[mID,folate$type=="TEL.AML1" ], 1, mean) 
tm <- esApply(folate[mID,folate$type=="T.Cell" ], 1, mean)
bam <- esApply(folate[mID,folate$type=="BCR.ABL" ], 1, mean)
pvals <- esApply(folate, 1, mykp)
outdat=data.frame(ID=mID,sym=msym,tam,bam,tm,om,pvals)#[om>500,]#[pvals<.05,] 
outdat  # dump out the table

counts=summary(outdat$sym)
noutdat=data.frame(NULL)
k=1
for (j in 1:length(counts)) 
   {
   kw=outdat$om[k:(k+counts[j]-1)]
   noutdat=rbind(noutdat,outdat[k+which(kw==max(kw))-1, ]    )
#   kw=outdat$pvals[k:(k+counts[j]-1)]
#   noutdat=rbind(noutdat,outdat[k+which(kw==min(kw))-1, ]    )
k=k+counts[j]
   }
k
row.names(noutdat)<-noutdat$sym
noutdat=noutdat[morrsym,]


write.table(noutdat,"table1Ross.csv",sep=",",row.names=F)
sID=as.character(noutdat$ID)
sfolate <- folate[sID, ]  # no doubles in sfolate 

gpeset=NULL
gpeset=exprs(sfolate)
row.names(gpeset)<-morrsym
DF=cbind(as.data.frame(t(gpeset)))


panel.hist <- function(x, ...)   {usr <- par("usr"); on.exit(par(usr))
 par(usr = c(usr[1:2], 0, 3) )
         h <- hist(x[cols==2], plot = FALSE)  # diagnosis
         breaks <- h$breaks; nB <- length(breaks)
         y <- h$counts; y <- y/max(y)
         lines(breaks[-nB], y, col=2, ...) 
         h <- hist(x[cols==3], plot = FALSE)  # diagnosis
         breaks <- h$breaks; nB <- length(breaks)
         y <- h$counts; y <- y/max(y)
         lines(breaks[-nB], y, col=3, ...) 
         h <- hist(x[cols==4], plot = FALSE)  # post treatment
         breaks <- h$breaks; nB <- length(breaks)
         y <- h$counts; y <- y/max(y)
         lines(breaks[-nB], y, col=4, ...) 
   	 par(usr = c(0,1, 0, 1.5) )
 	 wp=kruskal.test(x,g=pD$type)$p.value
         text(0.0,.93,"Kruskal-Wallis Test", font=3,col=1,cex = .7,pos=4)
  	 txt <- paste("P=",format(c(wp, 0.123456789), digits=4)[1], sep="")
          text(0.0,0.70, txt, col=1,cex = .7,pos=4)
             } 
panel.cor <- function(x, y,...)
     {
         usr <- par("usr"); on.exit(par(usr))
         par(usr = c(0, 1, 0, 1))
         r <- cor.test(x[cols==2], y[cols==2])
         txt <- paste("r=",format(c(r$estimate, 0.123456789), digits=2)[1]," P=",format(c(r$p.value, 0.123456789), digits=2)[1], sep="")
         text(0.0, 0.70, txt, col=2,cex = .7,pos=4)
         r <- cor.test(x[cols==3], y[cols==3])
         txt <- paste("r=",format(c(r$estimate, 0.123456789), digits=2)[1]," P=",format(c(r$p.value, 0.123456789), digits=2)[1], sep="")
         text(0.0, 0.30, txt, col=3,cex = .7,pos=4)
         r <- cor.test(x[cols==4], y[cols==4])
         txt <- paste("r=",format(c(r$estimate, 0.123456789), digits=2)[1]," P=",format(c(r$p.value, 0.123456789), digits=2)[1], sep="")
         text(0.0, 0.50, txt, col=4,cex = .7,pos=4)
     }

cols=1+c(1,3,2)[as.numeric(pD$type)]
pchs=c("b","B","T")[as.numeric(pD$type)]


pairs(DF, pch=pchs,col=cols,cex=.6,diag.panel=panel.hist,lower.panel=panel.cor,gap=0)
title(main="Yeoh et al: b=BCR-ABL B-cell (red); B=TEL-AML1 B-cell (blue);  T=T-cell (green)",line=3,cex.main = .7,   font.main= 3, col.main= "black")

dev.copy(pdf,file="fig4yeoh.pdf", width = 8, height = 8)
dev.off() # close the file device just opened, i.e. the dev.cur()


save(DF,file="yeohExprs.Rdata")


na=data.frame(exprs(sfolate));na
cnt=apply(cbind(apply(na[,1:15],1,median),apply(na[,16:35],1,median),apply(na[,36:49],1,median)),1,mean)#;cnt #this is the baseline for all perturbations
aa=na/(cnt%o%rep(1,dim(na)[1]));aa
aa=cbind(aa,control=rep(1,dim(gpeset)[1]))

rownames(aa)=morrsym

mi=summary(morr)
attach(mi)  # this gives rIDs

M=matrix(rep(1,dim(aa)[2]*length(rIDs)),nrow=length(rIDs))
rownames(M)<-rIDs
colnames(M)<-colnames(aa)
M[names(key),]
tmp=as.matrix(aa[key,])
rownames(tmp)<-names(key)
M[names(key),]=tmp
M

morr$species$EMTX$ic=0
out1=simulate(morr,seq(-20,0,1))
morr$species$EMTX$ic=1
out2=simulate(morr,0:30)
outs=data.frame(rbind(out1,out2))
morr$species$EMTX$ic=0

attach(outs)

par(mfrow=c(3,4))
plot(time,FH2b,type="l",xlab="Hours")
plot(time,FH2f,type="l",xlab="Hours")
plot(time,DHFRf,type="l",xlab="Hours")
plot(time,DHFRtot,type="l",xlab="Hours")
plot(time,CHOFH4,type="l",xlab="Hours")
plot(time,FH4,type="l",xlab="Hours")
plot(time,CH2FH4,type="l",xlab="Hours")
plot(time,CH3FH4,type="l",xlab="Hours")
plot(time,AICARsyn,type="l",xlab="Hours")
plot(time,MTR,type="l",xlab="Hours")
plot(time,TYMS,type="l",xlab="Hours")
#plot(time,EMTX,type="l",xlab="Hours")
plot(time,DHFReductase,type="l",xlab="Hours")
par(mfrow=c(1,1))
detach(outs)


nFluxes=length(rIDs)
# now make the big flux matrix. This takes time to run!!!!!
flux=matrix(rep(0,nFluxes*50),ncol=nFluxes,nrow=50)
conc=matrix(rep(0,nStates*50),ncol=nStates,nrow=50)
rownames(flux)<-c(paste("b",1:15,sep=""),paste("B",1:20,sep=""),paste("T",1:14,sep=""),"Control")
colnames(flux)<-rIDs
rownames(conc)<-c(paste("b",1:15,sep=""),paste("B",1:20,sep=""),paste("T",1:14,sep=""),"Control")
colnames(conc)<-names(y0)
flux
conc

for (patient in 1:50)
{
print(patient)
out1=simulate(morr,seq(-20,0,1),M[,patient])
out2=simulate(morr,0:30,M[,patient])
outs=data.frame(rbind(out1,out2))

conc[patient,]=as.numeric(outs[dim(outs)[1],2:(nStates+1)])
flux[patient,]=as.numeric(outs[dim(outs)[1],(nStates+2):(nStates+nFluxes+1)])

par(mfrow=c(3,1))
plot(FH2f~time,data=outs)
title(main=paste("patient ",patient))
plot(FH4~time,data=outs)
plot(CH2FH4~time,data=outs)
par(mfrow=c(1,1))
}

flux=data.frame(flux)
conc=data.frame(conc)
detach(mi)
save(flux,conc,file="FmorrYeohALL.Rdata")# save flux array since it takes much time to recompute
#  END big computation loop

#  Now do plotting and stats for the predicted fluxes
# load("FmorrRossBT.Rdata") # uncomment this if you saved the flux array 4 lines up 
attach(flux)
attach(conc)
#cols=1+as.numeric(pD$type)
#pchs=c("b","T","B")[as.numeric(pD$type)]
#DF=flux[1:49,c("MTHFR","MTHFD","TYMS")]
#DF
#
#pairs(DF, pch=pchs,col=cols,cex=.6,diag.panel=panel.hist,lower.panel=panel.cor,gap=0)
#title(main="Ross et al: b=BCR-ABL B-cell (red); B=TEL-AML1 B-cell (blue);  T=T-cell (green)",line=3,cex.main = .7,   font.main= 3, col.main= "black")
plot(TYMS,MTHFD/2,type="n",xlim=range(TYMS),ylim=range(MTHFD/2),xlab="dTMP Flux (uM/hr)",ylab="DNPS Flux (uM/hr)")
points(TYMS[folate$type=="TEL.AML1"],MTHFD[folate$type=="TEL.AML1"]/2,pch=7)
points(TYMS[folate$type=="T.Cell"],MTHFD[folate$type=="T.Cell"]/2,pch=2)
points(TYMS[folate$type=="BCR.ABL"],MTHFD[folate$type=="BCR.ABL"]/2,pch=1)
legend(locator(n=1),legend=c("BCR-ABL B-cell","TEL-AML1 B-cell","T-cell"), pch=c(1,7,2))
cor.test(TYMS[folate$type=="TEL.AML1"],MTHFD[folate$type=="TEL.AML1"])
cor.test(TYMS[folate$type=="T.Cell"],MTHFD[folate$type=="T.Cell"])
t.test(TYMS[folate$type=="TEL.AML1"],TYMS[folate$type=="T.Cell"])
t.test(MTHFD[folate$type=="TEL.AML1"],MTHFD[folate$type=="T.Cell"])
detach(flux)
detach(conc)

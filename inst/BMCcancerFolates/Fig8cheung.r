# This is Morrison's folate model driven by Vivian Cheung et al's Radiation time course response data (GDS479).
library(Biobase)
library(odesolve)
library(annotate)
library(hgu95av2)
library(SBMLR)  
setwd(file.path(.path.package("SBMLR"), "BMCcancerFolates")) #default dump site 
#setwd("C:/cwru/active/Morrison")  # set this to where figs should be dumped, with comment removed
morr=readSBMLR(file.path(.path.package("SBMLR"), "models/morrison.r"))  
library(cheungEset)
eset=cheung
pD=pData(eset);pD
morrsym=c('SHMT1','MTHFR','MTR','MTHFD1','GART','ATIC','TYMS','DHFR')
#morrsym=c('MTR','MTHFD1','GART','ATIC','TYMS','DHFR')
#key=c(MTHFR="MTHFR",MTR="MTR",SHMT="SHMT1",SHMTr="SHMT1",GARFT="GART",ATIC7="ATIC",MTHFD="MTHFD1",TYMS="TYMS",DHFReductase="DHFR",ATIC12="ATIC")
key=c(MTR="MTR",GARFT="GART",ATIC7="ATIC",MTHFD="MTHFD1",TYMS="TYMS",DHFReductase="DHFR",ATIC12="ATIC")

AffyID<- ls(env = hgu95av2SYMBOL)
lsym <- mget(AffyID, env = hgu95av2SYMBOL)
sym <- as.character(lsym)
names(sym)=names(lsym)
mID=names(sort(sym[is.element(sym,morrsym)]))
folate <- eset[mID, ]
msym=as.character(mget(mID, env = hgu95av2SYMBOL))
om <- esApply(folate[mID,], 1, mean) 
outdat=data.frame(ID=mID,sym=msym,om)#[om>500,] 
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
noutdat

row.names(noutdat)<-noutdat$sym
noutdat=noutdat[morrsym,]


noutdat
write.table(noutdat,"table1Cheung.csv",sep=",",row.names=F)
sID=as.character(noutdat$ID)
sfolate <- folate[sID, ]  # no doubles in sfolate 
na=data.frame(exprs(sfolate));na
#na[na<200]=200 # set a floor to avoid negative values
cnt=apply(na[,pD$dose==0],1,mean)#;cnt #this is the baseline for all perturbations
aa=na/(cnt%o%rep(1,dim(na)[2]));aa
aa=cbind(aa,ctrl=rep(1,dim(na)[1]))
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

M3= cbind(M[,"ctrl"],M[,pD$dose==3] )
M10=cbind(M[,"ctrl"],M[,pD$dose==10])
times=c(0,1,2,6,12,24)
finet=seq(0,24,.1)

# make up a dummy list of the right length
mods3=list(NULL)
mods10=list(NULL)
for (i in 1:length(rIDs))
{
mods3[[i]]=1
mods10[[i]]=1
}
names(mods3)<-rIDs
names(mods10)<-rIDs

for (i in 1:length(rIDs))
{
mods3[[i]]=approxfun(times,M3[rIDs[i],],method="linear",rule=2)
mods10[[i]]=approxfun(times,M10[rIDs[i],],method="linear",rule=2)
}

Mt=mods3
out1=simulate(morr,seq(-10,0,1),Mt)
out2=simulate(morr,0:30,Mt)
gray3=data.frame(rbind(out1,out2))
Mt=mods10  
out1=simulate(morr,seq(-10,0,1),Mt)
out2=simulate(morr,0:30,Mt)
gray10=data.frame(rbind(out1,out2))

graphics.off()  # clear off current figures
toPlot=c("TYMS","GARFT")
par(mfrow=c(1,2),mex=.7)
for (i in 1:length(toPlot))
{
plot(gray3[,"time"],gray3[,toPlot[i]],type="l",lty=1,ylab=toPlot[i],xlab="Hours",ylim=range(c(gray3[,toPlot[i]],gray10[,toPlot[i]])) ) #,xlim=c(0,24))
lines(gray10[,"time"],gray10[,toPlot[i]],lty=2)
if (i==2) legend(12,650,legend=c("3 gray","10 gray"),lty=c(1,2))
}

dev.copy(pdf,file="fig8cheungOut.pdf", width = 7.5, height = 7.5)
dev.off() # close the file device just opened, i.e. the dev.cur()

windows(xpos=-220,ypos=-50)  # place interpolation plots in a separate figure window down and to the left
toPlot=c("MTHFR","MTR","SHMT","MTHFD","GARFT","ATIC7","TYMS","DHFReductase")

par(mfrow=c(3,3),mex=.7)
for (i in 1:length(toPlot))
{
plot(times,M3[toPlot[i],],type="p",pch=2,ylab=toPlot[i],xlim=c(0,24),ylim=c(0,2))
lines(finet,mods3[[toPlot[i]]](finet))
points(times,M10[toPlot[i],],type="p",pch=3)
lines(finet,mods10[[toPlot[i]]](finet),lty=2)
if (i==1) legend(5,2,legend=c("3 gray","10 gray"),pch=c(2,3))
}
par(mfrow=c(1,1))
detach(mi)  

# use the other cheung script to get figure 7 with all of the genes!

#dev.copy(pdf,file="curCheungIn.pdf", width = 7.5, height = 7.5)  
#dev.off() # close the file device just opened, i.e. the dev.cur()




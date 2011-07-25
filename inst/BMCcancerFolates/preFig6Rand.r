# This is Morrison's folate model driven by random data
library(Biobase)
library(odesolve)
library(annotate)
library(hgu95av2)
library(SBMLR)  
setwd(file.path(.path.package("SBMLR"), "BMCcancerFolates")) #default dump site 
#setwd("C:/cwru/active/Morrison")  # set this to where figs should be dumped, with comment removed

morr=readSBMLR(file.path(.path.package("SBMLR"), "models/morrison.r"))  

morrsym=c('MTHFD1','GART','ATIC','TYMS','DHFR')
morrsym=c('SHMT1','MTHFR','MTR','MTHFD1','GART','ATIC','TYMS','DHFR')
key=c(GARFT="GART",ATIC7="ATIC",MTHFD="MTHFD1",TYMS="TYMS",DHFReductase="DHFR",ATIC12="ATIC")
key=c(MTHFR="MTHFR",MTR="MTR",SHMT="SHMT1",SHMTr="SHMT1",GARFT="GART",ATIC7="ATIC",MTHFD="MTHFD1",TYMS="TYMS",DHFReductase="DHFR",ATIC12="ATIC")

#npats=1000 # this is the real one in the paper, but it takes ~10 hours! ... I will look for faster ways shortly, perhaps via pysces. 
npats=10
aa=matrix(rnorm(npats*length(morrsym),mean=1,sd=.3),ncol=npats)
aa=cbind(aa,control=rep(1,length(morrsym)) )
rownames(aa)=morrsym
colnames(aa)<-c(paste("r",1:npats,sep=""),"control")

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


nFluxes=length(rIDs)
# now make the big flux matrix. This takes time to run!!!!!
flux=matrix(rep(0,nFluxes*(npats+1)),ncol=nFluxes,nrow=(npats+1))
conc=matrix(rep(0,nStates*(npats+1)),ncol=nStates,nrow=(npats+1))
rownames(flux)<-c(paste("r",1:npats,sep=""),"control")
colnames(flux)<-rIDs
rownames(conc)<-c(paste("r",1:npats,sep=""),"control")
colnames(conc)<-names(y0)
flux
conc

for (patient in 1:(npats+1))
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
save(flux,conc,file="FmorrRand.Rdata")# save flux array since it takes much time to recompute
#  END big computation loop

#  Now do plotting and stats for the predicted fluxes
# load("FmorrRossBT.Rdata") # uncomment this if you saved the flux array 4 lines up 
attach(flux)
attach(conc)
plot(TYMS,MTHFD/2,pch=1,xlim=range(TYMS),ylim=range(MTHFD/2),xlab="dTMP Flux (uM/hr)",ylab="DNPS Flux (uM/hr)")
detach(flux)
detach(conc)

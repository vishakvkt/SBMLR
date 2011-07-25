\name{simulate}
\alias{simulate}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Simulate a model of class SBML}
\description{
  This function simulates a model given the report times and optional modulators. It uses \code{lsoda} of the odesolve package.
}
\usage{
simulate(model, times, modulator=NULL,X0=NULL, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{model}{ The model object to be simulated. Initial conditions are passed through this object.}
  \item{times}{ The sequence of time points to be sampled and provided as rows of the output matrix.}
  \item{modulator}{ Null if there are no modulators (default), 
  a vector of numbers if there are steady state Vmax modulators,  
  and a list of interpolating functions if there are time course Vmax modulators. }
  \item{X0}{ Override model initial conditions in simulations, particularly piece-wise perturbation simulations.} 
  \item{...}{ For compatibility with \code{simulate} of the stats package.}
}
\details{
  This is a wrapper for lsoda. }
\value{ The data frame output that comes out of \code{lsoda}.}
\references{For the folate cycle example given below: Morrison PF, Allegra CJ: Folate cycle kinetics in human 
  breast cancer cells. JBiolChem 1989, 264(18):10552-10566.}
\author{ Tom Radivoyevitch}
\note{ Rules are implemented through time varying boundary conditions updated at each time point 
as a side effect within the (now internal) function \code{fderiv}. }

\examples{
##---- The following example performs a perturbation in PRPP from 5 to 50 uM in Curto et al.'s model.
library(SBMLR) 
library(odesolve)
curto=readSBML(file.path(system.file(package="SBMLR"), "models/curto.xml"))  
out1=simulate(curto,seq(-20,0,1))
curto$species$PRPP$ic=50
out2=simulate(curto,0:70)
outs=data.frame(rbind(out1,out2))
attach(outs)
par(mfrow=c(2,1))
plot(time,IMP,type="l")
plot(time,HX,type="l")
par(mfrow=c(1,1))
detach(outs)

# which should be the same plots as 
curto=readSBMLR(file.path(system.file(package="SBMLR"), "models/curto.r"))  
out1=simulate(curto,seq(-20,0,1))
curto$species$PRPP$ic=50
out2=simulate(curto,0:70)
outs=data.frame(rbind(out1,out2))
attach(outs)
par(mfrow=c(2,1))
plot(time,IMP,type="l")
plot(time,HX,type="l")
par(mfrow=c(1,1))
detach(outs)

##---- The following example uses fderiv to generate Morrison's folate system response to 1uM MTX 

morr=readSBMLR(file.path(system.file(package="SBMLR"), "models/morrison.r"))  
out1=simulate(morr,seq(-20,0,1))
morr$species$EMTX$ic=1
out2=simulate(morr,0:30)
outs=data.frame(rbind(out1,out2))
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
morr$species$EMTX$ic=0
}
\keyword{ math }% at least one, from doc/KEYWORDS

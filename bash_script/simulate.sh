#!/bin/bash
########################################
#                                      #
# SGE MPI job script for ECDF Cluster  #
#                                      #
# by ECDF System Team                  #
# ecdf-systems-team@lists.ed.ac.uk     #
#                                      #
########################################

# Grid Engine options

#$ -N SimulateModel				# jobname to check status
#$ -cwd					 	# to run in current working directory
#$ -l h_rt=00:60:00				# specifying jobtime since it is below 30 minutes
#$ -o o.dat
#$ -e e.dat
#$ -V 

# Initialise environment module
 
. /etc/profile.d/modules.sh

# Use OpenMPI and GNU compilers

module load R

N=$SGE_TASK_ID

TEMP=$R_LIBS
R_LIBS=/exports/home/s1061101/R/x86_64-redhat-linux-gnu-library/2.11

R --slave  <<EOF
#-v R_HOME = /exports/home/s1061101/R/x86_64-redhat-linux-gnu-library/2.11

num<-$N;								# number of jobs

library(XML)								# loading required packages
library(odesolve)
library(SBMLR)

load('model.dat')							# loading the model object
load('sobol.dat')							# loading the sobol matrix
source('sobol_sim.r')							# loading sobol function
mod <- model

#remain <- nrow(sets) %% num						# in case the world is going to end
#perNode <- floor(nrow(sets)/num)					# arbitrary calculation
									#first node gets a few extra sets if odd num of jobs
#if(num==1 && remain !=0) { out<- sobol_sim(mod,sets((num*perNode):(perNode+remain),])
#if(num==1 && remain ==0) { out<- sobol_sim(mod,sets(1:perNode,])}
#initial <- num * perNode
#out <- sobol_sim(mod, sets[(initial+1):((initial-1)+perNode), ])	# each node computes 10 sobol sets

out <-sobol_sim(mod, sets[(num-9):num, ])
save(file=paste(as.character(num),'.Rsim',sep=''),out)

done <- "process done"
save(file = paste(as.character(num),'.done',sep=''), done)		# written when job complete

if(FALSE) {

gcinfo(TRUE);
load("dat.Rdat");							# object deserialization and loading
rm(d_max,d_mean,d_min,d_nnint,c_max,c_mean,c_min,costs,c_sd,nnint,lf)

	datS<-dat[dat[['Total']]<1e5,];					# 1e3
	rm(dat)
gc()									# memory usage functions
	sens<-pcc(datS[1:8e4,parI:parN],datS[1:8e4,names[num]], rank = TRUE);
	str(sens)							# debug print
	prcc_l<-sens[['PRCC']];						# extract prcc
	colnames(prcc_l)<-names[num];					# set value according to names vector containing species names
	str(prcc_l)							# debug for sanity
save(file=paste(names[num],'.Rsens',sep=''),prcc_l);
q(save='no');
}

EOF
R_LIBS=$TEMP
date

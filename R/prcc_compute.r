# Wrapper function for sensitivity calculation
# dependencies - library(sensitivity)
#
# Author: vishak

"prcc_compute" <- function(SpeciesMatrix, sobolMatrix) {
	
	library(sensitivity)

	NumSets <- nrow(sobolMatrix)		#num of sobol sets
	nParams <- ncol(sobolMatrix)		#num of Parameter inputs
	nSpecies <- ncol(SpeciesMatrix)
	
	sens <- matrix(data=NA, nParams, nSpecies)			# matrix param_count x Area_curve_count which is the species count

	species_names <- colnames(SpeciesMatrix)

	frame_spec <- data.frame(SpeciesMatrix)
	frame_sobol <- data.frame(sobolMatrix)
		for(j in 1:nSpecies){

				senslist <- pcc(frame_sobol[1:NumSets,1:nParams],frame_spec[1:NumSets,species_names[j]], rank=TRUE)
#				sens<-pcc(datS[1:8e4,parI:parN],datS[1:8e4,names[num]], rank = TRUE);
				#str(senslist)
				#prcc_l<-senslist[['PRCC']];						# extract prcc
				#colnames(prcc_l)<-species_names[j];					# set value according to names vector containing species names
				#str(prcc_l)								# debug for sanity
				#print(prcc_l)
				print(senslist)
			}

		#senslist <- pcc(sobolvect, speciesvect, nboot = NumSets, rank=TRUE)
		#sens[,i] <- senslist[['PRCC']]
		#senslist <- pcc(sobolMatrix[1:NumSets, ], speciesvect, rank=TRUE)
	#}
#	rownames(sens) <- colnames(sobolMatrix)
	sens
}

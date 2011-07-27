# Wrapper function for sensitivity calculation
# dependencies - library(sensitivity)
#
# Author: vishak

"prcc_compute" <- function(SpeciesMatrix, sobolMatrix) {
	
	library(sensitivity)			#This wont execute if the library is already loaded

	NumSets <- nrow(sobolMatrix)		#num of sobol sets
	nParams <- ncol(sobolMatrix)		#num of Parameter inputs
	nSpecies <- ncol(SpeciesMatrix)
	
	sens <- list()				# matrix param_count x Area_curve_count which is the species count

	if(NumSets < nParams)			# If number of sobol sets less than parameter, cannot run sensitivity.
	{
		cat('Cannot run sensitivity when the number of sobol sets is less than the number of parameters. Simulate More points','\n')
		return()
	}
	species_names <- colnames(SpeciesMatrix)

	frame_spec <- data.frame(SpeciesMatrix)
	frame_sobol <- data.frame(sobolMatrix)
		for(j in 1:nSpecies){

				senslist <- pcc(frame_sobol[1:NumSets,1:nParams],frame_spec[1:NumSets,species_names[j]], rank=TRUE)
				prcc_l<-senslist[['PRCC']];						# extract prcc
				colnames(prcc_l) <- species_names[[j]];
				sens[[j]] <- prcc_l;
		}

	#colnames(sens) <- colnames(SpeciesMatrix)							# set matrix according to names vector containing species names
	#rownames(sens) <- colnames(sobolMatrix)								# set matrix row names as the parameter names
	sens
}

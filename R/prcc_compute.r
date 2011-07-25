# Wrapper function for sensitivity calculation
# dependencies - library(sensitivity)
#
# Author: vishak

"prcc_compute" <- function(SpeciesMatrix, sobolMatrix) {
	
	library(sensitivity)

	NumSets <- nrow(sobolMatrix)		#num of sobol sets
	
	nSpecies <- ncol(SpeciesMatrix)
	
	sens <- list()


	for(i in 1:NumSets) {
		#sens <- pcc(sobolMatrix[1:NumSets, ], SpeciesMatrix[1:NumSets, ], rank=TRUE)
		sobolvect <- sobolMatrix[i,]
		speciesvect <- SpeciesMatrix[i,]
		sens[[i]] <- pcc(sobolvect, speciesvect, rank=TRUE)
	}
	sens
}

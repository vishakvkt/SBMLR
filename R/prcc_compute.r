# Wrapper function for sensitivity calculation
# dependencies - library(sensitivity)
#
# Author: vishak

"prcc_compute" <- function(SpeciesMatrix, sobolMatrix) {
	
	library(sensitivity)

	NumSets <- nrow(sobolMatrix)
	
	nSpecies <- ncol(SpeciesMatrix)
	
	sens <- list()

	for(i in 1:nSpecies) {
		sens <- pcc(sobolMatrix[1:NumSets, ], SpeciesMatrix[1:NumSets, i], rank=TRUE)
	}
	sens
}

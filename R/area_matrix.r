
# Function to create the Species Matrix with rows equal to number of sobol sets, and columns equal to number of species
# This calculates the area under the curve by trapezoidal rule
#
# Author: vishak

"area_matrix" <- function(model, SimulatedModelSet) {

	nSets <- length(SimulatedModelSet)

	nSpecies <- length(model$species)
	SpeciesArea <- matrix(data=NA, nSets, nSpecies)
	snameVector <- NULL
	for(s in 1:nSets) {

	SimulateFrame <- data.frame(SimulatedModelSet[[s]])
	Integral <- 0
	
	time  <- length(SimulateFrame$time)
	
						# this will be n x m matrix where
						# n - Area under the curve for each sobol set
						# m - number of species Integral values
	for(i in 1:nSpecies) {

		SpeciesName <- model$species[[i]][["id"]]
		if(s==1)snameVector <- c(snameVector, SpeciesName)		#only need to do once
		nTime <- length(SimulateFrame[[SpeciesName]])			#sanity check.
	
		for(j in 2:nTime) {
			Integral <- Integral + ((SimulateFrame[[SpeciesName]][[j]] + SimulateFrame[[SpeciesName]][[j-1]]) / 2)
		}
		SpeciesArea[s, i] <- Integral
		Integral <- 0
	}
   }
	dimnames(SpeciesArea) <- list(NULL, snameVector) 		#setting the names of the columns
	SpeciesArea
}

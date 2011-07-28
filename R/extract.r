"extract" <- function(model, SimulatedModel) {

	SimulateFrame <- data.frame(SimulatedModel)
	nSpecies <- length(model$species)

	snameVector <- NULL
	time  <- length(SimulateFrame$time)
	SpeciesSim <- matrix(data=NA, time, nSpecies) # data, num_rows,num_columns

						# this will be m x n matrix where
						# n - number of species
						# m - number of timepoints
	for(i in 1:nSpecies) {
		SpeciesName <- model$species[[i]][["id"]]
		snameVector <- c(snameVector, SpeciesName)
		nTime <- length(SimulateFrame[[SpeciesName]])	#sanity check.
		cat('species name: ' , SpeciesName, '\n')
		for(j in 1:nTime) {
		SpeciesSim[j, i] <- SimulateFrame[[SpeciesName]][[j]]	
		}
	}
	dimnames(SpeciesSim) <- list(NULL, snameVector) 		#setting the names of the columns
	SpeciesSim
}

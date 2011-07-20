# Wrapper function for Reading in the sobol points and creating sobol data file
#
# Author: vishak

sobolset <- function(param_path, parameter_count=0, no_of_sets=0) {

	sets <- read.table(param_path, header=TRUE)			
	no_of_sets <- nrow(sets)
	parameter_count <- ncol(sets)

	save(file=paste(getwd(),'sobol.dat',sep='/'), sets)
}

# Wrapper function for simulating multiple sets of sobol points on a given model
#
# Author: vishak

sobol_sim <- function(model, sobolset)
{
	nParams <- length(model$globalParameters)
	results <- list()
	nSets <- nrow(sobolset)						#number of sets
	time <- seq(1,60,1)
	param_count <- 1

	for(i in 1:nSets) {
			while(param_count <= nParams) {
				#Substitute the values into the model$GlobalParameters[[param_count]]
				model$globalParameters[[param_count]] <- sobolset[i,param_count]
				param_count <- param_count + 1 	
			}
			param_count <- 1
			results[[i]] <- simulate(model, time)
	}
	results
}

# Function to create the Species Matrix with rows equal to number of sobol sets, and columns equal to number of species
# This calculates the area under the curve by trapezoidal rule
#
# Author: vishak

Areamatrix <- function(model, SimulatedModelSet) {

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

# Wrapper function for sensitivity calculation
# dependecies - library(sensitivity)
#
# Author: vishak

prcc_compute <- function(SpeciesMatrix, sobolMatrix) {

	NumSets <- nrow(sobolMatrix)
	
	nSpecies <- ncol(SpeciesMatrix)
	
	sens <- list()

	for(i in 1:nSpecies) {
		sens <- pcc(sobolMatrix[1:NumSets, ], SpeciesMatrix[1:NumSets, i], rank=TRUE)
	}
	sens
}




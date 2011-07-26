# Complete wrapper for the entire subset of functions
#
# Author: vishak

#if(!("package:SBMLR" %in% search())) library(SBMLR) 				# loading SBMLR package if it hasn't been already

#if(exists('sobolset')) { rm(sobolset) } 
#source('sobolset.r', local = FALSE);						# load matrix creation function. Remove Existing copies in memory/workspace
#if(exists('sobol_sim')) { rm(sobol_sim) }	
#source('sobol_sim.r', local = FALSE);						# load sobol simulator function. Remove Existing copies in memory/workspace

#Arguments
# path_to_sobol - path to the file containing the sobol sets.
# model_path - path to the model SBML file (.xml extension)
# time - vector depicting number of minutes to run the simulation eq. seq(1,60,1)

"load_simulate" <- function(path_to_sobol, model_path, time) {

	model <- readSBML(model_path)						#Create the model object
	print('model parsed successfully')

	if(path_to_sobol == 1)	{						#No set entered. Default model simulation requested.
	
		result <- simulate(model, time)
		result	
	}
	else {
	simulated_results <-list()
	sets <- sobolset(path_to_sobol)						#Create the sobol set matrix
	print('Sobol matrix created successfully')

	simulated_results <- sobol_sim(model, sets, time)			#Simulate the model for n sets
	extracted_results <- list()
		for(i in 1:nrow(sets)) {
			extracted_results[[i]] <- extract(model, simulated_results[[i]])
		}
	print('Simulation Completed successfully')
	extracted_results							#Return the simulation results
	}
}


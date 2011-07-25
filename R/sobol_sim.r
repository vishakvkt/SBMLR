# Wrapper function for simulating multiple sets of sobol points on a given model
#
# Author: vishak

# arguments - R model object
# sobolset - matrix containing sobol points. Each row is one set
# time - time series for simulation

"sobol_sim" <- function(model, sobolset, time)
{
	nParams <- length(model$globalParameters)
	results <- list()
	nSets <- nrow(sobolset)						#number of sets
	param_count <- 1
		cat("Simulation Begins", "\n")
		for(i in 1:nSets) {

			#model$globalParameters <- sobolset[i,]
			while(param_count <= nParams) {
				#Substitute the values into the models GlobalParameters
				model$globalParameters[[param_count]] <- sobolset[i,param_count]
				param_count <- param_count + 1 	
			}
			#DEBUG ---print(model$globalParameters)			#Currently substituted set
			param_count <- 1
			results[[i]] <- simulate(model, time)
			cat("Set ", i, "done.....", "\n")
	}
	results
}

# Wrapper function for Reading in the sobol points and creating sobol data file
#
# Author: vishak

#Arguments 
# param_path - path on disk containing the sobol sets
# parameter_count - Optional parameter if only a few parameters are being changed
# no_of_sets - Optional parameter if you want to only simulate the first 10 sets in a file containing 100 sets

"sobolset" <- function(param_path, parameter_count=0, no_of_sets=0) {


	sets <- read.table(param_path, header=TRUE)			

		if(no_of_sets==0) no_of_sets <- nrow(sets)	
	if(parameter_count==0) parameter_count <- ncol(sets)

	save(file=paste(getwd(),'sobol.dat',sep='/'), sets)
	sets
}

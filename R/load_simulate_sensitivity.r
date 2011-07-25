"load_simulate_sensitivity" <- function(path_to_sobol, model_path, time_series) {

	simulated_results <-list()
	sets <- sobolset(path_to_sobol)					#Create the sobol set matrix
	print('Sobol matrix created successfully')
	model <- readSBML(model_path)					#Create the model object
	print('model parsed successfully')
		
	simulated_results <- sobol_sim(model, sets, time_series)	#Simulate the model for n sets
	print('Simulation Completed successfully')

	species_area <- area_matrix(model, simulated_results)		#Use simulation results and calculate area under the curve for 										#each species.

	print('Area under the curve calculated successfully')

	sensitivity <- prcc_compute(species_area, sets)			#Calculate the partial rank correlation coefficient for all species
	print('Sensitivity Calculation successful.')

	sensitivity							# return prcc matrix
}

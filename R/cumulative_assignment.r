"cumulative_assignment" <- function(model, ExtractedModel)
{
		outs <- data.frame(ExtractedModel)
	
		nRules <- length(model$rules)
		nTime <- length(outs$time)							#getting length from time series. This is clean
		
		count <- 1
		n <- 1
		rules_values <-list(length(nRules))
		val <-NULL
			while(count <= nRules) 		{
				values <- model$rules[[count]]$inputs				#objects used in the law equation
			nInputs <- length(values)
			while(n <= nTime) {
				for(i in 1:nInputs){
				val <- c(val, outs[[values[[i]]]][[n]])				#Retriving values of each object for that time point
				}
			names(val) <- values							#Book-keeping to build object->value associative array
			concentration <- model$rules[[count]]$law(val)				#Evaluate expression
			rules_values[[model$rules[[count]]$idOutput]][[n]] <- concentration	#store value in list
			val <- NULL								#reset val expression
			n <- n + 1
			}
		n <- 1
		count <- count + 1
		}
		rules_values							#return list of objects containing the assignment object values over time.
}






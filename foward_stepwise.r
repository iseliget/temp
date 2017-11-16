library(dplyr)

stepwise_forward_BIC = function(data_object, response, included_covar) {
	all_covar = setdiff(colnames(data_object), response)
	included_covar = included_covar
	candidate_covar = setdiff(all_covar, included_covar)

	reg_formula = paste(response, "~", paste(included_covar, collapse=" + "))
	base_model_obj = do.call("lm", list(as.formula(reg_formula), data=as.name("data_object")))
	base_BIC = BIC(base_model_obj)
	print(paste("Base model has BIC ", base_BIC))

	for (i in 1:length(candidate_covar)) {
		print(paste("===================", "ITERATION ", i, " ========================="))
		print(paste("Iteration ", i))

		improvement = c()
		for (covar in candidate_covar) {
			current_covar = c(included_covar,covar)

			reg_formula = paste(response, "~", paste(current_covar, collapse=" + "))
			model_obj = do.call("lm", list(as.formula(reg_formula), data=as.name("data_object")))
			current_BIC = BIC(model_obj)

			improvement = c(improvement, base_BIC-current_BIC)
		}

		if (max(improvement) <= 0) {
			print("Iteration stops because no improvement in BIC can be made")
			break
		} else {
			names(improvement) = candidate_covar
			cand = names(improvement)[which.max(improvement)]
			candidate_covar = setdiff(candidate_covar,cand)
			included_covar = c(included_covar,cand)

			# update base BIC
			reg_formula = paste(response, "~", paste(included_covar, collapse=" + "))
			model_obj = do.call("lm", list(as.formula(reg_formula), data=as.name("data_object")))
			base_BIC = BIC(model_obj)
			print(paste(cand, " is included. BIC of new model is ",base_BIC ))
			print(summary(model_obj))

		}
	}
}

parkinsons = read.csv("parkinsons.csv", header=TRUE, sep=",")
parkinsons = dplyr::select(parkinsons, -motor_UPDRS)

stepwise_forward_BIC(data_object=parkinsons, 
	response="total_UPDRS", 
	included_covar=c("sex","age","test_time"))

Sys.setenv(LANG = "en")
library(dplyr)

stepwise_forward = function(data_object, response_var, included_covar) {
	# selection based on p-value will stop if p-value is not significant

	all_covariates = setdiff(colnames(data_object), response_var)
	included_covar = included_covar
	candidate_covar = setdiff(all_covariates, included_covar)

	for (i in 1:20) {
		smallest_p = 1

		for (covar in candidate_covar) {
			current_covar = c(included_covar,covar)
			rg_formula = paste(response_var, "~", paste(current_covar, collapse=" + "))
			model_obj = do.call("lm", list(as.formula(rg_formula), data=as.name("data_object")))
			coef = summary(model_obj)$coefficients
			current_p = coef[covar,4]

			if (current_p <= smallest_p) {
				smallest_p = current_p
				cand = covar
			}
		}

		if (smallest_p  > 0.05) {
			print(paste("FINAL MODEL: contains ", paste(included_covar, collapse=", ")))
			print("FINAL MODE SUMMARY: ")

			formula = paste(response_var, "~", paste(included_covar, collapse=" + "))
			final_model_obj = do.call("lm", list(as.formula(formula), data=as.name("data_object")))
			print(summary(final_model_obj))
			break		
		} else {
			included_covar = c(included_covar, cand)
			print(paste("STEP ",i, ": model contains ", paste(included_covar, collapse=", ")))
			candidate_covar = setdiff(candidate_covar, cand)
		}
	}

}

parkinsons = read.csv("parkinsons.csv", header=TRUE, sep=",")
parkinsons = dplyr::select(parkinsons, -motor_UPDRS)
stepwise_forward(data_object=parkinsons, response_var="total_UPDRS", included_covar=c("sex","age","test_time"))

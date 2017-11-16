library(dplyr)
set.seed(32) # just in case

parkinsons = read.csv("parkinsons.csv", header=TRUE, sep=",")

parkinsons = dplyr::select(parkinsons, -motor_UPDRS)
parkinsons_wo_sub = dplyr::select(parkinsons, -subject.) # "subject." should not be covariate
parkinsons_test = dplyr::select(parkinsons_wo_sub, -total_UPDRS) # data set without response

subject_index = unique(parkinsons$subject.)
length(subject_index)
# 42 unique subjects
# for 10-fold CV, each fold contains roughly 4 subjects

# generate fold list
fold = vector("list", 10)
for (i in 1:9) {
	t = sample(subject_index,4)
	fold[[i]] = t
	subject_index = setdiff(subject_index,t)
}
fold[[10]] = sample(subject_index,6)

# generate vec of covariate names
all_covar = setdiff(colnames(parkinsons_wo_sub), "total_UPDRS")
included_covar = c("sex","age","test_time")
candidate_covar = setdiff(all_covar, included_covar)

base_cv_error = 1e+60

for (i in 1:20) {
	print(paste("=================ITERATION ", i, "=================="))
	MSE_vec = c()

	for (covar in candidate_covar) {
		# set up formula for regression
		regression_formula = paste("total_UPDRS", "~", paste(c(included_covar,covar), collapse=" + "))
		print(paste("testing regression: ", regression_formula))

		MSE = 0
		for (j in 1:10) {
			foo = which(parkinsons$subject. %in% fold[[j]])
			model_obj = do.call("lm", list(as.formula(regression_formula), data=parkinsons_wo_sub[-foo,]))
			predicted_values = predict(model_obj,parkinsons_test[foo,])
			residual = predicted_values - parkinsons_wo_sub[foo, "total_UPDRS"]
			MSE = MSE + sum(residual^2)/length(residual)
		}

		print(paste("CV MSE is ", MSE/10))
		MSE_vec = c(MSE_vec, MSE/10)
	}

	if (base_cv_error - min(MSE_vec) <= 0) {
		print("Saw no improvement during this iteration. Iteration aborted.")
		break
	} else {
		base_cv_error = min(MSE_vec)
		names(MSE_vec) = candidate_covar
		cand = names(MSE_vec)[which.min(MSE_vec)]
		candidate_covar = setdiff(candidate_covar,cand) # remove chosen covariate from candidate list
		included_covar = c(included_covar, cand) # put chosen covariate into included list
		
		print(paste("Covariate ", cand, " was chosen. Cv MSE is ", min(MSE_vec)))
		print(paste("Best model so far: ", paste("total_UPDRS", "~", paste(included_covar, collapse=" + "))))

	}
}

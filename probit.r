library(MASS)
attach(birthwt)

# data preprocessings
birthwt$race <- factor(birthwt$race, label=c("white", "black", "other"))
birthwt$ptd <- factor(birthwt$ptl>0)
birthwt$ftv <- factor(birthwt$ftv)

levels(birthwt$ftv)[-(1:2)] <- "2+"

# setting up design matrix and response matrix
X<-model.matrix(~age + lwt + race + smoke + ptd + ht + ui + ftv, data=birthwt[c("age","lwt","race","smoke","ptd","ht","ui","ftv")])
Y<-as.matrix(birthwt$low)

# user-defined functions
eta = function(x, coeff) {
	return (x %*% coeff)
}

mu = function(x, link_func=default_link) {
	return (default_link(x))
}

default_link = function(x) {
	return (pnorm(x))
}

dnorm_deriv = function(x) {
	result = (-x)/(sqrt(2*pi)) * dnorm(x)
	return (result)
}

gradient = function(X,Y,coeff) {
	result = matrix(0, nrow=ncol(X), ncol=1)

	for (i in 1:nrow(X)) {
		first_const = as.numeric(Y[i,1]) - mu(eta(X[i,,drop=FALSE], coeff))
		second_const = 1/(mu(eta(X[i,,drop=FALSE], coeff))*(1 - mu(eta(X[i,,drop=FALSE], coeff))))
		third_const = dnorm(eta(X[i,,drop=FALSE], coeff))[1,1]

		result = result + (as.numeric(first_const*second_const*third_const))*t(X[i,,drop=FALSE])
	}

	return (result)
}

Hessian_first = function(X,Y,coeff) {
	result = matrix(0, nrow=ncol(X), ncol=ncol(X))

	for (i in 1:nrow(X)) {
		eta = eta(X[i,,drop=FALSE], coeff)

		const_term_num = (-1)*((dnorm(eta))^2)
		const_term_don = pnorm(eta)*(1-pnorm(eta))
		const_term = as.numeric(const_term_num/const_term_don)

		result = result + const_term * (t(X[i,,drop=FALSE]) %*% X[i,,drop=FALSE])
	}

	return (result)
}

Hessian_second = function(X,Y,coeff) {
	result = matrix(0, nrow=ncol(X), ncol=ncol(X))

	for (i in 1:nrow(X)) {
		eta = eta(X[i,,drop=FALSE], coeff)

		first_const = as.numeric(((Y[i,1]-pnorm(eta))*dnorm_deriv(eta))/(pnorm(eta)*(1-pnorm(eta))))
		second_const = as.numeric(((Y[i,1]-pnorm(eta))*(2*mu(eta)*dnorm(eta)^2 - dnorm(eta)^2))/(pnorm(eta)^2 * (1-pnorm(eta))^2))
		
		result = result + (first_const + second_const) * (t(X[i,,drop=FALSE]) %*% X[i,,drop=FALSE])
	}

	return (result)
}

log_likelihood = function(X,Y,coeff) {
	result = 0

	for (i in 1:nrow(X)) {
		first_const = as.numeric(Y[i,1,drop=FALSE]) * log(mu(eta(X[i,,drop=FALSE], coeff))/(1 - mu(eta(X[i,,drop=FALSE], coeff))))
		second_const = log(1 - mu(eta(X[i,,drop=FALSE], coeff)))

		result = result + first_const + second_const
	}

	return (result)
}

avg_log_likelihood = function(X,Y,coeff) {
	return (log_likelihood(X,Y,coeff)/nrow(X))
}

glm_probit = function(X,Y,max_iter=10000,init="zero",stop_cond="norm",fisher=TRUE) {
	# max_iter:		maximum number of iterations
	# stop_cond:	stopping condition. "norm" stands for difference in norm of coefficient.
	#				"likelihood" stands for difference in average log-likelihood
	# fisher:		if set TRUE, use Fisher scoring. if set FALSE, use regular Newton-Raphson

	coeff = matrix(0, nrow=ncol(X), ncol=1)
	if (init == "zero") {
		print("Initiall coefficients set to zeros")
	} else if (init == "ols") {
		ols_model = lm(Y~X+0)
		coeff = as.matrix(ols_model$coefficients)
		print("Initiall coefficients set to be OLS estimate")
	}
	average_log_likelihood = avg_log_likelihood(X,Y,coeff)
	counter = 0

	print(paste("Number of observations:", nrow(X)))
	print(paste("Number of covariates:", ncol(X)))

	if (fisher == TRUE) {
		print("Using Fisher scoring")
	} else {
		print("Using regular Newton-Raphson")
	}

	if (stop_cond == "norm") {
		print("Stopping condition: difference in consecutive norm of coefficient")
	} else if (stop_cond == "likelihood") {
		print("Stopping condition: difference in consecutive average log-likelihood")
	}

	for (i in 1:max_iter) {

		if (fisher == TRUE) {
			for (lambda in seq(1,0.1,-0.1)) {
				current_coeff = coeff - lambda * (solve(Hessian_first(X,Y,coeff)) %*% gradient(X,Y,coeff))
				if (avg_log_likelihood(X,Y,current_coeff) - avg_log_likelihood(X,Y,coeff) > 0) {
					break
				} else {
					print("back-track!")
				}

			}
		} else if (fisher == FALSE) {
			current_coeff = coeff - solve(Hessian_first(X,Y,coeff)+Hessian_second(X,Y,coeff)) %*% gradient(X,Y,coeff)
		}

		if ((stop_cond == "norm") & (sqrt(rowSums(t(current_coeff-coeff) %*% (current_coeff-coeff))) <= 0.00001)) {
			break
		} else if ((stop_cond == "likelihood") & (abs(avg_log_likelihood(X,Y,current_coeff) - average_log_likelihood) <= 0.00001)) {
			break
		} else {
			coeff = current_coeff
			average_log_likelihood = avg_log_likelihood(X,Y,coeff)
		}

		counter = counter + 1
	}

	print(paste("Total number of iterations:", counter))
	print("Point Estimation:")
	print(coeff)

	# calculate residual deviance
	residual_deviance = 0
	for (i in 1:nrow(X)) {
		residual_deviance = residual_deviance + (Y[i,1] * log(mu(eta(X[i,,drop=FALSE], coeff))/(1 -mu(eta(X[i,,drop=FALSE], coeff)))) + log(1-mu(eta(X[i,,drop=FALSE], coeff))))
	}
	residual_deviance = (-2)*residual_deviance
	print(paste("Residual Deviance:", residual_deviance, sep=""))

	# calculate null deviance
	deviance_of_null_model = 0
	for (i in 1:nrow(X)) {
		deviance_of_null_model = deviance_of_null_model + (Y[i] * log(mean(Y)/(1-mean(Y))) - log(1 - mean(Y)))
	}

	deviance_of_current_model = 0
	for (i in 1:nrow(X)) {
		deviance_of_current_model = deviance_of_current_model + (Y[i,1] * log(mu(eta(X[i,,drop=FALSE], coeff))/(1 -mu(eta(X[i,,drop=FALSE], coeff)))) + log(1-mu(eta(X[i,,drop=FALSE], coeff))))
	}

	null_deviance = 2*(deviance_of_null_model - deviance_of_current_model)
	print(paste("Null Deviance:", null_deviance))

	print("Standard Errors:")
	print(sqrt(diag((-1)*solve(Hessian_first(X,Y,current_coeff)+Hessian_second(X,Y,current_coeff)))))

	# return estimated coefficients
	return (current_coeff)
}


glm_probit(X,Y,init="ols",stop_cond="norm",fisher=FALSE)
glm_probit(X,Y,init="ols",stop_cond="norm",fisher=TRUE)
glm_probit(X,Y,init="ols",stop_cond="likelihood",fisher=FALSE)
glm_probit(X,Y,init="ols",stop_cond="likelihood",fisher=TRUE)
glm_probit(X,Y,init="zero",stop_cond="norm",fisher=FALSE)
glm_probit(X,Y,init="zero",stop_cond="norm",fisher=TRUE)
glm_probit(X,Y,init="zero",stop_cond="likelihood",fisher=FALSE)
glm_probit(X,Y,init="zero",stop_cond="likelihood",fisher=TRUE)

rm(list = ls())

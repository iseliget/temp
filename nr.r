gradient = function(X,Y,coef) {
	result = matrix(0, nrow=ncol(X), ncol=1)

	for (i in 1:nrow(X)) {

		first_const = as.numeric(Y[i,1]) - mu(X[i,,drop=FALSE], coef)
		second_const = 1/(mu(X[i,,drop=FALSE], coef)*(1 - mu(X[i,,drop=FALSE], coef)))
		third_const = pnorm(X[i,,drop=FALSE] %*% coef)[1,1]

		result = result + (first_const*second_const*third_const)*t(X[i,,drop=FALSE])

	}

	return (result)
}

Hessian_first = function(X,Y,coeff) {
	result = matrix(0, nrow=ncol(X), ncol=ncol(X))

	for (i in 1:nrow(X)) {
		eta = X[i,,drop=FALSE] %*% coeff

		const_term_num = (-1)*((dnorm(eta))^2)
		const_term_don = pnorm(eta)*(1-pnorm(eta))
		const_term = as.numeric(const_term_num/const_term_don)

		result = result + const_term * (t(X[i,,drop=FALSE]) %*% X[i,,drop=FALSE])
	}

	print("Hessian matrix:")
	print(result)
	return (result)
}

log_likelihood = function(X,Y,coeff) {
	result = 0

	for (i in 1:nrow(X)) {
		first_const = as.numeric(Y[i,1,drop=FALSE]) * log(mu(X[i,,drop=FALSE], coeff)/(1 - mu(X[i,,drop=FALSE], coeff)))
		second_const = log(1 - mu(X[i,,drop=FALSE], coeff))

		result = result + first_const + second_const
	}

	return (result)
}

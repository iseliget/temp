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

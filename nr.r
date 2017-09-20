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
		first_const = as.numeric(Y[i,1,drop=FALSE]) * log(mu(X[i,,drop=FALSE], coeff)/(1 - mu(X[i,,drop=FALSE], coeff)))
		second_const = log(1 - mu(X[i,,drop=FALSE], coeff))

		result = result + first_const + second_const
	}

	return (result)
}

Fisher_scoring = function(X,Y) {
	coeff = matrix(0, nrow=ncol(X), ncol=1)

	for (i in 1:10) {
		print(paste("=========================Fisher Scoring ITERATION ", i , "======================================="))
		coeff = coeff - solve(Hessian_first(X,Y,coeff)) %*% gradient(X,Y,coeff)
		print("coef after iteration")
		print(coeff)
		print("log-likelikhood after iteration")
		print(log_likelihood(X,Y,coeff))
	}

	return (coeff)
}

Newton_rapshon = function(X,Y) {
	coeff = matrix(0, nrow=ncol(X), ncol=1)

	for (i in 1:10) {
		print(paste("=========================Newton Rapshon ITERATION ", i , "======================================="))
		coeff = coeff - solve(Hessian_first(X,Y,coeff)+Hessian_second(X,Y,coeff)) %*% gradient(X,Y,coeff)
		print("coef after iteration")
		print(coeff)
		print("log-likelikhood after iteration")
		print(log_likelihood(X,Y,coeff))
	}

	return (coeff)
}

Fisher_scoring(X,Y)
Newton_rapshon(X,Y)

rm(list = ls())

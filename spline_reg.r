data <- read.table(file="LAozone.data", header=T, sep=",")

attach(data)

id <- order(dpg)

N = function(j,t,knots,K) {
	# t is the input
	# j is the index of basis

	if (j==1) {
		return (1)
	} else if (j==2) {
		return (t)
	} else if (t <= knots[j-2]) { # dealing with N_j where j >= 3
		return (0)
	} else if (t <= knots[K-1]) {
		return ((t - knots[j-2])^3 / (knots[K] - knots[j-2]))
	} else if (t <= knots[K]) {
		return ((t-knots[j-2])^3 / (knots[K] - knots[j-2]) - (t - knots[K-1])^3 / (knots[K] - knots[K-1]))
	} else if (t > knots[K]) {
		return ( ((t-knots[j-2])^3 - (t-knots[K])^3)/(knots[K]-knots[j-2]) - ((t-knots[K-1])^3 - (t-knots[K])^3)/(knots[K] - knots[K-1]) )
	}
}

generate_matrix_N = function(X,n,p) {
	# returns a n*p matrix where
	# n should be number of data points
	# p should be number of basis

	matrix_N = matrix(0,n,p)
	knots = sort(unique(X)) # this is a vector

	for (i in 1:n) {
		for (j in 1:p) {
			matrix_N[i,j] = N(j,t=X[i],knots=knots,K=length(knots))
		}
	}

	return (matrix_N)
}

generate_matrix_Sigma = function(knots) {
	# returns a p*p matrix where
	# p is the number of basis

	p = length(knots)
	K = length(knots)
	matrix_Sigma = matrix(0,p,p)

	for (i in 3:p) {
		for (j in i:p) {
			term_1 = (36/((knots[K] - knots[j-2])*(knots[K] - knots[i-2]))) * (-knots[i-2]*knots[j-2]^2/2 + knots[i-2]*knots[j-2]*knots[K] - knots[i-2]*knots[K]^2/2 + knots[j-2]^3/6 - knots[j-2]*knots[K]^2/2 + knots[K]^3/3)
			term_2 = 12*(knots[K] - knots[K-1])
			term_3 = (6/(knots[K] - knots[j-2])) * (knots[K] - knots[K-1]) * (3*knots[j-2] - knots[K-1] - 2*knots[K])
			term_4 = (6/(knots[K] - knots[i-2])) * (knots[K] - knots[K-1]) * (3*knots[i-2] - knots[K-1] - 2*knots[K])

			matrix_Sigma[i,j] = term_1 + term_2 + term_3 + term_4

		}
	}

	for (j in 3:p) {
		for (i in j:p) {
			matrix_Sigma[i,j] = matrix_Sigma[j,i]
		}
	}

	return (matrix_Sigma)
}

# load data
data <- read.table(file="LAozone.data", header=T, sep=",")
attach(data)

X = (dpg + 70)/180
y = ozone

# generate matrices
lambda = 1.6e-3
matrix_N = generate_matrix_N(X, length(X), length(unique(X)))
matrix_Sigma = generate_matrix_Sigma(knots=sort(unique(X)))

smooth_matrix = matrix_N %*% solve(t(matrix_N)%*%matrix_N + lambda*matrix_Sigma) %*% t(matrix_N)
predicted_y = smooth_matrix %*% y

result = cbind(dpg,predicted_y)
result = as.data.frame(result)
print(result[order(dpg),])

print("===========================================")

f = function(lambda) sum(diag(matrix_N %*% solve(t(matrix_N)%*%matrix_N + lambda*matrix_Sigma, tol=1.5e-60) %*% t(matrix_N))) - 8.001308
lambda = uniroot(f,lower=0, upper=100, tol=1.5e-20)$root
matrix_N = generate_matrix_N(X, length(X), length(unique(X)))
matrix_Sigma = generate_matrix_Sigma(knots=sort(unique(X)))

smooth_matrix = matrix_N %*% solve(t(matrix_N)%*%matrix_N + lambda*matrix_Sigma) %*% t(matrix_N)
predicted_y = smooth_matrix %*% y

result = cbind(dpg,predicted_y)
result = as.data.frame(result)
print(result[order(dpg),])

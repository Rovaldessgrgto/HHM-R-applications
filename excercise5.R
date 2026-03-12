pn2pw <- function(m, lambda, delta) {
	working_eta <- log(lambda)
	if(m==1) return (working_eta)
	
	working_tau <- log(delta[-1]/ delta[1])
	return(c(working_eta, working_tau))
}

pw2pn <- function(m,parpw) {
	lambda <- exp(parpw[1:m])
	if(m==1) return( list(lambda=lambda, delta=1))

	tau <- parpw[(m+1):(2*m-1)]
	delta <- c(1, exp(tau))/ (1+sum(exp(tau)))
	return(list(lambda= lambda,delta=delta))
}

mllk <- function(parpw,x,m){
	pn <- pw2pn(m, parpw)
	n <- length(x)
	density_matrix <- matrix(NA, nrow = n, ncol= m)
	
	for(i in 1:m) {
		density_matrix[,i] <- pn$delta[i] * dpois(x, lambda = pn$lambda[i])
	}

	mixture_density <- rowSums(density_matrix)
	if(any(mixture_density <= 0)) return(1e10)
	return(-sum(log(mixture_density)))
}

weeks <- 0:241
sales <- c(1, 6, 9, 18, 14, 8, 8, 1, 6, 7, 3, 3, 1, 3, 4, 12, 8, 10, 8, 2,
17, 15, 7, 12, 22, 10, 4, 7, 5, 0, 2, 5, 3, 4, 4, 7, 5, 6, 1, 3,
4, 5, 3, 7, 3, 0, 4, 5, 3, 3, 4, 4, 4, 4, 4, 3, 5, 5, 5, 7,
4, 0, 4, 3, 2, 6, 3, 8, 9, 6, 3, 4, 3, 3, 3, 3, 2, 1, 4, 5,
5, 2, 7, 5, 2, 3, 1, 3, 4, 6, 8, 8, 5, 7, 2, 4, 2, 7, 4, 15,
15, 12, 21, 20, 13, 9, 8, 0, 13, 9, 8, 0, 6, 2, 0, 3, 2, 4, 4, 6,
3, 2, 5, 5, 3, 2, 1, 1, 3, 1, 2, 6, 2, 7, 3, 2, 4, 1, 5, 6,
8, 14, 5, 3, 6, 5, 11, 4, 5, 9, 9, 7, 9, 8, 3, 4, 8, 6, 3, 5,
6, 3, 1, 7, 4, 9, 2, 6, 6, 4, 6, 6, 13, 7, 4, 8, 6, 4, 4, 4,
9, 2, 9, 2, 2, 2, 13, 13, 4, 5, 1, 4, 6, 5, 4, 2, 3, 10, 6, 15,
5, 9, 9, 7, 4, 4, 2, 4, 2, 3, 8, 15, 0, 0, 3, 4, 3, 4, 7, 5,
7, 6, 0, 6, 4, 14, 5, 1, 6, 5, 5, 4, 9, 4, 14, 2, 2, 1, 5, 2,
6, 4
)
x_data <- sales
cat("\n--- B() SINGLE POISSON (m=1) ---\n")
intial_pw_1 <- pn2pw(m=1, lambda=1, delta=1)
fit1 <- nlm(f= mllk, p= intial_pw_1, x= x_data, m=1)
print(pw2pn(1, fit1$estimate))
cat("Negative Log-Likelihood Score:", fit1$minimum, "\n")

cat("\n--- (A) 2-COMPONENT POISSON (m=2) ---\n")
initial_pw_2 <- pn2pw(m=2, lambda=c(7,17), delta=c(0.5,0.5))
fit2 <- nlm(f=mllk, p = initial_pw_2, x= x_data, m=2)
print(pw2pn(2, fit2$estimate))
cat("Negative Log-Likelihood Score:", fit2$minimum, "\n")

cat("\n--- (C) 3-COMPONENT POISSON (m=3) ---\n")
initial_pw_3 <- pn2pw(m=3,lambda=c(5,10,15), delta=c(0.33,0.33,0.34))
fit3 <- nlm(f = mllk, p = initial_pw_3, x = x_data, m = 3)
print(pw2pn(3, fit3$estimate))
cat("Negative Log-Likelihood Score:", fit3$minimum, "\n")

cat("\n--- (C) 4-COMPNENT POISSON (m=4) ---\n")
initial_pw_4 <- pn2pw(m=4,lambda=c(5,10,15,20), delta=c(0.25,0.25,0.25,0.25))
fit4 <- nlm(f= mllk, p= initial_pw_4, x = x_data, m=4)
print(pw2pn(4,fit4$estimate))
cat("Negative Log-likelihood Score:", fit4$minimum, "\n")

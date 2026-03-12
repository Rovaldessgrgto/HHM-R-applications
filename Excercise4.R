pn2pw <- function(m,mu,sigma,delta) {
	working_mu <- mu
	working_eta <- log(sigma)
	working_tau <- log(delta[-1]/delta[1])
	
	return(c(working_mu, working_eta, working_tau))
}

pw2pn <- function(m,parpw){
	mu <- parpw[1:m]
	sigma <- exp(parpw[(m+1):(2*m)])
	tau <- parpw[(2 * m+1):(3 * m - 1)]
	delta <- c(1, exp(tau)) /(1 + sum(exp(tau))) 
	
	return(list(mu=mu, sigma=sigma,delta=delta))
}

mllk <- function(parpw,x,m){
	pn <- pw2pn(m,parpw)
	n <-length(x)
	density_matrix <- matrix(NA, nrow = n, ncol = m)
	for(i in 1:m){
		density_matrix[,i] <- pn$delta[i] *dnorm(x, mean = pn$mu[i], sd = pn$sigma[i])
	}

	mixture_density <-rowSums(density_matrix)
	
	if(any(mixture_density <= 0)) return(1e10)

	score <- -sum(log(mixture_density))

	cat(sprintf("Guessing Means: [%.2f, %.2f] | Score (Lower is better): %.2f\n",
		pn$mu[1],pn$mu[2], score))
	
	return(score)
}

table_1_3_data <- c(-0.39, 0.06, 0.12, 0.48, 0.94, 1.01, 1.67, 1.68, 1.76, 1.80, 
                    2.44, 3.25, 3.72, 4.12, 4.28, 4.60, 4.92, 5.28, 5.53, 6.22)

m<-2
initial_mu <- c(1,4)
initial_sigma <-c(1,1)
initial_delta <- c(0.5,0.5)

initial_pw <- pn2pw(m, initial_mu, initial_sigma, intial_delta)

cat("Starting Optimization...\n")
fit<- nlm(f=mllk, p=initial_pw, x = table_1_3_data, m=m)

optimal_parameters <- pw2pn(m,fit$estimate)

cat("\n--- FINAL OPTIMAL PARAMETERS ---\n")
print(optimal_parameters)
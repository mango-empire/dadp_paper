library(rstan)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

set.seed(123)
Omega <- rbind(
    c(1, 0.3, 0.2),
    c(0.3, 1, 0.1),
    c(0.2, 0.1, 1)
)
sigma <- c(1, 2, 3)
Sigma <- diag(sigma) %*% Omega %*% diag(sigma)
N <- 200
x <- mvtnorm::rmvnorm(N, c(-1,5,15), Sigma)
beta <- c(1, 2.4, -6, 4)

y <- cbind(1,x) %*% beta + rnorm(N)
y <- as.numeric(y)


standata <- list(J = ncol(x), N=N, x = x, y = y)


stanmodel1 <- stan_model("linear.stan")


samples1 <- sampling(stanmodel1, data = standata,
                     iter = 5000, warmup = 1000, chains = 1)

library(coda)
codasamples1 <- do.call(mcmc.list, 
                        plyr::alply(rstan::extract(samples1, 
                                                   pars=c("sigma", "Omega[1,2]", "Omega[1,3]", "Omega[2,3]"), 
                                                   permuted=FALSE), 2, mcmc))

library(mcmc)

lupost_factory <- function(x, y) function(theta) {
  alpha <- theta[1]
  beta <- theta[2]
  mux <- theta[3]
  sx <- exp(theta[4])
  sy <- exp(theta[5])
  
  t1 <- dnorm(x, mean = mux, sd = exp(sx), log = TRUE)
  t2 <- dnorm(y, mean = alpha + beta * x, sd = exp(sy), log = TRUE)
  
  sum(t1) + sum(t2)
}

set.seed(1)
n <- 100
alpha <- 2
beta <- -3
x <- rnorm(n, 1, 3)
y <- alpha + beta*x + rnorm(n,0,5)
#----------

lupost_factory <- function(x, y) function(theta) {
  alpha <- theta[1]
  beta <- theta[2]
  sy <- exp(theta[3])
  t1 <- dnorm(y, mean = alpha + beta * x, sd = sy, log = TRUE)
  sum(t1)
}

post_smpl <- function(dmat, theta) {
  x <- dmat[,1]
  y <- dmat[,2]
  mout <- metrop(lupost_factory(x,y), 
         initial = theta, 
         scale = c(.5,.5,.1),
         nbatch = 1)
  c(mout$batch)
}

lik_smpl <- function(theta) {
  alpha <- theta[1]
  beta <- theta[2]
  sy <- exp(theta[3])
  x <- runif(1)
  y <- rnorm(1, mean = alpha + beta*x, sd = sy)
  c(x,y)
}

st_calc <- function(dmat) {
  x <- dmat[,1]
  y <- dmat[,2]
  beta <- cov(x,y)/var(x)
  alpha <- mean(y) - beta * mean(x)
  c(cov(x,y), var(x), alpha)
}

priv_mech_factory <- function(epsilon) {
  function(sdp, xt) {
    sum(VGAM::dlaplace(sdp - xt, 0, 3/epsilon, TRUE))
  }
}


set.seed(1)
epsilon <- 1
n <- 200
alpha <- 5
beta <- -3
x <- runif(n)
y <- alpha + beta*x + rnorm(n,0,3)
sdp <- st_calc(cbind(x,y)) + VGAM::rlaplace(3, 0, 3/epsilon)


dmod <- new_privacy(post_smpl = post_smpl,
                    lik_smpl = lik_smpl,
                    ll_priv_mech = priv_mech_factory(epsilon),
                    st_calc = st_calc,
                    add = TRUE,
                    npar = 3)

tmp <- mcmc_privacy(dmod,
                    sdp = sdp,
                    nobs = n,
                    init_par = c(1,1,1),
                    niter = 50,
                    chains = 1,
                    varnames = c("alpha", "beta", "sigma"))

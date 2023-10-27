lik_smpl <- function(theta) {
  alpha <- theta[1]
  beta <- theta[2]
  x <- runif(1)
  y <- rnorm(1, mean = alpha + beta * x, sd = 3) 
  c(x,y)
}

post_smpl <- function(dmat, theta) {
  x <- dmat[,1]
  y <- dmat[,2]
  xm <- cbind(1, x)
  Si <- (1/9) * t(xm) %*% xm + (1/100) * diag(2)
  mu <- (1/9) * solve(Si) %*% t(xm) %*% y
  MASS::mvrnorm(1, mu = mu, Sigma = solve(Si))
}

st_calc <- function(dmat) {
  x <- dmat[,1]
  y <- dmat[,2]
  n <- length(y) - cov(x,y)/var(x)
  s1 <- (n-1) * cov(x,y)
  s2 <- (n-1) * var(x)
  s3 <- mean(y) 
  s4 <- mean(x)
  c(s1, s2, s3, s4)
}

priv_mech_factory <- function(n, epsilon) {
  function(sdp, xt) {
    delta1 <- (1- 1/n)
    delta3 <- 1/n
    t1 <- VGAM::dlaplace(sdp[1] - xt[1], 0, 3 * delta1/epsilon, TRUE)
    t2 <- VGAM::dlaplace(sdp[2] - xt[2], 0, 3 * delta1/epsilon, TRUE)
    t3 <- VGAM::dlaplace(sdp[3] - xt[3], 0, 3 * delta3/epsilon, TRUE)
    t4 <- VGAM::dlaplace(sdp[4] - xt[4], 0, 3 * delta3/epsilon, TRUE)
    sum(c(t1,t2,t3,t4))
  }
}


set.seed(1)
epsilon <- 1
n <- 50
alpha <- 5
beta <- -3
x <- runif(n)
y <- alpha + beta*x + rnorm(n,0,3)
delta1 <- 1 - 1/n
delta3 <- 1/n
sdp <- st_calc(cbind(x,y))
sdp[1] <- sdp[1] + VGAM::rlaplace(1, 0, 3 * delta1/epsilon)
sdp[2] <- sdp[2] + VGAM::rlaplace(1, 0, 3 * delta1/epsilon)
sdp[3] <- sdp[3] + VGAM::rlaplace(1, 0, 3 * delta3/epsilon)
sdp[4] <- sdp[4] + VGAM::rlaplace(1, 0, 3 * delta3/epsilon)



dmod <- new_privacy(post_smpl = post_smpl,
                    lik_smpl = lik_smpl,
                    ll_priv_mech = priv_mech_factory(n, epsilon),
                    st_calc = st_calc,
                    add = FALSE,
                    npar = 2)

tmp <- mcmc_privacy(dmod,
               sdp = sdp,
               nobs = n,
               init_par = c(1,1),
               niter = 1000,
               chains = 1,
               varnames = c("alpha", "beta"))

summary(tmp)

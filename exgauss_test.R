library(TMB)
library(brms)
library(ggplot2)
# simulate data
n <- 500
dat <- data.frame(y = rnorm(n, 300, 50) + rexp(n, 1/80))

# try to recover via different packages

# roughly 30s per chain, converges easily though
# my branch is marginally more performant
m_b <- brm(y ~ 1, data = dat,
           family = exgaussian(),
           prior = c(prior('normal(300, 20)', class = 'Intercept'),
                     prior('normal(50, 10)', class = 'sigma'),
                     prior('normal(80, 10)', class = 'beta')),
           seed = 1)

data_list <- list(y = dat$y)
params <- list(mu = 300, sigma = 50, nu = 80)

compile('exgauss_test.cpp')
dyn.load(dynlib('exgauss_test'))
obj <- MakeADFun(data_list, params, DLL = 'exgauss_test', silent = TRUE)

system.time(opt <- nlminb(obj$par, obj$fn, obj$gr))


# profile
par(mfrow = c(1,3))
sapply(X = 1:3, 
       FUN = function(x) {v <- tmbprofile(obj, x, trace = FALSE)
                          print(confint(v))
plot(tmbprofile(obj, x, trace = FALSE))})

# simulate
sim <- replicate(1000, {
    simdata <- obj$simulate(par = opt$par, complete = TRUE)
    obj2 <- MakeADFun(simdata, params, DLL = 'exgauss_test', silent = TRUE)
    nlminb(obj2$par, obj2$fn, obj2$gr)$par
})

df <- data.frame(estimate=as.vector(sim), parameter=names(obj$par)[row(sim)])
ggplot(df, aes(x = estimate, fill = parameter)) + 
geom_density() + 
facet_wrap(~parameter, scales = 'free_x')
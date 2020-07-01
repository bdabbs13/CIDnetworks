# R6 Testing

library(CIDnetworks)

# cluster means
theta <- 2*pi/3 * 0:2
unit_centers <- cbind(cos(theta), sin(theta))
centers <- 5*unit_centers

samples <- matrix(aperm(replicate(12, apply(centers, 1, rmvnorm, n = 1, sigma = 0.5*diag(2))),
                        perm = 3:1), ncol = 2)

prob <- plogis(2 - sqrt(as.matrix(dist(samples))))
adjmat <- apply(prob, 1, sapply, rbinom, n = 1, size = 1)

# FULL GIBBS TEST
big_test_net = CIDnetwork$new(adjmat, components = list(INTERCEPT(), SBM(n.groups = 3)))
test_gibbs <- CIDGibbs$new(
  cid.network = big_test_net,
  draws = 1000,
  extend.max = 0
)

test_output = test_gibbs$gibbs.full()
plot(test_output$log.likelihood, type='l')

test_mean = test_gibbs$gibbs.mean(test_output)
outgroups <- test_mean$components[[2]]$membership
plot(samples, col = outgroups)

test_mean$value()

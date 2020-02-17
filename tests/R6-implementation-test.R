# R6 Testing

library(CIDnetworks)

# FULL GIBBS TEST
big_test_net = CIDnetwork$new(matrix(rbinom(900, 1, .7), nrow=30))
test_gibbs <- CIDGibbs$new(
  cid.network = big_test_net,
  draws = 300,
  thin = 10,
  extend.max = 0
)

test_output = test_gibbs$gibbs.full()
hist(test_output$component_output$intercept, breaks=30)
abline(v=qnorm(mean(big_test_net$outcome)), lwd=4, col="green")

plot(test_output$log.likelihood, type='l')

test_mean = test_gibbs$gibbs.mean(test_output)
test_mean$value()

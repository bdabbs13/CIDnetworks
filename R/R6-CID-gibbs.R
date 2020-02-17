CIDGibbsOutput <- R6Class(
  classname="CIDGibbsOutput",
  public = list(
    #fields
    total.draws = 0,
    log.likelihood = NULL,
    residual.variance = NULL,
    component_output = list(),
    initialize = function(
      cid.network,
      total.draws=0
    ){
      
      self$total.draws = total.draws
      self$log.likelihood = rep(NA, total.draws)
      self$residual.variance = rep(NA, total.draws)

      if (length(cid.network$components) > 0){
        self$component_output = unlist(lapply(cid.network$components, 
                                              function(cc) cc$create.output.list(self$total.draws)), 
                                       recursive=FALSE)
      }
    }
  )
)


CIDGibbs <- R6Class(
  classname = "CIDGibbs",
  public = list(
    # fields
    cid.network = NULL,
    report.interval = 100,
    draws = 100,
    burnin = -1,
    thin = 10,
    use.auto.burnin = FALSE,
    auto.burnin.steps = 200,
    make.random.start = TRUE,
    extend.max = 10,
    verbose = 2,
    
    ##### Methods
    initialize = function(
      cid.network,
      report.interval = 100,
      draws = 100,
      thin = 10,
      burnin = NULL,
      auto.burnin.steps = 200,
      make.random.start = TRUE,
      extend.max = 10,
      verbose = 2
    ){
      self$cid.network <- cid.network
      
      self$draws <- draws
      self$thin <- as.integer(thin)
      if (self$thin < 1)
        stop("thin must be an integer greater than zero.")

      
      # Setting up Burn In Parameters
      self$burnin = burnin
      if (is.null(self$burnin)){
        self$use.auto.burnin = TRUE
        self$auto.burnin.steps <- auto.burnin.steps
      }
      self$report.interval <- report.interval

      self$make.random.start <- make.random.start
      self$extend.max <- extend.max
      self$verbose <- verbose
      
      return(0)
    },

    gibbs.full = function(){

      gibbs.output.list <- CIDGibbsOutput$new(self$cid.network, self$draws)

      if (self$make.random.start)
        self$cid.network$random.start()
      
      if (self$use.auto.burnin){
        self$auto.burnin()
      }else{
        self$manual.burnin()
      }
      
      converged <- FALSE
      extend.iter <- -1
      while (!converged & (extend.iter < self$extend.max)){
        extend.iter = extend.iter + 1
        if (extend.iter > 0)
          message("Extending Chain...")
        
        for (current.draw in 1:(self$draws * self$thin)) {
          self$cid.network$draw()
          
          if (current.draw %% self$thin == 0){
            current.index = current.draw / self$thin
            self$cid.network$update.output.list(gibbs.output.list, current.index)
            if ((self$report.interval > 0) & ((current.index %% self$report.interval) == 0)) {
              if(self$verbose > 1) message("CID ", current.index)
            }
          }
        }
        converged = self$convergence.test(gibbs.output.list)
      }

      if (self$verbose > 0){
        if (converged){
          message("CID Converged")
        }else{
          message("CID Convergence not detected after maximum extensions")
        }
      }
    
      return(gibbs.output.list)
    },

    auto.burnin = function(){
      repeat {
        if(self$verbose > 1) message("Auto-Burning In")
        
        loglikes <- rep(NA, self$auto.burnin.steps)
        start.time <- proc.time()[3]
        for (kk in 1:self$auto.burnin.steps) {
          self$cid.network$draw()
          loglikes[kk] <- self$cid.network$log.likelihood
        }
        elapsed.time <- proc.time()[3] - start.time
        
        # This does a weighted average of the first half and last half
        # then takes the difference of those two averages. The result
        # gives more weight to draws at the beginning and end of the samples
        slopes <- sum(loglikes*(1:self$auto.burnin.steps - ((self$auto.burnin.steps + 1) / 2)))
        
        # If the likelihood has stopped increasing (on average) stop burn in
        if (slopes < 0 & !is.na(slopes)){
          break
        }
      }
      
      if(self$verbose > 0){
        message("CID Auto-Burned In. Estimated Seconds Remaining: ",
                round((elapsed.time/ self$auto.burnin.steps) * self$thin * self$draws))
      }
    },
    
    manual.burnin = function(){
      if (self$verbose > 1){
        message("CID Performing Manual Burn In: ", self$burnin)
      }
      start.time <- proc.time()[3]
      for (kk in 1:self$burnin) {
        self$cid.network$draw()
      }
      elapsed.time <- proc.time()[3] - start.time
      
      if (self$verbose > 0){
        message("CID Burned In. Estimated Seconds Remaining: ",
                round((elapsed.time/ self$burnin) * self$thin * self$draws))
        
      }

    },
    
    convergence.test = function(gibbs.output.list){
      
      log.lik <- gibbs.output.list$log.likelihood
      chain.1 <- log.lik[1:(self$draws/10)]
      chain.2 <- log.lik[(self$draws/2 + 1):self$draws]
      sd.est <- sqrt(var(chain.1)/length(chain.1) +
                       var(chain.2)/length(chain.2))
      
      # BD: We should update the test level to be dynamic
      return((mean(chain.2) - mean(chain.1)) <= qnorm(0.999) * sd.est)
    },
    
    # TODO: Test this functionality
    gibbs.mean = function(gibbs.output.list){
      return(self$cid.network$get.mean.CIDnetwork(gibbs.output.list))
    }
    
  )
)



# #Removes magic number later. This should be the length of(intercept, residual, loglik, cutoffs, int.correlation).
# non.comp.count = function() 5,




############################################################################
############################################################################
#####################  Functions to be Implemented  ########################
############################################################################
############################################################################
# 
# gibbs.summary = function(gibbs.out) {
#   switched <- gibbs.switcheroo(gibbs.out$results)
#   s.sum <- function(int1) {
#     c(min = round(min(int1), 3), max = round(max(int1), 3),
#       estimated.mean= round(mean(int1), 3),
#       estimated.sd = round(sd(int1), 3),
#       round(quantile(int1, c(0.025, 0.975)), 3))
#   }
#   out <- list()
#   
#   out$intercept <- s.sum(unlist(switched[[1]]))
#   
#   out$residual.variance <- s.sum(unlist(switched[[2]]))
#   
#   out$log.likelihood <- s.sum(unlist(switched[[3]]))
#   
#   if (class.outcome == "ordinal") if (ordinal.count > 2) {
#     o1 <- matrix(unlist(switched[[4]]), nrow = ordinal.count-2)
#     out$ordinal.cutoffs <- t(apply(o1, 1, s.sum))
#     rownames(out$ordinal.cutoffs) <- paste("Cutoff", 1:(ordinal.count-2), 2:(ordinal.count-1), sep = "-")
#   }else{
#     out$ordinal.cutoffs <- NULL
#   }
#   
#   if (length(components) > 0) for (cc in 1:length(components)) {
#     ncc <- non.comp.count()
#     out[[cc+ncc]] <- components[[cc]]$gibbs.summary(switched[[cc+ncc]])
#     names(out)[cc+non.comp.count()] <- class(components[[cc]])
#   }
#   
#   out$CID.object <- gibbs.out$CID.object
#   out <- structure(out, class = "summary.CID.Gibbs")
#   return(out)
# },
# 
# 
# print.gibbs.summary = function(gibbs.sum){
#   
#   message("Intercept:"); print(gibbs.sum$intercept[-(1:2)])
#   if (gibbs.sum$residual.variance[4] != 0) {  #SD
#     message("Residual variance:")
#     print(gibbs.sum$residual.variance)
#   }
#   
#   message("Log likelihood:"); print(gibbs.sum$log.likelihood[-(1:2)])
#   
#   if(!is.null(gibbs.sum$residual.variance[-(1:2)])){
#     if(gibbs.sum$residual.variance["estimated.sd"] != 0){
#       message("Residual variance:")
#       print(gibbs.sum$residual.variance)
#     }
#   }
#   
#   if (class.outcome == "ordinal") if (ordinal.count > 2) {
#     message("Ordinal cutoffs:")
#     print(gibbs.sum$ordinal.cutoffs)
#   }
#   
#   if(!is.null(gibbs.sum$ordinal.cutoffs)){
#     message("Ordinal cutoffs:")
#     print(gibbs.sum$ordinal.cutoffs)
#   }
#   
#   if (length(components) > 0){
#     cov.count <- 0
#     for (cc in 1:length(components)) {
#       ##		message("Component ", class(components[[cc]]),":")
#       ix <- which(names(gibbs.sum) == class(components[[cc]]))
#       if(length(ix) > 1){
#         cov.count <- cov.count + 1
#         ix <- ix[cov.count]
#       }
#       for(ii in 1:length(ix)){
#         components[[cc]]$print.gibbs.summary(gibbs.sum[[ix]])
#       }
#       
#       #		print(round(gibbs.sum$class(components[[cc]])), 3)
#     }
#   }
#   return()#invisible(gibbs.sum))
# },
# 
# gibbs.plot = function(gibbs.out, DIC = NULL, which.plots = 1:(length(components)+5), auto.layout = TRUE) {
#   
#   switched <- gibbs.switcheroo(gibbs.out)
#   
#   if (auto.layout) {
#     if (intercept.v < 0.0001) which.plots <- which.plots[which.plots != 1]
#     if (class.outcome != "gaussian") which.plots <- which.plots[which.plots != 2]
#     if (ordinal.count == 2) which.plots <- which.plots[which.plots != 4]
#     if (!reciprocity.present) which.plots <- which.plots[which.plots != 5]
#     cols <- ceiling(length(which.plots)/2)
#     par(mfrow = c(2, cols))
#   }
#   
#   ## 1: Grand Intercept
#   if (1 %in% which.plots & intercept.v >= 0.0001) plot.default(unlist(switched[[1]]),
#                                                                main = "Grand Intercept",
#                                                                xlab = "Iteration",
#                                                                ylab = "Intercept",
#                                                                type = "l")
#   
#   ## 2: Residual variance. Skipped if not gaussian.
#   if (2 %in% which.plots & class.outcome =="gaussian")
#     plot.default(unlist(switched[[2]]), main = "Residual Variance",
#                  xlab = "Iteration",
#                  ylab = "Residual Variance")
#   main.label <- "Log-likelihood"; if (!is.null(DIC)) {
#     main.label <- paste0(main.label, ": DIC = ",signif(DIC[1], 5))
#     if (length(DIC)>= 2) main.label <- paste0(main.label, "\nDeviance of Average = ",signif(DIC[2], 5))
#     if (length(DIC)>= 3) main.label <- paste0(main.label, "\nEffective Parameter Count = ",signif(DIC[3], 5))
#   }
#   
#   ## 3: Log-likelihood
#   if (3 %in% which.plots) {
#     plot.default(unlist(switched[[3]]), main = main.label,
#                  xlab = "Iteration",ylab = "Log Likelihood",
#                  type = "l")
#   }
#   
#   ## 4: Ordinal cutoffs.
#   if (4 %in% which.plots & class.outcome =="ordinal" & ordinal.count>2) {
#     draws <- length(unlist(switched[[1]]))
#     xx <- sort(rep(1:draws, ordinal.count-2))
#     plot.default(c(1, xx), c(0, unlist(switched[[4]])), col = c(0, rep(1:(ordinal.count-2), draws)),
#                  main = "Ordinal Cutoff Values",
#                  xlab = "Iteration",
#                  ylab = "Cutoff Value")
#     abline(h = 0, col = 8)
#   }
#   
#   ## 5: Reciprocity
#   if (5 %in% which.plots & reciprocity.present) {
#     main.label <- "Reciprocal Correlation"
#     plot.default(unlist(switched[[5]]), main = main.label,
#                  xlab = "Iteration",
#                  ylab = "Reciprocal Correlation")
#     
#     
#   }
#   
#   #6 to(5 + length(components)): components!
#   if (length(components) > 0) for (cc in 1:length(components)) if ((non.comp.count()+cc) %in% which.plots) components[[cc]]$gibbs.plot(switched[[cc+non.comp.count()]])
#   
# },
# 
# 
# DIC = function(gibbs.out, add.parts = FALSE) {
#   #all.values <- gibbs.value(gibbs.out)
#   #deviance.of.average <- -2*log.likelihood.by.value(apply(all.values, 1, mean))
#   deviance.of.average <- -2 * gibbs.mean(gibbs.out)$log.likelihood
#   average.deviance <- -2 * mean(unlist(gibbs.switcheroo(gibbs.out)$log.lik))
#   #2*apply(all.values, 2, log.likelihood.by.value))
#   output <- c(DIC = 2*average.deviance - deviance.of.average)
#   if (add.parts) output <- c(output,
#                              deviance.of.average = deviance.of.average,
#                              effective.parameters = average.deviance - deviance.of.average,
#                              average.deviance = average.deviance)
#   return(output)
# },
# 
# marginal.loglikelihood = function(gibbs.out) {
#   all.values <- gibbs.value(gibbs.out)
#   model.log.likelihoods <- apply(all.values, 2, log.likelihood.by.value)
#   1/mean(1/model.log.likelihoods)
# },
# 
# pseudo.CV.loglikelihood = function(gibbs.out) {
#   all.values <- gibbs.value(gibbs.out)
#   model.log.likelihoods <- apply(all.values, 2, log.likelihood.by.value, sumup = FALSE)
#   each.like <- log(1/apply(1/exp(model.log.likelihoods), 1, mean))
#   sum(each.like)
# }
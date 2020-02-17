#####################################################################################
#
# Gibbs Sampler collection for CID, given the collection of input terms.

#source("COV-reference.R"); source("SBM-reference.R")e; source("LSM-reference.R"); source("SR-reference.R"); library(Rcpp); library(mvtnorm); library(msm); sourceCpp("../src/cid.cpp"); source("CID-basefunctions.R");

library(R6)

# .onAttach <- function(...) {
#   packageStartupMessage("CIDnetworks v 0.8.0")
# }


# Components --------------------------------------------------------------

# All subclasses available to CID.
# note: calling LSMcid(...) works like this already

# unwrap.CID.Gibbs <- function(gibbs.out) list.output.to.matrices(gibbs.out)


CIDnetwork <- R6Class(
  classname = "CIDnetwork",
  public = list(
    # fields
    n.nodes = NULL,
    edge.list = NULL,
    edge.list.rows = NULL,
    outcome = NULL,
    node.names = NULL,
    is.directed = FALSE,
    class.outcome = "binary",
    intermediate.outcome = NULL,
    # BD_NOTE: I think we should move intercept out of the CIDnetwork and make it a component
    intercept = 0,
    intercept.m = 0,
    intercept.v = 1000,
    residual.variance = 1,
    residual.variance.ab = c(.001, .001),
    reciprocity.present = FALSE,
    reciprocal.match = NULL,
    log.likelihood = -Inf,
    components = list(),
    comp.values = NULL,
    
    ##### Methods
    initialize = function(
      sociomatrix,
      node.names,
      is.directed,
      outcome.is.continuous=FALSE,
      intercept=0,
      components = list(),
      # prior parameters
      intercept.m = 0,
      intercept.v = 0.000000000001,
      residual.variance.ab = c(0.001, 0.001),
      # include.reciprocity = FALSE, #No reciprocity unless asked.
      verbose = 2,
      reinit=FALSE
    ) {
      ## Initial checks. This is largely redundant right now, but that's OK.
      
      if (outcome.is.continuous){
        self$class.outcome <- "gaussian"
      }
      
      ### Loading Network
      if (nrow(sociomatrix) != ncol(sociomatrix)){
        stop("The provided sociomatrix is not square.")
      }
      self$n.nodes <- nrow(sociomatrix)

      # BD: We should add a test that edge.list and outcome are the same length
      if (missing(is.directed)){
        is.directed = any(sociomatrix[!is.na(sociomatrix)] != t(sociomatrix)[!is.na(t(sociomatrix))])
      }
      
      self$is.directed = is.directed
      if (self$is.directed){
        self$edge.list <- make.arc.list(self$n.nodes)
        self$outcome <- t(sociomatrix)[non.diag(self$n.nodes)]
      } else {
        self$edge.list <- make.edge.list(self$n.nodes)
        self$outcome <- sociomatrix[u.diag(self$n.nodes)]  #just to be clear.
      }
      self$comp.values <- matrix(0, nrow = nrow(self$edge.list))
      
      ## Inferring Node Names
      if (missing(node.names)){
        node.names = 1:self$n.nodes
      }
      if (length(node.names) != self$n.nodes){
        stop("Length of node.names differs from the number of nodes in network.")
      }
      self$node.names = node.names

      ### Setting Parameters
      self$intercept.m <- intercept.m
      self$intercept.v <- intercept.v
      self$residual.variance.ab <- residual.variance.ab

      self$intercept = intercept
      ### Loading Components  
      if (class(components) != "list")
        components <- list(components)
      
      # BD: Might want to add a check if you accidentally use intercept and SBM, etc.
      # BD: Leaving reinit for now.  It's to pass the nodes and edge list in
      if (reinit) {
        for (kk in seq_along(components.t)) {
          components[[kk]]$reinitialize(n.nodes = .self$n.nodes,
                                        edge.list = .self$edge.list,
                                        node.names = .self$node.names)
        }
      }

      self$components <- components
      self$update.intermediate.outcome()
    },
    
    get.sociomatrix = function(){
      n = self$n.nodes
      socio = array(NA, dim = c(n,n))
      df = data.frame(self$edge.list, self$outcome)
      if (self$is.directed){
        for(i in 1:dim(df)[1]){
          socio[df[i,1],df[i,2]] = df[i,3]
        }
      } else{
        for(i in 1:dim(df)[1]){
          socio[df[i,1],df[i,2]] = socio[df[i,2],df[i,1]] =df[i,3]
        }
      }
      return(socio)
    },
    
    update.output.list = function(output.list, draw){
      output.list$log.likelihood[draw] = self$log.likelihood
      output.list$residual.variance[draw] = self$residual.variance
      output.list$intercept[draw] = self$intercept
      if (length(self$components) > 0){
        for (kk in self$components){
          self$components[kk]$update.output.list(output.list, draw)
        }
      }
    },
    
    get.mean.CIDnetwork = function(gibbs.output.list){

      intercept.mean <- mean(gibbs.output.list$intercept)
      mean.components = list()
      if(length(components) > 0){
        for(cc in 1:length(components)){
          mean.components[[cc]] = self$components[[cc]]$get.mean.component(gibbs.output.list)
        }
      }
      
      CID.mean <- CIDnetwork$new(
        sociomatrix = self$get.sociomatrix(),
        node.names = self$node.names,
        is.directed = self$is.directed,
        outcome.is.continuous = self$class.outcome == "gaussian",
        intercept=intercept.mean,
        components = mean.components,
        # prior parameters
        intercept.m = self$intercept.m,
        intercept.v = self$intercept.v,
        residual.variance.ab = self$residual.variance.ab,
        verbose = self$verbose,
        reinit = FALSE
        )

      CID.mean$update.log.likelihood()
      CID.mean$update.comp.values()
      return(CID.mean)
    },
    # pieces = function(
    #   include.name = FALSE
    # ) {
    #   if (length(self$components)>0) {
    #     c(
    #       list(
    #         intercept = self$intercept,
    #         residual.variance = self$residual.variance,
    #         log.likelihood = self$log.likelihood,
    #         ordinal.cutoffs = self$ordinal.cutoffs,
    #         int.correlation = self$int.correlation
    #       ),
    #       lapply(self$components, function(cc) cc$pieces(include.name))
    #       )
    #   } else {
    #     list(
    #       intercept = self$intercept,
    #       residual.variance = self$residual.variance,
    #       log.likelihood = self$log.likelihood,
    #       ordinal.cutoffs = self$ordinal.cutoffs,
    #       int.correlation = self$int.correlation
    #     )
    #   }
    #   },
    
    value = function() {
      rowSums(self$comp.values) + self$intercept
    },
    
    rem.values = function(kk) {
      if (kk>0){
        self$value() - self$comp.values[, kk]
      }else {
        self$value() - self$intercept
      }
    },
    
    update.comp.values = function() {
      if (length(self$components)>0){
        self$comp.values <- sapply(self$components, function(cc) cc$value())
      }
    },
    
    ## Update Z_ij.
    update.intermediate.outcome = function() {
      # BD: Is this necessary?
      self$update.comp.values()
      current.value <- self$value()
      
      if(class.outcome == "gaussian"){
        self$intermediate.outcome <- outcome
      }else if(class.outcome == "binary"){
        observed.ones = self$outcome == 1 & !is.na(self$outcome)
        observed.zeros = self$outcome == 0 & !is.na(self$outcome)
        
        self$intermediate.outcome[observed.ones] = rtnorm(
          n=sum(observed.ones), 
          mean=current.value[observed.ones], 
          sd=1,
          lower=0
        )

        self$intermediate.outcome[observed.zeros] = rtnorm(
          n=sum(observed.zeros), 
          mean=current.value[observed.zeros], 
          sd=1,
          upper=0
        )
      }else{
        stop("class.outcome must be either binary or gaussian")
      }
      
      if(any(is.na(self$outcome))){
        missing.values = is.na(self$outcome)
        self$intermediate.outcome[missing.values] = rnorm(
          n=sum(missing.values),
          mean=current.value[missing.values],
          sd=1
        )
      }
    },
    
    draw.intercept = function() {
      outcome.residual <- self$intermediate.outcome - self$rem.values(0);
      posterior.variance <- 1 / (nrow(self$edge.list) / self$residual.variance + 1 / self$intercept.v)
      posterior.mean <- posterior.variance * (sum(outcome.residual) / self$residual.variance + self$intercept.m / self$intercept.v)
      self$intercept <- rnorm(1, posterior.mean, sqrt(posterior.variance))
    },
    
    log.likelihood.by.value = function(value,
                                       residual.variance,
                                       sumup = TRUE,
                                       use.intermediate = FALSE)
    {
      if (missing(value))
        value = self$value()
      if (missing(residual.variance))
        residual.variance = self$residual.variance

      not.missing.indicator = !is.na(self$outcome)      
      cm <- function(pp) matrix(c(1, pp, pp, 1), nrow = 2)
      output <- 0 * self$intermediate.outcome
      if (class.outcome == "gaussian" | use.intermediate) {
        outcomeresid <- self$intermediate.outcome - value
        output <- dnorm(outcomeresid[not.missing.indicator], 0, sqrt(residual.variance), log = TRUE)
      } else if (class.outcome == "binary") {
        observed.outcomes = self$outcome[not.missing.indicator]
        observed.values = value[not.missing.indicator]

        # This is a clever way to get probit probabilities.
        # The outcome is that for edges with value 1 we get 1 - pnorm
        # and for outcome with value 0 we get pnorm
        breaker.lower <- c(-Inf, 0)
        breaker.upper <- c(0, Inf)
        breakers.lower <- breaker.lower[observed.outcomes + 1]
        breakers.upper <- breaker.upper[observed.outcomes + 1]
  
        output <- log(
          pnorm(breakers.upper,
                observed.values,
                sqrt(residual.variance)) -
            pnorm(breakers.lower,
                  observed.values,
                  sqrt(residual.variance)))
      }
      if(sumup) output <- sum(output)

      return(output)
    },
    
    update.log.likelihood = function() {
      self$log.likelihood <- self$log.likelihood.by.value()
    },
    
    draw.variance = function() {
      outcomeresid <- self$intermediate.outcome - self$value()
  
      self$residual.variance <-
        1/rgamma(1,
                 residual.variance.ab[1] + nrow(edge.list)/2,
                 residual.variance.ab[2] + sum(outcomeresid^2)/2)
  
    },
  
    draw = function(verbose = FALSE) {
  
      if (class.outcome != "gaussian") self$update.intermediate.outcome()
  
      self$update.comp.values()
      self$draw.intercept()

      if (length(self$components)>0) for (kk in 1:length(self$components)) {
        self$update.comp.values()
        # BD NOTE: This seems like a messy way to pass around this information
        self$components[[kk]]$outcome <- self$intermediate.outcome - self$rem.values(kk)
        self$components[[kk]]$residual.variance <- self$residual.variance
        self$components[[kk]]$draw()
      }

      if (class.outcome == "gaussian") {
        self$update.comp.values()
        self$draw.variance(verbose)
      }
      
      self$update.comp.values()
      self$update.log.likelihood()
    },
  
    random.start = function() {
      self$intercept <- rnorm(1, 0, 1)
      if (length(self$components)>0) for (kk in 1:length(self$components)) self$components[[kk]]$random.start()
      if (class.outcome == "gaussian") draw.variance()
      self$update.log.likelihood()
    }
  )
)

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
      components = NULL,
      # prior parameters
      residual.variance.ab = c(0.001, 0.001),
      verbose = 2
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
        is.directed = !isSymmetric(sociomatrix)
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
      self$residual.variance.ab <- residual.variance.ab

      ### Loading Components  
      if (is.null(components))
        components <- list(INTERCEPT())
      if (class(components) != "list")
        components <- list(components)
      if (length(components) < 1)
        stop("ERROR: components list must have at least one item")
      
      for (kk in seq_along(components)) {
        self$components[[kk]] <- components[[kk]]$create.component(
          self$n.nodes, self$edge.list, self$node.names
          )
      }

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
      for (component in self$components){
        component$update.output.list(output.list, draw)
      }
    },
    
    get.mean.CIDnetwork = function(gibbs.output.list){

      mean.components = list()
      for(cc in 1:length(self$components)){
        mean.components[[cc]] = self$components[[cc]]$get.mean.component(gibbs.output.list)
      }
      
      CID.mean <- CIDnetwork$new(
        sociomatrix = self$get.sociomatrix(),
        node.names = self$node.names,
        is.directed = self$is.directed,
        outcome.is.continuous = self$class.outcome == "gaussian",
        components = mean.components,
        # prior parameters
        residual.variance.ab = self$residual.variance.ab,
        verbose = self$verbose
        )

      CID.mean$update.log.likelihood()
      CID.mean$update.comp.values()
      return(CID.mean)
    },
    
    value = function() {
      rowSums(self$comp.values)
    },
    
    rem.value = function(kk) {
      return(rowSums(self$comp.values[, -kk, drop=FALSE]))
    },
    
    update.comp.values = function() {
      self$comp.values <- sapply(self$components, function(cc) cc$value())
    },
    
    ## Update Z_ij.
    update.intermediate.outcome = function() {
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
  
      if (class.outcome != "gaussian")
        self$update.intermediate.outcome()

      for (kk in 1:length(self$components)) {
        # browser()
        residual.outcome <- self$intermediate.outcome - self$rem.value(kk)
        self$components[[kk]]$draw(residual.outcome, self$residual.variance)
        self$comp.values[,kk] = self$components[[kk]]$value()
      }

      if (class.outcome == "gaussian")
        self$draw.variance(verbose)

      self$update.log.likelihood()
    },
  
    random.start = function() {
      for (component in self$components){
        component$random.start()
      }
      if (class.outcome == "gaussian")
        draw.variance()
      self$update.log.likelihood()
    }
  )
)

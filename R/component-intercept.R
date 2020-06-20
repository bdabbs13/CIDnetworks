#' Specify intercept component for a CIDnetwork model
#'
#' Specify the prior mean and variance and starting value for the intercept
#' component of your CIDnetwork model. Returns an `InterceptParams` object that
#' stores the your supplied arguments. This output should be passed as part of a
#' `list` of other components to the `components` argument of the 
#' `CIDnetwork$new()` initialization method to define your model.
#'
#' @param intercept.m Numeric scalar specifying the prior mean (default = 0).
#' @param intercept.v Numeric scalar specifying the prior variance (must be > 0;
#'   default = 1000).
#' @param intercept Numeric scalar defining the initial value for sampling 
#'   (default = 0).
#'
#' @return An `InterceptParams` object that can be fed in a list of other 
#'   components to the `components` argument of `CIDnetwork$new()`.
#' @export
#'
#' TODO: refer to old examples for ideas
#' @examples
INTERCEPT <- function(intercept.m=0, intercept.v=1000, intercept=0) {
  return(InterceptParams$new(intercept, intercept.m, intercept.v))
}

InterceptParams <- R6Class(
  classname = "InterceptParams",
  inherit = BaseParams,
  public = list(
    intercept=0,
    intercept.m = 0,
    intercept.v = 1000,
    
    initialize = function(
      intercept=0,
      intercept.m = 0,
      intercept.v = 1000
    ){
      self$intercept = intercept
      self$intercept.m = intercept.m
      self$intercept.v = intercept.v
    },
    
    create.component = function(n.nodes, edge.list, node.names){
      return(InterceptComponent$new(edge.list, self))
    }
  )
)


InterceptComponent <- R6Class(
  classname = "InterceptComponent",
  inherit = BaseComponent,
  public = list(
    intercept=0,
    intercept.m = 0,
    intercept.v = 1000,
    
    initialize = function(edge.list, params){
      # Note this component doesn't actually need the edge.list
      self$intercept = params$intercept
      self$intercept.m = params$intercept.m
      self$intercept.v = params$intercept.v
      private$n.edges = nrow(edge.list)
    },
    
    random.start = function(){
      self$intercept = rnorm(1, self$intercept.m, self$intercept.v)
    },
    
    draw = function(outcome, residual.variance) {
      posterior.variance <- 1 / (private$n.edges / residual.variance + 1 / self$intercept.v)
      posterior.mean <- posterior.variance * (sum(outcome) / residual.variance + self$intercept.m / self$intercept.v)
      self$intercept <- rnorm(1, posterior.mean, sqrt(posterior.variance))
    },
    
    create.output.list = function(total.draws){
      return (list(
        intercept.vector = rep(NA, total.draws)
      ))
    },
    
    update.output.list = function(gibbs.output.list, draw){
      gibbs.output.list$component_output$intercept.vector[draw] = self$intercept
    },

    get.mean.component = function(gibbs.output.list){
      return (InterceptParams$new(
        intercept=mean(gibbs.output.list$component_output$intercept.vector)
        ))
    },
    
    value = function(){
      return(rep(self$intercept, private$n.edges))
    }
  ),
  
  private = list(
    n.edges = NA
  )
)
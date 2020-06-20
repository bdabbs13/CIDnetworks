#' Class for storing and transporting component parameters
#' 
#' An object used internally by the CIDnetwork object when initializing the
#' components of the model. Returned by user-facing component functions like
#' `INTERCEPT`, which allow the user to define prior parameters and intial values
#' for the Gibbs sampler. Each component has a `Params` class that inherits
#' from the `BaseParams` class.
#' @aliases InterceptParams
#' @export
BaseParams <- R6Class(
  classname = "BaseParams",
  public = list(
    #' Initialize Params object
    initialize = function(){},
    
    #' TODO: Maybe add some documentation here for the edification of future 
    #'   developers
    #'
    #' @param n.nodes 
    #' @param edge.list 
    #' @param node.names 
    #'
    #' @return
    #' @export
    #'
    #' @examples
    create.component = function(n.nodes, edge.list, node.names){
      return(BaseComponent$new(edge.list, self))
    }
  )
)

#' Class for defining component interface to CIDnetwork object
#' 
#' An object used internally by the CIDnetwork object, which interfaces with
#' each component in the model via a standardized set of functions. This class
#' allows developers to create plug-and-play model components
#' by cutomizing the expected functions. Each component has a `Component` class
#' that inherits from the `BaseComponent` class.
#' @aliases InterceptComponent
#' @export
BaseComponent <- R6Class(
  classname = "BaseComponent",
  public = list(
    
    initialize = function(edge.list, params){
      private$n.edges = nrow(edge.list)
    },
    
    random.start = function(){
      invisible(self)
    },
    
    draw = function(outcome, residual.variance) {
      invisible(self)      
    },
    
    create.output.list = function(total.draws){
      invisible(self)
    },
    
    update.output.list = function(gibbs.output.list, draw){
      invisible(self)      
    },
    
    get.mean.component = function(gibbs.output.list){
      invisible(self)
    },
    
    value = function(){
      return(rep(NA, private$n.edges))
    }
  ),
  
  private = list(
    n.edges = NA
  )
)
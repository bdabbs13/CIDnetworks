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

    #' Initialize Params fields and object.
    #'
    #' Handle any type- and shape-checking here. If the the shape depends on
    #' `n.nodes`, it will have to be handled in the `Component` class
    #' `initialize` function.
    #'
    #' @return Nothing; update the object fields.
    #' @export
    #'
    #' @examples
    initialize = function() {},
    
    #' Creates a component with the supplied parameters
    #' 
    #' The CIDnetwork object will call this function when initializing the model.
    #'
    #' @param n.nodes Number of nodes in network passed to CIDnetwork/gibbs
    #' @param edge.list Sender-receiver-outcome matrix created from adjacency
    #'   matrix passed to CIDnetwork/gibbs
    #' @param node.names Node names passed by user to CIDnetwork/gibbs or
    #'   default sequential node number
    #'
    #' @return Initialized component object with n.nodes, edge.list, node.names
    #'   specified by CIDnetwork, and params (self) defined by user.
    #' @export
    #'
    #' @examples
    create.component = function(n.nodes, edge.list, node.names){
      return(BaseComponent$new(n.nodes, edge.list, node.names, self))
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
    
    initialize = function(n.nodes, edge.list, node.names, params){
      private$n.edges = nrow(edge.list)
    },
    
    #' Initialize component parameters to random values
    #'
    #' Specify how initial values for component parameters should be randomly
    #' generated, if the user has not specified initial values in the the
    #' `Params` object.
    #' 
    #' @return Nothing; just update the fields in the component object.
    #' @export
    #'
    #' @examples
    random.start = function() {
      invisible(self)
    },
    
    #' Compute a single MCMC draw for the component
    #'
    #' Specify how the parameters are updated for each step of MCMC sampling.
    #' Utilize the residualized (intermediate) outcome that remains after
    #' accounting for the other components to determine the new values.
    #'
    #' @param outcome The vector of residualized (intermediate) outcome for each
    #'   edge.
    #' @param residual.variance The residual variance of the outcome unaccounted
    #'   for by the components.
    #'
    #' @return Nothing; just update the fields in the component object.
    #' @export
    #'
    #' @examples
    draw = function(outcome, residual.variance) {
      invisible(self)      
    },
    
    #' Initialize storage for MCMC draws
    #' 
    #' Specify how the shape of the MCMC sample chain for each parameter. In
    #' general, the container should be a matrix or array with the first
    #' dimension being for the sampling iteration.
    #'
    #' @param total.draws The number of expected draws for an MCMC run.
    #'
    #' @return A list with an entry for each of your parameters of interst.
    #' @export
    #'
    #' @examples
    create.output.list = function(total.draws){
      invisible(self)
    },
    
    #' Update the MCMC chain in CID.gibbs
    #' 
    #' Specify how the MCMC sample chain is updated with each draw. Especially
    #' important if your parameters are higher-dimensional (matrix/array).
    #'
    #' @param gibbs.output.list A list of MCMC samples that will be output by
    #'   CID.gibbs
    #' @param draw Integer identifying which iteration the gibbs sampler is at.
    #'   Use this to index the gibbs.output.list to update it with the most
    #'   recent draw for each of your parameters.
    #'
    #' @return Nothing; just update the gibbs.output.list
    #' @export
    #'
    #' @examples
    update.output.list = function(gibbs.output.list, draw){
      invisible(self)      
    },
    
    #' Get the mean of the component parameters over MCMC samples
    #'
    #' Specify how the mean of each of the component parameters is to be
    #' computed and then pass it to a new `Params` object.
    #'
    #' @param gibbs.output.list A list of MCMC samples that will be output by
    #'   CID.gibbs
    #'
    #' @return A `Params` object with fields set to the means of your parameters
    #'   of interest over the MCMC samples.
    #' @export
    #'
    #' @examples
    get.mean.component = function(gibbs.output.list){
      invisible(self)
    },
    
    #' Value contributed by the component to the (intermediate) outcome.
    #'
    #' The (intermediate, if ordinal) outcome can be thought of as being
    #' composed of the components of your model and residual variance. Use this
    #' function to tell CIDnetwork the value that your component contributes to
    #' outcome.
    #' @return Numeric vector with length equal to the number of edges in the
    #'   network
    #' @export
    #'
    #' @examples
    value = function(){
      return(rep(NA, private$n.edges))
    }
  ),
  
  private = list(
    n.edges = NA
  )
)
BaseParams <- R6Class(
  classname = "BaseParams",
  public = list(
    initialize = function(){},
    
    create.component = function(n.nodes, edge.list, node.names){
      return(BaseComponent$new(edge.list, self))
    }
  )
)


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
# Interface Definition ---------------------------------------------------
LSM <- function(
  dimension = 2,
  latent.space.pos.m = 0,
  latent.space.pos.v.ab = c(0.001, 0.001),
  latent.space.pos = NULL,
  latent.space.pos.v = 100,
  latent.space.target = NULL,
  inverted.model = FALSE,
  tune = 0.1
) {
  return(LSMParams$new(dimension, latent.space.pos.m, latent.space.pos.v.ab,
                       latent.space.pos, latent.space.pos.v,
                       latent.space.target, inverted.model, tune))
}

# Param Definition -------------------------------------------------------
LSMParams <- R6Class(
  classname = "LSMParams",
  inherit = BaseParams,
  public = list(
    dimension = 2,
    latent.space.pos.m = 0,
    latent.space.pos.v.ab = c(0.001, 0.001),
    latent.space.pos = NULL,
    latent.space.pos.v = 100,
    latent.space.target = NULL,
    inverted.model = FALSE,
    tune = 0.1,
    
    initialize = function(
      dimension = 2,
      latent.space.pos.m = 0,
      latent.space.pos.v.ab = c(0.001, 0.001),
      latent.space.pos = NULL,
      latent.space.pos.v = 100,
      latent.space.target = NULL,
      inverted.model = FALSE,
      tune = 0.1
    ) {
      self$dimension <- dimension
      self$latent.space.pos.m <- latent.space.pos.m
      
      if (length(latent.space.pos.v.ab) != 1) {
        stop("latent.space.pos.v.ab must have two elements c(alpha, beta) for",
             " each of the inv.gamma prior distribution parameters.")
      }
      self$latent.space.pos.v.ab <- latent.space.pos.v.ab
      
      # handle these in LSMComponent once we have n.nodes
      self$latent.space.pos <- latent.space.pos
      self$latent.space.target <- latent.space.target
      
      self$latent.space.pos.v <- latent.space.pos.v
      self$inverted.model <- inverted.model
      self$tune <- tune
    },
    
    create.component = function(n.nodes, edge.list, node.names) {
      return(LSMComponent$new(n.nodes, edge.list, node.names, self))
    }
  )
)

# Component Class Definition --------------------------------------------
LSMComponent <- R6Class(
  classname = "LSMComponent",
  inherit = BaseComponent,
  public = list(
    dimension = 2,
    n.nodes = NA,
    edge.list = NULL,
    node.names = NULL,
    latent.space.pos = NULL,
    latent.space.pos.m = 0,
    latent.space.pos.v = 100,
    latent.space.pos.v.ab = c(0.001, 0.001),
    latent.space.target = NULL,
    inverted.model = FALSE,
    tune = 0.1,
    
    initialize = function(n.nodes, edge.list, node.names, params) {
      self$dimension <- params$dimension
      self$latent.space.pos.m <- params$latent.space.pos.m
      self$latent.space.pos.v <- params$latent.space.pos.v
      self$latent.space.pos.v.ab <- params$latent.space.pos.v.ab
      self$inverted.model <- params$inverted.model
      self$tune <- params$tune
      
      if (!is.null(params$latent.space.pos)) {
        if (all(dim(params$latent.space.pos) != c(n.nodes, params$dimension))) {
          stop("latent.space.pos must be shape nrow = n.nodes, nrow = dimension")
        } else {
          if (is.null(params$latent.space.pos)) {
            self$latent.space.target <- params$latent.space.pos
          }
          self$latent.space.pos <- params$latent.space.pos
        }
      }
      
      if (!is.null(params$latent.space.target)) {
        if (all(dim(params$latent.space.target) != c(n.nodes, params$dimension))) {
          stop("latent.space.target must be shape nrow = n.nodes, nrow = dimension")
        } else {
          if (is.null(params$latent.space.pos)) {
            self$latent.space.pos <- params$latent.space.target
          }
          self$latent.space.target <- params$latent.space.target
        }
      }
      
      self$n.nodes <- n.nodes
      self$edge.list <- edge.list
      self$node.names <- node.names
      private$mult.factor <- 2 * self$inverted.model - 1
    },
    
    random.start = function() {
      "Generates a random latent positions to initiate the MCMC sampling. Places
      Them in the correct object field."
      if (!is.null(params$latent.space.pos)) {
        n <- self$dimension * self$n.nodes
        inits <- rnorm(n, self$latent.space.pos.m, sqrt(self$latent.space.pos.v))
        self$latent.space.pos <- matrix(inits, nrow = self$n.nodes)
        self$latent.space.target <- self$latent.space.pos
      }
    },
    
    log.likelihood = function(outcome,
                              residual.variance,
                              parameters,
                              edges = seq_len(nrow(self$edge.list))) {
      "Computes log-likelihood of outcome given with the mean being predicted by
    a regression on the latent-space distances."
      meanpart <- self$value.ext(parameters, edges)
      sum(dnorm(outcome[edges], meanpart, sqrt(residual.variance), log = TRUE))
    },
    
    draw = function(outcome, residual.variance) {
      "Computes every step for each MCMC draw. MH step for latent space
      positions, Gibbs update for variance."
      
      lsdim <- self$dimension
      
      latent.space.pos.hold <- self$latent.space.pos
      latent.space.pos.prop <- self$latent.space.pos
      for (dd in 1:self$n.nodes) {
        walk <- as.vector(rmvnorm(1, rep(0, lsdim), diag(self$tune^2, lsdim)))
        latent.space.pos.prop[dd, ] <- latent.space.pos.prop[dd, ] + walk
        
        edge.list.rows <- row.list.maker(self$edge.list)
        llike.prop <- self$log.likelihood(
          outcome = outcome,
          residual.variance = residual.variance,
          parameters = list(latent.space.pos.prop, private$mult.factor),
          edges = edge.list.rows[[dd]]
        )
        ldens.pos.prop <- dnorm(
          x = latent.space.pos.prop[dd, ],
          mean = self$latent.space.pos.m,
          sd = sqrt(self$latent.space.pos.v),
          log = TRUE
        )
        log.dens.prop <- llike.prop + sum(ldens.pos.prop)
        
        # compute log density of previously accepted draw
        llike.orig <- self$log.likelihood(
          outcome = outcome,
          residual.variance = residual.variance,
          parameters = list(latent.space.pos.hold, private$mult.factor),
          edges = edge.list.rows[[dd]]
        )
        ldens.pos.orig <- dnorm(
          x = latent.space.pos.hold[dd, ],
          mean = self$latent.space.pos.m,
          sd = sqrt(self$latent.space.pos.v),
          log = TRUE
        )
        log.dens.orig <- llike.orig + sum(ldens.pos.orig)
        
        # accepts or rejects proposal by comparing ratio to uniform(0, 1)
        if (log.dens.prop - log.dens.orig > -rexp(1)) {
          latent.space.pos.hold[dd, ] <- latent.space.pos.prop[dd, ]
        } else {
          # what does this do? latent.space.pos.prop is never used again
          latent.space.pos.prop[dd, ] <- latent.space.pos.hold[dd, ]
        }
      }
      
      ## Procrustean transformation using random start as a target matrix
      latent.space.pos.hold <- private$procrustean.post(
        latent.space.pos.hold,
        self$latent.space.target
      )
      rownames(latent.space.pos.hold) <- self$node.names
      
      # store accepted proposal/previously accepted draw as new draw
      self$latent.space.pos <- latent.space.pos.hold
      
      # Gibbs update draw from posterior of latent space position variance
      a.factor <- self$n.nodes * self$dimension / 2
      b.factor <- sum(self$latent.space.pos^2) / 2
      g.draw <- rgamma(n = 1,
                       shape = a.factor + self$latent.space.pos.v.ab[1],
                       rate = b.factor + self$latent.space.pos.v.ab[2]
      )
      self$latent.space.pos.v <- 1 / g.draw # stores in field
    },
    
    create.output.list = function(total.draws) {
      init <- replicate(total.draws, matrix(NA, self$dimension, self$n.nodes))
      init.dimperm <- aperm(init, perm = 3:1)
      return(list(
        latent.space.pos = init.dimperm,
        latent.space.pos.v = rep(NA, total.draws)
      ))
    },
    
    
    update.output.list = function(gibbs.output.list, draw) {
      gibbs.output.list$component_output$latent.space.pos[draw, , ] = self$latent.space.pos
      gibbs.output.list$component_output$latent.space.pos.v[draw] = self$latent.space.pos.v
    },
    
    get.mean.component = function(gibbs.output.list) {
      return(LSMParams$new(
        dimension = self$dimension,
        latent.space.pos = apply(gibbs.output.list$component_output$latent.space.pos, c(2, 3), mean),
        latent.space.pos.v = mean(gibbs.output.list$component_output$latent.space.pos.v),
        latent.space.target = self$latent.space.target,
        inverted.model = self$inverted.model
      ))
    },
    
    value = function() {
      "Returns latent distance between nodes multiplied by multiplicative factor."
      val <- private$mult.factor *
        edge.list.distance(self$latent.space.pos, self$edge.list)
      return(val)
    },
    
    value.ext = function(parameters,
                         edges = seq_len(nrow(self$edge.list))) {
      "Returns latent distance between nodes multiplied by multiplicative factor.
    parameters argument accepts list with latent.space.pos in first position,
    latent.space.v in second position, and multiplicative factor in third."
      # NOTES: this seems to be an error, parameters[[2]] from pieces() refers to
      #        latent.space.pos.v rather than mult.factor (which is in position 3)
      parameters[[2]] * edge.list.distance(parameters[[1]], rbind(self$edge.list[edges, ]))
    },    
    
    gibbs.summary = function(gibbs.out) {
      "Report mean latent space position across thinned gibbs draws."
      lsp.all <- gibbs.output.list$component_output$latent.space.pos
      # take mean across draws and reshape into matrix with nodes as rows
      output <- matrix(apply(lsp.all, c(2, 3), mean), nrow = n.nodes)
      rownames(output) <- node.names
      colnames(output) <- paste0("pos", 1:ncol(output))
      return(output)
    },
    
    print.gibbs.summary = function(gibbs.sum) {
      "Print summary without returning anything."
      # NOTE: Is this intended?
      message("Mean Latent Space Positions:")
      print(gibbs.sum)
      return()
    },
    
    gibbs.plot = function(gibbs.out, ...) {
      "Calls LSMcid plot method for mean gibbs draws."
      get.sum <- gibbs.summary(gibbs.out)
      main_title <- "Mean Latent Space Positions from Gibbs Sampler"
      private$plot(get.sum, main = main_title, ...)
    }
  ),
  # Private Methods and Attributes
  private = list(
    mult.factor = NA,
    plot = function(pos = self$latent.space.pos, ...) {
      "Plots the latent space positions. pos argument accepts latent space
       position matrix. Defaults to latent space position already in class."
      latent.space.plot(pos, labels = self$node.names, ...)
    },
    procrustean.post = function(latent.space.pos,
                                latent.space.target,
                                recenter = TRUE) {       
      if (recenter) {
        latent.space.pos <- scale(latent.space.pos, scale = FALSE)
        latent.space.target <- scale(latent.space.target, scale = FALSE)
      }
      
      projection = t(latent.space.target) %*% latent.space.pos
      ssZ = svd(projection)
      transformation = ssZ$v %*% t(ssZ$u)
      latent.space.pos = latent.space.pos %*% transformation
      
      return(latent.space.pos)
    }
  )
)

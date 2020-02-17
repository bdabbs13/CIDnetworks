#' @include CID-basefunctions.R plotting-pieces.R
NULL
# functions from CID-basefunctions.R used here:
#   make.edge.list, row.list.maker, edge.list.distance
# functions from plotting-pieces.R used here:
#   latent.space.plot

# !CTM!: WHAT DOES THIS STUFF DO UP HERE?
# library(Rcpp); library(mvtnorm); library(msm);
# sourceCpp ("../src/cid.cpp"); source("CID-basefunctions.R");
# Standard Latent Space Model: Reference Class

# Class Declaration -------------------------------------------------------

LSMcid <- setRefClass(
  Class = "LSMcid",
  fields = list(
    dimension = "numeric",
    latent.space.pos = "matrix",
    mult.factor = "numeric",
    latent.space.target = "matrix",
    latent.space.pos.m = "numeric", # mean
    latent.space.pos.v = "numeric", # variance
    latent.space.pos.v.ab = "numeric", # variance inverse gamma priors
    tune = "numeric",
    node.names = "character",
    n.nodes = "numeric",
    outcome = "numeric",
    edge.list = "matrix",
    residual.variance = "numeric",
    edge.list.rows = "list"
  )
)

# Methods ---------------------------------------------------------

LSMcid$methods(
  ## Internal Methods ---------------------------------------------
  initialize = function(dimension = 1,
                        n.nodes = 10,
                        edge.list = make.edge.list(n.nodes),
                        edge.list.rows = row.list.maker(edge.list),
                        residual.variance = 1,
                        outcome = numeric(0),
                        latent.space.pos = matrix(rnorm(dimension * n.nodes), nrow = n.nodes),
                        latent.space.target = matrix(rep(0, dimension * n.nodes), nrow = n.nodes),
                        latent.space.pos.m = 0,
                        latent.space.pos.v = 100,
                        latent.space.pos.v.ab = c(0.001, 0.001),
                        tune = 0.1,
                        inverted.model = FALSE,
                        generate = FALSE) {
    .self$n.nodes <<- n.nodes
    .self$edge.list <<- edge.list
    .self$edge.list.rows <<- edge.list.rows
    .self$residual.variance <<- residual.variance
    .self$node.names <<- as.character(1:n.nodes)
    .self$dimension <<- dimension
    .self$latent.space.pos <<- latent.space.pos
    .self$latent.space.pos.m <<- latent.space.pos.m
    .self$latent.space.pos.v <<- latent.space.pos.v
    .self$latent.space.pos.v.ab <<- latent.space.pos.v.ab
    .self$mult.factor <<- 2 * inverted.model - 1
    .self$tune <<- tune

    if (generate) {
      outcome <<- rnorm(nrow(edge.list), value(), sqrt(residual.variance))
    } else {
      outcome <<- outcome
    }
  },
  pieces = function(include.name = FALSE) {
    out <- list(
      latent.space.pos = latent.space.pos,
      latent.space.pos.v = latent.space.pos.v,
      mult.factor = mult.factor
    )
    class(out) <- "LSMout"
    out
  },
  value = function() {
    "Returns latent distance between nodes multiplied by multiplicative factor."
    mult.factor * edge.list.distance(latent.space.pos, edge.list)
  },
  value.ext = function(parameters = pieces(),
                       edges = seq_len(nrow(edge.list))) { # slightly slower.
    "Returns latent distance between nodes multiplied by multiplicative factor.
    parameters argument accepts list with latent.space.pos in first position,
    latent.space.v in second position, and multiplicative factor in third."
    # NOTES: this seems to be an error, parameters[[2]] from pieces() refers to
    #        latent.space.pos.v rather than mult.factor (which is in position 3)
    parameters[[2]] * edge.list.distance(parameters[[1]], rbind(edge.list[edges, ]))
  },
  random.start = function() {
    "Generates a random latent positions to initiate the MCMC sampling. Places
    Them in the correct object field."
    latent.space.pos <<- matrix(rnorm(dimension * n.nodes, 0, sqrt(latent.space.pos.v)), nrow = n.nodes)
    latent.space.target <<- latent.space.pos
  },
  log.likelihood = function(parameters = pieces(),
                            edges = 1:nrow(edge.list)) {
    "Computes log-likelihood of outcome given with the mean being predicted by
    a regression on the latent-space distances."
    meanpart <- value.ext(parameters, edges)
    sum(dnorm(outcome[edges], meanpart, sqrt(residual.variance), log = TRUE))
  },
  ## User Methods ------------------------------------------------
  show = function() {
    "Prints the latent space positions and multiplicative factor."
    message("t(latent.space.pos):")
    print(t(latent.space.pos))
    message("latent.space.pos.v:")
    print(latent.space.pos.v)
    
    message("mult.factor:")
    print(mult.factor)
  },
  plot = function(pos = latent.space.pos, ...) {
    "Plots the latent space positions. pos argument accepts latent space
    position matrix. Defaults to latent space position already in class."
    latent.space.plot(pos, labels = node.names, ...)
  },
  plot.network = function(color = outcome, ...) {
    "Plots the edge.list as a network. color argument is a numeric vector
    that controls the node color."
    image.netplot(edge.list, color, node.labels = node.names, ...)
  },
  draw = function(verbose = 0, mh.tune = tune, langevin.mode = FALSE) {
    "Computes every step for each MCMC draw. MH step for latent space positions,
    Gibbs update for variance."
    # NOTES: This could be more self-documenting by breaking it into more
    #        functions? It's also very slow looking. Not exploiting
    #        vectorization.
    
    # dimension of latent space positions
    # (why do we need a field called "dimension" then??)
    lsdim <- dim(latent.space.pos)[2]
    
    # initializing vector copies to store the previous positions and the
    # proposal.
    latent.space.pos.hold <- latent.space.pos
    latent.space.pos.prop <- latent.space.pos
    for (dd in 1:n.nodes) { # this is really slow, a lot of these steps can be vectorized
      ## proposal for latent space positions.
      walk <- as.vector(rmvnorm(1, rep(0, lsdim), diag(mh.tune^2, lsdim)))
      latent.space.pos.prop[dd, ] <- latent.space.pos.prop[dd, ] + walk
      
      # computes density of proposal
      llike.prop <- log.likelihood(
        parameters = list(latent.space.pos.prop, mult.factor),
        edges = edge.list.rows[[dd]]
      )
      ldens.pos.prop <- dnorm(
        x = nt.space.pos.prop[dd, ],
        mean = 0,
        sd = sqrt(latent.space.pos.v),
        log = TRUE
      )
      log.dens.prop <- llike.prop + sum(ldens.pos.prop)
      
      # compute log density of previously accepted draw
      llike.orig <- log.likelihood(
        parameters = list(latent.space.pos.hold, mult.factor),
        edges = edge.list.rows[[dd]]
      )
      ldens.pos.orig <- dnorm(
        x = latent.space.pos.hold[dd, ],
        mean = 0,
        sd = sqrt(latent.space.pos.v),
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
    
    ## Rotate back.
    ##        latent.space.pos.hold <- postprocess.latent.positions(latent.space.pos.hold)
    
    ## Procrustean transformation
    ## Using random start as a target matrix
    latent.space.pos.hold <- postprocess.latent.positions(
      latent.space.pos.hold,
      latent.space.target
    )
    rownames(latent.space.pos.hold) <- node.names
    
    # store accepted proposal/previously accepted draw as new draw
    latent.space.pos <<- latent.space.pos.hold
    
    # Gibbs update draw from posterior of latent space position variance
    # NOTE: Is this a gibbs update?
    a.factor <- nrow(latent.space.pos) * ncol(latent.space.pos) / 2
    b.factor <- sum(latent.space.pos^2) / 2
    g.draw <- rgamma(n = 1,
                     shape = a.factor + latent.space.pos.v.ab[1],
                     rate = b.factor + latent.space.pos.v.ab[2]
    )
    latent.space.pos.v <<- 1 / g.draw # stores in field
  },
  gibbs.full = function(report.interval = 0,
                        draws = 100,
                        burnin = 0,
                        thin = 1,
                        make.random.start = FALSE, ...) {
    "Performs each sampling step for specified number of draws and burn-in."
    out <- list()
    if (make.random.start) random.start()
    for (kk in 1:(draws * thin + burnin)) {
      # MHstep latent space positions and gibbs update for position variance
      draw(...)
      index <- (kk - burnin) / thin
      # NOTE: could use modulus operator for this
      if (kk > burnin & round(index) == index) {
        out[[index]] <- c(pieces(), list(log.likelihood = log.likelihood()))
        if (report.interval > 0) if (index %% report.interval == 0) message("LSM ", index)
      } else if (round(index) == index) {
        if (report.interval > 0) if (index %% report.interval == 0) message("LSM burnin ", index)
      }
    }
    return(out)
  },
  gibbs.value = function(gibbs.out) { 
    "Computes the value for each gibbs draw."
    sapply(gibbs.out, value.ext)
  },
  gibbs.summary = function(gibbs.out) {
    "Report mean latent space position across thinned gibbs draws."
    lsp.all <- sapply(gibbs.out, function(gg) gg$latent.space.pos)
    # take mean across draws and reshape into matrix with nodes as rows
    output <- matrix(apply(lsp.all, 1, mean), nrow = n.nodes)
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
  gibbs.mean = function(gibbs.out) {
    "Returns new LSM object with latent space positions equal to the mean."
    get.sum <- gibbs.summary(gibbs.out)
    lsm.out <- LSM(
      dimension = dimension, 
      n.nodes = n.nodes,
      edge.list = edge.list,
      edge.list.rows = edge.list.rows,
      residual.variance = residual.variance,
      outcome = outcome,
      latent.space.pos = get.sum, # mean of gibbs draws for latent pos
      latent.space.pos.m = latent.space.pos.m,
      latent.space.pos.v = latent.space.pos.v,
      latent.space.pos.v.ab = latent.space.pos.v.ab,
      tune = tune,
      inverted.model = mult.factor == 1 # this seems to be wrong?, should be -1?
    )
    return(lsm.out)
  },
  gibbs.plot = function(gibbs.out, ...) {
    "Calls LSMcid plot method for mean gibbs draws."
    get.sum <- gibbs.summary(gibbs.out)
    plot(get.sum, main = "Mean Latent Space Positions from Gibbs Sampler", ...)
  },
  gibbs.node.colors = function(gibbs.out) {
    "Colors the nodes from the gibbs draws."
    rep("#DDDDFF", n.nodes)
  },
  reinitialize = function(n.nodes = NULL,
                          edge.list = NULL, node.names = NULL) {
    "Reinitializes the object, clearing previous data, with an option to change
    the input data."
    if (!is.null(n.nodes)) n.nodes <<- n.nodes
    if (!is.null(edge.list)) {
      edge.list <<- edge.list
      edge.list.rows <<- row.list.maker(edge.list)
    }
    if (nrow(latent.space.pos) != .self$n.nodes) {
      message("Reinitializing LSM Positions")
      latent.space.pos <<- matrix(rnorm(dimension * n.nodes, 0, sqrt(latent.space.pos.v)), nrow = n.nodes)
      #   adjust.lsp()
    }
    
    if (ncol(latent.space.pos) > .self$n.nodes) {
      warning("Latent space has more dimensions than nodes. Not impossible, but probably counter-productive.")
    }
    
    if (!is.null(node.names)) {
      if (length(node.names) == .self$n.nodes) node.names <<- node.names
    } else {
      node.names <<- as.character(1:.self$n.nodes)
    }
  }
)

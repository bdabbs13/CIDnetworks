# Interface Definition ---------------------------------------------------
#' Specify SBM component for a CIDnetwork model
#'
#' @description
#' Specify the prior mean and variance for the block-block tie matrix, the
#' multinomial prior on the membership vector, and their starting values for the
#' SBM component of your CIDnetwork model. Returns an `SBMParams` object that
#' stores the your supplied arguments. This output should be passed as part of a
#' `list` of other components to the `components` argument of the
#' `CIDnetwork$new()` initialization method to define your model.
#' 
#' @details
#' In a stochastic block model (SBM), the probability of observing a network tie
#' between two nodes depends on which groups (called blocks) the nodes are
#' members of, and the relationship between the groups. Nodes in the same block
#' are more likely to have a tie than nodes in different block; and some times
#' certain pairs of blocks are more likely to have nodes that share ties than
#' other pairs.
#' 
#' The model is written in notation as
#' \deqn{Y_{ij} \sim Bernoulli(g_i^T B g_j)}
#' \deqn{g_i \sim Multinomial(p_i)}
#' \deqn{B_{\ell m} \sim Beta(\mu_{\ell m}, \nu_{\ell m})}
#' 
#' This package uses the mean-variance parameterization to specify the prior on
#' the block-block tie matrix.
#' 
#' The variables in this model correspond to the following arguments for this
#' function
#' \itemize{
#'   \item \eqn{g_i} is the group/block membership vector. You can pass the
#'   initial value in the `membership` argument.
#'   \item \eqn{p_i} is the prior for the group/block membership vector which is
#'   the `membership.a` argument. In this package, the prior can be specified
#'   for each node.
#'   \item \eqn{B_{\ell m}} is the block-block tie matrix. You can pass the
#'   initial value in the `block.matrix` argument.
#'   \item \eqn{\mu_{\ell m}} corresponds to the `block.matrix.m` argument.
#'   \item \eqn{\nu_{\ell m}} corresponds to the `block.matrix.v` argument.
#' }
#' 
#' The block.matrix priors can be used to incorporate information about how
#' likely inter-group ties are.
#'   
#' @param n.groups Integer scalar specifying the number of groups/blocks.
#'   (default = 2)
#' @param block.matrix.m Numeric matrix of size `G X G` specifying the mean of
#'   the beta prior distribution for each block-block tie probability. `G` is
#'   the number of groups/blocks. (default = matrix of all 0s)
#' @param block.matrix.v Numeric matrix of size `G X G` specifying the variance
#'   of the beta prior distribution for each block-block tie probability. `G` is
#'   the number of groups/blocks. (default = matrix of all 1000s; "weak" prior,
#'   nearly uniform)
#' @param membership.a Numeric matrix of size `N X G` specifying the prior
#'   probability for each node to be a member of each group. `N` is the number
#'   of nodes and `G` is the number of groups/blocks.  Rows do not need to sum
#'   to 1, as MCMC takes care of normalization. However, it is highly
#'   recommended to specify rows that sum to 1 if not using a uniform
#'   probability (like the default) so that you are fully aware of the
#'   normalized prior probability you are specifying. (default = matrix of all
#'   1s)
#' @param symmetric.b Logical scalar, if `TRUE`, will force `block.matrix` to be
#'   symmetric. Ensures that group sending-receiving probabilities are
#'   identical. (default = TRUE)
#' @param stong.block Logical scalar, if `TRUE`, will force `block.matrix` to
#'   have larger diagonal entries than off-diagonal. Ensures that intra-group
#'   tie probabilities are strong than ties with nodes in other groups. (default
#'   = FALSE)
#' @param block.matrix Numeric matrix of size `G X G` specifying the initial
#'   block-block tie probability matrix. `G` is the number of groups/blocks.
#'   (default = NULL, which randomly generates the matrix from the prior)
#' @param membership Numeric vector of length `N` specifying the initial value
#'   for the group/block membership vector. `N` is the number of nodes. (default
#'   = NULL, which randomly generates the vector from the prior)
#'
#' @return An `SBMParams` object that can be fed in a list of other components
#'   to the `components` argument of `CIDnetwork$new()`.
#' @export
#'
#' TODO: refer to old examples for ideas
#' @examples
SBM <- function(
  n.groups = 2,
  block.matrix.m = matrix(0, nrow = n.groups, ncol = n.groups),
  block.matrix.v = matrix(10000, nrow = n.groups, ncol = n.groups),
  membership.a = NULL,
  symmetric.b = TRUE,
  strong.block = FALSE,
  block.matrix = NULL,
  membership = NULL
) {
  return(SBMParams$new(n.groups, block.matrix.m, block.matrix.v,
                       membership.a, symmetric.b, strong.block,
                       block.matrix, membership))
}

# Param Definition -------------------------------------------------------
SBMParams <- R6Class(
  classname = "SBMParams",
  inherit = BaseParams,
  public = list(
    n.groups = 2,
    block.matrix.m = NULL,
    block.matrix.v = NULL,
    membership.a = NULL,
    symmetric.b = TRUE,
    strong.block = FALSE,
    block.matrix = NULL,
    membership = NULL,
    
    initialize = function(
      n.groups = 2,
      block.matrix.m = matrix(0, nrow = n.groups, ncol = n.groups),
      block.matrix.v = matrix(10000, nrow = n.groups, ncol = n.groups),
      membership.a = NULL,
      symmetric.b = TRUE,
      strong.block = FALSE,
      block.matrix = NULL,
      membership = NULL
    ) {
      self$n.groups <- n.groups
      
      if (!all(dim(block.matrix.m) == c(n.groups, n.groups))) {
        stop("block.matrix.m must be a square matrix with nrow = ncol = n.groups")
      }
      self$block.matrix.m <- block.matrix.m
      
      if (!all(dim(block.matrix.v) == c(n.groups, n.groups))) {
        stop("block.matrix.v must be a square matrix with nrow = ncol = n.groups")
      }
      self$block.matrix.v <- block.matrix.v
      
      # have to deal with shape in the Component class when we know n.nodes
      self$membership.a <- membership.a
      self$symmetric.b <- symmetric.b
      self$strong.block <- strong.block
      
      if (!is.null(block.matrix)) {
        if (!all(dim(block.matrix) == c(n.groups, n.groups))) {
          stop("block.matrix must be a square matrix with nrow = ncol = n.groups")
        } else if (symmetric.b & !isSymmetric(block.matrix)) {
          stop("with symmetric.b = TRUE, block.matrix initial value must be symmetric")
        }
      }
      self$block.matrix <- block.matrix
      
      # have to deal with shape in the Component class when we know n.nodes
      self$membership <- membership
    },
    
    create.component = function(n.nodes, edge.list, node.names) {
      return(SBMComponent$new(n.nodes, edge.list, node.names, self))
    }
  )
)

# Component Class Definition --------------------------------------------
SBMComponent <- R6Class(
  classname = "SBMComponent",
  inherit = BaseComponent,
  public = list(
    n.groups = 2,
    block.matrix.m = NULL,
    block.matrix.v = NULL,
    membership.a = NULL,
    symmetric.b = TRUE,
    strong.block = FALSE,
    block.matrix = NULL,
    membership = NULL,
    n.nodes = NA,
    edge.list = NULL,
    edge.list.rows = NULL,
    node.names = NULL,
    
    initialize = function(n.nodes, edge.list, node.names, params) {
      self$n.groups <- params$n.groups
      self$block.matrix.m <- params$block.matrix.m
      self$block.matrix.v <- params$block.matrix.v
      
      if (is.null(params$membership.a)) {
        params$membership.a = matrix(1, nrow = n.nodes, ncol = params$n.groups)
      } else if (!all(dim(params$membership.a) == c(n.nodes, params$n.groups))) {
        stop("membership.a must be matrix with shape nrow = n.nodes, ncol = n.groups")
      }
      self$membership.a <- params$membership.a
      self$symmetric.b <- params$symmetric.b
      self$strong.block <- params$strong.block
      self$block.matrix <- params$block.matrix
      
      if (!is.null(params$membership) & length(params$membership) != n.nodes) {
        stop("membership must be a numeric vector with length = n.nodes")
      } 
      self$membership <- params$membership
      
      self$n.nodes <- n.nodes
      self$edge.list <- edge.list
      self$edge.list.rows <- row.list.maker(edge.list)
      self$node.names <- node.names
    },
    
    
    random.start = function() {
      "Initializes group membership matrix and block-block tie probability matrix"
      # QUESTION: does not use membership.a prior probability?
      if (length(unique(self$membership)) != self$n.groups) { # TRUE when self$membership = NULL
        if (!is.null(self$membership)) {
          message("Resampling group membership. Some classes had 0 members.")
        }
        membership <- sample(self$n.groups, self$n.nodes, replace = TRUE) 
        while (length(unique(membership)) != self$n.groups) {
          message("Resampling group membership. Some classes had 0 members.")
          membership <- sample(self$n.groups, self$n.nodes, replace = TRUE)
        }
        self$membership <- membership
      } 
      
      # QUESTION: does not use membership.a prior parameters, isn't beta distributed?
      if (is.null(self$block.matrix)) {
        block.matrix <- matrix(rnorm(self$n.groups * self$n.groups, 0, 1), nrow = self$n.groups)
      }
      if (self$strong.block) {
        pivots <- sapply(1:n.groups, function(kk) {
          max(c(block.matrix[kk, -kk], block.matrix[-kk, kk]))
        })
        diag(block.matrix) <- pivots + rexp(self$n.groups)
      }
      self$block.matrix <- block.matrix
    },
    
    log.likelihood = function(outcome,
                              residual.variance,
                              parameters,
                              edges = seq_len(nrow(self$edge.list))) {
      "Computes log-likelihood of residualized outcome after accounting for 
      other components."
      meanpart <- self$value.ext(parameters, edges)
      sum(dnorm(outcome[edges], meanpart, sqrt(residual.variance), log = TRUE))
    },
    
    draw = function(outcome, residual.variance) {
      "Computes every step for each MCMC draw."
      if (length(outcome) != nrow(self$edge.list)) {
        stop("SBM: outcome and edge.list have different lengths.")
      }
      
      # Drawing membership vector
      b.memb <- self$membership
      
      # TODO: Should probably be a private method
      get_log_posterior <- function(gg, node) {
        b.memb[node] <- gg
        ## BD: We don't need to calculate the likelihood for all edges, just the ones
        sender_node <- b.memb[self$edge.list[self$edge.list.rows[[node]], 1]]
        receiver_node <- b.memb[self$edge.list[self$edge.list.rows[[node]], 2]] - 1
        piece <- self$block.matrix[sender_node + receiver_node * self$n.groups]
        
        log_prior <- log(self$membership.a[node, gg])
        log_likelihood <- sum(dnorm(outcome[self$edge.list.rows[[node]]],
                                    piece, sqrt(residual.variance),
                                    log = TRUE
        ))
        log_posterior <- log_prior + log_likelihood
        
        return(log_posterior)
      }
      
      for (node in sample(1:self$n.nodes)) {
        ## Note 2014-12-05: If a move empties a class, disallow it. -AT
        if (any(b.memb[-node] == b.memb[node])) {
          log.pp.vec <- sapply(1:self$n.groups, get_log_posterior, node = node)
          log.pp.vec <- log.pp.vec - max(log.pp.vec) # Handling Underflow
          b.memb[node] <- sample(1:self$n.groups, 1, prob = exp(log.pp.vec))
        }
      }
      self$membership <- b.memb
      
      # Drawing block-block tie probs.
      membership.pairs <- cbind(
        b.memb[self$edge.list[, 1]],
        b.memb[self$edge.list[, 2]]
      )
      for (ss in 1:self$n.groups) {
        for (rr in 1:self$n.groups) {
          if (!self$symmetric.b | (self$symmetric.b & ss <= rr)) {
            if (self$symmetric.b) {
              picks <- which((b.memb[self$edge.list[, 1]] == ss &
                                b.memb[self$edge.list[, 2]] == rr) |
                               (b.memb[self$edge.list[, 2]] == ss &
                                  b.memb[self$edge.list[, 1]] == rr))
            } else {
              picks <- which(b.memb[self$edge.list[, 1]] == ss &
                               b.memb[self$edge.list[, 2]] == rr)
            }
            
            if (length(picks) > 0) {
              prec.b <- ((length(picks) / residual.variance) +
                           (1 / self$block.matrix.v[ss, rr]))
              var.b <- 1 / prec.b
              
              likelihood_piece <- sum(outcome[picks]) / residual.variance
              prior_piece <- (self$block.matrix.m[ss, rr] / self$block.matrix.v[ss, rr])
              mean.b <- var.b * (likelihood_piece + prior_piece)
            } else {
              var.b <- 0.5^2
              mean.b <- 0
            }
          }
        }
        
        if (!self$strong.block) {
          output <- rnorm(1, mean.b, sqrt(var.b))
        } else {
          if (ss == rr) {
            pivot <- max(c(
              self$block.matrix[ss, -ss],
              self$block.matrix[-ss, ss]
            ))
            output <- rtnorm(1, mean.b, sqrt(var.b),
                             lower = pivot
            )
          } else {
            pivot <- min(c(
              self$block.matrix[ss, ss],
              self$block.matrix[rr, rr]
            ))
            output <- rtnorm(1, mean.b, sqrt(var.b),
                             upper = pivot
            )
          }
        }
        
        self$block.matrix[ss, rr] <- output
        if (self$symmetric.b) {
          self$block.matrix[rr, ss] <- output
        }
      }
      private$rotate()
    },
    
    create.output.list = function(total.draws) {
      init <- replicate(total.draws, matrix(NA, self$n.groups, self$n.groups))
      init.dimperm <- aperm(init, perm = 3:1)
      return(list(
        block.matrix = init.dimperm,
        membership = matrix(NA, nrow = total.draws, ncol = self$n.nodes)
      ))
    },
    
    update.output.list = function(gibbs.output.list, draw) {
      gibbs.output.list$component_output$block.matrix[draw, , ] <- self$block.matrix
      gibbs.output.list$component_output$membership[draw, ] <- self$membership
    },
    
    get.mean.component = function(gibbs.output.list) {
      return(SBMParams$new(
        n.groups = self$n.groups,
        block.matrix.m = self$block.matrix.m,
        block.matrix.v = self$block.matrix.v,
        membership.a = self$membership.a,
        symmetric.b = self$symmetric.b,
        strong.block = self$strong.block,
        block.matrix = apply(gibbs.output.list$component_output$block.matrix, c(2, 3), mean),
        membership = apply(gibbs.output.list$component_output$membership, 2, mean)
      ))
    },
    
    value = function() {
      self$block.matrix[self$membership[self$edge.list[, 1]] +
                          dim(self$block.matrix)[1] * (self$membership[self$edge.list[, 2]] - 1)]
    },
    
    value.ext = function(parameters,
                         edges = seq_len(nrow(self$edge.list))) {
      sbm.matrix <- parameters[[1]]
      sbm.matrix[parameters[[2]][self$edge.list[edges, 1]] +
                   dim(sbm.matrix)[1] * (parameters[[2]][self$edge.list[edges, 2]] - 1)]
    },
    
    gibbs.summary = function(gibbs.out) {
      membs1 <- {
        d1 <- sapply(gibbs.out,
          function(gg) {number.to.vector(gg$membership, nrow(gg$block.matrix))}
        )
        matrix(apply(d1, 1, mean), ncol = n.nodes)
      }
      colnames(membs1) <- self$node.names
      this.block.matrix <- matrix(
        apply(
          sapply(gibbs.out,
                 function(gg) { c(gg$block.matrix)}),
          1, mean),
        nrow = n.groups
      )
      modal.membership <- apply(membs1, 2, which.max)
      return(list(
        membership = membs1,
        modal.membership = modal.membership,
        block.matrix = this.block.matrix
      ))
    },
    
    print.gibbs.summary = function(gibbs.sum) {
      message("Probabilistic block memberships:")
      print(gibbs.sum$membership)
      
      message("Modal block memberships:")
      print(gibbs.sum$modal.membership)
      
      message("Block value matrix:")
      print(gibbs.sum$block.matrix)
      
      return()
    },
    
    gibbs.plot = function(gibbs.out, ...) {
      get.sum <- gibbs.summary(gibbs.out)
      block.membership.plot(
        get.sum$membership,
        get.sum$block.matrix,
        node.labels = self$node.names,
        main = "SBM Summary from Gibbs Sampler",
        col = gibbs.node.colors(gibbs.out),
        ...
      )
    }
  ),
  # Private Methods and Attributes
  private = list(
    mult.factor = NA,
    
    rotate = function() {
      rotation <- SBM.ID.rotation(self$membership, self$n.groups)
      self$membership <- rotation[self$membership]
      self$block.matrix <- SBM.rotate.block(self$block.matrix, rotation)
    },
    
    gibbs.node.colors = function(gibbs.out, colors = (1:self$n.groups) + 1) {
      get.sum <- gibbs.summary(gibbs.out)
      return(colors[get.sum$modal.membership])
    }
  )
)

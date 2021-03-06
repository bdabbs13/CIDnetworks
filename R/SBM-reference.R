### Single-membership Stochastic Block Model: Reference Class
### input: ID.labels for a number of nodes, chosen to be the "dominant" ones.
### output: the permutation to switch labels to a simpler convention.

SBMcid <-
    setRefClass(
        "SBMcid",
        fields = list(
            n.groups="numeric",

            block.matrix="matrix",
            block.matrix.m="matrix",
            block.matrix.v="matrix",

            membership="integer",
            membership.a="matrix",

            symmetric.b="logical",
            strong.block="logical", # 2014-12-08, ACT -- should the diagonal always be greater?
            group.pairs="matrix",

            ## inherited from main. Must fix later, but OK for now.
            node.names="character",
            n.nodes="numeric",
            outcome="numeric",
            edge.list="matrix",
            residual.variance="numeric",
            edge.list.rows="list"
        ),

        methods=list(

            initialize = function (

                                   n.groups=1,

                                   n.nodes=10,
                                   edge.list=make.edge.list(n.nodes),
                                   edge.list.rows=row.list.maker(edge.list),
                                   residual.variance=1,
                                   outcome=numeric(0),

                                   block.matrix=matrix(0, nrow=n.groups, ncol=n.groups),
                                   block.matrix.m=matrix(0, nrow=n.groups, ncol=n.groups),
                                   block.matrix.v=matrix(10000, nrow=n.groups, ncol=n.groups),

                                   membership=sample(n.groups, n.nodes, replace=TRUE),
                                   symmetric.b=TRUE,
                                   strong.block=FALSE,

                                   membership.a=matrix(1, nrow=n.nodes, ncol=n.groups),

                                   generate=FALSE

                                   ) {


                .self$n.nodes <<- n.nodes
                .self$edge.list <<- edge.list
                .self$edge.list.rows <<- edge.list.rows
                .self$node.names <<- as.character(1:.self$n.nodes)

                .self$n.groups <<- n.groups

                .self$block.matrix <<- block.matrix
                .self$block.matrix.m <<- block.matrix.m
                .self$block.matrix.v <<- block.matrix.v
                .self$membership <<- as.integer(membership)

                .self$membership.a <<- membership.a
                .self$residual.variance <<- residual.variance

                .self$group.pairs <<- makeEdgeListSelfies(n.groups)

                .self$symmetric.b <<- symmetric.b
                .self$strong.block <<- strong.block


                if (symmetric.b) {
                    b.block <- .self$block.matrix
                    b.block[u.diag(.self$n.groups)] <- b.block[l.diag(.self$n.groups)]
                    .self$block.matrix <<- as.matrix(b.block)
                }

                if (generate) .self$generate() else .self$outcome <<- outcome

            },

            reinitialize = function (n.nodes=NULL,

                                     edge.list=NULL, node.names=NULL) {

                if (!is.null(n.nodes)) n.nodes <<- n.nodes
                if (!is.null(edge.list)) {
                    edge.list <<- edge.list
                    edge.list.rows <<- row.list.maker(edge.list)
                }
                if (!is.null(node.names)) {
                    if (length(node.names) == .self$n.nodes) node.names <<- node.names
                } else node.names <<- as.character(1:.self$n.nodes)


                if (n.groups > n.nodes) {
                    warning("SBM: Resetting number of groups to one less than the number of nodes.")
                    n.groups <<- n.nodes - 1
                    block.matrix <<- matrix(0, nrow=n.groups,ncol=n.groups)

                    membership <<- sample(n.groups, n.nodes, replace=TRUE)
                    membership.a <<- matrix(1, nrow=n.nodes, ncol=n.groups)
                }

                if (length(membership) != n.nodes) {
                    message ("Reinitializing SBM Memberships")
                    membership <<- sample(n.groups, n.nodes, replace=TRUE)
                }
                if(!identical(dim(membership.a), c(n.nodes,n.groups))){
                    membership.a <<- matrix(1, nrow=n.nodes, ncol=n.groups)
                }

                while (length(unique(membership)) != n.groups) {
                    message ("reinitialize: Group membership omits classes.")
                    membership <<- sample(n.groups, n.nodes, replace=TRUE)
                }


            },

            pieces = function (include.name=TRUE) {
                out <- list (block.matrix=block.matrix, membership=membership)
                class(out) <- "SBMout"
                out
            },

            show = function () {
                message("block.matrix:"); print(block.matrix)
                message("membership:"); print(membership)
            },
            plot = function (memb=membership, block=block.matrix, ...) {
                single.membership.plot (memb, block, node.labels=node.names, ...)
            },
            plot.network = function (color=outcome, ...) {
                image.netplot (edge.list, color, node.labels=node.names, ...)
            },



            value = function () {
                block.matrix[membership[edge.list[,1]] +
                             dim(block.matrix)[1]*(membership[edge.list[,2]]-1)]
            },
            value.ext = function (parameters=pieces(), edges=1:nrow(edge.list)) { #slightly slower.
                sbm.matrix <- parameters[[1]]
                sbm.matrix[parameters[[2]][edge.list[edges,1]] +
                           dim(sbm.matrix)[1]*(parameters[[2]][edge.list[edges,2]]-1)]
            },

            generate = function () {
                outcome <<- rnorm(nrow(edge.list), value(), sqrt(residual.variance))
            },

            log.likelihood = function(parameters=pieces(), edges=1:nrow(edge.list)) {
                meanpart <- value.ext (parameters, edges)
                sum(dnorm(outcome[edges], meanpart, sqrt(residual.variance), log=TRUE))
            },


            random.start = function () {
                membership <<- sample(n.groups, n.nodes, replace=TRUE)
                while (length(unique(membership)) != n.groups) {
                    message ("reinitialize: Group membership omits classes.")
                    membership <<- sample(n.groups, n.nodes, replace=TRUE)
                }

                block.matrix <<- matrix(rnorm(n.groups*n.groups, 0, 1), nrow=n.groups)
                if (strong.block) {
                    pivots <- sapply(1:n.groups, function(kk) max(c(block.matrix[kk, -kk],
                                                                    block.matrix[-kk, kk])))
                    diag(block.matrix) <<- pivots + rexp(n.groups)
                }

            },

            rotate = function () {
                rotation <- SBM.ID.rotation(membership, n.groups)
                membership <<- rotation[membership]
                block.matrix <<- SBM.rotate.block(block.matrix, rotation)
            },

            draw = function (verbose=0) {

                if (length(outcome) != nrow(edge.list))
                    stop ("SBM: outcome and edge.list have different lengths.")

                b.memb <- membership

                if (verbose>1) print(b.memb)

                get_log_posterior <- function(gg, node){
                    b.memb[node] <- gg
                    ## BD: We don't need to calcula te the likelihood for all edges, just the ones
                    sender_node <-  b.memb[edge.list[edge.list.rows[[node]], 1]]
                    receiver_node <- b.memb[edge.list[edge.list.rows[[node]], 2]] - 1
                    piece <- block.matrix[sender_node + receiver_node * n.groups]

                    log_prior <- log(membership.a[node,gg])
                    log_likelihood <- sum(dnorm(outcome[edge.list.rows[[node]]],
                                                piece, sqrt(residual.variance), log=TRUE))
                    log_posterior <- log_prior + log_likelihood

                    return(log_posterior)

                }

                for (node in sample(1:n.nodes)){
                    ## Note 2014-12-05: If a move empties a class, disallow it. -AT
                    if (any(b.memb[-node] == b.memb[node])){
                        log.pp.vec <- sapply(1:n.groups, get_log_posterior, node=node)
                        log.pp.vec <- log.pp.vec - max(log.pp.vec) # Handling Underflow
                        b.memb[node] <- sample (1:n.groups, 1, prob=exp(log.pp.vec))
                    }
                }
                if (verbose>1) print(b.memb)
                membership <<- b.memb

                ## draw block probs.
                membership.pairs <- cbind(b.memb[edge.list[,1]],
                                          b.memb[edge.list[,2]])
                for (ss in 1:n.groups){
                    for (rr in 1:n.groups){
                        if (!symmetric.b | (symmetric.b & ss <= rr)) {
                            if (symmetric.b) {
                                picks <- which((b.memb[edge.list[,1]] == ss &
                                                b.memb[edge.list[,2]] == rr) |
                                               (b.memb[edge.list[,2]] == ss &
                                                b.memb[edge.list[,1]] == rr))
                            } else {
                                picks <- which(b.memb[edge.list[,1]] == ss &
                                               b.memb[edge.list[,2]] == rr)
                            }

                            if (length(picks) > 0) {
                                prec.b <- ((length(picks) / residual.variance) +
                                           (1 / block.matrix.v[ss,rr]))
                                var.b <- 1 / prec.b

                                likelihood_piece <- sum(outcome[picks]) / residual.variance
                                prior_piece <- (block.matrix.m[ss,rr] / block.matrix.v[ss,rr])
                                mean.b <- var.b * (likelihood_piece + prior_piece)
                                ## mean.b <- var.b*(sum(outcome[picks]) / residual.variance + block.matrix.m[ss,rr]/block.matrix.v[ss,rr])

                            }else{
                                var.b <- 0.5^2; mean.b <- 0
                            }
                        }
                    }

                    if (!strong.block) {
                        output <- rnorm(1, mean.b, sqrt(var.b))
                    } else {
                        if (ss == rr) {
                            pivot <- max (c(block.matrix[ss, -ss],
                                            block.matrix[-ss, ss]))
                            output <- rtnorm(1, mean.b, sqrt(var.b),
                                             lower=pivot)
                        } else {
                            pivot <- min (c(block.matrix[ss, ss],
                                            block.matrix[rr, rr]))
                            output <- rtnorm(1, mean.b, sqrt(var.b),
                                             upper=pivot)
                        }
                    }

                    block.matrix[ss,rr] <<- output
                    if (symmetric.b) block.matrix[rr,ss] <<- output
                }
                rotate()
            },

            gibbs.full = function (report.interval=0, draws=100, burnin=0, thin=1,
                                   make.random.start=FALSE) {
                out <- list()
                if (make.random.start) random.start()
                for (kk in 1:(draws*thin+burnin)) {
                    draw();
                    index <- (kk-burnin)/thin
                    if (kk > burnin & round(index)==index) {
                        out[[index]] <- c(pieces(), list(log.likelihood=log.likelihood()))
                        if ((report.interval > 0) && (index %% report.interval == 0)){
                            message("SBM ",index)
                        }
                    }
                }
                return(out)
            },

            gibbs.value = function (gibbs.out) sapply(gibbs.out, function(gg) {
                value.ext (gg)
            }),

            gibbs.mean = function(gibbs.out){
                get.sum <- gibbs.summary(gibbs.out)

                return(SBM(n.groups=n.groups,
                           n.nodes=n.nodes,
                           edge.list=edge.list,
                           edge.list.rows=edge.list.rows,
                           residual.variance=residual.variance,
                           outcome=outcome,

                           block.matrix=get.sum$block.matrix,
                           block.matrix.m=block.matrix.m,
                           block.matrix.v=block.matrix.v,

                           membership=get.sum$modal.membership,
                           symmetric.b=symmetric.b,
                           strong.block=strong.block,

                           membership.a=membership.a))

            },

            gibbs.summary = function (gibbs.out) {
                membs1 <- {
                    d1 <- sapply(gibbs.out,
                                 function(gg){
                                     number.to.vector(gg$membership, nrow(gg$block.matrix))
                                 })
                    matrix(apply(d1, 1, mean),  ncol=n.nodes)
                }
                colnames(membs1) <- node.names
                this.block.matrix <- matrix(apply(sapply(gibbs.out,
                                                         function(gg){
                                                             c(gg$block.matrix)
                                                         }),
                                                  1, mean),
                                            nrow=n.groups)
                modal.membership <- apply(membs1, 2, which.max)
                return(list(membership=membs1,
                            modal.membership=modal.membership,
                            block.matrix=this.block.matrix))
            },
            print.gibbs.summary = function (gibbs.sum) {
                message ("Probabilistic block memberships:")
                print (gibbs.sum$membership)

                message ("Modal block memberships:")
                print (gibbs.sum$modal.membership)

                message ("Block value matrix:")
                print (gibbs.sum$block.matrix)

                return()
            },

            gibbs.node.order = function (gibbs.out) {
                get.sum <- gibbs.summary(gibbs.out)
            },

            gibbs.plot = function (gibbs.out, ...) {
                get.sum <- gibbs.summary(gibbs.out)
                block.membership.plot (get.sum$membership, get.sum$block.matrix,
                                       node.labels=node.names,
                                       main = "SBM Summary from Gibbs Sampler", ...)
            },

            gibbs.node.colors = function (gibbs.out, colors=(1:n.groups) + 1) {
                get.sum <- gibbs.summary(gibbs.out)
                return(colors[get.sum$modal.membership])
            }




        )
    )


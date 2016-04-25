#' Marxan reserve design with Gurobi
#'
#' Solve the Marxan reserve design optimization problem using Integer Quadratic
#' Programming (IQP) as implemented by the commercial optimization software
#' Gurobi.
#'
#' @details Marxan is the world's most widely used systematic conservation
#'   planning software. The optimization problem that Marxan solves is to find
#'   the set of planning units that minimizes the overall cost and boundary
#'   length of a reserve network, while meeting a set of representation targets
#'   for the conservation features. The cost is often either the area of the
#'   planning units or the opportunity cost of foregone commericial activities
#'   (e.g. loggin or agriculture). Minimizing the boundary length of the reserve
#'   network produces more compact, less fragmented spatial configurations.
#'   Finally, the representation targets ensure that each species is adequately
#'   represented in the reserve network.
#'
#'   Marxan solves this optimization problem using simulated annealing, a
#'   stochastic heuristic for approximating global optima of functions. However,
#'   the Marxan problem can be formualted as a Integer Quadratic Program (IQP)
#'   for which exact algorithms exist. This function uses the R interface to the
#'   commercial optimization software Gurobi to solve reserve design problems
#'   either exactly or to within some specified distance to optimality.
#'
#' @param pu pu RasterLayer or SpatialPolygonsDataFrame object; the planning
#'   units to use in the reserve design exercise and their corresponding cost.
#'   If a RasterLayer object is provided the cost should be stored in the cell
#'   values and if a SpatialPolygonsDataFrame object is provided it should have
#'   an attribute field named \code{cost}.
#' @param features RasterStack object; the distribution of conservation
#'   features. If \code{pu} is a RasterStack object then \code{features} should
#'   be defined on the same template. If \code{pu} is a SpatialPolygonsDataFrame
#'   \code{features} will be summarize over the polygons using
#'   \code{\link{summarize_features}}.
#' @param targets numeric; representation targets either as percent of total
#'   representation level when \code{target_type = "percent"} (the default), or
#'   absolute targets when \code{target_type = "absolute"}.
#' @param blm numeric; the boundary length modifier (BLM). The BLM sets the
#'   importance of minimizing the perimeter of the reserve network relative to
#'   minimizing the cost.
#' @param boundary matrix or data.frame; matrix of pairwise boundary lengths
#'   between planning units. Can be a sparse matrix from the slam package or a
#'   normal base matrix object. Alternatively, a data frame representation of
#'   this matrix with three variables: first planning unit index
#'   (\code{boundary$id1}), second planning unit index (\code{boundary$id2}),
#'   and corresponding boundary length (\code{boundary$boundary}). Note that
#'   unlike the corresponding Marxan input file (\code{boundary.dat}, which only
#'   has components with \code{id1 <= id2}, this data frame should be complete,
#'   i.e. have all components even it redundent. If this matrix is not provided
#'   it will be calculated based on the planning units defined by \code{pu}.
#' @param rij numeric matrix; matrix of representation level of conservation
#'   features (rows) within planning units (columns). Can be a sparse matrix
#'   from the slam package or a normal base matrix object. Alternatively, a data
#'   frame representation of this matrix with three variables: feature index
#'   (\code{rij$feature}), planning unit index (\code{rij$pu}), and
#'   corresponding representation level (\code{rij$amount}). If this matrix is
#'   not provided it will be calculated based on the planning units and
#'   RasterStack of conservation feature distributions.
#' @param edge_factor numeric; factor to scale lengths of edge boundaries that
#'   aren't shared with another planning unit. Set this to 0 if you don't wish
#'   to include these boundaries as diagonal components of the matrix. This
#'   scale factor can also be used in cases where edge boundaries are much
#'   longer than other shared boundaries, such as when edges are composed of
#'   complicated coast lines. In such situations, the long boundary lengths can
#'   bias selection against these edge planning units, and it may be desirable
#'   to scale the length to smaller values.
#' @param target_type character; whether the provided \code{targets} should be
#'   interpreted as percent-based (\code{"percent"}) or absolute
#'   (\code{"absolute"}).
#' @param gap numeric; proportional gap to optimality in returned solution. For
#'   example, a value of 0.005 (the default) will result in Gurobi stopping when
#'   the solution is within 5 percent of optimality.
#' @param time_limit numeric; time in seconds to run the Gurobi Optimizer. This
#'   will cut Gurobi off before is reaches an optimal solution.
#'
#' @return This function returns the results object created by Gurobi. This
#'   object is a named list, the most important elements of which are:
#'
#'   \itemize{
#'     \item{\code{objval}: The value of the objective function at the returned
#'     solution.}
#'     \item{\code{x}: Vector of decision variables for best solution found.}
#'   }
#' @export
#' @examples
#' # raster template
#' e <- raster::extent(0, 100, 0, 100)
#' r <- raster::raster(e, nrows = 100, ncols = 100, vals = 1)
#'
#' # generate 9 feature distributions with different rarities and scales
#' f <- mapply(function(x, y, r) gaussian_field(r = r, range = x, prop = y),
#'             rep(c(5, 15, 25), each = 3),
#'             rep(c(0.1, 0.25, 0.5), times = 3),
#'             MoreArgs = list(r = r))
#' f <- raster::stack(f)
#' f <- setNames(f, letters[1:raster::nlayers(f)])
#' # genrate cost layer
#' cost <- gaussian_field(r, 20, mean = 1000, variance = 500)
#' cost <- setNames(cost, "cost")
#' # solve to within 5% of optimality (the default)
#' results <- gurobi_marxan(cost, f, targets = 0.3, blm = 100)
#' # objective function value
#' results$objval
#'
#' # planning units can also be supplied as an SpatialPolygonsDataFrame object
#' # with cost stored as an attribute (pu$cost). Typically the function takes
#' # longer to execute with polygons because calculating boundaries and
#' # summarizing features over planning units is less efficient.
#' pu_spdf <- raster::rasterToPolygons(cost)
#' results_spdf <- gurobi_marxan(pu_spdf, f, targets = 0.3, blm = 100)
#' # cost
#' results_spdf$objval
gurobi_marxan <- function(pu, features, targets, blm, boundary, rij,
                          edge_factor, target_type, gap, time_limit) {
  UseMethod("gurobi_marxan")
}

#' @export
gurobi_marxan.Raster <- function(pu, features, targets, blm, boundary, rij,
                                 edge_factor = 1,
                                 target_type = c("percent", "absolute"),
                                 gap = 0.005, time_limit = Inf) {
  # assertions
  assert_that(requireNamespace("gurobi", quietly = TRUE),
              raster::nlayers(pu) == 1,
              assertthat::is.number(blm), blm >= 0,
              assertthat::is.number(edge_factor),
              assertthat::is.number(gap), gap >= 0,
              assertthat::is.number(time_limit), time_limit >= 0)

  # if blm = 0 this reduces to a linear minimum set coverage problem
  if (isTRUE(all.equal(blm, 0))) {
    if (missing(rij)) {
      results <- gurobi_minsetcover(pu = pu, features = features,
                                    targets = targets,
                                    target_type = target_type,
                                    gap = gap,
                                    time_limit = time_limit)
    } else {
      results <- gurobi_minsetcover(pu = pu, rij = rij,
                                    targets = targets,
                                    target_type = target_type,
                                    gap = gap,
                                    time_limit = time_limit)
    }
    return(results)
  }

  # representation matrix rij
  if (missing(rij)) {
    # if not provided, calculate it
    assert_that(inherits(features, "Raster"),
               raster::compareRaster(pu, features))
    rij <- slam::as.simple_triplet_matrix(t(unname(features[])))
  } else {
    # ensure that rij is a matrix, sparse matrix, or data frame
    assert_that(inherits(rij, c("matrix", "simple_triplet_matrix",
                                "data.frame")))
    if (is.matrix(rij)) {
      rij <- slam::as.simple_triplet_matrix(unname(rij))
    } else if (is.data.frame(rij)) {
      rij <- df_to_matrix(rij, nrow = raster::nlayers(features),
                          ncol = raster::ncell(pu),
                          vars = c("feature", "pu", "amount"))
    }
    # number of columns should be equal to number of planning units
    assert_that(rij$ncol == raster::ncell(pu))
  }

  # set proportional targets or check absolute targets
  target_type <- match.arg(target_type)
  if (target_type == "percent") {
    # convert percent-based targets to absolute targets
    targets <- set_targets(slam::row_sums(rij), targets)
  } else {
    # check that all targets are attainable
    assert_that(length(targets) == rij$nrow,
                all(targets <= slam::row_sums(rij)))
  }

  # planning unit boundaries
  if (missing(boundary)) {
    # calculate boundaries if not provided
    boundary <- calculate_boundaries(pu, edge_factor = edge_factor)
  } else {
    # ensure that boundary is a matrix or sparse matrix
    assert_that(inherits(boundary, c("matrix", "simple_triplet_matrix",
                                     "data.frame")))
    if (is.matrix(boundary)) {
      boundary <- slam::as.simple_triplet_matrix(unname(boundary))
    } else if (is.data.frame(boundary)) {
      boundary <- df_to_matrix(boundary,
                               nrow = raster::ncell(pu),
                               ncol = raster::ncell(pu),
                               vars = c("id1", "id2", "boundary"))
    }
    # number of columns and rows should be equal to number of planning units
    assert_that(boundary$nrow == raster::ncell(pu),
                boundary$ncol == raster::ncell(pu))
  }
  # matrix diagonal, external/edge boundaries
  d <- diag(as.matrix(boundary))
  # remove diagonal from main boundary matrix
  boundary <- boundary - slam::simple_triplet_diag_matrix(d)

  # linear component of objective function
  obj <- unname(pu[])
  obj <- obj + blm * slam::row_sums(boundary)
  # if edge_factor != 0 add external/edge boundaries to objective function
  if (!isTRUE(all.equal(edge_factor, 0))) {
    obj <- obj + blm * d
  }

  # construct model
  model <- list()
  # goal is to minimize objective function
  model$modelsense <- "min"
  # binary decision variables
  model$vtype <- "B"
  # objective function
  model$obj <- obj
  model$Q <- -blm * boundary
  # structural constraints
  model$A <- rij
  model$rhs <- unname(targets)
  model$sense <- rep(">=", length(targets))

  # set the parameters that control the algorithm
  # MIPgap defines the % gap to optimality
  params <- list(Presolve = 2, MIPGap = gap)
  if (is.finite(time_limit)) {
    params$TimeLimit = time_limit
  }

  # gurobi is very RAM hungry, remove any un-needed objects
  rm(pu, features, targets, obj, rij, boundary, d)

  # solve
  gurobi::gurobi(model, params)
}

#' @export
gurobi_marxan.SpatialPolygons <- function(pu, features, targets, blm,
                                      boundary, rij,
                                      edge_factor = 1,
                                      target_type = c("percent", "absolute"),
                                      gap = 0.005, time_limit = Inf) {

  # assertions
  assert_that(requireNamespace("gurobi", quietly = TRUE),
              "cost" %in% names(pu),
              assertthat::is.number(blm), blm >= 0,
              assertthat::is.number(edge_factor),
              assertthat::is.number(gap), gap >= 0,
              assertthat::is.number(time_limit), time_limit >= 0)

  # if blm = 0 this reduces to a linear minimum set coverage problem
  if (isTRUE(all.equal(blm, 0))) {
    if (missing(rij)) {
      results <- gurobi_minsetcover(pu = pu, features = features,
                                    targets = targets,
                                    target_type = target_type,
                                    gap = gap,
                                    time_limit = time_limit)
    } else {
      results <- gurobi_minsetcover(pu = pu, rij = rij,
                                    targets = targets,
                                    target_type = target_type,
                                    gap = gap,
                                    time_limit = time_limit)
    }
    return(results)
  }

  # representation matrix rij
  if (missing(rij)) {
    # if not provided, calculate it
    assert_that(inherits(features, "Raster"))
    rij <- summarize_features(pu, features)
  } else {
    # ensure that rij is a matrix or sparse matrix
    assert_that(inherits(rij, c("matrix", "simple_triplet_matrix",
                                "data.frame")))
    if (is.matrix(rij)) {
      rij <- slam::as.simple_triplet_matrix(unname(rij))
    } else if (is.data.frame(rij)) {
      rij <- df_to_matrix(rij, nrow = raster::nlayers(features),
                          ncol = length(pu),
                          vars = c("feature", "pu", "amount"))
    }
    # number of columns should be equal to number of planning units
    assert_that(rij$ncol == length(pu))
  }

  # set proportional targets or check absolute targets
  target_type <- match.arg(target_type)
  if (target_type == "percent") {
    # convert percent-based targets to absolute targets
    targets <- set_targets(slam::row_sums(rij), targets)
  } else {
    # check that all targets are attainable
    assert_that(length(targets) == rij$nrow,
                all(targets <= slam::row_sums(rij)))
  }

  # planning unit boundaries
  if (missing(boundary)) {
    # calculate boundaries if not provided
    boundary <- calculate_boundaries(pu, edge_factor = edge_factor)
  } else {
    # ensure that boundary is a matrix or sparse matrix
    assert_that(inherits(boundary, c("matrix", "simple_triplet_matrix",
                                     "data.frame")))
    if (is.matrix(boundary)) {
      boundary <- slam::as.simple_triplet_matrix(unname(boundary))
    } else if (is.data.frame(boundary)) {
      boundary <- df_to_matrix(boundary,
                               nrow = length(pu),
                               ncol = length(pu),
                               vars = c("id1", "id2", "boundary"))
    }
    # number of columns and rows should be equal to number of planning units
    assert_that(boundary$nrow == length(pu),
                boundary$ncol == length(pu))
  }
  # matrix diagonal, external/edge boundaries
  d <- diag(as.matrix(boundary))
  # remove diagonal from main boundary matrix
  boundary <- boundary - slam::simple_triplet_diag_matrix(d)

  # linear component of objective function
  obj <- pu$cost
  obj <- obj + blm * slam::row_sums(boundary)
  # if edge_factor != 0 add external/edge boundaries to objective function
  if (!isTRUE(all.equal(edge_factor, 0))) {
    obj <- obj + blm * d
  }

  # construct model
  model <- list()
  # goal is to minimize objective function
  model$modelsense <- "min"
  # binary decision variables
  model$vtype <- "B"
  # objective function
  model$obj <- obj
  model$Q <- -blm * boundary
  # structural constraints
  model$A <- rij
  model$rhs <- unname(targets)
  model$sense <- rep(">=", length(targets))

  # set the parameters that control the algorithm
  # MIPgap defines the % gap to optimality
  params <- list(Presolve = 2, MIPGap = gap)
  if (is.finite(time_limit)) {
    params$TimeLimit = time_limit
  }

  # gurobi is very RAM hungry, remove any un-needed objects
  rm(pu, features, targets, obj, rij, boundary, d)

  # solve
  gurobi::gurobi(model, params)
}

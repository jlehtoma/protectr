% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/gurobi-minsetcover.r
\name{gurobi_minsetcover}
\alias{gurobi_minsetcover}
\title{Minimum set cover reserve design with Gurobi}
\usage{
gurobi_minsetcover(pu, features, targets, rij, target_type, gap, time_limit)
}
\arguments{
\item{pu}{pu RasterLayer or SpatialPolygonsDataFrame object; the planning
units to use in the reserve design exercise and their corresponding cost.
If a RasterLayer object is provided the cost should be stored in the cell
values and if a SpatialPolygonsDataFrame object is provided it should have
an attribute field named \code{cost}.}

\item{features}{RasterStack object; the distribution of conservation
features. If \code{pu} is a RasterStack object then \code{features} should
be defined on the same template. If \code{pu} is a SpatialPolygonsDataFrame
\code{features} will be summarize over the polygons using
\code{\link{summarize_features}}.}

\item{targets}{numeric; representation targets either as percent of total
representation level when \code{target_type = "percent"} (the default), or
absolute targets when \code{target_type = "absolute"}.}

\item{rij}{numeric matrix; matrix of representation level of conservation
features (rows) within planning units (columns). Can be a sparse matrix
from the slam package or a normal base matrix object. Alternatively, a data
frame representation of this matrix with three variables: feature index
(\code{rij$feature}), planning unit index (\code{rij$pu}), and
corresponding representation level (\code{rij$amount}). If this matrix is
not provided it will be calculated based on the planning units and
RasterStack of conservation feature distributions.}

\item{target_type}{character; whether the provided \code{targets} should be
interpreted as percent-based (\code{"percent"}) or absolute
(\code{"absolute"}).}

\item{gap}{numeric; proportional gap to optimality in returned solution. For
example, a value of 0.005 (the default) will result in Gurobi stopping when
the solution is within 5 percent of optimality.}

\item{time_limit}{numeric; time in seconds to run the Gurobi Optimizer. This
will cut Gurobi off before is reaches an optimal solution.}
}
\value{
This function returns the results object created by Gurobi. This
  object is a named list, the most important elements of which are:

  \itemize{
    \item{\code{objval}: The value of the objective function at the returned
    solution.}
    \item{\code{x}: Vector of decision variables for best solution found.}
  }
}
\description{
Solve the minimum set cover reserve design optimization problem using Integer
Linear Programming (ILP) as implemented by the commercial optimization
software Gurobi.
}
\details{
In the context of systematic reserve design, the minimum set cover
  problem seeks to find the set of planning units that minimizes the overall
  cost of a reserve network, while meeting a set of representation targets
  for the conservation features. The cost is often either the area of the
  planning units or the opportunity cost of foregone commericial activities
  (e.g. loggin or agriculture). The representation targets ensure that each
  species is adequately represented in the reserve network. This is a
  simplified version of the Marxan objective function that doesn't account
  for the boundary length of the resulting reserve.

  Marxan solves this optimization problem using simulated annealing, a
  stochastic heuristic for approximating global optima of functions. However,
  this problem can be formualted as a Integer Linar Program (ILP) for which
  exact algorithms exist. This function uses the R interface to the
  commercial optimization software Gurobi to solve reserve design problems
  either exactly or to within some specified distance to optimality.
}
\examples{
# raster template
e <- raster::extent(0, 100, 0, 100)
r <- raster::raster(e, nrows = 100, ncols = 100, vals = 1)

# generate 9 feature distributions with different
f <- mapply(function(x, y, r) gaussian_field(r = r, range = x, prop = y),
            rep(c(5, 15, 25), each = 3),
            rep(c(0.1, 0.25, 0.5), times = 3),
            MoreArgs = list(r = r))
f <- raster::stack(f)
f <- setNames(f, letters[1:raster::nlayers(f)])
# genrate cost layer
cost <- gaussian_field(r, 20, mean = 1000, variance = 500)
cost <- setNames(cost, "cost")
# solve to within 5\% of optimality (the default)
results <- gurobi_minsetcover(cost, f, targets = 0.3)
# cost
results$objval

# planning units can also be supplied as an SpatialPolygonsDataFrame object
# with cost stored as an attribute (pu$cost). Typically the function takes
# longer to execute with polygons because summarizing features over planning
# units is less efficient.
pu_spdf <- raster::rasterToPolygons(cost)
results_spdf <- gurobi_minsetcover(pu_spdf, f, targets = 0.3)
# cost
results_spdf$objval
}


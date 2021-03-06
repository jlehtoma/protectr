% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/gurobi-maxcoverage.r
\name{gurobi_maxcoverage}
\alias{gurobi_maxcoverage}
\title{Maximum coverage reserve design with Gurobi}
\usage{
gurobi_maxcoverage(pu, features, budget, rij, gap, time_limit)
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

\item{budget}{numeric; budget for reserve.}

\item{rij}{numeric matrix; matrix of representation level of conservation
features (rows) within planning units (columns). Can be a sparse matrix
from the slam package or a normal base matrix object. Alternatively, a data
frame representation of this matrix with three variables: feature index
(\code{rij$feature}), planning unit index (\code{rij$pu}), and
corresponding representation level (\code{rij$amount}). If this matrix is
not provided it will be calculated based on the planning units and
RasterStack of conservation feature distributions.}

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
Solve the maximum coverage reserve design optimization problem using Integer
Linear Programming (ILP) as implemented by the commercial optimization
software Gurobi.
}
\details{
In the context of systematic reserve design, the maximum coverage
  problem seeks to find the set of planning units that minimizes the overall
  level of representation across a suite of conservation features, while
  keeping cost within a fixed budget. The cost is often either the area of
  the planning units or the opportunity cost of foregone commericial
  activities (e.g. loggin or agriculture). This is essentially the opposite
  of what the conservation planning software Marxan does.

  This problem can be formualted as a Integer Linar Program (ILP) for which
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
# set budget to 20\% of total cost of all cells
b <- 0.2 * raster::cellStats(cost, "sum")
# solve to within 5\% of optimality (the default)
results <- gurobi_maxcoverage(cost, f, budget = b)
# cost
results$objval

# planning units can also be supplied as an SpatialPolygonsDataFrame object
# with cost stored as an attribute (pu$cost). Typically the function takes
# longer to execute with polygons because summarizing features over planning
# units is less efficient.
pu_spdf <- raster::rasterToPolygons(cost)
results_spdf <- gurobi_maxcoverage(pu_spdf, f, budget = b)
# cost
results_spdf$objval
}


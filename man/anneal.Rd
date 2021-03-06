% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/anneal.r
\name{anneal}
\alias{anneal}
\title{Optimization via simulated annealing}
\usage{
anneal(x_init, objective, neighbour, t_init, alpha, n = 10000, ntemp = n/10,
  use_optim = TRUE)
}
\arguments{
\item{x_init}{vector; initial values for the decision variables.}

\item{objective}{function; the objective function that is to be optimized.
This function should only have a single argument: a vector of decision
variables. If any other data or parameters are required to calculate the
objective function, they should be must be stored within the enclosing
environment of the function. See the section on \bold{closures} for further
details.}

\item{neighbour}{function; neighbour selection function, i.e. given a set of
decision variable values this function should return a neighbouring state
to consider. This function should only have a single argument: a vector of
decision variables. If any other data or parameters are required to
calculate the objective function, they should be must be stored within the
enclosing environment of the function. See the section on \bold{closures}
forfurther details.}

\item{t_init}{numeric; starting temperature for simulated annealing.}

\item{alpha}{numeric; cooling factor to specify annealing schedule (see
Details). Ignored unless \code{use_optim = FALSE}.}

\item{n}{integer; the number of annealing iterations, which is also the
number of times the objective function will be evaluated.}

\item{ntemp}{integer; the number of temperature decreases. Must divide evenly
into the number of annealing iterations (\code{n}). See Details.}

\item{use_optim}{logical; whether to use \code{\link[stats]{optim}}, which
implements simulated annealing in C, or a pure R version implemented in
this package. The pure R version is significantly slower and is provided
mainly for educational purposes. Also, unlike \code{\link[stats]{optim}}
this pure R implementation uses the same annealing schedule as the
systematic conservation planning software Marxan (see Details).}
}
\value{
A vector specifying the state (i.e. values of the decision variables)
  that minimize the objective function. Since simulated annealing is
  stochastic, each time it is run a different set of values will be returned.
  Furthermore, there is no guarantee that the returned state corresponds to
  the global minimum, is may be a local minimum.
}
\description{
Minimize an objective function over a set of decision variables using
simulated annealing, a stochastic heuristic for approximating the global
optimum of a function.
}
\details{
This function can act as a wrapper for the C-based
function \code{\link[stats]{optim}} in the stats package (\code{use_optim =
TRUE}) or use a pure-R implementation of simulated annealing that closely
matches the one used by the popular conservation planning software Marxan
(see Details). The latter implementation is several order of magnitude
slower, and is included primarly for the purpose of emulating Marxan's
functionality exactly for educational purposes.

Simulated annealing is a flexible optimization heurstic, which can be applied
to complex functions with many local optima and functions of many discrete or
continuous variables. However, simulated annealing does not find the true
optimal solutions, rather it finds a near optimal solution and, since it is
stochastic, each time simulated annealing is run, a different near optimal
solution is returned.

Starting from some initial state, simulated annealing explores the state
space of the decision variables iteratively and stochastically by:

\enumerate{
  \item Randomly picking a \bold{neighbouring state} from a
  well-defined set of candidates, and calculating the value of the objective
  function at this new state.
  \item If this change reduces the objective
  function it is accepted. If the change increases the objective function, it
  is accepted according to an \bold{acceptance probability}, which depends on
  the magnitude of the increase in objective function and a global
  \bold{temperature} parameter.
  \item As the process progresses, the temperature parameter gradually
  decreases according to an \bold{annealing schedule}. As a result, changes
  that increase the objective function become less likely to be accepted
}

Initially having some non-zero acceptance probability for changes that
increase the objective function, reduces the likelihood of being caught in a
local minimum. This \bold{acceptance probability} is given by:

\deqn{min[1, exp(-\Delta/T_i)]}

where \eqn{\Delta} is the change in objective function and \eqn{T_i} is the
temperature after \eqn{i} decreases. The \bold{annealing schedule} determines
the rate at which temperature decreases, and alternate implementations of
simulated annealing use different functional forms for the schedule. In
particular, \code{\link[stats]{optim}} in the stats package uses

\deqn{T_i = \alpha^i * T_0}

while Marxan (and this function when \code{use_optim = FALSE}) uses

\deqn{T_i = T_0 / log(i + e)}

where \eqn{T_0} is the initial temperature and \eqn{\alpha} is a cooling
factor, both of which are user-defined parameters. These parameters should
be provided to this function via the \code{anneal_params} argument.
}
\section{Systematic conservation planning}{


  Simulated annealing is widely used to solve optimization problems in
  conservation planning, including by Marxan, the most widely used systematic
  reserve design tool. In conservation planning, knowning the exactly optimal
  reserve configuration is often not required, which makes simluated
  annealing appropriate. Instead, simulated annealing is performed many times
  to generate a suite of possible reserve networks. Thus, multiple options
  can be presented to decision makers, and the importance of particular
  planning units can be measured by the selection frequency (i.e. the number
  of solutions containing that planning unit).

  In the context of a Marxan-like systematic reserve design exercise, the
  decision variables are a set of binary values specifying whether each
  planning is selected or not. Marxan is exponentially faster than this
  simple R implementation of simulated annealing. However, if you wish to
  emulate Marxan use this function with \code{use_optim = FALSE} and use
  \code{\link{adaptive_schedule}} to set the initial temperature and cooling
  factor adaptively.
}

\section{Closures}{


  To keep this simulated annealing implementation as flexible as possible,
  the only argument the objective and neighbour selection functions take is
  the vector of decision variables. In general, other data or parameters will
  be required to evaluate these functions, but this additional information
  must be stored within the enclosing environment of the function. This can
  be accomplished using \bold{closures}, functions created by other
  functions, which have access to all variables in the parent function's
  environment. This will become more clear when I define these closure below.
  The Advanced R book by Hadley Wickham has great coverage of this topic,
  particularly the sections on function environments and closures. See the
  examples for a demonstration.
}
\examples{
# travelling salesman problem from optim documentation:
# what is the shortest distance one needs to travel to visit all
# cities in a set, given a matrix of distances between cities

# matrix of distances between cities from the stats package
euro_dist <- as.matrix(eurodist)
# the decision variable specifies the order the cities are visited in
# assume initial order is alphabetic
x_init <- c(1:nrow(euro_dist), 1)

# objective function generator
# distance matrix stored in enclosing environment
generate_ts_objective <- function(dist_mat) {
  function(x) {
  x <- embed(x, 2)
  sum(dist_mat[cbind(x[,2], x[,1])])
  }
}
objective <- generate_ts_objective(euro_dist)
objective(x_init)

# neighbour function generator, switches order of two cities
# distance matrix stored in enclosing environment
generate_ts_neighbour <- function(dist_mat) {
  function(x) {
    idx <- seq(2, nrow(dist_mat) - 1)
    change_points <- sample(idx, size = 2, replace = FALSE)
    tmp <- x[change_points[1]]
    x[change_points[1]] <- x[change_points[2]]
    x[change_points[2]] <- tmp
    x
  }
}
neighbour <- generate_ts_neighbour(euro_dist)
neighbour(x_init)

set.seed(123)
# solve using C-implementation in optim
x_optim <- anneal(x_init, objective, neighbour, n = 30000, t_init = 2000)
# shortest distance found
objective(x_optim)

# solve using pure R implementation
x_purer <- anneal(x_init, objective, neighbour, n = 30000,
                  t_init = 2000, alpha = 0.999, use_optim = FALSE)
# shortest distance found
objective(x_purer)
}


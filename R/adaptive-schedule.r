#' Adaptively choose parameters for annealing schedule
#'
#' In simulated annealing, the acceptance probability for a candidate solution
#' is dependent on a temperature parameter. As the heuristic progresses, this
#' temperature decreases according to an annealing schedule. In the Marxan
#' implementation of simulated annealing, the temperature decays exponentially
#' at a rate depending on a given initial and final temperature. The optimal
#' values for the annealing schedule parameters depend on the objective
#' function, and Marxan can chose these parameters adaptively. Following Marxan,
#' this function iteratively explores the phase space to find maximum and
#' minimum size of "bad" changes (i.e. those that increased the objective
#' function). The initial temperature is set to the maximum bad change and the
#' final temperature is set to 10% of the minimum bad change. The rate of
#' cooling can then be inferred from the initial and final temperatures.
#'
#' @details The cooling schedule used by Marxan is:
#'
#'   \deqn{T_i = \alpha^i * T_initial}
#'
#'   where \eqn{T_i} is the temperature after \eqn{i} decreases, \eqn{T_0} is
#'   the initial temperature, and \eqn{\alpha} is a cooling factor. \eqn{\alpha}
#'   is related to the initial and final temperatures by \code{exp(log(t_final /
#'   t_initial) / ntemp)}, where \code{ntemp} is the number of temperature
#'   decreases in the simulated annealing run.
#'
#' @param objective function; the objective function that is to be optimized.
#'   This function should only have a single argument: a vector of decision
#'   variables. If any other data or parameters are required to calculate the
#'   objective function, they should be must be stored within the enclosing
#'   environment of the function. See the section on \bold{closures} for further
#'   details.
#' @param neighbour function; neighbour selection function, i.e. given a set of
#'   decision variable values this function should return a neighbouring state
#'   to consider. This function should only have a single argument: a vector of
#'   decision variables. If any other data or parameters are required to
#'   calculate the objective function, they should be must be stored within the
#'   enclosing environment of the function. See the section on \bold{closures}
#'   forfurther details.
#' @param x_init vector; initial values for the decision variables.
#' @param n integer; number of iterations. Often chosen to be n / 100 where n is
#'   the numnber of annealing iterations.
#' @param ntemp integer; number of temperature decrases that will be used for
#'   the simulated annealing run. This parameter is required to calculate the
#'   cooling factor. By default it is equal to the argument \code{n}, however,
#'   it should be set to the number of temperature decrases that will be used
#'   for simulated annealing.
#'
#' @return Named numeric vector with two elements: the initial temperature
#'   (\code{t_initial}) and the final temperature (\code{t_final}). The cooling
#'   factor can then be calculated as \code{exp(log(t_final / t_initial) /
#'   ntemp)}, where \code{ntemp} is the number of temperature decreases in the
#'   simulated annealing run.
#' @export
#' @examples
#' NULL
adaptive_schedule <- function(objective, neighbour, x_init, n, ntemp = n) {
  # assertions
  assert_that(is.function(objective), is.function(neighbour),
              is.vector(x_init),
              assertthat::is.count(n),
              assertthat::is.count(ntemp))

  epsilon <- 1e-10
  scores <- numeric(n)
  x <- x_init
  scores[1] <- objective(x)
  for (i in 2:n) {
    x <- neighbour(x)
    scores[i] <- objective(x)
  }
  delta <- diff(scores)
  delta <- delta[delta > epsilon]
  # annealing schedule parameters
  t_initial <- max(delta)
  t_final <- 0.1 * min(delta)
  alpha <- exp(log(t_final / t_initial) / ntemp)
  c(t_initial = t_initial, t_final = t_final, alpha = alpha)
}

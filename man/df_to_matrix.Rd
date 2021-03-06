% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils.r
\name{df_to_matrix}
\alias{df_to_matrix}
\title{Convert data frame to matrix}
\usage{
df_to_matrix(x, nrow, ncol, sparse = TRUE, add_lower = FALSE,
  vars = c("i", "j", "v"))
}
\arguments{
\item{x}{data.frame; with variables \code{x$i} and \code{x$j} for the row and
column indices, respectively, and \code{x$v} for the matrix components.
Alternative names for these variables can be specified by \code{vars}.}

\item{nrow, ncol}{integer; the number of rows and columns, respectively.
Defaults to the maximum row and column indices, respectively.}

\item{sparse}{logical; whether to return a sparse matrix or a normal base
matrix object.}

\item{add_lower}{logical; if the data frame represents a symmetric matrix
where only the upper triangle is provided explicitly, set \code{add_lower}
to \code{TRUE} to add in these values explicitly to the matrix.}

\item{vars}{character vector; three element character vector with names of
row, column, and matrix component variables repsectively.}
}
\value{
A \code{\link[slam]{simple_triplet_matrix}} sparse matrix from the
  slam package, or a normal base matrix object.
}
\description{
Convert a data frame with row index, column index, and matrix components into
a corresponding matrix. By default returns a
\code{\link[slam]{simple_triplet_matrix}} sparse matrix from the slam
package.
}
\examples{
df <- data.frame(i = c(1, 1, 3), j = c(1, 3, 3), v = c(10, 20, 30))
df_to_matrix(df)
df_to_matrix(df, sparse = FALSE)
df_to_matrix(df, sparse = FALSE, add_lower = TRUE)

df <- data.frame(a = c(1, 1), b = c(1, 2), c = c(5, 25))
df_to_matrix(df, sparse = FALSE, vars = c("a", "b", "c"))
}


% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/hypertraps.R
\name{plotHypercube.influences}
\alias{plotHypercube.influences}
\title{Visualize pairwise influences (matrix)}
\usage{
plotHypercube.influences(
  my.post,
  featurenames = TRUE,
  use.regularised = FALSE,
  use.final = FALSE,
  reorder = FALSE,
  upper.right = FALSE,
  cv.thresh = Inf,
  red.green = FALSE
)
}
\arguments{
\item{my.post}{A model fit returned from HyperTraPS}

\item{featurenames}{Either logical (TRUE = use names from model fit), or a character vector of feature names. Default TRUE.}

\item{use.regularised}{Logical, whether to use the regularised (by parameter pruning) parameterisation. Defaults to FALSE.}

\item{use.final}{Logical, whether to use only the final parameterisation. Makes sense for point estimates (e.g. from simulated annealing), not for Bayesian posteriors. Defaults to FALSE.}

\item{reorder}{Logical, whether to reorder features by mean acquisition. Defaults to FALSE.}

\item{upper.right}{Logical, whether to arrange features in an upper-right direction. Defaults to FALSE.}

\item{cv.thresh}{Numeric. Upper threshold of CV for interactions to be plotted. Defaults to Inf.}

\item{red.green}{Logical. Adjust colors for red-green color-blindness. 
Defaults to FALSE.}
}
\value{
a ggplot
}
\description{
Plot pairwise influences between features in the L^2 picture, as a matrix.
}
\examples{
observations <- matrix(c(0,0,0,
                         0,0,1,
                         0,1,1,
                         1,1,1), byrow=TRUE, ncol=3)
fitted.cube <- HyperTraPS(observations)
plotHypercube.bubbles(fitted.cube)
}

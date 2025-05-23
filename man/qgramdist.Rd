% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/hypertraps.R
\name{qgramdist}
\alias{qgramdist}
\title{Calculate q-gram distance}
\usage{
qgramdist(my.post.1, my.post.2)
}
\arguments{
\item{my.post.1}{A fitted hypercube returned from [HyperTraPS()].}

\item{my.post.2}{A fitted hypercube returned from [HyperTraPS()].}
}
\value{
a list containing the distance between the two cubes.
}
\description{
Construct the (probability-weighted) q-gram distance between two hypercubes.
}
\examples{
observations.1 <- matrix(c(0,0,0,
                           0,0,1,
                           0,1,1,
                           1,1,1), byrow=TRUE, ncol=3)
fitted.cube.1 <- HyperTraPS(observations.1)
observations.2 <- matrix(c(0,0,0,
                           1,0,0,
                           1,1,0,
                           1,1,1), byrow=TRUE, ncol=3)
fitted.cube.2 <- HyperTraPS(observations.2)
qgramdist(fitted.cube.1, fitted.cube.2)
}

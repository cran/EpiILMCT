\name{SpatialData}
\docType{data}
\alias{SpatialData}
\title{Simulated epidemic data set from a distance-based SIR ILM}
\description{

This is a simulated epidemic data set of population size = 50 individuals that was generated using simple distance-based SIR individual-level model (ILM) with power-law distance kernel. The infectivity rate is given by:

\deqn{\lambda_{jt} = \alpha \sum_{i \in I_{t}}{d_{ij}^{-\beta}}}{\lambda_{jt} = \alpha \sum_{i in I_t}{d_{ij}^{-\beta}}}

The infectious period is assumed to follow exponential distribution with rate \eqn{\delta}. The epidemic was simulated with the following parameter values: \eqn{\alpha = 1.5}{\alpha = 1.5}, \eqn{\beta = 3} and \eqn{\delta = 2}. 

The data set file contains the XY coordinates of individuals, the Euclidean distance matrix, and the simulated epidemic data set.
}

\usage{SpatialData}
\format{
A list of the simulated data:    
   \describe{
  \item{loc}{
The XY coordinates of individuals.
  }	
  \item{dis}{
The Euclidean distance matrix.
  }	
  \item{epi}{
The simulated epidemic with four columns as: the id numbers of individuals, removal times, infectious periods, and infection times.
  }	
}
}
\keyword{datasets}
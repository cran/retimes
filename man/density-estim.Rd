\name{density-estim}
\Rdversion{1.1}
\alias{histestim}
\alias{distestim}
\alias{kernestim}

\title{Density Estimators}

\description{
Functions to identify the peak of a probability density distribution using a Gaussian kernel
estimator (\code{kernestim}), a fixed-width histogram estimator (\code{histestim}), and a
distance density estimator (\code{distestim}).
}

\usage{
kernestim(x, smoothing = NULL)
histestim(x, smoothing = NULL)
distestim(x, smoothing = NULL)
}

\arguments{
    \item{x}{Vector of data.}
    \item{smoothing}{
        Smothing parameter, also named bandwidth. If the argument is defined as \code{NULL}, the
        smoothing parameter will be calculated from the data (see details).
    }
}

\details{
The functions \code{kernestim} and \code{histestim} implement two well-known methods to estimate
the density of an empirical disribution (see Van Zandt, 2000, for their application in RT analysis).
The functions can be used to find the value corresponding to the peak of a empirical distribution.

The function \code{kernestim} cernter a Gaussian density over each observation and identifies
the value with the greater density.

The function \code{histestim} divides the data in bins, starting from the lower to the higher value
of data. The function searches the bin with the higher data frequency. The peak of the distribution
is identified calculating the mean of the data into the bin.

The function \code{distestim} is an experimental method which mixes the two previous techniques.
Around each data point (pivot), an interval \eqn{[x_{i}-h/2, x_{i}+h/2]} is builded, where \eqn{h}
is the smoothing parameter. The function searches the interval (bin) with the higher data frequency.
The output value is the weighted average of the values into the selected bin, in which each
observation is weighted on the basis of the distance from the pivot. If bins with equal densities
are found, the bin presenting the smallest deviance from the pivot is chosen.

For the Gaussian kernel estimator, the smoothing parameter is calculated using the Silverman's
method (Silverman, 1986). Differently, using histogram and distance estimators, the smoothing
paramete is calculated as: \eqn{(Q_{0.975}-Q_{0.025}) / \sqrt{n}}, where \eqn{Q_{p}} are the
quantiles for \eqn{\alpha = 0.05} and \eqn{n} is the sample dimension.
}

\value{The value at which the peak of density is located.}

\references{
Silverman, B. W. (1986).
Density estimation for statistics and data analysis.
London: Chapman & Hall.

Parzen, E. (1962).
On estimation of a probability density function and  mode.
Annals of Mathematical Statistics, 33(3), 1065-1076.

Van Zandt, T. (2000).
How to fit a response time distribution.
Psychonomic Bulletin & Review, 7(3), 424-465.
}

\author{Davide Massidda \email{davide.massidda@humandata.it}}

\examples{
x <- rexgauss(1000, mu=500, sigma=100, tau=250)
k <- kernestim(x); k
h <- histestim(x); h
d <- distestim(x); d
plot(density(x))
segments(k,0,k,1,col="red")
segments(h,0,h,1,col="blue")
segments(d,0,d,1,col="green")
}

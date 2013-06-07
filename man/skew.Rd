\name{skew}
\Rdversion{1.1}
\alias{skew}

\title{Skewness Index}

\description{Calculate the skewness index of an empirical distribution.}

\usage{skew(x)}

\arguments{\item{x}{Vector of data.}}

\details{
The skewness measures the asymmetry of a probability distribution. If a positive skewness
identifies a greater right tail, a negative skewness is symptom of a greater left tail of
the distribution.

The skewness is calculated as: \eqn{k_{3} / k_{2}^{1.5}}, where \eqn{k_{2}} and \eqn{k_{3}} are
respectively the second and the third moment of the distribution.
}

\value{The skewness index.}

\author{Davide Massidda \email{davide.massidda@humandata.it}}

\examples{
# Vector of reaction times from Heathcote (1996):
RT <- c(474.688, 506.445, 524.081, 530.672, 530.869,
        566.984, 582.311, 582.940, 603.574, 792.358)
skew(RT)
}

\name{exgauss}
\Rdversion{1.1}
\alias{exgauss}
\alias{dexgauss}
\alias{rexgauss}
\alias{mexgauss}

\title{Ex-Gaussian Distribution}

\description{Density function, parameter estimation by moments and random variable
generation for the ex-Gaussian distribution.}

\usage{
dexgauss(q, mu = 0, sigma = 1, tau = 1)
mexgauss(x, n = length(x))
rexgauss(n, mu = 0, sigma = 1, tau = 1, positive = TRUE)
}

\arguments{
    \item{x,q}{Vector of reaction times.}
    \item{n}{Number of observations.}
    \item{mu,sigma,tau}{Parameters of the ex-Gaussian distribution.}
    \item{positive}{Logical: imposes the generation of a vector containing only positive values.}
}

\details{
The function \code{mexgauss} implements the method of moments to estimate the parameters of an
ex-Gaussian distribution as described by Heathcote (1996).
}

\references{
Heathcote, A. (1996).
RTSYS: A DOS application for the analysis of reaction time data.
Behavior Research Methods, Instruments, & Computers, 28(3), 427-445.

Luce, R. D. (1986).
Response times: their role in inferring elementary mental organization.
New York: Oxford University Press.
}

\author{Davide Massidda \email{davide.massidda@humandata.it}}

\seealso{ \code{\link[retimes]{timefit}} }
\examples{
# Vector of reaction times from Heathcote (1996):
RT <- c(474.688, 506.445, 524.081, 530.672, 530.869,
        566.984, 582.311, 582.940, 603.574, 792.358)
mexgauss(RT)
}

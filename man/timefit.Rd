\name{timefit}
\alias{timefit}
\alias{timefit-class}
\alias{show,timefit-method}
\alias{plot,timefit-method}
\alias{logLik,timefit-method}
\alias{AIC,timefit-method}
\alias{BIC,timefit-method}
\Rdversion{1.1}

\title{Analyzing Reaction Times}
\description{
Parameter estimation for reaction time distributions using the maximum likelihood method.
Currently, the function works only for the ex-Gaussian probability distribution.
}

\usage{
timefit(x, iter = 0, size = length(x), replace = TRUE,
        plot = FALSE, start = NULL, ...)
}

\arguments{
    \item{x}{Vector of reaction times, in milliseconds or seconds.}
    \item{iter}{
        if \code{iter > 0}, the estimation is performed by bootstrap resampling. For each iteration,
        a sub-sample with a positive skewness will be estracted from the original sample, and
        distribution parameters will be estimated on resampled data.
    }
    \item{size}{Sample dimension to use in resampling, when the bootstrap is active.}
    \item{replace}{Logical: specifies if the resampling must be performed with replacement.}
    \item{plot}{
        Logical: if \code{TRUE}, it shows the data distribution (continuous line) and the
        estimated distribution (dashed line). Furthermore, when the bootstrap is used, the
        histogram of each parameter distribution is shown.
    }
    \item{start}{
        Vector containing starting values (in order: \eqn{mu}, \eqn{sigma} and \eqn{tau}) for
        the optimization routine. As default, the method of moments is used to estimate initial
        parameters. If the bootstrap is active, the argument will be ignored (starting values
        will be calculated on each bootstrapped sample using the method of moments).
    }
    \item{...}{
        Further arguments for the optimization routine. Arguments will be passed to the \code{optim}
        function: see the \code{optim} documentation for further details.
    }
}

\details{
    The function \code{timefit} estimates ex-Gaussian parameters by maximum likelihood, using the
    \code{optim} function. As default, the Simplex method (Nelder-Mead) is applied to find the
    minimum of the objective function.
    
    The function implements a bootstrap approach to identify distribution parameters. The bootstrap
    is useful especially working on small samples. The implemented method identifies \eqn{\mu} and
    \eqn{\sigma} using a gaussian kernel estimator (see Van Zandt, 2000, for details). Since small
    samples often present problems of convergence for \eqn{\sigma}, for this parameter only the
    values included a theoretical plausible range are considered. This range is defined by:
    
    \deqn{ \sqrt{\frac{min(x-M)^2}{n-1}} \leq \sigma \leq S }
    
    where M and S are respectively the mean and the standard deviation of data. When \eqn{\mu} and
    \eqn{\sigma} are identified, \eqn{\tau} is chosen within the bootstrapped values, according to
    the maximum likelihood criterion.
}

\note{The methods \code{plot}, \code{logLik}, \code{AIC} and \code{BIC} are available.}

\value{
    An object of class \code{timefit}. The operator \code{@} should be used to extract the values
    from the slots of the output object.
    \itemize{
    \item @x: vector of data.
    \item @samples: bootstrapped samples.
    \item @par: estimated parameters.
    \item @bootPar: bootstrapped parameters.
    \item @bootValue: log-likelihood of bootstrapped parameters.
    \item @sigmaValid: not rejected (\code{TRUE}) and rejected
                       (\code{FALSE}) values for \eqn{\sigma}.
    \item @start: starting values.
    \item @logLik: log-likelihood of estimated parameters.
    \item @AIC: Akaike Information Criterion.
    \item @BIC: Bayesian Information Criterion.
    }
}

\references{
Cousineau, D., Brown, S. & Heathcote, A. (2004).
Fitting distributions using maximum likelihood: Methods and packages.
Behavior Research Methods, Instruments, & Computers, 36(4), 742-756.

Heathcote, A. (1996).
RTSYS: A DOS application for the analysis of reaction time data.
Behavior Research Methods, Instruments, & Computers, 28(3), 427–445.

Lacouture, Y., & Cousineau, D. (2008).
How to use MATLAB to fit the ex-Gaussian and other probability functions
to a distribution of response times.
Tutorials in Quantitative Methods for Psychology, 4(1), 35-45.

Luce, R. D. (1986).
Response times: their role in inferring elementary mental organization.
New York: Oxford University Press.

Ratcliff, R. (1978).
A theory of memory retrieval. Psychological Review, 85(2), 59-108.

Ratcliff, R. (1979).
Group Reaction Time Distributions and an Analysis of Distribution Statistics.
Psychological Bulletin, 86(3), 446-461.

Ratcliff, R. & Murdock Jr., B. B. (1976).
Retrieval Processes in Recognition Memory.
Psychological Review, 83(3), 190-214.

Van Zandt, T. (2000).
How to fit a response time distribution.
Psychonomic Bulletin & Review, 7(3), 424-465.

Van Zandt, T.(2002).
Analysis of response time distributions.
In J. T. Wixted (Vol. Ed.) & H. Pashler (Series Ed.):
Stevens' Handbook of Experimental Psychology (3rd Edition),
Volume 4: Methodology in Experimental Psychology (pp. 461-516). New York: Wiley Press.
}

\author{Davide Massidda \email{davide.massidda@humandata.it}}

\seealso{
    \code{\link[retimes]{exgauss}},
    \code{\link[retimes]{density-estim}},
    \code{\link[retimes]{skew}}
}

\examples{\dontrun{

# Vector of reaction times from Heathcote (1996):
RT <- c(474.688, 506.445, 524.081, 530.672, 530.869,
        566.984, 582.311, 582.940, 603.574, 792.358)
RTfit <- timefit(RT); RTfit

# Boostrap tesing with simulated data (m=500, s=100, t=250):

x1 <- c(451.536,958.612,563.538,505.735,1266.825,
        860.663,457.707,268.679,587.303,669.594)
fit1a <- timefit(x1); fit1a
fit1b <- timefit(x1, iter=1000); fit1b

x2 <- c(532.474,525.185,1499.471,480.732,631.752,
        672.419,322.341,571.356,428.832,680.848)
fit2a <- timefit(x2, plot=TRUE); fit2a
fit2b <- timefit(x2, iter=1000, plot=TRUE); fit2b

}
}

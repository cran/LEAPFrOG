\name{constrOptim2}
\alias{constrOptim2}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
constrOptim2}
\description{
Called by LEAPFrOG and LEAPFrOG EM. Allows non-linearly constrained optimisation with hessian matrix returned
}
\usage{
constrOptim2(theta, f, grad, ui, ci, mu = 1e-04, control = list(), method = if (is.null(grad)) "Nelder-Mead" else "BFGS", outer.iterations = 100, outer.eps = 1e-05, hessian = FALSE, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{theta}{
%%     ~~Describe \code{theta} here~~
}
  \item{f}{
%%     ~~Describe \code{f} here~~
}
  \item{grad}{
%%     ~~Describe \code{grad} here~~
}
  \item{ui}{
%%     ~~Describe \code{ui} here~~
}
  \item{ci}{
%%     ~~Describe \code{ci} here~~
}
  \item{mu}{
%%     ~~Describe \code{mu} here~~
}
  \item{control}{
%%     ~~Describe \code{control} here~~
}
  \item{method}{
%%     ~~Describe \code{method} here~~
}
  \item{outer.iterations}{
%%     ~~Describe \code{outer.iterations} here~~
}
  \item{outer.eps}{
%%     ~~Describe \code{outer.eps} here~~
}
  \item{hessian}{
%%     ~~Describe \code{hessian} here~~
}
  \item{\dots}{
%%     ~~Describe \code{\dots} here~~
}
}
\author{
Ravi Varadhan
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (theta, f, grad, ui, ci, mu = 1e-04, control = list(), 
    method = if (is.null(grad)) "Nelder-Mead" else "BFGS", outer.iterations = 100, 
    outer.eps = 1e-05, hessian = FALSE, ...) 
{
    if (!is.null(control$fnscale) && control$fnscale < 0) 
        mu <- -mu
    R <- function(theta, theta.old, ...) {
        ui.theta <- ui \%*\% theta
        gi <- ui.theta - ci
        if (any(gi < 0)) 
            return(NaN)
        gi.old <- ui \%*\% theta.old - ci
        bar <- sum(gi.old * log(gi) - ui.theta)
        if (!is.finite(bar)) 
            bar <- -Inf
        f(theta, ...) - mu * bar
    }
    dR <- function(theta, theta.old, ...) {
        ui.theta <- ui \%*\% theta
        gi <- drop(ui.theta - ci)
        gi.old <- drop(ui \%*\% theta.old - ci)
        dbar <- colSums(ui * gi.old/gi - ui)
        grad(theta, ...) - mu * dbar
    }
    if (any(ui \%*\% theta - ci <= 0)) 
        stop("initial value not feasible")
    obj <- f(theta, ...)
    r <- R(theta, theta, ...)
    for (i in 1L:outer.iterations) {
        obj.old <- obj
        r.old <- r
        theta.old <- theta
        fun <- function(theta, ...) {
            R(theta, theta.old, ...)
        }
        gradient <- function(theta, ...) {
            dR(theta, theta.old, ...)
        }
        a <- optim(theta.old, fun, gradient, control = control, 
            method = method, hessian = hessian, ...)
        r <- a$value
        if (is.finite(r) && is.finite(r.old) && abs(r - r.old)/(outer.eps + 
            abs(r - r.old)) < outer.eps) 
            break
        theta <- a$par
        obj <- f(theta, ...)
        if (obj > obj.old * sign(mu)) 
            break
    }
    if (i == outer.iterations) {
        a$convergence <- 7
        a$message <- "Barrier algorithm ran out of iterations and did not converge"
    }
    if (mu > 0 && obj > obj.old) {
        a$convergence <- 11
        a$message <- paste("Objective function increased at outer iteration", 
            i)
    }
    if (mu < 0 && obj < obj.old) {
        a$convergence <- 11
        a$message <- paste("Objective function decreased at outer iteration", 
            i)
    }
    a$outer.iterations <- i
    a$barrier.value <- a$value
    a$value <- f(a$par, ...)
    a$barrier.value <- a$barrier.value - a$value
    a
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line


\name{best.logsym}
\alias{best.logsym}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Best log-symmetric model for the data
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
This function finds the best log-symmetric model for the data and presents the respective maximum likelihood estimates, as well as the
}
\usage{
best.logsym(y)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{y}{Vector of observations.
%%     ~~Describe \code{x} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
best.logsym gives the best log-symmetric model according to the Akaike information criterion (AIC). Supported families include BS, Normal, Student, Powerexp, Slash and Hyperbolic, which correspond to Birnbaum-Saunders, normal, Student-t, power-exponential, slash and symmetric hyperbolic distributions, respectively.
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
\item{eta }{Estimate of eta.}
\item{phi }{Estimate of phi.}
\item{xi }{Estimate of xi.}
\item{AIC }{Value of Akaike information criterion.}
\item{BIC }{Value of Bayesian information criterion.}
\item{family }{Family member obtained as the best model.}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
Leiva, V., Santos, R A., Saulo, H., and Marchant, C. (2020). Bootstrap control charts based on log-symmetric
distributions. Under review.
}
\author{Helton Saulo <heltonsaulo@gmail.com> and Rafael do Santos.
%%  ~~who you are~~
}
%\note{
%%  ~~further notes~~
%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

%\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
%}
\examples{
#====================================================
## data: Failure stresses (in GPa) for bundles of
## carbon fibres - mm150

## finds the best log-symmetric model
res150 <- best.logsym(y=mm150)
res150


#====================================================
## data: Failure stresses (in GPa) for bundles of
## carbon fibres - mm300

## finds the best log-symmetric model
res300 <- best.logsym(y=mm300)
res300
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line

\name{chart.logsym1}
\alias{chart.logsym1}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Log-symmetric bootstrap control chart - Phase I
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
This function provides the log-symmetric bootstrap control chart in phase I of control chart application.
}
\usage{
chart.logsym1(data, p, n, m, gamma, boot)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{Vector of observations.}
  \item{p}{Percentile of interest.}
  \item{n}{Number of samples.}
  \item{m}{Sample size of each sample.}
  \item{gamma}{False alarm rate (FAR).}
  \item{boot}{Number of bootstrap replications to construct the control chart.}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
chart.logsym1 constructs the log-symmetric bootstrap control chart in phase I of control chart application. The data are read as a vector and transformed into a matrix of n rows and m columns (the matrix is filled by rows). Supported families include BS, Normal, Student, Powerexp, Slash and Hyperbolic, which correspond to Birnbaum-Saunders, normal, Student-t, power-exponential, slash and symmetric hyperbolic distributions, respectively.
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
\item{family }{Member of the of log-symmetric family.}
\item{eta }{Estimate of eta.}
\item{phi }{Estimate of phi.}
\item{xi }{Estimate of xi.}
\item{lcl }{Lower control limit.}
\item{cl }{Center line.}
\item{ucl }{Upper control limit.}
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
## data: Table X of Qual. Reliab. Engng. Int. 2008, 24:585-600

## Phase I
resI <- chart.logsym1(data=WEInControl, p=0.01, n=10, m=5, gamma=0.0027, boot=10000)
resI

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line

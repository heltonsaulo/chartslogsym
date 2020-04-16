\name{chart.logsym.paper}
\alias{chart.logsym.paper}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Log-symmetric bootstrap control chart (paper)
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
This function provides the log-symmetric bootstrap control chart as shown in the application of Leiva et al. (2020).
}
\usage{
chart.logsym.paper(dataic, dataofc, paric, family, p, n, m, gamma, boot)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dataic}{In-control vector of observations.}
  \item{dataofc}{Out-of-control vector of observations.}
  \item{paric}{Vector of estimated parameters.}
  \item{family}{Log-symmetric family member.}
  \item{p}{Percentile of interest.}
  \item{n}{Number of samples.}
  \item{m}{Sample size of each sample.}
  \item{gamma}{False alarm rate (FAR).}
  \item{boot}{Number of bootstrap replications to construct the control chart.}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
chart.logsym.paper constructs the log-symmetric bootstrap control chart as shown in the application of Leiva et al. (2020). The data (dataic and dataofc) are read as vectors and transformed into matrices of n rows and m columns (the matrices are filled by rows). Supported families include BS, Normal, Student, Powerexp, Slash and Hyperbolic, which correspond to Birnbaum-Saunders, normal, Student-t, power-exponential, slash and symmetric hyperbolic distributions, respectively.

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
%% ...family = familyhat, eta = etahat, phi = phihat, xi = xihat, lcl = lcl, cl = cl, ucl = ucl)
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

## best model
res150 <- best.logsym(y=mm150)

## control chart
chart.logsym.paper(dataic=StudentInControl, dataofc=StudentOutofControl,
                   paric=c(res150$eta,res150$phi,res150$xi), family=res150$family,
                   p=0.01, n=20, m=5, gamma=0.01, boot= 10000)


#====================================================
## data: Failure stresses (in GPa) for bundles of
## carbon fibres - mm300

## best model
res300 <- best.logsym(y=mm300)

## control chart
chart.logsym.paper(dataic=SlashInControl, dataofc=SlashOutofControl,
                   paric=c(res300$eta,res300$phi,res300$xi), family=res300$family,
                   p=0.01, n=20, m=5, gamma=0.01, boot= 10000)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line

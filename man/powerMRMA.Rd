\name{powerMRMA}
\alias{powerMRMA}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{powerMRMA}
\description{
Power function to examine the mediated path from gene to outcome through an intermediate phenotype
}
\usage{
powerMRMA(plot.name = "powerMRMAplot", methodnames = c("MR.Classical", "MR.Egger", "MR.IVW", "MR.Median", "MA.Imai", "MA.4Way"), n = 1000, n.sim = 500, MeasurementError = T, alpha.level = 0.05, legend.include = T, color = T, seed = 1, nSNP = 4, MAF = c(0.2, 0.2, 0.2, 0.2), varM = 1, varY = 1, gamma0 = 0, gammaX = c(0.15, 0.15, 0.15, 0.15), gammaU = 0, beta0 = 0, betaM = c(0.15, 0.25), betaX = c(0, 0, 0, 0), betaI = c(0, 0, 0, 0), betaU = 0, muU = 0, varU = 1, muME = 0, varME = 1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{plot.name}{
%%     ~~Describe \code{plot.name} here~~
}
  \item{methodnames}{
%%     ~~Describe \code{methodnames} here~~
}
  \item{n}{
%%     ~~Describe \code{n} here~~
}
  \item{n.sim}{
%%     ~~Describe \code{n.sim} here~~
}
  \item{MeasurementError}{
%%     ~~Describe \code{MeasurementError} here~~
}
  \item{alpha.level}{
%%     ~~Describe \code{alpha.level} here~~
}
  \item{legend.include}{
%%     ~~Describe \code{legend.include} here~~
}
  \item{color}{
%%     ~~Describe \code{color} here~~
}
  \item{seed}{
%%     ~~Describe \code{seed} here~~
}
  \item{nSNP}{
%%     ~~Describe \code{nSNP} here~~
}
  \item{MAF}{
%%     ~~Describe \code{MAF} here~~
}
  \item{varM}{
%%     ~~Describe \code{varM} here~~
}
  \item{varY}{
%%     ~~Describe \code{varY} here~~
}
  \item{gamma0}{
%%     ~~Describe \code{gamma0} here~~
}
  \item{gammaX}{
%%     ~~Describe \code{gammaX} here~~
}
  \item{gammaU}{
%%     ~~Describe \code{gammaU} here~~
}
  \item{beta0}{
%%     ~~Describe \code{beta0} here~~
}
  \item{betaM}{
%%     ~~Describe \code{betaM} here~~
}
  \item{betaX}{
%%     ~~Describe \code{betaX} here~~
}
  \item{betaI}{
%%     ~~Describe \code{betaI} here~~
}
  \item{betaU}{
%%     ~~Describe \code{betaU} here~~
}
  \item{muU}{
%%     ~~Describe \code{muU} here~~
}
  \item{varU}{
%%     ~~Describe \code{varU} here~~
}
  \item{muME}{
%%     ~~Describe \code{muME} here~~
}
  \item{varME}{
%%     ~~Describe \code{varME} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ mediation }% use one of  RShowDoc("KEYWORDS")
\keyword{ power }% __ONLY ONE__ keyword per line

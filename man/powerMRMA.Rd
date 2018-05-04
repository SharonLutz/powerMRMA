\name{powerMRMA}
\alias{powerMRMA}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{powerMRMA}
\description{
Power function to examine the mediated path from gene to outcome through an intermediate phenotype
}
\usage{
powerMRMA(plot.name = "powerMRMAplot", MethodNames = c("MR.Classical", "MR.Egger", "MR.IVW", "MR.Median", "MA.Imai", "MA.4Way"), 
n = 1000, n.sim = 500, MeasurementError = T, alpha.level = 0.05, legend.include = T, color = T, SEED = 1, nSNP = 4, 
MAF = c(0.2, 0.2, 0.2, 0.2), varM = 1, varY = 1, gamma0 = 0, gammaX = c(0.15, 0.15, 0.15, 0.15), gammaU = 0, beta0 = 0, 
betaM = c(0,0.20, 0.25), betaX = c(0, 0, 0, 0), betaI = c(0, 0, 0, 0), betaU = 0, muU = 0, varU = 1, muME = 0, varME = 1)
}

\arguments{
  \item{plot.name}{ is the name of the plot
}
  \item{MethodNames}{ is the vector of methods to run
}
  \item{n}{ is the sample size, (default = 1000)
}
  \item{n.sim}{ is the number of simulations to run (default = 500)
}
  \item{MeasurementError}{
is generated on the mediator if T, is not generated if F
}
  \item{alpha.level}{
is the default alpha level
}
  \item{legend.include}{
is T for adding a legend to the plot, is F for not adding a legend
}
  \item{color}{
is T for producing a color plot, is F for producing a plot in grayscale
}
  \item{SEED}{
sets the seed for the function
}
  \item{nSNP}{
is the number of SNPS to generate
}
  \item{MAF}{
is the minor allele frequency
}
  \item{varM}{
is the variance of the mediator M
}
  \item{varY}{
is the variance of the outcome Y
}
  \item{gamma0}{
is the intercept for M
}
  \item{gammaX}{
is a vector of the associations of X with M
}
  \item{gammaU}{
is the unmeasured confounding on M
}
  \item{beta0}{
is the intercept for Y
}
  \item{betaM}{
is the association of M with Y 
}
  \item{betaX}{
is a vector of the direct effects of X on Y 
}
  \item{betaI}{
is a vector of the interaction effects between X and M on Y 
}
  \item{betaU}{
is the unmeasured confounding on Y 
}
  \item{muU}{
is the mean for U
}
  \item{varU}{
is the variance for U
}
  \item{muME}{
is the mean for ME
}
  \item{varME}{
is the varaince for ME
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
a matrix of the power of each method to detect an effect of the mediator M on the ouctome Y given that the first SNP is associated with the mediator (i.e. the indirect path) and corresponding plot
}
\references{
Thwing A, Ghosh D, Hokanson JE, Lutz SM. (2018) Mediated Paths 
in Genetic Association Studies: A Comparison of Mendelian Randomization 
and Mediation Analysis Approaches. (Target Journal).
}
\author{

}
\note{
%%  ~~further notes~~
}

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{

}

\keyword{ mediation }% use one of  RShowDoc("KEYWORDS")
\keyword{ power }% __ONLY ONE__ keyword per line

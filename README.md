# powerMRMA
Power package to examine the mediated path from gene to outcome through an intermediate phenotype

#### Installation
```
install.packages("devtools") # devtools must be installed first
install.packages("mediation")
install.packages("MendelianRandomization")

devtools::install_github("SharonLutz/powerMRMA")
```

#### Input

The desired methods to be run can be input using methodnames. The default is that all six methods will run, "MR.Classical", "MR.Egger", "MR.IVW", "MR.Median", "MA.Imai", and "MA.4Way".

The number of SNPs to be used as instrumental variables is set using nSNP, (default = 4). The MAF of these SNPS is set using MAF, and should be input as a vector (i.e. MAF = c(0.2,0.2,0.2,0.2)).

M is generated such that E\[M\] = ....

\Gamma0 is set with gamma0 (default = 0), \gammaX is set with gammaX and must be input as a vector whose length matches the number of SNPS (nSNP). \gammaU is set with gammaU (default = 0), a non-zero value of gammaU will generate unmeasured confounding. The variance of M is set with varM (default = 1).

Y is generated such that E\[Y\] = ....

\beta0 is set with beta0 (default = 0). \betaX is set with betaX, and indicates a direct effect from a SNP to the outcome,  it must be input as a vector whose length matches the number of SNPS (nSNP). The default is no direct effect from any SNP to the outcome, betaX = c(0,0,0,0). \betaU is set with betaU (default = 0), a non-zero value of betaU will generate unmeasured confounding. The association of M with Y is set using betaM. This must be input as a vector with length >= 2 and these values will be the x-axis of the output power plot. The variance of Y is set with varY (default = 1).

U is generated such that U = rnorm(\muU,\varU) where the mean can be set with muU (default = 0) and the variance can be set with varU (default = 1).

An interaction between a SNP_i and M on Y can be generated by setting index i of the betaI vector to be nonzero. BetaI must be a vector of length 4 (default is no interaction, betaI = c(0,0,0,0)).

Measurement error is generated such that M* = M + rnorm(\muME,\varME), and M* is used in analysis. muME is set using muME (default = 0) and varME is set using varME (default = 1). Measurement error can be included by setting MeasurementError to TRUE, setting MeasurementError to FALSE generates no measurement error of the mediator.

The name of the corresponding plot can be set using plot.name, the default is "powerMRMAplot". Legend.include can be set to true (default) if the user wants a legend included, or false if the user does not want the legend. color can be set to true (default) if the user wants the plot in color or false if the user wants the plot in gray scale.

The sample size can be set using n, the default sample size is 1000. The number of simulations can be set using n.sim, the default is 500. User should note that this function is slow. The seed can be set using seed (default = 1).

The default alpha level is 0.05, this can be changed using alpha.level. The mediation analysis and MR Classical methods which evaulate one SNP at a time are evaluated using an alpha level of alpha.level/ # of SNPs. The MR Egger, MR IVW and MR Median methods are evaluated using alpha.level.




#### Example
The code below evaluates the power of all six methods to detect an indirect effect from a SNP X to outcome Y given that a mediated path is simulated. This example detects the power of these six metohds four SNPs with a MAF of 0.2 are generated. The power is evaluated at two levels of association between M and Y, 0.15 and 0.25. There is no measurement error of the mediator or unmeasured confounding of the mediator outcome generated. There is no direct effect from any SNP X to the outcome Y or interaction between any X and M on Y generated. 

The classical approach to MR and mediation analysis methods detect the mediated pathway from one SNP to the outcome through the mediator M while the MR Egger, MR IVW and MR Median approach detect the mediated pathway from the same SNP to the outcome through using all four SNPS as instrumental variables for M. This code runs 100 simulations of a sample size of 1000.
```
library(powerMRMA)
?powerMRMA # For details on this function

powerMRMA (plot.name = "powerMRMAplot",methodnames = c("MR.Classical","MR.Egger","MR.IVW","MR.Median","MA.Imai","MA.4Way")
,n = 1000,n.sim=100,MeasurementError=F, nSNP = 4, MAF=c(0.2,0.2,0.2,0.2), gammaX = c(.15,.15,.15,.15), betaM = c(0.15,0.25))

```
#### Output
For this example analysis we get the following matrix of the power of each method to detect an indirect effect, and corresponding plot:

#### Reference


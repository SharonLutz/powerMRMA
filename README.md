# powerMRMA
Power package to examine the mediated path from gene to outcome through an intermediate phenotype. This R package compares Mendelian Randomization (MR) and mediation analysis approaches to detect the path from the mediator to the outcome given at least one SNP serves as an instrumental variable for the mediator.

#### Installation
```
install.packages("devtools") # devtools must be installed first
install.packages("mediation")
install.packages("MendelianRandomization")

devtools::install_github("SharonLutz/powerMRMA")
```

#### Input
nSNP is the number of SNPs generated from a binomial distribution for n subjects (input n) for a given minor allele frequency (input vector MAF).

For the SNPs Xi and unmeasured confounder U, the mediator M is generated from a normal distirbution with the variance (input varM) and the mean as follows:

E\[M\] = \gamma0 + sum \gammaX * Xi + \gammaU * U

All of these values are inputted by the user (i.e. the intercept gamma0, the genetic effect size as a vector gammaX, and the effect of the unmeasured confouder U as gammU).

The outcome Y is generated from a normal distribution with the variance (input varY) and the mean as follows:

E\[Y\] = beta0 + sum \betaX Xi + sum betaI* Xi* M + betaM * M + betaU* U 

All of these values are inputted by the user (i.e. the intercept beta0, the direct effect from the SNPs X to the outcome Y as a vector betaX, the effect of the unmeasured confouder U as gammU, the interaction effect of the SNPs with the mediator on the outcome Y as a vector betaI, and the effect of the mediator directly on the oucome as betaM).

The user can also generate measurment error. After the SNPs X, mediator M, and outcome Y are generated, then the powerMRMA package compare the power and type 1 error rate of the following 6 methods to detect the path from M to Y (i.e. betaM) given that at least one SNP serves as an instrumental variable for the mediator.

MethodNames denotes the possible methods used. (i.e. MethodNames= c("MR.Classical", "MR.Egger", "MR.IVW", "MR.Median", "MA.Imai", "MA.4Way")) where 

MR.Classical is the classical approach to MR.
Davey Smith, G., & Hemani, G. (2014). Mendelian randomization: genetic anchors for causal inference in epidemiological studies. Human Molecular Genetics, 23(R1), 89-98. 

MR.Egger is the Egger Regression approach to MR.
Bowden J., Davey Smith G., & Burgess S. (2015). Mendelian Randomization with invalid instruments: effect estimation and bias detection through Egger regression. International Journal of Epidemiology, 44(2), 512-525. 

MR.IVW is the Inverse Variant Weighted approach to MR.
Burgess, S., Butterworth, A., & Thompson, S. G. (2013). Mendelian Randomization Analysis With Multiple Genetic Variants Using Summarized Data. Genetic Epidemiology, 37(7), 658-665.

MR. Median is the Median Weighted approach to MR.
Bowden, J., Davey Smith, G., Haycock, P. C., & Burgess, S. (2016). Consistent Estimation in Mendelian Randomization with Some Invalid Instruments Using a Weighted Median Estimator. Genetic Epidemiology, 40(4), 304-314. 

MA.Imai is the Imai et al. approach to mediation analysis.
Imai, K., Keele, L., & Tingley, D. (2010). A general approach to causal mediation analysis. Psychological methods, 15(4), 309-334.

MA.4way is the 4 way decompoisition to mediation analysis.
VanderWeele, T. J. (2014). A unification of mediation and interaction: a four-way decomposition. Epidemiology (Cambridge, Mass.), 25(5), 749-761. 

Please use ?powerMRMA to see the man page which gives full details for all of the input parameters.

#### Example
This example displays the power of all six methods, using four SNPs as instrumental variables, each with a MAF of 0.2. Two levels of association between M and Y are evaluated, 0.15 and 0.25. There is no measurement error of the mediator or unmeasured confounding of the mediator outcome generated. There is no direct effect from any SNP X to the outcome Y or interaction between any X and M on Y generated. This code runs 100 simulations of a sample size of 1000.
```
library(powerMRMA)
?powerMRMA # For details on this function

powerMRMA (plot.name = "powerMRMAplot",methodnames = c("MR.Classical","MR.Egger","MR.IVW","MR.Median","MA.Imai","MA.4Way")
,n = 1000,n.sim=100,MeasurementError=F, nSNP = 4, MAF=c(0.2,0.2,0.2,0.2), gammaX = c(.15,.15,.15,.15), betaM = c(0, 0.2, 0.25))

```

#### Output
For this example analysis, we get the following matrix of the power of each method to detect an indirect effect, and corresponding plot:

```
     MR.Classical MR.Egger MR.IVW MR.Median MA.Imai MA.4Way
[1,]         0.00     0.00   0.00      0.00     0.0     0.0
[2,]         0.01     0.05   0.26      0.06     0.5     0.4
[1] "The plot is saved in your home directory."
```
<img src="https://github.com/SharonLutz/powerMRMA/blob/master/powerMRMAplot.png" width="600">

#### Reference
Thwing A, Ghosh D, Hokanson JE, Lutz SM. (2018) Mediated Paths in Genetic Association Studies: A Comparison of Mendelian Randomization and Mediation Analysis Approaches. (Target Journal).


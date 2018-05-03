# powerMRMA
Power package to examine the mediated path from gene to outcome through an intermediate phenotype. This R package compares Mendelian Randomization (MR) and mediation analysis approaches to detect the path from the mediator to the outcome given at least one SNP serves as an instrumental variable for the mediator.

# Installation
```
install.packages("devtools") # devtools must be installed first
install.packages("mediation")
install.packages("MendelianRandomization")

devtools::install_github("SharonLutz/powerMRMA")
```

# Input
nSNP is the number of SNPs generated from a binomial distribution for n subjects (input n) for a given minor allele frequency (input vector MAF).

For the SNPs Xi and unmeasured confounder U, the mediator M is generated from a normal distirbution with the variance (input varM) and the mean as follows:

```
E\[M\] = \gamma0 + sum \gammaX * Xi + \gammaU * U
```

All of these values are inputted by the user (i.e. the intercept gamma0, the genetic effect size as a vector gammaX, and the effect of the unmeasured confouder U as gammU).

The outcome Y is generated from a normal distribution with the variance (input varY) and the mean as follows:

```
E\[Y\] = beta0 + sum \betaX Xi + sum betaI* Xi* M + betaM * M + betaU* U 
```

All of these values are inputted by the user (i.e. the intercept beta0, the direct effect from the SNPs X to the outcome Y as a vector betaX, the effect of the unmeasured confouder U as gammU, the interaction effect of the SNPs with the mediator on the outcome Y as a vector betaI, and the effect of the mediator directly on the oucome as betaM).

After the SNPs X, mediator M, and outcome Y are generated, then the powerMRMA package compares the power and type 1 error rate of the following 6 methods to detect the path from M to Y (i.e. betaM) given that at least one SNP serves as an instrumental variable for the mediator. MethodNames denotes the possible methods used. (i.e. MethodNames= c("MR.Classical", "MR.Egger", "MR.IVW", "MR.Median", "MA.Imai", "MA.4Way")) where the citations are given at the end of this page.

The user can also generate measurment error. Please use ?powerMRMA to see the man page which gives full details for all of the input parameters.

# Example
Here, we consider 4 SNPs (nSNP=4) with MAF of 0.2. We vary the direct effect of the mediator M on the ouctome Y (i.e. betaM=c(0.15,0.25)). We simulate no measurement error of the mediator, no unmeasured confounding of the mediator-outcome relationship, no direct effect from any SNP X to the outcome Y, and no interaction between any SNP X and meditor M on the outcome Y. This code runs 200 simulations for n subjects with n=1000.

```
library(powerMRMA)
?powerMRMA # For details on this function

powerMRMA (plot.name = "powerMRMAplot",methodnames = c("MR.Classical","MR.Egger","MR.IVW","MR.Median","MA.Imai","MA.4Way")
,n = 1000,n.sim=100,MeasurementError=F, nSNP = 4, MAF=c(0.2,0.2,0.2,0.2), gammaX = c(.15,.15,.15,.15), betaM = c(0, 0.2, 0.25))

```

# Output
For this example, we get the following matrix of the type 1 error rate (row 1 with betaM=0) and power (row 2 with betaM=0.25 and row 3 with betaM=0.4) of each method to detect an effect of the mediator M on the ouctome Y given that the first SNP is associated with the mediator (i.e. the indirect path) and corresponding plot:

```
     MR.Classical MR.Egger MR.IVW MR.Median MA.Imai MA.4Way
[1,]         0.00     0.00   0.00      0.00     0.0     0.0
[2,]         0.01     0.05   0.26      0.06     0.5     0.4
[1] "The plot is saved in your home directory."
```
<img src="https://github.com/SharonLutz/powerMRMA/blob/master/powerMRMAplot.png" width="600">

# References
The power analysis used here is detailed in the following manuscript: <br/>
```
Thwing A, Ghosh D, Hokanson JE, Lutz SM. (2018) Mediated Paths 
in Genetic Association Studies: A Comparison of Mendelian Randomization 
and Mediation Analysis Approaches. (Target Journal).
```

MR.Classical is the classical approach to MR.<br/>
```
Davey Smith, G., & Hemani, G. (2014). Mendelian randomization: 
genetic anchors for causal inference in epidemiological studies. 
Human Molecular Genetics, 23(R1), 89-98. 
```

MR.Egger is the Egger Regression approach to MR.<br/>
```
Bowden J., Davey Smith G., & Burgess S. (2015). Mendelian Randomization 
with invalid instruments: effect estimation and bias detection through 
Egger regression. International Journal of Epidemiology, 44(2), 512-525. 
```

MR.IVW is the Inverse Variant Weighted approach to MR.<br/>
```
Burgess, S., Butterworth, A., & Thompson, S. G. (2013). Mendelian 
Randomization Analysis With Multiple Genetic Variants Using Summarized 
Data. Genetic Epidemiology, 37(7), 658-665.
```

MR. Median is the Median Weighted approach to MR.<br/>
```
Bowden, J., Davey Smith, G., Haycock, P. C., & Burgess, S. (2016). Consistent 
Estimation in Mendelian Randomization with Some Invalid Instruments Using a 
Weighted Median Estimator. Genetic Epidemiology, 40(4), 304-314. 
```

MA.Imai is the Imai et al. approach to mediation analysis.<br/>
```
Imai, K., Keele, L., & Tingley, D. (2010). A general approach to causal mediation 
analysis. Psychological methods, 15(4), 309-334.
```

MA.4way is the 4 way decompoisition to mediation analysis.<br/>
```
VanderWeele, T. J. (2014). A unification of mediation and interaction: a 
four-way decomposition. Epidemiology (Cambridge, Mass.), 25(5), 749-761. 
```

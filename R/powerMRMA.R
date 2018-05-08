powerMRMA <-function(plot.name = "powerMRMAplot",MethodNames = c("MR.Classical","MR.Egger","MR.IVW","MR.Median","MA.Imai","MA.4Way"),n = 1000,n.sim=500,
                      MeasurementError=T, alpha.level=0.05, legend.include=T, color=T, SEED=1, nSNP = 4, MAF=c(0.2,0.2,0.2,0.2), varM = 1, varY = 1, gamma0 = 0,
                      gammaX = c(.15,.15,.15,.15), gammaU = 0, beta0 = 0, betaM = c(0,0.2,0.25), betaX = c(0,0,0,0), betaI = c(0,0,0,0), betaU = 0, 
                      muU = 0, varU = 1, muME = 0, varME = 1){
  
  library(MendelianRandomization)
  library(mediation)
  
  set.seed(SEED)
  
  possibleMethods <- c("MR.Classical","MR.Egger","MR.IVW","MR.Median","MA.Imai","MA.4Way")
  
  # Length of methods to be used for color/ line type/ symbol vectors
  numMethods <- length(MethodNames)
  
  colorVec <- c("black" ,"red" , "green3" , "blue" , "cyan" , "magenta")
  
  # Set grayscale if color = F, set to color if color = T
  if(!typeof(color)=="logical"){stop("Error: color must be T/True or F/False")}
  if(color){
    colorVec <- c(colorVec[1:numMethods])
  }else{
    colorVec <- gray.colors(numMethods)
  }
  
  pchv <- c(1:numMethods)
  ltyv <- rep_len(c(2:numMethods),length.out = length(MethodNames))
  
  # Check that the methods are entered correctly
  for(mm in 1:length(MethodNames)){
    if(!toupper(MethodNames[mm])%in%toupper(possibleMethods)){
      stop(paste("Error: ",MethodNames[mm]," is not a valid method name. See ?powerMRMA.",sep=""))
    }
  }
  
  ### Reorder and rename methods to be same order every time
  tempMethods <- rep(NA,length(possibleMethods))
  
  for(jj in 1:length(MethodNames)){
    if(toupper(MethodNames[jj])=="MR.CLASSICAL"){tempMethods[1]<-"MR.Classical"}
    if(toupper(MethodNames[jj])=="MR.EGGER"){tempMethods[2]<-"MR.Egger"}
    if(toupper(MethodNames[jj])=="MR.IVW" ){tempMethods[3]<-"MR.IVW"}
    if(toupper(MethodNames[jj])=="MR.MEDIAN" ){tempMethods[4]<-"MR.Median"}
    if(toupper(MethodNames[jj])=="MA.IMAI"){tempMethods[5]<-"MA.Imai"}
    if(toupper(MethodNames[jj])=="MA.4WAY"){tempMethods[6]<-"MA.4Way"}
  }
  
  MethodNames <- tempMethods[!is.na(tempMethods)]
  ####
  
  # Error Checks for type/ length
  if(nSNP!=length(MAF)){stop("Error: nSNP must equal length(MAF).")}
  if(nSNP!=length(gammaX)){stop("Error: nSNP must equal length(gammaX).")}
  if(nSNP!=length(betaX)){stop("Error: nSNP must equal length(betaX).")}
  if(nSNP!=length(betaI)){stop("Error: nSNP must equal length(betaI).")}
  if(length(muU)>1){stop("Error: muU must not be a vector")}
  if(length(varU)>1){stop("Error: varU must not be a vector")}
  if(length(varM)>1){stop("Error: varM must not be a vector")}
  if(length(varY)>1){stop("Error: varY must not be a vector")}
  if(length(muME)>1){stop("Error: muME must not be a vector")}
  if(length(varME)>1){stop("Error: varME must not be a vector")}
  if(length(unique(betaM))<2){stop("Error: there must be at least two distinct values of betaM")}
  if(!typeof(MeasurementError)=="logical"){stop("Error: MeasurementError must be T/True or F/False")}
  if(!typeof(legend.include)=="logical"){stop("Error: legend.include must be T/True or F/False")}
  
  # Error checks for possible values
  if(!(varU>0)){stop("Error: varU must be > 0")}
  if(!(varM>0)){stop("Error: varM must be > 0")}
  if(!(varY>0)){stop("Error: varY must be > 0")}
  if(!(varME>0)){stop("Error: varME must be > 0")}
  if(min(MAF)<0 | max(MAF)>1){stop("Error: all MAF values must be between 0 and 1")}
  if(min(alpha.level)<0 | max(alpha.level)>1){stop("Error: alpha level must be between 0 and 1")}
  if(floor(n)!=ceiling(n)){stop("Error: n must be an integer.")}
  if(n<0 | n==0){stop("Error: n must be greater than zero.")}
  if(floor(n.sim)!=ceiling(n.sim)){stop("Error: n.sim must be an integer.")}
  if(n.sim<0 | n.sim==0){stop("Error: n.sim must be greater than zero.")}
  
  # Create Results Matrix
  mat_results <- matrix(0,nrow=length(betaM),ncol=6)
  colnames(mat_results) <- c("MR.Classical","MR.Egger","MR.IVW","MR.Median","MA.Imai","MA.4Way")
  
  for(ind in 1:n.sim){
    print(paste("sim ",ind," of ",n.sim,sep = ""))
    
    mat <- matrix(0,nrow=length(betaM),ncol=6)
    colnames(mat) <- c("MR.Classical","MR.Egger","MR.IVW","MR.Median","MA.Imai","MA.4Way")
    
    for(betaM.ind in 1:length(betaM)){
      ## Data Generation
      matX <- matrix(0,nrow=n,ncol=nSNP)
      for(xx in 1:nSNP){
        matX[,xx] <- rbinom(n,2,MAF[xx]) # Generate x, where column xx of matX has MAF from index xx of MAF.vec
      }
      
      U <- rnorm(n,muU,varU) # Generate U
      
      gammaX.total <- matX %*% gammaX # Create gammaX as the sum of gammaXi*Xi for each i in 1:n
      muM <- gamma0 + gammaX.total + gammaU*U # Set the seed and generate, make sure the rows Y2 <- rnorm(2)

      M <- rnorm(n,muM,varM) # Generate M with mean muM and variance varM
      
      betaX.total <- matX %*% betaX
      betaI.total <- matX %*% betaI
      
      muY <- beta0 + betaX.total + betaM[betaM.ind]*M + betaI.total + betaU*U # 
      
      Y <- rnorm(n,muY,varY) # Generate Y with mean muY and variance varY
      
      if(MeasurementError){M <- M + rnorm(n,muME,varME)} # Add measurement error
      
      ################# Run the methods
      
      betaXMR <- c()
      betaXseMR <- c()
      betaXMRy <- c()
      betaXseMRy <- c()
      # Mendelian Randomization
      for(snp.ind in 1:ncol(matX)){
        fit <- lm(M~matX[,snp.ind])
        # Add the beta estimate to a betaXMR vector
        betaXMR <- c(betaXMR,summary(fit)$coefficients[2,1])
        # Add the SE to a betaXMRse vector
        betaXseMR <- c(betaXseMR,summary(fit)$coefficients[2,2])
        
        # betaY are the beta-coefficients from regression analyses of the outcome on each genetic
        # variant in turn, and betaYse are the standard errors
        fitY <- lm(Y~matX[,snp.ind])
        betaXMRy <- c(betaXMRy,summary(fitY)$coefficients[2,1])
        betaXseMRy <- c(betaXseMRy,summary(fitY)$coefficients[2,2])
      }
      
      if("MR.Classical" %in% MethodNames){
        fitY <- lm(Y~matX[,1])
        if(summary(fitY)$coefficients[2,4] < (alpha.level/nSNP)){ mat[betaM.ind,"MR.Classical"] <- mat[betaM.ind,"MR.Classical"] + 1 }
      }
      
      # Variants aren't correlated so don't need to include a correlation matrix
      mr.input <- mr_input(bx = betaXMR, bxse = betaXseMR, by = betaXMRy, byse = betaXseMRy)
      # Run all MR methods
      if("MR.Median" %in% MethodNames){ 
        mr.median.p <- mr_median(mr.input,seed = NA)$Pvalue
        if(mr.median.p < alpha.level){ mat[betaM.ind,"MR.Median"] <- mat[betaM.ind,"MR.Median"] + 1 }
      }
      if("MR.Egger" %in% MethodNames){
        mr.egger.p <- mr_egger(mr.input,robust = F,penalized = F)$Pvalue.Est
        if(mr.egger.p < alpha.level){ mat[betaM.ind,"MR.Egger"] <- mat[betaM.ind,"MR.Egger"] + 1 }
      }
      if("MR.IVW" %in% MethodNames){
        mr.ivw.p <- mr_ivw(mr.input,robust = F,penalized = F)$Pvalue
        if(mr.ivw.p < alpha.level){ mat[betaM.ind,"MR.IVW"] <- mat[betaM.ind,"MR.IVW"] + 1 }
      }
      
      if("MA.Imai" %in% MethodNames){
        nSimImai <- 1000
        SNPinterest <- matX[,snp.ind]
        # Interaction between x1 and M on Y
        med.fit <- lm(M ~ SNPinterest)
        out.fit <- lm(Y ~ M + SNPinterest + M*SNPinterest)
        med.out <- mediate(med.fit, out.fit, treat = "SNPinterest", mediator = "M", sims = nSimImai, conf.level = 1-(alpha.level/nSNP))
        pindirect <- summary(med.out)$d.avg.p # Gets the p-value (of ACME(average))
        if(pindirect<(alpha.level/nSNP)){mat[betaM.ind,"MA.Imai"] <- mat[betaM.ind,"MA.Imai"]+1}
      }
      
      if("MA.4Way" %in% MethodNames){
        SNPinterest <- matX[,snp.ind]
        medV.out <- mediationV(dataS=cbind(SNPinterest,M,Y),yvar="Y", mvar="M", avar="SNPinterest", yreg="linear", mreg="linear",int=TRUE,output="full",boot=T,alpha_level = alpha.level/nSNP)
        medV.out <- as.matrix(medV.out)
        if((as.numeric(medV.out[3,5])<(alpha.level/nSNP))&(!is.na(medV.out[3,5]))){
          mat[betaM.ind,"MA.4Way"]<-mat[betaM.ind,"MA.4Way"]+1}
      }
      
    } # betaM loop end
    
    # Add the current matrix (mat) to the ongoing matrix (mat_results)
    mat_results <- mat_results + mat
    
  } # End of nSim loop
  
  # All of the simulations have been run, divide the matrix by the number of simulations run
  mat_results <- mat_results/n.sim
  
  # Output the table for the user
  print(mat_results)
  
  # Plot the results
  pdf(paste("~/",plot.name,".pdf",sep=""))
  
  plot(-1,-1,xlim=c(min(betaM),max(betaM)),ylim=c(0,1),ylab="Power",xlab="beta M")
  
    for(mp in 1:length(MethodNames)){
      if(MethodNames[mp]=="MR.Classical"){lines(betaM,mat_results[,"MR.Classical"],pch = pchv[mp],col = colorVec[mp],type="b",lty = ltyv[mp])}
      if(MethodNames[mp]=="MR.Egger"){lines(betaM,mat_results[,"MR.Egger"],pch = pchv[mp],col = colorVec[mp],type="b",lty = ltyv[mp])}
      if(MethodNames[mp]=="MR.IVW"){lines(betaM,mat_results[,"MR.IVW"],pch = pchv[mp],col = colorVec[mp],type="b",lty = ltyv[mp])}
      if(MethodNames[mp]=="MR.Median"){lines(betaM,mat_results[,"MR.Median"],pch = pchv[mp],col = colorVec[mp],type="b",lty = ltyv[mp])}
      if(MethodNames[mp]=="MA.Imai"){lines(betaM,mat_results[,"MA.Imai"],pch = pchv[mp],col = colorVec[mp],type="b",lty = ltyv[mp])}
      if(MethodNames[mp]=="MA.4Way"){lines(betaM,mat_results[,"MA.4Way"],pch = pchv[mp],col=colorVec[mp],type="b",lty = ltyv[mp])}
    }
  
  # Create a legend only out of the included methods
  leg.text<-rep(0,length(MethodNames))
  for(ff in 1:length(MethodNames)){
    if(MethodNames[ff]=="MR.Classical"){leg.text[ff]<-c("MR.Classical")}
    if(MethodNames[ff]=="MR.Egger"){leg.text[ff]<-c("MR.Egger")}
    if(MethodNames[ff]=="MR.IVW"){leg.text[ff]<-c("MR.IVW")}
    if(MethodNames[ff]=="MR.Median"){leg.text[ff]<-c("MR.Median")}
    if(MethodNames[ff]=="MA.Imai"){leg.text[ff]<-c("MA.Imai")}
    if(MethodNames[ff]=="MA.4Way"){leg.text[ff]<-c("MA.4Way")}
  }
  
  if(legend.include==TRUE){legend("topleft",legend = leg.text,col=colorVec,lty=ltyv,pch=pchv)}
  
  dev.off()
  
  print("The plot is saved in your home directory.")
  
}

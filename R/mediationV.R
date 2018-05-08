mediationV <-
function(dataS,yvar="y",avar="x",mvar="m",cvar='',a0=0,a1=1,m=0,nc=0,yreg="linear",mreg="linear",
           int=FALSE,casecontrol=FALSE,output="full",c='',boot=TRUE,alpha_level=0.05){
    numboot = 1000
    alphalev = alpha_level
    dataS = as.data.frame(dataS)
    na.matrix <- NULL
    # Data House Keeping
    if(cvar[1]!=''){data1 <- cbind(dataS[avar],dataS[mvar],dataS[yvar],dataS[cvar])}
    if(cvar[1]==''){data1 <- cbind(dataS[avar],dataS[mvar],dataS[yvar])}
    # Adds an int column that is avar times mvar
    if(int==TRUE){ 
      data1["int"] <- dataS[avar]*dataS[mvar]
      na.matrix <- matrix(NA,6,5,dimnames = list(c("cde","pnde","pnie","tnde","tnie","total effect"),
                                                 c("Estimate","Standard Error", "95% CI lower", "95% CI upper", "p-value")))
    }
    if(int!=TRUE){
      na.matrix <- matrix(NA,3,5,dimnames = list(c("cde=nde","nie","total effect"),
                                                 c("Estimate","Standard Error", "95% CI lower", "95% CI upper", "p-value")))
    }
    
    # Might need to move this to below data housekeeping section:
    numCoefficients = ncol(data1)
    # print(numCoefficients)
    
    # print("Made it here")
    if(cvar[1]!='' & casecontrol==FALSE | cvar[1]!='' & casecontrol==''){
      cvars=cvar; i=1
      data2 <- c()
      colnames <- c()
      meanVector <- c()
      for(i in 1:length(cvars)){
        covariateIndex <- cvars[i]
        currentMean <- mean(unlist(data1[covariateIndex]))
        # Assign current mean to data2i
        assign(paste("data2",i,sep=""),currentMean)
        meanName <- paste(cvars[i],"_Mean",sep="")
        meanVector <- c(meanVector,meanName)
        assign(paste("data2new",i,sep=""),currentMean)
        data2 <- c(data2,get(paste("data2new",i,sep="")))
        i=i+1
      }
      vb <- as.data.frame(data2)
      data2 <- t(vb)
      
      colnames(data2) <- meanVector
      rownames(data2) <- NULL
      if(c!=''){
        cval <- c ; i <- 1
        for(m in length(cval)){
          newname <- paste("cval",m)
          data2[newname] <- mean(data2[m])
          m = m+1
        }
      }
    } # End of if(cvar!='' & casecontrol==FALSE | cvar!='' & casecontrol=='')
    # print(data2)
    
    if(cvar[1] != '' & casecontrol==TRUE){
      cvars=cvar; i=1
      data3 <- c()
      for(i in 1:length(cvars)){
        nums <- which(data1[yvar]==0)
        newset <- c()
        for(s in nums){
          newset <- rbind(newset,data1[s,1:ncol(data1)])
        }
        # Create a data frame of means
        means <- c()
        for(s in 1:length(newset)){
          tempmean <- mean(unlist(newset[s]))
          means <- c(means, tempmean)
        } 
        assign(paste("data2",i,sep=""),means[i])
        temp <- get(paste("data2",i,sep=""))
        temp <- as.matrix(temp)
        colnames(temp) <- c("mean")
        rownames(temp) <- NULL
        assign(paste("data2new",i,sep=""),temp)
        data3 <- c(data3,get(paste("data2new",i,sep="")))
        i=i+1
      }
      vb <- as.data.frame(data3)
      data3 <- t(vb)
      data2 <- data3
      colnames(data2) <- rep("mean",length(data2))
      rownames(data2) <- NULL
      if(c!=''){
        cval <- c ; i <- 1
        for(m in length(cval)){
          newname <- paste("cval",m)
          data2[newname] <- mean(data2[m])
          m = m+1
        }
      }
    }  # End of if(cvar != '' & casecontrol=TRUE)
    
    # Bootstrap Procedure
    #***************   BOOTSTRAP PROCEDURE   ***************
    if (boot!='' & boot!= FALSE){
      if(boot==TRUE){nboot=numboot}
      if(boot!=TRUE){nboot=boot}
      
      #***************** bootstrap samples   ******************
      # To create b bootstrap replications
      numrows <- nrow(data1)
      total = nboot*numrows
      sampleCol <- numeric(total)
      iCol <- numeric(total)
      points <- matrix(nrow=total,ncol=ncol(data1))
      for(samp in 1:nboot){
        temp <- data1[sample(1:numrows,replace = T),]
        temp <- as.matrix(temp)
        points[(((samp-1)*numrows)+1):(samp*numrows),]=temp
      } # end of sample 1 to nboot
      rownames(points) <- NULL
      iCol <- rep(seq_len(numrows),nboot)
      sampleCol <- rep(1:nboot,each=numrows)
      newData <- cbind(sampleCol,iCol,points)
      newData <- as.data.frame(newData)
      
      if(cvar[1]==''){
        if(int){colnames(newData) <- c("sample","i",avar,mvar,yvar,"int")}
        if(!int){colnames(newData) <- c("sample","i",avar,mvar,yvar)}
      }
      if(cvar[1]!=''){
        if(int){colnames(newData) <- c("sample","i",avar,mvar,yvar,cvar,"int")}
        if(!int){colnames(newData) <- c("sample","i",avar,mvar,yvar,cvar)}
      }
      data1 <- newData
      # Sets data1t to be dataset containing all of the columns of sample t
      var=0
      for(t in 1:nboot){
        temp <- data1[(1+nrow(dataS)*var):(nrow(dataS)+nrow(dataS)*var),1:ncol(data1)]
        assign(paste("data1",t,sep=""),temp)
        var = var + 1
      } # end of do t=1:nboot
    } # End of if (boot!='' & boot!= FALSE)
    
    #***************** regression-for bootstrap ************************
    if(boot!='' & boot!=FALSE){
      for(t in 1:nboot){
        #print(paste("t2",t,"in",nboot))
        #***************************************************************
        if(yreg=="linear"){
          #***************************************************************
          if(int==FALSE & cvar[1]!=''){
            d = paste("data1",t,sep="")
            tempdata <- cbind(get(d)[avar],get(d)[mvar],get(d)[cvar])
            fit <- lm(unlist(get(d)[yvar]) ~ .,data=tempdata)
            if(nrow(summary(fit)$coefficients) != numCoefficients)
            {
              stop("Error (mediationV): There is not enough variability. Increase n or MAF.")
              # print(paste("int is ",int,sep=""))
              # print(paste("numCoefficients is ",numCoefficients,sep=""))
              # print(summary(fit))
              # return(na.matrix)
            }
            assign(paste("out1",t,sep=""),rbind(fit$coefficients,vcov(fit)))
          }
          if(int==FALSE & cvar[1]==''){
            d = paste("data1",t,sep="")
            fit <- lm(unlist(get(d)[yvar]) ~ unlist(get(d)[avar]) + unlist(get(d)[mvar]))
            if(nrow(summary(fit)$coefficients) != numCoefficients)
            {
              stop("Error (mediationV): There is not enough variability. Increase n or MAF.")
              # print(paste("int is ",int,sep=""))
              # print(paste("numCoefficients is ",numCoefficients,sep=""))
              # print(summary(fit))
              # return(na.matrix)
            }
            assign(paste("out1",t,sep=""),rbind(fit$coefficients,vcov(fit)))
          }
          if(int==TRUE & cvar[1]!=''){
            d = paste("data1",t,sep="")
            tempdata <- cbind(get(d)[avar],get(d)[mvar],get(d)[cvar],get(d)["int"])
            fit <- lm(unlist(get(d)[yvar]) ~ .,data=tempdata)
            if(nrow(summary(fit)$coefficients) != numCoefficients)
            {
              stop("Error (mediationV): There is not enough variability. Increase n or MAF.")
              # print(paste("int is ",int,sep=""))
              # print(paste("numCoefficients is ",numCoefficients,sep=""))
              # print(summary(fit))
              # return(na.matrix)
            }
            assign(paste("out1",t,sep=""),rbind(fit$coefficients,vcov(fit)))
          }
          if(int==TRUE & cvar[1]==''){
            d = paste("data1",t,sep="")
            fit <- lm(unlist(get(d)[yvar]) ~ unlist(get(d)[avar]) + unlist(get(d)[mvar]) +unlist(get(d)["int"]))
            if(nrow(summary(fit)$coefficients) != numCoefficients)
            {
              stop("Error (mediationV): There is not enough variability. Increase n or MAF.")
              # print(paste("int is ",int,sep=""))
              # print(paste("numCoefficients is ",numCoefficients,sep=""))
              # print(summary(fit))
              # return(na.matrix)
            }
            assign(paste("out1",t,sep=""),rbind(fit$coefficients,vcov(fit)))
          }
        } # End of if(yreg=="linear")
        
        #*****************************************************************************
        if (yreg=="logistic"  | yreg=="loglinear" | yreg=="poisson" | yreg=="negbin"){
          #*****************************************************************************
          if(int==TRUE & cvar[1]!=''){
            if (yreg=="logistic"){
              d = paste("data1",t,sep="")
              tempdata <- cbind(get(d)[avar],get(d)[mvar],get(d)[cvar],get(d)["int"])
              fit <- glm(as.factor(unlist(get(d)[yvar])) ~ .,data=tempdata,family=binomial)

              if(nrow(summary(fit)$coefficients) != numCoefficients)
              {
                stop("Error (mediationV): There is not enough variability. Increase n or MAF.")
                # print(paste("int is ",int,sep=""))
                # print(paste("numCoefficients is ",numCoefficients,sep=""))
                # print(summary(fit))
                # print("Line 213")
                # return(na.matrix)
              }
              assign(paste("out1",t,sep=""),rbind(fit$coefficients,vcov(fit)))
            } # End of if (yreg=="logistic")
            if(yreg=="loglinear"){
              d = paste("data1",t,sep="")
              tempdata <- cbind(get(d)[avar],get(d)[mvar],get(d)[cvar],get(d)["int"])
              fit <- glm(unlist(get(d)[yvar]) ~ .,family=poisson(link = "log"),data=tempdata)

              if(nrow(summary(fit)$coefficients) != numCoefficients)
              {
                stop("Error (mediationV): There is not enough variability. Increase n or MAF.")
                # print(paste("int is ",int,sep=""))
                # print(paste("numCoefficients is ",numCoefficients,sep=""))
                # print(summary(fit))
                # return(na.matrix)
              }
              cov = vcov(fit)
              gmparms <- summary(fit)$coefficients
              gmparms <- as.data.frame(gmparms)
              par <- gmparms$Estimate
              assign(paste("out1",t,sep=""),rbind(par,cov))
            } # End of if(yreg=="loglinear")
            if(yreg=="poisson"){
              d = paste("data1",t,sep="")
              tempdata <- cbind(get(d)[avar],get(d)[mvar],get(d)[cvar],get(d)["int"])
              fit <- glm(unlist(get(d)[yvar]) ~ .,family=poisson(link = "log"),data=tempdata)

              if(nrow(summary(fit)$coefficients) != numCoefficients)
              {
                stop("Error (mediationV): There is not enough variability. Increase n or MAF.")
                # print(paste("int is ",int,sep=""))
                # print(paste("numCoefficients is ",numCoefficients,sep=""))
                # print(summary(fit))
                # return(na.matrix)
              }
              cov = vcov(fit)
              gmparms <- summary(fit)$coefficients
              gmparms <- as.data.frame(gmparms)
              par <- gmparms$Estimate
              assign(paste("out1",t,sep=""),rbind(par,cov))
            } # end of &yreg=poisson
            if(yreg=="negbin"){
              d = paste("data1",t,sep="")
              tempdata <- cbind(get(d)[avar],get(d)[mvar],get(d)[cvar],get(d)["int"])
              fit <- glm.nb(unlist(get(d)[yvar]) ~ .,data=tempdata)

              if(nrow(summary(fit)$coefficients) != numCoefficients)
              {
                stop("Error (mediationV): There is not enough variability. Increase n or MAF.")
                # print(paste("int is ",int,sep=""))
                # print(paste("numCoefficients is ",numCoefficients,sep=""))
                # print(summary(fit))
                # return(na.matrix)
              }
              cov = vcov(fit)
              gmparms <- summary(fit)$coefficients
              gmparms <- as.data.frame(gmparms)
              par <- gmparms$Estimate
              assign(paste("out1",t,sep=""),rbind(par,cov))
            } # end of &yreg=negbin
          } # End of if(int==TRUE & cvar!='') # Completely accounted for
          if(int==TRUE & cvar[1]==''){
            if(yreg=="logistic"){
              d = paste("data1",t,sep="")
              fit <- glm(as.factor(unlist(get(d)[yvar])) ~ unlist(get(d)[avar]) + unlist(get(d)[mvar]) 
                         + unlist(get(d)["int"]),family=binomial)
              if(nrow(summary(fit)$coefficients) != numCoefficients)
              {
                stop("Error (mediationV): There is not enough variability. Increase n or MAF.")
                # print(paste("int is ",int,sep=""))
                # print(paste("numCoefficients is ",numCoefficients,sep=""))
                # print(summary(fit))
                # return(na.matrix)
              }
              assign(paste("out1",t,sep=""),rbind(fit$coefficients,vcov(fit)))
            } # end of if(yreg="logistic")
            if(yreg=="loglinear"){
              d = paste("data1",t,sep="")
              fit <- glm(unlist(get(d)[yvar]) ~ unlist(get(d)[avar]) + unlist(get(d)[mvar]) + unlist(get(d)["int"])
                         ,family=poisson(link = "log"))
              if(nrow(summary(fit)$coefficients) != numCoefficients)
              {
                stop("Error (mediationV): There is not enough variability. Increase n or MAF.")
                # print(paste("int is ",int,sep=""))
                # print(paste("numCoefficients is ",numCoefficients,sep=""))
                # print(summary(fit))
                # return(na.matrix)
              }
              cov = vcov(fit)
              gmparms <- summary(fit)$coefficients
              gmparms <- as.data.frame(gmparms)
              par <- gmparms$Estimate
              assign(paste("out1",t,sep=""),rbind(par,cov))
            } # end of if(yreg=""loglinear"")
            if(yreg=="poisson"){
              d = paste("data1",t,sep="")
              fit <- glm(unlist(get(d)[yvar]) ~ unlist(get(d)[avar]) + unlist(get(d)[mvar]) 
                         + unlist(get(d)["int"]),family=poisson(link = "log"))
              if(nrow(summary(fit)$coefficients) != numCoefficients)
              {
                stop("Error (mediationV): There is not enough variability. Increase n or MAF.")
                # print(paste("int is ",int,sep=""))
                # print(paste("numCoefficients is ",numCoefficients,sep=""))
                # print(summary(fit))
                # return(na.matrix)
              }
              cov = vcov(fit)
              gmparms <- summary(fit)$coefficients
              gmparms <- as.data.frame(gmparms)
              par <- gmparms$Estimate
              assign(paste("out1",t,sep=""),rbind(par,cov))
            } # end of yreg=poisson
            if(yreg=="negbin"){
              d = paste("data1",t,sep="")
              fit <- glm.nb(unlist(get(d)[yvar]) ~ unlist(get(d)[avar]) + unlist(get(d)[mvar]) 
                            + unlist(get(d)["int"]))
              if(nrow(summary(fit)$coefficients) != numCoefficients)
              {
                stop("Error (mediationV): There is not enough variability. Increase n or MAF.")
                # print(paste("int is ",int,sep=""))
                # print(paste("numCoefficients is ",numCoefficients,sep=""))
                # print(summary(fit))
                # return(na.matrix)
              }
              cov = vcov(fit)
              gmparms <- summary(fit)$coefficients
              gmparms <- as.data.frame(gmparms)
              par <- gmparms$Estimate
              assign(paste("out1",t,sep=""),rbind(par,cov))
            } # End of if(yreg=="negbin")
          } # End of if(int==TRUE & cvar=='')
          
          if(int==FALSE & cvar[1]!=''){
            if(yreg=="logistic"){
              d = paste("data1",t,sep="")
              tempdata <- cbind(get(d)[avar],get(d)[mvar],get(d)[cvar])
              fit <- glm(as.factor(unlist(get(d)[yvar])) ~ .,data=tempdata,family=binomial)

              if(nrow(summary(fit)$coefficients) != numCoefficients)
              {
                stop("Error (mediationV): There is not enough variability. Increase n or MAF.")
                # print(paste("int is ",int,sep=""))
                # print(paste("numCoefficients is ",numCoefficients,sep=""))
                # print(summary(fit))
                # return(na.matrix)
              }
              assign(paste("out1",t,sep=""),rbind(fit$coefficients,vcov(fit)))
            } # end of if(yreg=="logistic")
            if(yreg=="loglinear"){
              d = paste("data1",t,sep="")
              tempdata <- cbind(get(d)[avar],get(d)[mvar],get(d)[cvar])
              fit <- glm(unlist(get(d)[yvar]) ~ .,family=poisson(link = "log"),data=tempdata)
              
              if(nrow(summary(fit)$coefficients) != numCoefficients)
              {
                stop("Error (mediationV): There is not enough variability. Increase n or MAF.")
                # print(paste("int is ",int,sep=""))
                # print(paste("numCoefficients is ",numCoefficients,sep=""))
                # print(summary(fit))
                # return(na.matrix)
              }
              cov = vcov(fit)
              gmparms <- summary(fit)$coefficients
              gmparms <- as.data.frame(gmparms)
              par <- gmparms$Estimate
              assign(paste("out1",t,sep=""),rbind(par,cov))
            } # end of if(yreg=="loglinear")
            if(yreg=="poisson"){
              d = paste("data1",t,sep="")
              tempdata <- cbind(get(d)[avar],get(d)[mvar],get(d)[cvar])
              fit <- glm(unlist(get(d)[yvar]) ~ .,family=poisson(link = "log"),data=tempdata)
              
              if(nrow(summary(fit)$coefficients) != numCoefficients)
              {
                stop("Error (mediationV): There is not enough variability. Increase n or MAF.")
                # print(paste("int is ",int,sep=""))
                # print(paste("numCoefficients is ",numCoefficients,sep=""))
                # print(summary(fit))
                # return(na.matrix)
              }
              cov = vcov(fit)
              gmparms <- summary(fit)$coefficients
              gmparms <- as.data.frame(gmparms)
              par <- gmparms$Estimate
              assign(paste("out1",t,sep=""),rbind(par,cov))
            } # end of if &yreg=poisson
            if(yreg=="negbin"){
              d = paste("data1",t,sep="")
              tempdata <- cbind(get(d)[avar],get(d)[mvar],get(d)[cvar])
              fit <- glm.nb(unlist(get(d)[yvar]) ~ .,data=tempdata)
              
              if(nrow(summary(fit)$coefficients) != numCoefficients)
              {
                stop("Error (mediationV): There is not enough variability. Increase n or MAF.")
                # print(paste("int is ",int,sep=""))
                # print(paste("numCoefficients is ",numCoefficients,sep=""))
                # print(summary(fit))
                # return(na.matrix)
              }
              cov = vcov(fit)
              gmparms <- summary(fit)$coefficients
              gmparms <- as.data.frame(gmparms)
              par <- gmparms$Estimate
              assign(paste("out1",t,sep=""),rbind(par,cov))
            }# end of if(yreg=="negbin")
          } # end of &int=FALSE & &cvar^=
          
          if(int==FALSE & cvar[1]==''){
            if(yreg=="logistic"){
              d = paste("data1",t,sep="")
              fit <- glm(as.factor(unlist(get(d)[yvar])) ~ unlist(get(d)[avar]) + unlist(get(d)[mvar]),
                         family=binomial)
              if(nrow(summary(fit)$coefficients) != numCoefficients)
              {
                stop("Error (mediationV): There is not enough variability. Increase n or MAF.")
                # print(paste("int is ",int,sep=""))
                # print(paste("numCoefficients is ",numCoefficients,sep=""))
                # print(summary(fit))
                # return(na.matrix)
              }
              assign(paste("out1",t,sep=""),rbind(fit$coefficients,vcov(fit)))
            } # end of if(yreg=="logistic")
            if(yreg=="loglinear"){
              d = paste("data1",t,sep="")
              fit <- glm(unlist(get(d)[yvar]) ~ unlist(get(d)[avar]) + unlist(get(d)[mvar]),family=poisson(link = "log"))
              if(nrow(summary(fit)$coefficients) != numCoefficients)
              {
                stop("Error (mediationV): There is not enough variability. Increase n or MAF.")
                # print(paste("int is ",int,sep=""))
                # print(paste("numCoefficients is ",numCoefficients,sep=""))
                # print(summary(fit))
                # return(na.matrix)
              }
              cov = vcov(fit)
              gmparms <- summary(fit)$coefficients
              gmparms <- as.data.frame(gmparms)
              par <- gmparms$Estimate
              assign(paste("out1",t,sep=""),rbind(par,cov))
            } # end of if(yreg="loglinear")
            if(yreg=="poisson"){
              d = paste("data1",t,sep="")
              fit <- glm(unlist(get(d)[yvar]) ~ unlist(get(d)[avar]) + unlist(get(d)[mvar]) ,family=poisson(link = "log"))
              if(nrow(summary(fit)$coefficients) != numCoefficients)
              {
                stop("Error (mediationV): There is not enough variability. Increase n or MAF.")
                # print(paste("int is ",int,sep=""))
                # print(paste("numCoefficients is ",numCoefficients,sep=""))
                # print(summary(fit))
                # return(na.matrix)
              }
              cov = vcov(fit)
              gmparms <- summary(fit)$coefficients
              gmparms <- as.data.frame(gmparms)
              par <- gmparms$Estimate
              assign(paste("out1",t,sep=""),rbind(par,cov))
            } # end of if(yreg=poisson)
            if(yreg=="negbin"){
              d = paste("data1",t,sep="")
              fit <- glm.nb(unlist(get(d)[yvar]) ~ unlist(get(d)[avar]) + unlist(get(d)[mvar]))
              if(nrow(summary(fit)$coefficients) != numCoefficients)
              {
                stop("Error (mediationV): There is not enough variability. Increase n or MAF.")
                # print(paste("int is ",int,sep=""))
                # print(paste("numCoefficients is ",numCoefficients,sep=""))
                # print(summary(fit))
                # return(na.matrix)
              }
              cov = vcov(fit)
              gmparms <- summary(fit)$coefficients
              gmparms <- as.data.frame(gmparms)
              par <- gmparms$Estimate
              assign(paste("out1",t,sep=""),rbind(par,cov))
            } # end of if(yreg="negbin")
          } # end of if int==FALSE & cvar==''
          
        } # End of if (yreg=="logistic"  | yreg=="loglinear" | yreg==poisson | yreg=="negbin") # Accounted for
        
        #************************************************************
        if(mreg=="linear" & yreg=="linear"){
          if(cvar[1]!=''){
            if(casecontrol!=TRUE){
              d = paste("data1",t,sep="")
              tempdata <- cbind(get(d)[avar],get(d)[cvar])
              fit <- lm(unlist(get(d)[mvar]) ~ .,data=tempdata)
              if(nrow(summary(fit)$coefficients) != (numCoefficients - (int*1 +1)))
              {
                stop("Error (mediationV): There is not enough variability. Increase n or MAF.")
                # print(paste("int is ",int,sep=""))
                # print(paste("numCoefficients is ",numCoefficients,sep=""))
                # print(summary(fit))
                # return(na.matrix)
              }
              assign(paste("out2",t,sep=""),rbind(fit$coefficients,vcov(fit)))
              # if(t==100){print(get(paste("out2",t,sep="")))}
            }
            if(casecontrol==TRUE){
              d = paste("data1",t,sep="")
              nums <- which(get(d)[yvar]==0)
              newset <- c()
              for(s in nums){
                newset <- rbind(newset,get(d)[s,1:ncol(get(d))])
              }
              tempdata <- cbind(newset[avar],newset[cvar])
              fit <- lm(unlist(newset[mvar]) ~ .,data=tempdata)
              if(nrow(summary(fit)$coefficients) != (numCoefficients - (int*1 +1)))
              {
                stop("Error (mediationV): There is not enough variability. Increase n or MAF.")
                # print(paste("int is ",int,sep=""))
                # print(paste("numCoefficients is ",numCoefficients,sep=""))
                # print(summary(fit))
                # return(na.matrix)
              }
              assign(paste("out2",t,sep=""),rbind(fit$coefficients,vcov(fit)))
            }
          } # end of if &cvar^= then do;
          if(cvar[1]==''){
            if(casecontrol!=TRUE){
              d = paste("data1",t,sep="")
              fit <- lm(unlist(get(d)[mvar]) ~ unlist(get(d)[avar]))
              if(nrow(summary(fit)$coefficients) != (numCoefficients - (int*1 +1)))
              {
                stop("Error (mediationV): There is not enough variability. Increase n or MAF.")
                # print(paste("int is ",int,sep=""))
                # print(paste("numCoefficients is ",numCoefficients,sep=""))
                # print(summary(fit))
                # return(na.matrix)
              }
              assign(paste("out2",t,sep=""),rbind(fit$coefficients,vcov(fit)))
            }
            if(casecontrol==TRUE){
              d = paste("data1",t,sep="")
              nums <- which(get(d)[yvar]==0)
              newset <- c()
              for(s in nums){
                newset <- rbind(newset,get(d)[s,1:ncol(get(d))])
              }
              fit <- lm(unlist(newset[mvar]) ~ unlist(newset[avar]))
              if(nrow(summary(fit)$coefficients) != (numCoefficients - (int*1 +1)))
              {
                stop("Error (mediationV): There is not enough variability. Increase n or MAF.")
                # print(paste("int is ",int,sep=""))
                # print(paste("numCoefficients is ",numCoefficients,sep=""))
                # print(summary(fit))
                # return(na.matrix)
              }
              assign(paste("out2",t,sep=""),rbind(fit$coefficients,vcov(fit)))
            }
          } # end of if  &cvar= then do;
        } # end of if  &mreg="linear" & &yreg="linear" then do
        
        #*************************************************************
        if(mreg=="linear" & yreg!="linear"){
          if(cvar[1]!=''){
            if(casecontrol!=TRUE){
              d = paste("data1",t,sep="")
              tempdata <- cbind(get(d)[avar],get(d)[cvar])
              fit <- lm(unlist(get(d)[mvar]) ~ .,data=tempdata)
              if(nrow(summary(fit)$coefficients) != (numCoefficients - (int*1 +1)))
              {
                stop("Error (mediationV): There is not enough variability. Increase n or MAF.")
                # print(paste("int is ",int,sep=""))
                # print(paste("numCoefficients is ",numCoefficients,sep=""))
                # print(summary(fit))
                # return(na.matrix)
              }
              temp <- rbind(fit$coefficients,vcov(fit))
              RMSE <- rep(summary(fit)$sigma,nrow(temp))
              assign(paste("out2",t,sep=""),cbind(RMSE,temp))
            }
            if(casecontrol==TRUE){
              d = paste("data1",t,sep="")
              nums <- which(get(d)[yvar]==0)
              newset <- c()
              for(s in nums){
                newset <- rbind(newset,get(d)[s,1:ncol(get(d))])
              }
              tempdata <- cbind(newset[avar],newset[cvar])
              fit <- lm(unlist(newset[mvar]) ~ .,data=tempdata)
              if(nrow(summary(fit)$coefficients) != (numCoefficients - (int*1 +1)))
              {
                stop("Error (mediationV): There is not enough variability. Increase n or MAF.")
                # print(paste("int is ",int,sep=""))
                # print(paste("numCoefficients is ",numCoefficients,sep=""))
                # print(summary(fit))
                # return(na.matrix)
              }
              temp <- rbind(fit$coefficients,vcov(fit))
              RMSE <- rep(summary(fit)$sigma,nrow(temp))
              assign(paste("out2",t,sep=""),cbind(RMSE,temp))
            }
          } # end of if &cvar^= then do;
          if(cvar[1]==''){
            if(casecontrol!=TRUE){
              d = paste("data1",t,sep="")
              fit <- lm(unlist(get(d)[mvar]) ~ unlist(get(d)[avar]))
              if(nrow(summary(fit)$coefficients) != (numCoefficients - (int*1 +1)))
              {
                stop("Error (mediationV): There is not enough variability. Increase n or MAF.")
                # print(paste("int is ",int,sep=""))
                # print(paste("numCoefficients is ",numCoefficients,sep=""))
                # print(summary(fit))
                # return(na.matrix)
              }
              temp <- rbind(fit$coefficients,vcov(fit))
              RMSE <- rep(summary(fit)$sigma,nrow(temp))
              assign(paste("out2",t,sep=""),cbind(RMSE,temp))
            }
            if(casecontrol==TRUE){
              d = paste("data1",t,sep="")
              nums <- which(get(d)[yvar]==0)
              newset <- c()
              for(s in nums){
                newset <- rbind(newset,get(d)[s,1:ncol(get(d))])
              }
              fit <- lm(unlist(newset[mvar]) ~ unlist(newset[avar]))
              if(nrow(summary(fit)$coefficients) != (numCoefficients - (int*1 +1)))
              {
                stop("Error (mediationV): There is not enough variability. Increase n or MAF.")
                # print(paste("int is ",int,sep=""))
                # print(paste("numCoefficients is ",numCoefficients,sep=""))
                # print(summary(fit))
                # return(na.matrix)
              }
              temp <- rbind(fit$coefficients,vcov(fit))
              RMSE <- rep(summary(fit)$sigma,nrow(temp))
              assign(paste("out2",t,sep=""),cbind(RMSE,temp))
            }
          } # end of if &cvar= then do
        } # end of if &mreg="linear" & &yreg^="linear"
        
        #************************************************************
        if(mreg=="logistic"){
          if(cvar[1]!=''){
            if(casecontrol!=TRUE){
              d = paste("data1",t,sep="")
              tempdata <- cbind(get(d)[avar],get(d)[cvar])
              fit <- glm(as.factor(unlist(get(d)[mvar])) ~ .,data=tempdata,family=binomial)

              if(nrow(summary(fit)$coefficients) != (numCoefficients - (int*1 +1)))
              {
                stop("Error (mediationV): There is not enough variability. Increase n or MAF.")
                # print(paste("int is ",int,sep=""))
                # print(paste("numCoefficients is ",numCoefficients,sep=""))
                # print(summary(fit))
                # return(na.matrix)
              }
              assign(paste("out2",t,sep=""),rbind(fit$coefficients,vcov(fit)))
            }
            if(casecontrol==TRUE){
              d = paste("data1",t,sep="")
              nums <- which(get(d)[yvar]==0)
              newset <- c()
              for(s in nums){
                newset <- rbind(newset,get(d)[s,1:ncol(get(d))])
              }
              tempdata <- cbind(newset[avar],newset[cvar])
              fit <- glm(as.factor(unlist(newset[mvar])) ~ .,data=tempdata,family=binomial)

              if(nrow(summary(fit)$coefficients) != (numCoefficients - (int*1 +1)))
              {
                stop("Error (mediationV): There is not enough variability. Increase n or MAF.")
                # print(paste("int is ",int,sep=""))
                # print(paste("numCoefficients is ",numCoefficients,sep=""))
                # print(summary(fit))
                # return(na.matrix)
              }
              assign(paste("out2",t,sep=""),rbind(fit$coefficients,vcov(fit)))
            }
          } # end of &cvar^= then do
          if(cvar[1]==''){
            if(casecontrol!=TRUE){
              d = paste("data1",t,sep="")
              fit <- glm(as.factor(unlist(get(d)[mvar])) ~ unlist(get(d)[avar]),family=binomial)
              if(nrow(summary(fit)$coefficients) != (numCoefficients - (int*1 +1)))
              {
                stop("Error (mediationV): There is not enough variability. Increase n or MAF.")
                # print(paste("int is ",int,sep=""))
                # print(paste("numCoefficients is ",numCoefficients,sep=""))
                # print(summary(fit))
                # return(na.matrix)
              }
              assign(paste("out2",t,sep=""),rbind(fit$coefficients,vcov(fit)))
            }
            if(casecontrol==TRUE){
              d = paste("data1",t,sep="")
              nums <- which(get(d)[yvar]==0)
              newset <- c()
              for(s in nums){
                newset <- rbind(newset,get(d)[s,1:ncol(get(d))])
              }
              fit <- glm(as.factor(unlist(newset[mvar])) ~ unlist(newset[avar]),family=binomial)
              if(nrow(summary(fit)$coefficients) != (numCoefficients - (int*1 +1)))
              {
                stop("Error (mediationV): There is not enough variability. Increase n or MAF.")
                # print(paste("int is ",int,sep=""))
                # print(paste("numCoefficients is ",numCoefficients,sep=""))
                # print(summary(fit))
                # return(na.matrix)
              }
              assign(paste("out2",t,sep=""),rbind(fit$coefficients,vcov(fit)))
            }
          } # end of if(cvar= ) then do
        } # end of if(mreg="logistic")
        
      } # end of for(t in 1:nboot)
      #***************** regression-for bootstrap 	END **************
      
      #***************** causal effects for bootstrap  ***********************
      # Create objects in which we save the bootstrap samples of causal effects
      if(mreg=="linear" & int==FALSE){
        bootsample = matrix(0,nrow=nboot,ncol=3)
      }
      if(mreg=="linear" & int==TRUE){
        if(cvar[1]!='' & output=="full"){
          bootsample = matrix(0,nrow=nboot,ncol=12)
        }
        if(cvar[1]=='' | (cvar[1]!='' & output!="full")){
          bootsample = matrix(0,nrow=nboot,ncol=6)
        }
      }
      if(mreg=="logistic"){
        if(cvar[1]!='' & output=="full"){
          bootsample = matrix(0,nrow=nboot,ncol=12)
        }
        if(int==FALSE & cvar[1]=='' & yreg=="linear"){
          bootsample = matrix(0,nrow=nboot,ncol=3)
        }
        if((int==TRUE & cvar[1]=='') | (cvar[1]!='' & output!="full") | (int==FALSE & cvar[1]=='' & yreg!="linear" )){
          bootsample = matrix(0,nrow=nboot,ncol=6)
        }
      }
      ######## */compute the causal effects*/
      if((mreg=="linear" & int==FALSE) | (yreg=="linear" & mreg=="logistic" & int==FALSE & cvar[1]=='')){
        for(t in 1:nboot){
          vb = get(paste("out2",t,sep=""))
          if(yreg=="linear"){
            beta0 = vb[1,1]
            beta1 = vb[1,2]
          }
          if(yreg!="linear"){
            beta1=vb[1,3]
          }
          vb = get(paste("out1",t,sep=""))
          theta1=vb[1,2]
          theta2=vb[1,3]
          ##*/cde and nde*
          if(yreg=="linear" & mreg=="logistic"){
            bootsample[t,1]=(theta1)*(a1-a0)
            #*/nie*/;
            bootsample[t,2]=(theta2)*(exp(beta0+beta1*a1)/(1+exp(beta0+beta1*a1))-exp(beta0+beta1*a0)/(1+exp(beta0+beta1*a0)))
            #*/te*/;
            bootsample[t,3]=bootsample[t,1]+bootsample[t,2]
          }
          if(yreg=="linear" & mreg=="linear" & int==FALSE){
            bootsample[t,1]=((theta1)*(a1-a0))
            #*/nie*/;
            bootsample[t,2]=((theta2*beta1)*(a1-a0))
            #*/te*/;
            bootsample[t,3]=((theta1+theta2*beta1)*(a1-a0))
          }
          if(yreg!="linear" & mreg=="linear" ){
            bootsample[t,1]=exp((theta1)*(a1-a0))
            bootsample[t,2]=exp((theta2*beta1)*(a1-a0))
            bootsample[t,3]=bootsample[t,1]*bootsample[t,2]
          }
        } # end of for (t in 1:nboot)
        bootdata = bootsample
        colnames(bootdata) = c("boot1", "boot2", "boot3")
        bootdata <- as.data.frame(bootdata)
        
      } # End of if((mreg=="linear" & int==FALSE) | (yreg=="linear" & mreg=="logistic"...
      
      if((int==TRUE)|(mreg=="logistic" & int==FALSE & cvar[1]!='') |
         (yreg!="linear" & mreg=="logistic" & int==FALSE)){
        if(cvar[1]!=''){
          vb = data2
          if(c==''){
            cmean=vb[1,1:ncol(vb)]
          }
          if(c!=''){
            cmean=vb[1,1:(ncol(vb)-nc)]
            c=vb[1,(ncol(vb)-nc+1):ncol(vb)] 
          }
        }
        if(cvar[1]!=''){
          for(t in 1:nboot){
            vb = get(paste("out1",t,sep=""))
            theta1=vb[1,2];
            theta2=vb[1,3];
            if(int==TRUE){
              theta3=vb[1,ncol(vb)]
            }
            
            vb = get(paste("out2",t,sep=""))
            
            if((yreg=="linear" & mreg=="linear") | (mreg=="logistic")){
              beta0=vb[1,1]
              beta1=vb[1,2]
              beta2= vb[1,3:ncol(vb)]
            }
            if((yreg!="linear" & mreg=="linear")){
              s2=vb[1,1]
              s2=s2**2
              beta0=vb[1,2]
              beta1=vb[1,3]
              beta2= vb[1,4:ncol(vb)]
              tsq=(theta3**2)
              rm=s2
              asq=(a1**2)
              a1sq=(a0**2)
            }
            if((yreg=="linear" & mreg=="linear" & int==TRUE)){
              # print(cmean)
              #*/MARGINAL CDE*/;
              bootsample[t,1]=(theta1+theta3*m)*(a1-a0)
              #*/MARGINAL NDE*/;
              bootsample[t,2]=(theta1+theta3*beta0+theta3*beta1*a0+(theta3*beta2%*%t(t(cmean))))*(a1-a0)
              #*/MARGINAL NIE*/;
              bootsample[t,3]=(theta2*beta1+theta3*beta1*a0)*(a1-a0)
              #*/ MARGINAL TNDE*/;
              bootsample[t,4]=(theta1+theta3*beta0+theta3*beta1*a1+(theta3*beta2%*%t(t(cmean))))*(a1-a0)
              #*/ MARGINAL TNIE*/;
              bootsample[t,5]=(theta2*beta1+theta3*beta1*a1)*(a1-a0)
              #*/te marginal*/;
              bootsample[t,6]=(theta1+theta3*beta0+theta3*beta1*a0+(theta3*beta2%*%t(t(cmean)))+theta2*beta1+theta3*beta1*a1)*(a1-a0)
              if(c!=''){
                #*/CONDITIONAL CDE*/;
                bootsample[t,7]=(theta1)*(a1-a0)+(theta3*(m))*(a1-a0)
                #*/CONDITIONAL NDE*/;
                bootsample[t,8]=(theta1+theta3*beta0+theta3*beta1*a0+(theta3*beta2*t(c)))*(a1-a0)
                #*/CONDITIONAL NIE*/;
                bootsample[t,9]=(theta2*beta1+theta3*beta1*a0)*(a1-a0)
                #*/CONDITIONAL TNDE*/;
                bootsample[t,10]=(theta1+theta3*beta0+theta3*beta1*a1+(theta3*beta2*t(c)))*(a1-a0)
                #*/ CONDITIONAL TNIE*/;
                bootsample[t,11]=(theta2*beta1+theta3*beta1*a1)*(a1-a0)
                #*/te conditional*/;	
                bootsample[t,12]=(theta1+theta3*beta0+theta3*beta1*a0+(theta3*beta2*t(c))+theta2*beta1+theta3*beta1*a1)*(a1-a0)
              } # end of if(c!='') 
            } # end of if((yreg=="linear" & mreg=="linear" & int==TRUE))
            if((yreg=="linear" & mreg=="logistic" & int==TRUE)){
              #*/MARGINAL CDE*/;
              bootsample[t,1]=(theta1+theta3*m)*(a1-a0)
              #*/MARGINAL NDE*/;
              bootsample[t,2]=(theta1+theta3*exp(beta0+beta1*a0+sum(beta2*t(cmean)))/(1+exp(beta0+beta1*a0+sum(beta2*t(cmean)))))*(a1-a0)
              #*/MARGINAL NIE*/;
              bootsample[t,3]=(theta2+theta3*a0)*(exp(beta0+beta1*a1+sum(beta2*t(cmean)))/(1+exp(beta0+beta1*a1+sum(beta2*t(cmean))))-
                                                    exp(beta0+beta1*a0+sum(beta2*t(cmean)))/(1+exp(beta0+beta1*a0+sum(beta2*t(cmean)))))
              #*/ MARGINAL TNDE*/;
              bootsample[t,4]=(theta1+theta3*exp(beta0+beta1*a1+sum(beta2*t(cmean)))/(1+exp(beta0+beta1*a1+sum(beta2*t(cmean)))))*(a1-a0)
              #*/ MARGINAL TNIE*/;
              bootsample[t,5]=(theta2+theta3*a1)*(exp(beta0+beta1*a1+sum(beta2*t(cmean)))/(1+exp(beta0+beta1*a1+sum(beta2*t(cmean))))-
                                                    exp(beta0+beta1*a0+sum(beta2*t(cmean)))/(1+exp(beta0+beta1*a0+sum(beta2*t(cmean)))))
              #*/te marginal*/;
              bootsample[t,6]=bootsample[t,2]+bootsample[t,5]
              if(c!=''){
                #*/CONDITIONAL CDE*/;
                bootsample[t,7]=(theta1)*(a1-a0)+(theta3*(m))*(a1-a0)
                #*/CONDITIONAL NDE*/;
                bootsample[t,8]=(theta1+theta3*exp(beta0+beta1*a0+sum(beta2*t(c)))/(1+exp(beta0+beta1*a0+sum(beta2*t(c)))))*(a1-a0)
                #*/CONDITIONAL NIE*/;
                bootsample[t,9]=(theta2+theta3*a0)*(exp(beta0+beta1*a1+sum(beta2*t(c)))/(1+exp(beta0+beta1*a1+sum(beta2*t(c))))-
                                                      exp(beta0+beta1*a0+sum(beta2*t(c)))/(1+exp(beta0+beta1*a0+sum(beta2*t(c)))))
                #*/CONDITIONAL TNDE*/;
                bootsample[t,10]=(theta1+theta3*exp(beta0+beta1*a1+sum(beta2*t(c)))/(1+exp(beta0+beta1*a1+sum(beta2*t(c)))))*(a1-a0)
                #*/ CONDITIONAL TNIE*/;
                bootsample[t,11]=(theta2+theta3*a1)*(exp(beta0+beta1*a1+sum(beta2*t(c)))/(1+exp(beta0+beta1*a1+sum(beta2*t(c))))-
                                                       exp(beta0+beta1*a0+sum(beta2*t(c)))/(1+exp(beta0+beta1*a0+sum(beta2*t(c)))))
                #*/te conditional*/;	
                bootsample[t,12]=bootsample[t,8]+bootsample[t,11]
              } # end of if(c!='')
            } # end of if((yreg=="linear" & mreg=="logistic" & int==TRUE))
            if((yreg=="linear" & mreg=="logistic" & int==FALSE)){
              #*/MARGINAL CDE#*/
              bootsample[t,1]=(theta1)*(a1-a0)
              #*/MARGINAL NDE#*/
              bootsample[t,2]=(theta1)*(a1-a0)
              #*/MARGINAL NIE#*/
              bootsample[t,3]=(theta2)*(exp(beta0+beta1*a1+sum(beta2*t(cmean)))/(1+exp(beta0+beta1*a1+sum(beta2*t(cmean))))-
                                          exp(beta0+beta1*a0+sum(beta2*t(cmean)))/(1+exp(beta0+beta1*a0+sum(beta2*t(cmean)))))
              #*/ MARGINAL TNDE#*/
              bootsample[t,4]=(theta1)*(a1-a0)
              #*/ MARGINAL TNIE#*/
              bootsample[t,5]=(theta2)*(exp(beta0+beta1*a1+sum(beta2*t(cmean)))/(1+exp(beta0+beta1*a1+sum(beta2*t(cmean))))-
                                          exp(beta0+beta1*a0+sum(beta2*t(cmean)))/(1+exp(beta0+beta1*a0+sum(beta2*t(cmean)))))
              #*/te marginal#*/
              bootsample[t,6]=bootsample[t,2]+bootsample[t,5]
              if(c!=''){
                #*/CONDITIONAL CDE#*/
                bootsample[t,7]=(theta1)*(a1-a0)
                #*/CONDITIONAL NDE#*/
                bootsample[t,8]=(theta1)*(a1-a0)
                #*/CONDITIONAL NIE#*/
                bootsample[t,9]=(theta2)*(exp(beta0+beta1*a1+sum(beta2*t(c)))/(1+exp(beta0+beta1*a1+sum(beta2*t(c))))-
                                            exp(beta0+beta1*a0+sum(beta2*t(c)))/(1+exp(beta0+beta1*a0+sum(beta2*t(c)))))
                #*/CONDITIONAL TNDE#*/
                bootsample[t,10]=(theta1)*(a1-a0)
                #*/ CONDITIONAL TNIE#*/
                bootsample[t,11]=(theta2)*(exp(beta0+beta1*a1+sum(beta2*t(c)))/(1+exp(beta0+beta1*a1+sum(beta2*t(c))))-
                                             exp(beta0+beta1*a0+sum(beta2*t(c)))/(1+exp(beta0+beta1*a0+sum(beta2*t(c)))))
                #*/te conditional#*/	
                bootsample[t,12]=bootsample[t,8]+bootsample[t,11]
              } # end of if(c!='')
            } # end of if((yreg=="linear" & mreg=="logistic" & int==FALSE))
            if((yreg!="linear" & mreg=="linear" & int==TRUE)){
              #*/MARGINAL CDE#*/
              x6=(theta1+theta3*m)*(a1-a0)
              bootsample[t,1]=exp(x6)
              #*/MARGINAL NDE#*/
              x7=(theta1+theta3*beta0+theta3*beta1*a0+sum(theta3*beta2*t(cmean))+theta3*theta2*rm)*(a1-a0)+1/2*tsq*rm*(asq-a1sq)
              bootsample[t,2]=exp(x7)
              #*/MARGINAL NIE#*/
              x8=(theta2*beta1+theta3*beta1*a0)*(a1-a0)
              bootsample[t,3]=exp(x8)
              #*/ MARGINAL TNDE#*/
              x9=(theta1+theta3*beta0+theta3*beta1*a1+sum(theta3*beta2*t(cmean))+theta3*theta2*rm)*(a1-a0)+1/2*tsq*rm*(asq-a1sq)
              bootsample[t,4]=exp(x9)
              #*/ MARGINAL TNIE#*/
              x10=(theta2*beta1+theta3*beta1*a1)*(a1-a0)
              bootsample[t,5]=exp(x10)
              #*/te#*/
              bootsample[t,6]=bootsample[t,2]*bootsample[t,5]
              if(c!=''){
                #*/CONDITIONAL CDE#*/
                bootsample[t,7]=exp((theta1+theta3*m)*(a1-a0))
                #*/CONDITIONAL NDE#*/	  
                bootsample[t,8]=exp((theta1+theta3*beta0+theta3*beta1*a0+sum(theta3*beta2*t(c))
                                     +theta3*theta2*rm)*(a1-a0)+(1/2)*tsq*rm*(asq-a1sq))
                #*/CONDITIONAL NIE#*/
                x3=(theta2*beta1+theta3*beta1*a0)*(a1-a0)
                bootsample[t,9]=exp(x3)
                #*/CONDITIONAL TNDE#*/
                x4=(theta1+theta3*beta0+theta3*beta1*a1+sum(theta3*beta2*t(c))+theta3*theta2*rm)*(a1-a0)+(1/2)*tsq*rm*(asq-a1sq)
                bootsample[t,10]=exp(x4)
                #*/ CONDITIONAL TNIE#*/
                x5=(theta2*beta1+theta3*beta1*a1)*(a1-a0)
                bootsample[t,11]=exp(x5)
                #*/te#*/
                bootsample[t,12]=bootsample[t,8]*bootsample[t,11]
              } # end of if(c!='')
            } # end of if((yreg!="linear" & mreg=="linear" & int==TRUE))
            if((yreg!="linear" & mreg=="logistic" & int==FALSE)){
              #*/MARGINAL CDE#*/
              x6=(theta1)*(a1-a0)
              bootsample[t,1]=exp(x6)
              #*/MARGINAL NDE#*/
              bootsample[t,2]=exp((theta1)*(a1-a0))*(1+exp(theta2+beta0+beta1*a0+sum(beta2*t(cmean))))/
                (1+exp(theta2+beta0+beta1*a0+sum(beta2*t(cmean))))
              #*/MARGINAL NIE#*/
              bootsample[t,3]=((1+exp(beta0+beta1*a0+sum(beta2*t(cmean))))*(1+exp(theta2+beta0+beta1*a1+sum(beta2*t(cmean)))))/
                ((1+exp(beta0+beta1*a1+sum(beta2*t(cmean))))*(1+exp(theta2+beta0+beta1*a0+sum(beta2*t(cmean)))))
              #*/ MARGINAL TNDE#*/
              bootsample[t,4]=exp((theta1)*(a1-a0))*(1+exp(theta2+beta0+beta1*a1+sum(beta2*t(cmean))))/
                (1+exp(theta2+beta0+beta1*a1+sum(beta2*t(cmean))))
              #*/ MARGINAL TNIE#*/
              bootsample[t,5]=((1+exp(beta0+beta1*a0+sum(beta2*t(cmean))))*(1+exp(theta2+beta0+beta1*a1+sum(beta2*t(cmean)))))/
                ((1+exp(beta0+beta1*a1+sum(beta2*t(cmean))))*(1+exp(theta2+beta0+beta1*a0+sum(beta2*t(cmean)))))
              bootsample[t,6]=bootsample[t,2]*bootsample[t,5]
              if(c!=''){
                #*/CONDITIONAL CDE#*/
                x1=exp(theta1*(a1-a0))
                bootsample[t,7]=x1
                #*/CONDITIONAL NDE#*/
                bootsample[t,8]=exp((theta1)*(a1-a0))*(1+exp(theta2+beta0+beta1*a0+sum(beta2*t(c))))/
                  (1+exp(theta2+beta0+beta1*a0+sum(beta2*t(c))))
                #*/CONDITIONAL NIE#*/
                bootsample[t,9]=((1+exp(beta0+beta1*a0+beta2*t(c)))*(1+exp(theta2+beta0+beta1*a1+sum(beta2*t(c)))))/
                  ((1+exp(beta0+beta1*a1+sum(beta2*t(c))))*(1+exp(theta2+beta0+beta1*a0+sum(beta2*t(c)))))
                #*/CONDITIONAL TNDE#*/
                bootsample[t,10]=exp((theta1)*(a1-a0))*(1+exp(theta2+beta0+beta1*a1+sum(beta2*t(c))))/
                  (1+exp(theta2+beta0+beta1*a1+sum(beta2*t(c))))
                #*/ CONDITIONAL TNIE#*/
                bootsample[t,11]=((1+exp(beta0+beta1*a0+sum(beta2*t(c))))*(1+exp(theta2+beta0+beta1*a1+sum(beta2*t(c)))))/
                  ((1+exp(beta0+beta1*a1+sum(beta2*t(c))))*(1+exp(theta2+beta0+beta1*a0+sum(beta2*t(c)))))
                bootsample[t,12]=bootsample[t,8]*bootsample[t,11]
              } # end of if(c!='')
            }# end of if((yreg!="linear" & mreg=="logistic" & int==FALSE))
            if((yreg!="linear" & mreg=="logistic" & int==TRUE)){
              #*/MARGINAL CDE#*/
              x6=(theta1+theta3*m)*(a1-a0)
              bootsample[t,1]=exp(x6)
              #*/MARGINAL NDE#*/
              bootsample[t,2]=exp(theta1*(a1-a0))*(1+exp(theta2+theta3*a1+beta0+beta1*a0+sum(beta2*t(cmean))))/
                (1+exp(theta2+theta3*a0+beta0+beta1*a0+sum(beta2*t(cmean))))
              #*/MARGINAL NIE#*/
              bootsample[t,3]=((1+exp(beta0+beta1*a0+sum(beta2*t(cmean))))*(1+exp(theta2+theta3*a0+beta0+beta1*a1+sum(beta2*t(cmean)))))/
                ((1+exp(beta0+beta1*a1+sum(beta2*t(cmean))))*(1+exp(theta2+theta3*a0+beta0+beta1*a0+sum(beta2*t(cmean)))))
              #*/ MARGINAL TNDE#*/
              bootsample[t,4]=exp(theta1*(a1-a0))*(1+exp(theta2+theta3*a1+beta0+beta1*a1+sum(beta2*t(cmean))))/
                (1+exp(theta2+theta3*a0+beta0+beta1*a1+sum(beta2*t(cmean))))
              #*/ MARGINAL TNIE#*/
              bootsample[t,5]=((1+exp(beta0+beta1*a0+sum(beta2*t(cmean))))*
                                 (1+exp(theta2+theta3*a1+beta0+beta1*a1+sum(beta2*t(cmean)))))/
                ((1+exp(beta0+beta1*a1+sum(beta2*t(cmean))))*(1+exp(theta2+theta3*a1+beta0+beta1*a0+sum(beta2*t(cmean)))))
              bootsample[t,6]=bootsample[t,2]*bootsample[t,5]
              if(c!=''){
                #*/CONDITIONAL CDE#*/
                x1=exp((theta1+theta3*m)*(a1-a0))
                bootsample[t,7]=x1
                #*/CONDITIONAL NDE#*/
                bootsample[t,8]=exp(theta1*(a1-a0))*(1+exp(theta2+theta3*a1+beta0+beta1*a0+sum(beta2*t(c))))/
                  (1+exp(theta2+theta3*a0+beta0+beta1*a0+sum(beta2*t(c))))
                #*/CONDITIONAL NIE#*/
                bootsample[t,9]=((1+exp(beta0+beta1*a0+sum(beta2*t(c))))*(1+exp(theta2+theta3*a0+beta0+beta1*a1+sum(beta2*t(c)))))/
                  ((1+exp(beta0+beta1*a1+sum(beta2*t(c))))*(1+exp(theta2+theta3*a0+beta0+beta1*a0+sum(beta2*t(c)))))
                #*/CONDITIONAL TNDE#*/
                bootsample[t,10]=exp(theta1*(a1-a0))*
                  (1+exp(theta2+theta3*a1+beta0+beta1*a1+sum(beta2*t(c))))/
                  (1+exp(theta2+theta3*a0+beta0+beta1*a1+sum(beta2*t(c))))
                #*/ CONDITIONAL TNIE#*/
                bootsample[t,11]=((1+exp(beta0+beta1*a0+sum(beta2*t(c))))*(1+exp(theta2+theta3*a1+beta0+beta1*a1+sum(beta2*t(c)))))/
                  ((1+exp(beta0+beta1*a1+sum(beta2*t(c))))*(1+exp(theta2+theta3*a1+beta0+beta1*a0+sum(beta2*t(c)))))
                bootsample[t,12]=bootsample[t,8]*bootsample[t,11]
              } # end of if(c!='')
            } # end of if((yreg!="linear" & mreg=="logistic" & int==TRUE))
          } # end of for(t in 1:nboot)
          if(c!=''){
            bootdata = bootsample
            colnames(bootdata) <- c("boot1","boot2","boot3","boot4","boot5","boot6","boot7","boot8","boot9","boot10","boot11","boot12")
            bootdata <- as.data.frame(bootdata)
          }
          if(c==''){
            bootdata = bootsample[,1:6]
            colnames(bootdata) <- c("boot1","boot2","boot3","boot4","boot5","boot6")
            bootdata <- as.data.frame(bootdata)
          }
        } # end of if(cvar!='')
        if(cvar[1]==''){
          for(t in 1:nboot){
            vb = get(paste("out1",t,sep=""))
            NVB1= nrow(vb)
            V1=vb[2:NVB1,]
            theta1=vb[1,2]
            theta2=vb[1,3]
            if(int==TRUE){
              theta3=vb[1,4] 
            }
            vb = get(paste("out2",t,sep=""))
            NVB2= nrow(vb)
            V2=vb[2:NVB2,]
            if((yreg=="linear" & mreg=="linear") | (mreg=="logistic")){
              beta0=vb[1,1]
              beta1=vb[1,2]
            }
            if(yreg!="linear" & mreg=="linear"){
              s2=vb[1,1]
              s2=s2**2
              beta0=vb[1,2]
              beta1=vb[1,3]
              tsq=(theta3**2)
              rm=s2
              asq=(a1**2)
              a1sq=(a0**2)
            }
            if(yreg=="linear" & mreg=="linear" & int==TRUE){
              #*/CONDITIONAL=MARGINAL CDE#*/
              bootsample[t,1]=(theta1)*(a1-a0)+(theta3*(m))*(a1-a0)
              #*/CONDITIONAL=MARGINAL NDE#*/
              bootsample[t,2]=(theta1+theta3*beta0+theta3*beta1*a0)*(a1-a0)
              #*/CONDITIONAL=MARGINAL NIE#*/
              bootsample[t,3]=(theta2*beta1+theta3*beta1*a0)*(a1-a0)
              #*/CONDITIONAL=MARGINAL TNDE#*/
              bootsample[t,4]=(theta1+theta3*beta0+theta3*beta1*a1)*(a1-a0)
              #*/ CONDITIONAL=MARGINAL TNIE#*/
              bootsample[t,5]=(theta2*beta1+theta3*beta1*a1)*(a1-a0)
              #*/te#*/
              bootsample[t,6]=(theta1+theta3*beta0+theta3*beta1*a0+theta2*beta1+theta3*beta1*a1)*(a1-a0)
            }
            if(yreg=="linear" & mreg=="logistic" & int==TRUE){
              bootsample[t,1]=(theta1+theta3*m)*(a1-a0)
              #*/CONDITIONAL=MARGINAL NDE#*/
              bootsample[t,2]=(theta1+theta3*exp(beta0+beta1*a0)/(1+exp(beta0+beta1*a0)))*(a1-a0)
              #*/ CONDITIONAL=MARGINAL TNIE#*/
              bootsample[t,3]=(theta2+theta3*a0)*(exp(beta0+beta1*a1)/(1+exp(beta0+beta1*a1))-exp(beta0+beta1*a0)/(1+exp(beta0+beta1*a0)))
              #*/CONDITIONAL=MARGINAL TNDE#*/
              bootsample[t,4]=(theta1+theta3*exp(beta0+beta1*a1)/(1+exp(beta0+beta1*a1)))*(a1-a0)
              #*/ CONDITIONAL=MARGINAL TNIE#*/
              bootsample[t,5]=(theta2+theta3*a1)*(exp(beta0+beta1*a1)/(1+exp(beta0+beta1*a1))-exp(beta0+beta1*a0)/(1+exp(beta0+beta1*a0)))
              #*/te#*/
              bootsample[t,6]=bootsample[t,2]+bootsample[t,5]
            }
            if(yreg!="linear" & mreg=="linear" & int==TRUE){
              #*/MARGINAL=CONDITIONAL CDE#*/
              x1=(theta1+theta3*m)*(a1-a0)
              bootsample[t,1]=exp(x1)
              #*/MARGINAL=CONDITIONAL NDE#*/
              x2=(theta1+theta3*beta0+theta3*beta1*a0+theta3*theta2*rm)*(a1-a0)+(1/2)*tsq*rm*(asq-a1sq)
              bootsample[t,2]=exp(x2)
              #*/MARGINAL=CONDITIONAL NIE#*/
              x3=(theta2*beta1+theta3*beta1*a0)*(a1-a0)
              bootsample[t,3]=exp(x3)
              #*/MARGINAL=CONDITIONAL TNDE#*/
              x4=(theta1+theta3*beta0+theta3*beta1*a1+theta3*theta2*rm)*(a1-a0)+(1/2)*tsq*rm*(asq-a1sq)
              bootsample[t,4]=exp(x4)
              #*/ MARGINAL=CONDITIONAL TNIE#*/
              x5=(theta2*beta1+theta3*beta1*a1)*(a1-a0)
              bootsample[t,5]=exp(x5)
              #*/ MARGINAL=CONDITIONAL TE#*/
              bootsample[t,6]=bootsample[t,2]*bootsample[t,5]
            }
            if(yreg!="linear" & mreg=="logistic" & int==FALSE){
              #*/MARGINAL=CONDITIONAL CDE#*/
              x1=exp(theta1*(a1-a0))
              bootsample[t,1]=x1
              #*/MARGINAL=CONDITIONAL NDE#*/
              bootsample[t,2]=exp((theta1)*(a1-a0))*(1+exp(theta2+beta0+beta1*a0))/(1+exp(theta2+beta0+beta1*a0))
              #*/MARGINAL=CONDITIONAL NIE#*/
              bootsample[t,3]=((1+exp(beta0+beta1*a0))*(1+exp(theta2+beta0+beta1*a1)))/
                ((1+exp(beta0+beta1*a1))*(1+exp(theta2+beta0+beta1*a0)))
              #*/MARGINAL=CONDITIONAL TNDE#*/
              bootsample[t,4]=exp((theta1)*(a1-a0))*(1+exp(theta2+beta0+beta1*a1))/(1+exp(theta2+beta0+beta1*a1))
              #*/ MARGINAL=CONDITIONAL TNIE#*/
              bootsample[t,5]=((1+exp(beta0+beta1*a0))*(1+exp(theta2+beta0+beta1*a1)))/
                ((1+exp(beta0+beta1*a1))*(1+exp(theta2+beta0+beta1*a0)))
              bootsample[t,6]=bootsample[t,2]*bootsample[t,5]
            }
            if(yreg!="linear" & mreg=="logistic" & int==TRUE){
              #*/MARGINAL CDE#*/
              x6=(theta1+theta3*m)*(a1-a0)
              bootsample[t,1]=exp(x6)
              #*/MARGINAL NDE#*/
              bootsample[t,2]=exp(theta1*(a1-a0))*(1+exp(theta2+theta3*a1+beta0+beta1*a0))/(1+exp(theta2+theta3*a0+beta0+beta1*a0))
              #*/MARGINAL NIE#*/
              bootsample[t,3]=((1+exp(beta0+beta1*a0))*(1+exp(theta2+theta3*a0+beta0+beta1*a1)))/
                ((1+exp(beta0+beta1*a1))*(1+exp(theta2+theta3*a0+beta0+beta1*a0)))
              #*/ MARGINAL TNDE#*/
              bootsample[t,4]=exp(theta1*(a1-a0))*(1+exp(theta2+theta3*a1+beta0+beta1*a1))/(1+exp(theta2+theta3*a0+beta0+beta1*a1))
              #*/ MARGINAL TNIE#*/
              bootsample[t,5]=((1+exp(beta0+beta1*a0))*(1+exp(theta2+theta3*a1+beta0+beta1*a1)))/
                ((1+exp(beta0+beta1*a1))*(1+exp(theta2+theta3*a1+beta0+beta1*a0)))
              bootsample[t,6]=bootsample[t,2]*bootsample[t,5]
            }
          } # finished end to for(t in 1:nboot)
          
          bootdata = bootsample
          colnames(bootdata) <- c("boot1","boot2","boot3","boot4","boot5","boot6")
          bootdata <- as.data.frame(bootdata)
          
        } # finished end of if(cvar=='')
      } # finished goes to: if((int==TRUE)|(mreg=="logistic" & int==FALSE & cvar!='')|...
      #***************** causal effects for bootstrap END  *************************
      
      #***************** causal effects, standard errors and confidence intervals from bootstrap *************************
      #*/ no int 
      if((mreg=="linear" & int==FALSE )|(yreg=="linear" & mreg=="logistic" & int==FALSE & cvar[1]=='')){
        #*effects*
        effect = matrix(1,nrow=1,ncol=3)
        for(j in 1:3){
          effect[,j]=sum((bootdata[,j]))/nboot
        }
        colnames(effect) <- c("effect1","effect2","effect3")
        se = matrix(1,nrow=3,ncol=1)
        square = matrix(1,nrow=nboot,ncol=3)
        
        # *standard errors*
        
        for(j in 1:3){
          for(t in 1:nboot){
            square[t,j]=((bootdata[t,j])-effect[,j])^2
          }
          se[j,]=(sqrt(sum((square[,j]))))/sqrt(nboot)
        }
        y=se
        
        #* p value
        pvalue = matrix(1,nrow=1,ncol=3)
        for(j in 1:3){
          pvalue[,j] = 2*min(1-abs(pnorm((effect[,j])/(se[j,]))),abs(pnorm((effect[,j])/(se[j,]))))
        }
        
        #*Percentile confidence intervals*
        a1 = (alphalev/2)*100
        a2 = (1 - alphalev/2)*100
        ##################*#*#**#**#**#*##**#*#*#**#*#*#
        for(j in 1:3){
          temp <- paste("boot",j,sep="")
          e <- mean(bootdata[[temp]])
          cl <- quantile(bootdata[[temp]],a1/100,na.rm=TRUE)
          cu <- quantile(bootdata[[temp]],a2/100,na.rm=TRUE)
          temp <- cbind(e,cl,cu)
          rownames(temp) <- NULL
          colnames(temp) <- c(paste("effect",j,sep=""),paste("p_cil",j,sep=""),paste("p_ciu",j,sep=""))
          assign(paste("pmethod",j,sep=""),temp)
        }
        #***#*#*#**#*#*#**#**#*#**#*#**#*##**#**#**#*###*##**#*#*#*#*
        for(j in 1:3){
          vb = get(paste("pmethod",j,sep=""))
          assign(paste("cil",j,sep=""),vb[1,2])
          assign(paste("ciu",j,sep=""),vb[1,3])
        }
        cil=c(cil1,cil2,cil3)
        ciu=c(ciu1,ciu2,ciu3)
        ci = cbind(cil,ciu)
      } # end of if (mreg=="linear" & int==FALSE )| (yreg=="linear" & mreg=="logistic" ...
      
      #*/effects, standard errors, confidence intervals and p-value int*/
      #*/ other
      if((int==TRUE)|(mreg=="logistic" & int==FALSE & cvar[1]!='')|
         (yreg!="linear" & mreg=="logistic" & int==FALSE)){
        if(c!='' & cvar[1]!=''){
          effect=matrix(1,nrow=1,ncol=12)
          for(j in 1:12){
            effect[,j]=sum((bootdata[,j]))/nboot
          }
          colnames(effect) <- c("effect1","effect2","effect3","effect4","effect5",
                                "effect6","effect7","effect8","effect9","effect10","effect11","effect12")
          se=matrix(1,nrow=12,ncol=1)
          square=matrix(1,nrow=nboot,ncol=12)
          #*standard errors*
          for(j in 1:12){
            for(t in 1:nboot){
              #print(paste("t4",t,"in",nboot))
              square[t,j]=((bootdata[t,j])-effect[,j])^2
            }
            se[j,]=(sqrt(sum((square[,j]))))/sqrt(nboot);
          }
          y=se
          
          #* P-value
          # print("Line 1101")
          # print(se)
          # print(effect)
          pvalue = matrix(1,nrow=1,ncol=12)
          for(j in 1:12){
            pvalue[,j] = 2*min(1-abs(pnorm((effect[,j])/(se[j,]))),abs(pnorm((effect[,j])/(se[j,]))))
          }
          
        } # end of if(c!='' & cvar!='')
        if(cvar[1]=='' | (cvar[1]!='' & c=='')){
          effect=matrix(1,nrow=1,ncol=6)
          for(j in 1:6){
            effect[,j]=sum((bootdata[,j]))/nboot
          }
          colnames(effect) <- c("effect1","effect2","effect3","effect4","effect5","effect6")
          se=matrix(1,nrow=6,ncol=1)
          square=matrix(1,nrow=nboot,ncol=6)
          #*standard errors*
          for(j in 1:6){
            for(t in 1:nboot){
              square[t,j]=((bootdata[t,j])-effect[,j])^2
            }
            se[j,]=(sqrt(sum((square[,j]))))/sqrt(nboot)
          }
          y=se
        } # end of if(cvar=='' | (cvar!='' & c==''))
        
        #* p value
        pvalue = matrix(1,nrow=1,ncol=6)
        for(j in 1:6){
          pvalue[,j] = 2*min(1-abs(pnorm((effect[,j])/(se[j,]))),abs(pnorm((effect[,j])/(se[j,]))))
        }
        
        #*Percentile confidence intervals*
        a1 = (alphalev/2)*100
        a2 = (1 - alphalev/2)*100
        if(c!=''){
          for(j in 1:12){
            temp <- paste("boot",j,sep="")
            e <- mean(bootdata[[temp]])
            cl <- quantile(bootdata[[temp]],a1/100,na.rm=TRUE)
            cu <- quantile(bootdata[[temp]],a2/100,na.rm=TRUE)
            temp <- cbind(e,cl,cu)
            rownames(temp) <- NULL
            colnames(temp) <- c(paste("effect",j,sep=""),paste("p_cil",j,sep=""),paste("p_ciu",j,sep=""))
            assign(paste("pmethod",j,sep=""),temp)
          }
        }
        if((cvar[1]!='' & c=='') | cvar[1]==''){
          for(j in 1:6){
            temp <- paste("boot",j,sep="")
            e <- mean(bootdata[[temp]])
            cl <- quantile(bootdata[[temp]],a1/100,na.rm=TRUE)
            cu <- quantile(bootdata[[temp]],a2/100,na.rm=TRUE)
            temp <- cbind(e,cl,cu)
            rownames(temp) <- NULL
            colnames(temp) <- c(paste("effect",j,sep=""),paste("p_cil",j,sep=""),paste("p_ciu",j,sep=""))
            assign(paste("pmethod",j,sep=""),temp)
          }
        }
        if(c!='' & cvar[1]!=''){
          for(j in 1:12){
            vb = get(paste("pmethod",j,sep=""))
            assign(paste("cil",j,sep=""),vb[1,2])
            assign(paste("ciu",j,sep=""),vb[1,3])
          }
          cil=c(cil1,cil2,cil3,cil4,cil5,cil6,cil7,cil8,cil9,cil10,cil11,cil12)
          ciu=c(ciu1,ciu2,ciu3,ciu4,ciu5,ciu6,ciu7,ciu8,ciu9,ciu10,ciu11,ciu12)
          ci = cbind(cil,ciu)
        }
        if(cvar[1]=='' | (cvar[1]!='' & c=='')){
          for(j in 1:6){
            vb = get(paste("pmethod",j,sep=""))
            assign(paste("cil",j,sep=""),vb[1,2])
            assign(paste("ciu",j,sep=""),vb[1,3])
          }
          cil=c(cil1,cil2,cil3,cil4,cil5,cil6)
          ciu=c(ciu1,ciu2,ciu3,ciu4,ciu5,ciu6)
          ci = cbind(cil,ciu)
        }
      } # End of if((int==TRUE)|(mreg=="logistic" & int==FALSE & cvar!='')...
      
    } # end of if(boot!='' & boot!=FALSE)
    #***************************   BOOTSTRAP PROCEDURE -END-  ************************************************
    
    # Clear the output window
    # cat("\014")
    
    #************* regression to print ***********************************************************************
    #*********************************************************************************************************
    if(yreg=="linear"){
      #*******************************************************************************************************
      if(int==FALSE & cvar[1]!=''){
        tempdata <- cbind(data1[avar],data1[mvar],data1[cvar])
        fit <- lm(unlist(data1[yvar]) ~ .,data=tempdata)
        if(nrow(summary(fit)$coefficients) != numCoefficients)
        {
          # print(paste("int is ",int,sep=""))
          print(paste("numCoefficients is ",numCoefficients,sep=""))
          print(summary(fit))
          return(na.matrix)
        }
        out1 <- rbind(fit$coefficients,vcov(fit))
      }
      if(int==FALSE & cvar[1]==''){
        fit <- lm(unlist(data1[yvar]) ~ unlist(data1[avar]) + unlist(data1[mvar]))
        if(nrow(summary(fit)$coefficients) != numCoefficients)
        {
          # print(paste("int is ",int,sep=""))
          print(paste("numCoefficients is ",numCoefficients,sep=""))
          print(summary(fit))
          return(na.matrix)
        }
        out1 <- rbind(fit$coefficients,vcov(fit))
      }
      if(int==TRUE & cvar[1]!=''){
        tempdata <- cbind(data1[avar],data1[mvar],data1[cvar],data1["int"])
        fit <- lm(unlist(data1[yvar]) ~ .,data=tempdata)
        if(nrow(summary(fit)$coefficients) != numCoefficients)
        {
          # print(paste("int is ",int,sep=""))
          print(paste("numCoefficients is ",numCoefficients,sep=""))
          print(summary(fit))
          return(na.matrix)
        }
        out1 <- rbind(fit$coefficients,vcov(fit))
      }
      if(int==TRUE & cvar[1]==''){
        fit <- lm(unlist(data1[yvar]) ~ unlist(data1[avar]) + unlist(data1[mvar]) + unlist(data1["int"]))
        if(nrow(summary(fit)$coefficients) != numCoefficients)
        {
          # print(paste("int is ",int,sep=""))
          print(paste("numCoefficients is ",numCoefficients,sep=""))
          print(summary(fit))
          return(na.matrix)
        }
        out1 <- rbind(fit$coefficients,vcov(fit))
      }
    } # end of if(yreg=="linear")
    #************************************************************************************************************************
    if(yreg=="logistic" | yreg=="loglinear" | yreg=="poisson" | yreg=="negbin"){
      if(int==TRUE & cvar[1]!=''){
        if(yreg=="logistic"){
          tempdata <- cbind(data1[avar],data1[mvar],data1[cvar],data1["int"])
          fit <- glm(as.factor(unlist(data1[yvar])) ~ .,data=tempdata,family=binomial)

          if(nrow(summary(fit)$coefficients) != numCoefficients)
          {
            # print(paste("int is ",int,sep=""))
            print(paste("numCoefficients is ",numCoefficients,sep=""))
            print(summary(fit))
            return(na.matrix)
          }
          out1 <- rbind(fit$coefficients,vcov(fit))
        }
        if(yreg=="loglinear"){
          tempdata <- cbind(data1[avar],data1[mvar],data1[cvar],data1["int"])
          fit <- glm(unlist(data1[yvar]) ~ .,family=poisson(link = "log"),data=tempdata)
          
          if(nrow(summary(fit)$coefficients) != numCoefficients)
          {
            # print(paste("int is ",int,sep=""))
            print(paste("numCoefficients is ",numCoefficients,sep=""))
            print(summary(fit))
            return(na.matrix)
          }
          cov = vcov(fit)
          gmparms <- summary(fit)$coefficients
          gmparms <- as.data.frame(gmparms)
          par <- gmparms$Estimate
          out1 <- rbind(par,cov)
        }
        if(yreg=="poisson"){
          tempdata <- cbind(data1[avar],data1[mvar],data1[cvar],data1["int"])
          fit <- glm(unlist(data1[yvar]) ~ .,family=poisson(link = "log"),data=tempdata)

          if(nrow(summary(fit)$coefficients) != numCoefficients)
          {
            # print(paste("int is ",int,sep=""))
            print(paste("numCoefficients is ",numCoefficients,sep=""))
            print(summary(fit))
            return(na.matrix)
          }
          cov = vcov(fit)
          gmparms <- summary(fit)$coefficients
          gmparms <- as.data.frame(gmparms)
          par <- gmparms$Estimate
          out1 <- rbind(par,cov)
        }
        if(yreg=="negbin"){
          tempdata <- cbind(data1[avar],data1[mvar],data1[cvar],data1["int"])
          fit <- glm.nb(unlist(data1[yvar]) ~ .,data=tempdata)

          if(nrow(summary(fit)$coefficients) != numCoefficients)
          {
            # print(paste("int is ",int,sep=""))
            print(paste("numCoefficients is ",numCoefficients,sep=""))
            print(summary(fit))
            return(na.matrix)
          }
          cov = vcov(fit)
          gmparms <- summary(fit)$coefficients
          gmparms <- as.data.frame(gmparms)
          par <- gmparms$Estimate
          out1 <- rbind(par,cov)
        }
      }
      if(int==TRUE & cvar[1]==''){
        if(yreg=="logistic"){
          fit <- glm(as.factor(unlist(data1[yvar])) ~ unlist(data1[avar]) + unlist(data1[mvar]) 
                     + unlist(data1["int"]),family=binomial)
          if(nrow(summary(fit)$coefficients) != numCoefficients)
          {
            # print(paste("int is ",int,sep=""))
            print(paste("numCoefficients is ",numCoefficients,sep=""))
            print(summary(fit))
            return(na.matrix)
          }
          out1 <- rbind(fit$coefficients,vcov(fit))
        }
        if(yreg=="loglinear"){
          fit <- glm(unlist(data1[yvar]) ~ unlist(data1[avar]) + unlist(data1[mvar]) 
                     + unlist(data1["int"]),family=poisson(link = "log"))
          if(nrow(summary(fit)$coefficients) != numCoefficients)
          {
            # print(paste("int is ",int,sep=""))
            print(paste("numCoefficients is ",numCoefficients,sep=""))
            print(summary(fit))
            return(na.matrix)
          }
          cov = vcov(fit)
          gmparms <- summary(fit)$coefficients
          gmparms <- as.data.frame(gmparms)
          par <- gmparms$Estimate
          out1 <- rbind(par,cov)
        }
        if(yreg=="poisson"){
          fit <- glm(unlist(data1[yvar]) ~ unlist(data1[avar]) + unlist(data1[mvar]) 
                     + unlist(data1["int"]),family=poisson(link = "log"))
          if(nrow(summary(fit)$coefficients) != numCoefficients)
          {
            # print(paste("int is ",int,sep=""))
            print(paste("numCoefficients is ",numCoefficients,sep=""))
            print(summary(fit))
            return(na.matrix)
          }
          cov = vcov(fit)
          gmparms <- summary(fit)$coefficients
          gmparms <- as.data.frame(gmparms)
          par <- gmparms$Estimate
          out1 <- rbind(par,cov)
        }
        if(yreg=="negbin"){
          fit <- glm.nb(unlist(data1[yvar]) ~ unlist(data1[avar]) + unlist(data1[mvar]) + unlist(data1["int"]))
          if(nrow(summary(fit)$coefficients) != numCoefficients)
          {
            # print(paste("int is ",int,sep=""))
            print(paste("numCoefficients is ",numCoefficients,sep=""))
            print(summary(fit))
            return(na.matrix)
          }
          cov = vcov(fit)
          gmparms <- summary(fit)$coefficients
          gmparms <- as.data.frame(gmparms)
          par <- gmparms$Estimate
          out1 <- rbind(par,cov)
        }
      }
      if(int==FALSE & cvar[1]!=''){
        if(yreg=="logistic"){
          tempdata <- cbind(data1[avar],data1[mvar],data1[cvar])
          fit <- glm(as.factor(unlist(data1[yvar])) ~ .,data=tempdata,family=binomial)
          
          if(nrow(summary(fit)$coefficients) != numCoefficients)
          {
            # print(paste("int is ",int,sep=""))
            print(paste("numCoefficients is ",numCoefficients,sep=""))
            print(summary(fit))
            return(na.matrix)
          }
          out1 <- rbind(fit$coefficients,vcov(fit))
        }
        if(yreg=="loglinear"){
          tempdata <- cbind(data1[avar],data1[mvar],data1[cvar])
          fit <- glm(unlist(data1[yvar]) ~ .,family=poisson(link = "log"),data=tempdata)
          
          if(nrow(summary(fit)$coefficients) != numCoefficients)
          {
            # print(paste("int is ",int,sep=""))
            print(paste("numCoefficients is ",numCoefficients,sep=""))
            print(summary(fit))
            return(na.matrix)
          }
          cov = vcov(fit)
          gmparms <- summary(fit)$coefficients
          gmparms <- as.data.frame(gmparms)
          par <- gmparms$Estimate
          out1 <- rbind(par,cov)
        }
        if(yreg=="poisson"){
          tempdata <- cbind(data1[avar],data1[mvar],data1[cvar])
          fit <- glm(unlist(data1[yvar]) ~ .,family=poisson(link = "log"),data=tempdata)

          if(nrow(summary(fit)$coefficients) != numCoefficients)
          {
            # print(paste("int is ",int,sep=""))
            print(paste("numCoefficients is ",numCoefficients,sep=""))
            print(summary(fit))
            return(na.matrix)
          }
          cov = vcov(fit)
          gmparms <- summary(fit)$coefficients
          gmparms <- as.data.frame(gmparms)
          par <- gmparms$Estimate
          out1 <- rbind(par,cov)
        }
        if(yreg=="negbin"){
          tempdata <- cbind(data1[avar],data1[mvar],data1[cvar])
          fit <- glm.nb(unlist(data1[yvar]) ~ .,data=tempdata)

          if(nrow(summary(fit)$coefficients) != numCoefficients)
          {
            # print(paste("int is ",int,sep=""))
            print(paste("numCoefficients is ",numCoefficients,sep=""))
            print(summary(fit))
            return(na.matrix)
          }
          cov = vcov(fit)
          gmparms <- summary(fit)$coefficients
          gmparms <- as.data.frame(gmparms)
          par <- gmparms$Estimate
          out1 <- rbind(par,cov)
        }
      }
      if(int==FALSE & cvar[1]==''){
        if(yreg=="logistic"){
          fit <- glm(as.factor(unlist(data1[yvar])) ~ unlist(data1[avar]) + unlist(data1[mvar]),family=binomial)
          if(nrow(summary(fit)$coefficients) != numCoefficients)
          {
            # print(paste("int is ",int,sep=""))
            print(paste("numCoefficients is ",numCoefficients,sep=""))
            print(summary(fit))
            return(na.matrix)
          }
          out1 <- rbind(fit$coefficients,vcov(fit))
        }
        if(yreg=="loglinear"){
          fit <- glm(unlist(data1[yvar]) ~ unlist(data1[avar]) + unlist(data1[mvar]),family=poisson(link = "log"))
          if(nrow(summary(fit)$coefficients) != numCoefficients)
          {
            # print(paste("int is ",int,sep=""))
            print(paste("numCoefficients is ",numCoefficients,sep=""))
            print(summary(fit))
            return(na.matrix)
          }
          cov = vcov(fit)
          gmparms <- summary(fit)$coefficients
          gmparms <- as.data.frame(gmparms)
          par <- gmparms$Estimate
          out1 <- rbind(par,cov)
        }
        if(yreg=="poisson"){
          fit <- glm(unlist(data1[yvar]) ~ unlist(data1[avar]) + unlist(data1[mvar]),family=poisson(link = "log"))
          if(nrow(summary(fit)$coefficients) != numCoefficients)
          {
            # print(paste("int is ",int,sep=""))
            print(paste("numCoefficients is ",numCoefficients,sep=""))
            print(summary(fit))
            return(na.matrix)
          }
          cov = vcov(fit)
          gmparms <- summary(fit)$coefficients
          gmparms <- as.data.frame(gmparms)
          par <- gmparms$Estimate
          out1 <- rbind(par,cov)
        }
        if(yreg=="negbin"){
          fit <- glm.nb(unlist(data1[yvar]) ~ unlist(data1[avar]) + unlist(data1[mvar]))
          if(nrow(summary(fit)$coefficients) != numCoefficients)
          {
            # print(paste("int is ",int,sep=""))
            print(paste("numCoefficients is ",numCoefficients,sep=""))
            print(summary(fit))
            return(na.matrix)
          }
          cov = vcov(fit)
          gmparms <- summary(fit)$coefficients
          gmparms <- as.data.frame(gmparms)
          par <- gmparms$Estimate
          out1 <- rbind(par,cov)
        }
      }
    } # end of if(yreg=="logistic"  | yreg=="loglinear"  |yreg==poisson | yreg=="negbin")
    #**************************************************************************************************************
    if(mreg=="linear" & yreg=="linear"){
      if(cvar[1]!=''){
        if(casecontrol!=TRUE){
          tempdata <- cbind(data1[avar],data1[cvar])
          fit <- lm(unlist(data1[mvar]) ~ .,data=tempdata)
          if(nrow(summary(fit)$coefficients) != (numCoefficients - (int*1 +1)))
          {
            # print(paste("int is ",int,sep=""))
            print(paste("numCoefficients is ",numCoefficients,sep=""))
            print(summary(fit))
            return(na.matrix)
          }
          out2 <- rbind(fit$coefficients,vcov(fit))
        }
        if(casecontrol==TRUE){
          nums <- which(data1[yvar]==0)
          newset <- c()
          for(s in nums){
            newset <- rbind(newset,data1[s,1:ncol(data1)])
          }
          tempdata <- cbind(newset[avar],newset[cvar])
          fit <- lm(unlist(newset[mvar]) ~ .,data=tempdata)
          if(nrow(summary(fit)$coefficients) != (numCoefficients - (int*1 +1)))
          {
            # print(paste("int is ",int,sep=""))
            print(paste("numCoefficients is ",numCoefficients,sep=""))
            print(summary(fit))
            return(na.matrix)
          }
          out2 <- rbind(fit$coefficients,vcov(fit))
        }
      }
      if(cvar[1]==''){
        if(casecontrol!=TRUE){
          fit <- lm(unlist(data1[mvar]) ~ unlist(data1[avar]))
          if(nrow(summary(fit)$coefficients) != (numCoefficients - (int*1 +1)))
          {
            # print(paste("int is ",int,sep=""))
            print(paste("numCoefficients is ",numCoefficients,sep=""))
            print(summary(fit))
            return(na.matrix)
          }
          out2 <- rbind(fit$coefficients,vcov(fit))
        }
        if(casecontrol==TRUE){
          nums <- which(data1[yvar]==0)
          newset <- c()
          for(s in nums){
            newset <- rbind(newset,data1[s,1:ncol(data1)])
          }
          fit <- lm(unlist(newset[mvar]) ~ unlist(newset[avar]))
          if(nrow(summary(fit)$coefficients) != (numCoefficients - (int*1 +1)))
          {
            # print(paste("int is ",int,sep=""))
            print(paste("numCoefficients is ",numCoefficients,sep=""))
            print(summary(fit))
            return(na.matrix)
          }
          out2 <- rbind(fit$coefficients,vcov(fit))
        }
      }
    } # end of if(mreg=="linear" & yreg=="linear")
    #************************************************************************************************************************;
    if(mreg=="linear" & yreg!="linear")
    {
      #************************************************************************************************************************;
      if(cvar[1]!=''){
        if(casecontrol!=TRUE){
          tempdata <- cbind(data1[avar],data1[cvar])
          fit <- lm(unlist(data1[mvar]) ~ .,data=tempdata)
          if(nrow(summary(fit)$coefficients) != (numCoefficients - (int*1 +1)))
          {
            # print(paste("int is ",int,sep=""))
            print(paste("numCoefficients is ",numCoefficients,sep=""))
            print(summary(fit))
            return(na.matrix)
          }
          temp <- rbind(fit$coefficients,vcov(fit))
          RMSE <- rep(summary(fit)$sigma,nrow(temp))
          out2 <- cbind(RMSE,temp)
        }
        if(casecontrol==TRUE){
          nums <- which(data1[yvar]==0)
          newset <- c()
          for(s in nums){
            newset <- rbind(newset,data1[s,1:ncol(data1)])
          }
          tempdata <- cbind(newset[avar],newset[cvar])
          fit <- lm(unlist(newset[mvar]) ~ .,data=tempdata)
          if(nrow(summary(fit)$coefficients) != (numCoefficients - (int*1 +1)))
          {
            # print(paste("int is ",int,sep=""))
            print(paste("numCoefficients is ",numCoefficients,sep=""))
            print(summary(fit))
            return(na.matrix)
          }
          temp <- rbind(fit$coefficients,vcov(fit))
          RMSE <- rep(summary(fit)$sigma,nrow(temp))
          out2 <- cbind(RMSE,temp)
        }
      }
      if(cvar[1]==''){
        if(casecontrol!=TRUE){
          fit <- lm(unlist(data1[mvar]) ~ unlist(data1[avar]))
          if(nrow(summary(fit)$coefficients) != (numCoefficients - (int*1 +1)))
          {
            # print(paste("int is ",int,sep=""))
            print(paste("numCoefficients is ",numCoefficients,sep=""))
            print(summary(fit))
            return(na.matrix)
          }
          temp <- rbind(fit$coefficients,vcov(fit))
          RMSE <- rep(summary(fit)$sigma,nrow(temp))
          out2 <- cbind(RMSE,temp)
        }
        if(casecontrol==TRUE){
          nums <- which(data1[yvar]==0)
          newset <- c()
          for(s in nums){
            newset <- rbind(newset,data1[s,1:ncol(data1)])
          }
          fit <- lm(unlist(newset[mvar]) ~ unlist(newset[avar]))
          if(nrow(summary(fit)$coefficients) != (numCoefficients - (int*1 +1)))
          {
            # print(paste("int is ",int,sep=""))
            print(paste("numCoefficients is ",numCoefficients,sep=""))
            print(summary(fit))
            return(na.matrix)
          }
          temp <- rbind(fit$coefficients,vcov(fit))
          RMSE <- rep(summary(fit)$sigma,nrow(temp))
          out2 <- cbind(RMSE,temp)
        }
      }
    } # end of if(mreg=="linear" & yreg!="linear")
    #***********************************************************************************************
    if(mreg=="logistic"){
      if(cvar[1]!=''){
        if(casecontrol!=TRUE){
          tempdata <- cbind(data1[avar],data1[cvar])
          fit <- glm(as.factor(unlist(data1[mvar])) ~ .,data=tempdata,family=binomial)

          if(nrow(summary(fit)$coefficients) != (numCoefficients - (int*1 +1)))
          {
            # print(paste("int is ",int,sep=""))
            print(paste("numCoefficients is ",numCoefficients,sep=""))
            print(summary(fit))
            return(na.matrix)
          }
          out2 <- rbind(fit$coefficients,vcov(fit))
        }
        if(casecontrol==TRUE){
          nums <- which(data1[yvar]==0)
          newset <- c()
          for(s in nums){
            newset <- rbind(newset,data1[s,1:ncol(data1)])
          }
          tempdata <- cbind(newset[avar],newset[cvar])
          fit <- glm(as.factor(unlist(newset[mvar])) ~ .,data=tempdata,family=binomial)

          if(nrow(summary(fit)$coefficients) != (numCoefficients - (int*1 +1)))
          {
            # print(paste("int is ",int,sep=""))
            print(paste("numCoefficients is ",numCoefficients,sep=""))
            print(summary(fit))
            return(na.matrix)
          }
          out2 <- rbind(fit$coefficients,vcov(fit))
        }
      }
      if(cvar[1]==''){
        if(casecontrol!=TRUE){
          fit <- glm(as.factor(unlist(data1[mvar])) ~ unlist(data1[avar]),family=binomial)
          if(nrow(summary(fit)$coefficients) != (numCoefficients - (int*1 +1)))
          {
            # print(paste("int is ",int,sep=""))
            print(paste("numCoefficients is ",numCoefficients,sep=""))
            print(summary(fit))
            return(na.matrix)
          }
          out2 <- rbind(fit$coefficients,vcov(fit))
        }
        if(casecontrol==TRUE){
          nums <- which(data1[yvar]==0)
          newset <- c()
          for(s in nums){
            newset <- rbind(newset,data1[s,1:ncol(data1)])
          }
          fit <- glm(as.factor(unlist(newset[mvar])) ~ unlist(newset[avar]),family=binomial)
          if(nrow(summary(fit)$coefficients) != (numCoefficients - (int*1 +1)))
          {
            # print(paste("int is ",int,sep=""))
            print(paste("numCoefficients is ",numCoefficients,sep=""))
            print(summary(fit))
            return(na.matrix)
          }
          out2 <- rbind(fit$coefficients,vcov(fit))
        }
      }
    } # end of if(mreg=="logistic")
    # *******************************************************************
    
    # ***************************   OUTPUT    *************************** 
    
    #*/"linear" "linear" no int 
    if((mreg=="linear" & int==FALSE )|(yreg=="linear" & mreg=="logistic" & int==FALSE & cvar[1]=='')){
      # pvalue <- as.matrix(pvalue)
      pvalue <- t(pvalue)
      if((boot=='' | boot==FALSE)){
        # pvalue <- data.matrix(pvalue)
        cil <- data.matrix(cil)
        ciu <- data.matrix(ciu)
      }
      if((boot!='' & boot!=FALSE)){
        # ci <- data.matrix(ci)
        x=cbind(t(effect),se,ci)
        cname1=c("Estimate","s.e.","95% CI lower","95% CI upper")
      }
      x3 = x
      colnames(x3) <- cname1
      name = c("cde=nde","nie","total effect")
      name=t(t(name))
      cname2=c("Effect")
      x4 = name
      colnames(x4) <- cname2
    } # End of if((mreg=="linear" & int==FALSE )|(yreg=="linear" & .... )
    
    #*/"linear" "linear" int 
    if((int==TRUE)|(mreg=="logistic" & int==FALSE & cvar[1]!='')|(yreg!="linear" & mreg=="logistic" & int==FALSE)){
      # pvalue <- as.matrix(pvalue)
      pvalue <- t(pvalue)
      if(output=="full" & c!='' & cvar[1]!=''){
        name = c("marginal cde" , "marginal pnde","marginal pnie","marginal tnde","marginal tnie","marginal total effect","conditional cde", 
                 "conditional pnde","conditional pnie","conditional tnde","conditional tnie","conditional total effect")
        name=t(t(name))
        cname2=("effect")
        x4 = name
        colnames(x4) <- cname2
        if((boot!='' & boot!=FALSE)){
          x=cbind(t(effect),se,ci)
          cname1 = c("Estimate","s.e.","95% CI lower","95% CI upper")
        }
        x3 = x
        colnames(x3) <- cname1
      }
      if((output=="full" & cvar[1]=='') | (output=="full" & cvar[1]!='' & c=='')){
        name = c("cde", "pnde","pnie","tnde","tnie","total effect")
        name = t(t(name))
        cname2 = c("effect")
        x4 = name
        colnames(x4) <- cname2
        if((boot!='' & boot!=FALSE)){
          x = cbind(t(effect),se,ci)
          cname1 = c("Estimate","s.e.","95% CI lower","95% CI upper")
        }
        x3 = x
        colnames(x3) <- cname1
      }
    } # *end int
    
    x5 <- cbind(x4,x3)
    if (ncol(x5)==5){
      x5 <- cbind(x5,pvalue)
    }
    
    perc.name.low <- paste((1-alphalev)*100,"% CI lower",sep = "")
    perc.name.upp <- paste((1-alphalev)*100,"% CI upper",sep = "")

    colnames(x5) <- c("Effect","Estimate","Standard Error",perc.name.low,perc.name.upp,"p-value")
    if((int==TRUE)|(mreg=="logistic" & int==FALSE & cvar[1]!='')|(yreg!="linear" & mreg=="logistic" & int==FALSE)){
      rownames(x5) <- c("cde","pnde","pnie","tnde","tnie","total effect")
    }
    if((mreg=="linear" & int==FALSE )|(yreg=="linear" & mreg=="logistic" & int==FALSE & cvar[1]=='')){
      rownames(x5) <- c("cde=nde","nie","total effect")
    }
    
    if(cvar!=""){
      if(int==TRUE & yreg!="linear"){
        if(mreg=="linear"){
          x5[c("pnde","tnde","total effect"),] <- NA
        }
        if(mreg=="logistic"){
          x5[c("pnde","pnie","tnde","tnie","total effect"),] <- NA
        }
      }
    }

    x5 <- as.data.frame(x5)
    x5$Effect <- NULL
    return(x5)
    
    # ***************************   OUTPUT -END-   *************************** 
    
  }

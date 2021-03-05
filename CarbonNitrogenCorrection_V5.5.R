###Please note that this code only deals with carbon labeling experiment.
###Please use other version of correction codes for other experiments.
#Please make sure you have the packages below installed.### 
#Use install.packages("") to install packages###

require(gsubfn)
require(nnls)
require(dplyr)
require(gdata)
require(xlsx)
require(CHNOSZ)
require(stringr)

CarbonNaturalAbundace <- c(0.9893, 0.0107)
HydrogenNaturalAbundace <- c(0.999885, 0.000115)
NitrogenNaturalAbundace <- c(0.99636, 0.00364)
OxygenNaturalAbundace <- c(0.99757, 0.00038, 0.00205)
SulfurNaturalAbundace <- c(0.95, 0.0075, 0.0425)
SiliconNaturalAbundace <- c(0.92223, 0.04685, 0.03092)
ChlorineNaturalAbundance<- c(0.7576,0.2424)
BromineNaturalAbundance<- c(0.5069,0.4931)

###Please make sure these parameters are accurate.#
C13Purity <- 0.99
N15Purity <- 0.99
Resolution <- 70000
ResDefAt <- 200
ReportPoolSize <- TRUE
#For Exactive, the Resolution is 100000, defined at Mw 200#

#Please make sure the paths below are correct. 
#Make sure to use / in the paths

InputFile <- "CN_Compound_Test.xlsx"
InputSheetName <- "Sheet1"
MetaboliteList <- read.csv("Metabolite Formula and Charge Info.csv", header = TRUE, check.names = FALSE,stringsAsFactors=FALSE)

InputDF <- read.xlsx(InputFile, InputSheetName, header = TRUE, check.names=FALSE, stringsAsFactors=FALSE)
InputDF[,1] <- as.character(InputDF[,1])
OutputMatrix <- matrix(0, ncol=(ncol(InputDF)-5), nrow=0)
OutputPercentageMatrix <- matrix(0, ncol=(ncol(InputDF)-5), nrow=0)
OutputPoolBefore <- matrix(0, ncol=(ncol(InputDF)-5), nrow=0)
OutputPoolAfter <- matrix(0, ncol=(ncol(InputDF)-5), nrow=0)
# OutputEnrichment <- matrix(0, ncol=(ncol(InputDF)-5), nrow=0)
OutputCompound <- NULL
OutputLabel <- NULL
OutputPoolCompound <- NULL
names(InputDF)[1] <- "Compound"

if.not.null <- function(x) if(!is.null(x)) x else 0

for (i in 1:nrow(InputDF)) {
  if(startsWith(InputDF[i,1], "Unknown")) {InputDF[i,1] <- InputDF[i-1,1]}
}

Erh <- function(x) {
  n <- length(x)-1
  return(sum(x*c(0:n))/n)
}

#########################

IsotopeCorrection <- function(formula, datamatrix, label,Autofill = F, Charge = -1) {
  if (Charge == 0){
    print("Charge cannot be 0, replaced with -1")
    Charge <- -1
  }
  AtomNumber <- rep(0,9)
  names(AtomNumber) <- c("C","H","N","O","P","S","Si","Cl","Br")
  AtomicComposition <- makeup(formula)
  # adjust the H# according to the charge
  AtomicComposition["H"]<-AtomicComposition["H"]+Charge
  
  for (i in names(AtomicComposition)) {
    AtomNumber[i] <- AtomicComposition[i]
  }
  
  MZ <- sum(AtomNumber*c(12,1,14,16,31,32,28,35.5,80))/abs(Charge)
  Mass.Limit <- 1.66*MZ^1.5/Resolution/sqrt(ResDefAt)
  
  
  Mass.Difference<-(label$`13C#`)*1.00628+(label$`15N#`)*0.99703
  for (i in 1:(length(Mass.Difference)-1)){
    for (j in (i+1):length(Mass.Difference)){
      if (Mass.Limit>abs(Mass.Difference[i]-Mass.Difference[j])){
        print(paste(compound,": 13C",label[i,1],"-15N",label[i,2]," is indistinguishable with 13C",label[j,1],"-15N",label[j,2]," under provided resolution, check your data again.",sep = ""))
      }
    }
  }
  
  CompleteLabel<-expand.grid(C13=c(0:AtomNumber["C"]),N15=c(0:AtomNumber["N"]))
  ExpMatrix <- matrix(0,ncol=ncol(datamatrix),nrow=nrow(CompleteLabel))
  LabelIndex <-c()
  for (i in 1:nrow(label)){
    for (j in 1:nrow(CompleteLabel)){
      if (prod(CompleteLabel[j,]==label[i,])){
        ExpMatrix[j,]<-datamatrix[i,]
        LabelIndex<-append(LabelIndex,j)
      }
    }
  }
  CorrectedMatrix <- matrix(0, ncol=ncol(datamatrix), nrow=nrow(CompleteLabel))
  
  if((AtomNumber["C"]+1)*(AtomNumber["N"]+1) < length(label)){
    if (length(compound)!=0){
      print(paste(compound,":the number of labeling species exceeded the maximum allowed number of carbon and Nitrogen isotopomers, check the input file."))
    } else{
      print(paste(formula,":the number of labeling species exceeded the maximum allowed number of carbon and Nitrogen isotopomers, check the input file."))
    }
    return(CorrectedMatrix)
  } else{
    
    #Construct CNMatrix
    CNMatrix <- diag((AtomNumber["C"]+1)*(AtomNumber["N"]+1))
    CMatrix <- diag((AtomNumber["C"]+1))
    NMatrix <- diag((AtomNumber["N"]+1))
    for(i in 1:(AtomNumber["C"]+1)){
      CMatrix[,i] <- sapply(0:AtomNumber["C"], function(x) dbinom(x-i+1, AtomNumber["C"]-i+1 , CarbonNaturalAbundace[2]))
    }
    for(i in 1:(AtomNumber["N"]+1)){
      NMatrix[,i] <- sapply(0:AtomNumber["N"], function(x) dbinom(x-i+1, AtomNumber["N"]-i+1 , NitrogenNaturalAbundace[2]))
    }
    for(m in 1:(AtomNumber["N"]+1)){
      for (n in 1:m){
        CNMatrix[((m-1)*(AtomNumber["C"]+1)+1):(m*(AtomNumber["C"]+1)),((n-1)*(AtomNumber["C"]+1)+1):(n*(AtomNumber["C"]+1))] <- CMatrix*NMatrix[m,n]
      }
    }

    #Construct PurityMatrix
    PurityMatrix <- diag((AtomNumber["C"]+1)*(AtomNumber["N"]+1))
    CPurityMatrix <- diag((AtomNumber["C"]+1))
    NPurityMatrix <- diag((AtomNumber["N"]+1))
    for(i in 1:(AtomNumber["C"]+1)){
      CPurityMatrix[i,] <- sapply(0:AtomNumber["C"], function(x) dbinom(x-i+1, x , 1-C13Purity))
    }
    for(i in 1:(AtomNumber["N"]+1)){
      NPurityMatrix[i,] <- sapply(0:AtomNumber["N"], function(x) dbinom(x-i+1, x , 1-N15Purity))
    }
    for(n in 1:(AtomNumber["N"]+1)){
      for (m in 1:n){
        PurityMatrix[((m-1)*(AtomNumber["C"]+1)+1):(m*(AtomNumber["C"]+1)),((n-1)*(AtomNumber["C"]+1)+1):(n*(AtomNumber["C"]+1))] <- CPurityMatrix*NPurityMatrix[m,n]
      }
    }

    #Construct NonTracerMatrix  
    Isotope.Combinations <- expand.grid(C13=c(0:AtomNumber["C"]),N15=c(0:AtomNumber["N"]),H2=c(0:AtomNumber["H"]),O17=c(0:AtomNumber["O"]),
                                        O18=c(0:AtomNumber["O"]),S33=c(0:AtomNumber["S"]),S34=c(0:AtomNumber["S"]),
                                        Si29=c(0:AtomNumber["Si"]),Si30=c(0:AtomNumber["Si"]),Cl37=c(0:AtomNumber["Cl"]),
                                        Br81=c(0:AtomNumber["Br"]))
  
    Isotope.Combinations1 <- Isotope.Combinations %>% mutate(NonTracerMass=H2+O17+O18*2+S33+S34*2+Si29+Si30*2+Cl37*2+Br81*2) %>%
      filter((O17+O18)<=AtomNumber["O"] & (S33+S34)<=AtomNumber["S"] & (Si29+Si30)<=AtomNumber["Si"] & NonTracerMass<=(AtomNumber["C"]+AtomNumber["N"])) %>% filter((NonTracerMass==C13+N15))%>%
      mutate(MassDiff=1.00628*H2+1.00422*O17+2.00425*O18+0.99939*S33+1.99580*S34+0.99957*Si29+1.99684*Si30+1.99705*Cl37+1.99796*Br81
             -1.00335*C13-0.99703*N15) %>%
      filter(abs(MassDiff/Charge)<Mass.Limit)
    
    Isotope.Combinations2 <- Isotope.Combinations %>% mutate(NonTracerMass=H2+O17+O18*2+S33+S34*2+Si29+Si30*2+Cl37*2+Br81*2) %>%
      filter((O17+O18)<=AtomNumber["O"] & (S33+S34)<=AtomNumber["S"] & (Si29+Si30)<=AtomNumber["Si"] & NonTracerMass<=(AtomNumber["C"]+AtomNumber["N"])) %>% filter((NonTracerMass==abs(C13-N15))&(C13*N15*NonTracerMass>0))%>%
      mutate(MassDiff=1.00628*H2+1.00422*O17+2.00425*O18+0.99939*S33+1.99580*S34+0.99957*Si29+1.99684*Si30+1.99705*Cl37+1.99796*Br81
             -abs(1.00335*C13-0.99703*N15)) %>%
      filter(abs(MassDiff/Charge)<Mass.Limit)

    NonTracerMatrix <- matrix(0,ncol=(AtomNumber["C"]+1)*(AtomNumber["N"]+1),nrow=(AtomNumber["C"]+1)*(AtomNumber["N"]+1))
    for(m in 1:nrow(NonTracerMatrix)){
      Cindex<-CompleteLabel[m,1]
      Nindex<-CompleteLabel[m,2]
      Current.Combinations<-Isotope.Combinations1%>%filter(C13==Cindex & N15==Nindex)
      if (nrow(Current.Combinations)!=0){
        for (i in 1:nrow(Current.Combinations)){
          p <- dbinom(Current.Combinations[i,3], AtomNumber["H"], HydrogenNaturalAbundace[2])*
            dmultinom(unlist(c(AtomNumber["O"]-sum(Current.Combinations[i,4:5]),Current.Combinations[i,4:5])), AtomNumber["O"], OxygenNaturalAbundace)*
            dmultinom(unlist(c(AtomNumber["S"]-sum(Current.Combinations[i,6:7]),Current.Combinations[i,6:7])), AtomNumber["S"], SulfurNaturalAbundace)*
            dmultinom(unlist(c(AtomNumber["Si"]-sum(Current.Combinations[i,8:9]),Current.Combinations[i,8:9])), AtomNumber["Si"], SiliconNaturalAbundace)*
            dbinom(Current.Combinations[i,10], AtomNumber["Cl"], ChlorineNaturalAbundance[2])*
            dbinom(Current.Combinations[i,11], AtomNumber["Br"], BromineNaturalAbundance[2])
          for (n in m:nrow(NonTracerMatrix)){
            if ((n-m)%%(AtomNumber["C"]+1)<=(n-1)%%(AtomNumber["C"]+1)){
              NonTracerMatrix[n,n-m+1] <- NonTracerMatrix[n,n-m+1]+p
            }
          }
        }
      }
    }
    #Correct the situations such as 2H+15N -> 13C2
    Current.Combinations<-Isotope.Combinations2
    if (nrow(Current.Combinations)!=0){
      for (i in 1:nrow(Current.Combinations)){
        Cindex<- Current.Combinations[i,1]
        Nindex<- Current.Combinations[i,2]
        p <- dbinom(Current.Combinations[i,3], AtomNumber["H"], HydrogenNaturalAbundace[2])*
          dmultinom(unlist(c(AtomNumber["O"]-sum(Current.Combinations[i,4:5]),Current.Combinations[i,4:5])), AtomNumber["O"], OxygenNaturalAbundace)*
          dmultinom(unlist(c(AtomNumber["S"]-sum(Current.Combinations[i,6:7]),Current.Combinations[i,6:7])), AtomNumber["S"], SulfurNaturalAbundace)*
          dmultinom(unlist(c(AtomNumber["Si"]-sum(Current.Combinations[i,8:9]),Current.Combinations[i,8:9])), AtomNumber["Si"], SiliconNaturalAbundace)*
          dbinom(Current.Combinations[i,10], AtomNumber["Cl"], ChlorineNaturalAbundance[2])*
          dbinom(Current.Combinations[i,11], AtomNumber["Br"], BromineNaturalAbundance[2])
        if (Cindex<Nindex){
          n<-which(CompleteLabel$C13==Cindex&CompleteLabel$N15==0)
          m<-which(CompleteLabel$C13==0&CompleteLabel$N15==Nindex)
          for (j in 0:(AtomNumber["C"]+1-n)){
            for (h in 0: (AtomNumber["N"]-m%/%(AtomNumber["C"]+1))){
              NonTracerMatrix[m+j+h*(AtomNumber["C"]+1),n+j+h*(AtomNumber["C"]+1)] <- NonTracerMatrix[m+j+h*(AtomNumber["C"]+1),n+j+h*(AtomNumber["C"]+1)]+p
            }
          }
        } else {
          n<-which(CompleteLabel$C13==0&CompleteLabel$N15==Nindex)
          m<-which(CompleteLabel$C13==Cindex&CompleteLabel$N15==0)
          for (j in 0:(AtomNumber["C"]+1-m)){
            for (h in 0: (AtomNumber["N"]-n%/%(AtomNumber["C"]+1))){
            NonTracerMatrix[m+j+h*(AtomNumber["C"]+1),n+j+h*(AtomNumber["C"]+1)] <- NonTracerMatrix[m+j+h*(AtomNumber["C"]+1),n+j+h*(AtomNumber["C"]+1)]+p
            }
          }
        }
      }
    }
    #remove matrix rows and columns that autofilled in the input file
    if (Autofill==F){
      LabelIndex<-sort(LabelIndex)
      NonTracerMatrix<-NonTracerMatrix[LabelIndex,LabelIndex]
      CNMatrix<-CNMatrix[LabelIndex,LabelIndex]
      PurityMatrix<-PurityMatrix[LabelIndex,LabelIndex]
      ExpMatrix<-ExpMatrix[LabelIndex,]
      CorrectedMatrix<-CorrectedMatrix[LabelIndex,]
      CompleteLabel<-CompleteLabel[LabelIndex,]
    }
    
    for(i in 1:ncol(CorrectedMatrix)) {
      CorrectedMatrix[,i] <- coef(nnls(NonTracerMatrix %*% 
                                         CNMatrix %*% PurityMatrix, ExpMatrix[,i]))
    }
    CorrectedMatrix<-cbind(CompleteLabel,CorrectedMatrix)
    return(CorrectedMatrix)
  }
}

for (i in unique(InputDF$Compound)) {
  Formula=MetaboliteList[MetaboliteList$compound==i,2]
  charge=MetaboliteList[MetaboliteList$compound==i,3]
  compound <- i
  if(length(Formula)==0) {
    print(paste("The formula of",i,"is unknown",sep=" "))
    break
  }
  CurrentMetabolite <- filter(InputDF, Compound==i)
  DataMatrix <- as.matrix(CurrentMetabolite[,6:(ncol(CurrentMetabolite))])
  Corrected.raw <- IsotopeCorrection(Formula, DataMatrix, CurrentMetabolite[,3:4],Charge = charge, Autofill = F)
  Corrected<-NULL
  #Delete unexpected rows
  for (m in 1:nrow(Corrected.raw)){
    for (n in 1:nrow(CurrentMetabolite)){
      if(prod(Corrected.raw[m,1:2]==CurrentMetabolite[n,3:4]) & CurrentMetabolite[n,5] ){
        Corrected<- rbind(Corrected,Corrected.raw[m,])
      }
    }
  }
  CorrectedPercentage <- scale(Corrected[,3:ncol(Corrected)],scale=colSums(Corrected)[3:ncol(Corrected)],center=FALSE)
  OutputMatrix <- rbind(OutputMatrix, Corrected[,3:ncol(Corrected)])
  OutputPercentageMatrix <- rbind(OutputPercentageMatrix, CorrectedPercentage)
  OutputPoolBefore <- rbind(OutputPoolBefore, colSums(DataMatrix))
  OutputPoolAfter <- rbind(OutputPoolAfter, colSums(Corrected[,3:ncol(Corrected)]))
  OutputCompound <- append(OutputCompound, rep(i,nrow(Corrected[,3:ncol(Corrected)])))
  OutputLabel <- rbind(OutputLabel, Corrected[,1:2])
  OutputPoolCompound <- append(OutputPoolCompound, i)
} 

OutputDF <- data.frame(OutputCompound, OutputLabel, OutputMatrix)
OutputPercentageDF <- data.frame(OutputCompound, OutputLabel, OutputPercentageMatrix)
OutputPoolBeforeDF <- data.frame(OutputPoolCompound, OutputPoolBefore)
OutputPoolAfterDF <- data.frame(OutputPoolCompound, OutputPoolAfter)
names(OutputDF) <- c("Compound", "C13","N15", names(InputDF)[6:(length(names(InputDF)))])
names(OutputPercentageDF) <- names(OutputDF)
names(OutputPoolBeforeDF) <- c("Compound", names(InputDF)[6:(length(names(InputDF)))])
names(OutputPoolAfterDF) <- c("Compound", names(InputDF)[6:(length(names(InputDF)))])

write.xlsx2(OutputDF, file=InputFile, sheetName = "Corrected", row.names=FALSE, append=TRUE)
write.xlsx2(OutputPercentageDF, file=InputFile, sheetName = "Normalized", row.names=FALSE, append=TRUE)
write.xlsx2(OutputPoolAfterDF, file=InputFile, sheetName = "Pool Size", row.names=FALSE, append=TRUE)

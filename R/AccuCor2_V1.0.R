#' Natural Abundance correction for 13C-2H tracer labeling data
#'
#' @param formula string describing the chemical formula of the metabolite (in neutral form)
#' @param datamatrix matrix that contains the observed mass fractions
#' @param label matrix that contains the number of isotopes (C followed by H) associated with the observed mass fractions
#' @param Resolution For Exactive, the Resolution is 100000, defined at Mw 200
#' @param Autofill default: FALSE. If the intensity of an isotopologue is not provided, should 0 be assumed?
#' @param CarbonNaturalAbundace
#' @param HydrogenNaturalAbundace
#' @param NitrogenNaturalAbundace
#' @param SulfurNaturalAbundace
#' @param SiliconNaturalAbundace
#' @param ChlorineNaturalAbundance
#' @param BromineNaturalAbundance
#' @param C13Purity default:0.99.The isotopic purity for C13 tracer.
#' @param H2N15Purity default:0.99. The isotopic purity for H2/N15 tracer. 
#' @param ResDefAt Resolution defined at (in Mw), e.g. 200 Mw
#' @param ReportPoolSize default: TRUE
#' @import xlsx
#' @import gsubfn
#' @import nnls
#' @import dplyr
#' @import gdata
#' @import CHNOSZ
#' @import stringr
#' @return matrix of labeling pattern after the correction
CH_Correction <- function(formula, datamatrix, label, Resolution, 
                          Autofill = F, 
                          Charge = -1,
                          CarbonNaturalAbundace = c(0.9893, 0.0107),
                          HydrogenNaturalAbundace = c(0.999885, 0.000115),
                          NitrogenNaturalAbundace = c(0.99636, 0.00364),
                          OxygenNaturalAbundace = c(0.99757, 0.00038, 0.00205),
                          SulfurNaturalAbundace = c(0.95, 0.0075, 0.0425),
                          SiliconNaturalAbundace = c(0.92223, 0.04685, 0.03092),
                          ChlorineNaturalAbundance = c(0.7576,0.2424),
                          BromineNaturalAbundance = c(0.5069,0.4931),
                          C13Purity = 0.99,
                          H2N15Purity = 0.99,
                          ResDefAt = 200) {
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
  
  
  Mass.Difference<-(label$`13C#`)*1.00628+(label$`2H#`)*1.00335
  for (i in 1:(length(Mass.Difference)-1)){
    for (j in (i+1):length(Mass.Difference)){
      if (Mass.Limit>abs((Mass.Difference[i]-Mass.Difference[j])/Charge)){
        print(paste(compound,": 13C",label[i,1],"-2H",label[i,2]," is indistinguishable with 13C",label[j,1],"-2H",label[j,2]," under provided resolution, check your data again.",sep = ""))
      }
    }
  }
  
  CompleteLabel<-expand.grid(C13=c(0:AtomNumber["C"]),H2=c(0:AtomNumber["H"]))
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
  
  if((AtomNumber["C"]+1)*(AtomNumber["H"]+1) < length(label)){
    if (length(compound)!=0){
      print(paste(compound,":the number of labeling species exceeded the maximum allowed number of carbon and hydrogen isotopomers, check the input file."))
    } else{
      print(paste(formula,":the number of labeling exceeded the number of carbon atoms, check the input file."))
    }
    return(CorrectedMatrix)
  } else{
    
    #Construct CHMatrix
    CHMatrix <- diag((AtomNumber["C"]+1)*(AtomNumber["H"]+1))
    CMatrix <- diag((AtomNumber["C"]+1))
    HMatrix <- diag((AtomNumber["H"]+1))
    for(i in 1:(AtomNumber["C"]+1)){
      CMatrix[,i] <- sapply(0:AtomNumber["C"], function(x) dbinom(x-i+1, AtomNumber["C"]-i+1 , CarbonNaturalAbundace[2]))
    }
    for(i in 1:(AtomNumber["H"]+1)){
      HMatrix[,i] <- sapply(0:AtomNumber["H"], function(x) dbinom(x-i+1, AtomNumber["H"]-i+1 , HydrogenNaturalAbundace[2]))
    }
    for(m in 1:(AtomNumber["H"]+1)){
      for (n in 1:m){
        CHMatrix[((m-1)*(AtomNumber["C"]+1)+1):(m*(AtomNumber["C"]+1)),((n-1)*(AtomNumber["C"]+1)+1):(n*(AtomNumber["C"]+1))] <- CMatrix*HMatrix[m,n]
      }
    }
    
    #Construct PurityMatrix
    PurityMatrix <- diag((AtomNumber["C"]+1)*(AtomNumber["H"]+1))
    CPurityMatrix <- diag((AtomNumber["C"]+1))
    HPurityMatrix <- diag((AtomNumber["H"]+1))
    for(i in 1:(AtomNumber["C"]+1)){
      CPurityMatrix[i,] <- sapply(0:AtomNumber["C"], function(x) dbinom(x-i+1, x , 1-C13Purity))
    }
    for(i in 1:(AtomNumber["H"]+1)){
      HPurityMatrix[i,] <- sapply(0:AtomNumber["H"], function(x) dbinom(x-i+1, x , 1-H2N15Purity))
    }
    for(n in 1:(AtomNumber["H"]+1)){
      for (m in 1:n){
        PurityMatrix[((m-1)*(AtomNumber["C"]+1)+1):(m*(AtomNumber["C"]+1)),((n-1)*(AtomNumber["C"]+1)+1):(n*(AtomNumber["C"]+1))] <- CPurityMatrix*HPurityMatrix[m,n]
      }
    }
    
    #Construct NonTracerMatrix  
    Isotope.Combinations <- expand.grid(C13=c(0:AtomNumber["C"]),H2=c(0:AtomNumber["H"]),N15=c(0:AtomNumber["N"]),O17=c(0:AtomNumber["O"]),
                                        O18=c(0:AtomNumber["O"]),S33=c(0:AtomNumber["S"]),S34=c(0:AtomNumber["S"]),
                                        Si29=c(0:AtomNumber["Si"]),Si30=c(0:AtomNumber["Si"]),Cl37=c(0:AtomNumber["Cl"]),
                                        Br81=c(0:AtomNumber["Br"]))
    
    Isotope.Combinations1 <- Isotope.Combinations %>% mutate(NonTracerMass=N15+O17+O18*2+S33+S34*2+Si29+Si30*2+Cl37*2+Br81*2) %>%
      filter((O17+O18)<=AtomNumber["O"] & (S33+S34)<=AtomNumber["S"] & (Si29+Si30)<=AtomNumber["Si"] & NonTracerMass<=(AtomNumber["C"]+AtomNumber["H"])) %>% filter((NonTracerMass==C13+H2))%>%
      mutate(MassDiff=0.99703*N15+1.00422*O17+2.00425*O18+0.99939*S33+1.99580*S34+0.99957*Si29+1.99684*Si30+1.99705*Cl37+1.99796*Br81
             -1.00335*C13-1.00627*H2) %>%
      filter(abs(MassDiff/Charge)<Mass.Limit)
    
    Isotope.Combinations2 <- Isotope.Combinations %>% mutate(NonTracerMass=N15+O17+O18*2+S33+S34*2+Si29+Si30*2+Cl37*2+Br81*2) %>%
      filter((O17+O18)<=AtomNumber["O"] & (S33+S34)<=AtomNumber["S"] & (Si29+Si30)<=AtomNumber["Si"] & NonTracerMass<=(AtomNumber["C"]+AtomNumber["H"])) %>% filter((NonTracerMass==abs(C13-H2))&(C13*H2*NonTracerMass>0))%>%
      mutate(MassDiff=0.99703*N15+1.00422*O17+2.00425*O18+0.99939*S33+1.99580*S34+0.99957*Si29+1.99684*Si30+1.99705*Cl37+1.99796*Br81
             -abs(1.00335*C13-1.00627*H2)) %>%
      filter(abs(MassDiff/Charge)<Mass.Limit)
    
    NonTracerMatrix <- matrix(0,ncol=(AtomNumber["C"]+1)*(AtomNumber["H"]+1),nrow=(AtomNumber["C"]+1)*(AtomNumber["H"]+1))
    for(m in 1:nrow(NonTracerMatrix)){
      Cindex<-CompleteLabel[m,1]
      Hindex<-CompleteLabel[m,2]
      Current.Combinations<-Isotope.Combinations1%>%filter(C13==Cindex & H2==Hindex)
      if (nrow(Current.Combinations)!=0){
        for (i in 1:nrow(Current.Combinations)){
          p <- dbinom(Current.Combinations[i,3], AtomNumber["N"], NitrogenNaturalAbundace[2])*
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
        Hindex<- Current.Combinations[i,2]
        p <- dbinom(Current.Combinations[i,3], AtomNumber["N"], NitrogenNaturalAbundace[2])*
          dmultinom(unlist(c(AtomNumber["O"]-sum(Current.Combinations[i,4:5]),Current.Combinations[i,4:5])), AtomNumber["O"], OxygenNaturalAbundace)*
          dmultinom(unlist(c(AtomNumber["S"]-sum(Current.Combinations[i,6:7]),Current.Combinations[i,6:7])), AtomNumber["S"], SulfurNaturalAbundace)*
          dmultinom(unlist(c(AtomNumber["Si"]-sum(Current.Combinations[i,8:9]),Current.Combinations[i,8:9])), AtomNumber["Si"], SiliconNaturalAbundace)*
          dbinom(Current.Combinations[i,10], AtomNumber["Cl"], ChlorineNaturalAbundance[2])*
          dbinom(Current.Combinations[i,11], AtomNumber["Br"], BromineNaturalAbundance[2])
        if (Cindex<Hindex){
          n<-which(CompleteLabel$C13==Cindex&CompleteLabel$H2==0)
          m<-which(CompleteLabel$C13==0&CompleteLabel$H2==Hindex)
          for (j in 0:(AtomNumber["C"]+1-n)){
            for (h in 0: (AtomNumber["H"]-m%/%(AtomNumber["C"]+1))){
              NonTracerMatrix[m+j+h*(AtomNumber["C"]+1),n+j+h*(AtomNumber["C"]+1)] <- NonTracerMatrix[m+j+h*(AtomNumber["C"]+1),n+j+h*(AtomNumber["C"]+1)]+p
            }
          }
        } else {
          n<-which(CompleteLabel$C13==0&CompleteLabel$H2==Hindex)
          m<-which(CompleteLabel$C13==Cindex&CompleteLabel$H2==0)
          for (j in 0:(AtomNumber["C"]+1-m)){
            for (h in 0: (AtomNumber["H"]-n%/%(AtomNumber["C"]+1))){
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
      CHMatrix<-CHMatrix[LabelIndex,LabelIndex]
      PurityMatrix<-PurityMatrix[LabelIndex,LabelIndex]
      ExpMatrix<-ExpMatrix[LabelIndex,]
      CorrectedMatrix<-CorrectedMatrix[LabelIndex,]
      CompleteLabel<-CompleteLabel[LabelIndex,]
    }
    
    for(i in 1:ncol(CorrectedMatrix)) {
      CorrectedMatrix[,i] <- coef(nnls(NonTracerMatrix %*% 
                                         CHMatrix %*% PurityMatrix, ExpMatrix[,i]))
    }
    CorrectedMatrix<-cbind(CompleteLabel,CorrectedMatrix)
    return(CorrectedMatrix)
  }
}


#' Natural Abundance correction for 13C-15N tracer labeling data
#'
#' @param formula string describing the chemical formula of the metabolite (in neutral form)
#' @param datamatrix matrix that contains the observed mass fractions
#' @param label matrix that contains the number of isotopes (C followed by N) associated with the observed mass fractions
#' @param Resolution For Exactive, the Resolution is 100000, defined at Mw 200
#' @param Autofill default: FALSE. If the intensity of an isotopologue is not provided, should 0 be assumed?
#' @param CarbonNaturalAbundace
#' @param HydrogenNaturalAbundace
#' @param NitrogenNaturalAbundace
#' @param SulfurNaturalAbundace
#' @param SiliconNaturalAbundace
#' @param ChlorineNaturalAbundance
#' @param BromineNaturalAbundance
#' @param C13Purity default:0.99.The isotopic purity for C13 tracer.
#' @param H2N15Purity default:0.99. The isotopic purity for H2/N15 tracer. 
#' @param ResDefAt Resolution defined at (in Mw), e.g. 200 Mw
#' @param ReportPoolSize default: TRUE
#' @import xlsx
#' @import gsubfn
#' @import nnls
#' @import dplyr
#' @import gdata
#' @import CHNOSZ
#' @import stringr
#' @return matrix of labeling pattern after the correction
CN_Correction <- function(formula, datamatrix, label, Resolution, 
                          Autofill = F, 
                          Charge = -1,
                          CarbonNaturalAbundace = c(0.9893, 0.0107),
                          HydrogenNaturalAbundace = c(0.999885, 0.000115),
                          NitrogenNaturalAbundace = c(0.99636, 0.00364),
                          OxygenNaturalAbundace = c(0.99757, 0.00038, 0.00205),
                          SulfurNaturalAbundace = c(0.95, 0.0075, 0.0425),
                          SiliconNaturalAbundace = c(0.92223, 0.04685, 0.03092),
                          ChlorineNaturalAbundance = c(0.7576,0.2424),
                          BromineNaturalAbundance = c(0.5069,0.4931),
                          C13Purity = 0.99,
                          H2N15Purity = 0.99,
                          ResDefAt = 200) {
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
      NPurityMatrix[i,] <- sapply(0:AtomNumber["N"], function(x) dbinom(x-i+1, x , 1-H2N15Purity))
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


#' Natural Abundance correction for 13C-15N or 13C-2H tracer labeling data
#'
#' @param InputFile String representing the name of the Input file
#' @param InputSheetName String representing the name of excel sheet that contains the input data
#' @param MetaboliteListName String representing the name of the excel database the name, formula and charge information of each metabolite
#' @param Isotopes String that specify the type of tracer isotopes, "CN" or "CH"
#' @param Resolution For Exactive, the Resolution is 100000, defined at Mw 200
#' @param Autofill default: FALSE. If the intensity of an isotopologue is not provided, should 0 be assumed?
#' @param CarbonNaturalAbundace
#' @param HydrogenNaturalAbundace
#' @param NitrogenNaturalAbundace
#' @param SulfurNaturalAbundace
#' @param SiliconNaturalAbundace
#' @param ChlorineNaturalAbundance
#' @param BromineNaturalAbundance
#' @param C13Purity default:0.99.The isotopic purity for C13 tracer.
#' @param H2N15Purity default:0.99. The isotopic purity for H2/N15 tracer. 
#' @param ResDefAt Resolution defined at (in Mw), e.g. 200 Mw
#' @param ReportPoolSize default: TRUE
#' @import xlsx
#' @import gsubfn
#' @import nnls
#' @import dplyr
#' @import gdata
#' @import CHNOSZ
#' @import stringr
#' @return New excel sheet named: 'Corrected', 'Normalized', 'Pool size' added to the original excel file.
#' @export

dual_correction <- function(InputFile,
                          InputSheetName,
                          MetaboliteListName,
                          Isotopes,
                          Resolution = 100000,
                          Autofill = F,
                          CarbonNaturalAbundace = c(0.9893, 0.0107),
                          HydrogenNaturalAbundace = c(0.999885, 0.000115),
                          NitrogenNaturalAbundace = c(0.99636, 0.00364),
                          OxygenNaturalAbundace = c(0.99757, 0.00038, 0.00205),
                          SulfurNaturalAbundace = c(0.95, 0.0075, 0.0425),
                          SiliconNaturalAbundace = c(0.92223, 0.04685, 0.03092),
                          ChlorineNaturalAbundance = c(0.7576,0.2424),
                          BromineNaturalAbundance = c(0.5069,0.4931),
                          C13Purity = 0.99,
                          H2N15Purity = 0.99,
                          ResDefAt = 200,
                          ReportPoolSize = TRUE){
  
  MetaboliteList <- read.csv(MetaboliteListName, header = TRUE, check.names = FALSE,stringsAsFactors=FALSE)
  InputDF <- read.xlsx(InputFile, InputSheetName, header = TRUE, check.names=FALSE, stringsAsFactors=FALSE)
  InputDF[,1] <- as.character(InputDF[,1])
  OutputMatrix <- matrix(0, ncol=(ncol(InputDF)-5), nrow=0)
  OutputPercentageMatrix <- matrix(0, ncol=(ncol(InputDF)-5), nrow=0)
  OutputPoolBefore <- matrix(0, ncol=(ncol(InputDF)-5), nrow=0)
  OutputPoolAfter <- matrix(0, ncol=(ncol(InputDF)-5), nrow=0)
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
  
  for (i in unique(InputDF$Compound)) {
    Formula=MetaboliteList[MetaboliteList$compound==i,2]
    Charge=MetaboliteList[MetaboliteList$compound==i,3]
    compound <- i
    if(length(Formula)==0) {
      print(paste("The formula of",i,"is unknown",sep=" "))
      break
    }
    CurrentMetabolite <- filter(InputDF, Compound==i)
    DataMatrix <- as.matrix(CurrentMetabolite[,6:(ncol(CurrentMetabolite))])
    if (Isotopes=="CH"){
      Corrected.raw <- CH_Correction(Formula, DataMatrix, CurrentMetabolite[,3:4],Resolution=Resolution,Charge = Charge, Autofill = Autofill, 
                                     CarbonNaturalAbundace = CarbonNaturalAbundace,
                                     HydrogenNaturalAbundace = HydrogenNaturalAbundace,
                                     NitrogenNaturalAbundace = NitrogenNaturalAbundace,
                                     OxygenNaturalAbundace = OxygenNaturalAbundace,
                                     SulfurNaturalAbundace = SulfurNaturalAbundace,
                                     SiliconNaturalAbundace = SiliconNaturalAbundace,
                                     ChlorineNaturalAbundance = ChlorineNaturalAbundance,
                                     BromineNaturalAbundance = BromineNaturalAbundance,
                                     C13Purity = C13Purity,
                                     H2N15Purity = H2N15Purity,
                                     ResDefAt = ResDefAt)
    } else if (Isotopes=="CN"){
      Corrected.raw <- CN_Correction(Formula, DataMatrix, CurrentMetabolite[,3:4],Resolution=Resolution,Charge = Charge, Autofill = Autofill)
    } else {
      print("The isotopes must be CH or CN, please check the input again")
      stop()
    }
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
  if (Isotopes=="CH"){
    names(OutputDF) <- c("Compound", "C13","H2", names(InputDF)[6:(length(names(InputDF)))])
  } else if(Isotopes=="CN"){
    names(OutputDF) <- c("Compound", "C13","N15", names(InputDF)[6:(length(names(InputDF)))])
  }
  names(OutputPercentageDF) <- names(OutputDF)
  names(OutputPoolBeforeDF) <- c("Compound", names(InputDF)[6:(length(names(InputDF)))])
  names(OutputPoolAfterDF) <- c("Compound", names(InputDF)[6:(length(names(InputDF)))])
  write.xlsx2(OutputDF, file=InputFile, sheetName = "Corrected", row.names=FALSE, append=TRUE)
  write.xlsx2(OutputPercentageDF, file=InputFile, sheetName = "Normalized", row.names=FALSE, append=TRUE)
  write.xlsx2(OutputPoolAfterDF, file=InputFile, sheetName = "Pool Size", row.names=FALSE, append=TRUE)
}



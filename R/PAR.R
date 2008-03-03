PAR <- function(D,E,C=NULL, model){
                
AtR               <- AdjAR(D,E,C,model)[,ncol(E)+1] 
bincomE           <- bincombinations(ncol(E))   
colnames(bincomE) <- colnames(E)

PAR <- c()

######### running through all exposure factors ##########

for(i in 1:(ncol(E))){  
    SAR  <- 0
    B    <- subset(bincomE, bincomE[,i]==1)
    AR_1 <- AtR[which(bincomE[,i]==1)]
    AR_2 <- AtR[which(bincomE[,i]==0)]
    SAR  <- AR_1 - AR_2    
    
    if(is.matrix(B)){
       sum_B <- rowSums(B)
    } else   
    {sum_B <- B}  
    
    ######### weights of the SAR ########
      
    WF <- (factorial(ncol(E)-(sum_B)))*(factorial(sum_B-1)) 
    WF <-  WF/(factorial(ncol(E)))  
    
    PAR[i] <- WF%*%SAR    
}

PAR           <- matrix(PAR, ncol=1)
rownames(PAR) <- colnames(E)
return(PAR)

}
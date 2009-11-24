AdjAR <- function(D,E,C=NULL,model){

all_cases <- length(D)

missings1 <- which(is.na(D))
missings2 <- which(is.na(E), arr.ind=TRUE)[,1]

if(is.null(C)==FALSE){
	missings3 <- which(is.na(C), arr.ind=TRUE)[,1]
}
if(is.null(C)==FALSE){
	missings <- unique(c(missings1, missings2, missings3))
} else{
				missings <- unique(c(missings1, missings2))
}

if(length(missings ==0)){

		D <- D[-missings]
		E <- E[-missings,]

		if(is.null(C)==FALSE){
			C <- C[-missings]
		}
		
	 warning(paste(length(missings), "out of", all_cases, "cases were removed because of missing values! "), call.=FALSE)
}

#source('Class_AR.r')

if(is.numeric(D)==TRUE){
		n <- length(D)
		name_D <- as.character(formula(model))[2] 
		D_m <- as.matrix(D)
		colnames(D_m) <- name_D
}
if(is.data.frame(D)){
		n   <- length(D[,1])
		D_m <- D
}

################## plausibility check of the entered data ######################
                                                                     

if(nrow(E)!=n) 
   stop(paste(sQuote("E"),"is not of length ",n))
if(!(is.null(C))){     
   if(nrow(C)!=n)
      stop(paste(sQuote("C"),"is not of length", n))
   if((sum(is.na(C)))>0) {
      stop(paste(sQuote("C"),"contains missing values"))}      
}

for(i in 1:(ncol(E))){
    if((length(levels(as.factor(E[,i]))))>2){
       stop(paste(sQuote("E"),"is not dichotomous"))}
}

if((length(levels(as.factor(D_m[,1]))))>2) {
   stop(paste(sQuote("D"),"is not dichotomous"))}

if((sum(is.na(E)))>0) {
   stop(paste(sQuote("E"),"contains missing values"))}
if((sum(is.na(D_m[,1])))>0) {
   stop(paste(sQuote("D"),"contains missing values"))}

######### call of the estimation of the (adjusted) attributable risk ###########

bincomE           <- bincombinations(ncol(E))                                      
colnames(bincomE) <- colnames(E)
AtR               <- c()                        
AtR[1]            <- 0   
AtR2              <- c()                        
AtR2[1]           <- 0                     
                     
for(i in 2:(nrow(bincomE))){                       
    conf        <- which(bincomE[i,]==0)                                   
    if(is.null(C)){
       help1  <- AR_woC(D_m,E,model,bincomE,conf)
       AtR[i] <- help1$Attributable_Risk_without_Confounding 
    }
    else{
       help2  <- AR_wC(D_m,E,C,model,bincomE,conf)  
       AtR[i] <- help2$Attributable_Risk_with_Confounding  
    }
}


if(is.matrix(C)){   
for(i in 1:(ncol(C))){
          a <- 0  
          if((length(levels(as.factor(C[,i]))))==2){
              a <- a + 1
          }
}
if(a > 0){   
               for(i in 1:(ncol(C))){
                   if(is.factor(C[,i]))
                      C[,i] <- as.numeric(levels(C[,i]))[C[,i]]
               }

              sample_EC          <- as.data.frame(cbind(C,E))
              bincomEC           <- bincombinations(ncol(sample_EC))
              colnames(bincomEC) <- colnames(sample_EC)
              for(i in 2:(nrow(bincomEC))){ 
                  conf    <- which(bincomEC[i,]==0)
                  help3   <- AR_woC(D_m,sample_EC,model,bincomEC,conf)
                  AtR2[i] <- help3$Attributable_Risk_without_Confounding
              }
          }
if(a==0){
				
}
}  
  

AtR <- cbind(bincomE, AtR)
rownames(AtR) <- rep("", nrow(AtR))
if(is.null(C)){
		combinations <- matrix()
		a <- list(AttRisk = AtR, OddsRatio = help1$Odds_Ratios, CondProb = help1$Conditional_Probability)
		b <- new("AR", AttRisk = a$AttRisk, OddsRatio=a$OddsRatio, CondProb=a$CondProb, combinations=combinations ) 
		return(b)
}
else{
	if(a>0){
			a <- list(AttRisk = AtR, Odds_Ratio_plot = help3$Odds_Ratios, CondProb_plot = help3$Conditional_Probability, combinations=bincomEC)
			b <- new("AR", AttRisk = a$AttRisk, OddsRatio=a$Odds_Ratio_plot, CondProb=a$CondProb_plot, combinations=a$combinations)
			return(b)
	} else{
	  OR <- matrix()
		CP <- numeric()
		comb <- matrix()  
		b <- new("AR", AttRisk = AtR, OddsRatio=OR , CondProb=CP, combinations=comb)
		return(b)
	}
}


}
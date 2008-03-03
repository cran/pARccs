AdjAR <- function(D,E,C=NULL,model){

n <- length(D) 

################## plausibility check of the entered data ######################
                                                                     
if(is.matrix(D))
   stop(paste(sQuote("D"),"is not a vector"))
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

if((length(levels(as.factor(D))))>2) {
   stop(paste(sQuote("D"),"is not dichotomous"))}

if((sum(is.na(E)))>0) {
   stop(paste(sQuote("E"),"contains missing values"))}
if((sum(is.na(D)))>0) {
   stop(paste(sQuote("D"),"contains missing values"))}

######### call of the estimation of the (adjusted) attributable risk ###########

bincomE           <- bincombinations(ncol(E))                                      
colnames(bincomE) <- colnames(E)
AtR               <- c()                        
AtR[1]            <- 0                       
                     
for(i in 2:(nrow(bincomE))){                         
    conf        <- which(bincomE[i,]==0)                                   
    if(is.null(C)){ 
       AtR[i] <- AR_woC(D,E,model,bincomE,conf) 
    }
    else {
       AtR[i] <- AR_wC(D,E,C,model,bincomE,conf)  
    }
}   

AtR <- cbind(bincomE, AtR)
return(AtR)
}
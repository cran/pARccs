AR_woC <- function(D,E,model,bincomE,conf=NULL){

########################### logistic regression #########################################

sample_DE  <- as.data.frame(cbind(D,E))
model.out  <- glm(model,family=binomial,data=sample_DE)       
coeff      <- coef(model.out)
formul     <- as.formula(paste("~", paste(as.character(formula(model)[3]))))

############################ model matrix ###############################################

M          <- model.matrix(formul,as.data.frame(bincomE))
M[,1]      <- 0      
M[,conf+1] <- 0   

if(length(conf)>0){
   splm <- strsplit(colnames(M),":")
   for(j in 1:length(colnames(M))){
       if(length(splm[[j]])>1){
          hb <- 0
          for(k in 1:length(splm[[j]])){
              for(i in 1:length(names(conf))){
                  if(names(conf)[i]==splm[[j]][k]) hb <- hb+1                                             
              }
          }
          if(hb==k) M[,j] <- 0
       }
   } 
}

############## attributable risk adjusted for the rest of exposure factors ############## 

OR        <- c()     
sum_terms <- c()    
pebd      <- c()   
p         <- c()

for(l in 1:nrow(bincomE)){
        a       <- t(matrix(bincomE[l,],ncol=nrow(E[D==1,]),nrow=ncol(bincomE)))
        q       <- abs(E[D==1,] - as.numeric(a))
        rq      <- rowSums(q)
        pebd[l] <- (sum(ifelse(rq==0,1,0)))/(nrow(sample_DE[D==1,]))
}

ln_OR     <- M%*%coeff
OR        <- exp(ln_OR)
sum_terms <- pebd/OR 
AtR_oC    <- 1-sum(sum_terms)  
return(AtR_oC)
}
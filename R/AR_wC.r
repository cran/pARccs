AR_wC <- function(D,E,C=NULL,model, bincomE, conf=NULL){

cn <- colnames(C)
en <- colnames(E)
gn <- c(en,cn)

pos_dum    <- c()
pos_dum[1] <- 1

if(ncol(C)==1){
conf_nam   <- colnames(C)
# C          <- as.vector(C)
C          <- as.factor(C[,1])
pos_dum[2] <- length(levels(C))+1
} else{
       C <- as.data.frame(C)
       for(i in 1:(ncol(C))){
           C[,i]        <- as.factor(C[,i])     
           pos_dum[i+1] <- length(levels(C[,i]))+pos_dum[i]                                  
       }
}
  
EC           <- data.frame(E,C)
colnames(EC) <- gn

########################### logistic regression ################################

sample_DEC          <- data.frame(D, EC)
model.out           <- glm(formula(model), family=binomial, data=sample_DEC)
coeff               <- coef(model.out)
coeff[is.na(coeff)] <- 0
new_formula         <- as.formula(paste("~",paste(names(coeff)[-1],collapse="+")))
             
####################### confounder as dummy variables ##########################

C_dummy   <- as.dummy(C, drop=TRUE, sep="")
C_dummy   <- ifelse(C_dummy==TRUE,1,0)
pos_dum_v <- pos_dum[-length(pos_dum)]
C_dummy   <- C_dummy[,-pos_dum_v] 

conf_names <- c()
if(ncol(as.matrix(C))==1){
   for(k in 1:ncol(as.matrix(C_dummy))){
       conf_names[k] <- names(coeff)[ncol(E)+1+k]
   }
   C_dummy <- as.matrix(C_dummy)
   colnames(C_dummy) <- conf_names  
}
 
if(is.vector(C_dummy)){
   C_dummy           <- matrix(C_dummy, ncol=1)
   colnames(C_dummy) <- c(names(coeff)[ncol(E)+2])
} 

EC_dummy <- cbind(E,C_dummy)                                    

bincomC           <- bincombinations(ncol(C_dummy))                                                    
colnames(bincomC) <- colnames(C_dummy)

for(j in 2:length(pos_dum)){
    pos_dum[j] <- pos_dum[j]-(j-1)
    a          <- (pos_dum[j-1]):(pos_dum[j]-1)
    if(length(a)>1) bincomC[rowSums(bincomC[,a])>1,] <- NA                          
}
bincomC <- na.omit(bincomC)

strat_AR <- c() 
 
############## running through the strata of the dummy confounder ##############


for(s in 1:nrow(bincomC)){        
		he               <- matrix(rep(bincomC[s,],nrow(bincomE)),nrow=nrow(bincomE), byrow=TRUE)
    colnames(he)     <- colnames(C_dummy)
    vm               <- cbind(bincomE, he) 
    dum_conf         <- c(conf,(ncol(E)+1):(ncol(vm)))                                      
    names(dum_conf)  <- colnames(vm)[dum_conf]
    
    ####################### model matrix ############################
                  
    M                <- model.matrix(new_formula, as.data.frame(vm))
    M[,1]            <- 0              
    M[,(dum_conf+1)] <- 0   

    if(length(dum_conf)>0){
       splm <- strsplit(colnames(M),":")
       for(j in 1:length(colnames(M))){
           if(length(splm[[j]])>1){
              hb <- 0
              for(k in 1:length(splm[[j]])){
                  for(i in 1:length(names(dum_conf))){
                      if(names(dum_conf)[i]==splm[[j]][k]) hb <- hb+1                                        
                  }
              }
              if(hb==k) M[,j] <- 0
           }
       } 
    }
    
    ######### attributable risk adjusted to the confounders #########
    ############# and to the rest of exposure factors ###############
              
    OR        <- c()    
    sum_terms <- c()   
    pebd      <- c()
		penod     <- c()    
    p         <- c()
    
    for(l in 1:nrow(vm)){
        if(sum(D[,1])==0) {pebd[l] <- 0} else {
           a  <- t(matrix(vm[l,],ncol=ifelse((is.vector(EC_dummy[D[,1]==1,])),1, nrow(EC_dummy[D[,1]==1,])), nrow=ncol(vm)))
           q  <- abs((EC_dummy[D[,1]==1,]) - as.numeric(a))
           a2 <- t(matrix(vm[l,],ncol=ifelse((is.vector(EC_dummy[D[,1]==0,])),1, nrow(EC_dummy[D[,1]==0,])), nrow=ncol(vm)))
           q2 <- abs(EC_dummy[D[,1]==0,] - as.numeric(a2))
           if(is.vector(q)){
                 rq  <- sum(q)
                 rq2 <- sum(q2)
           }else {
              rq   <- rowSums(q)
              rq2  <- rowSums(q2)
           }
           pebd[l]  <- (sum(ifelse(rq==0,1,0)))/(nrow(sample_DEC[sample_DEC[,1]==1,]))
           penod[l] <- (sum(ifelse(rq2==0,1,0)))/(nrow(sample_DEC[sample_DEC[,1]==0,]))
        }                        
    }
    
    if(s==1){
    		penod_null <- penod
    }
    
       ln_OR          <- M%*%coeff
       OR             <- exp(ln_OR)
       if(s == 1){
       		OR_null <- OR
       }
       sum_terms      <- pebd/OR     
       strat_AR[s]    <- sum(sum_terms) 
}

AtR_mc <- 1-sum(strat_AR) 
return(list(Attributable_Risk_with_Confounding = AtR_mc))
}



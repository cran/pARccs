Boot_CI <- function(D,E,C=NULL, model,  method=c("PAR", "AR"), stepwise=FALSE, scope=NULL, nboot=1000, alpha=0.025, original, type="perc", strat_boot=TRUE){

name_D <- as.character(formula(model))[2] 

D_m <- as.matrix(D)
colnames(D_m) <- name_D

################## plausibility check of the entered data ######################

if((alpha<0)||(alpha>1)) 
   stop(paste(sQuote("alpha"),"do not range between 0 and 1"))
if(nboot<=0)
   stop(paste(sQuote("nboot"),"have to be greater than 0"))

######################## nonparametric bootstrap ###############################

if(method=="PAR"){
   the_boot <- matrix(0,nrow=ncol(E), ncol=nboot) 
}
if(method=="AR"){
   the_boot <- matrix(0,nrow=nrow(bincombinations(ncol(E))), ncol=nboot)  
}   

for(i in 1:(nboot)){
    if(strat_boot==TRUE){
    ################### stratification ############################ 
       numbers_ca    <- which(D_m[,1]==1)
       num_sample_ca <- sample(numbers_ca, length(numbers_ca), replace=TRUE)
       numbers_co    <- which(D_m[,1]==0)
       num_sample_co <- sample(numbers_co, length(numbers_co), replace=TRUE)
       num_sample    <- c(num_sample_ca, num_sample_co)
    } 
    else{
    #################### without stratification ################### 
    numbers    <- c(1:(length(D_m[,1])))
    num_sample <- sample(numbers, length(D_m[,1]), replace=TRUE)      
    }
    E_sam           <- as.matrix(E[num_sample,])
    colnames(E_sam) <- colnames(E)
    D_sam           <- as.vector(D_m[num_sample,])
    if(stepwise==TRUE){
    ############ choose a model in a stepwise algorithm ###########
       if(!(is.null(C))){
          C_sam           <- as.matrix(C[num_sample,])
          colnames(C_sam) <- colnames(C)           
          for(j in 1:(ncol(C_sam))){
              if(!is.factor(C_sam[,j]))
              C_sam[,j] <- as.factor(C_sam[,j])                           
          }
       sample_DE    <- data.frame(D_sam, E_sam, C_sam)
       assign("sample_DE", sample_DE, env=.GlobalEnv)
       help_model   <- glm(model,family=binomial, data=sample_DE)
       sample_model <- step(help_model, scope, trace=FALSE)
        if(method=="PAR"){
          the_boot[,i] <- getPAR(PAR(D_sam,E_sam, C_sam, sample_model))
       }
       if(method=="AR"){
          this_AR <- AdjAR(D,E_sam, C_sam, sample_model)
          AR <- getAR(this_AR)
          the_boot[,i] <- AR[,ncol(E)+1]
          
       }
       } 
       else {
             data_DE      <- data.frame(D_sam,E_sam)
             assign("data_DE", data_DE, env=.GlobalEnv)
             help_model   <- glm(formula(model),family=binomial, data=data_DE)  
             sample_model <- step(help_model, scope, trace=FALSE)
						 if(method=="PAR"){
          			the_boot[,i] <- getPAR(PAR(D=D_sam,E=E_sam,model=sample_model))
       				}
       				if(method=="AR"){
      				  this_AR <- AdjAR(D=D,E=E_sam,model=sample_model)
                AR <- getAR(this_AR)
                the_boot[,i] <- AR[,ncol(E)+1]       
       				}  
       
       }
    }
    
    
    #################### use of the entered model ################# 
    else{ 
    if(!(is.null(C))){
       C_sam           <- as.matrix(C[num_sample,])
       colnames(C_sam) <- colnames(C)         
       for(j in 1:(ncol(C_sam))){
           if(!is.factor(C_sam[,j]))
              C_sam[,j] <- as.factor(C_sam[,j])                           
       }
    if(method=="PAR"){
          the_boot[,i] <- getPAR(PAR(D_sam,E_sam, C_sam, model))
       }
       if(method=="AR"){
          this_AR <- AdjAR(D,E_sam,C_sam,model)
          AR <- getAR(this_AR)
          the_boot[,i] <- AR[,ncol(E)+1]
          
       }
    } else {
    if(method=="PAR"){
          			the_boot[,i] <- getPAR(PAR(D=D_sam,E=E_sam,C=NULL, model=model))
       				}
       				if(method=="AR"){
       				   this_AR <- AdjAR(D=D,E=E_sam,model=model)
                 AR <- getAR(this_AR)
    	           the_boot[,i] <- AR[,ncol(E)+1]        
       				}  
    }
    }
}


assign("the_boot", the_boot, envir=.GlobalEnv)
 
for(j in 1:nrow(the_boot)){
    the_boot[j,] <- sort(the_boot[j,])
}

############# percentile method - (1-2*alpha) confidence interval ##############

if((type=="perc")||(type=="both")){
   lower_num <- alpha*nboot        
   upper_num <- (1-alpha)*nboot
 
if(!((lower_num-floor(lower_num))==0)){
   help_k    <- floor((nboot+1)*alpha)
   lower_num <- help_k 
   upper_num <- nboot+1-help_k
}
 
if(lower_num<1) {lower_num=1
   warning("nboot * alpha < 1")}

if(upper_num>nboot) {upper_num <- nboot}
  
lower   <- the_boot[,lower_num]
upper   <- the_boot[,upper_num] 

CI_perc           <- as.matrix(cbind(lower,upper))
colnames(CI_perc) <- c("lower", "upper")

if(method=="PAR"){
  rownames(CI_perc) <- colnames(E)
}
if(method=="AR"){
   h <- bincombinations(ncol(E))
   colnames(h) <- colnames(E)
   rownames(h) <- rep("", nrow(h))
   CI_perc <- cbind(h, CI_perc)
}
 
}
  
 
################ BCa method - (1-2*alpha)confidence interval ###################

   ############## estimation of z0 - bias correction ##############
 
if((type=="bca")||(type=="both")){
   s_z0 <- c()
   s_z0[1] <- 0
   if(method=="PAR"){
       	begin <- 1:nrow(the_boot)
       }
       if(method=="AR"){
       		begin <- 2:nrow(the_boot)
       }
   for(j in begin){
			 if(method=="PAR"){
       		prop    <- ifelse(the_boot[j,]< getPAR(original)[j,1],1,0) 
       }
       if(method=="AR"){
       		prop    <- ifelse(the_boot[j,]< getAR(original)[j,ncol(E)+1],1,0) 
       }
                                        
       s_z0[j] <- qnorm(mean(prop))
   }
   if(is.infinite(s_z0[2])) 
      stop(paste(sQuote("nboot"),"is to small"))
	
   ############## estimation of a - acceleration ##################

if(method=="PAR"){
   the_acc <- matrix(0,nrow=ncol(E), ncol=length(D_m[,1])) 
}

if(method=="AR"){
   the_acc <- matrix(0,nrow=nrow(bincombinations(ncol(E))), ncol=length(D_m[,1]))  
}  
 
E_orig    <- E
D_orig    <- D_m

for(k in 1:(length(D_orig[,1]))){
    E               <- as.data.frame(E_orig[-k,])
    colnames(E)     <- colnames(E_orig)
    D               <- as.data.frame(D_orig[-k,])
    colnames(D)     <- name_D
    if(!(is.null(C))){
       C_acc           <- C[-k,]
       colnames(C_acc) <- colnames(C)
       sample_data     <- data.frame(D,E,C_acc)
       sample_model    <- glm(model, family=binomial, data=sample_data)
      if(method=="PAR"){
          the_acc[,k]     <- getPAR(PAR(D=D,E=E, C=C_acc, model=sample_model))
       }
        if(method=="AR"){
					this_acc <- AdjAR(D=D,E=E, C=C_acc, model=sample_model)
          AR <- getAR(this_acc)     
          the_acc[,k] <- AR[,ncol(E)+1]
       }
    } 
    else{
       sample_DE    <- data.frame(D,E)
       sample_model <- glm(model, family=binomial, data=sample_DE)
      if(method=="PAR"){
          the_acc[,k]  <- getPAR(PAR(D=D,E=E,C=NULL, model=sample_model))
       }
       if(method=="AR"){
          this_acc  <- AdjAR(D=D,E=E,model=sample_model)
          AR <- getAR(this_acc) 
          the_acc[,k] <- AR[,ncol(E)+1]
       }
    }
}
                               
D <- D_orig
E <- E_orig

acc_mean <- rowMeans(the_acc)
help_ma  <- (-1)*(the_acc)+acc_mean
help_a_z <- rowSums((help_ma)^3)
help_a_n <- 6*((rowSums((help_ma)^2))^(1/2))
acc      <-  (help_a_z)/(help_a_n)
acc[1] <- 0
 
   #################### bca confidence interval ####################
 
z_alpha     <- qnorm(alpha)
z_min_alpha <- qnorm(1-alpha)
 
z0_mod_1    <- s_z0 + z_alpha
z0_mod_2    <- s_z0 + z_min_alpha
 
help_alpha1 <- c()
help_alpha2 <- c()
alpha1      <- c()
alpha2      <- c()

       if(method=="PAR"){
          n_values <- ncol(E) 
       }
       if(method=="AR"){
          n_values <- nrow(bincombinations(ncol(E)))
       }
 
for(l in 1:(n_values)){
    help_alpha1[l] <- s_z0[l]+((z0_mod_1[l])/(1-acc[l]*z0_mod_2[l]))
    help_alpha2[l] <- s_z0[l]+((z0_mod_2[l])/(1-acc[l]*z0_mod_2[l]))
    alpha1[l]      <- pnorm(help_alpha1[l])
    alpha2[l]      <- pnorm(help_alpha2[l])     
}


lower_num <- alpha1*nboot
upper_num <- alpha2*nboot

assign("lower_num", lower_num, env=.GlobalEnv)
assign("upper_num", upper_num, env=.GlobalEnv)
 
lower <- c()
upper <- c()

for(p in 1:length(lower_num)){
    
    if((floor(lower_num[p]))==0)
       lower[p]=NA
    else{
    
    if(!((lower_num[p]-floor(lower_num[p]))==0)){
       f0l <- the_boot[p,floor(lower_num[p])] 
       f1l <- the_boot[p,ceiling(lower_num[p])]
       x0l <- floor(lower_num[p])
       x1l <- ceiling(lower_num[p])
       lower[p] <- f0l + ((f1l-f0l)/(x1l-x0l))*(lower_num[p]-x0l)
    } 
    else{
    lower[p] <- the_boot[p,lower_num[p]]}
     
    }
    if(!((upper_num[p]-floor(upper_num[p]))==0)){
       f0u <- the_boot[p,floor(upper_num[p])] 
       f1u <- the_boot[p,ceiling(upper_num[p])]
       x0u <- floor(upper_num[p])
       x1u <- ceiling(upper_num[p])
       upper[p] <- f0u + ((f1u-f0u)/(x1u-x0u))*(upper_num[p]-x0u)
    }
    else{
    upper[p] <- the_boot[p,upper_num[p]]
    }
}

CI_bca           <- as.matrix(cbind(lower,upper))

if(method=="PAR"){
  rownames(CI_bca) <- colnames(E)
}
if(method=="AR"){
   h <- bincombinations(ncol(E))
   colnames(h) <- colnames(E)
   rownames(h) <- rep("", nrow(h))
   CI_bca <- cbind(h, CI_bca)
}
}                                                                 
  
if(type=="perc"){
   empty <- matrix()  
   CI <- new("Boot", CI_perc=CI_perc, CI_bca= empty)
} else{
     if(type=="bca"){
  		empty <- matrix()  
   		CI <- new("Boot", CI_perc=empty, CI_bca= CI_bca)
	 } else{ 
	 if(type=="both"){
				 CI <- new("Boot", CI_perc=CI_perc, CI_bca= CI_bca)
   }
  }
}

return(CI)

}






    




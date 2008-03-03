Boot_CI <- function(D,E,C=NULL, model, stepwise=FALSE, scope=NULL, nboot=250, alpha=0.025, original, type=c("perc","bca","both"), strat_boot=TRUE){

################## plausibility check of the entered data ######################

if((alpha<0)||(alpha>1)) 
   stop(paste(sQuote("alpha"),"do not range between 0 and 1"))
if(nboot<=0)
   stop(paste(sQuote("nboot"),"have to be greater than 0"))

######################## nonparametric bootstrap ###############################

PAR_boot <- matrix(0,nrow=ncol(E), ncol=nboot)   

for(i in 1:(nboot)){
    if(strat_boot==TRUE){
    ################### stratification ############################ 
       numbers_ca    <- which(D==1)
       num_sample_ca <- sample(numbers_ca, length(numbers_ca), replace=TRUE)
       numbers_co    <- which(D==0)
       num_sample_co <- sample(numbers_co, length(numbers_co), replace=TRUE)
       num_sample    <- c(num_sample_ca, num_sample_co)
    } 
    else{
    #################### without stratification ################### 
    numbers    <- c(1:(length(D)))
    num_sample <- sample(numbers, length(D), replace=TRUE)      
    }
    D_orig          <- D
    E_sam           <- as.matrix(E[num_sample,])
    colnames(E_sam) <- colnames(E)
    D               <- as.vector(D_orig[num_sample])
    if(stepwise==TRUE){
    ############ choose a model in a stepwise algorithm ###########
       if(!(is.null(C))){
          C_sam           <- as.matrix(C[num_sample,])
          colnames(C_sam) <- colnames(C)           
          for(j in 1:(ncol(C_sam))){
              if(!is.factor(C_sam[,j]))
              C_sam[,j] <- as.factor(C_sam[,j])                           
          }
       new_model    <- update(formula(model), D~.)
       sample_DE    <- data.frame(D, E_sam, C_sam)
       assign("sample_DE", sample_DE, env=.GlobalEnv)
       help_model   <- glm(new_model,family=binomial, data=sample_DE)
       sample_model <- step(help_model, scope, trace=FALSE)
       PAR_boot[,i] <- PAR(D,E_sam, C_sam, sample_model)[,1]
       } 
       else {
             new_model    <- update(formula(model), D~.)
             data_DE      <- data.frame(D,E_sam)
             assign("data_DE", data_DE, env=.GlobalEnv)
             help_model   <- glm(formula(new_model),family=binomial, data=data_DE)  
             sample_model <- step(help_model, scope, trace=FALSE)  
             PAR_boot[,i] <- PAR(D=D,E=E_sam,model=sample_model)[,1]
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
    PAR_boot[,i] <- PAR(D,E_sam,C_sam,model)[,1]} else {
    PAR_boot[,i] <- PAR(D=D,E=E_sam,model=model)[,1]}
    }
}

D <- D_orig

for(j in 1:nrow(PAR_boot)){
    PAR_boot[j,] <- sort(PAR_boot[j,])
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
  
lower   <- PAR_boot[,lower_num]
upper   <- PAR_boot[,upper_num]

CI_perc           <- as.matrix(cbind(lower,upper))
rownames(CI_perc) <- colnames(E)
colnames(CI_perc) <- c("lower", "upper")
 
}
  
 
################ BCa method - (1-2*alpha)confidence interval ###################

   ############## estimation of z0 - bias correction ##############
 
if((type=="bca")||(type=="both")){
   s_z0 <- c()
   for(j in 1:nrow(PAR_boot)){
       prop    <- ifelse(PAR_boot[j,]< original[j,1],1,0)
       s_z0[j] <- qnorm(mean(prop))
   }
   if(is.infinite(s_z0[1])) 
      stop(paste(sQuote("nboot"),"is to small"))
 
   ############## estimation of a - acceleration ##################
 
PAR_acc   <- matrix(0,nrow=ncol(E), ncol=length(D))
E_orig    <- E
D_orig    <- D


for(k in 1:(length(D_orig))){
    E               <- as.matrix(E_orig[-k,])
    colnames(E)     <- colnames(E_orig)
    D               <- as.vector(D_orig[-k])
    if(!(is.null(C))){
       C_acc           <- as.matrix(C[-k,])
       colnames(C_acc) <- colnames(C)
       sample_data     <- as.data.frame(cbind(D,E,C_acc))
       up_model  <- update(formula(model), D~.)
       assign("up_model", up_model, env=.GlobalEnv)
       sample_model    <- glm(up_model, family=binomial, data=sample_data)
       PAR_acc[,k]     <- PAR(D=D,E=E, C=C_acc, model=sample_model)[,1]
    } 
    else{
       sample_DE    <- as.data.frame(cbind(D,E))
       up_model  <- update(formula(model), D~.)
       assign("up_model", up_model, env=.GlobalEnv)
       sample_model <- glm(up_model, family=binomial, data=sample_DE)
       PAR_acc[,k]  <- PAR(D=D,E=E,model=sample_model)[,1]
    }
}
D <- D_orig
E <- E_orig


acc_mean <- rowMeans(PAR_acc)
help_ma  <- (-1)*(PAR_acc)+acc_mean
help_a_z <- rowSums((help_ma)^3)
help_a_n <- 6*((rowSums((help_ma)^2))^(1/2))
acc      <-  (help_a_z)/(help_a_n)
 
   #################### bca confidence interval ####################
 
z_alpha     <- qnorm(alpha)
z_min_alpha <- qnorm(1-alpha)
 
z0_mod_1    <- s_z0 + z_alpha
z0_mod_2    <- s_z0 + z_min_alpha
 
help_alpha1 <- c()
help_alpha2 <- c()
alpha1      <- c()
alpha2      <- c()
 
for(l in 1:(ncol(E))){
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
       f0l <- PAR_boot[p,floor(lower_num[p])] 
       f1l <- PAR_boot[p,ceiling(lower_num[p])]
       x0l <- floor(lower_num[p])
       x1l <- ceiling(lower_num[p])
       lower[p] <- f0l + ((f1l-f0l)/(x1l-x0l))*(lower_num[p]-x0l)
    } 
    else{
    lower[p] <- PAR_boot[p,lower_num[p]]}
     
    }
    if(!((upper_num[p]-floor(upper_num[p]))==0)){
       f0u <- PAR_boot[p,floor(upper_num[p])] 
       f1u <- PAR_boot[p,ceiling(upper_num[p])]
       x0u <- floor(upper_num[p])
       x1u <- ceiling(upper_num[p])
       upper[p] <- f0u + ((f1u-f0u)/(x1u-x0u))*(upper_num[p]-x0u)
    }
    else{
    upper[p] <- PAR_boot[p,upper_num[p]]
    }
}

CI_bca           <- as.matrix(cbind(lower,upper))
rownames(CI_bca) <- colnames(E)
}                                                                 
  
if(type=="perc"){  
   return(CI_perc)} else{if(type=="bca"){
   return(CI_bca)} else{ if(type=="both"){
   return(list(Perzentil=CI_perc, BCa=CI_bca))
   }
  }
}

}






    




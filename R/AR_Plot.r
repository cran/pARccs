
AR_Plot <- function(AR, Confounding = FALSE, type=c("onlyplot", "addplot"), rplot_path=getwd(), title_plot="Depicting excess risk of disease", legend_plot=TRUE){ 

assign("AR", AR, envir = .GlobalEnv)

if((ncol(getCombi(AR))>=3) | (ncol(getAR(AR))>=3)) legend_plot=FALSE

if(Confounding == TRUE){

  if(sum(is.na(getOR(AR)))>0){
      stop(paste("Confounders have more than two categories. \n Only Confounders with two categories can used for AR plots."))  
  }
  # Position labels for plot #
  
  sum_CondProb <- cumsum(getCProb(AR))
  pos_text <- sum_CondProb- (getCProb(AR)/2)
  assign("sum_CondProb", sum_CondProb, env=.GlobalEnv)
  
  which_val <- seq(1, nrow(getOR(AR)), by=2^(ncol(getAR(AR)-1)))
  seg_a <- getOR(AR)[which_val]
  risk_null <- c()
  for(k in 1:length(seg_a)){
      help1 <- rep(seg_a[k],2^(ncol(getAR(AR))-1)) 
      risk_null <- c(risk_null, help1)          
  }
  seg_b <- getCProb(AR)[which_val]               
  
  # Create dataset to plot #
  
  help1 <- getOR(AR)-risk_null
  help2a <- ifelse(help1>0, help1, 0)
  help2b <- ifelse(help1<0, help1, 0)
  
  neg <- 0
  for(i in 1:length(help1)){
      if(help1[i] < 0){
         neg <- neg + 1
      }                    
  }
   
  help3 <- as.matrix(cbind(risk_null, help2a, help2b))
  prob_zero <- which(getCProb(AR)==0)
  colnames(help3) <- c("1","OR-1", "Neg")
  plot_data <- t(help3)
  if(length(prob_zero)>0){
     pos_text <- pos_text[-prob_zero]
  }
  assign("pos_text", pos_text, env=.GlobalEnv)
  assign("plot_data", plot_data, env=.GlobalEnv)
  
  # External File for plot with labels to include, only up to two factors (one exposure and one confounder) #
  
  if((ncol(getCombi(AR)))<=2){
  
    nexp <- ncol(getAR(AR))-1
   
    Name_vis <- c()
    help2 <- getCombi(AR)
    Exp <- matrix(NA, nrow=nrow(help2), ncol=2*(ncol(help2))-1)
    
    for(j in 1:nrow(help2)){
        for(k in seq(1, 2*(ncol(help2))-1, by = 2)){
            if(k==1){
               if(help2[j,k]==0){
                if(k<(2*(ncol(help2))-1)){
                   Exp[j,k] <- paste("expression(paste(bar(C)[",k,"],")
                   Exp[j,k+1] <- paste("symbol('\\307')")
                }
                else{
                   Exp[j,k] <- paste("expression(paste(bar(C)[",k,"]))")
                }
               }
            else {
                 if(k<(2*(ncol(help2))-1)){
                   Exp[j,k] <- paste("expression(paste(C[",k,"],")
                     Exp[j,k+1] <- paste("symbol('\\307')")
                 }
                 else{
                   Exp[j,k] <- paste("expression(paste(C[",k,"]))") }
            }
          }
          else{
             if(help2[j,ceiling(k/2)]==0){
                if(k<(2*(ncol(help2))-1)){
                   if(k > (ncol(getCombi(AR)) - nexp)){
                      Exp[j,k] <- paste(", ' ', bar(E)[",(ceiling(k/2)-(ncol(getCombi(AR)) - nexp)),"],")
                   } 
                   else{
                        Exp[j,k] <- paste(", ' ', bar(C)[",ceiling(k/2),"],")
                   }
                   Exp[j,k-1] <- paste("symbol('\\307')")
                }
                else{
                  if(k > (ncol(getCombi(AR)) - nexp)){
                      Exp[j,k] <- paste(", ' ', bar(E)[",(ceiling(k/2)-(ncol(getCombi(AR)) - nexp)),"]))")
                   }
                   else{ 
                        Exp[j,k] <- paste(", ' ', bar(C)[",ceiling(k/2),"]))")
                   }
                   Exp[j,k-1] <- paste("symbol('\\307')")
                }
            }
            else {
                 if(k<(2*(ncol(help2))-1)){
                   if(k > (ncol(getCombi(AR)) - nexp)){
                       Exp[j,k] <- paste(", ' ', E[",(ceiling(k/2)-(ncol(getCombi(AR)) - nexp)),"],") 
                    }
                    else{
                         Exp[j,k] <- paste(", ' ', C[",ceiling(k/2),"],")
                    }
                    Exp[j,k-1] <- paste("symbol('\\307')")
                 }
                 else{
                    if(k > (ncol(getCombi(AR)) - nexp)){
                       Exp[j,k] <- paste(", ' ', E[",(ceiling(k/2)-(ncol(getCombi(AR)) - nexp)),"]))")
                    }
                    else{
                         Exp[j,k] <- paste(", ' ', C[",ceiling(k/2),"]))")
                    }
                    Exp[j,k-1] <- paste("symbol('\\307')")
                 }
            }
                
          }
        }
    }
    
    
    for(i in 1:nrow(Exp)){
        a <- ""
        for (j in 1:ncol(Exp)){
             a <- paste(a, Exp[i,j])
        }
        Name_vis[i] <- a
    }
    
    if(length(prob_zero)>0){
		   Name_vis <- Name_vis[-prob_zero]
    }
    
    a <- Name_vis[1]
    for(i in 2:length(Name_vis)){
        a <- paste(a, Name_vis[i], sep=",")
    }
    
    text_plot <- a
    
    exposures_names <- c()
    for(i in 1:(ncol(getCombi(AR)))){
        if(i > (ncol(getCombi(AR)) - nexp)){
           exposures_names[i] <- paste("expression(paste(E[",(i-(ncol(getCombi(AR)) - nexp)),"],' ... ",colnames(getCombi(AR))[i],"'))")
        }
        else{
        exposures_names[i] <- paste("expression(paste(C[",i,"],' ... ",colnames(getCombi(AR))[i],"'))")
        }
    }
       
    if(type=="onlyplot"){
      
      file_plot1 <- paste(rplot_path,"/plot1.r", sep="")    
      zz <- file(file_plot1, open = "w")
      cat("library(gplots)",file=zz, sep="\n")
      cat("x11()",file=zz, sep="\n")
      if(legend_plot==TRUE){
        cat("y_upper <- ceiling(max(getOR(AR)))", file=zz, sep="\n")
      }else{ cat("y_upper <- 0.5*ceiling(max(getOR(AR))/0.5)+0.2", file=zz, sep="\n")}
      cat("mp <- barplot2(plot_data, beside = FALSE, space = 0, angle =c(45, 135, 45)  , ylim = c(0,y_upper ) ,density = c(0, 30, 30), col=c('white',rep('grey',3)), width = c(getCProb(AR)), ylab='Odds Ratio', xlab='Exposure Distribution', xlim = c(0,1), col.axis='white')", file=zz, sep="\n")
      cat("axis(1,cex.axis=1.0,cex.lab=1.0)",file=zz, sep="\n")
      cat("axis(2,cex.axis=1.0,cex.lab=1.0)",file=zz, sep="\n")
      cat("segments(sum_CondProb,0, sum_CondProb, 0.5*ceiling(max(getOR(AR))/0.5),lty=3)", file=zz, sep="\n")
      cat("mtext(side = 1, at = c(pos_text) , cex = 0.9, line = -2, text =c(",file=zz, sep="")
      cat(text_plot,file=zz, sep="")
      cat("), col = 'black')",file=zz, sep="\n")
      cat("title(main='", file=zz, sep="")
      cat(title_plot, file=zz, sep="")
			cat("')", file=zz, sep="\n")
			if(legend_plot==TRUE){
	      cat("smartlegend(x='left', y.intersp = 1.5, y='top',c(", file=zz, sep="")
	      for(i in 1:(ncol(getCombi(AR)))){
	          cat(exposures_names[i],file=zz, sep = "")
	          if(i != (ncol(getCombi(AR)))){ 
	             cat(",", file=zz, sep = "")
	          }  
	      }
	      cat("))",file=zz, sep="\n")
      }else{
			  cat("x11(width=0.393700787+(max(nchar(colnames(getCombi(AR))))*0.1)/0.393700787, height=(ncol(getCombi(AR))*0.25)/0.393700787+ 0.393700787)", file=zz, sep="\n")
				cat("par(mar = c(0,0,0,0))", file=zz, sep="\n")
        cat("x=0.05", file=zz, sep="\n")
        cat("y=0.05", file=zz, sep="\n")
        cat("plot(x,y, axes= FALSE, col= 'white', ylab='', xlab='')", file=zz, sep="\n") 
        cat("smartlegend(x='center', y.intersp = 1.5, y='center',c(", file=zz, sep="")
	      for(i in 1:(ncol(getCombi(AR)))){
	          cat(exposures_names[i],file=zz, sep = "")
	          if(i != (ncol(getCombi(AR)))){ 
	             cat(",", file=zz, sep = "")
	          }
				} 
				cat("))",file=zz, sep="\n") 
			
			}
      close(zz)
      
      source(file_plot1)        
    }
    
    if(type=="addplot"){
    
      x11()
			l <- layout(matrix(c(1,1,2,2), 2,2, byrow = T), c(3,3), c(3,1,3))
      
      file_plot1 <- paste(rplot_path,"/plot1.r", sep="")
      zz <- file(file_plot1, open = "w")
      cat("library(gplots)",file=zz, sep="\n")
      cat("par(mar = c(4.5,4.5,5,3))", file=zz, sep="\n")
      if(legend_plot==TRUE){
        cat("y_upper <- ceiling(max(getOR(AR)))", file=zz, sep="\n")
      }else{ cat("y_upper <- 0.5*ceiling(max(getOR(AR))/0.5)+0.2", file=zz, sep="\n")}
      cat("mp <- barplot2(plot_data, beside = FALSE, space = 0, angle =c(45, 135, 45)  , ylim = c(0, y_upper) ,density = c(20, 30, 30), col=c('white', rep('grey',3)), width = c(getCProb(AR)), ylab='Odds Ratio', xlab='Exposure Distribution', xlim = c(0,1), col.axis='white')", file=zz, sep="\n")
      cat("axis(1,cex.axis=1.0,cex.lab=1.0)",file=zz, sep="\n")
      cat("axis(2,cex.axis=1.0,cex.lab=1.0)",file=zz, sep="\n")
      cat("segments(sum_CondProb,0, sum_CondProb, 0.5*ceiling(max(getOR(AR))/0.5) ,lty=3)", file=zz, sep="\n")
      cat("mtext(side = 1, at = c(pos_text) , cex = 0.9, line = -2, text =c(",file=zz, sep="\n")
      cat(text_plot,file=zz, sep="\n")
      cat("), col = 'black')",file=zz, sep="\n")
      cat("title(main='", file=zz, sep="")
      cat(title_plot, file=zz, sep="")
			cat("') ", file=zz, sep="\n")
			if(legend_plot==TRUE){
		      cat("smartlegend(x='left', y.intersp = 1.5, y='top',c(", file=zz, sep="\n")
		      for(i in 1:(ncol(getCombi(AR)))){
		          cat(exposures_names[i],file=zz, sep = "")
		          if(i != (ncol(getCombi(AR)))){ 
		             cat(",", file=zz, sep = "")
		          }  
		      }
		      cat("))",file=zz, sep="\n")
      }else{
        file_plot1a <- paste(rplot_path,"/plot1a.r", sep="")
        zza <- file(file_plot1a, open = "w")
        cat("x11(width=0.393700787+(max(nchar(colnames(getCombi(AR))))*0.1)/0.393700787, height=(ncol(getCombi(AR))*0.25)/0.393700787+ 0.393700787)", file=zza, sep="\n")
				cat("par(mar = c(0,0,0,0))", file=zza, sep="\n")
        cat("x=0.05", file=zza, sep="\n")
        cat("y=0.05", file=zza, sep="\n")
        cat("plot(x,y, axes= FALSE, col= 'white', ylab='', xlab='')", file=zza, sep="\n") 
        cat("smartlegend(x='center', y.intersp = 1.5, y='center',c(", file=zza, sep="")
	      for(i in 1:(ncol(getCombi(AR)))){
	          cat(exposures_names[i],file=zza, sep = "")
	          if(i != (ncol(getCombi(AR)))){ 
	             cat(",", file=zza, sep = "")
	          }
			   }
				 cat("))",file=zza, sep="\n")  
				 close(zza)
			
			}
      close(zz)
      
        if(neg ==0){
         text_a <- "The left-hatched area in proportion to the total area of all rectangles is equal to the \\n combined attributable risk of all exposure factors adjusted to the confounder(s) \\n = "
         text_b <- "\\n The graphic shows the situation that all exposure factors were eliminated totally \\n from the population." 
      }
      else{
         text_a <- "The left-hatched area minus the right-hatched area in proportion to the total area \\n of all rectangles (less the right-hatched area) is equal to the combined attributable \\n risk of all exposure factors adjusted to the confounder(s) ="
         text_b <- "\\n The graphic shows the situation that all exposure factors were eliminated \\n totally from the population." 
      }
      
      file_plot3 <- paste(rplot_path,"/plot3.r", sep="")    
      zz3 <- file(file_plot3, open="w")
      cat("par(mar = c(5,5,0,0))", file=zz3, sep="\n")
      cat("x=4", file=zz3, sep="\n")
      cat("y=4", file=zz3, sep="\n")
      cat("plot(x,y, axes= FALSE, col= 'white', ylab='', xlab='')", file=zz3, sep="\n")
      cat("mtext(side=1, cex=0.9, adj = 0,  line=0, text=paste('", file=zz3, sep="")
      cat(text_a, file=zz3, sep="")
      cat(round(getAR(AR)[nrow(getAR(AR)),ncol(getAR(AR))],digits=4), file=zz3, sep="")
      cat(".",file=zz3, sep="")
      cat(text_b, file=zz3, sep = "")
      cat("'))",file=zz3, sep="\n")
      close(zz3)
      
      source(file_plot1)
      source(file_plot3)
      source(file_plot1a)
    
    }
  }
  else
  {
  
  # Plots with lables outside of the bars (for more than 2 exposure factors)
  
  alphabet <- c("A", "B", "C", "D", "E", "F", "G", "H", "I", "J", "K", "L", "M","N", "O", "P", "Q", "R", "S", "T", "U", "V", "W", "X", "Y", "Z")
  combinations <- c(seq(1,nrow(getCombi(AR)), by=1))
  label <- alphabet[combinations]
  
  Name_vis <- c()
  help2 <- getCombi(AR)
  Exp <- matrix(NA, nrow=nrow(help2), ncol=2*(ncol(help2))-1)
  
  nexp <- ncol(getAR(AR))-1
   
  for(j in 1:nrow(help2)){
      for(k in seq(1, 2*(ncol(help2))-1, by = 2)){
          if(k==1){
          if(help2[j,k]==0){
              if(k<(2*(ncol(help2))-1)){
                 Exp[j,k] <- paste("expression(paste('",label[j],"... ', bar(C)[",k,"],")
                 Exp[j,k+1] <- paste("symbol('\\307')")
              }
              else{
                 Exp[j,k] <- paste("expression(paste('",label[j],"... ', bar(C)[",k,"]))")
              }
          }
          else {
               if(k<(2*(ncol(help2))-1)){
                 Exp[j,k] <- paste("expression(paste('",label[j]," ... ', C[",k,"],")
                   Exp[j,k+1] <- paste("symbol('\\307')")
               }
               else{
                 Exp[j,k] <- paste("expression(paste('",label[j]," ... ', C[",k,"]))") }
        }
        }
        else{
           if(help2[j,ceiling(k/2)]==0){
              if(k<(2*(ncol(help2))-1)){
                 if(ceiling(k/2) > (ncol(getCombi(AR)) - nexp)){
                    Exp[j,k] <- paste(", ' ', bar(E)[",(ceiling(k/2)-(ncol(getCombi(AR)) - nexp)),"],")
                 }
                 else{
                      Exp[j,k] <- paste(", ' ', bar(C)[",ceiling(k/2),"],")
                 }
                 
                 Exp[j,k-1] <- paste("symbol('\\307')")
              }
              else{
                 if(ceiling(k/2) > (ncol(getCombi(AR)) - nexp)){
                    Exp[j,k] <- paste(", ' ', bar(E)[",(ceiling(k/2)-(ncol(getCombi(AR)) - nexp)),"]))")
                 }
                 else{
                      Exp[j,k] <- paste(", ' ', bar(C)[",ceiling(k/2),"]))")
                 }
                 Exp[j,k-1] <- paste("symbol('\\307')")
              }
          }
          else {
               if(k<(2*(ncol(help2))-1)){
                   if(ceiling(k/2)> (ncol(getCombi(AR)) - nexp)){
                     Exp[j,k] <- paste(", ' ', E[",(ceiling(k/2)-(ncol(getCombi(AR)) - nexp)),"],")
                  }
                  else{
                      Exp[j,k] <- paste(", ' ', C[",ceiling(k/2),"],")
                  }
                  Exp[j,k-1] <- paste("symbol('\\307')")
               }
               else{
                  if(ceiling(k/2)> (ncol(getCombi(AR)) - nexp)){
                     Exp[j,k] <- paste(", ' ', E[",(ceiling(k/2)-(ncol(getCombi(AR)) - nexp)),"]))")
                  }
                  else{
                       Exp[j,k] <- paste(", ' ', C[",ceiling(k/2),"]))")
                  }
                  Exp[j,k-1] <- paste("symbol('\\307')")
               }
          }
              
       }
    }
  }
  

  
  for(i in 1:nrow(Exp)){
      a <- ""
      for (j in 1:ncol(Exp)){
           a <- paste(a, Exp[i,j])
      }
  Name_vis[i] <- a
  }
  
  if(length(prob_zero)>0){
	   label <- label[-prob_zero]
  }
  assign("label", label, env=.GlobalEnv)
  
  a <- Name_vis[1]
  for(i in 2:length(Name_vis)){
      a <- paste(a, Name_vis[i], sep=",")
  }
  text_plot <- a
  
  exposures_names <- c()
    for(i in 1:(ncol(getCombi(AR)))){
        if(i > (ncol(getCombi(AR)) - nexp)){
           exposures_names[i] <- paste("expression(paste(E[",(i-(ncol(getCombi(AR)) - nexp)),"],' ... ",colnames(getCombi(AR))[i],"'))")
        }
        else{
        exposures_names[i] <- paste("expression(paste(C[",i,"],' ... ",colnames(getCombi(AR))[i],"'))")
        }
    }
  
  if(type=="onlyplot"){
  
    file_plot1 <- paste(rplot_path,"/plot1.r", sep="") 
    zz <- file(file_plot1, open = "w")
    cat("library(gplots)",file=zz, sep="\n")
    cat("x11()",file=zz, sep="\n")
    cat("par(mar = c(4.5,4.5,4.5,4.5))", file=zz, sep="\n")
    if(legend_plot==TRUE){
        cat("y_upper <- ceiling(max(getOR(AR)))", file=zz, sep="\n")
    }else{ cat("y_upper <- 0.5*ceiling(max(getOR(AR))/0.5)+0.2", file=zz, sep="\n")}
    cat("mp <- barplot2(plot_data, beside = FALSE, space = 0, angle =c(45, 135, 45)  , ylim = c(0, y_upper) ,density = c(0, 30, 30), col=c('white',rep('grey',3)), width = c(getCProb(AR)), ylab='Odds Ratio', xlab='Exposure Distribution', xlim = c(0,1), col.axis='white')", file=zz, sep="\n")
    cat("axis(1,cex.axis=1.0,cex.lab=1.0)",file=zz, sep="\n")
    cat("axis(2,cex.axis=1.0,cex.lab=1.0)",file=zz, sep="\n")
    cat("segments(sum_CondProb,0, sum_CondProb, 0.5*ceiling(max(getOR(AR))/0.5) ,lty=3)", file=zz, sep="\n")
    cat("mtext(side = 1, at = c(pos_text) , cex = 1, line = -2, text =c(label), col='black')",file=zz, sep="\n")
    cat("title(main='", file=zz, sep="")
		cat(title_plot, file=zz, sep="")
		cat("') ", file=zz, sep="/n")
    close(zz)
    
    source(file_plot1)
    
    if(legend_plot==FALSE){
    file_plot2 <- paste(rplot_path,"/plot2.r", sep="")
	  zz2 <- file(file_plot2, open="w")
	  cat("x11(width=(ncol(getCombi(AR))/0.393700787)*0.5, height=((nrow(getCombi(AR))*0.2)/0.393700787))", file=zz2, sep="\n")
		cat("par(mar = c(2,1,4,1))", file=zz2, sep="\n")
	  cat("x=4", file=zz2, sep="\n")
	  cat("y=4", file=zz2, sep="\n")
	  cat("plot(x,y, axes= FALSE, col= 'white', ylab='', xlab='')", file=zz2, sep="\n")
	  cat("title(main='Legend')", file=zz2, sep="\n")
	  cat("smartlegend(x='center', y.intersp = 1.5, cex=1, y='top',c(", file=zz2, sep ="")
	  cat(text_plot, file=zz2, sep = "")
	  cat("))", file=zz2, sep="\n")
    cat("x11(width=0.393700787+(max(nchar(colnames(getCombi(AR))))*0.1)/0.393700787, height=(ncol(getCombi(AR))*0.25)/0.393700787+ 0.393700787)", file=zz2, sep="\n")
		cat("par(mar = c(0,0,0,0))", file=zz2, sep="\n")
    cat("x=0.05", file=zz2, sep="\n")
    cat("y=0.05", file=zz2, sep="\n")
    cat("plot(x,y, axes= FALSE, col= 'white', ylab='', xlab='')", file=zz2, sep="\n") 
    cat("smartlegend(x='center', y.intersp = 1.5, y='center',c(", file=zz2, sep="")
    for(i in 1:(ncol(getCombi(AR)))){
	        cat(exposures_names[i],file=zz2, sep = "")
	        if(i != (ncol(getCombi(AR)))){ 
	           cat(",", file=zz2, sep = "")
	        }  
    }
    cat("))",file=zz2, sep="\n")
    close(zz2)
    source(file_plot2)
		}
    
    
  
  }
   
  if(type=="addplot"){
  
		x11() 
	  if(legend_plot==TRUE){
	    l <- layout(matrix(c(rep(1, 2), seq(2,3),rep(4, 2)), 3,2, byrow = T), c(2,4,1,4,1), c(1,6,2))
	  }else{   l <- layout(matrix(c(rep(1, 2), rep(2,2),rep(3, 2)), 3,2, byrow = T), c(2,4,1), c(1,6,2))}
    
    file_plot1 <- paste(rplot_path,"/plot1.r", sep="")
    zz <- file(file_plot1, open = "w")
    cat("library(gplots)",file=zz, sep="\n")
    cat("par(mar = c(4.5,4,0,2))", file=zz, sep="\n")
    if(legend_plot==TRUE){
        cat("y_upper <- ceiling(max(getOR(AR)))", file=zz, sep="\n")
    }else{ cat("y_upper <- 0.5*ceiling(max(getOR(AR))/0.5)+0.2", file=zz, sep="\n")}
    cat("mp <- barplot2(plot_data, beside = FALSE, space = 0, angle =c(45, 135, 45)  , ylim = c(0, y_upper) ,density = c(0, 30, 30), col=c('white',rep('grey',3)), width = c(getCProb(AR)), ylab='Odds Ratio', xlab='Exposure Distribution', xlim = c(0,1), col.axis='white')", file=zz, sep="\n")
    cat("axis(1,cex.axis=1.0,cex.lab=1.0)",file=zz, sep="\n")
    cat("axis(2,cex.axis=1.0,cex.lab=1.0)",file=zz, sep="\n")
    cat("segments(sum_CondProb,0, sum_CondProb, 0.5*ceiling(max(getOR(AR))/0.5) ,lty=3)", file=zz, sep="\n")
    cat("mtext(side = 1, at = c(pos_text) , cex = 1, line = -2, text =c(label), col='black')",file=zz, sep="\n")
    # cat("title(main='Depicting excess risk of disease') ", file=zz, sep="/n")
    close(zz)
    
	    file_plot2 <- paste(rplot_path,"/plot2.r", sep="")
	    zz2 <- file(file_plot2, open="w")
	    cat("par(mar = c(2,1,4,1))", file=zz2, sep="\n")
	    cat("x=4", file=zz2, sep="\n")
	    cat("y=4", file=zz2, sep="\n")
	    cat("plot(x,y, axes= FALSE, col= 'white', ylab='', xlab='')", file=zz2, sep="\n")
	    cat("title(main='Legend')", file=zz2, sep="\n")
	    if(legend_plot==FALSE){
	     cat("smartlegend(x='center',cex=1, y.intersp = 1.5, y='top',c(", file=zz2, sep ="")}else{
	    cat("smartlegend(x='center',cex=1.2, y.intersp = 1.5, y='top',c(", file=zz2, sep ="")
	    }
	    cat(text_plot, file=zz2, sep = "")
	    cat("))", file=zz2, sep="\n")
	    if(legend_plot==FALSE){
			  cat("x11(width=0.393700787+(max(nchar(colnames(getCombi(AR))))*0.1)/0.393700787, height=(ncol(getCombi(AR))*0.25)/0.393700787+0.393700787)", file=zz2, sep="\n")
				cat("par(mar = c(0,0,0,0))", file=zz2, sep="\n")
        cat("x=0.05", file=zz2, sep="\n")
        cat("y=0.05", file=zz2, sep="\n")
        cat("plot(x,y, axes= FALSE, col= 'white', ylab='', xlab='')", file=zz2, sep="\n") 
        cat("smartlegend(x='center', y.intersp = 1.5, y='center',c(", file=zz2, sep="")
      }else{ cat("smartlegend(x='center', y.intersp = 1.5, cex=1.2, y='bottom',c(", file=zz2, sep="")}
	    for(i in 1:(ncol(getCombi(AR)))){
	        cat(exposures_names[i],file=zz2, sep = "")
	        if(i != (ncol(getCombi(AR)))){ 
	           cat(",", file=zz2, sep = "")
	        }  
	    }
	    cat("))",file=zz2, sep="\n")
	    close(zz2)
    
    
    
    if(neg ==0){
         text_a <- "The left-hatched area in proportion to the total area of all rectangles is equal to the \\n combined attributable risk of all exposure factors adjusted to the confounder(s) = "
         text_b <- "\\n The graphic shows the situation that all exposure factors were eliminated totally \\n from the population." 
      }
      else{
         text_a <- "The left-hatched area minus the right-hatched area in proportion to the total area of \\n all rectangles (less the right-hatched area) is equal to the combined attributable risk \\n of all exposure factors adjusted to the confounder(s) ="
         text_b <- "\\n The graphic shows the situation that all exposure factors were eliminated \\n totally from the population." 
      }
    
    file_plot5 <- paste(rplot_path,"/plot5.r", sep="")  
    zz5 <- file(file_plot5, open="w") 
    cat("par(mar = c(4,5,1,3))", file=zz5, sep="\n")
    cat("x=4", file=zz5, sep="\n")
    cat("y=4", file=zz5, sep="\n")
    cat("plot(x,y, axes= FALSE, col= 'white', ylab='', xlab='')", file=zz5, sep="\n")   
    cat("mtext(side=1, cex=0.9, adj = 0,  line=0, text=paste('", file=zz5, sep="")
    cat(text_a, file=zz5, sep="")
    cat(round(getAR(AR)[nrow(getAR(AR)),ncol(getAR(AR))],digits=4), file=zz5, sep="")
    cat(".",file=zz5, sep="")
    cat(text_b, file=zz5, sep = "")
    cat("'))",file=zz5, sep="\n")
    close(zz5)
    
    file_plot4 <- paste(rplot_path,"/plot4.r", sep="")
    zz4 <- file(file_plot4, open="w")
    cat("par(mar = c(0,3,5,3))", file=zz4, sep="\n")
    cat("x=0.01", file=zz4, sep="\n")
    cat("y=0.01", file=zz4, sep="\n")
    cat("plot(x,y, axes= FALSE, col= 'white', ylab='', xlab='')", file=zz4, sep="\n")
    cat("title(main='", file=zz4, sep="")
		cat(title_plot, file=zz4, sep="")
		cat("')", file=zz4, sep="")
    close(zz4)
    
    source(file_plot4)
    if(legend_plot==TRUE){
       source(file_plot2)
    }
    source(file_plot1)
    source(file_plot5)
    if(legend_plot==FALSE){
			 x11(width=(ncol(getCombi(AR))/0.393700787)*0.5, height=((nrow(getCombi(AR))*0.2)/0.393700787))
			 source(file_plot2)
    }
    
  
  }
 }
}

if(Confounding == FALSE){
  
  # Position labels for plot #
  
  sum_CondProb <- cumsum(getCProb(AR))
  pos_text <- sum_CondProb- (getCProb(AR)/2)
  assign("sum_CondProb", sum_CondProb, env=.GlobalEnv)
    
  
  # Create dataset to plot #
  
  help1 <- getOR(AR)-1 
  
  neg <- 0
  for(i in 1:length(help1)){
      if(help1[i] < 0){
         neg <- neg + 1
      }
  }
  
  help2a <- ifelse(help1>0, help1, 0)
  help2b <- ifelse(help1<0, help1, 0)
    
  help2c <- rep(1,length(getOR(AR)))
  help3 <- as.matrix(cbind(help2c, help2a, help2b))
  prob_zero <- which(getCProb(AR)==0)
  colnames(help3) <- c("1","OR-1", "Neg")
  plot_data <- t(help3)
  if(length(prob_zero)>0){
     pos_text <- pos_text[-prob_zero]
  }
  assign("plot_data", plot_data, env=.GlobalEnv)
  assign("pos_text", pos_text, env=.GlobalEnv)
																							 
  # External File for plot with labels to include, only up to two exposure factors #
  
  if((ncol(getAR(AR))-1)<=2){
  
    Name_vis <- c()
    help2 <- getAR(AR)
    Exp <- matrix(NA, nrow=nrow(help2), ncol=2*(ncol(help2)-1)-1)
    
    for(j in 1:nrow(help2)){
        for(k in seq(1, 2*(ncol(help2)-1)-1, by = 2)){
            if(k==1){
               if(help2[j,k]==0){
                if(k<(2*(ncol(help2)-1)-1)){
                   Exp[j,k] <- paste("expression(paste(bar(E)[",k,"],")
                   Exp[j,k+1] <- paste("symbol('\\307')")
                }
                else{
                   Exp[j,k] <- paste("expression(paste(bar(E)[",k,"]))")
                }
               }
            else {
                 if(k<(2*(ncol(help2)-1)-1)){
                   Exp[j,k] <- paste("expression(paste(E[",k,"],")
                     Exp[j,k+1] <- paste("symbol('\\307')")
                 }
                 else{
                   Exp[j,k] <- paste("expression(paste(E[",k,"]))") }
            }
          }
          else{
             if(help2[j,ceiling(k/2)]==0){
                if(k<(2*(ncol(help2)-1)-1)){
                   Exp[j,k] <- paste(", ' ', bar(E)[",ceiling(k/2),"],")
                   Exp[j,k-1] <- paste("symbol('\\307')")
                }
                else{
                   Exp[j,k] <- paste(", ' ', bar(E)[",ceiling(k/2),"]))")
                   Exp[j,k-1] <- paste("symbol('\\307')")
                }
            }
            else {
                 if(k<(2*(ncol(help2)-1)-1)){
                    Exp[j,k] <- paste(", ' ', E[",ceiling(k/2),"],")
                    Exp[j,k-1] <- paste("symbol('\\307')")
                 }
                 else{
                    Exp[j,k] <- paste(", ' ', E[",ceiling(k/2),"]))")
                    Exp[j,k-1] <- paste("symbol('\\307')")
                 }
            }
                
          }
        }
    }
    
    
    for(i in 1:nrow(Exp)){
        a <- ""
        for (j in 1:ncol(Exp)){
             a <- paste(a, Exp[i,j])
        }
        Name_vis[i] <- a
    }
    
    
    if(length(prob_zero)>0){
		   Name_vis <- Name_vis[-prob_zero]
    }
    
    a <- Name_vis[1]
    for(i in 2:length(Name_vis)){
        a <- paste(a, Name_vis[i], sep=",")
    }
    text_plot <- a
    
    exposures_names <- c()
    for(i in 1:(ncol(getAR(AR))-1)){
        exposures_names[i] <- paste("expression(paste(E[",i,"],' ... ",colnames(getAR(AR))[-ncol(getAR(AR))][i],"'))")
    }
    
    if(type=="onlyplot"){
      file_plot1 <- paste(rplot_path,"/plot1.r", sep="") 
      zz <- file(file_plot1, open = "w")
      cat("library(gplots)",file=zz, sep="\n")
      cat("x11()",file=zz, sep="\n")
      if(legend_plot==TRUE){
        cat("y_upper <- ceiling(max(getOR(AR)))", file=zz, sep="\n")
      }else{ cat("y_upper <- 0.5*ceiling(max(getOR(AR))/0.5)+0.2", file=zz, sep="\n")}
      cat("mp <- barplot2(plot_data, beside = FALSE, space = 0, angle =c(45, 135, 45)  , ylim = c(0, y_upper) ,density = c(0, 30, 30), col=c('white',rep('grey',3)), width = c(getCProb(AR)), ylab='Odds Ratio', xlab='Exposure Distribution', xlim = c(0,1), col.axis='white')", file=zz, sep="\n")
      cat("axis(1,cex.axis=1.0,cex.lab=1.0)",file=zz, sep="\n")
      cat("axis(2,cex.axis=1.0,cex.lab=1.0)",file=zz, sep="\n")
      cat("segments(sum_CondProb,0, sum_CondProb, 0.5*ceiling(max(getOR(AR))/0.5) ,lty=3)", file=zz, sep="\n")
      cat("mtext(side = 1, at = c(pos_text) , cex = 0.9, line = -2, text =c(",file=zz, sep="")
      cat(text_plot,file=zz, sep="")
      cat("), col = 'black')",file=zz, sep="\n")
      cat("title(main='", file=zz, sep="")
      cat(title_plot, file=zz, sep="")
			cat("') ", file=zz, sep="\n")
			if(legend_plot==TRUE){
		      cat("smartlegend(x='left', y.intersp = 1.5, y='top',c(", file=zz, sep="")
		      for(i in 1:(ncol(getAR(AR))-1)){
		          cat(exposures_names[i],file=zz, sep = "")
		          if(i != (ncol(getAR(AR))-1)){ 
		             cat(",", file=zz, sep = "")
		          }  
		      }
		      cat("))",file=zz, sep="\n")
      }else{
			  cat("x11(width=0.393700787+(max(nchar(colnames(getAR(AR))))*0.1)/0.393700787, height=(ncol(getAR(AR))*0.25)/0.393700787+ 0.393700787)", file=zz, sep="\n")
				cat("par(mar = c(0,0,0,0))", file=zz, sep="\n")
        cat("x=0.05", file=zz, sep="\n")
        cat("y=0.05", file=zz, sep="\n")
        cat("plot(x,y, axes= FALSE, col= 'white', ylab='', xlab='')", file=zz, sep="\n") 
        cat("smartlegend(x='center', y.intersp = 1.5, y='center',c(", file=zz, sep="")
	      for(i in 1:(ncol(getAR(AR))-1)){
	          cat(exposures_names[i],file=zz, sep = "")
	          if(i != (ncol(getAR(AR))-1)){ 
	             cat(",", file=zz, sep = "")
	          }  
				 }
				 cat("))",file=zz, sep="\n")
			}
      
      close(zz)
      
      source(file_plot1)
    }
    
    if(type=="addplot"){
    
			x11()
      l <- layout(matrix(c(1,1,2,2), 2,2, byrow = T), c(3,3), c(3,1,3))
      
      file_plot1 <- paste(rplot_path,"/plot1.r", sep="")
      zz <- file(file_plot1, open = "w")
      cat("library(gplots)",file=zz, sep="\n")
      cat("par(mar = c(4.5,4.5,5,3))", file=zz, sep="\n")
      if(legend_plot==TRUE){
        cat("y_upper <- ceiling(max(getOR(AR)))", file=zz, sep="\n")
      }else{ cat("y_upper <- 0.5*ceiling(max(getOR(AR))/0.5)+0.2", file=zz, sep="\n")}
      cat("mp <- barplot2(plot_data, beside = FALSE, space = 0, angle =c(45, 135, 45)  , ylim = c(0, y_upper) ,density = c(0, 30, 30), col=c('white', rep('grey',3)), width = c(getCProb(AR)), ylab='Odds Ratio', xlab='Exposure Distribution', xlim = c(0,1), col.axis='white')", file=zz, sep="\n")
      cat("axis(1,cex.axis=1.0,cex.lab=1.0)",file=zz, sep="\n")
      cat("axis(2,cex.axis=1.0,cex.lab=1.0)",file=zz, sep="\n")
  	  cat("segments(sum_CondProb,0, sum_CondProb, 0.5*ceiling(max(getOR(AR))/0.5) ,lty=3)", file=zz, sep="\n")
      cat("mtext(side = 1, at = c(pos_text) , cex = 0.9, line = -2, text =c(",file=zz, sep="")
      cat(text_plot,file=zz, sep="")
      cat("), col = 'black')",file=zz, sep="\n")
      cat("title(main='", file=zz, sep="")
      cat(title_plot, file=zz, sep="")
			cat("') ", file=zz, sep="\n")
			if(legend_plot==TRUE){
		      cat("smartlegend(x='left', y.intersp = 1.5, y='top',c(", file=zz, sep="\n")
		      for(i in 1:(ncol(getAR(AR))-1)){
		          cat(exposures_names[i],file=zz, sep = "")
		          if(i != (ncol(getAR(AR))-1)){ 
		             cat(",", file=zz, sep = "")
		          }  
		      }
		      cat("))",file=zz, sep="\n")
      }else{
        file_plot2 <- paste(rplot_path,"/plot2.r", sep="")
        zz2 <- file(file_plot2, open = "w") 
			  cat("x11(width=0.393700787+(max(nchar(colnames(getAR(AR))))*0.1)/0.393700787, height=(ncol(getAR(AR))*0.25)/0.393700787+ 0.393700787)", file=zz2, sep="\n")
				cat("par(mar = c(0,0,0,0))", file=zz2, sep="\n")
        cat("x=0.05", file=zz2, sep="\n")
        cat("y=0.05", file=zz2, sep="\n")
        cat("plot(x,y, axes= FALSE, col= 'white', ylab='', xlab='')", file=zz2, sep="\n") 
        cat("smartlegend(x='center', y.intersp = 1.5, y='center',c(", file=zz2, sep="")
	      for(i in 1:(ncol(getAR(AR))-1)){
	          cat(exposures_names[i],file=zz2, sep = "")
	          if(i != (ncol(getAR(AR))-1)){ 
	             cat(",", file=zz2, sep = "")
	          }  
				 }
				 cat("))",file=zz2, sep="\n")
				 close(zz2)
			}
      close(zz)
      
      if(neg ==0){
         text_a <- "The left-hatched area in proportion to the total area of all rectangles is equal to the \\n combined attributable risk = "
         text_b <- "\\n The graphic shows the situation that all exposure factors were eliminated totally \\n from the population." 
      }
      else{
         text_a <- "The left-hatched area minus the right-hatched area in proportion to the total area of \\n all rectangles (less the right-hatched area) is equal to the combined attributable \\n risk ="
         text_b <- "The graphic shows the situation that all exposure factors were eliminated \\n totally from the population." 
      }
      
      file_plot3 <- paste(rplot_path,"/plot3.r", sep="")
      zz3 <- file(file_plot3, open="w")
      cat("par(mar = c(5,5,0,0))", file=zz3, sep="\n")
      cat("x=4", file=zz3, sep="\n")
      cat("y=4", file=zz3, sep="\n")
      cat("plot(x,y, axes= FALSE, col= 'white', ylab='', xlab='')", file=zz3, sep="\n")
      cat("mtext(side=1, cex=0.9, adj = 0,  line=0, text=paste('", file=zz3, sep="")
      cat(text_a, file=zz3, sep="")
      cat(round(getAR(AR)[nrow(getAR(AR)),ncol(getAR(AR))],digits=4), file=zz3, sep="")
      cat(".",file=zz3, sep="")
      cat(text_b, file=zz3, sep = "")
      cat("'))",file=zz3, sep="\n")
      close(zz3)
      
      source(file_plot1)
      source(file_plot3)
      if(legend_plot==FALSE){
				 source(file_plot2)
      }
      
    
    }
  }
  else
  {
  
  # Plots with lables outside of the bars (for more than 2 exposure factors
  
  alphabet <- c("A", "B", "C", "D", "E", "F", "G", "H", "I", "J", "K", "L", "M","N", "O", "P", "Q", "R", "S", "T", "U", "V", "W", "X", "Y", "Z")
  combinations <- c(seq(1,nrow(getAR(AR)), by=1))
  label <- alphabet[combinations]
  
  Name_vis <- c()
  help2 <- getAR(AR)
  Exp <- matrix(NA, nrow=nrow(help2), ncol=2*(ncol(help2)-1)-1)
  
  for(j in 1:nrow(help2)){
      for(k in seq(1, 2*(ncol(help2)-1)-1, by = 2)){
          if(k==1){
          if(help2[j,k]==0){
              if(k<(2*(ncol(help2)-1)-1)){
                 Exp[j,k] <- paste("expression(paste('",label[j],"... ', bar(E)[",k,"],")
                 Exp[j,k+1] <- paste("symbol('\\307')")
              }
              else{
                 Exp[j,k] <- paste("expression(paste('",label[j],"... ', bar(E)[",k,"]))")
              }
          }
          else {
               if(k<(2*(ncol(help2)-1)-1)){
                 Exp[j,k] <- paste("expression(paste('",label[j]," ... ', E[",k,"],")
                   Exp[j,k+1] <- paste("symbol('\\307')")
               }
               else{
                 Exp[j,k] <- paste("expression(paste('",label[j]," ... ', E[",k,"]))") }
        }
        }
        else{
           if(help2[j,ceiling(k/2)]==0){
              if(k<(2*(ncol(help2)-1)-1)){
                 Exp[j,k] <- paste(", ' ', bar(E)[",ceiling(k/2),"],")
                 Exp[j,k-1] <- paste("symbol('\\307')")
              }
              else{
                 Exp[j,k] <- paste(", ' ', bar(E)[",ceiling(k/2),"]))")
                 Exp[j,k-1] <- paste("symbol('\\307')")
              }
          }
          else {
               if(k<(2*(ncol(help2)-1)-1)){
                  Exp[j,k] <- paste(", ' ', E[",ceiling(k/2),"],")
                  Exp[j,k-1] <- paste("symbol('\\307')")
               }
               else{
                  Exp[j,k] <- paste(", ' ', E[",ceiling(k/2),"]))")
                  Exp[j,k-1] <- paste("symbol('\\307')")
               }
          }
              
       }
    }
  }
  
  
  for(i in 1:nrow(Exp)){
      a <- ""
      for (j in 1:ncol(Exp)){
           a <- paste(a, Exp[i,j])
      }
  Name_vis[i] <- a
  }
  
  if(length(prob_zero)>0){
	   label <- label[-prob_zero]
  }
  assign("label", label, env=.GlobalEnv)
    
  a <- Name_vis[1]
  for(i in 2:length(Name_vis)){
      a <- paste(a, Name_vis[i], sep=",")
  }
  text_plot <- a
  
  exposures_names <- c()
  for(i in 1:(ncol(getAR(AR))-1)){
      exposures_names[i] <- paste("expression(paste(E[",i,"],' ... ",colnames(getAR(AR))[-ncol(getAR(AR))][i],"'))")
  }
  
  if(type=="onlyplot"){
    file_plot1 <- paste(rplot_path,"/plot1.r", sep="")
    zz <- file(file_plot1, open = "w")
    cat("library(gplots)",file=zz, sep="\n")
    cat("x11()",file=zz, sep="\n")
    cat("par(mar = c(4.5,4.5,4.5,4.5))", file=zz, sep="\n")
    if(legend_plot==TRUE){
        cat("y_upper <- ceiling(max(getOR(AR)))", file=zz, sep="\n")
    }else{ cat("y_upper <- 0.5*ceiling(max(getOR(AR))/0.5)+0.2", file=zz, sep="\n")}
    cat("mp <- barplot2(plot_data, beside = FALSE, space = 0, angle =c(45, 135, 45)  , ylim = c(0, y_upper) ,density = c(0, 30, 30), col=c('white',rep('grey',3)), width = c(getCProb(AR)), ylab='Odds Ratio', xlab='Exposure Distribution', xlim = c(0,1), col.axis='white')", file=zz, sep="\n")
    cat("axis(1,cex.axis=1.0,cex.lab=1.0)",file=zz, sep="\n")
    cat("axis(2,cex.axis=1.0,cex.lab=1.0)",file=zz, sep="\n")
    cat("segments(sum_CondProb,0, sum_CondProb, 0.5*ceiling(max(getOR(AR))/0.5) ,lty=3)", file=zz, sep="\n")
    cat("mtext(side = 1, at = c(pos_text) , cex = 1, line = -2, text =c(label), col='black')",file=zz, sep="\n")
    cat("title(main='", file=zz, sep="")
    cat(title_plot, file=zz, sep="")
		cat("') ", file=zz, sep="/n")
    close(zz)
    
    source(file_plot1)
    
   if(legend_plot==FALSE){
    file_plot2 <- paste(rplot_path,"/plot2.r", sep="")
	  zz2 <- file(file_plot2, open="w")
	  cat("x11(width=(ncol(getAR(AR))/0.393700787)*0.5, height=((nrow(getAR(AR))*0.2)/0.393700787))", file=zz2, sep="\n")
		cat("par(mar = c(2,1,4,1))", file=zz2, sep="\n")
	  cat("x=4", file=zz2, sep="\n")
	  cat("y=4", file=zz2, sep="\n")
	  cat("plot(x,y, axes= FALSE, col= 'white', ylab='', xlab='')", file=zz2, sep="\n")
	  cat("title(main='Legend')", file=zz2, sep="\n")
	  cat("smartlegend(x='center', y.intersp = 1.5, cex=1, y='top',c(", file=zz2, sep ="")
	  cat(text_plot, file=zz2, sep = "")
	  cat("))", file=zz2, sep="\n")
    cat("x11(width=0.393700787+(max(nchar(colnames(getAR(AR))))*0.1)/0.393700787, height=(ncol(getAR(AR))*0.25)/0.393700787+ 0.393700787)", file=zz2, sep="\n")
		cat("par(mar = c(0,0,0,0))", file=zz2, sep="\n")
    cat("x=0.05", file=zz2, sep="\n")
    cat("y=0.05", file=zz2, sep="\n")
    cat("plot(x,y, axes= FALSE, col= 'white', ylab='', xlab='')", file=zz2, sep="\n") 
    cat("smartlegend(x='center', y.intersp = 1.5, y='center',c(", file=zz2, sep="")
    for(i in 1:(ncol(getAR(AR))-1)){
	        cat(exposures_names[i],file=zz2, sep = "")
	        if(i != (ncol(getAR(AR))-1)){ 
	           cat(",", file=zz2, sep = "")
	        }  
    }
    cat("))",file=zz2, sep="\n")
    close(zz2)
     
    source(file_plot2)
    }
   
  
  }
   
  if(type=="addplot"){
    
		x11() 
    if(legend_plot==TRUE){ 
	    l <- layout(matrix(c(rep(1, 2), seq(2,3),rep(4, 2)), 3,2, byrow = T), c(2,4,1,4,1), c(1,6,2))
    }else{   l <- layout(matrix(c(rep(1, 2), rep(2,2),rep(3, 2)), 3,2, byrow = T), c(2,4,1), c(1,6,2))}
    
    file_plot1 <- paste(rplot_path,"/plot1.r", sep="")
    zz <- file(file_plot1, open = "w")
    cat("library(gplots)",file=zz, sep="\n")
    cat("par(mar = c(4.5,4,0,2))", file=zz, sep="\n")
    if(legend_plot==TRUE){
        cat("y_upper <- ceiling(max(getOR(AR)))", file=zz, sep="\n")
    }else{ cat("y_upper <- 0.5*ceiling(max(getOR(AR))/0.5)+0.2", file=zz, sep="\n")}
    cat("mp <- barplot2(plot_data, beside = FALSE, space = 0, angle =c(45, 135,45)  , ylim = c(0, y_upper) ,density = c(0, 30,30), col=c('white', rep('grey',3)), width = c(getCProb(AR)), ylab='Odds Ratio', xlab='Exposure Distribution', xlim = c(0,1), col.axis='white')", file=zz, sep="\n")
    cat("axis(1,cex.axis=1.0,cex.lab=1.0)",file=zz, sep="\n")
    cat("axis(2,cex.axis=1.0,cex.lab=1.0)",file=zz, sep="\n")
    cat("segments(sum_CondProb,0, sum_CondProb, 0.5*ceiling(max(getOR(AR))/0.5) ,lty=3)", file=zz, sep="\n")
    cat("mtext(side = 1, at = c(pos_text) , cex = 1, line = -2, text =c(label), col='black')",file=zz, sep="\n")
    #cat("title(main='Depicting excess risk of disease') ", file=zz, sep="/n")
    close(zz)
    
	    file_plot2 <- paste(rplot_path,"/plot2.r", sep="")
	    zz2 <- file(file_plot2, open="w")  
	    cat("par(mar = c(2,1,4,1))", file=zz2, sep="\n")
	    cat("x=4", file=zz2, sep="\n")
	    cat("y=4", file=zz2, sep="\n")
	    cat("plot(x,y, axes= FALSE, col= 'white', ylab='', xlab='')", file=zz2, sep="\n")
	    cat("title(main='Legend') ", file=zz2, sep="\n")
	     if(legend_plot==FALSE){
  	      cat("smartlegend(x='center', y.intersp = 1.5, cex=1, y='center',c(", file=zz2, sep ="")
	     }else{ cat("smartlegend(x='center', y.intersp = 1.5, cex=1.2, y='top',c(", file=zz2, sep ="")}
	    cat(text_plot, file=zz2, sep = "")
	    cat("))", file=zz2, sep="\n")  
	    cat("x11(width=0.393700787+(max(nchar(colnames(getAR(AR))))*0.1)/0.393700787, height=(ncol(getAR(AR))*0.25)/0.393700787+ 0.393700787)", file=zz2, sep="\n")
	    cat("par(mar = c(0,0,0,0))", file=zz2, sep="\n")
      cat("x=0.05", file=zz2, sep="\n")
      cat("y=0.05", file=zz2, sep="\n")
      cat("plot(x,y, axes= FALSE, col= 'white', ylab='', xlab='')", file=zz2, sep="\n") 
      if(legend_plot==FALSE){
         cat("smartlegend(x='center', y.intersp = 1.5, cex=1, y='center',c(", file=zz2, sep="")
      }else{
	    cat("smartlegend(x='center', y.intersp = 1.5, cex=1.2, y='bottom',c(", file=zz2, sep="")
	    }
	    for(i in 1:(ncol(getAR(AR))-1)){
	        cat(exposures_names[i],file=zz2, sep = "")
	        if(i != (ncol(getAR(AR))-1)){ 
	           cat(",", file=zz2, sep = "")
	        }  
	    }
	    cat("))",file=zz2, sep="\n")
	    close(zz2)
    
    
    if(neg ==0){
         text_a <- "The left-hatched area in proportion to the total area of all rectangles is equal to the \\n combined attributable risk = "
         text_b <- "\\n The graphic shows the situation that all exposure factors were eliminated totally \\n from the population." 
      }
      else{
         text_a <- "The left-hatched area minus the right-hatched area in proportion to the total area of \\n all rectangles (less the right-hatched area) is equal to the combined attributable \\n risk ="
         text_b <- "The graphic shows the situation that all exposure factors were eliminated \\n totally from the population." 
      }
    
    file_plot3 <- paste(rplot_path,"/plot3.r", sep="")      
    zz3 <- file(file_plot3, open="w")
    cat("par(mar = c(4,5,1,3))", file=zz3, sep="\n")
    cat("x=4", file=zz3, sep="\n")
    cat("y=4", file=zz3, sep="\n")
    cat("plot(x,y, axes= FALSE, col= 'white', ylab='', xlab='')", file=zz3, sep="\n")
    cat("mtext(side=1, cex=0.9, adj = 0,  line=0, text=paste('", file=zz3, sep="")
    cat(text_a, file=zz3, sep="")
    cat(round(getAR(AR)[nrow(getAR(AR)),ncol(getAR(AR))],digits=4), file=zz3, sep="")
    cat(".",file=zz3, sep="")
    cat(text_b, file=zz3, sep = "")
    cat("'))",file=zz3, sep="\n")
    close(zz3)
    
    file_plot4 <- paste(rplot_path,"/plot4.r", sep="")
    zz4 <- file(file_plot4, open="w")
    cat("par(mar = c(0,3,5,3))", file=zz4, sep="\n")
    cat("x=0.01", file=zz4, sep="\n")
    cat("y=0.01", file=zz4, sep="\n")
    cat("plot(x,y, axes= FALSE, col= 'white', ylab='', xlab='')", file=zz4, sep="\n")
    cat("title(main='", file=zz4, sep="")
		cat(title_plot, file=zz4, sep="")
		cat("')", file=zz4, sep="")
    close(zz4)

    source(file_plot4)
		if(legend_plot==TRUE){ 
       source(file_plot2)
	  }   
    source(file_plot1)
    source(file_plot3)
    if(legend_plot==FALSE){ 
       x11(width=(ncol(getAR(AR))/0.393700787)*0.5, height=((nrow(getAR(AR))*0.2)/0.393700787))
			 source(file_plot2)
	  } 
    }
  }
 }
}

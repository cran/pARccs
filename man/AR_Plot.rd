\name{AR_Plot}
\alias{AR_Plot}
\title{Visualizing of attributable risks from case-control data}
\description{
	 Visualize attributable risks in situations with  single and multiple exposure factors including adjusted attributable risks in form of modified scaled Venn diagrams.  
}
\usage{
AR_Plot(AR, Confounding = FALSE, type=c("onlyplot", "addplot"), 
        rplot_path=getwd(), title_plot="Depicting excess risk of disease", 
				legend_plot=TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{AR}{an object of class "AR" containing the estimated attributable risks and also the odds ratios and conditional propabilities needed for plotting}
  \item{Confounding}{Is "TRUE" if confounders are included, otherwise (and default) it is "FALSE"}
  \item{type}{A plotting area only with die diagram you will get with choosing "onlyplot". A short description of the plot you will get additionally to the diagram in choosing "addplot".}
  \item{rplot_path}{Some temporary files for plotting will be created in this path. At default, it is the working directory.}
  \item{title_plot}{Title of the plot. Default title is "Depicting excess risk of disease".}
  \item{legend_plot}{Is "TRUE" if the legend should be intergrated in the plot (default). Otherwise, the legend is printed separately.}
}

\details{
   
  All of the exposure factors and confounders have to be dichotomous.
	
	The diagram has the odds ratio (as estimator for the relative risk) along the vertical axis and the distribution of exposure among the   controls (as approximation  for the distribution of the exposure) at the horizontal axis. 
	For more then one exposure, the needed values are estimated in the resutling combinations of all exposure factors and confounders.
	The combinations will be declared in legends.
	
	An additional plot for legends (coding of the combinations) is shown  if the analysis was conducted for more than 2 factors.
	So therefor "legend\_plot" is always "FALSE".
	   
}

\note{

   This plot is only recommended if you have not more than 3 variables (confounders and exposure factors together).
   
}

 
\value{
  
}
\references{            
              Eide, G.; Heuch I. (2001)
               Attributable fractions: fundamental concepts and their visualization
               \emph{Statistical Methods in Medical Research} \bold{10}, 159-93
               
}

\author{ Christiane Raemsch }


\seealso{ \code{\link{AdjAR}}} 
\examples{

set.seed(27092009)
dicho        <- c(0,1)
two          <- c(0,1,2)
cc_state     <- sample(dicho, 100, replace=TRUE, prob=c(0.6, 0.4))
exposure1    <- sample(dicho, 100, replace=TRUE)
exposure2    <- sample(dicho, 100, replace=TRUE, prob=c(0.2,0.8))
exposure3    <- sample(dicho, 100, replace=TRUE, prob=c(0.4,0.6))
exposure4    <- sample(dicho, 100, replace=TRUE)
confounder1  <- sample(dicho, 100, replace=TRUE, prob=c(0.3, 0.7))
confounder2  <- sample(dicho, 100, replace=TRUE)

#### 1 Exposure ####

relation     <- as.formula(cc_state~exposure1)
data_exp     <- cbind(exposure1)
AR_exposures <- AdjAR(D=cc_state, E=data_exp, model=relation)
AR_exposures_plot <-AR_Plot(AR_exposures, type="onlyplot", Confounding=FALSE, 
title_plot="Example 1", legend_plot=FALSE)
AR_exposures_plot <-AR_Plot(AR_exposures, type="addplot", Confounding=FALSE, 
title_plot="Example 1", legend_plot=TRUE)

#### 1 Exposure and 1 Confounder #### 

relation     <- as.formula(cc_state~exposure1+confounder1)
data_exp     <- cbind(exposure1)
conf         <- cbind(confounder1)
AR_exposures <- AdjAR(D=cc_state,C=conf, E=data_exp, model=relation)
AR_exposures_plot <-AR_Plot(AR_exposures, type="onlyplot", Confounding=TRUE)
AR_exposures_plot <-AR_Plot(AR_exposures, type="addplot", Confounding=TRUE, 
legend_plot=FALSE)

#### 1 Exposure and 2 Confounders ####

relation     <- as.formula(cc_state~exposure1+confounder1+confounder2)
data_exp     <- cbind(exposure1)
conf         <- cbind(confounder1, confounder2)
AR_exposures <- AdjAR(D=cc_state,C=conf, E=data_exp, model=relation)
AR_exposures_plot <-AR_Plot(AR_exposures, type="onlyplot", Confounding=TRUE)
AR_exposures_plot <-AR_Plot(AR_exposures, type="addplot", Confounding=TRUE)

#### 2 Exposures ####  

relation     <- as.formula(cc_state~exposure1+exposure2)
data_exp     <- cbind(exposure1,exposure2)
AR_exposures <- AdjAR(D=cc_state, E=data_exp, model=relation)
AR_exposures_plot <-AR_Plot(AR_exposures, type="onlyplot", Confounding=FALSE, 
legend_plot=TRUE)
AR_exposures_plot <-AR_Plot(AR_exposures, type="addplot", Confounding=FALSE)

#### 2 Exposures and 1 Confounder ####    

relation     <- as.formula(cc_state~exposure1+exposure2+confounder1)
data_exp     <- cbind(exposure1,exposure2)
conf         <- cbind(confounder1)
AR_exposures <- AdjAR(D=cc_state, C=conf, E=data_exp, model=relation)
AR_exposures_plot <-AR_Plot(AR_exposures, type="onlyplot", Confounding=TRUE)
AR_exposures_plot <-AR_Plot(AR_exposures, type="addplot", Confounding=TRUE)

#### 2 Exposures and 2 Confounder ####   

relation     <- as.formula(cc_state~exposure1+exposure2+confounder1+confounder2)
data_exp     <- cbind(exposure1,exposure2)
conf         <- cbind(confounder1, confounder2)
AR_exposures <- AdjAR(D=cc_state, C=conf, E=data_exp, model=relation)
AR_exposures_plot <-AR_Plot(AR_exposures, type="onlyplot", Confounding=TRUE)
AR_exposures_plot <-AR_Plot(AR_exposures, type="addplot", Confounding=TRUE)

#### 3 Expsoures ####    

relation     <- as.formula(cc_state~exposure2+exposure1+exposure4)
data_exp     <- cbind(exposure2,exposure1, exposure4)
AR_exposures <- AdjAR(D=cc_state, E=data_exp, model=relation)
AR_exposures_plot <-AR_Plot(AR_exposures, type="onlyplot", Confounding=FALSE)
AR_exposures_plot <-AR_Plot(AR_exposures, type="addplot", Confounding=FALSE)

}


\keyword{ manip }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line

\name{Boot_CI}
\alias{Boot_CI}
\alias{getCIperc}
\alias{getCIbca}

\title{ Bootstrap confidence intervals for (partial) attributable risks (AR and PAR) from case-control data }
\description{
  With \code{Boot_CI} you can determine confidence intervals for (partial) atributable 
  risks from case-control data.
  Therefor the nonparametric bootstrap is used with whose bootstrap replications eighter 
  percentile confidence intervals or BCa confidence intervals are developed (or both, if you want to).
}
\usage{
Boot_CI(D, E, C = NULL, model, method=c("PAR", "AR"), stepwise = FALSE, 
        scope = NULL, nboot = 1000, alpha = 0.025, original, 
        type = "perc", strat_boot = TRUE)
}

\arguments{
  \item{D}{ a vector which holds the case-control state ("1" = case, "0"=control)}
  \item{E}{ a matrix of the exposure factor/s (all of them have to be dichotomous!)}
  \item{C}{ a matrix of the confounder/s (all of them have to be categorical!)}
  \item{model}{ a model formula or an object of class "\code{glm}" }
  \item{method}{ defines, if confidence intervals for attributable risks (AR) or partial attributable risks (PAR) should be estimated}
  \item{stepwise}{ a logical value indicating whether a stepwise-selected model should 
                   be used in the computation, default is \code{FALSE} }
  \item{scope}{ a description of the variables which should be taken into account in the 
                stepwise selection (upper model) and which variables are necessarily part 
                of the model (lower model)}
  \item{nboot}{ number of (bootstrap-)replication, default is \code{250} }
  \item{alpha}{ left- and right-hand error (default is \code{0.025}), so you will get 
                a \eqn{100\cdot(1-2\alpha)}{100*(1-2*alpha)}\% confidence interval}
  \item{original}{ a vector of the computed partial attributable risks from the original data}
  \item{type}{ a description of the type of confidence intervals which should be computed, 
               \code{"perc"} stands for the percentile confidence interval, \code{"bca"} for 
               the BCa confidence interval. You should
               choose \code{"both"} if you want to have calculated both types of confidence intervals. 
               \code{type="perc"} is the default. }
  \item{strat_boot}{ a logical value indicating whether a stratified or a non-stratified 
                     bootstrap should be executed, default is \code{TRUE} }
}
\details{

The computation of the (partial) attributable risks from the data set does not take place in this function.
You have to estimate them separately and pass the results (from classes "AR"/"PAR") through \code{original} to the function 
\code{Boot_CI}. 

To generate the bootstrap sample in every replication step one may use eighter the stratified or the non-
stratified method. If \code{strat_boot=TRUE} the sampling occurs separately from case-data and control-data, 
otherwise the sampling occurs from the complete data set.
 
If \code{stepwise=TRUE} the logistic regression model fitting the data from the bootstrap sample is 
choosen in a stepwise algorithm by the AIC. Therefor the argument \code{scope} is needed (look 
\code{?step} for more information). Note, that at least the main effects of the exposure factors 
(and confounders) have to be part of the lower model, so that the stepwise algorithm is only used to 
identify the most significant interactions. 
If \code{stepwise=FALSE} the formula of the argument \code{model} 
is used to build a model fitting the data of the bootstrap sample.

The bootstrap replications for the partial attributable risks are used to build confidence intervals
(as default 95\% confidence intervals are computed).
Therefor two methods are implemented: the percentile method (\code{type="perc"}) and the bias-corrected 
and accelerated (BCa) method (\code{type="bca"}). In conjunction with the choice between these two methods
you should take note of the great computational effort by using the BCa method.


}

\note{
 Also if there are only a single exposure factor/confounder you have to enter a matrix, so
   this will be a matrix with only one column.
   
 The names of the variables (outcome, exposure factors, confounders) in the argument \code{model} 
  have to be identically to the (column-)names of the entered data.  
  Furthermore all given exposure factors and confounders have to be part of the argument \code{model}.
  
  It is also important that the given variables in \code{D}, \code{E} and \code{C} are not defined
  as factors.
  
  Validity of the interval estimation can only be taken for granted for data with simple random sampling, stratified
   random sampling or frequency-matching of controls.
	 
	 Please note, if nboot is too small, BCa interval estimation may fail (for lower, upper or both limits). In this case, the limit is given with "NA".   
}


\value{
  \code{Boot_CI} returns a object of class "Boot". To extract the percentile confidence interval use function \code{getCIperc}, for the BCa confidence interval use \code{getCIbca}. 
}
\references{ Efron, B.; Tibshirani, R. (1986)
              Bootstrap methods for standard errors, confidence intervals, and other measure 
              of statistical accuracy
              \emph{Statistical Science} \bold{1}, 54-75
              
             Efron, B.; Tibshirani, R. (1993)
              \emph{An Introduction to the Bootstrap}
              Chapman \& Hall (Monographs on Statistics and Applied Probability 57)
                   
}
\author{ Christiane Raemsch }

\seealso{ \code{\link{PAR}} }
\examples{

###### Computation of BCa confidence intervals #######
###### for the PAR if there are no confounders #######

set.seed(2007)
dicho           <- c(0,1)
cc_state        <- sample(dicho, 100, replace=TRUE)
exposure1       <- sample(dicho, 100, replace=TRUE, prob=c(0.7, 0.3))
exposure2       <- sample(dicho, 100, replace=TRUE, prob=c(0.4, 0.6))
relation        <- as.formula(cc_state~exposure1+exposure2)
data_exp        <- cbind(exposure1, exposure2)
PAR_exposures   <- PAR(cc_state, data_exp, model=relation)
CI_95           <- Boot_CI(D=cc_state, E=data_exp, model=relation, method="PAR", 
                           nboot=70,original=PAR_exposures, type="bca")
CI_95

###### Computation of percentile confidence intervals #######
###### for the AR if there are confounders            #######

set.seed(2008)
dicho           <- c(0,1)
cc_state        <- sample(dicho, 100, replace=TRUE)
exposure1       <- sample(dicho, 100, replace=TRUE, prob=c(0.7, 0.3))
exposure2       <- sample(dicho, 100, replace=TRUE, prob=c(0.4, 0.6))
cat_confounder  <- c(0,1,2,3)
confounder1     <- sample(cat_confounder, 100, replace=TRUE)
relation        <- as.formula(cc_state~exposure1+exposure2+confounder1)
data_exp        <- cbind(exposure1, exposure2)
conf            <- matrix(confounder1, ncol=1)
colnames(conf)  <- c("confounder1")
AR_exposures    <- AdjAR(cc_state, data_exp, conf, model=relation)
CI_95           <- Boot_CI(D=cc_state, E=data_exp, C=conf, model=relation, 
                           method="AR", nboot=70,original=AR_exposures)
CI_95

}
\keyword{ manip }




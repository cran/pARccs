\name{AR}
\alias{AR_woC}
\alias{AR_wC}
\title{Estimating (adjusted) attributable risks from case-control data}
\description{
  With the functions \code{AR_wC} and \code{AR_woC} the estimation of  the 
  attributable risks (AR) from case-control data is realized.
  
  From \code{AR_woC} you get the ARs for the exposure factors of primary interest
  adjusted to the rest of the exposure factors,
  the resulting ARs for the exposure factors of primary interest
  from \code{AR_wC} are additionally adjusted to the given confounders.
}
\usage{
AR_woC(D, E, model, bincomE, conf = NULL)

AR_wC (D, E, C = NULL, model, bincomE, conf = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{D}{a matrix which holds the case-control state ("1" = case, "0"=control) }
  \item{E}{a matrix of the exposure factor/s (all of them have to be dichotomous!) }
  \item{C}{a matrix of the confounder/s (all of them have to be categorical!)}
  \item{model}{ a model formula or an object of class "\code{glm}"}
  \item{bincomE }{a matrix which contains all binary combinations of the exposures in \code{E} }
  \item{conf}{ a vector which holds the corresponding number of column of the exposure factors 
  which should act as an additional confounder; the default is \code{NULL}, which means no exposure factor
  acts an as additional confounder
  }
}
\details{
   With \code{AR_wC} the (joint) attributable risk for the exposure factor(s) of primary interest, 
   which are not mentioned in vector \code{conf}, is estimated. It is adjusted to the rest of exposure 
   factor/s in \code{E} (these are defined by \code{conf}).
   If you want an additional adjustment to secondary confounders, use function \code{AR_woC} with
   indicating the confounder/s in \code{C}.
   
   If \code{conf=NULL} the joint attributable risk of all given exposure factors is estimated 
   (as the case may be with adjustment to the confounder/s).
   
   For the adjusted estimation regression models are used, here it is a logistic regression model. Through 
   this model the needed Odds Ratio (OR) is estimated.
   The argument \code{model} can be eighter of the form \code{D~terms}, where \code{terms} is a series of 
   terms out of the exposure factors and confounders, or an object of class "\code{glm}". 
   (In the process of model fitting with \code{glm} you have to choose \code{family=binomial} to get a 
   logistic regression model.)\cr
   All given exposure factors and confounders have to be part of the argument \code{model}.
   The names of the variables (outcome, exposure factors, confounders) in the argument \code{model} have to 
   be identical to the (col-)names of the entered data.
   Also the colnames of \code{bincomE} have to be identical to the colnames of \code{E}.
   
   To get the matrix \code{bincomE} you may use the function \code{bincombinations()}\cr
   (use \code{help(bincombinations, package=e1071)} for further information).
   
}

\note{ 
   Also if there are only a single exposure factor/confounder you have to enter a matrix, so
   this will be a matrix with only one column.
   
   It is also important that the given variables in \code{D}, \code{E} and \code{C} are not defined
   as factors.
   
   Validity of the estimation can only be taken for granted for data with simple random sampling, stratified
   random sampling or frequency-matching of controls.    
  
   Here the (adjusted) attributable risk for only one defined (binary) combination of the exposure factors 
   is estimated. To get the (adjusted) attributable risks for every possible (binary) combination of all 
   given exposure factors use function \code{\link{AdjAR}}.
 }
 
\value{
   \code{AR_woC} returns a list which contains the (joint) attributable risk of one (or more) exposure 
   factor(s) adjusted to the rest of the exposure factors, the odds ratio and the conditional probability for every combination of the exposures in \code{E} (last ones are needed if plotting the attributable risks) .

   \code{AR_wC} returns a single value which is the (joint) attributable risk of one (or more) exposure 
   factor(s) adjusted to the rest of the exposure factors and to the given confounders.
}
\references{
              Levin, M. (1953)
               The occurrence of lung cancer in man
               \emph{Acta Unio Internationalis Contra Cancrum} \bold{9}, 531-41

              Bruzzi, P.; Green S.; Byard, D. \emph{et al.} (1985)
               Estimating the population attributable risk for multiple risk factors using case-control data
               \emph{American Journal of Epidemiology} \bold{122}, 904-14

              Benichou, J. (1991)
               Methods of adjustment for estimating the attributable risk in case-control studies: a review
               \emph{Statistics in Medicine} \bold{10}, 1753-73
               }

\author{ Christiane Raemsch }


\seealso{ \code{\link{AdjAR}}} 
\examples{

##### use of function 'AR_woC':       #####
##### attributable risk for exposure2 #####
##### adjusted for exposure 1         #####

set.seed(2007)
dicho            <- c(0,1)
cc_state         <- as.matrix(sample(dicho, 100, replace=TRUE))
exposure1        <- sample(dicho, 100, replace=TRUE, prob=c(0.7, 0.3))
exposure2        <- sample(dicho, 100, replace=TRUE, prob=c(0.4, 0.6))
relation         <- as.formula(cc_state~exposure1+exposure2)
data_exp         <- cbind(exposure1, exposure2)
bincom           <- bincombinations(2)
colnames(bincom) <- colnames(data_exp)
AR_exposure2     <- AR_woC(cc_state, data_exp, relation, bincom, c(2))
AR_exposure2


##### use of function 'AR_wC':               #####
##### joint attributable risk for exposure1  #####
##### and exposure2 adjusted for confounder1 #####

set.seed(2008)
dicho            <- c(0,1)
cc_state         <- as.matrix(sample(dicho, 100, replace=TRUE))
exposure1        <- sample(dicho, 100, replace=TRUE, prob=c(0.7, 0.3))
exposure2        <- sample(dicho, 100, replace=TRUE, prob=c(0.4, 0.6))
cat_confounder   <- c(0,1,2,3)
confounder1      <- sample(cat_confounder, 100, replace=TRUE)
rel_mod          <- glm(cc_state~exposure1+exposure2+confounder1, 
                        family=binomial)
data_exp         <- cbind(exposure1, exposure2)
conf             <- matrix(confounder1, ncol=1)
colnames(conf)   <- c("confounder1")
bincom           <- bincombinations(2)
colnames(bincom) <- colnames(data_exp)
AR_exposure1_2   <- AR_wC(cc_state, data_exp, conf, rel_mod, bincom)
AR_exposure1_2 

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ manip }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line

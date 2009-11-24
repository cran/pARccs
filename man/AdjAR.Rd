\name{AdjAR}
\alias{AdjAR}
\alias{getAR}
\alias{getOR}
\alias{getCombi}
\alias{getCProb}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Calling the estimation of the (adjusted) attributable risks from case-control data}
\description{
   \code{AdjAR} realizes the estimation of (adjusted) attributable risk (AR) 
   from case-control data via logistic regression by calling the adequate 
   function which holds the computation}
\usage{
AdjAR(D, E, C = NULL, model)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{D}{ a vector which holds the case-control state ("1" = case, "0"=control)}
  \item{E}{ a matrix of the exposure factor/s (all of them have to be dichotomous!) }
  \item{C}{ a matrix of the confounder/s (all of them have to be categorical!) }
  \item{model }{ a model formula or an object of class "\code{glm}"}
}
\details{
   Depending from the entered data \code{AdjAR} accesses to two additional functions:
   \code{AR_woC} is selected if there is no variable which is only a confounder, 
   expressed as \code{C=NULL}.
   \code{AR_wC} is selected if there are also variables which only act as confounders, 
   that means \code{C} is a matrix. 
   
   See \code{\link{AR_woC}} and \code{\link{AR_wC}} for further information about the computation.
}

\note{
   
	 \code{AdjAR} handles with missing data. Cases with missing values are removed and a warning message is printed with the information how many cases were removed.
	 
	 Also if there are only a single exposure factor/confounder you have to enter a matrix, so
   this will be a matrix with only one column.
   
   It is also important that the given variables in \code{D}, \code{E} and \code{C} are not defined
   as factors.
   
   The names of the variables (outcome, exposure factor/s, confounder/s) in the argument \code{model} 
   have to be identical to the (column-)names of the entered data. 
   Furthermore all given exposure factors and confounders have to be part of the argument \code{model}.
   
   Validity of the estimation can only be taken for granted for data with simple random sampling, stratified
   random sampling or frequency-matching of controls.    

}

\value{
   \code{AdjAR} returns a object of class "AR" containing the attributable risk for every possible (binary) 
   combination of the exposure factors in \code{E}.
   
   If  \code{C=NULL} these are only adjusted to the rest of the exposure factors 
   (which are not part of the interested combination). 
   If there are given confounders in \code{C}
   the attributable risks are additionally adjusted to them.
}

\references{  Levin, M. (1953)
               The occurrence of lung cancer in man 
               \emph{Acta Unio Internationalis Contra Cancrum} \bold{9}, 531-41

              Bruzzi, P.; Green S.; Byard, D. \emph{et al.} (1985)
               Estimating the population attributable risk for multiple risk factors 
               using case-control data
               \emph{American Journal of Epidemiology} \bold{122}, 904-14 
               
              Benichou, J. (1991)
               Methods of adjustment for estimating the attributable risk in 
               case-control studies: a review
               \emph{Statistics in Medicine} \bold{10}, 1753-73      }
                
\author{Christiane Raemsch}

\seealso{ \code{\link{AR_woC}}, \code{\link{AR_wC}}}
\examples{


##### Computation of the AR for every combination of two #####
##### exposure factors if there are no confounders       #####

set.seed(2007)
dicho        <- c(0,1)
cc_state     <- sample(dicho, 100, replace=TRUE)
exposure1    <- sample(dicho, 100, replace=TRUE, prob=c(0.7, 0.3))
exposure2    <- sample(dicho, 100, replace=TRUE, prob=c(0.4, 0.6))
relation     <- as.formula(cc_state~exposure1+exposure2)
data_exp     <- cbind(exposure1, exposure2)
AR_exposures <- AdjAR(D=cc_state, E=data_exp, model=relation)
AR_exposures


##### Computation of the AR for every combination of two #####
##### exposure factors with adjustment to confounder1    #####

set.seed(2008)
cc_state        <- sample(dicho, 100, replace=TRUE)
exposure1       <- sample(dicho, 100, replace=TRUE, prob=c(0.9, 0.1))
exposure2       <- sample(dicho, 100, replace=TRUE, prob=c(0.6, 0.4))
cat_confounder  <- c(0,1,2,3)
confounder1     <- sample(cat_confounder, 100, replace=TRUE)
data_exp        <- cbind(exposure1, exposure2)
conf            <- matrix(confounder1, ncol=1)
colnames(conf)  <- c("confounder1")
rel_mod         <- glm(cc_state~exposure1+exposure2+confounder1,
                       family=binomial)
AR_exposures    <- AdjAR(cc_state, data_exp, conf, rel_mod)
AR_exposures



}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ manip } % data manipulation
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line










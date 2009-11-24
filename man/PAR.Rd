\name{PAR}
\alias{PAR}
\alias{getPAR}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ Estimating partial attributable risks from case-control data }

\description{ \code{PAR} estimates the partial attributable risks (PAR) for multiple exposure factors. 
The underlying data have to arise from a case-control-study. 
}     

\usage{
PAR(D, E, C = NULL, model)
}
%- maybe also 'usage' for other objects documented here.

\arguments{
  \item{D}{ a vector which holds the case-control state ("1" = case, "0"=control)}
  \item{E}{ a matrix of the exposure factor/s (all of them have to be dichotomous!)}
  \item{C}{ a matrix of the confounder/s (all of them have to be categorical!)}
  \item{model }{ a model formula or an object of class "\code{glm}". }
}
\details{
  For the estimation of the PAR the partitioning technique through the interstep of the 
  sequential attributable risks by Eide and Gefeller (1995) is used. 
  
  It is assumed that all exposure factors are equally ranking and that there are no
  equally or hierarchically structured classes of exposure factors. 
  
  The needed (adjusted) attributable risks are estimated within the function with access to
  the function \code{\link{AdjAR}} (look there for further information). 
}

\note{

  Also if there are only a single exposure factor/confounder you have to enter a matrix, so
   this will be a matrix with only one column.
   
  It is also important that the given variables in \code{D}, \code{E} and \code{C} are not defined
  as factors.
    
  The names of the variables (outcome, exposure factors, confounders) in the argument \code{model} 
  have to be identically to the (column-)names of the entered data.  
  Furthermore all given exposure factors and confounders have to be part of the argument \code{model}.
  
  Validity of the estimation can only be taken for granted for data with simple random sampling, stratified
   random sampling or frequency-matching of controls.    
  
  To simplify the computation a compressed definition for the PAR is used (see Wille and Gefeller (1996) 
  for detailed information). 
}
\value{

  \code{PAR} returns a object of class "PAR". The partial attributable risk for 
  each exposure factor is given. To get the PAR/s for further calculations use function \code{getPAR}.
  
  
}
\references{ Eide, G.; Gefeller, O. (1995)
              Sequential and average attributable fractions as aids in the selection of 
              preventive strategies
               \emph{Journal of Clinical Epidemiology} \bold{48}, 645-55
             
             Wille, L.; Gefeller, O. (1996)
               Partitioning the disease risk among several exposure factors: a computational
               solution to an epidemiological problem
               \emph{Advances in Statistical Software} \bold{5}, 249-56 
               
             Land, M.; Vogel, C.; Gefeller, O. (2001)
              Partitioning methods for multifactorial risk attribution
               \emph{Statistical Methods in Medical Research} \bold{10}, 217-30   
            
             }
\author{ Christiane Raemsch }

\seealso{ \code{\link{AdjAR}}}
\examples{

#### partial attributable risks of exposure1    #####
#### and exposure2 if there are no confounders  #####

set.seed(2007)
dicho            <- c(0,1)
cc_state         <- sample(dicho, 100, replace=TRUE)
exposure1        <- sample(dicho, 100, replace=TRUE, prob=c(0.7, 0.3))
exposure2        <- sample(dicho, 100, replace=TRUE, prob=c(0.4, 0.6))
relation         <- as.formula(cc_state~exposure1+exposure2)
data_exp         <- cbind(exposure1, exposure2)
PAR_exposures    <- PAR(D=cc_state, E=data_exp, model=relation)
PAR_exposures

#### partial attributable risks of exposure1 and     #####
#### exposure2 with taking into account confounder1  #####

set.seed(2008)
dicho            <- c(0,1)
cc_state         <- sample(dicho, 100, replace=TRUE)
exposure1        <- sample(dicho, 100, replace=TRUE, prob=c(0.9, 0.1))
exposure2        <- sample(dicho, 100, replace=TRUE, prob=c(0.6, 0.4))
cat_confounder   <- c(0,1,2,3)
confounder1      <- sample(cat_confounder, 100, replace=TRUE)
data_exp         <- cbind(exposure1, exposure2)
conf             <- matrix(confounder1, ncol=1)
colnames(conf)   <- c("confounder1")
rel_mod          <- glm(cc_state~exposure1+exposure2+confounder1,
                        family=binomial)
PAR_exposures    <- PAR(cc_state, data_exp, conf, rel_mod)
PAR_exposures

}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{manip}

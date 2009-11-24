\name{AR-class}
\docType{class}
\alias{AR-class}
\alias{getAR,AR-method}
\alias{getCombi,AR-method}
\alias{getCProb,AR-method}
\alias{getOR,AR-method}
\alias{print,AR-method}
\alias{show,AR-method}

\title{Class "AR"}
\description{Presentation of attributable risks}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("AR", ...)}.
}
\section{Slots}{
	 \describe{
    \item{\code{AttRisk}:}{Attributable Risks for every combination of expsoure factors = Object of class \code{"matrix"}}
    \item{\code{OddsRatio}:}{Odds Ratio for every combination of expsoure factors = Object of class \code{"matrix"} }
    \item{\code{CondProb}:}{Conditional probabilities for every combination of expsoure factors = Object of class \code{"numeric"} }
    \item{\code{combinations}:}{every combination of expsoure factors = Object of class \code{"matrix"} }
  }
}
\section{Methods}{
  \describe{
    \item{getAR}{\code{signature(object = "AR")}: get a matrix with attributable risks }
    \item{getCombi}{\code{signature(object = "AR")}: get all combinations of exposure factors}
    \item{getCProb}{\code{signature(object = "AR")}: get conditional probabilities  }
    \item{getOR}{\code{signature(object = "AR")}: get odds ratios }
    \item{print}{\code{signature(x = "AR")}}
    \item{show}{\code{signature(object = "AR")}}
	 }
}

\note{ The slots OddsRatio, CondProb and combinations are only filled if confounders have not more than 2 categories. }

\author{Christiane Raemsch }


\seealso{
	Class \code{\linkS4class{PAR}} 
}
\examples{
showClass("AR")
}
\keyword{classes}

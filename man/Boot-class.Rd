\name{Boot-class}
\docType{class}
\alias{Boot-class}
\alias{getCIbca,Boot-method}
\alias{getCIperc,Boot-method}
\alias{print,Boot-method}
\alias{show,Boot-method}

\title{Class "Boot"}
\description{Presentation of bootstrap confidence intervals}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("Boot", ...)}.
}
\section{Slots}{
	 \describe{
    \item{\code{CI_perc}:}{Percentile confidence intervals = Object of class \code{"matrix"}}
    \item{\code{CI_bca}:}{BCa confidence intervals = Object of class \code{"matrix"}}
  }
}
\section{Methods}{
  \describe{
    \item{getCIbca}{\code{signature(object = "Boot")}: get BCa confidence intervals }
    \item{getCIperc}{\code{signature(object = "Boot")}: get percentile confidence intervals }
    \item{print}{\code{signature(x = "Boot")}}
    \item{show}{\code{signature(object = "Boot")}}
	 }
}

\author{Christiane Raemsch}

\seealso{
Class \code{\linkS4class{PAR}},  
Class \code{\linkS4class{AR}} 
}
\examples{
showClass("Boot")
}
\keyword{classes}

\name{PAR-class}
\docType{class}
\alias{PAR-class}
\alias{getPAR,PAR-method}
\alias{print,PAR-method}
\alias{show,PAR-method}

\title{Class "PAR"}
\description{Presentation of partial attributable risks}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("PAR", ...)}.
}
\section{Slots}{
	 \describe{
    \item{\code{PAR}:}{Object of class \code{"matrix"}}
  }
}
\section{Methods}{
  \describe{
    \item{getPAR}{\code{signature(object = "PAR")}:  get a matrix with the partial attributable risks }
    \item{print}{\code{signature(x = "PAR")} }
    \item{show}{\code{signature(object = "PAR")}}
	 }
}

\author{Christiane Raemsch}

\seealso{
Class	\code{\linkS4class{AR}} 
}
\examples{
showClass("PAR")
}
\keyword{classes}

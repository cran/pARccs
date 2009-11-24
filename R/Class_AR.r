########## Definition Class AR ########

setClass(Class="AR", representation(AttRisk="matrix", OddsRatio="matrix",CondProb="numeric", combinations="matrix"))

setMethod ("print","AR",
 function(x,...){
 cat("\n\n")
 cat("Attributable Risks\n")
 cat("********************")
 cat("\n\n")
 print (x@AttRisk)
 cat("\n")
 }
)

setMethod ("show","AR",
 function(object){
 cat("\n\n")
 cat("Attributable Risks \n=")
 cat("********************")
 cat("\n\n")
 print (object@AttRisk)
 cat("\n")
 }
)

setGeneric("getAR",function(object){standardGeneric ("getAR")})

setMethod("getAR","AR",
  function(object){
  return(object@AttRisk)
  }
)

setGeneric("getOR",function(object){standardGeneric ("getOR")})

setMethod("getOR","AR",
  function(object){
  return(object@OddsRatio)
  }
)

setGeneric("getCProb",function(object){standardGeneric ("getCProb")})

setMethod("getCProb","AR",
  function(object){
  return(object@CondProb)
  }
)

setGeneric("getCombi",function(object){standardGeneric ("getCombi")})

setMethod("getCombi","AR",
  function(object){
  return(object@combinations)
  }
)

########### Definition Class PAR ################

setClass(Class="PAR", representation(PAR="matrix"))

setMethod ("print","PAR",
 function(x,...){
 cat("\n\n")
 cat("Partial Attributable Risk/s\n")
 cat("********************")
 cat("\n\n")
 print (x@PAR)
 cat("\n")
 }
)

setMethod ("show","PAR",
 function(object){
 cat("\n\n")
 cat("Partial Attributable Risk/s \n=")
 cat("********************")
 cat("\n\n")
 print (object@PAR)
 cat("\n")
 }
)

setGeneric("getPAR",function(object){standardGeneric ("getPAR")})

setMethod("getPAR","PAR",
  function(object){
  return(object@PAR)
  }
)

########### Definition Class Boot ################

setClass(Class="Boot", representation(CI_perc="matrix", CI_bca="matrix"))

setMethod ("print","Boot",
 function(x,...){
 cat("\n\n")
 cat("Bootstrap Confidence Intervall/s\n")
 cat("**********************************")
 cat("\n\n")
 cat("Percentile Confidence Intervall:\n")
 cat("\n")
 print (x@CI_perc)
 cat("\n")
 cat("BCa Confidence Intervall:\n")
 cat("\n")
 print (x@CI_bca)
 cat("\n")
 }
)

setMethod ("show","Boot",
 function(object){
 cat("\n\n")
 cat("Bootstrap Confidence Intervall/s\n")
 cat("**********************************")
 cat("\n\n")
 cat("Percentile Confidence Intervall:\n")
 cat("\n")
 print (object@CI_perc)
 cat("\n")
 cat("BCa Confidence Intervall:\n")
 cat("\n")
 print (object@CI_bca)
 cat("\n")
 }
)

setGeneric("getCIperc",function(object){standardGeneric ("getCIperc")})

setMethod("getCIperc","Boot",
  function(object){
  return(object@CI_perc)
  }
)

setGeneric("getCIbca",function(object){standardGeneric ("getCIbca")})

setMethod("getCIbca","Boot",
  function(object){
  return(object@CI_bca)
  }
)
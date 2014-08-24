\name{constants-package}
\alias{constants-package}
\alias{constants}
\docType{package}
\title{
Constants
}
\description{
Defines some physical constants and functions commonly used in optics.
}
\details{
\tabular{ll}{
Package: \tab constants\cr
Type: \tab Package\cr
Version: \tab 1.0\cr
Date: \tab 2007-12-29\cr
License: \tab GPL\cr
}
Constants is a data frame with some physical constants, plus units and name in comment. Gold, silver and glass permittivities, together with spline interpolation functions. luxpop data for other compounds. N and k / epsilon conversion both ways. energy/ wavelength conversion. Fano and Lorentzian functions.
}
\author{
baptiste Auguie

Maintainer: baptiste Auguie <ba208@ex.ac.uk>
}
\references{
 NIST 
}
\keyword{ physical constants optical }
\seealso{
 \code{\link[<pkg>:<pkg>-package]{<pkg>}} 
}
\examples{
# library(constants)
# data(Constants)
# x<-seq(from=-2,to=3,by=.1)
# p<-c(0,1,1,1,1)
# y<-cbind(fano(p,x),lorentz(p[1:4],x))
# y<-y%*%diag(1/apply(y,2,"max")) # normalize to 1
# matplot(x,y)
# legend("topleft",c("fano","lorentz"),pch=c("1","2"),col=c(1,2))
# Constants$cel
# comment(Constants$cel)
# 
# data(luxpop) # loads of names...
# fit.material("Cr",all=FALSE,range=c(0.5,0.8))
# 
# ## or directly from luxpop
# 
# baseUrl="http://www.luxpop.com/Material/"
# material="Cr.nk"
# readLines(url(paste(baseUrl,material,sep=""), open = "", blocking = TRUE, encoding = getOption("encoding")))->file
# grep(";",file,value=TRUE)->comment
# comment
# read.table(textConnection(file),skip=length(comment))->values
# head(values)
# names(values)<-c("l","n","k")
# attach(values)
# nk2epsilon(list(n=n,k=k))->valuesEpsilon
# valuesEpsilon$lambda<-l/1e4 # Angstrom to um
# detach(values)
# attach(valuesEpsilon)
# 
# graphics.off()
# par(mfrow=c(2,2),mar=c(5, 5, 1, 1),bty="n",pty="m")
# cutePlot(lambda,epsr)
# cutePlot(lambda,epsi,yl=expression(paste(epsilon[i])))
# cutePlot(lambda,values$n,yl=expression(paste("n")))
# cutePlot(lambda,values$k,yl=expression(paste("k")))




}

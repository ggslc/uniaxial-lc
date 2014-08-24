\name{plot.material}
\alias{plot.material}
\alias{fit.material}
\alias{material2plot}

\title{spline fit and plotting material permittivities}
\description{
plots raw data, and fit, returns invisible raw data, and fine interpolation
}
\usage{
fit.material(material = AuJChristy, all = FALSE, draw= TRUE, range.lambda = c(0.2, 2), fine.n = 200)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{material}{ string }
  \item{all}{all range? }
  \item{draw}{plot? }
  \item{range.lambda}{ range }
  \item{fine.n}{ number of points in spline fit }
}
\details{
  
}
\value{
  list
  \item{raw.data }{lambda , epsilon, nk }
  \item{fine.sample }{lambda , epsilon, nk }
  ...
}
\references{ luxpop }
\author{ baptiste auguie }
\note{ 
}
\seealso{}
\examples{
# library(constants)
# data(luxpop)
# plot.material() -> test
# data(AuJChristy)
# plot.material("AuJChristy",range.lambda = c(0.4, 1), draw=FALSE) -> gold
# 
# subset<- sort(sample(length(gold$fine.sample$lambda), 20))
# small.sample <-as.data.frame(gold$fine.sample)[subset,]
# str(small.sample)
# with(small.sample,plot(lambda,Re(epsilon)))

\dontrun{
require(constants)
require(baptMisc)
require(ggplot2)
# clr()
data(AuJChristy) # load some permittivity values
gold <- fit.material(AuJChristy, range=c(0.4, 1), draw=F)

p0.basic <- c(1.5, 1.3, 1.54)

model.basic <- function(p, lambda=gold$fine.sample$lambda*1e-6){
	epsilon <- drude(lambda=lambda, interband=FALSE, basic=TRUE, surf=FALSE, 
		omega.p = p[1]*1e+16, 
		gamma.p = p[2]*1e+14,
		eps.inf=p[3])
		data.frame(lambda=lambda*1e6, epsilon=epsilon, nk=sqrt(epsilon))
}

obj.basic <- function(p, restrict=c(0.7, 1.2)){
 data.restrict <- subset(gold$fine.sample, isIn(restrict, lambda)  )
	model <- model.basic(p, lambda=data.restrict$lambda*1e-6)
	sos.real <- sum( (Re(model$epsilon) - Re(data.restrict$epsilon))^2 ) / max(abs(Re(data.restrict$epsilon)))
	sos.imag <- sum( (Im(model$epsilon) - Im(data.restrict$epsilon))^2 ) / max(abs(Im(data.restrict$epsilon)))
	sos <-  sos.real + sos.imag
	sos
}

res.basic <- optim(p0.basic, obj.basic)
fit.basic <- model.basic(res.basic$par)

mdrude <- material2plot(fit.basic)
gold$qplot+layer(data=mdrude$meltdf, mapping = aes(colour = variable), geom = c("line"), 
                stat = "identity")
}
}

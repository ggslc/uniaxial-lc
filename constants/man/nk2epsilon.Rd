\name{nk2epsilon}
\alias{nk2epsilon}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ refractive index to permittivity  }
\description{
  Converts n and k (list) to epsilon (list) 
}
\usage{
nk2epsilon(nk)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{nk}{ list(n=, k=)}
}
\details{
 
}
\value{
list(epsr,epsi)

}
\references{   }
\author{ baptiste }
\note{ 
}
\seealso{ reverse of: \code{\link{epsilon2nk}},  }
\examples{

## The function is currently defined as
function(nk)
{
epsi<-Im((nk$n+1i*nk$k)^2)
epsr<-Re((nk$n+1i*nk$k)^2)
list(epsr=epsr,epsi=epsi)
	}
}


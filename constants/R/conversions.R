`epsilon2nk`<-
function(epsilon)
{
nk <- sqrt(epsilon)
}

nm2um <- 1e-3
um2nm <- 1e3

`L2eV` <- 
function(lambda)
{
#with(Constants, h * cel / ee / lambda)
cel<-2.99792458e8 
h<-6.62606896e-34 
ee<-1.602176487e-19 
h * cel / ee / lambda -> eV
eV
	}

`eV2L` <- function(eV)
{
## with(constants::Constants, h * cel / ee / eV)

cel <- 2.99792458e8 
h <- 6.62606896e-34 
ee <- 1.602176487e-19 

h * cel / ee / eV

}

`L2w` <- 
function(lambda)
{
  ##with(Constants, 2*pi * cel / lambda)
  cel <- 2.99792458e8
  2*pi * cel / lambda 
}


`w2L` <- 
function(omega)
{
## with(Constants, 2*pi * cel / omega )
cel <- 2.99792458e8 
2*pi * cel / omega 
}

t2eV <- function(tau){
##   with(Constants, h/(pi*tau*ee))
  h <- 6.62606896e-34 
  ee <- 1.602176487e-19 
  h/(pi*tau*ee)
}


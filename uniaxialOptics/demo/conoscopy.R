require(uniaxialOptics)
#compute a conoscopy figure of a uniaxial, homeotropic
#layer between ITO coated glass
#shows the difference bewteen Ko/Berreman and Lien/Jones
#methods. Ko/Berreman includes thin film inteference.
#simple stack, one birefringent, near homeotropic layer between
#ito coated glass

depth <- c(0,.03,30,.03,0) * 1e-6
glass <- 1.52^2 + 0.0i
ito <- 3.8 + 0.08i
eo <- c(glass,ito,glass,ito,glass)
ee <- c(glass,ito,(sqrt(Re(glass)) + 0.2)^2 + 0i,ito,glass)
tilt <- c(0,0,0,0,0)
twist <- rep(0,length(tilt))

angle <- seq(0,pi/3,l=100)
azimuth <- rep(pi/4,length(angle))
lambda <- c(600e-9,400e-9)

numap <- .55/1.52 # numerical aperture, in glass, of ~50X lens
n <- 96 #horizontal and vertical resolution




system.time(rko <- conorect.stack(lambda,depth,tilt,twist,eo,ee,n,n, 0,
                      numap,az0=pi/4,method="ko"))




system.time(rlien <- conorect.stack(lambda,depth,tilt,twist,eo,ee,n,n, 0,
                      numap,az0=pi/4,method="lien"))







#define gamma corrected image plotting function 
imp <- function(r,l,bright,...){

  s <- seq(0,1,l=64)^(1/2.2)
  b <- col2rgb(bright)
  col <- rgb(b[1]*s,b[2]*s,b[3]*s,maxColorValue=255)
  
  image((r$crossed[,,l]),col=col,axes=FALSE,...)
}

par(mfrow=c(2,2))
imp(rko,1,bright="red",main=paste(lambda[1]*1e+9," nm, ko"))
imp(rko,2,bright="blue",main=paste(lambda[2]*1e+9," nm, ko"))
imp(rlien,1,bright="red",main=paste(lambda[1]*1e+9," nm, lien"))
imp(rlien,2,bright="blue",main=paste(lambda[2]*1e+9," nm, lien"))

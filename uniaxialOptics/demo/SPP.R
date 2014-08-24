#surface plasma polariton/ Fabry-Perot
#angle scan calculation
#using Berreman method and dsipersion data from Baptiste
#Augie's constants package.
#System is 50nm Ag,5um 5CB, 50nm Ag, bounded by Silica.

require(uniaxialOptics)
require(constants)
data(AgPalik)

a <- 85 * pi/180
po <- seq(a/2,a,l=500)
az <- rep(0,length(po))

#select a few wavlengths
sq <- (AgPalik$lambda > 450e-9 & AgPalik$lambda < 650e-9)
nl <- length(AgPalik$lambda[sq])
                             
d <- c(0,50e-9,5e-6,50e-9,0)

#express permittivity as matrix
esi <- rep(3.24 + 0i,nl)
eo5CB <- 1.52^2 + 0.0i
ee5CB <- 1.68^2 + 0.00i
eAg <- AgPalik$epsr[sq] + 1i * AgPalik$epsi[sq]

eo <- rbind(esi,eAg,eo5CB,eAg,esi)
ee <- eo
ee[3,] <- ee5CB


of <- function(tilt){

  nt <- length(tilt)
  optics.stack(AgPalik$lambda[sq] ,az,po,
               c(0,50e-9,rep(5e-6/nt,nt),50e-9,0),
                 c(0,0,tilt,0,0),
                 c(0,0,rep(0,nt),0,0),
                 rbind(esi,eAg,
                       rbind(matrix(rep(eo5CB,nl*nt),nrow=nt)),
                       eAg,esi),
                 rbind(esi,eAg,
                       rbind(matrix(rep(ee5CB,nl*nt),nrow=nt)),
                       eAg,esi))
}

o <- of(pi/2)

dt <- pi/6
o2 <- of(c(pi/2,pi/2-dt,pi/2)) # change mid-cell - barely changes SPP
o3 <- of(c(pi/2-dt,pi/2,pi/2)) # change at lower surface - massive change to SPP
#

p <- 3
r <- rev(rainbow(length(o$lambda)+p))[-(1:p)]

matplot(o$angle*180/pi,
        o$Rpp + rep(seq(0,4,l=nl),each=length(o$angle)),
        type='l',lty=1,main="Fabry-Perot modes and SPPs",
        col=r,xlab=expression("Incident angle, " * alpha * " / degrees"),
        ylab=expression("Reflectivity, " * R[pp] * " (offset)"))


matlines(o2$angle*180/pi,
        o2$Rpp + rep(seq(0,4,l=nl),each=length(o$angle)),
        type='l',lty=2,col=r)

matlines(o3$angle*180/pi,
        o3$Rpp + rep(seq(0,4,l=nl),each=length(o$angle)),
        type='l',lty=3,col=r,lwd=1)


legend(x="bottomleft",lwd=1,leg=o$lambda*1e9,
       col=r,
       title=expression(lambda * phantom(0) / phantom(0) * nm)
       ,bg="white")


legend(x="bottomright",lwd=1,
       leg=c("planar","tilted at miplane","tilted at surface"),
       lty=1:3,title=expression(theta(z)),
       ,bg="white")

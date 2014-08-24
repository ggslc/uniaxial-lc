require(uniaxialOptics)

#plot conoscopy figures for some typical LC cells

#5CB LC films between ITO coated glass
glass <- 1.52^2 + 0.0i
ito <- 3.8 + 0.08i 
nz <- 48   #divide lc into sublayers
d <- 30e-6 #lc layer thickness
eolc <- 1.52^2 + 0i
eelc <- 1.68^2 + 0i
di <- c(0,30e-9)
ei <- c(glass,ito)
numap <- asin(.55/1.52) #50X refracted by glass

#CCD plot function
ccdf <- function(lc,main) {
  ccd <- conorect.lc(633e-9,lc,eolc, eelc, di, ei, rev(di), rev(ei),
              nx=96,ny=96,pmin=0,pmax=numap)
  image((ccd$crossed[,,1]),col=grey.colors(64,start=0,end=1),axes=FALSE,
        main=main)
  
}


par(mfrow=c(2,3))

#han cell
#note that the lc profiles can be computed by
#the fucntion lc_switching_seqeunce in package
#Ericksen-Leslie, but are done on the cheap here.

lch <- list(depthz=c(0.5,rep(1,nz-2),0.5) * d/nz,
              tilt=matrix(seq(0,pi/2,l=nz),ncol=1),
              twist=matrix(seq(0,0,l=nz),ncol=1))
ccdf(lch,"HAN cell" )

#planar cell
lcp <- list(depthz=c(0.5,rep(1,nz-2),0.5) * d/nz,
              tilt=matrix(seq(pi/2,pi/2,l=nz),ncol=1),
              twist=matrix(seq(0,0,l=nz),ncol=1))
ccdf(lcp,"Planar homogeneous cell")

#pi/16 slab cell
lctn <- list(depthz=c(0.5,rep(1,nz-2),0.5) * d/nz,
              tilt=matrix(seq(pi/16,pi/16,l=nz),ncol=1),
              twist=matrix(seq(0,0,l=nz),ncol=1))
ccdf(lctn,expression(pi/16 * " slab cell"))


#homeotropic  cell
lcht <- list(depthz=c(0.5,rep(1,nz-2),0.5) * d/nz,
              tilt=matrix(seq(0,0,l=nz),ncol=1),
              twist=matrix(seq(0,0,l=nz),ncol=1))
ccdf(lcht,"Homeotropic cell")

#pi/32midplane symettric cell
lcs <- list(depthz=c(0.5,rep(1,nz-2),0.5) * d/nz,
              tilt=matrix( ifelse(1:nz > nz/2, pi/32, -pi/32),ncol=1),
              twist=matrix(seq(0,0,l=nz),ncol=1))
ccdf(lcs,expression(pi/16 * " midplane symmetric cell"))


#tn cell
lctn <- list(depthz=c(0.5,rep(1,nz-2),0.5) * d/nz,
              tilt=matrix(seq(pi/2,pi/2,l=nz),ncol=1),
              twist=matrix(seq(0,pi/2,l=nz),ncol=1))
ccdf(lctn,"twisted cell")

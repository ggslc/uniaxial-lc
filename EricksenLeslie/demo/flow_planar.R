require(EricksenLeslie)

lcf <- function(pg,az,depth,
                k,alpha,epsilon) {

  #returns
  #-------
  #list including director profile
  #for a pressure gradient along the x axis
  #director is planar homogeneous, same rubbing drection at both surfaces
  
  #arguments
  #---------
  #pg:    vector, arbitrary length, pressure gradient, Pa/m
  #az:    azimuthal angle between pg and director
  #depth: layer depth
  #k,alpha,epsilon: lc parameters

  #notes
  #-----
  #Finding steady-state solutions for
  #planar cells is unstable when the
  #flow causes a tilt of more than a couple of degrees.
  #To get round this, the calculations are time-dependent,
  #so the first element of pg should be low, or zero
  
  aza <- 0
  
  lc_switching_sequence(nz = 48, times = seq(0,10,l=length(pg)),
                        voltages = rep(0,length(pg)),
                        freq = rep(1,length(pg)), depth = depth,
                        k = k , epsilon = epsilon,
                        tilt_easy_angle = c(pi/2,pi/2),
                        twist_easy_angle=c(az,az),
                        pgx = pg,
                        alpha = alpha, shear = TRUE,
                        niter=20,ntol=1e-6,ttol=1e-5)
  
}


lcv <- function(U,az,depth,k,alpha,epsilon){

  #returns
  #-------
  #list including director profile
  #for a flow rate U along the x axis
  #director is planar homogeneous, same rubbing drection at both surfaces
  
  #arguments
  #---------
  #U:    vector, arbitray length, volume flux (per square meter)
  #az:    azimuthal angle between pg and director
  #depth: layer depth
  #k,alpha,epsilon: lc parameters

  #notes
  #-----
  #unless the Ericksen-Leslie equations are linearised, it
  #the flux must be specified indirectly, by
  #finding a opressure gradient along x such
  #that integral(u(z, pg) dz) = U

  f <- function(pg,U){
    #difference between computed and desired flux
    #for a pressure gradient pg
    lc <- lcf(c(0,pg),az,depth,k,alpha,epsilon)
    sum(lc$flowx[,2] * lc$depthz) - U
  }

  p <- function(U,etac){
    #pressure gradient for given U 
    if (U != 0){
      #maximum pressure gradient
      pgmax <- -2 * 12 * etac * U / depth^3
      uniroot(function(pg)f(pg,U),c(min(0,pgmax),max(0,pgmax)))$root
    } else {
      0
    }
  }
  
  #made up effective shear viscosity etac for
  #limits of pressure gradiant 
  eta2 <-  0.5 * (alpha[2] + 2 * alpha[3] + alpha[4] + alpha[5])   
  etac <- 3 * eta2
  
  pg <- NULL
  for (Ue in U){
    pg <- c(pg,p(Ue,etac))
  }

  #lazy... but convenient - recompute
  #all the lc profiles over pg
  lcf(pg,az,depth,k,alpha,epsilon)
  
  
}

#define matrix plotting function for convenience
mp <- function(x,y,...){
matplot(x,y,type='l',lty=1,
col=colorRampPalette(c("blue","orange"),space="Lab")(ncol(y)),...)
}

#parameters for 5CB
k <- c(.62,.39,.89) * 1e-11
alpha <- c(-0.006,-0.0824,-0.0036,0.0652,0.0640,0)
epsilon <- c(7,11.5)

#20 microns lc layer
depth <- 20e-6
#cell width, m
dy <- 3e-3

#flux, m^2 / s

cfac <- 1e-6 * 1e-3 / dy / 3600
#       mul to l | l to m^3 | mean in y | hours to seconds

par(mfrow=c(2,3))

#flow parallel to rubbing direction
#high flow rates required to see much at all
#flow rate in ul/h
urate <- seq(0,100, l = 10)
U <- urate * cfac
lc <- lcv(U,0,depth,k,alpha,epsilon)

#plot the tilt profiles
mp(lc$z * 1e6,lc$tilt * 180 / pi,
   xlab=expression(z * " / " * mu * m),
   ylab=expression(theta * " / " * degrees))
   

#plot the maximum tilt against flow rate
plot(urate,apply(lc$tilt,2,max)*180/pi,type='b',
     xlab=expression('Requested flux, ' * U * " / " *   mu * l * h^{-1}),
     ylab=expression(theta[max] * " / " * degrees))
                                        
#plot computed flux against requested flux 
plot(urate,apply(lc$flowx * lc$depthz,2,sum)/cfac,type='b',
     xlab=expression('Requested flux, ' * U * " / " *   mu * l * h^{-1}),
     ylab=expression('Computed flux, ' * U[c] * " / " *   mu * l * h^{-1}))


#flow at pi/4 to rubbing direction
#lower flow rates needed - director twists out quickly
urate <- seq(0,40, l = 10)
U <- urate * cfac

lc <- lcv(U,pi/4,depth,k,alpha,epsilon)

#plot the tilt profiles
mp(lc$z * 1e6,lc$tilt * 180 / pi,
   xlab=expression(z * " / " * mu * m),
   ylab=expression(theta * " / " * degrees))
   
#plot the twist profiles
mp(lc$z * 1e6,lc$twist * 180 / pi,
   xlab=expression(z * " / " * mu * m),
   ylab=expression(theta * " / " * degrees))


#plot the minimum twist against flow rate
plot(urate,apply(lc$twist,2,min)*180/pi,type='b',
     xlab=expression('Requested flux, ' * U * " / " *   mu * l * h^{-1}),
     ylab=expression(phi[max] * " / " * degrees))#

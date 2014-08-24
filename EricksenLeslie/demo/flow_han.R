require(EricksenLeslie)

lcf <- function(pg,az,depth,
                k,alpha,epsilon,v=0) {

  #returns
  #-------
  #list including director profile
  #for a pressure gradient along the x axis
  #HAN cell
  
  #arguments
  #---------
  #pg:    vector, arbitrary length, pressure gradient, Pa/m
  #az:    azimuthal angle between pg and director
  #depth: layer depth
  #k,alpha,epsilon: lc parameters

 
  
  aza <- 0
  
  lc_switching_sequence(nz = 48, times = seq(0,0,l=length(pg)),
                        voltages = rep(v,length(pg)),
                        freq = rep(1,length(pg)), depth = depth,
                        k = k , epsilon = epsilon,
                        tilt_easy_angle = c(0,pi/2),
                        twist_easy_angle=c(az,az),
                        pgx = pg,
                        alpha = alpha, shear = TRUE,
                        niter=20,ntol=1e-6,ttol=1e-5)
  
}


lcv <- function(U,az,depth,k,alpha,epsilon,v=0){

  #returns
  #-------
  #list including director profile
  #for a flow rate U along the x axis
  #HAN cell
  
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
    lc <- lcf(pg,az,depth,k,alpha,epsilon,v=v)
    sum(lc$flowx * lc$depthz) - U
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
  lcf(pg,az,depth,k,alpha,epsilon,v=v)
  
  
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

#flow parallel or anti-parallel to rubbing direction
#flow rate in ul/h
urate <- seq(-5,10, l = 10)
U <- urate * cfac
lc <- lcv(U,0,depth,k,alpha,epsilon,v=0)

#plot the tilt profiles
mp(lc$z * 1e6,lc$tilt * 180 / pi,
   xlab=expression(z * " / " * mu * m),
   ylab=expression(theta * " / " * degrees))
   
                                        
#plot computed flux against requested flux 
plot(urate,apply(lc$flowx * lc$depthz,2,sum)/cfac,type='b',
     xlab=expression('Requested flux, ' * U * " / " *   mu * l * h^{-1}),
     ylab=expression('Computed flux, ' * U[c] * " / " *   mu * l * h^{-1}))


#applied field
lc <- lcv(U,0,depth,k,alpha,epsilon,v=5)

#plot the tilt profiles
mp(lc$z * 1e6,lc$tilt * 180 / pi,
   xlab=expression(z * " / " * mu * m),
   ylab=expression(theta * " / " * degrees))

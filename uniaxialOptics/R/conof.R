#  conof.R: some functions to compute conoscopy
#  figures for anisotropic multilayer stacks
#  copyright (C) 2007  Steph Cornford
#
#  This program is free software; you can redistribute it and/or
#  modify it under the terms of the GNU General Public License
#  as published by the Free Software Foundation; either version 2
#  of the License, or (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  You should have received a copy of the GNU General Public License
#  along with this program; if not, write to the Free Software
#  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.


conopix.stack <- function(lambda,depth,tilt,twist,eo,ee,x,y,az0=pi/4,method="ko"){

  #compute the conoscopic intensity at pixel pairs (x,y)
  #lambda - wavelength of light
  #depth - vector, n long, of layer depths, d[1]=d[n]  = 0
  #tilt,twist - vector, n long, of tilt and twist angles
  #eo,ee - ordinary, extraodinary permitivities 
  #        (complex, perhaps n x m matrices)
  #x,y   - pixels to compute
  #az0   - inclination of polarizer transmission axis to stack y axis
  
  r <- sqrt(x^2 + y^2)
  po <- ifelse(r < 1, asin(r), 1)
  az <- ifelse(x > 0 & y > 0, atan(y/x) ,
                ifelse(x < 0, pi + atan(y/x),
                       ifelse(x > 0 & y < 0, 2 * pi + atan(y/x),
                              pi/2 * sign(y))))
  
 	
  o <- optics.stack(lambda,az,po,depth,tilt,twist,eo,ee,method)
  
  p <- cos(az - az0)
  s <- -sin(az - az0)

  crossed <- Mod((o$tpp - o$tss) * p * s - o$tps * p * p + o$tsp * s * s)^2
  
  dim(crossed) <- c(length(x),length(lambda))
  
  list(x=x,y=y,crossed=crossed,lambda=lambda)
  

}





conorect.stack <- function(lambda,depth,tilt,twist,eo,ee,nx,ny,
                     pmin,pmax,az0=pi/4,method="ko"){

  #conoscopic image as a function of x,y (ie CCD image)
  #regions outside [pmin,pmax] are filled with zeros

  x <- seq(sin(-pmax),sin(pmax),l=nx)
  y <- seq(sin(-pmax),sin(pmax),l=ny)
  
  xx <- rep(x,each=length(y))
  yy <- rep(y,length(x))
  
  s <- sqrt(xx^2 + yy^2) <= sin(pmax) & sqrt(xx^2 + yy^2) >= sin(pmin)
  
  cotmp <- conopix.stack(lambda,depth,tilt,twist,eo,ee,xx,yy,az0,method)
  #although it seems ineffcient to compute pixels beyond s
  #when they will be overwritten with zeros, it is quicker
  #than unpacking the compressed data returned by
  #cotmp <- conopix.stack(lambda,depth,tilt,twist,eo,ee,xx[s],yy[s],az0,method)
  
  co <- list(x=xx,y=yy,crossed=NULL)
  
  for (j in 1:dim(cotmp$crossed)[2]){
    co$crossed <- c(co$crossed,ifelse(s,cotmp$crossed[,j],1))
     
  }
  
  dim(co$crossed) <- c(nx,ny,dim(cotmp$crossed)[2])

  co$pmax <- pmax
  co$pmin <- pmin
  co$nx <- nx
  co$ny <- ny
  co

}




conopix.lc <- function(lambda,lc,eolc,eelc,db,eb,dt,et,x,y,az0=pi/4,i=1,method="ko"){

  #compute conocopy pixels for time i of lc
  s <- .stack.lc(lc,eolc,eelc,db,eb,dt,et,i)
  conopix.stack(lambda,s$depth,s$tilt,s$twist,s$eo,s$ee,x,y,az0,method)
  
}


conorect.lc <- function(lambda,lc,eolc,eelc,db,eb,dt,et,nx,ny,pmin,pmax,az0=pi/4,i=1,method="ko"){

  #compute conocopy rectangle(CCD image) for time i of lc
  s <- .stack.lc(lc,eolc,eelc,db,eb,dt,et,i)
  conorect.stack(lambda,s$depth,s$tilt,s$twist,s$eo,s$ee,nx,ny,pmin,pmax,az0,method)
  
}








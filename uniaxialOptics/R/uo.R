#  uo.R: uniaxial optics calculation
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


.ko4x4 <- function(lambda,azimuth,angle,h,tilt,twist,eo,ee){
  
  nlayer <- length(h)
  nangle <- length(angle)

  
  d <- .Fortran("ainter",nlayer=as.integer(nlayer),nangle=as.integer(nangle)
                ,eo=as.complex(eo),ee=as.complex(ee)
                ,twist=as.numeric(twist),tilt=as.numeric(tilt),h=h,lambda=lambda
                ,coeffs=matrix(0,8,nangle),zcoeffs=matrix(0+0i,8,nangle)
                ,angle=as.numeric(angle),aziumuth=as.numeric(azimuth)
                ,PACKAGE="uniaxialOptics")

  i <- 1
  for (sig in c("Rpp","Rps","Rsp","Rss","Tpp","Tps","Tsp","Tss"))
    {
      d[[sig]] <- d$coeffs[i,]
      i <- i + 1
    }
  i <- 1
  for (sig in c("rpp","rps","rsp","rss","tpp","tps","tsp","tss"))
    {
      d[[sig]] <- d$zcoeffs[i,]
      i <- i + 1
    } 

  d$coeffs <- NULL
  d$zcoeffs <- NULL
  d
  
}

.lien2x2 <- function(lambda,azimuth,angle,h,tilt,twist,eo,ee){

  nlayer <- length(h)
  nangle <- length(angle)
  
   d <- .C("lien_2x2_optics",
           Tpp = numeric(nangle),
           Tps = numeric(nangle),
           Tsp = numeric(nangle),
           Tss = numeric(nangle),
           tpp = complex(nangle),
           tps = complex(nangle),
           tsp = complex(nangle),
           tss = complex(nangle),
           nangle = as.integer(nangle),
           angle = as.numeric(angle),
           azimuth = as.numeric(azimuth),
           lambda = as.numeric(lambda),
           nlayer = as.integer(nlayer),
           h = as.numeric(h),
           tilt = as.numeric(tilt),
           twist = as.numeric(twist),
           eo = as.complex(eo),
           eo = as.complex(ee),
           PACKAGE="uniaxialOptics")
           

  d$Rpp = rep(0,nangle)
  d$Rps = d$Rpp 
  d$Rsp = d$Rpp
  d$Rss = d$Rpp
  d$rpp = rep(0 + 0i,nangle)
  d$rps = d$rpp
  d$rsp = d$rpp
  d$rss = d$rpp

  d
}


optics.stack <- function(lambda,azimuth,angle,depth,tilt,twist,eo,ee, method="ko" )
{
  #compute the reflection and transmission
  #coefficients for a stratified, uniaxial system,
  #for each combination of lambda and the doublet (azimuth,angle)

  if (!is.numeric(angle)) stop ("angle must be numeric vector")
  if (!is.numeric(depth)) stop ("depth must be numeric vector")
  if (!is.complex(eo)) stop ("epsilon_1 must be complex vector")
  if (!is.complex(ee)) stop ("epsilon_3 must be complex vector")
  if (!is.numeric(tilt)) stop ("tilt must be numeric vector")
  if (!is.numeric(twist)) stop ("twist must be numeric vector")
  if (!is.numeric(lambda)) stop ("lambda must be numeric vector")
  if (!is.numeric(azimuth)) stop ("lambda must be numeric vector")

  if (!identical(length(depth),length(tilt))) stop ("length(tilt) != length(depth)")
  if (length(depth) != length(twist)) stop ("length(twist) != length(depth)") 

 
  if (length(angle) != length(azimuth)) stop ("length(angle) != length(azimuth)") 

  if (length(depth) != length(eo)
      & length(depth)  * length(lambda) != length(eo) )
    stop ("eo / depth mismatch") 

  if (length(depth) != length(ee)
     & length(depth)  * length(lambda) != length(ee))
    stop ("ee/ depth mismatch")  

  r <- list(outer=0,lambda=lambda,azimuth=azimuth,angle=angle) 
  siga <-  c("Rpp","Rps","Rsp","Rss","Tpp","Tps","Tsp","Tss",
             "rpp","rps","rsp","rss","tpp","tps","tsp","tss")
  
  first <- TRUE

  of <- if (method == "ko") {
    twist <- pi/2 - twist
    azimuth <- pi/2 - azimuth
    .ko4x4
  } else {
    if (method == "lien") {
      .lien2x2
    } else {
      stop (paste("unknown method",method)) 
    }
  }

  
  for (l in 1:length(lambda)) {
  #loop over wavelengths
    
    d <- of(lambda[l],azimuth,angle,depth,tilt,twist,
            if (is.matrix(eo)) {eo[,l]} else {eo},
            if (is.matrix(ee)) {ee[,l]} else {ee})


    if (first){
      first <- FALSE
      for (sig in siga)
        r[[sig]] <- d[[sig]]
    } else {
      for (sig in siga)
        r[[sig]] <- append(r[[sig]],d[[sig]])
    }
  }
    
    for (sig in siga){
      dim(r[[sig]]) <-  c(length(r$angle),length(r$lambda))
    }

r   

}

.stack.lc <- function(lc,eolc,eelc,db,eb,dt,et,i=1){
  #build a multilayer stack from lc at time j

  nz <- dim(lc$tilt)[1]
  
  ef <- function(eb,et,elc,nz){
    if (is.matrix(eb)){
      if (!is.matrix(et) | !is.matrix(elc)) stop ("dispersion must be applied to all epsilon or none")
      rbind(eb,matrix(rep(elc,nz),ncol=ncol(elc),byrow=TRUE),et)
    } else {
      if (is.matrix(eb) | is.matrix(elc)) stop ("dispersion must be applied to all epsilon or none")
      c(eb,rep(elc,nz),et)
    }
  }

  list(depth=c(db,lc$depthz,dt),
       tilt=c(rep(0.0,length(db)),
         (2*lc$tilt[1,i] + lc$tilt[2,i])/3,
         lc$tilt[2:(nz-1),i],
         (2*lc$tilt[nz,i] + lc$tilt[nz-1,i])/3,
         rep(0.0,length(dt))),
       twist=c(rep(0.0,length(db)),
         (2*lc$twist[1,i] + lc$twist[2,i])/3,
         lc$twist[2:(nz-1),i],
         (2*lc$twist[nz,i] + lc$twist[nz-1,i])/3,
         rep(0.0,length(dt))),
       eo=ef(eb,et,eolc,nz),
       ee=ef(eb,et,eelc,nz))
  
  
}


optics.lc <- function(lambda,azimuth,angle,lc,eolc,eelc,db,eb,dt,et,method="ko"){

  #compute transmission and reflection coefficients for a series of lc states
  #for each combination of lambda and the doublet (angle,azimuth)
  #The list lc must contain (at least) matrices tilt and twist, and vector depthz.
  #Each column of tilt (and twist) defines director profile, with each
  #element giving the value at the center of a lyer of depth depthz.
  #The first and last elements give the value at the boundary.
  #Just such output is produced by the lc_switching_sequence in the EricksenLeslie
  #package.

  i <- 1

   siga <-  c("Rpp","Rps","Rsp","Rss","Tpp","Tps","Tsp","Tss",
              "rpp","rps","rsp","rss","tpp","tps","tsp","tss")
  

  first <- TRUE
  nz <- nrow(lc$tilt)
  
  for (time in 1:ncol(lc$tilt)){

    s <- .stack.lc(lc,eolc,eelc,db,eb,dt,et,time)
    
    rt <- optics.stack(lambda,azimuth,angle,
                          depth=s$depth,tilt=s$tilt,
                          twist=s$twist,eo=s$eo,
                          ee=s$ee,method)
    
    
    if (first){
      r <- rt
      r$outer <- 1:(ncol(lc$tilt))
      first <- FALSE
    } else {
      for (sig in siga){
        r[[sig]] <- cbind(r[[sig]],rt[[sig]])
      }
      
    }
  }

  r
}
  
  
  





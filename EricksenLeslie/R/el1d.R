#
#  el1d.R: numerical solutions to Ericksen-Leslie equations for 1D
#  liquid crystal films: R interface
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


meisowicz <-function (alpha){

  list(gamma=c(alpha[3]-alpha[2],
               alpha[3]+alpha[2]),
       eta=c(0.5*(alpha[4] + alpha[5] - alpha[2]),
         	0.5 * (alpha[2] + 2.0 * alpha[3] + alpha[4] + alpha[5]),
         	0.5 * alpha[4],
         	alpha[1]))
}



lc_switching_sequence <- function
(nz, times, voltages, frequencies, depth, k, alpha, epsilon,
 tilt_easy_angle, tilt_anchoring=c(0,0),
 twist_easy_angle=c(0,0),twist_anchoring=c(0,0),
 flexo=c(0,0), pgx=rep(0,length(times)),pgy=rep(0,length(times)),
 rho0=0,mu=rep(0,length(rho0)),
 chi=rep(0,length(rho0)), charge=rep(0,length(rho0)), sigma0=rep(0,length(rho0)*2),
 shear=TRUE, shear_y = shear, initial=NULL,
 ntol=1e-5,niter=12,ttol=2e-3,tmin=1e-6,ftime=1.0)
{
  # compute a sequence of lcd states.
  # the vector times specifies the times at which states
  # are required. If times[i] = time[i-1] then state[i]
  # will be a steady-state solution. Otherwise, state[i]
  # will be a state evolved from state[i-1]

  if (nz < 3){
    stop("nz < 3")
  }
  
  if (is.null(times)){
    stop("NULL times") 
  }
  
  nt <- length(times)
  
  if (length(voltages) != nt){
    stop("length(voltages) != length(times)")
  }
  
  if (length(frequencies) != nt){
    stop("length(frequencies) != length(times)")
  }

  if (length(pgx) != nt){
    stop("length(pgx) != length(times)")
  }

    if (length(pgy) != nt){
    stop("length(pgy) != length(times)")
  }
  
  if (depth <= 0){
    stop("depth <= 0") 
  }
      
  if (length(k) != 3){
    stop("length(k) != 3")
  }

  if (! identical(rep(TRUE,3),k > 0)){
    stop("k < 0")
  }

  if (length(alpha) != 6){
    stop("length(alpha) != 6")
  }

  #TODO check alpha

  #append null boundary epsilon if needed
  if (length(epsilon) == 2){
    epsilon = c(epsilon,0,0)
  }
  
  if (length(epsilon) != 4){
    stop("length(epsilon) != 4")
  }

  if ((epsilon[1] + epsilon[2]) < 0){
    stop("delta epsilon < 0")
  }


  if (any(rho0 < 0)) {
    stop("negative rho0")
  }

  if (any(sigma0 < 0)) {
    stop("negative sigma0")
  }
  
  nsp <- length(rho0)
  spo <- 1:nsp
  sp <- function(v,each=1) v
  
  if (length(mu) < nsp){
    stop("length(mu) < length(rho0[rho0>0]")
  }
  if (length(chi) < nsp){
    stop("length(chi) < length(rho0[rho0>0]")
  }
    
  if (length(charge) < nsp){
    stop("length(charge) < length(rho0[rho0>0]")
  }
    
  if (length(sigma0) < 2*nsp){
    stop("length(sigma0) < 2 * length(rho0[rho0>0]")
  }

  if (is.null(initial)){

    #default initial guesses
    initial <- list(tilt = matrix( if (epsilon[2] > 0){0} else {pi/2} ,nz,nt),
                    twist = matrix(0,nz,nt),volt_ac=matrix(0,nz,nt),volt_dc=matrix(0,nz,nt),
                    flowx=matrix(0,nz,nt),flowy=matrix(0,nz,nt))  
                      
                    
  }

  #check size of initial guess matrices
  if ( any(dim(initial$tilt) != c(nz,nt))
      | any(dim(initial$twist) != c(nz,nt))
      | any(dim(initial$volt_ac) != c(nz,nt))
      | any(dim(initial$volt_dc) != c(nz,nt))
      | any(dim(initial$flowx) != c(nz,nt))
      | any(dim(initial$flowx) != c(nz,nt))){

    stop("if an initial lc profiles is given then all its elements (tilt,twist,volt_ac,volt_dc,flowx,flowy) must by nz * nt matrices")
  }
      
  .C("lcd_dynamics",
     nz=as.integer(nz),
     nt=as.integer(nt),
     times=as.double(times),
     voltages=as.double(voltages),
     frequencies=as.double(frequencies),
     depth=as.double(depth),
     k=as.double(k),
     alpha=as.double(alpha),
     shear=as.logical(shear),shear_y =as.logical(shear_y),
     epsilon=as.double(epsilon),
     flexo=as.double(flexo),
     tilt_anchoring=as.double(tilt_anchoring),
     tilt_easy_angle=as.double(tilt_easy_angle),
     twist_anchoring=as.double(twist_anchoring),
     twist_easy_angle=as.double(twist_easy_angle),
     pgx=as.double(pgx),pgy=as.double(pgy),
     nsp=as.integer(nsp),
     rho0=as.double(sp(rho0)),
     mu=as.double(sp(mu)),
     chi=as.double(sp(chi)),
     charge=as.double(sp(charge)),sigma0=as.double(matrix(sigma0,nrow=2)[,spo]),
     z=numeric(nz),
     depthz=numeric(nz),
     tilt=initial$tilt,
     twist=initial$twist,
     volt_ac=initial$volt_ac,
     volt_dc=initial$volt_dc,
     flowx=initial$flowx,
     flowy=initial$flowy,
     rho=array(0,c(nz,nt,nsp)),sigma=array(0,c(2,nt,nsp)),
     niter=as.integer(niter),ntol=as.double(ntol), tmin=as.double(tmin),
     ttol=as.double(ttol),ftime=as.double(ftime),iters=integer(nt),
     PACKAGE="EricksenLeslie")
}




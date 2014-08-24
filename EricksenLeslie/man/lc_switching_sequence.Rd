\name{lc_switching_sequence}
\alias{lc_switching_sequence}
\title{LC switching simulation}
\description{Compute consecutive LC state over a switching sequence.
The Ericksen-Leslie equations are solved given elastic constants,
viscosities, dielectric constants, and sequence of times/applied voltages/frequency
of applied voltages. Shear flow is included in dynamical calculations unless
otherwise specified. Optionally, the simulation can include
$n$ species of charged impurities, charges bound to either surface,
the flexoelectric effect, and a seqeunce of pressure gradients which
drive flow. A full description is available in S L Cornford 2008, PhD thesis, University
of Exeter.}
\usage{ lc_switching_sequence (
nz, times, voltages, frequencies, depth, k, alpha, epsilon,
 tilt_easy_angle, tilt_anchoring=c(0,0),
 twist_easy_angle=c(0,0),twist_anchoring=c(0,0),
 flexo=c(0,0), pgx=rep(0,length(times)),pgy=rep(0,length(times)),
 rho0=0,mu=rep(0,length(rho0)),
 chi=rep(0,length(rho0)), charge=rep(0,length(rho0)), sigma0=rep(0,length(rho0)*2),
 shear=TRUE, shear_y=shear, initial=NULL,
 ntol=1e-5,niter=12,ttol=2e-3,tmin=1e-6,ftime=1.0,iters=integer(nt))
}
\arguments{
\item{nz}{number of sublayers to consider}
\item{times}{vector of times (in seconds). If two consecutive times are equals,
a steady-state calculation is computed for the second.}
\item{voltages}{vector of RMS voltages at each time}
\item{frequencies}{vector of frequencies of the applied voltage at each time.}
\item{depth}{cell depth}
\item{k}{3 elastic constants}
\item{alpha}{6 Leslie viscosities}
\item{epsilon}{$\epsilon _\perp$ and $\delta \epsilon$}
\item{tilt_easy_angle}{tilt easy angle at $z=0$ and $z= d$}
\item{tilt_anchoring}{tilt anchoring strength. Defaults to 0 (infinite anchoring)}
\item{twist_easy_angle}{twist easy angle at $z=0$ and $z=d$. Defaults to 0}
\item{twist_anchoring}{twist anchoring strength. Defaults to 0 (infinite anchoring)}
\item{flexo}{flexo-electric coeffients $(e_s - e_b)$ and $(e_s + e_b)$. Defaults to 0}
\item{pgx}{vector of pressure gradients parallel to the x-axis at each time. Defaults to 0}
\item{pgy}{vector of pressure gradients parallel to the y-axis at each time. Defaults to 0}
\item{rho0}{vector of mean ion concentrations, ordered by species. Defaults to 0}
\item{mu}{vector of ion mobilities, ordered by species. Defaults to 0}
\item{chi}{vector of ion diffusivities, ordered by species. Defaults to 0}
\item{charge}{charge of ion species / e, , ordered by species. Defaults to 0}
\item{sigma0}{initial surface charge density * cell depth , ordered by boundary,species}
\item{shear}{consider shear flow? Defaults to TRUE}
\item{shear_y}{consider shear flow along y? Defaults to shear}
\item{initial}{allows an initial guesses to be set for lc profiles, held
  in a list of nz x nt matrices (tilt,twist,volt\_ac,volt\_dc,flowx,flowy)}
\item{ntol}{tolerance for Newton iterations. Defaults to 1e-5 }
\item{niter}{Maximium number of Newton iterations. Defaults to 12}
\item{ttol}{tolerance for time steps. Defaults to 2e-3}
\item{tmin}{minumum time step. Defaults to 1e-5}
\item{ftime}{weight for current profile in time integrals. default is 1,
  for implicit scheme. 0.5 gives Crank-Nicholson scheme}
\item{iters}{on return, filled with the numebr of newton iterations
  required for each solution}
}          
\seealso{}
\examples{

#define matrix plotting function for convenience
mp <- function(x,y,...){
matplot(x,y,type='l',lty=1,
col=colorRampPalette(c("blue","orange"),space="Lab")(ncol(y)))
}

#ZLI-4788
k <- c(13.7,10,20)*1e-12
alpha <- c(0.01,-0.25,0.05,0.2,0,0)
epsilon <- c(10.2,-4.7)

#3 micron HAN cell
d <- 3e-6
te <- c(0,pi/2)
nz <- 48

#ac static profiles
v <- c(0,1,3,5,7)
f <- rep(1e+5,length(v))
t <- rep(0,length(v))
lc <- lc_switching_sequence (nz, t, v, f, d, k, alpha, epsilon, te, shear=FALSE,flexo=c(5e-11,0))

mp(lc$z,lc$tilt,xlab="z",ylab=expression(theta(z,t)))
mp(lc$z,lc$volt_ac,xlab="z",ylab=expression(v(z,t)))

#ac 3v switch-on
n <- 8
v <- c(0,rep(3,n))
f <- rep(1e+5,n+1)
t <- c(0,2^(1:n - 1) * 1e-3)
lc <- lc_switching_sequence (nz, t, v, f, d, k, alpha, epsilon, te, flexo=c(5e-11,0))

mp(lc$z,lc$tilt,xlab="z",ylab=expression(theta(z,t)))
mp(lc$z,lc$flowx,xlab="z",ylab=expression(u(z,t)))


#dc 3v switch-on, with ions

kTq <- 290 * 1.3806503e-23 / 1.602176462e-19
nt <- 8

tp <- 2^(0:nt)*1e-3
t<-c(0,tp,tp + max(tp),tp+2*max(tp),tp+3*max(tp),tp*max(tp)+4*max(tp))
v<-c(rep(0,length(tp)+1),rep(-3,length(tp)),rep(3,length(tp)),rep(0,2*length(tp)))
f <- rep(0,length(v))
mu<-c(1,1)*1e-10

lcp <- lc_switching_sequence (nz, t, v, f, d, k, alpha, epsilon, te,flexo=c(0,0),
  shear=FALSE,rho0=c(25,30),sigma0=c(5,0,0,0),mu=mu,chi=kTq*mu,charge=c(-1,+1))

lcm <- lc_switching_sequence (nz, t, - v, f, d, k, alpha, epsilon, te,flexo=c(0,0),
  shear=FALSE,rho0=c(25,30),sigma0=c(5,0,0,0), mu=mu,chi=kTq*mu,charge=c(-1,+1))


plot(lcp$times,lcp$tilt[nz/2,],type='l',col="red")
lines(lcm$times,lcm$tilt[nz/2,],type='l',col="blue")

plot(lcm$times,lcm$rho[1,,1],type='l',col="red")
lines(lcp$times,lcp$rho[1,,1],type='l',col="blue")

plot(lcm$times,lcm$rho[1,,2],type='l',col="red")
lines(lcp$times,lcp$rho[1,,2],type='l',col="blue")





}
\keyword{}

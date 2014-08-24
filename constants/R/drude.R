##' @param lambda
##' @param omega
##' @param omega.p
##' @param gamma.p
##' @param basic
##' @param surf
##' @param interband
##' @param eps.inf ionbackground
##' @param nu.f
##' @param leff surface scattering
##' @return ...
##' @examples

##' library(constants)
##' data(AuJC)
##' 
##' gold.raw <- predict.material(AuJC, range=c(0.7, 1.2)) # original values for fit
##' gold <- predict.material(AuJC, range=c(0.3, 1.2), n=500) # interpolated values in wider range
##' 
##' model <- function(p, lambda=gold$wavelength){
##'   epsilon <- drude(lambda=lambda*1e-6, # in meters
##'                    interband=FALSE, basic=TRUE, surf=FALSE,
##'                    omega.p = p[1]*1e+16, 
##'                    gamma.p = p[2]*1e+14,
##'                    eps.inf=p[3])
##' 
##'   data.frame(wavelength=lambda, epsilon=epsilon)
##' }
##' 
##' 
##' 
##' objective <- function(p){
##'   
##' 	model <- model(p, lambda=gold.raw$wavelength)$epsilon
##' 	sos.real <- sum( (Re(model) - Re(gold.raw$epsilon))^2 ) / sum(Re(gold.raw$epsilon)^2)
##' 	sos.imag <- sum( (Im(model) - Im(gold.raw$epsilon))^2 ) / sum(Im(gold.raw$epsilon)^2)
##' 	sos <-  sos.real + sos.imag
##' 	sos
##' }
##' 
##' 
##' p0 <- c(1, 1, 1)
##' guess <- model(p0)
##' 
##' res <- optim(p0, objective)
##' fit <- model(res$par)
##' 
##' splitComplex <- function(d){
##' 
##'   transform(d, real=Re(epsilon), imag=Im(epsilon))
##' }
##' 
##' comparison <- lapply(list(raw=gold.raw,
##'                           data=gold,
##'                           guess=guess,
##'                           fit=fit), splitComplex)
##' \dontrun{
##' library(ggplot2)
##' 
##' m <- melt(comparison[-1], meas=c("real", "imag"))
##' 
##' str(m)
##' 
##' ggplot(m, aes(wavelength, value)) +
##'   facet_grid(variable~., scales="free_y")+
##'   geom_path(aes(colour=L1)) +
##'   geom_point(data=melt(comparison[1], meas=c("real", "imag"))) +
##'   ylab(expression(epsilon)) +
##'   xlab(expression(wavelength*" / "*mu*m)) +
##'   scale_colour_discrete("") +
##'   theme_bw()
##' }



`drude` <-
function (lambda=NULL, omega=NULL, omega.p = 1.3172e+16, gamma.p = 1.2991e+14, 
			basic=FALSE, surf = FALSE, interband = FALSE, 
			eps.inf = 1.54,
			A = 1, nu.f = 1.4 * 1e+06, leff = 1,
			A.i = c(1.27, 1.1),
			phi.i = c(-pi/4, -pi/4), 
			lambda.i = c(470, 325)*1e-9, 
			gamma.i = c(1900, 1060)*1e-9) 
{
	if(is.null(omega) & is.null(lambda)) stop("you need to provide omega or lambda!")
	if(!is.null(omega) & !is.null(lambda)) stop("you need to provide omega or lambda!")
	if(is.null(omega)) omega <- L2w(lambda)
	
gamma <- gamma.p
interbands <- 0.0

if (surf) {
        correction <- A * nu.f/leff
        gamma <- gamma.p + correction
}
if (basic){
	drude.C <- eps.inf -  omega.p^2 / (omega^2 +  (0+1i) * omega * gamma)
} else {
	drude.C <- 1 - omega.p^2 / (omega^2 +  (0+1i) * omega * gamma)
}

if (interband){
if(is.null(lambda)) lambda <- w2L(omega) # easier to work in wavelength for this
# gamma.l <- w2L(gamma.p) # 143*1e-9 
# lambda.p <- w2L(omega.p) # 14500*1e-9 
interbands <- rowSums(mapply(function(.A, .phi, .lambda, .gamma)	
	{
		.A / .lambda * (exp(1i*.phi) / (1/.lambda - 1/ lambda - 1i/.gamma) +
						exp(-1i*.phi) / (1/.lambda + 1/ lambda + 1i/.gamma) )
	}, .A = A.i, .phi = phi.i, .lambda = lambda.i, .gamma = gamma.i) )
}

    drude.C + interbands
}
##' @param p
##' @return ...

`drudeFit` <-
function (p) 
{
	if(length(p) == 2){
    drude(omega, omega.p = p[1], gamma.0 = p[2], basic=TRUE)
}
else{
    drude(omega, eps.inf = p[1], eps.0 = p[2], omega.p = p[3], gamma.0 = p[4])
	}
}

#include <complex.h>
#include <math.h>
#include <stdlib.h>
 #include <stdio.h>
void _lien_2x2_layer(complex *g, complex *s,double tilt, double twist, complex eo,
	   complex ee, double h, double kx, double k0)
{

  /* extended Jones matrices for a single layer */
  /* g - output : phase differences for ordinary,extraodrinary waves */
  /* s - output : amplitude division matrix */


  double ny = sin(tilt) * cos(twist);
  double nx = sin(tilt) * sin(twist);
  double nz = cos(tilt);
  
  complex de = ee-eo;
 
  //  complex exx = eo + de * nx * nx;
  complex exy = de * nx * ny;
  complex exz = de * nx * nz;
    
  complex eyy = eo + de * ny * ny;
  complex eyz = de * ny * nz;

  complex ezz = eo + (ee - eo) * nz * nz;

  complex kx2 = kx * kx;

  /*z-component of ordinary wave*/
  complex kzo2 = eo - kx2;
  complex kzo = csqrt(kzo2);
  /*z-component of extraordinary wave*/
  complex kze = - exz / ezz * kx +
    csqrt(eo*ee) /ezz * csqrt((ezz - (1.0 - de/ee * ny * ny )*kx2));
  complex kze2 = kze * kze;

  /* phase differences */
  complex kh = k0 * h;
  
  g[0] = cexp(I * kh * kzo);
  g[1] = g[2] = 0.0 + 0*I;
  g[3] = cexp(I * kh * kze);



  double eps = 1e-6;
  if (fabs(sin(tilt)) < eps){
    s[0] = 0; s[1] = 1; s[2] = 1; s[3] = 0; 
  } else {
    if (fabs(sin(twist)) < eps){
      s[0] = 0; s[1] = -1; s[2] = 1; s[3] = 0; 
    } else {
      s[0] = 1;
      s[3] = 1;
      s[2] = ((kx2 - ezz)*exy  + (kx * kzo + exz) * eyz)/
	((kzo2 + kx2 - eyy)*(kx2 - ezz) - eyz * eyz);

      s[1] = ((kx2 - ezz)*exy  + (kx * kze + exz) * eyz)/
	((kze2 + kx2 - eyy)*(kx2 - ezz) - eyz*eyz);
      s[1] = 1.0/s[1];
      
    }
  }
 
  
}

#define CMMEL(i,j) a[i]*b[j] + a[i+1]*b[j+2] 

void _lien_2x2_cmm(complex *a, complex *b, complex *r){
  /* product of two 2x2 matrices */
  r[0] = CMMEL(0,0);
  r[1] = CMMEL(0,1);
  r[2] = CMMEL(2,0);
  r[3] = CMMEL(2,1);
}

void _lien_2x2_cinv(complex *a, complex *r){
  /* inverse of 2x2 matrix*/
  complex d = a[0] * a[3] - a[1] * a[2];
  if (creal(d) != 0 || cimag(d) !=0){
    r[0] = a[3] / d;
    r[1] = -a[1] / d;
    r[2] = -a[2] / d;
    r[3] = a[0] / d;
  } else {

    r[0] = r[1] = r[2] = r[3] = 0;
    
  }
}

void _lien_2x2_stack(complex *s, int n, double *tilt, double *twist, 
		      complex *eo,complex *ee, double *h, double kx, double k0)
{
  /* extended jones matrix for stack */
  s[0] = 1;
  s[1] = 0;
  s[2] = 0;
  s[3] = 1; 

  complex u[4],v[4]; /* workspace matrices */
  complex ls[4],lg[4]; /* layer division and propagation matrices */

  int i;
  for (i = n-1; i >= 0 ; i--){
    _lien_2x2_layer(lg,ls,tilt[i],twist[i],eo[i],ee[i],h[i],kx,k0); 
    _lien_2x2_cinv(ls,u);
    _lien_2x2_cmm(u,s,v);
    _lien_2x2_cmm(lg,v,u);
    _lien_2x2_cmm(ls,u,s);
  }
}


#define CABS2(z) pow(creal(z),2) + pow(cimag(z),2)

void lien_2x2_optics(double *Tpp, double *Tps, double *Tsp, double *Tss,
		     complex *tpp, complex *tps, complex *tsp, complex *tss,
		     int *nangle, double *angle, double *azimuth, double *lambda,
		     int *nlayer, double *h, double *tilt, double *twist, 
		     complex *eo, complex *ee)
{

  /* compute transmission intensity and amplitude coefficients */
  /* for a multilayer anisotropic stack */
  /* assumes that eo[0] = eo[nlayer] and eo[0] is real */
  
  double k0 = 2 * M_PI / *lambda;
  double *twist_a;
  twist_a = malloc(sizeof(double)*(*nlayer));
  complex s[4];

  int i,j;
  for (i = 0; i < *nangle; i++){

    double kx = sin(angle[i]) * sqrt(creal(eo[0]));
    
    for (j = 0; j < *nlayer; j++){
      twist_a[j] = twist[j] + azimuth[i];
      }

    _lien_2x2_stack(s, *nlayer, tilt, twist_a, eo,ee, h, kx, k0);

    double ca = cos(angle[i]);
    tpp[i] = s[0];
    tps[i] = s[2] * ca;
    tsp[i] = s[1] / ca;
    tss[i] = s[3];

    Tpp[i] = CABS2(tpp[i]);
    Tps[i] = CABS2(tps[i]);
    Tsp[i] = CABS2(tsp[i]);
    Tss[i] = CABS2(tss[i]);

  }
  
  free(twist_a);

}



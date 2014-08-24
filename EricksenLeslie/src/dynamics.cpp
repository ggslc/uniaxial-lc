/*
  lcphysics1d: simulation of one dimensional liquid crystal cells
  copyright (C) 2007  Steph Cornford

  This program is free software; you can redistribute it and/or
  modify it under the terms of the GNU General Public License
  as published by the Free Software Foundation; either version 2
  of the License, or (at your option) any later version.

  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with this program; if not, write to the Free Software
  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.
*/


//dynamics.cpp : implements 1D Leslie-Ericksen
//     equations for splay-bend geometries
//#include <gsl/gsl_blas.h>
#include "dynamics.hh"
#include <math.h>
#include <iostream>

#define XMIN 2.225074e-308 

//#define UPWIND
//#define UPJACO
//#define EXPON

//#define POW

//#include <gsl/gsl_blas.h>

// inline functions to return values
// of tilt,volt,flowx at current and prior
// timesteps

inline 
double Power(double x, int n)
{
  return pow(x,n);
  //return x;
}

inline 
double flexo(State* s, int index){
  return s->parameters()->flexo[index];
}

inline 
double pgx(State* s){
  return s->parameters()->pgx;
}


inline 
double pgy(State* s){
  return s->parameters()->pgy;
}


inline 
double Sin(double x)
{
  return sin(x);
  //return x;
}

inline 
double Cos(double x)
{
  return cos(x);
}

inline double sgn(double x){
  if (x > 0){
    return 1.0;
  } else {
    return -1.0;
  }
  
}


inline 
double tilt(State* s, int i, int j){
  if ((j == 0) || !s->previous()) 
    return s->tilt[i];
  return tilt(s->previous(),i,j+1);
}

inline 
double twist(State* s, int i, int j){
  if ((j == 0) || !s->previous()) 
    return s->twist[i];
  return twist(s->previous(),i,j+1);
}

inline 
double volt_ac(State* s, int i, int j){
  if ((j == 0) || !s->previous()) {
    //double v =  s->volt_ac[i];
    return s->volt_ac[i];
    
  }
  return volt_ac(s->previous(),i,j+1);
}

inline 
double volt_dc(State* s, int i, int j){
  if ((j == 0) || !s->previous()) 
    return s->volt_dc[i];
  return volt_dc(s->previous(),i,j+1);
}





inline 
double flowx(State* s, int i, int j){
  if  ((j == 0) || !s->previous()) 
    return   s->flowx[i];
  return flowx(s->previous(),i,j+1);
}
inline 
double flowy(State* s, int i, int j){
  if  ((j == 0) || !s->previous()) 
    return  s->flowy[i];
  return flowy(s->previous(),i,j+1);
}

inline 
double rho_p(State* s, int i, int j){
  if  ((j == 0) || !s->previous()) 
    return s->rho_p[i];
  return rho_p(s->previous(),i,j+1);
}

inline 
double rho_m(State* s, int i, int j){
  if  ((j == 0) || !s->previous()) 
    return s->rho_m[i];
  return rho_m(s->previous(),i,j+1);
}


inline 
double rho(State* s, size_t sp, int i, int j){
  if  ((j == 0) || !s->previous()) 
    return s->rho[sp][i];
  return rho(s->previous(),sp,i,j+1);
}

inline double  sigma(State* s, size_t sp, int i,int j){
  if  ((j == 0) || !s->previous()) 
    return s->sigma[sp][i];
  return sigma(s->previous(),sp,i,j+1);
}


inline double surface_charge(State*s, int i, int j){
  if  ((j == 0) || !s->previous()) {
    double q = 0;
    for (size_t sp = 0; sp < s->sigma.size(); sp++){
      q += s->parameters()->charge[sp] * s->sigma[sp][i];
      
    }
    return q;
  }
  return surface_charge(s->previous(),i,j+1);
}


inline double charge(State* s, int i, int j){
  if  ((j == 0) || !s->previous()) {
    double q = 0;
    for (size_t sp = 0; sp < s->rho.size(); sp++)
      q += s->parameters()->charge[sp] * s->rho[sp][i];
    return q;
  }
  return charge(s->previous(),i,j+1);
}

inline 
double sign(State* s, size_t sp){
  return (s->parameters()->charge[sp] > 0)?1:-1;
}

inline 
double charge(State* s, size_t sp){
  return s->parameters()->charge[sp];
}


inline 
double mu(State* s, size_t sp){
  return s->parameters()->mu[sp];
}

inline 
double chi(State* s, size_t sp){
  return s->parameters()->chi[sp];
}

inline 
double ftime(State* s){
  return s->parameters()->ftime;
}

inline 
double alpha(State* s, int index){
  return s->parameters()->alpha[index];
}

inline 
double K(State* s, int index){
  return s->parameters()->K[index];
}

inline 
double w(State *s, int index){
  return s->parameters()->w[index]; 
}; 

//inline 
//double mu(State *s, int index){
// return s->parameters()->mu; 
//}; 


inline 
double ea(State* s){
  return  s->parameters()->ea;
}

inline 
double eo(State* s){
  return s->parameters()->eo;
}

inline 
double dt(State* s){
  return s->dt();
}

inline 
double dz(State* s){
  return s->dz();
}

inline
double power_law(double p){
  // 5th power law of patankar
  return max(0.0,pow(1-0.1*fabs(p),5));
}


inline double boundary_volt_ac(State* s, int boundary){
  return (s->parameters()->boundary_volt_ac[boundary]) ;
}

inline double boundary_volt_dc(State* s, int boundary){
  return s->parameters()->boundary_volt_dc[boundary] ;
}

inline double field_dc_centre(State* s, size_t i){
  //dc field at the centre of cell i (multipled by dz(s))

  if ((i < s->volt_dc.size() - 1) && i > 0) {
    return (- s->volt_dc[i+1] + s->volt_dc[i-1]) * .5
      // contribution to field from surface charges
      + s->parameters()->d * dz(s) * (surface_charge(s,0,0) - surface_charge(s,1,0)) 
      / (2.0 * E0 * (eo(s) + ea(s) * pow(cos(tilt(s,i,0)),2)));  
  }
  return 0;

}


inline double field_dc_edge(State* s, size_t i){
  //dc field at the east edge of cell i (multipled by dz(s))

  if (i < s->volt_dc.size() - 1) {
    return (- s->volt_dc[i+1] + s->volt_dc[i]) 
      // contribution to field from surface charges
      + s->parameters()->d * dz(s) * (surface_charge(s,0,0) - surface_charge(s,1,0)) 
      / (2.0 * E0 * (eo(s) + ea(s) * pow(0.5 * (cos(tilt(s,i,0)) + cos(tilt(s,i,0))),2)));  
  }
  return 0;

}





inline double boundary_tilt(State* s, int boundary){
  return s->parameters()->boundary_tilt[boundary];
}


inline double boundary_twist(State* s, int boundary){
  return s->parameters()->boundary_twist[boundary];
}

inline double boundary_flowx(State* s, int boundary){
  return s->parameters()->boundary_flowx[boundary];
}

inline double boundary_flowy(State* s, int boundary){
  return s->parameters()->boundary_flowy[boundary];
}
inline double dissipation_m(State*s, int i){
  return alpha(s,3) * pow(sin (tilt(s,i,0)),2) 
    - alpha(s,2) * pow(cos(tilt(s,i,0)),2);
}


inline double vdamping(State *s){
  //under-relaxation factor, for voltage equations
  
   double f = fabs(flexo(s,1)/ ea(s));
   return (f < 1.0)?1.0:exp(-0.5 * (f-1.0));
}

// dc electric potential equation (includes terms from flexoelectricity, and ion terms)


#define TINY 1e-8

double Volt_DCEquation :: difference(State *a, State *b)
{
  double diff = 0;
  double sum = 0;
  
  //double d = vdamping(s);

   for (size_t k = 0; k  < b->n_lc_layers(); k++) {
     diff += pow(volt_dc(b,k,0) - volt_dc(a,k,0) ,2);
     sum +=  pow(volt_dc(b,k,0) + volt_dc(a,k,0) ,2);
    
   }	
   if (sum < TINY) sum = 1.0;

   return diff/sum;

}


void Volt_DCEquation :: update(VEC dv, State*s)
{
 
  // 
  double damping = vdamping(s);
  
  for (size_t i = 0; i < s->n_lc_layers(); i++){
    s->volt_dc[i] -= damping * (*dv)[i];
  }
  
  
}

void Volt_DCEquation :: evaluate(VEC rhs, VEC ju, VEC jl, VEC jd, State* s){

  int n = s->n_lc_layers();
  double fluxp = 0, fluxm = 0;
  for (int i = 1; i < n-1; i++){

    //residual
    if (i > 1) {
      fluxm = fluxp;
    }
    else {
      fluxm = (
	       (flexo(s,1)*(Sin(2*tilt(s,-1 + i,0)) + Sin(2*tilt(s,i,0)))*(tilt(s,-1 + i,0) - tilt(s,i,0))) +       
	       Power(Cos(tilt(s,-1 + i,0)) + Cos(tilt(s,i,0)),2)*ea(s)*(volt_dc(s,-1 + i,0) - volt_dc(s,i,0)) +       
	       4*eo(s)*(   volt_dc(s,-1 + i,0) - volt_dc(s,i,0)))/4.
	;
    }



    fluxp = (
	     // sign error ....- 
	     (flexo(s,1)*(Sin(2*tilt(s,i,0)) + Sin(2*tilt(s,1 + i,0)))*(+ tilt(s,i,0) - tilt(s,1 + i,0))) + 	      
	     Power(Cos(tilt(s,i,0)) + Cos(tilt(s,1 + i,0)),2)*ea(s)*(volt_dc(s,i,0) - volt_dc(s,1 + i,0)) +	      
	     4*eo(s)*(volt_dc(s,i,0) - volt_dc(s,1 + i,0)))/4.
      ;



   
    
    
    double src =  -1.0/E0 * charge(s,i,0)*Power(dz(s),2) + fluxp - fluxm;
   
   

    (*rhs)[i] = src;  
    
    (*jd)[i] =   (Power(Cos(tilt(s,-1 + i,0)) + Cos(tilt(s,i,0)),2)*ea(s) 
	       + Power(Cos(tilt(s,i,0)) + Cos(tilt(s,1 + i,0)),2)*ea(s) + 8*eo(s))/4.   
      ;
    

    if (i > 0) (*ju)[i-1] = 
		 (-(Power(Cos(tilt(s,-1 + i,0)) + Cos(tilt(s,i,0)),2)*ea(s)) 
		  - 4*eo(s))/4.
			      ;
				   


    if (i < n-1) (*jl)[i] =
		   (-(Power(Cos(tilt(s,i,0)) + Cos(tilt(s,1 + i,0)),2)*ea(s)) - 4*eo(s))/4.
		   ;

  }

  // boundary conditions
 //  double vsig = 0.0;
//   for (int i = 0; i < n; i++)
//     {
//       vsig += s->parameters()->d * dz(s) * (surface_charge(s,0,0) - surface_charge(s,1,0)) 
// 	/ (2.0 * E0 * (eo(s) + ea(s) * pow(cos(tilt(s,i,0)),2)));
//     }



  //residual
  (*rhs)[0] = BIG*(volt_dc(s,0,0) - boundary_volt_dc(s,0));
  (*rhs)[n-1] = BIG*(volt_dc(s,n-1,0) - boundary_volt_dc(s,1));

  //diagonal
  (*jd)[0] = BIG;
  (*jd)[n-1] = BIG;
  
//off diagonal
  (*ju)[n-2] = 0;
  (*jl)[0] = 0;
  
}


// ac electric potential equation (no flexoelectric term)

void Volt_ACEquation :: update(VEC dv, State*s)
{

  //double damping = vdamping(s);
     for (size_t i = 0; i < s->n_lc_layers(); i++){
       s->volt_ac[i] -=  (*dv)[i];
}
}

double Volt_ACEquation :: difference(State *a, State *b)
{
  double diff = 0;
  double sum = 0;

   for (size_t k = 0; k  < b->n_lc_layers(); k++) {
     diff += pow(volt_ac(b,k,0) - volt_ac(a,k,0) ,2);
     sum +=  pow(volt_ac(b,k,0) + volt_ac(a,k,0) ,2);
    
   }	
   if (sum < TINY) sum = 1.0;

   return diff/sum;

}


void Volt_ACEquation :: evaluate(VEC rhs, VEC ju, VEC jl, VEC jd, State* s){

  int n = s->n_lc_layers();


  double fluxm = 0, fluxp = 0;

  for (int i = 1; i < n-1; i++){

    //residual

      if (i > 1) {
       fluxm = fluxp;
     }
     else {
    fluxm = (
	     Power(Cos(tilt(s,-1 + i,0)) + Cos(tilt(s,i,0)),2)*ea(s)*(volt_ac(s,-1 + i,0) - volt_ac(s,i,0)) + 
	     4*eo(s)*(volt_ac(s,-1 + i,0) - volt_ac(s,i,0)))/4.
      ;
     }

    fluxp = (
	     Power(Cos(tilt(s,i,0)) + Cos(tilt(s,1 + i,0)),2)*ea(s)*(volt_ac(s,i,0) - volt_ac(s,1 + i,0)) + 
	     4*eo(s)*(volt_ac(s,i,0) - volt_ac(s,1 + i,0)))/4.
      ;

    double src =  fluxp - fluxm;
   
    (*rhs)[i] = src;  
    
   
    (*jd)[i] = (Power(Cos(tilt(s,-1 + i,0)) + Cos(tilt(s,i,0)),2)*ea(s) 
	    + Power(Cos(tilt(s,i,0)) + Cos(tilt(s,1 + i,0)),2)*ea(s) + 8*eo(s))/4.
      ;


    if (i > 0) (*ju)[i-1] =
		 (-(Power(Cos(tilt(s,-1 + i,0)) + Cos(tilt(s,i,0)),2)*ea(s)) - 4*eo(s))/4.
		 ;
				   

    
    if (i < n-1) (*jl)[i] = 
		   (-(Power(Cos(tilt(s,i,0)) + Cos(tilt(s,1 + i,0)),2)*ea(s)) - 4*eo(s))/4.
		   ;
    

  }

  // boundary conditions

  //residual
  (*rhs)[0] = BIG*(volt_ac(s,0,0) - boundary_volt_ac(s,0));
  (*rhs)[n-1] = BIG*(volt_ac(s,n-1,0) - boundary_volt_ac(s,1));
  
  //diagonal
  (*jd)[0] = BIG;
  (*jd)[n-1] = BIG;

  //off diagonal
  (*ju)[n-2] = 0;
  (*jl)[0] = 0;
  
}





// ion transport




void IonEquation :: update(VEC drho, State*s)
{
  

  for (size_t i = 0; i < s->n_lc_layers(); i++){
    s->rho[species][i]  -= 1.0 * (*drho)[i];
    s->rho[species][i] = max(XMIN,s->rho[species][i]);
  }
  
}

double IonEquation :: difference(State *a, State *b)
{
  double diff = 0;
  double sum = 0;

   for (size_t k = 0; k  < b->n_lc_layers(); k++) {
     diff += pow(rho(b,species,k,0) - rho(a,species,k,0) ,2);
     sum +=  pow(rho(b,species,k,0) + rho(a,species,k,0) ,2);
    
   }	
   if (sum < TINY) sum = 1.0;

   return diff/sum;

}


void IonEquation :: evaluate(VEC rhs, VEC jl, VEC ju, VEC jd, State* s){
  
  int n = s->n_lc_layers();
  double fluxp = 0, fluxm = 0, ep = 0, em = 0;
  double Pe = 0 , Pw = 0;


  // equation and  boundary condition are evaluated together



  for (int i = 0; i < n ; i++){
    fluxm = fluxp;


    em = ep;
    Pw = Pe;

    if (i < n-1) {
      // convection coefficient
      ep =  mu(s,species) * charge(s,species)/fabs(charge(s,species)) * (field_dc_edge(s,i));
 
#ifdef UPWIND
      //upwind scheme   
      fluxp = -  chi(s,species) * (rho(s,species,1 + i,0) - rho(s,species,i,0));
      
      if (ep > 0){
      	fluxp += ep * rho(s,species,i,0);
      } else {
      	fluxp += ep * rho(s,species,i+1,0);
      }



      //exponential scheme
#else
#define PECLET 1e-16
      Pe = ep /  chi(s,species);
      if ( fabs(Pe) > PECLET){
	//cout << Pe << endl;
	
	//	fluxp = ep * (exp(Pe) * (rho(s,species,i,0) - rho(s,species,i+1,0))/ exp(2*Pe - 1.0)
	//	      + (exp(Pe) - 1.0) * (exp(Pe) * rho(s,species,i,0) + rho(s,species,i+1,0))/
	//	      (exp(2 * Pe) - 1.0));

	fluxp = ep * (exp(Pe) * rho(s,species,i,0) -  rho(s,species,i+1,0))/(exp(Pe) - 1.);

      } else {
	
	fluxp = chi(s,species) * (rho(s,species,i,0) - rho(s,species,i+1,0));
	  //+ 0.5 * ep *  (rho(s,species,i,0) + rho(s,species,i+1,0));

      }
      
#endif

    } else {
      fluxp = 0;
      ep = 0;
    }

    //residual

    double src = (fluxp - fluxm);

    if (s->time_dependent()){
      src +=  (Power(dz(s),2)*(-rho(s,species,i,-1) + rho(s,species,i,0)))/dt(s);
    } 

   

    (*rhs)[i] = src ;
    

    bool cond = true;//(fabs(boundary_volt_dc(s,0) - boundary_volt_dc(s,1)) > 0.1);
    
    double d = 0; 
    double du = 0;
    double dl = 0;

 

    if ( i > 0){

#ifdef UPWIND

      //upwind upwind
      d += chi(s,species);
      dl += -chi(s,species);
#ifdef UPJACO  
      //upwind term
      if (cond){
	if (em > 0)
	  dl +=  - em ;
	else 
	  d  += -em ;
      }
#endif

#else
      //exponential
      if (cond && 
	  (fabs(Pw) > PECLET)){
	 dl +=  - em * (exp(Pw))/(exp(Pw) - 1.);
	 d +=    em / (exp(Pw) - 1.);  
	   
	   
       } else {
	 d +=  chi(s,species);// + 0.5 * em;
	   dl += - chi(s,species);// + 0.5 * em;
       }

#endif

    }
     if (i < n-1){

#ifdef UPWIND

      //upwind
      d += chi(s,species);
      du += -chi(s,species);
     
#ifdef UPJACO
      //upwind term
      if (cond){
	if (ep > 0)
	  d += ep ;
	else 
	  du += ep ;
      }
#endif
#else
             //exponential
      if (cond && 
	   (fabs(Pe) > PECLET)){
	d +=  ep * (exp(Pe))/(exp(Pe) - 1.);
	du += - ep / (exp(Pe) - 1.);
      } else {
	d += chi(s,species);// + 0.5 * ep;
	du += -  chi(s,species);// + 0.5 * ep;
      }
#endif
      
     }
     
     
     //}
     
     
     if (s->time_dependent()){
       d +=    pow(dz(s),2)/dt(s);

     } else {
      
       
	 d += 1e-12;
       
     }
    

    //diagonal of jacobian
    (*jd)[i] = d ;
    //subdiagonal df(i)/dx(i-1)
    if (i > 0) (*jl)[i-1] = dl;
    //superdiagonal  df(i)/dx(i+1)
    if (i < n-1) (*ju)[i] = du;
  }
  
}





// flowx equation
void FlowxEquation :: update(VEC d, State*s)
{
  for (size_t i = 0; i < s->n_lc_layers(); i++){
    s->flowx[i] -= (*d)[i];
  }
}


double FlowxEquation :: difference(State *a, State *b)
{
  double diff = 0;
  double sum = 0;

   for (size_t k = 0; k  < b->n_lc_layers(); k++) {
     diff += pow(flowx(b,k,0) - flowx(a,k,0) ,2);
     sum +=  pow(flowx(b,k,0) + flowx(a,k,0) ,2);
    
   }	
   if (sum < TINY) sum = 1.0;

   return diff/sum;

}


double FlowxEquation :: constant(int i, State* s) {
  
  return 0.;
    
}

inline double g(State* s, double theta){
  return 0.5 * (alpha(s,4) 
		+(alpha(s,5) - alpha(s,2)) * pow(cos(theta),2)
		+(alpha(s,3) + alpha(s,6)) * pow(sin(theta),2)
		+ 2.0 * alpha(s,1) * pow(sin(theta)*cos(theta),2));
		
}



inline double m(State* s, double theta){
  return alpha(s,3) * pow (sin(theta),2) 
    - alpha(s,2) * pow (cos(theta),2); 
}


inline double gxux(State *s, double theta, double phi){
  //coefficient of dux/dz in x-momentum equation
  return 0.5 * ((2.0 * alpha(s,1) *pow(cos(theta),2) + alpha(s,6) + alpha(s,3)) 
		* pow(cos(phi),2) * pow(sin(theta),2)
		+ (alpha(s,5) - alpha(s,2)) * pow(cos(theta),2) + alpha(s,4));
	  
}

inline double gxuy(State *s, double theta, double phi){
  //coefficient of duy/dz in x-momentum equation
    return 0.5 * ((2.0 * alpha(s,1) *pow(cos(theta),2) + alpha(s,6) + alpha(s,3)) 
  	* 0.5 * sin(2.0 * phi) * pow(sin(theta),2));

  //return 0.5 * alpha(s,3) * pow(sin(theta),2) * sin(phi)*cos(phi);
  

}

inline double gxp(State *s, double theta, double phi){
  //coefficient of dphi/dt in x-momentum equation
  return - alpha(s,2) * sin(phi) * 0.5 * sin(2.0*theta);

}

inline double gxt(State *s, double theta, double phi){
  //coefficient dtheta/t in x-momentum equation
  return (- alpha(s,3) * pow (sin(theta),2) 
	  + alpha(s,2) * pow (cos(theta),2)) 
   * cos(phi);
}


void FlowxEquation :: evaluate (VEC rhs, VEC ju, VEC jl, VEC jd, State* s){
  // navier-stokes equation for x component of flow

  int n = s->n_lc_layers();
  
  
  for (int i = 1; i < n-1; i++){
    
    double tilte = 0.5 * (tilt(s,i,0) + tilt(s,1 + i,0));
    double tiltw = 0.5 * (tilt(s,i-1,0) + tilt(s,i,0));
    double twiste = 0.5 * (twist(s,i,0));// + twist(s,1 + i,0));
  double twistw = 0.5 * (twist(s,i-1,0));// + twist(s,i,0));
    
    double guxe = gxux(s,tilte,twiste);
    double guxw = gxux(s,tiltw,twistw);
    
    double src =  - pow(dz(s),2) * pgx(s) +
      (flowx(s,i-1,0) - flowx(s,i,0))*guxw - (flowx(s,i,0) - flowx(s,i+1,0))*guxe;
    
    //terms in duy/dz
    src += gxuy(s,tiltw,twistw) * (flowy(s,i-1,0) - flowy(s,i,0))
      - gxuy(s,tilte,twiste) * (flowy(s,i,0) - flowy(s,i+1,0));
   

    if (s->time_dependent()) {
      //time-dependent terms

      double gxte = gxt(s,tilt(s,i+1,0),twist(s,i+1,0));
      double gxtw = gxt(s,tilt(s,i-1,0),twist(s,i-1,0));

      src +=  
	 dz(s) * 0.5 / dt(s) * ( - gxte *  (tilt(s,i+1,-1) - tilt(s,i+1,0))
				+  gxtw *  (tilt(s,i-1,-1) - tilt(s,i-1,0))
				-  gxp(s,tilt(s,i+1,0),twist(s,i+1,0)) 
				* (twist(s,i+1,-1) - twist(s,i+1,0))
				+  gxp(s,tilt(s,i-1,0),twist(s,i-1,0)) 
				* (twist(s,i-1,-1) - twist(s,i-1,0)));
     
    }

    (*rhs)[i] = src;
    
    //diagonal of jacobian
    (*jd)[i] =  - guxe - guxw;
    
    //subdiagonal df(i)/dx(i-1)
    if (i > 0) (*ju)[i-1] =  guxw;
    
    //superdiagonal  df(i)/dx(i+1)
    if (i < n-1) (*jl)[i] =  guxe;
    
  }
  
  // boundary conditions
  
  //residual
  (*rhs)[0]  = BIG*(flowx(s,0,0) - boundary_flowx(s,0));
  (*rhs)[n-1] = BIG*(flowx(s,n-1,0) - boundary_flowx(s,1));
  
  //diagonal
  (*jd)[0] = BIG;
  (*jd)[n-1] = BIG;
  
}

// flowx equation
void FlowxEquationXZ :: update(VEC d, State*s)
{
  for (size_t i = 0; i < s->n_lc_layers(); i++){
    s->flowx[i] -= (*d)[i];
  }
}


double FlowxEquationXZ :: difference(State *a, State *b)
{
  double diff = 0;
  double sum = 0;

   for (size_t k = 0; k  < b->n_lc_layers(); k++) {
     diff += pow(flowx(b,k,0) - flowx(a,k,0) ,2);
     sum +=  pow(flowx(b,k,0) + flowx(a,k,0) ,2);
    
   }	
   if (sum < TINY) sum = 1.0;

   return diff/sum;

}


double FlowxEquationXZ :: constant(int i, State* s) {
  
  return 0.;
    
}



inline double gxuxxz(State *s, double theta){
  //coefficient of dux/dz in x-momentum equation
  return 0.5 * ((2.0 * alpha(s,1) *pow(cos(theta),2) + alpha(s,6) + alpha(s,3)) 
		* pow(sin(theta),2)
		+ (alpha(s,5) - alpha(s,2)) * pow(cos(theta),2) + alpha(s,4));
}


inline double gxtxz(State *s, double theta){
  //coefficient dtheta/t in x-momentum equation
  return (- alpha(s,3) * pow (sin(theta),2) 
	  + alpha(s,2) * pow (cos(theta),2));
}


void FlowxEquationXZ :: evaluate (VEC rhs, VEC ju, VEC jl, VEC jd, State* s){
  // navier-stokes equation for x component of flow

  int n = s->n_lc_layers();
  
  
  for (int i = 1; i < n-1; i++){
    
    double tilte = 0.5 * (tilt(s,i,0) + tilt(s,1 + i,0));
    double tiltw = 0.5 * (tilt(s,i-1,0) + tilt(s,i,0));
    
    double guxe = gxuxxz(s,tilte);
    double guxw = gxuxxz(s,tiltw);
    
    double src =  - pow(dz(s),2) * pgx(s) +
      (flowx(s,i-1,0) - flowx(s,i,0))*guxw 
      - (flowx(s,i,0) - flowx(s,i+1,0))*guxe;
    
    

    if (s->time_dependent()) {
      //time-dependent terms

      double gxte = gxtxz(s,tilt(s,i+1,0));
      double gxtw = gxtxz(s,tilt(s,i-1,0));

      src +=  
	 dz(s) * 0.5 / dt(s) * ( - gxte *  (tilt(s,i+1,-1) - tilt(s,i+1,0))
				+  gxtw *  (tilt(s,i-1,-1) - tilt(s,i-1,0))
				 );
     
    }

    (*rhs)[i] = src;
    
    //diagonal of jacobian
    (*jd)[i] =  - guxe - guxw;
    
    //subdiagonal df(i)/dx(i-1)
    if (i > 0) (*ju)[i-1] =  guxw;
    
    //superdiagonal  df(i)/dx(i+1)
    if (i < n-1) (*jl)[i] =  guxe;
    
  }
  
  // boundary conditions
  
  //residual
  (*rhs)[0]  = BIG*(flowx(s,0,0) - boundary_flowx(s,0));
  (*rhs)[n-1] = BIG*(flowx(s,n-1,0) - boundary_flowx(s,1));
  
  //diagonal
  (*jd)[0] = BIG;
  (*jd)[n-1] = BIG;
  
}



// flowy equation
void FlowyEquation :: update(VEC d, State*s)
{
  for (size_t i = 0; i < s->n_lc_layers(); i++){
    s->flowy[i] -= (*d)[i];
  }
}


double FlowyEquation :: difference(State *a, State *b)
{
  double diff = 0;
  double sum = 0;

   for (size_t k = 0; k  < b->n_lc_layers(); k++) {
     diff += pow(flowy(b,k,0) - flowy(a,k,0) ,2);
     sum +=  pow(flowy(b,k,0) + flowy(a,k,0) ,2);
    
   }	
   if (sum < TINY) sum = 1.0;

   return diff/sum;

}


inline double gyuy(State *s, double theta, double phi){
  //coefficient of duy/dz in y-momentum equation
  return 0.5 * ((2.0 * alpha(s,1) *pow(cos(theta),2) + alpha(s,6) + alpha(s,3)) 
		* pow(sin(phi),2) * pow(sin(theta),2)
		+ (alpha(s,5) - alpha(s,2)) * pow(cos(theta),2) + alpha(s,4));
	  
}

inline double gyux(State *s, double theta, double phi){
  //coefficient of dux/dz in y-momentum equation
   return 0.5 * ((2.0 * alpha(s,1) *pow(cos(theta),2) + alpha(s,6) + alpha(s,3)) 
  		* 0.5 * sin(2.0 * phi) * pow(sin(theta),2));  

   //return 0.5 * alpha(s,3) * pow(sin(theta),2) * sin(phi)*cos(phi);

}

inline double gyp(State *s, double theta, double phi){
  //coefficient of dphi/dt in y-momentum equation
  return + alpha(s,2) * cos(phi) * 0.5 * sin(2.0*theta);

}

inline double gyt(State *s, double theta, double phi){
  //coefficient dtheta/t in y-momentum equation
  return (- alpha(s,3) * pow (sin(theta),2) 
	  + alpha(s,2) * pow (cos(theta),2)) 
    * sin(phi);
}



void FlowyEquation :: evaluate (VEC rhs, VEC ju, VEC jl, VEC jd, State* s){
  // navier-stokes equation for y component of flow

 int n = s->n_lc_layers();
  
  
  for (int i = 1; i < n-1; i++){
    
    double tilte = 0.5 * (tilt(s,i,0) + tilt(s,1 + i,0));
    double tiltw = 0.5 * (tilt(s,i-1,0) + tilt(s,i,0));
    double twiste = 0.5 * (twist(s,i,0)+ twist(s,1 + i,0));
    double twistw = 0.5 * (twist(s,i-1,0) + twist(s,i,0));
    
    double guye = gyuy(s,tilte,twiste);
    double guyw = gyuy(s,tiltw,twistw);
    
    double src =  - pow(dz(s),2) *  pgy(s) + 
      (flowy(s,i-1,0) - flowy(s,i,0))*guyw - (flowy(s,i,0) - flowy(s,i+1,0))*guye;
    
    //terms in dux/dz
    src += 1.0*( gyux(s,tiltw,twistw) * (flowx(s,i-1,0) - flowx(s,i,0))
		  - gyux(s,tilte,twiste) * (flowx(s,i,0) - flowx(s,i+1,0)));
   

    if (s->time_dependent()) {
      //time-dependent terms
      double gyte = gyt(s,tilt(s,i+1,0),twist(s,i+1,0));
      double gytw = gyt(s,tilt(s,i-1,0),twist(s,i-1,0));

      src +=  
	dz(s) * 0.5 / dt(s) * ( - gyte * (tilt(s,i+1,-1) - tilt(s,i+1,0))
				+  gytw * (tilt(s,i-1,-1) - tilt(s,i-1,0))
				- gyp(s,tilt(s,i+1,0),twist(s,i+1,0)) 
				* (twist(s,i+1,-1) - twist(s,i+1,0))
				+  gyp(s,tilt(s,i-1,0),twist(s,i-1,0)) 
				* (twist(s,i-1,-1) - twist(s,i-1,0)));
      
    }

    (*rhs)[i] = src;
    
    //diagonal of jacobian
    (*jd)[i] =  - guye - guyw;
    
    //subdiagonal df(i)/dx(i-1)
    if (i > 0) (*ju)[i-1] =  guyw;
    
    //superdiagonal  df(i)/dx(i+1)
    if (i < n-1) (*jl)[i] =  guye;
    
  }
  

  // boundary conditions

  //residual
  (*rhs)[0]  = BIG*(flowy(s,0,0) - boundary_flowy(s,0));
  (*rhs)[n-1] = BIG*(flowy(s,n-1,0) - boundary_flowy(s,1));
  
  //diagonal
  (*jd)[0] = BIG;
  (*jd)[n-1] = BIG;

}



//twist equation

void TwistEquation :: update (VEC dtwist, State*s)
{
  for (size_t i = 0; i < s->n_lc_layers(); i++){
        s->twist[i] -= (*dtwist)[i];
   }
}


double TwistEquation :: difference(State *a, State *b)
{
  double diff = 0;
  double sum = 0;

   for (size_t k = 0; k  < b->n_lc_layers(); k++) {
     diff += pow((twist(b,k,0) - twist(a,k,0)) 
		 * 0.5 * (sin(tilt(b,k,0)) + sin(tilt(a,k,0))) ,2);
     sum +=  pow(twist(b,k,0) + twist(a,k,0) ,2);
    
   }	
   if (sum < TINY) sum = 1.0;

   return diff/sum;

}

inline double gtw(State *s, double theta){
  double s2 = pow(sin(theta),2);
  return K(s,3)/K(s,1) * s2 * ((K(s,2)/K(s,3) - 1) * s2 + 1); 
}
inline double dgtw(State *s, double theta){
  double ss = sin(theta);
  double cc = cos(theta);
  return 2 * K(s,3)/K(s,1) *  ss * cc * (2 * ss * ss * (K(s,2)/K(s,3) - 1) + 1);
}


void TwistEquation :: evaluate (VEC rhs, VEC ju, VEC jl, VEC jd, State* s){
  // conservation of angular momentum: twist component 

  int n = s->n_lc_layers();


  double a2 = alpha(s,2)/K(s,1);
  //double a3 = alpha(s,3)/K(s,1);
  //double gamma = a3 - a2;
  double gammap = alpha(s,3) - alpha(s,2);

  for (int i = 1; i < n-1; i++){

    double g =  gtw(s,tilt(s,i,0));
    double dg = dgtw(s,tilt(s,i,0));
    double e = 1e-9;
 
    g = max(g,pow(e,2));

    double q =  0.25 * dg/g * (tilt(s,i+1,0) - tilt(s,i-1,0));
    double dux = dz(s) * 0.5 * (flowx(s,i+1,0) - flowx(s,i-1,0));
    double duy = dz(s) * 0.5 * (flowy(s,i+1,0) - flowy(s,i-1,0));
    double fc =  0.5 *  sin(2. * tilt(s,i,0)) * a2 / g;
    double f =   fc *  (duy * cos(twist(s,i,0)) - dux * sin(twist(s,i,0)));
    
    //if (abs(sin(tilt(s,i,0))) < 1e-7) f = 0.;

	
    (*rhs)[i]  = q * (twist(s,i+1,0) - twist(s,i-1,0))
      +  (twist(s,i+1,0) + twist(s,i-1,0) - 2 * twist(s,i,0)) - f;
	
    double twe = twist(s,i+1,0);
    double tww = twist(s,i-1,0);
    double twp = twist(s,i,0);

    //double aemw = asin(sin(twe)*cos(tww) - cos(twe)*sin(tww));
    //double aemp = asin(sin(twe)*cos(twp) - cos(twe)*sin(twp));
    //double awmp = asin(sin(tww)*cos(twp) - cos(tww)*sin(twp));


    //double aemw = asin(sin(twe-tww));
    double aemp = asin(sin(twe-twp));
    double awmp = asin(sin(tww-twp));
    
    
    //"gets stuck" formula (slow)
    (*rhs)[i]  = q * (twe - tww) +  (aemp + awmp) - f;
    //(*rhs)[i]  = q * aemw +  (aemp + awmp) - f;
    //"no over-tilt" formula (quick)
    //(*rhs)[i]  = q * aemw + (twe + tww - 2 * twp)  - f;

    //"no over-tilt" formula (quick)
    //(*rhs)[i]  = q * (twe - tww) +  (twe + tww - 2 * twp) - f;


    (*jd)[i] =  - 2. 
      + min(0.0,fc * (- duy * sin(twist(s,i,0)) - dux * cos(twist(s,i,0)))) ;
    (*ju)[i-1] =  1. - q ;
    (*jl)[i] =  1. + q ; 


    if (s->time_dependent()){
      //max(pow(sin(tilt(s,i,0)),2),pow(e,2))
      double uf  = pow(dz(s),2)/dt(s) 
	* gammap / ( pow(sin(tilt(s,i,0)),2) * (K(s,2) - K(s,3)) + K(s,3));
	//* max(0.*pow(e,2),pow(sin(tilt(s,i,0)),2))  * gamma/g;
      (*rhs)[i] += uf * (twist(s,i,-1) - twist(s,i,0));
      (*jd)[i]  -= uf;
    }

    

  }

  //boundary conditions (just strong anchoring for now)

  //residual
  //if (twist(s,n-1,0) < twist(s,n-1,0) ){
    (*rhs)[0] = BIG*(twist(s,0,0) - boundary_twist(s,0));
  //} else {
    // (*rhs)[0] = BIG*(twist(s,0,0) - (boundary_twist(s,0) + 2. * M_PI) );
    //}
    //(*rhs)[0] = BIG*(twist(s,0,0) - 1.99*M_PI);

  (*rhs)[n-1] = BIG*(twist(s,n-1,0) - boundary_twist(s,1));

  //diagonal
  (*jd)[0] = BIG;
  (*jd)[n-1] = BIG;
  
//off diagonal
  (*ju)[n-2] = 0;
  (*jl)[0] = 0;

}





// tilt equation

void TiltEquation :: update(VEC dtilt, State*s)
{
  for (size_t i = 0; i < s->n_lc_layers(); i++){
    s->tilt[i] -=  (*dtilt)[i];
  }
}

 double TiltEquation :: difference(State *a, State *b)
{
  double diff = 0;
  double sum = 0;

   for (size_t k = 0; k  < b->n_lc_layers(); k++) {
     diff += pow(tilt(b,k,0) - tilt(a,k,0) ,2);
     sum +=  pow(tilt(b,k,0) + tilt(a,k,0) ,2);
    
   }	
   if (sum < TINY) sum = 1.0;

   return diff/sum;

}

double TiltEquation::steadyrhs(State *s, size_t i,int j){
  
  double kk = K(s,3)/K(s,1);
  double ee = E0/K(s,1);
  double ff =  (E0 * flexo(s,1))/K(s,1);
  double a2 = alpha(s,2)/K(s,1);
  double a3 = alpha(s,3)/K(s,1);
  //double gamma = a3 - a2;

 double r =   
   (cos(twist(s,i,j)) 
    * (4*a2*Power(Cos(tilt(s,i,j)),2)*dz(s)*(flowx(s,-1 + i,j) - flowx(s,1 + i,j)) + 
       4*a3*dz(s)*(-flowx(s,-1 + i,j) + flowx(s,1 + i,j))*Power(Sin(tilt(s,i,j)),2)) 
    +  sin(twist(s,i,j)) 
   *(4*a2*Power(Cos(tilt(s,i,j)),2)*dz(s)*(flowy(s,-1 + i,j) - flowy(s,1 + i,j)) + 
     4*a3*dz(s)*(-flowy(s,-1 + i,j) + flowy(s,1 + i,j))*Power(Sin(tilt(s,i,j)),2)))
   +
   Sin(2*tilt(s,i,j))*Power(tilt(s,-1 + i,j) - tilt(s,1 + i,j),2) - 
   kk*Sin(2*tilt(s,i,j))*Power(tilt(s,-1 + i,j) - tilt(s,1 + i,j),2) + 
   8*Power(Cos(tilt(s,i,j)),2)*kk*(tilt(s,-1 + i,j) - 2*tilt(s,i,j) + tilt(s,1 + i,j)) + 
   8*Power(Sin(tilt(s,i,j)),2)*(tilt(s,-1 + i,j) - 2*tilt(s,i,j) + tilt(s,1 + i,j)) 
   -ee*ea(s)*Sin(2*tilt(s,i,j))*Power(volt_ac(s,-1 + i,j) - volt_ac(s,1 + i,j),2)
   -ee*ea(s)*Sin(2*tilt(s,i,j))*Power(2.*field_dc_centre(s,i),2)
   ;

    //flexo terms
    r += - 4* ff *Sin(2*tilt(s,i,0)) 
      * (-field_dc_edge(s,i)  + field_dc_edge(s,i-1)) ;
   
    //twist terms;
    double twe = twist(s,i+1,0);
    double tww = twist(s,i-1,0);
    double aemw = asin(sin(twe-tww));

    //r += -8.0 * 0.5 * dgtw(s,tilt(s,i,0)) 
    //  * pow(0.5 * (twist(s,i+1,0) - twist(s,i-1,0)),2);
    r += -8.0 * 0.5 * dgtw(s,tilt(s,i,0)) 
      * pow(0.5 * (aemw),2);
    return r;

}


void TiltEquation :: evaluate (VEC rhs, VEC ju, VEC jl, VEC jd, State* s){
  // euler-lagrange equation for 
  // tilt in cell i

  int n = s->n_lc_layers();

  double kk = K(s,3)/K(s,1);
  double ee = E0/K(s,1);
  double ff =  (E0 * flexo(s,1))/K(s,1);
  double a2 = alpha(s,2)/K(s,1);
  double a3 = alpha(s,3)/K(s,1);
  double gamma = a3 - a2;

  
  for (int i = 1; i < n-1; i++){

   
      
    // residual

    double src = steadyrhs(s,i,0);
    //    double ftime(s) = 0.5;

    //double stf = 0.0;

    if (s->time_dependent()) {
      // time dependent term
      src *= ftime(s);
      src += (1-ftime(s)) * steadyrhs(s,i,-1); 
      src += (8*gamma*Power(dz(s),2)*(tilt(s,i,-1) - tilt(s,i,0)))/dt(s);

    }

    (*rhs)[i] = src;


    //diagonal of jacobian df(i)/dx(i)
    //time independent terms

    double diff = - 8* ff * Cos(2*tilt(s,i,0)) *(-field_dc_edge(s,i) + field_dc_edge(s,i-1));
    double d_flexo = min(0.0,diff)     ;

    double d  =
      
      2*( d_flexo - 
	  8*Power(Cos(tilt(s,i,0)),2)*kk - 8 * Power(Sin(tilt(s,i,0)),2)  - 
	  2*a2*dz(s)*(flowx(s,-1 + i,0) - flowx(s,1 + i,0))*Sin(2*tilt(s,i,0)) - 
	  2*a3*dz(s)*(flowx(s,-1 + i,0) - flowx(s,1 + i,0))*Sin(2*tilt(s,i,0)) +
	  Cos(2*tilt(s,i,0))*Power(tilt(s,-1 + i,0) - tilt(s,1 + i,0),2) - 
	  Cos(2*tilt(s,i,0))*kk*Power(tilt(s,-1 + i,0) - tilt(s,1 + i,0),2) + 
	  8*Cos(tilt(s,i,0))*Sin(tilt(s,i,0))*
	  (tilt(s,-1 + i,0) - 2*tilt(s,i,0) + tilt(s,1 + i,0)) - 
	  8*Cos(tilt(s,i,0))*kk*Sin(tilt(s,i,0))*
	  (tilt(s,-1 + i,0) - 2*tilt(s,i,0) + tilt(s,1 + i,0)) - 
	  ee*Cos(2*tilt(s,i,0))*ea(s)*Power(volt_ac(s,-1 + i,0) - volt_ac(s,1 + i,0),2) 
	  - ee*Cos(2*tilt(s,i,0))*ea(s)*Power(2.0 * field_dc_centre(s,i),2))
      ;
    
    (*jd)[i] = d;


    

    if (s->time_dependent()) {
      // add time dependent terms
      (*jd)[i] *= ftime(s);
      (*jd)[i] +=  -8*(gamma)*Power(dz(s),2)/dt(s);
    }   

    

    //sub-diagonal df(i)/dx(i-1)
    if (i > 0) (*ju)[i-1] =
		 8*Power(Cos(tilt(s,i,0)),2)*kk + 8*Power(Sin(tilt(s,i,0)),2) + 
		 2*Sin(2*tilt(s,i,0))*(tilt(s,-1 + i,0) - tilt(s,1 + i,0)) - 
		 2*kk*Sin(2*tilt(s,i,0))*(tilt(s,-1 + i,0) - tilt(s,1 + i,0))
		 ;	
     
    //super-diagonal  df(i)/dx(i+1)
    if (i < n-1) (*jl)[i] = 
		   8*Power(Cos(tilt(s,i,0)),2)*kk + 8*Power(Sin(tilt(s,i,0)),2) - 
		   2*Sin(2*tilt(s,i,0))*(tilt(s,-1 + i,0) - tilt(s,1 + i,0)) + 
		   2*kk*Sin(2*tilt(s,i,0))*(tilt(s,-1 + i,0) - tilt(s,1 + i,0))
		   ;

    if (s->time_dependent()){
      (*ju)[i-1] *= ftime(s);
      (*jl)[i] *= ftime(s);

    }
    
  }

  // boundary conditions (still need to include surface viscosity)
  
  // z = 0 surface
  if (w(s,0) == 0) {
    // strong anchoring 
    (*rhs)[0] = n*BIG*(tilt(s,0,0) - boundary_tilt(s,0));
    (*jd)[0]  = n*BIG;
    (*jl)[0] = 0 * 1.0;
  }
  else {

    //weak anchoring
    
    if (fabs(tilt(s,0,0)-boundary_tilt(s,0)) < M_PI/4 - 0.1){

#define SVISC 0

      double src = (1.0 + (K(s,3) - K(s,1)) * pow(cos(tilt(s,0,0)),2)) * (-tilt(s,0,0) + tilt(s,1,0))
	- 0.5 * E0 * flexo(s,1)/K(s,1) * sin (2.0 * tilt(s,0,0)) * field_dc_edge(s,0)
	- 0.5 * dz(s) * w(s,0)/K(s,1) * sin(2.0*(tilt(s,0,0)-boundary_tilt(s,0)));

      if (s->time_dependent()){
      	src  -= SVISC  * dz(s) / (dt(s) * K(s,1))  * (tilt(s,0,0) - tilt(s,0,-1));
	}


      (*rhs)[0] = src;

      double d =  - (1.0 + (K(s,3) - K(s,1)) * pow(cos(tilt(s,0,0)),2))
	- 0.5 * E0 * flexo(s,1)/K(s,1) * 2.0 * cos (2.0 * tilt(s,0,0)) * field_dc_edge(s,0)
	- 0.5 * dz(s) * w(s,0)/K(s,1) * 2.0 * cos(2.0*(tilt(s,0,0)-boundary_tilt(s,0)));

      if (s->time_dependent()){
	d +=  - SVISC * dz(s) / (dt(s) * K(s,1));
      }

      (*jd)[0] = d;
      (*jl)[0] = 1.0 + (K(s,3) - K(s,1)) * pow(cos(tilt(s,0,0)),2);

    } else {
      (*rhs)[0] =
	(-tilt(s,0,0) + tilt(s,1,0)) - dz(s)
	*Power(tilt(s,0,0)-boundary_tilt(s,0),2)*w(s,0)/K(s,1)
	;
	 
	 
      (*jd)[0] =
	- 1.0
	- 2*(tilt(s,0,0)-boundary_tilt(s,0))*dz(s)
	*w(s,0)/K(s,1)
	;

      (*jl)[0] = 0 * 1.0;
	 
   }

  }

  // z = 1 surface
  if (w(s,1) == 0) {
    // strong anchoring
    (*rhs)[n-1] = n*BIG*(tilt(s,n-1,0) - boundary_tilt(s,1));
    (*jd)[n-1] = n*BIG;
    (*ju)[n-2] = 0;
  }
  else {

    int i = n-1;

    if (fabs(tilt(s,i,0)-boundary_tilt(s,1)) <  M_PI/4 - 0.1){
      

      if (false){
      (*rhs)[i] = 
	(-tilt(s,i-1,0) + tilt(s,i,0)) - dz(s)
	*Power(Sin(tilt(s,i,0)-boundary_tilt(s,1)),2)*w(s,1)/K(s,1)
	;
      (*jd)[i] =
	+ 1.0
	- 2*Cos(tilt(s,i,0)-boundary_tilt(s,1))*dz(s)
	* Sin(tilt(s,i,0)-boundary_tilt(s,1))*w(s,1)/K(s,1)
	;
      (*ju)[i-1] = -1.0;
      (*jl)[i-1] = 0.0;
      }
      else {

      double src = (1.0 + (K(s,3) - K(s,1)) * pow(cos(tilt(s,i,0)),2)) * (- tilt(s,i-1,0) + tilt(s,i,0))
	- 0.5 * E0 * flexo(s,1)/K(s,1) * sin (2.0 * tilt(s,i,0)) * field_dc_edge(s,i)
	+ 0.5 * dz(s) * w(s,1)/K(s,1) * sin(2.0*(tilt(s,i,0)-boundary_tilt(s,1)));

      if (s->time_dependent()){
      	src  -= SVISC  * dz(s) / (dt(s) * K(s,1))  * (tilt(s,i,0) - tilt(s,i,-1));
      }


      (*rhs)[i] = src;

      double d =    (1.0 + (K(s,3) - K(s,1)) * pow(cos(tilt(s,i,0)),2))
	- 0.5 * E0 * flexo(s,1)/K(s,1) * 2.0 * cos (2.0 * tilt(s,i,0)) * field_dc_edge(s,i)
	+ 0.5 * dz(s) * w(s,1)/K(s,1) * 2.0 * cos(2.0*(tilt(s,i,0)-boundary_tilt(s,1)));

      if (s->time_dependent()){
	d +=  - SVISC * dz(s) / (dt(s) * K(s,1));
      }

      (*jd)[i] = d;
      (*jl)[i-1] = 0;
      (*ju)[i-1] = - (1.0 + (K(s,3) - K(s,1)) * pow(cos(tilt(s,i,0)),2));

      }
    }
    else {
      (*rhs)[i] = 
	(-tilt(s,i-1,0) + tilt(s,i,0)) - dz(s)
	*Power(tilt(s,i,0)-boundary_tilt(s,1),2)*w(s,1)/K(s,1)	      
	;  
      (*jd)[i] =
	+ 1.0
	- 2*(tilt(s,i,0)-boundary_tilt(s,1))*dz(s)
	*w(s,1)/K(s,1)
	; 
      (*ju)[i-1] = -1.0;

    }

    

  }

}



// tilt equation, director & flow confined to xz plane

void TiltEquationXZ :: update(VEC dtilt, State*s)
{
  double dtmax = M_PI/4;
  double a = 1.0;
 //  for (size_t i = 0; i < s->n_lc_layers(); i++){
//     double dt = fabs((*dtilt)[i]);
//     a = (a > dtmax/dt)?dtmax/dt:a;
//   }
  

  for (size_t i = 0; i < s->n_lc_layers(); i++){
    s->tilt[i] -= a * (*dtilt)[i];
  }
}

 double TiltEquationXZ :: difference(State *a, State *b)
{
  double diff = 0;
  double sum = 0;

   for (size_t k = 0; k  < b->n_lc_layers(); k++) {
     diff += pow(tilt(b,k,0) - tilt(a,k,0) ,2);
     sum +=  pow(tilt(b,k,0) + tilt(a,k,0) ,2);
    
   }	
   if (sum < TINY) sum = 1.0;

   return diff/sum;

}

inline double TiltEquationXZ :: steadyrhs(State *s, size_t i,int j){
  
  double kk = K(s,3)/K(s,1);
  double ee = E0/K(s,1);
  double ff =  (E0 * flexo(s,1))/K(s,1);
  double a2 = alpha(s,2)/K(s,1);
  double a3 = alpha(s,3)/K(s,1);
  //double gamma = a3 - a2;

 double r =   
   
   (4*a2*Power(Cos(tilt(s,i,j)),2)*dz(s)*
    (flowx(s,-1 + i,j) - flowx(s,1 + i,j)) + 
    4*a3*dz(s)*(-flowx(s,-1 + i,j) + flowx(s,1 + i,j))
    *Power(Sin(tilt(s,i,j)),2)) 
   +
   Sin(2*tilt(s,i,j))*Power(tilt(s,-1 + i,j) - tilt(s,1 + i,j),2) - 
   kk*Sin(2*tilt(s,i,j))*Power(tilt(s,-1 + i,j) - tilt(s,1 + i,j),2) + 
   8*Power(Cos(tilt(s,i,j)),2)*kk*(tilt(s,-1 + i,j) - 2*tilt(s,i,j) + tilt(s,1 + i,j)) + 
   8*Power(Sin(tilt(s,i,j)),2)*(tilt(s,-1 + i,j) - 2*tilt(s,i,j) + tilt(s,1 + i,j)) 
   -ee*ea(s)*Sin(2*tilt(s,i,j))*Power(volt_ac(s,-1 + i,j) - volt_ac(s,1 + i,j),2)
   -ee*ea(s)*Sin(2*tilt(s,i,j))*Power(2.*field_dc_centre(s,i),2)
   ;

    //flexo terms
    r += - 4* ff *Sin(2*tilt(s,i,0)) 
      * (-field_dc_edge(s,i)  + field_dc_edge(s,i-1)) ;
    return r;

}


void TiltEquationXZ :: evaluate (VEC rhs, VEC ju, VEC jl, VEC jd, State* s){
  // euler-lagrange equation for 
  // tilt in cell i

  int n = s->n_lc_layers();

  double kk = K(s,3)/K(s,1);
  double ee = E0/K(s,1);
  double ff =  (E0 * flexo(s,1))/K(s,1);
  double a2 = alpha(s,2)/K(s,1);
  double a3 = alpha(s,3)/K(s,1);
  double gamma = a3 - a2;

  
  for (int i = 1; i < n-1; i++){

   
      
    // residual

    double src = steadyrhs(s,i,0);
    //    double ftime(s) = 0.5;

    //double stf = 0.0;

    if (s->time_dependent()) {
      // time dependent term
      src *= ftime(s);
      src += (1-ftime(s)) * steadyrhs(s,i,-1); 
      src += (8*gamma*Power(dz(s),2)*(tilt(s,i,-1) - tilt(s,i,0)))/dt(s);

    }

    (*rhs)[i] = src;


    //diagonal of jacobian df(i)/dx(i)
    //time independent terms

    double diff = - 8* ff * Cos(2*tilt(s,i,0)) *(-field_dc_edge(s,i) + field_dc_edge(s,i-1));
    double d_flexo = min(0.0,diff)     ;

    double d  =
      
      2*( d_flexo - 
	  8*Power(Cos(tilt(s,i,0)),2)*kk - 8 * Power(Sin(tilt(s,i,0)),2)  
	  - 2*(a2+a3)*dz(s)*(flowx(s,-1 + i,0) - flowx(s,1 + i,0))*Sin(2*tilt(s,i,0)) 
	  +
	  Cos(2*tilt(s,i,0))*Power(tilt(s,-1 + i,0) - tilt(s,1 + i,0),2) - 
	  Cos(2*tilt(s,i,0))*kk*Power(tilt(s,-1 + i,0) - tilt(s,1 + i,0),2) + 
	  8*Cos(tilt(s,i,0))*Sin(tilt(s,i,0))*
	  (tilt(s,-1 + i,0) - 2*tilt(s,i,0) + tilt(s,1 + i,0)) - 
	  8*Cos(tilt(s,i,0))*kk*Sin(tilt(s,i,0))*
	  (tilt(s,-1 + i,0) - 2*tilt(s,i,0) + tilt(s,1 + i,0)) - 
	  ee*Cos(2*tilt(s,i,0))*ea(s)*Power(volt_ac(s,-1 + i,0) - volt_ac(s,1 + i,0),2) 
	  - ee*Cos(2*tilt(s,i,0))*ea(s)*Power(2.0 * field_dc_centre(s,i),2))
      ;
    
    (*jd)[i] = d;


    

    if (s->time_dependent()) {
      // add time dependent terms
      (*jd)[i] *= ftime(s);
      (*jd)[i] +=  -8*(gamma)*Power(dz(s),2)/dt(s);
    }   

    

    //sub-diagonal df(i)/dx(i-1)
    if (i > 0) (*ju)[i-1] =
		 8*Power(Cos(tilt(s,i,0)),2)*kk + 8*Power(Sin(tilt(s,i,0)),2) + 
		 2*Sin(2*tilt(s,i,0))*(tilt(s,-1 + i,0) - tilt(s,1 + i,0)) - 
		 2*kk*Sin(2*tilt(s,i,0))*(tilt(s,-1 + i,0) - tilt(s,1 + i,0))
		 ;	
     
    //super-diagonal  df(i)/dx(i+1)
    if (i < n-1) (*jl)[i] = 
		   8*Power(Cos(tilt(s,i,0)),2)*kk + 8*Power(Sin(tilt(s,i,0)),2) - 
		   2*Sin(2*tilt(s,i,0))*(tilt(s,-1 + i,0) - tilt(s,1 + i,0)) + 
		   2*kk*Sin(2*tilt(s,i,0))*(tilt(s,-1 + i,0) - tilt(s,1 + i,0))
		   ;

    if (s->time_dependent()){
      (*ju)[i-1] *= ftime(s);
      (*jl)[i] *= ftime(s);

    }
    
  }

  // boundary conditions (still need to include surface viscosity)
  
  // z = 0 surface
  if (w(s,0) == 0) {
    // strong anchoring 
    (*rhs)[0] = n*BIG*(tilt(s,0,0) - boundary_tilt(s,0));
    (*jd)[0]  = n*BIG;
    (*jl)[0] = 0 * 1.0;
  }
  else {

    //weak anchoring

    if (fabs(tilt(s,0,0)-boundary_tilt(s,0)) < M_PI/4 - 0.1){

#define SVISC 0

      double src = (1.0 + (K(s,3) - K(s,1)) * pow(cos(tilt(s,0,0)),2)) * (-tilt(s,0,0) + tilt(s,1,0))
	- 0.5 * E0 * flexo(s,1)/K(s,1) * sin (2.0 * tilt(s,0,0)) * field_dc_edge(s,0)
	- 0.5 * dz(s) * w(s,0)/K(s,1) * sin(2.0*(tilt(s,0,0)-boundary_tilt(s,0)));

      if (s->time_dependent()){
      	src  -= SVISC  * dz(s) / (dt(s) * K(s,1))  * (tilt(s,0,0) - tilt(s,0,-1));
	}


      (*rhs)[0] = src;

      double d =  - (1.0 + (K(s,3) - K(s,1)) * pow(cos(tilt(s,0,0)),2))
	- 0.5 * E0 * flexo(s,1)/K(s,1) * 2.0 * cos (2.0 * tilt(s,0,0)) * field_dc_edge(s,0)
	- 0.5 * dz(s) * w(s,0)/K(s,1) * 2.0 * cos(2.0*(tilt(s,0,0)-boundary_tilt(s,0)));

      if (s->time_dependent()){
	d +=  - SVISC * dz(s) / (dt(s) * K(s,1));
      }

      (*jd)[0] = d;
      (*jl)[0] = 1.0 + (K(s,3) - K(s,1)) * pow(cos(tilt(s,0,0)),2);

    } else {
      (*rhs)[0] =
	(-tilt(s,0,0) + tilt(s,1,0)) - dz(s)
	*Power(tilt(s,0,0)-boundary_tilt(s,0),2)*w(s,0)/K(s,1)
	;
	 
	 
      (*jd)[0] =
	- 1.0
	- 2*(tilt(s,0,0)-boundary_tilt(s,0))*dz(s)
	*w(s,0)/K(s,1)
	;

      (*jl)[0] = 0 * 1.0;
	 
   }

  }

  // z = 1 surface
  if (w(s,1) == 0) {
    // strong anchoring
    (*rhs)[n-1] = n*BIG*(tilt(s,n-1,0) - boundary_tilt(s,1));
    (*jd)[n-1] = n*BIG;
    (*ju)[n-2] = 0;
  }
  else {

    int i = n-1;

    if (fabs(tilt(s,i,0)-boundary_tilt(s,1)) <  M_PI/4 - 0.1){
      

      if (false){
      (*rhs)[i] = 
	(-tilt(s,i-1,0) + tilt(s,i,0)) - dz(s)
	*Power(Sin(tilt(s,i,0)-boundary_tilt(s,1)),2)*w(s,1)/K(s,1)
	;
      (*jd)[i] =
	+ 1.0
	- 2*Cos(tilt(s,i,0)-boundary_tilt(s,1))*dz(s)
	* Sin(tilt(s,i,0)-boundary_tilt(s,1))*w(s,1)/K(s,1)
	;
      (*ju)[i-1] = -1.0;
      (*jl)[i-1] = 0.0;
      }
      else {

      double src = (1.0 + (K(s,3) - K(s,1)) * pow(cos(tilt(s,i,0)),2)) * (- tilt(s,i-1,0) + tilt(s,i,0))
	- 0.5 * E0 * flexo(s,1)/K(s,1) * sin (2.0 * tilt(s,i,0)) * field_dc_edge(s,i)
	+ 0.5 * dz(s) * w(s,1)/K(s,1) * sin(2.0*(tilt(s,i,0)-boundary_tilt(s,1)));

      if (s->time_dependent()){
      	src  -= SVISC  * dz(s) / (dt(s) * K(s,1))  * (tilt(s,i,0) - tilt(s,i,-1));
      }


      (*rhs)[i] = src;

      double d =    (1.0 + (K(s,3) - K(s,1)) * pow(cos(tilt(s,i,0)),2))
	- 0.5 * E0 * flexo(s,1)/K(s,1) * 2.0 * cos (2.0 * tilt(s,i,0)) * field_dc_edge(s,i)
	+ 0.5 * dz(s) * w(s,1)/K(s,1) * 2.0 * cos(2.0*(tilt(s,i,0)-boundary_tilt(s,1)));

      if (s->time_dependent()){
	d +=  - SVISC * dz(s) / (dt(s) * K(s,1));
      }

      (*jd)[i] = d;
      (*jl)[i-1] = 0;
      (*ju)[i-1] = - (1.0 + (K(s,3) - K(s,1)) * pow(cos(tilt(s,i,0)),2));

      }
    }
    else {
      (*rhs)[i] = 
	(-tilt(s,i-1,0) + tilt(s,i,0)) - dz(s)
	*Power(tilt(s,i,0)-boundary_tilt(s,1),2)*w(s,1)/K(s,1)	      
	;  
      (*jd)[i] =
	+ 1.0
	- 2*(tilt(s,i,0)-boundary_tilt(s,1))*dz(s)
	*w(s,1)/K(s,1)
	; 
      (*ju)[i-1] = -1.0;

    }

    

  }

}


State* fixed_step_compute_state(vector<Equation*> e,
				Parameters* p, Geometry* g,  State* ps,
				double dt,	double newton_tol,  int newton_iter)
{
  // compute s(t + dt) given fixed dt and ps = s(t)
  TriDiagWorkSpace tdws(g->n_cells()); 
  State *s = new State(ps,dt,p,g);
  semi_implicit_solve(e,s,newton_tol, newton_iter, &tdws);
  return s;

}


bool test_difference(State *a, State *b, vector<Equation *> e, double eps){
  // compare states a and b...
  bool same = true;
  for (size_t j = 0; j < e.size(); j++) {
    same &= (e[j]->difference(a,b) < eps);     
  }
  return same;
}


bool test_difference(State *a, State *b,double eps){
  //measure of difference between states...
  // really needs a tidy up, say give each Equation a
  // difference(a,b,eps) 

  double diff_tilt = 0;
  double sum_tilt = 0;

	
  for (size_t k = 0; k  < b->n_lc_layers(); k++) {
    diff_tilt += pow(b->tilt[k] - a->tilt[k],2);
    sum_tilt +=  pow(b->tilt[k] + a->tilt[k],2);
       
  }
	
  if (sum_tilt < eps) sum_tilt = 1.0;

  double diff_flow = 0;
  double sum_flow = 0;

	
  for (size_t k = 0; k  < b->n_lc_layers(); k++) {
    diff_flow += pow(b->flowx[k] - a->flowx[k],2);
    sum_flow +=  pow(b->flowx[k] + a->flowx[k],2);
       
  }
  //diff_flow = 0;
  if (sum_flow < 1e-8) sum_flow = 1.0;

  double diff_volt_ac = 0;
  double sum_volt_ac = 0;

  for (size_t k = 0; k  < b->n_lc_layers(); k++) {
    diff_volt_ac += pow(volt_ac(b,k,0) - volt_ac(a,k,0) ,2);
    sum_volt_ac +=  pow(volt_ac(b,k,0) + volt_ac(a,k,0) ,2 );
       
  }	

  if (sum_volt_ac < 1e-8) sum_volt_ac = 1.0;

  double diff_volt_dc = 0;
  double sum_volt_dc = 0;

  for (size_t k = 0; k  < b->n_lc_layers(); k++) {
    diff_volt_dc += pow(volt_dc(b,k,0) - volt_dc(a,k,0) ,2);
    sum_volt_dc +=  pow(volt_dc(b,k,0) + volt_dc(a,k,0) ,2 );
       
  }	

  if (sum_volt_dc < 1e-8) sum_volt_dc = 1.0;

  

  double diff_charge = 0;
  double sum_charge = 0;
	
  for (size_t k = 0; k  < b->n_lc_layers(); k++) {
    diff_charge += pow(charge(b,k,0) - charge(a,k,0) ,2);
    sum_charge +=  pow(charge(b,k,0) + charge(a,k,0) ,2 );
	  
  }	
	
  if (sum_charge < 1e-8) sum_charge = 1.0;
	
  //diff_charge=0;

	
  return  (diff_tilt/sum_tilt < eps)  
    && (diff_volt_ac/sum_volt_ac < eps)
    && (diff_volt_dc/sum_volt_dc < eps)
    && (diff_flow/sum_flow < eps) 
    && (diff_charge/sum_charge < eps);
	
}

State* refine_step_compute_state(vector<Equation*> e,
				 Parameters* p, Geometry* g,  State* s, State* ps,
				 double dt,	double newton_tol,  int newton_iter, 
				 double time_tol, double min_timestep){

  // refine estimate of the state s(t + dt), where previous
  // is ps(t)
  
  // the time mesh is refined recusively.

  //bool s_local = false;
  if (s == 0) 
    s = fixed_step_compute_state(e,p,g,ps,dt,newton_tol,newton_iter);
  else
    s = new State(s);

  if (dt < 2.0 * min_timestep) 
    return s;

  State *sa = fixed_step_compute_state(e,p,g,ps,dt/2.0,newton_tol,newton_iter);
  State *sb = fixed_step_compute_state(e,p,g,sa,dt - dt/2.0,newton_tol,newton_iter);
    
  if (test_difference(s,sb,e,time_tol)){
    delete(sa);
    delete(s);
    sb->set_previous(ps);
    return sb;
  }

  
  State *sc = refine_step_compute_state
    (e,p,g,sb,
     refine_step_compute_state
     (e,p,g,sa,ps,dt/2.0,newton_tol,newton_iter,time_tol,min_timestep),
     dt - dt/2.0,newton_tol,newton_iter,time_tol,min_timestep);

  delete(sa);
  delete(sb);
  delete(s);
  sc->set_previous(ps);
  return sc;


}


State* compute_state(vector<Equation*> e,
		     Parameters* p, Geometry* g,  State* previous,
		     double dt,	double newton_tol,  int newton_iter, 
		     double time_tol, double min_timestep)
{
  // solve the a time-dependent  or equilibrium 
  // equations using Newton's method.
  // return a state representing the solution

  // if  time-dependent (ie dt > 0 is specified) 
  // the state is determined by succesively refining the time step 

  // accepts:

 
  // int newton_tol: tolerance for Newton iterations, test 
  //                 is perfromed on the step-size for all
  //                 equations

  // int newton_iter: the maximum number of newton iterations
  //                  to be attempted

  // time_tol       : tolerance for time resolution
  // min_timestep   : minimum size for the succesively refined
  //                  timestep


 
  TriDiagWorkSpace tdws(g->n_cells()); 

  
  State* s = 0;

  if (dt > 0){

#define RECURSIVE
     
#ifdef RECURSIVE
    s = refine_step_compute_state(e,p, g,  0, previous,
				  dt,	newton_tol,  newton_iter, time_tol,  min_timestep);
#else
    bool resolved = false;
    vector<State*> states_a, states_b;

    int steps = 1;

    while ((!resolved) && (dt > min_timestep)) {
                  
      states_a.push_back(previous);

      for (int j = 0; j < steps; j++) {

	State* st = new State(states_a.back(),dt,p,g);
	states_a.push_back(st);
	semi_implicit_solve(e,st,newton_tol, newton_iter, &tdws);

      }


      // compare the final state of the new sequence
      // with the the final state of the old sequence
      // (after the first iteration, when no comparison 
      // exists
      if (states_b.size() > 0) {

	State state_a = *states_a.back();
	State state_b = *states_b.back();
	
	double diff_tilt = 0;
	double sum_tilt = 0;

	
	for (int k = 0; k  < state_b.n_lc_layers(); k++) {
	  diff_tilt += pow(state_b.tilt[k] - state_a.tilt[k],2);
	  sum_tilt +=  pow(state_b.tilt[k] + state_a.tilt[k],2);
       
	}
	
	if (sum_tilt < time_tol) sum_tilt = 1.0;

	double diff_flow = 0;
	double sum_flow = 0;

	
	for (int k = 0; k  < state_b.n_lc_layers(); k++) {
	  diff_flow += pow(state_b.flowx[k] - state_a.flowx[k],2);
	  sum_flow +=  pow(state_b.flowx[k] + state_a.flowx[k],2);
       
	}
	//diff_flow = 0;
	if (sum_flow < 1e-8) sum_flow = 1.0;

	double diff_volt = 0;
	double sum_volt = 0;

	for (int k = 0; k  < state_b.n_lc_layers(); k++) {
	  diff_volt += pow(volt(&state_b,k,0) - volt(&state_a,k,0) ,2);
	  sum_volt +=  pow(volt(&state_b,k,0) + volt(&state_a,k,0) ,2 );
       
	}	

	if (sum_volt < 1e-8) sum_volt = 1.0;


	double diff_charge = 0;
	double sum_charge = 0;
	
	for (int k = 0; k  < state_b.n_lc_layers(); k++) {
	  diff_charge += pow(charge(&state_b,k,0) - charge(&state_a,k,0) ,2);
	  sum_charge +=  pow(charge(&state_b,k,0) + charge(&state_a,k,0) ,2 );
	  
	}	
	
	if (sum_charge < 1e-8) sum_charge = 1.0;
	
	//diff_charge=0;

	
	resolved =   (((diff_tilt/sum_tilt < time_tol)  
		       && ((diff_volt/sum_volt < time_tol)
			   && (diff_flow/sum_flow < time_tol))) 
		      && (diff_charge/sum_charge < time_tol));
	  
	//cerr << "diff = " << diff_tilt/sum_tilt << endl;
	//resolved = true;
	

      }
      dt = 0.5*dt;
      steps *= 2;
      //free memory in sequence b (expect first pointer, which is previous)
      for (int i = 1; i < states_b.size(); i++) delete(states_b[i]);
	     
      // copy the new sequence to the old
      states_b.assign(states_a.begin(),states_a.end());
      states_a.resize(0);
    }

    s = states_b.back();
    s->set_previous(previous);

    // free all but first and last member of sequence b
    for (int i = 1; i < states_b.size() - 1; i++){
      delete(states_b[i]);
    }
#endif
  }
  else {

    // steady state solution.
    //State *s;

    if (previous == 0)
      s = new State(g->n_cells(),p,g);
    else
      s = new State(previous,0,p,g);
    //bool converged = 
    semi_implicit_solve(e,s,newton_tol, newton_iter, &tdws);
    

  }
    
  for (size_t i = 1; i < s->n_lc_layers(); i++){
    if (s->twist[i] - s->twist[i-1] > (M_PI * 2)){
      //s->twist[i] -=  2 * M_PI;
    }
  }

  
  return s;
   
}


inline bool ion_model(Parameters *p){
  // does the parameter set specify an ion model?
  return ((p->rho0.size() > 0) && (p->rho0[0] > 0));
 
}


vector<State*> compute_states 
(vector<double> const& times, vector<double> const& voltages, vector<double> const& frequencies,
 vector<double> const& pgx, vector<double> const& pgy, vector<State*> const& initial_states,
 bool shear, bool shear_y,  Parameters *p, Geometry *g, double newton_tol, int newton_iter, 
 double time_tol, double min_timestep)
{
  // given a vectors of times & voltages. , 
  // solve the Leslie-Ericksen equations for a 
  // series of lcd states.


  // vector<State*> states : input/ouput
  //   
  //    on entry, each element of states points
  //    to a state object which will define
  //    lcd parameters, geometry, and a previous
  //    lcd state.

  // set up equations

  vector<State*> states;
  vector<Equation*> e;

  
  // need to include an ac voltage equation only if
  // non-zero frequency is applied
  bool ac = false;
  size_t i = 0;
  while (!ac && (i < times.size())){
    ac |= !(frequencies[i++] == 0.0);
  }
	 
  if (ac)
    e.push_back(new Volt_ACEquation());

  // need to include a dc voltage equation if flexoelectricity
  // or ions are considered
  bool dc = false;

  if (ion_model(p)) {

    // add an ion transport equation for each non-zero density species
    for (size_t sp = 0; sp < p->rho0.size(); sp++){
      if (p->rho0[sp] > 0){
	e.push_back(new IonEquation(sp));
      }
    }
  }

  if (((abs(p->flexo[1]) > 1e-6) || (ion_model(p)))){
      e.push_back(new Volt_DCEquation());
      dc = true;
  }

  
 
 


  if (shear){
    if (p->boundary_twist[0] != p->boundary_twist[1]){
      // both momentum equations needed
         e.push_back(new FlowxEquation());
	 if (shear_y){
	   //unless flagged otherwise
	   e.push_back(new FlowyEquation());
	 } 
    } else {
      if (abs(cos((p->boundary_twist[0]))) > 1e-6){
	e.push_back(new FlowxEquationXZ());
      }
      if (abs(sin((p->boundary_twist[0]))) > 1e-6){
	if (shear_y){
	  e.push_back(new FlowyEquation());
	}
      }      
    }

  }

  if ((p->boundary_twist[0] != p->boundary_twist[1]) ||
      (p->boundary_twist[0] != 0)){
     e.push_back(new TiltEquation());
    e.push_back(new TwistEquation());
   
  } else {
    e.push_back(new TiltEquationXZ());
  }

 

  

   
  for (size_t i = 0; i < times.size(); i++) {

    Parameters* pi = p;

    if (dc && (frequencies[i] == 0.0)){
      pi->boundary_volt_ac[0] = 0.0; 
      pi->boundary_volt_dc[0] = voltages[i];
    }
    else{
      pi->boundary_volt_dc[0] = 0.0; 
      pi->boundary_volt_ac[0] = voltages[i]; 
    }

    pi->boundary_volt_ac[1] = 0.0;
    pi->boundary_volt_dc[1] = 0.0;
   
    pi->pgx = pgx[i];
    pi->pgy = pgy[i];

    if ((i > 0) && (times[i] > times[i-1])){

      State* s = compute_state
	(e, pi, g, states.back() , times[i] - times[i-1] ,
	 newton_tol,  newton_iter , time_tol , min_timestep);

     
      states.push_back(s);
      
    } 
    else
      {
	State* s = compute_state
	  (e, pi, g, initial_states[i] , 0.0 , newton_tol,  newton_iter , 0.0 , 0);
	states.push_back(s);
      }

    

  }

  for (size_t i = 0; i < e.size(); i++){
    delete(e.at(i));
  }
 
  return states;

}



//constructor for NewtonResult
NewtonResult::NewtonResult(){residual = 0;stepsize = 0;}

//constructor for TriDiagWorkSpace
TriDiagWorkSpace::TriDiagWorkSpace(int n){
  x = new valarray<double>(n);
  r = new valarray<double>(n);
  jd = new valarray<double>(n);
  jl = new valarray<double>(n-1);
  ju = new valarray<double>(n-1);



}

//destructor for TriDiagWorkSpace
TriDiagWorkSpace::~TriDiagWorkSpace(){
  delete(x);
  delete(r);
  delete(jd);
  delete(ju);
  delete(jl);
  
}


inline double norm2(VEC x){
  double d = 0;
  for (size_t i = 0; i < x->size(); i++){
    d += pow((*x)[i],2);
  }
  return sqrt(d);
}

int solve_tridiag (VEC b,VEC c,VEC a,VEC d, VEC x){

  //thomas alogorithm
  int n = x->size();

  for (int i = 0; i < n - 1; i++){
    double q = (*a)[i] / (*b)[i];

    (*b)[i+1] = (*b)[i+1] - q * (*c)[i];
    (*d)[i+1] = (*d)[i+1] - q * (*d)[i];
  }

  (*x)[n-1] = (*d)[n-1]/(*b)[n-1];

  for (int i = n-2; i >= 0; i--){
    (*x)[i] = ((*d)[i] - (*c)[i]*(*x)[i+1])/(*b)[i];

  }


  return 0;

}

// int gsl_linalg_solve_tridiag (VEC b,VEC c,VEC a,VEC d,VEC x){

//   //thomas alogorithm

//   int n = x->size;

//   for (int i = 0; i < n - 1; i++){
//     double q = gsl_vector_get(a,i) / gsl_vector_get(b,i);

//     cout << i 
// 	 << " " << gsl_vector_get(a,i) 
// 	 << " " << gsl_vector_get(b,i) 
// 	 << " " << q << endl;
//     gsl_vector_set(b,i+1,gsl_vector_get(b,i+1) - q * gsl_vector_get(c,i));
//     gsl_vector_set(d,i+1,gsl_vector_get(d,i+1) - q * gsl_vector_get(d,i));
//   }

//   gsl_vector_set(x,n-1,gsl_vector_get(d,n-1)/
// 		 gsl_vector_get(b,n-1));

//   for (int i = n-2; i >= 0; i--){
//     gsl_vector_set(x,i,(gsl_vector_get(d,i) -
// 			gsl_vector_get(c,i)*gsl_vector_get(x,i+1))
// 		   /gsl_vector_get(b,i));
//      cout << i 
// 	  << " " << gsl_vector_get(x,i) << endl;

//   }


//   return 0;

// }


NewtonResult tri_newton(int offset,int n, Equation* eq,State* s,
			TriDiagWorkSpace* tdws)

{
  //  perform a newton step toward
  //  solving f[i](x[j]) = 0, assuming
  //  a tri-diagonal jacobian df[i]/dx[j]

  NewtonResult r;
  
  // compute f(x), jacobian and rms residual

  

  eq->evaluate(tdws->r,tdws->jl,tdws->ju,tdws->jd,s)  ;
  r.residual = norm2(tdws->r)/tdws->r->size();


  // compute newton step by solving jx = rhs
    if (solve_tridiag (tdws->jd,tdws->ju,tdws->jl,
  			tdws->r,tdws->x) == 0){
  
      
      r.stepsize = norm2(tdws->x)/tdws->x->size();
      eq->update(tdws->x,s);

      eq->residual.push_back(r.residual);
      eq->stepsize.push_back(r.stepsize);

      

      }

  return r;

}

bool semi_implicit_solve(vector<Equation*> e, State* s,
			 double tol, int niter, 
			 TriDiagWorkSpace* tdws)
{
  // semi-implict Newton scheme for a sequence of equations
  // whose Jacobian is tri-diagonal
  // each equation is updated in turn by a single newton step
  
  // vector<Equation> e  : the sequence of equations
  // State*s             : state upon which equations act
  // double tol          : desired tolerance, measured 
  //                       as the norm of the update size
  // int niter           : the maximun number of iterations


  int n = s->n_lc_layers();
  bool converged = false;
  int i = 0;
 
  while ((i++ < 2) || ((!converged) && (i < niter))){
    s->iters = i+1;
    converged = true;
    for (size_t j = 0; j < e.size(); j++) {
      NewtonResult r = tri_newton(0,n,e[j],s,tdws);
     
      converged &= r.stepsize < tol;
      //cout << i << " " << e[j]->name() << " " << r.residual << " " << r.stepsize << endl;

    }  
  }
  
  return converged;
}


/*
  lc_dynamics_r.cpp: C wrapper to C++ code for R access,
  simulation of one dimensional liquid crystal cells
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

// interface functions to allow C++/C  functions
// in dynamics.cpp etc to be called from R

#include "dynamics.hh"
#include <iostream>


#define NEWTON_TOL 1.0e-5
#define TIME_TOL 2.0e-3
#define NEWTON_ITER 12
#define MIN_TIMESTEP 1.0e-6

Parameters init_parameters( double *d, double *K, double *alpha, double *epsilon, double *flexo,
			    double *w, double *boundary_tilt,  double *wtw, double *boundary_twist, 
			    int  *nsp, double *rho0, double *mu, double *chi, 
			    double *charge, double *sigma0, double *ftime)
{

  Parameters p;
  
  p.d = *d;
  //p.pg = *pg;
  p.ftime = *ftime;

  for (int i = 0; i < 6; i++){
    p.alpha[i+1] = alpha[i];
  
  }
  
  if (p.alpha[5] == 0.0){
      //comapatibility kludge: before this code
      //did twist, it was usual to pass
      //alpha[4] = alpha[4] + alpha[5] and alpha[5] = 0
      //as only the sum appeared in the governing equations.
      p.alpha[5] = 0.5 * p.alpha[4];
      p.alpha[4] = 0.5 * p.alpha[4];
    }
  //Parodi relation
  p.alpha[6] = p.alpha[2] + p.alpha[3] + p.alpha[5];


  for (int i = 0; i < 3; i++){
    p.K[i+1] = K[i];
  } 
  
  p.eo = epsilon[0];
  p.ea = epsilon[1];
  p.flexo[1] = flexo[0] /E0 ;
  p.flexo[2] = flexo[1] /E0 ;
  
  for (size_t sp = 0; sp < (size_t)*nsp; sp++){
    p.add_ions(mu[sp],chi[sp],charge[sp],rho0[sp],&sigma0[2*sp]);
  }
  
  
  for (int i = 0; i < 2; i++){
    p.w[i] = w[i];
    p.wtw[i] = wtw[i];
    p.boundary_tilt[i] = boundary_tilt[i];
    p.boundary_twist[i] = boundary_twist[i];
    p.boundary_volt_ac[i] = 0;
    p.boundary_volt_dc[i] = 0;
    p.boundary_flowx[i] = 0.0;
    p.boundary_flowy[i] = 0.0;
  }

  return p;
};



extern "C" {

  void lcd_state ( int *nz, double *dt, double *volt, double *freq,
		   double *d, double *K, double *alpha, int *shear, 
		   double *epsilon, double *flexo, double *w, double *boundary_tilt,
		   double *wtw, double *boundary_twist,double *pgx, double *pgy,
		   int *nsp, double *rho0, double *mu, double *chi, double *charge, double *sigma0,
		   double *z, double *layer_depth, double *tilt, double *twist, 
		   double *volt_ac, double *volt_dc, double *flowx, double *flowy, double *rho, double *sigma,
		   int *newton_iter, double *newton_tol, double *min_timestep, 
		   double *time_tol, double *ftime,
		   double* residual, double* stepsize)
		   
  {
    
    // purpose: compute an lcd state given an lcd state as an initial estimate
    
    // accepts : nz - input - number of grid points to consider (including boundaries)
    //           dt - input - time difference between initial state and desired state
    //                        if dt == 0 a steady state is computed


    Geometry g(*nz,*d);
    Parameters p = init_parameters(d,K,alpha,epsilon,flexo, w, boundary_tilt, wtw, boundary_twist,
				   nsp, rho0, mu, chi, charge, sigma0, ftime);

    p.boundary_volt_ac[0] = *volt;
    p.pgx = *pgx;
    p.pgy = *pgy;


    // create the previous state
    State prev (*nz,&p,&g);
    
    for (size_t j = 0; j < (size_t)*nz; j++){
      prev.tilt[j] = tilt[j];
      prev.volt_ac[j] = volt_ac[j];
      prev.volt_dc[j] = volt_dc[j];
      prev.flowx[j] = flowx[j];
       
      for (size_t sp = 0; sp < (size_t)*nsp; sp++){
	size_t l = j + sp * (*nz);
	prev.rho[sp][j] = rho[l];
      }      
    }
    
    vector<Equation*> e;
    //TODO include more equations!
    e.push_back(new Volt_ACEquation());
    e.push_back(new TiltEquation());
    
    
    State *s = compute_state
      (e, &p, &g, &prev , *dt ,
       *newton_tol,  *newton_iter , *time_tol , *min_timestep);
    
    
    for (size_t j = 0; j < (size_t)*nz; j++){
      tilt[j] = s->tilt[j];
      volt_ac[j] = s->volt_ac[j];
      volt_dc[j] = s->volt_dc[j];
      flowx[j] = s->flowx[j];
      flowy[j] = s->flowy[j];
      for (size_t sp = 0; sp < (size_t)*nsp; sp++){
	size_t l = j + sp * (*nz);
	rho[l] = s->rho[sp][j];
      }      
    } 
    

    delete(s);

#define COPYR(x) if (x != 0) for (size_t i = 0; i < e.size(); i++) for (size_t j = 0; j < e.at(i)->x.size(); j++) x[i*(*newton_iter) + j] = e.at(i)->x[j];  

    //copy residuals 
    COPYR(residual);
    COPYR(stepsize);
      for (size_t i = 0; i < e.size(); i++){
	delete(e.at(i));
      }

  }


  void lcd_dynamics
  ( int *nz, int *nt, double *times, double *voltages, double *frequencies,
    double *d, double *K, double *alpha, int *shear, int* shear_y, double *epsilon, double *flexo,
    double *w, double *boundary_tilt, double *wtw, double *boundary_twist, 
    double *pgx, double *pgy, int *nsp, double *rho0, double *mu, 
    double *chi, double *charge, double *sigma0, double *z, double *layer_depth, double *tilt, 
    double *twist, double *volt_ac, double *volt_dc, double *flowx, double *flowy, double *rho, double *sigma,
    int *newton_iter, double *newton_tol, double *min_timestep, double *time_tol, double *ftime, int *iters)
  {
    // purpose : compute a (time-dependent) sequence of
    // lcd state given a set of lcd parameters and a 
    // sequence of {time,voltage,frequency} tuples

    // accepts:
    // nz - input - 
    //    number of grid points to consider (including boundaries)
    // offset - input
    //    offset into ouptut arrays that lcd starts at
    // nt - input - 
    //    number of time steps
    // times - input 
    //    array of times[nt]. if time[i] - time[i-1] = 0 then solution[i]
    //    will be a steady-state solution.     
    // voltages[nt],frequencies[nt] - input
    //    voltage and frequency sequences
    // z - output - z values at which state variables are computed
    // tilt,twist,volt_ac,volt_dc,flowx - input/output: state variables, array length[nz*nt]
    // tilt[i + j * nz] = tilt(z_i,t_j). On input, these are the initial guesses for steady
    // state solutions 


    vector<double> _times;
    vector<double> _voltages;
    vector<double> _frequencies;
    vector<double> _pgx;
    vector<double> _pgy;
  
    for (int i = 0; i < *nt; i++){
    
      _times.push_back(times[i]);
      _voltages.push_back(voltages[i]);
      _frequencies.push_back(frequencies[i]);
      _pgx.push_back(pgx[i]);
      _pgy.push_back(pgy[i]);


    }

    Geometry g(*nz,*d);
    Parameters p = init_parameters(d,K,alpha,epsilon,flexo, w, boundary_tilt, wtw, boundary_twist,
				   nsp, rho0, mu, chi, charge, sigma0, ftime);

    vector<State*> initial_states;
    for (int i = 0; i < *nt; i++){
      State *s = new State(*nz,&p,&g);
      
      //copy bulk properties (ion densities not yet supported)
      for (int j = 0; j < *nz; j++){
	size_t k = j + i*(*nz); 

	//lc fields
	s->tilt[j] = tilt[k];
	s->twist[j] = twist[k];
	s->volt_ac[j] = volt_ac[k] ;
	s->volt_dc[j] = volt_dc[k];
	s->flowx[j] = flowx[k];
	s->flowy[j] = flowy[k];
      }

      initial_states.push_back(s);
    }


    vector<State*> states = compute_states
      (_times, _voltages, _frequencies, _pgx, _pgy, initial_states,
       *shear > 0, *shear_y > 0,  &p, &g, *newton_tol, *newton_iter, *time_tol, *min_timestep);

    for (int j = 0 ; j < *nz; j++){
      z[j] = states[0]->geometry()->z_centre[j];
      layer_depth[j] = states[0]->geometry()->depth[j];
    }

    for (int i = 0; i < *nt; i++){

      //destroy initial states
      delete(initial_states[i]);

      //copy surface ion densities
      for (size_t sp = 0; sp < (size_t)*nsp; sp++){
	size_t j = 2 * i + sp * (*nt) * (2);
	sigma[j] = states[i]->sigma[sp][0];
	sigma[j+1] = states[i]->sigma[sp][1];
      }
      
      //copy bulk properties
      for (int j = 0; j < *nz; j++){
	size_t k = j + i*(*nz); 

	//lc fields
	tilt[k] = states[i]->tilt[j];
	twist[k] = states[i]->twist[j];
	volt_ac[k] = states[i]->volt_ac[j];
	volt_dc[k] = states[i]->volt_dc[j];
	flowx[k] = states[i]->flowx[j];
	flowy[k] = states[i]->flowy[j];

	for (size_t sp = 0; sp < (size_t)*nsp; sp++){
	  //copy bulk ion densities
	  size_t l = k + sp * (*nt) * (*nz);
	  rho[l] =  states[i]->rho[sp][j];
	}      
      }
    }

     

    for (int i = 0; i < *nt; i++){
      iters[i] = states[i]->iters;
      delete(states[i]);
    }

  }


}

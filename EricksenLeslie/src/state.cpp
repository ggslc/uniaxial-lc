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

#include "state.hh"
#include <math.h>

#include <iostream>


#define XMAX 1.797693e+308 
#define E0 8.854187817e-12


Parameters::Parameters(){
  for (int i = 0; i < 4; i++) K[i] = 1e-12; // elastic constants
  for (int i = 0; i < 7; i++) alpha[i] = 0; // leslie viscosities
  for (int i = 0; i < 3; i++) flexo[i] =  0; //flexo-electric coefficents
 
  eo = ea = 10; // dielectric constants
  


  boundary_tilt[0] =  boundary_tilt[1] = 0;//(easy) tilt at boundaries
  boundary_twist[0] =  boundary_twist[1] = 0;//(easy) twist at boundaries
  w[0] = w[1] = wtw[0] = wtw[1] = 0.0; // anchoring strength at boundaries
  
  boundary_volt_ac[0] = boundary_volt_ac[1] = 0 ; //voltage at boundaries
  boundary_volt_dc[0] = boundary_volt_dc[1] = 0 ; //voltage at boundaries

  boundary_flowx[0] = boundary_flowx[1] = 0; //flow at boundaries
  boundary_flowy[0] = boundary_flowy[1] = 0; //flow at boundaries

  d=0; //total cell depth, finite volume depth
  pgx = pgy =  0; // pressure gradient
}


Geometry::Geometry(size_t n_cells, double depth){
  _n_cells = n_cells;
  double dz = depth/double((n_cells-1));
  z_centre.resize(n_cells,0);
  this->depth.resize(n_cells,0);

  z_centre[0] = 0.0;
  this->depth[0] = 0.5*dz;
  for (size_t i = 1; i < n_cells; i++){
    z_centre[i] = z_centre[i-1] + dz;
    this->depth[i] = dz;
  }
  this->depth[n_cells - 1] = 0.5*dz;
}

void State::init(size_t n_lc_layers, Parameters* parameters,
	    Geometry* geometry)
{
  // initialisation common to steady state and time
  // dependent states
  _n_lc_layers = n_lc_layers;
  _parameters = parameters;
  _geometry = geometry;
  
}

State::State(size_t n_lc_layers, Parameters* parameters,
	     Geometry* geometry)
{
  // construct a steady state object with no previous state
  _time_dependent = false;
  _t = 0;

  if (parameters->ea < 0) tilt.resize(n_lc_layers,M_PI/2);
  else tilt.resize(n_lc_layers,0);

  if (parameters->boundary_tilt[1] > M_PI/2){
    for (size_t i = (size_t) (n_lc_layers/2); i < n_lc_layers; i++){
      tilt[i] = M_PI;
    }
  }

  twist.resize(n_lc_layers,parameters->boundary_twist[0]);
  volt_ac.resize(n_lc_layers,0);
  volt_dc.resize(n_lc_layers,0);
  //try a linear voltage...
  for (size_t i = 0; i < n_lc_layers; i++){
    volt_dc[i] = (parameters->boundary_volt_dc[0] * double(n_lc_layers - 1 - i)
		  + parameters->boundary_volt_dc[1] * double(i))/double(n_lc_layers - 1);
  
    volt_ac[i] = (parameters->boundary_volt_ac[0] * double(n_lc_layers - 1 - i)
		  + parameters->boundary_volt_ac[1] * double(i))/double(n_lc_layers - 1); 
  }


  flowx.resize(n_lc_layers,0);
  flowy.resize(n_lc_layers,0);
  //rho_p.resize(n_lc_layers,parameters->rho_p_0);
  //rho_m.resize(n_lc_layers,parameters->rho_m_0);

  // double cv = parameters->boundary_volt_dc[1] 
  //  - parameters->boundary_volt_dc[0];

  double efield = (- parameters->boundary_volt_dc[1] 
		   + parameters->boundary_volt_dc[0]);


  for (size_t sp = 0; sp < parameters->rho0.size(); sp++){
   
 
    vector<double> sg(2,parameters->sigma0[2 * sp]);
    sg[1] = parameters->sigma0[2 * sp+1];

    efield += pow(parameters->d,2) * parameters->charge[sp] * (sg[0] - sg[1])
     /(parameters->eo * E0);

    sigma.push_back(sg);

    


  }



  
#define BIG 1.0e+200
  for (size_t sp = 0; sp < parameters->rho0.size(); sp++){

    vector<double> r(n_lc_layers,parameters->rho0[sp]);
    if (true && (fabs(efield) > 1e-16) && (r[0] > 0)){

      double l = parameters->charge[sp] * efield * parameters->mu[sp]/parameters->chi[sp];
      double a = l * parameters->rho0[sp] / (min(exp(l),XMAX) - 1.0);

      for (size_t i = 0; i < n_lc_layers; i++){
      	double z = double(i)/double(n_lc_layers-1);  
      	r[i] = a * min(exp(l * z),XMAX);
      }

    }
    rho.push_back(r);
  }

  init(n_lc_layers, parameters, geometry);
}

State::State(State *s)
{
  // copy state from pointer to state
  *this = *s;
}


State::State(State* previous, double dt, Parameters* parameters,
	     Geometry* geometry){

  // construct a time valued state object from a previous state object

 
    _time_dependent = (dt >  0);
    _previous = previous;

    
  // copy flow and tilt values from previous state vectors into current
  // state vector, set voltage to zero (it is effectively always
  // steady state and low voltage states sit outside the radius
  // of convergence for high-voltage states)

  tilt.assign(previous->tilt.begin(),previous->tilt.end()); 

  //temp...
  //volt_ac.assign(previous->volt_ac.begin(),previous->volt_ac.end());
 
  flowx.assign(previous->flowx.begin(),previous->flowx.end());
  flowy.assign(previous->flowy.begin(),previous->flowy.end());
  twist.assign(previous->twist.begin(),previous->twist.end());

 //  for (int i = 0; i < twist.size(); i++){
//     if (fabs(tilt[i]) < 0.05)
//       if (twist[i] > M_PI)
// 	twist[i] = 2.0 * M_PI - twist[i];
//       else if (twist[i] <= -M_PI)
// 	twist[i] = - 2.0 * M_PI - twist[i];
//   }

  //twist.resize(previous->n_lc_layers(),0);



  //  rho_p.assign(previous->rho_p.begin(),previous->rho_p.end());
  // rho_m.assign(previous->rho_m.begin(),previous->rho_m.end());

  rho.assign(previous->rho.begin(),previous->rho.end());
  sigma.assign(previous->sigma.begin(),previous->sigma.end());

  //for (size_t sp = 0; sp < previous->rho.size(); sp++){
  //   vector<double> r = previous->rho[sp]
  //  rho[sp].assign(previous->rho[sp].begin(),previous->rho[sp].end());
  //}

  volt_ac.resize(previous->n_lc_layers(),0);
  volt_dc.resize(previous->n_lc_layers(),0);
  _t = previous->t() + dt;
    
 

  init(previous->n_lc_layers(), parameters, geometry);
  
 }


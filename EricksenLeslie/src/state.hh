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

#ifndef STATE_H
#define STATE_H

//define Geometry, Parameters, and State classes
//both the dynamics and optics calculations
//utilise these classes for data storage

#include <vector>


using namespace std;


class Geometry {
  // defines a (1D) regular grid. the intention
  // is to provide a basis for irregular grids at
  // some point...
  int _n_cells;
 public:
  vector<double> z_centre;
  vector<double> depth;
  Geometry(size_t n_cells, double depth); //constructor
  size_t n_cells(){return _n_cells;};
  
};


class Miescowicz
{
  // struct to store Miescowicz viscosities
 public:
  double gamma1, eta1, eta2, eta3, eta12;
  double alpha(size_t i);
};

inline double Miescowicz::alpha(size_t i){

  switch (i) 
    {
    case 1:
      return eta12;
    case 2:
      return (eta1-eta2-gamma1)/2.;
    case 3:
      return (eta1-eta2+gamma1)/2.; 
    case 4:
      return 2.0 * eta3;
    case 5:
      return (eta1 + 3*eta2 - 4*eta3 - gamma1)/2.;
    case 6:
      return (3*eta1 + eta2 - 4*eta3 - gamma1)/2.;
    default:
      return 0;
    }
}


class Parameters  {
  // parameters for the liquid crystal
 public: 
  double K[4]; // elastic constants
  double alpha[7]; // leslie viscosities
  double flexo[3]; //flexo-electric coefficents
  double eo,ea; // dielectric constants
  vector<double> mu,chi,charge,rho0,sigma0; 
  //ion mobility, diffusion, bulk charge carrier density, 
  //surface charge carrier density
  double boundary_tilt[2],boundary_twist[2]; //(easy) tilt,twist at boundaries
  double w[2],wtw[2]; // anchoring tilt,twist strength at boundaries
  double boundary_volt_ac[2],boundary_volt_dc[2] ; //voltage at boundaries
  double boundary_flowx[2],boundary_flowy[2]; //flow at boundaries
  double d; //total cell depth
  double pgx,pgy; //pressure gradient in x,y direction
  double ftime; // weighting factor for time dependence in equations  
  // 0 = explicit, 1 - implicit, 0.5 = crank-nicholson

  void set_miescowicz(Miescowicz m);
  Miescowicz get_miescowicz();

  Parameters();

  void add_ions(double _mu, double _chi, double _charge,
		double _rho,double *_sigma){
    mu.push_back(_mu);
    chi.push_back(_chi);
    charge.push_back(_charge);
    rho0.push_back(_rho);
    sigma0.push_back(_sigma[0]);
    sigma0.push_back(_sigma[1]);
  };


};


inline void Parameters::set_miescowicz
(Miescowicz m)
{
  for (size_t i = 1; i <=6 ; i++)
    {
      alpha[i] = m.alpha(i);
    }
}

inline Miescowicz Parameters::get_miescowicz()
{
  Miescowicz m;
  m.eta1 = 0.5 * (alpha[3] + alpha[4] + alpha[6]);
  m.eta2 = 0.5 * (-alpha[2] + alpha[4] + alpha[5]);
  m.eta3 = 0.5 * alpha[4];
  m.gamma1 = alpha[3] - alpha[2];
  m.eta12 = alpha[1];
  return m;
}



class State {  
 private:

  //private data
  size_t _n_lc_layers;
  size_t _left_lc_layer;
  double _t;
  bool _time_dependent;
  Geometry* _geometry;
  Parameters* _parameters;
  State* _previous; //pointer to previous state

  //private methods
  void init(size_t n_lc_layers, Parameters* parameters,
	    Geometry* geometry);


 public:
  //public data
  vector<double> tilt,twist,flowx,flowy,volt_ac,volt_dc,
    rho_m,rho_p; //state variables

  vector<vector<double> > rho,sigma;

  size_t iters; // no of iterations taken to compute state;
  

  //public methods
  State(size_t n_lc_layers, Parameters* parameters, 
	Geometry* geometry); //constructor
  State(State* previous, double dt, Parameters* parameters, 
	Geometry* geometry );
  State(State *s);
  size_t n_lc_layers(){return _n_lc_layers;};
  double t() {return _t;};
  double dt() {return(_time_dependent)?(t() - previous()->t()):0.0;};
  double dz() {return _parameters->d/double(_n_lc_layers);};

  Parameters* parameters() {return _parameters;};
  void set_parameters(Parameters* p){_parameters = p;};

  State* previous() {return _previous;};
  void set_previous(State* p){_previous = p;};

  Geometry* geometry(){return _geometry;};
  void set_geometry(Geometry* g){_geometry = g;};


  bool time_dependent(){return _time_dependent;};

 
 
};
#endif

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

#ifndef DYNAMICS_H
#define DYNAMICS_H
 
//stdlib includes
#include <vector>
#include <string>


#define E0 8.854187817e-12
#define BIG 1


#define VDB vector<double>
//#define VEC gsl_vector*
#include <valarray>
#define VEC valarray<double> *

#include "state.hh"

using namespace std;




class Equation {
  // abstract class specifies the structure of
  // equations 
 public:
  virtual double difference(State *a, State *b)=0;
  virtual void evaluate(VEC res, VEC ju, VEC jl, VEC jd, State* s)=0;
  virtual void update(VEC dx, State*s)=0;
  virtual string name()=0;

  vector<double> residual, stepsize;

};


class TiltEquation : public Equation {
  // specialisiation to tilt equation
  double steadyrhs(State *s,size_t i, int j);
 public:
  double difference(State *a, State *b);
  void evaluate(VEC res, VEC ju, VEC jl, VEC jd, State* s);
  void update(VEC dtilt, State*s);
  string name(){return "tilt equation";};
};


class TiltEquationXZ : public Equation {
  // specialisiation to tilt equation when director
  // and flow are confined to xz plane (ie no twist)
  double steadyrhs(State *s,size_t i, int j);
 public:
  double difference(State *a, State *b);
  void evaluate(VEC res, VEC ju, VEC jl, VEC jd, State* s);
  void update(VEC dtilt, State*s);
  string name(){return "tilt equation (no twist)";};
};


class TwistEquation : public Equation {
  // specialisiation to tilt equation
 public:
  double difference(State *a, State *b);
  void evaluate(VEC res, VEC ju, VEC jl, VEC jd, State* s);
  void update(VEC dtwist, State*s);
  string name(){return "twist equation";};
};

class FlowxEquation : public Equation {
  // specialisation for x component of Navier Stokes equation
 public:
  double difference(State *a, State *b);
  void evaluate(VEC res, VEC ju, VEC jl, VEC jd, State* s);
  void update(VEC dflowx, State*s);  
  double constant(int i,State* s);
  string name(){return "flow x equation";};
};

class FlowxEquationXZ : public Equation {
  // specialisation for x component of Navier Stokes equation
  // when all the director and flow fields are 
  // confined to the xz plane (ie no twist)
 public:
  double difference(State *a, State *b);
  void evaluate(VEC res, VEC ju, VEC jl, VEC jd, State* s);
  void update(VEC dflowx, State*s);  
  double constant(int i,State* s);
  string name(){return "flow x equation (no twist)";};
};


class FlowyEquation : public Equation {
  // specialisation for x component of Navier Stokes equation
 public:
  double difference(State *a, State *b);
  void evaluate(VEC res, VEC ju, VEC jl, VEC jd, State* s);
  void update(VEC dflowy, State*s);  
  string name(){return "flow y equation";};
};


class Volt_ACEquation :  public Equation {
  // specialisation for Poisson equation
 public:
  double difference(State *a, State *b);
  void evaluate(VEC res, VEC ju, VEC jl, VEC jd, State* s);
  void update(VEC dvolt_ac, State*s);  
  string name(){return "ac potential equation";};
};

class Volt_DCEquation :  public Equation {
  // specialisation for Poisson equation
 public:
  double difference(State *a, State *b);
  void evaluate(VEC res, VEC ju, VEC jl, VEC jd, State* s);
  void update(VEC dvolt_dc, State*s);  
  string name(){return "dc potential equation";};
};


class IonEquation : public Equation {
  // specialisation for ion transport equation
  size_t species;
 public:
  IonEquation(size_t sp){species = sp;};
  double difference(State *a, State *b);
  void evaluate(VEC res, VEC ju, VEC jl, VEC jd, State* s);
  void update(VEC drhoe, State*s);
   string name(){return "ion transport equation";};
};


class TriDiagWorkSpace{
 public:
  // gsl_vector *x,*r,*jd,*ju,*jl;
  valarray<double> *x,*r,*jd,*ju,*jl;
  TriDiagWorkSpace(int n);
  ~TriDiagWorkSpace();
};

class NewtonResult {
 public:
  NewtonResult();
  float residual, stepsize;
  size_t iters;
};

double norm2(VEC v);

NewtonResult tri_newton
(int offset,int n, Equation* eq, State* s, TriDiagWorkSpace* tdws) ;

bool semi_implicit_solve
(vector<Equation*> e, State* s,double tol, int niter, 
 TriDiagWorkSpace* tdws);

State* compute_state
(vector<Equation*> e, Parameters* p, Geometry* g, State* previous,
 double time, double newton_tol, int newton_iter, 
 double time_tol, double min_timestep);

vector<State*> compute_states 
(vector<double> const& times, vector<double>  const& voltages,  vector<double>  const& frequencies,
 vector<double>  const& pgx, vector<double>  const& pgy, vector<State*> const& initial_states,
 bool shear, bool shear_y, Parameters* p, Geometry* g, double newton_tol, int newton_iter, 
 double time_tol, double min_timestep);
#endif

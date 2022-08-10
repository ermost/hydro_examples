/*
 * =====================================================================================
 *
 *       Filename:  main.cc
 *
 *    Description:  1D Advection test problem
 *
 *        Version:  1.0
 *        Created:  19.06.2018 13:51:38
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Elias Roland Most (ERM), emost@itp.uni-frankfurt.de
 *   Organization:  Goethe University Frankfurt
 *
 * =====================================================================================
 */


#include "../../advance.hh"
#include "../../advect.hh"
#include "../../boundaries.hh"
#include "../../output.hh"
#include "../../reconstruct.hh"
#include "../../riemann.hh"
#include "../../storage.hh"
#include "./system.hh"
#include "./initial_data.hh"


static constexpr int this_ndim = 1;
const double global_cfl = 0.2;
const double global_final_time = 10;
const int global_ngz = 4; 
const int global_ext = 1024; 
const double global_dx = 1./global_ext; 
const double global_origin_x = 0.;
const int output_iteration = 4000000;

static constexpr bool HO = true; // High order scheme (4th) 

using this_system_t = SimpleAdvectionSystem<1,double>;

using this_grid_t = SimpleCartesianGrid<this_ndim>;
using this_storage = SimpleStorage<this_grid_t,this_system_t,double>;

using this_reconstruct_t = WenoZ_Reconstruct<false,double>;

using this_riemann_t = HLL_RiemannSolver<this_system_t>;
using this_hrsc_t = McCorquodale_FV<HO,this_system_t, double,this_grid_t,
      this_reconstruct_t, this_riemann_t>;


using initial_data_t = Advected_Wave1D;

static this_grid_t global_grid;


int main(){

  std::array<int,this_ndim> ngz {{global_ngz}};
  std::array<int,this_ndim> extent {{global_ext + 2*global_ngz}};
  std::array<double,this_ndim> dx{{global_dx}};
  std::array<double,this_ndim> origin{{global_origin_x - global_dx*global_ngz}};

  global_grid = this_grid_t(std::move(extent), std::move(ngz), std::move(dx), std::move(origin));


  this_storage U(global_grid);


  this_hrsc_t HRSC(U.grid);

  //Select a time stepper:

  auto tstepper = TimeStepperRK4<this_system_t,this_storage>(U.grid);
  tstepper.set_timestep_from_CFL(global_cfl);


  auto evolve_func = [&] (auto & U, auto &Uout) {

    //First boundaries
    
    PeriodicBC<this_ndim>(U);

    //Second fluxes

    HRSC.advect(U, Uout);

  };

  //setup initial data
  initial_data_t::initial_data(U);

  //Output initial data
  PeriodicBC<this_ndim>(U);
  output(U,tstepper.get_iteration());

  //Initial data is pointwise, now switch to volumes
  HRSC.switch_to_volume(U);
  PeriodicBC<this_ndim>(U);


  while(tstepper.get_current_time() < global_final_time){
    tstepper.advance(U,evolve_func);

    if(tstepper.get_iteration() % output_iteration == 0){
      PeriodicBC<this_ndim>(U);
      output(U,tstepper.get_iteration());
    }
  };

  //If not yet at right time, adjust!
  auto t_remain = global_final_time - tstepper.get_current_time();

  if(t_remain>0){
    tstepper.set_timestep(t_remain);
    tstepper.advance(U,evolve_func);
  };


  //Output final data
  PeriodicBC<this_ndim>(U);
  output(U,tstepper.get_iteration());

  return 0;

};





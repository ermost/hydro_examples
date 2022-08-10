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


#include "../../../advance.hh"
#include "../../../advect.hh"
#include "../../../boundaries.hh"
#include "../../../output.hh"
#include "../../../reconstruct.hh"
#include "../../../riemann.hh"
#include "../../../storage.hh"
#include "./../system.hh"
#include "./initial_data.hh"

#include<chrono>


static constexpr int this_ndim = 2;
const double global_cfl = 0.3;
const double global_final_time = 3.;
const int global_ngz = 4; 
const int global_ext_x = 128; 
const int global_ext_y = 256; 
const double global_dx = 1./global_ext_x; 
const double global_dy = 2./global_ext_y; 
const double global_origin_x = -0.5;
const double global_origin_y = -1.0;
const int output_iteration = 10000;


const int info_it = 10;

static constexpr bool HO = false; // High order scheme (4th) 

using this_system_t = SRHD<this_ndim,double,SimpleGammaLaw<double>>;

using this_grid_t = SimpleCartesianGrid<this_ndim>;
using this_storage = SimpleStorage<this_grid_t,this_system_t,double>;

using this_reconstruct_t = WenoZ_Reconstruct<false,double>;

using this_riemann_t = HLL_RiemannSolver<this_system_t,false>;
using this_hrsc_t = McCorquodale_FV<HO,this_system_t, double,this_grid_t,
      this_reconstruct_t, this_riemann_t>;


using initial_data_t = KelvinHelmholtz2D;

static this_grid_t global_grid;


int main(){

  SimpleGammaLaw<double>::Gamma =4./3.;

  std::array<int,this_ndim> ngz {{global_ngz, global_ngz}};
  std::array<int,this_ndim> extent {{global_ext_x + 2*global_ngz, 
    global_ext_y + 2*global_ngz}};
  std::array<double,this_ndim> dx{{global_dx, global_dy}};
  std::array<double,this_ndim> origin{{global_origin_x - global_dx*global_ngz,
    global_origin_y - global_dy*global_ngz }};

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

//    this_system_t::fill_aux(U);
    HRSC.advect(U, Uout);

  };

  //setup initial data
  initial_data_t::initial_data(U);

  //Output initial data
  PeriodicBC<this_ndim>(U);
  output(U,tstepper.get_iteration());

  this_system_t::switch_to_cons(U);

  //Initial data is pointwise, now switch to volumes
  HRSC.switch_to_volume(U);
  PeriodicBC<this_ndim>(U);

  std::chrono::time_point<std::chrono::high_resolution_clock> start, end;

    start = std::chrono::high_resolution_clock::now();
 


  while(tstepper.get_current_time() < global_final_time){
    tstepper.advance(U,evolve_func);

    if(tstepper.get_iteration() % info_it == 0){
      end = std::chrono::high_resolution_clock::now();

      double elapsed_seconds = std::chrono::duration_cast<std::chrono::milliseconds>
                             (end-start).count();

      auto cell_up = U.ndof*U.nsystem*info_it/elapsed_seconds;
      
      std::cout << "Current iteration: " << tstepper.get_iteration()
       << "\t" << "time: "<< tstepper.get_current_time() << "   " << "Cell up per second "<< cell_up*1e+3 <<"  elapsed time[s]: "<< elapsed_seconds*1.e-3 <<std::endl;

      start = std::chrono::high_resolution_clock::now();
    }

    if(tstepper.get_iteration() % output_iteration == 0){
      
      std::cout << "Iteration: " << tstepper.get_iteration() << "  Writing output" << std::endl;

      HRSC.switch_to_point(U);
      this_system_t::switch_to_prims(U);
      PeriodicBC<this_ndim>(U);
      output(U,tstepper.get_iteration());
      this_system_t::switch_to_cons(U);
      HRSC.switch_to_volume(U);
      PeriodicBC<this_ndim>(U);
    }
  };

  //If not yet at right time, adjust!
  auto t_remain = global_final_time - tstepper.get_current_time();

  if(t_remain>0){
    tstepper.set_timestep(t_remain);
    tstepper.advance(U,evolve_func);
  };


  //Output final data
  HRSC.switch_to_point(U);
  PeriodicBC<this_ndim>(U);
  output(U,tstepper.get_iteration());

  return 0;

};





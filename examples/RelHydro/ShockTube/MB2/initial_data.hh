/*
 * =====================================================================================
 *
 *       Filename:  initial_data.hh
 *
 *    Description:  Initial data
 *
 *        Version:  1.0
 *        Created:  19.06.2018 22:57:47
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Elias Roland Most (ERM), emost@itp.uni-frankfurt.de
 *   Organization:  Goethe University Frankfurt
 *
 * =====================================================================================
 */


#pragma once

struct SodsProblem{

  template<typename Tstorage>
  static inline void initial_data(Tstorage &U){
    
    // This explicitly assumes 1D for indices!
    #pragma omp parallel for
    for(int ijk = 0; ijk < U.ndof; ++ijk){
      auto const x = U.grid.get_coords(0,ijk);

      if(x < 0.){
      	U[ijk + U.ndof*0] = 1.;
      	U[ijk + U.ndof*1] = 10.;
      	U[ijk + U.ndof*2] = -0.6;
      	U[ijk + U.ndof*3] = 0.;
      	U[ijk + U.ndof*4] = 0.;
      }
      else{
      	U[ijk + U.ndof*0] = 10.;
      	U[ijk + U.ndof*1] = 20.;
      	U[ijk + U.ndof*2] = 0.5;
      	U[ijk + U.ndof*3] = 0.0;
      	U[ijk + U.ndof*4] = 0.;
      }

      auto const v2 = U[ijk + U.ndof*2]*U[ijk + U.ndof*2];
      auto const lorentz_sq = 1./(1.-v2);
      auto const lorentz = std::sqrt(lorentz_sq);

      U[ijk + U.ndof*2] *= lorentz;

      	U[ijk + U.ndof*1] /= SimpleGammaLaw<double>::Gamma - 1.;
    }
  };
};


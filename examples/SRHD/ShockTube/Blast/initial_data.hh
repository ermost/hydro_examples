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

struct BlastProblem{

  template<typename Tstorage>
  static inline void initial_data(Tstorage &U){
    
    // This explicitly assumes 1D for indices!
    #pragma omp parallel for
    for(int ijk = 0; ijk < U.ndof; ++ijk){
      auto const x = U.grid.get_coords(0,ijk);

      if(x < 0.){
      	U[ijk + U.ndof*0] = 1.;
      	U[ijk + U.ndof*1] = 1000.;
      	U[ijk + U.ndof*2] = 0.;
      	U[ijk + U.ndof*3] = 0.;
      	U[ijk + U.ndof*4] = 0.;
      }
      else{
      	U[ijk + U.ndof*0] = 1.;
      	U[ijk + U.ndof*1] = 0.01;
      	U[ijk + U.ndof*2] = 0.;
      	U[ijk + U.ndof*3] = 0.;
      	U[ijk + U.ndof*4] = 0.;
      }

      	U[ijk + U.ndof*1] /= SimpleGammaLaw<double>::Gamma - 1.;
    }
  };
};


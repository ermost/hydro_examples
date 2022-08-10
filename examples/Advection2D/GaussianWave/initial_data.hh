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

struct Advected_Wave2D{

  static constexpr double x0 = 0.5;
  static constexpr double y0 = 0.5;

  template<typename Tstorage>
  static inline void initial_data(Tstorage &U){
    
    // This explicitly assumes 1D for indices!
    #pragma omp parallel for
    for(int j = 0; j < U.grid.extent[1]; ++j)
#pragma omp simd
    for(int i = 0; i < U.grid.extent[0]; ++i){

      int ijk = i+ U.grid.extent[0]*j;

      auto const x = U.grid.get_coords(0,i) -x0;
      auto const y = U.grid.get_coords(1,j) -y0;

    //Implement u = exp (-256 r^2)

//      U[ijk] = ...

    }
  };
};



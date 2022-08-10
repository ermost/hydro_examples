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

struct Burgers_Cosine1D{

  template<typename Tstorage>
  static inline void initial_data(Tstorage &U){
    //Implement u = cos(2 pi x)
    
    // This explicitly assumes 1D for indices!
    #pragma omp parallel for simd
    for(int ijk = 0; ijk < U.ndof; ++ijk){
      auto const x = U.grid.get_coords(0,ijk);

      // Implement initial data
      ///U[ijk] = ...;

    }
  };
};


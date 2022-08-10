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


struct Slotted_Cylinder2D{

  static constexpr double x0 = 0.5;
  static constexpr double y0 = 0.75;

  static constexpr double W = 0.05;
  static constexpr double H = 0.25;
  static constexpr double R = 0.15;



  template<typename Tstorage>
  static inline void initial_data(Tstorage &U){
    //Implement u = exp (-2 cos(2 pi x))
    
    // This explicitly assumes 1D for indices!
    #pragma omp parallel for
    for(int j = 0; j < U.grid.extent[1]; ++j)
#pragma omp simd
    for(int i = 0; i < U.grid.extent[0]; ++i){

      int ijk = i+ U.grid.extent[0]*j;

      auto const x = U.grid.get_coords(0,i) -x0;
      auto const y = U.grid.get_coords(1,j) -y0;

      auto const r = std::sqrt(x*x + y*y);

      U[ijk] = 1;

      if(r > R) U[ijk] = 0;

      if((std::abs(x) < 0.5*W) && (y +R >0) && (y+R <H))
	U[ijk] = 0; 

    }
  };
};



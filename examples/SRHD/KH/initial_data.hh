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

struct KelvinHelmholtz2D{

  static constexpr double A = 0.1;
  static constexpr double a = 0.01;
  static constexpr double sigma = 0.1;
  static constexpr double Vshear = 0.5;
  static constexpr double rho0 = 0.505;
  static constexpr double rho1 = 0.495;


  template<typename Tstorage>
  static inline void initial_data(Tstorage &U){
    
    // This explicitly assumes 1D for indices!
    #pragma omp parallel for
    for(int j = 0; j < U.grid.extent[1]; ++j)
    for(int i = 0; i < U.grid.extent[0]; ++i){

      int ijk = i+ U.grid.extent[0]*j;

      auto const x = U.grid.get_coords(0,i);
      auto const y = U.grid.get_coords(1,j);

      if(y <= 0.){
      	U[ijk + U.ndof*0] = rho0 - rho1*std::tanh((y+0.5)/a);
      	U[ijk + U.ndof*1] = 1.;
      	U[ijk + U.ndof*2] = -Vshear*std::tanh((y+0.5)/a);
      	U[ijk + U.ndof*3] = -A*Vshear*std::sin(2.*M_PI*x)*std::exp(-(y+0.5)*(y+0.5)/sigma/sigma);
      	U[ijk + U.ndof*4] = 0.;
      }
      else{
      	U[ijk + U.ndof*0] = rho0 + rho1*std::tanh((y-0.5)/a);
      	U[ijk + U.ndof*1] = 1.;
      	U[ijk + U.ndof*2] = Vshear*std::tanh((y-0.5)/a);
      	U[ijk + U.ndof*3] = A*Vshear*std::sin(2.*M_PI*x)*std::exp(-(y-0.5)*(y-0.5)/sigma/sigma);
      	U[ijk + U.ndof*4] = 0.;
      }

      auto const v2 =      U[ijk + U.ndof*2]*U[ijk + U.ndof*2]
			  +U[ijk + U.ndof*3]*U[ijk + U.ndof*3]
			  +U[ijk + U.ndof*4]*U[ijk + U.ndof*4];

      auto const lorentz = std::sqrt(1./(1.-v2));

      U[ijk + U.ndof*2] *= lorentz;
      U[ijk + U.ndof*3] *= lorentz;
      U[ijk + U.ndof*4] *= lorentz;

      	U[ijk + U.ndof*1] /= SimpleGammaLaw<double>::Gamma - 1.;

//	if(std::abs(y) < 0.1) std::cout  << U[ijk + U.ndof*0] << std::endl;
    }
  };
};


/*
 * =====================================================================================
 *
 *       Filename:  system.hh
 *
 *    Description:  Systems of advection type equations
 *
 *        Version:  1.0
 *        Created:  11.06.2018 17:58:25
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Elias Roland Most (ERM), emost@itp.uni-frankfurt.de
 *   Organization:  Goethe University Frankfurt
 *
 * =====================================================================================
 */

template <int Tndim = 1, typename T = double> class SimpleAdvectionSystem {

public:
  static constexpr bool needs_c2p = false;
  static constexpr int ndim = Tndim;
  static constexpr int num_vars = 1;
  static constexpr int num_aux = 2;

  static constexpr T uL = 2.*M_PI;


  template <int dir,typename Td>
  static inline decltype(auto) compute_flux(std::array<Td, num_vars + num_aux> &U) {

    // Compute fluxes for all components, here we have only one!

//    for(int nv=0; nv <  num_vars + num_aux; ++nv)
//      std::cout << U[nv] << "\t";
//    std::cout << std::endl;
    return std::array<Td, num_vars>{{U[0]*U[num_vars+dir]}};
  };

  template <typename Tstorage>
  static inline void add_source(Tstorage &storage, Tstorage &source) {

    // This system has no source
    //      std::memset(source,    0, Tstorage::ndof*sizeof(T))

    return;
  };

  template <typename Tstorage>
  static inline void fill_aux(Tstorage &U) {


#pragma omp parallel for
      for (int j = 0; j < U.grid.extent[1]; ++j){
	auto const y = U.grid.get_coords(1,j) - 0.5*U.grid.dx[1];
#pragma omp simd
        for (int i = 0; i < U.grid.extent[0]; ++i) {
	  auto const x = U.grid.get_coords(0,i) - 0.5*U.grid.dx[0];

            U.aux[i + U.grid.extent[0] * (j + U.grid.extent[1] * 0)]
	      = - ( y- 0.5)*uL;
            U.aux[i + U.grid.extent[0] * (j + U.grid.extent[1] * 1)]
	      =  ( x- 0.5)*uL;

        } // for i
      }// for j

    return;
  };

  template <typename Tstorage>
  static inline void switch_to_cons(Tstorage &storage) {
    storage.is_primitive = false;
    return; // No C2P needed
  };

  template <typename Td>
  static inline void switch_to_cons_single(std::array<Td, num_vars+num_aux> &U) {
    return; // No C2P needed
  };

  template <typename Td>
  static inline void switch_to_prims_single(std::array<Td, num_vars+num_aux> &U) {
    return; // No C2P needed
  };

  template <typename Tstorage>
  static inline void switch_to_prims(Tstorage &storage) {
    storage.is_primitive = false;
    return; // No C2P needed
  };

  template <int dir, typename Td>
  static inline decltype(auto)
  compute_max_characteristics(std::array<Td, num_vars+num_aux> &U) {
    auto a = std::abs(U[num_vars+dir]);
    return std::array<T, 2>{{ a, -a  }};
  };
};

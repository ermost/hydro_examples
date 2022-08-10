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
  static constexpr int num_aux = 0;

  static constexpr T uL = 1.; // CONSTANT ADVECTION SPEED

  template <typename Tstorage, typename Fstorage>
  static inline void compute_flux(Tstorage &storage, Fstorage &flux) {

    // Compute fluxes for all components, here we have only one!

#pragma omp simd
    for (int ijk = 0; ijk < Tstorage::ndof; ++ijk) {
      // IMPLEMENT FLUX FORMULA
      //flux.U[ijk] = 
      
    };

    return;
  };

  template <int dir,typename Td>
  static inline decltype(auto) compute_flux(std::array<Td, num_vars> &U) {

    // Compute fluxes for all components, here we have only one!

    return std::array<Td, num_vars>{{
    	 // IMPLEMENT FLUX FOR ALL COMPONENTS (just one)
    	
    }};
  };

  template <typename Tstorage>
  static inline void add_source(Tstorage &storage, Tstorage &source) {

    // This system has no source
    //      std::memset(source,    0, Tstorage::ndof*sizeof(T))

    return;
  };

  template <typename Tstorage>
  static inline void fill_aux(Tstorage &storage) {

    // This system has no source
    //      std::memset(source,    0, Tstorage::ndof*sizeof(T))

    return;
  };

  template <typename Tstorage>
  static inline void switch_to_cons(Tstorage &storage) {
    storage.is_primitive = false;
    return; // No C2P needed
  };

  template <typename Td>
  static inline void switch_to_cons_single(std::array<Td, num_vars> &U) {
    return; // No C2P needed
  };

  template <typename Td>
  static inline void switch_to_prims_single(std::array<Td, num_vars> &U) {
    return; // No C2P needed
  };

  template <typename Tstorage>
  static inline void switch_to_prims(Tstorage &storage) {
    storage.is_primitive = false;
    return; // No C2P needed
  };

  template <int dir,typename Td>
  static inline decltype(auto)
  compute_max_characteristics(std::array<Td, num_vars> &U) {
    // IMPLEMENT CHARACTERISTICS:
    // +-c 
    return std::array<T, 2>{{
       // +-c ..
    }};
  };
};

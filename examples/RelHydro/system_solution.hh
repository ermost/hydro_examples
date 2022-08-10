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

#pragma once

#include "./eos.hh"

template <int Tndim = 1 , typename T = double, typename EOS_t=SimpleGammaLaw<T>> class SRHD {

public:
  // PRIMITVE VARIABLES ARE:
  // \rho, \varepsilon, u_i
  
  enum { RHOB = 0, EPS, WVX, WVY, WVZ, NUM_VARS };
  enum { RHOSTAR = 0, TAUENERGY, STX, STY, STZ };

  static constexpr bool needs_c2p = true;
  static constexpr int ndim = Tndim;
  static constexpr int num_vars = NUM_VARS;
  static constexpr int num_aux = 0;
  static T rho_atmo;
  static T eps_atmo;
  static constexpr double c2p_tol = 1.e-10;

  template <typename Tstorage, typename Fstorage>
  static inline void compute_flux(Tstorage &storage, Fstorage &flux) {

    // Compute fluxes for all components, here we have only one!
    assert(!"Not implemented");
    return;
  };

  template <int dir, typename Td>
  static inline decltype(auto) compute_flux(std::array<Td, num_vars> &U) {

    // Compute fluxes for all components, here we have only one!

    auto const z2 = U[WVX] * U[WVX] + U[WVY] * U[WVY] + U[WVZ] * U[WVZ];
    auto const lorentz = std::sqrt(1. + z2);

    // Need EOS call here!!
    T eps_tot;
    auto press = EOS_t::press_eps__rho_eps_th(eps_tot, U[RHOB], U[EPS]);

    auto const rhoh = U[RHOB] + (press + eps_tot);

    auto const stx = rhoh * U[WVX];
    auto const sty = rhoh * U[WVY];
    auto const stz = rhoh * U[WVZ];

//    auto const taud = rhoh*lorentz - U[RHOB];
//    This should also work without double precision for small velocities
    auto const taud = (press + eps_tot) * lorentz + U[RHOB]*z2/(lorentz+1.);

    auto result = std::array<Td, num_vars>{
        {U[RHOB] * U[WVX + dir], (taud)*U[WVX + dir], stx * U[WVX + dir],
         sty * U[WVX + dir], stz * U[WVX + dir]}};

    result[STX + dir] += press;

    return result;
  };

  template <typename Tstorage>
  static inline void add_source(Tstorage &storage, Tstorage &source) {

    // This system has no source
    //      std::memset(source,    0, Tstorage::ndof*sizeof(T))

    return;
  };

  template <typename Tstorage> static inline void fill_aux(Tstorage &UL) {

    return;
  };

  template <typename Td>
  static inline void switch_to_cons_single(std::array<Td, num_vars> &U) {

    auto const z2 = U[WVX] * U[WVX] + U[WVY] * U[WVY] + U[WVZ] * U[WVZ];
    auto const lorentz2 = 1. + z2;
    auto const lorentz = std::sqrt(lorentz2);

    // Need EOS call here!!
    T eps_tot;
    auto press = EOS_t::press_eps__rho_eps_th(eps_tot, U[RHOB], U[EPS]);

    U[RHOSTAR] *= lorentz;
    auto const rhohW = U[RHOSTAR] + (press + eps_tot) * lorentz;

    U[STX] *= rhohW;
    U[STY] *= rhohW;
    U[STZ] *= rhohW;

//    U[TAUENERGY] = rhohW * lorentz - press -  U[RHOSTAR];
//    This should also work without double precision for small velocities
   U[TAUENERGY]  = press*z2 + eps_tot * lorentz2 + U[RHOSTAR]*z2/(lorentz+1.);

    return;
  };

  template <typename Td>
  static inline void switch_to_prims_single(std::array<Td, num_vars> &U) {

//    for(auto & x : U) std::cout << x << " ;  ";
//    std::cout << std::endl;

    if (U[RHOSTAR] < rho_atmo) {
      U[RHOB] = rho_atmo;
      U[EPS] = eps_atmo;
      U[WVX] = 0;
      U[WVY] = 0;
      U[WVZ] = 0;
      return;
    }

//    T press = 0.;
    T epst;
    T press = EOS_t::press_eps__rho_eps_th(epst, U[RHOB], U[EPS]);
    T press_prev = 1.e99;
    T press_pp = 1.e109;

    std::array<Td, num_vars> Uout;

    auto const rhostari = 1./U[RHOSTAR];

    auto const Snorm = std::sqrt(U[STX] * U[STX] + U[STY] * U[STY]+ U[STZ] * U[STZ])*rhostari;

    int nn = 0;
    T lorentz = 1.;
    T hWi =  1.;

    bool print = false;
    while ((std::abs(press - press_prev) > c2p_tol * press) && (nn < 1000)) {
//      if(nn>500 && print== false){  print=true; nn=0; press= 0; press_prev = 1.e99;};



      auto const hW = (U[TAUENERGY] + press) * rhostari + 1.;
      hWi = 1. / hW;


      auto const v = std::min( std::abs(Snorm*hWi), 0.9999499987499375);
      auto const v2 = v*v;

      auto const lorentz2i = std::abs(1.-(v2));
      auto const lorentzi = std::sqrt(lorentz2i);
      lorentz = 1./lorentzi;


      Uout[RHOB] = U[RHOSTAR] * lorentzi;

//      auto eps = Uout[RHOB]*(hW * lorentzi - 1.) - press;
      auto eps = U[TAUENERGY] * lorentz2i - (U[RHOSTAR] /(lorentz+1.) + press)*v2;
      eps =  std::max(eps, Uout[RHOB]*eps_atmo);
      Uout[EPS] = EOS_t::eps_th__eps_rho(eps, Uout[RHOB]);


      press_pp = press_prev;
      press_prev = press;

      // Need EOS call here!!
      press = EOS_t::press_eps__rho_eps_th(eps, Uout[RHOB], Uout[EPS]);
      press = std::abs(press); //std::max(-(U[RHOSTAR]+U[TAUENERGY]), press);


      //Aitken accelerator
      auto const R = (press - press_prev)/(press_prev - press_pp);
      auto const Paitken = std::abs(press_prev + (press-press_prev)/(1.-R));

     if(std::abs(R) < 1. && nn > 2){
	      press_pp = press_prev;
	      press_prev = press;
	      press=Paitken;
     };


      nn++;
    };

    assert(nn < 1000); // Abort if C2P fails!

     auto const conv =  hWi * rhostari * lorentz;

      U[RHOB] =Uout[RHOB];
      U[EPS] = Uout[EPS] ;
      U[WVX] = U[STX] * conv;
      U[WVY] = U[STY] * conv;
      U[WVZ] = U[STZ] * conv;


    return; 
  };

  template <typename Tstorage> static inline void switch_to_cons(Tstorage &UL) {
    UL.is_primitive = false;
    if (ndim == 1) {
      //#pragma omp simd
      for (int i = 0; i < UL.grid.extent[0]; ++i) {

        std::array<typename Tstorage::data_t,
                   num_vars + num_aux>
            Ul;

        for (int nv = 0; nv < num_vars; ++nv) {
          Ul[nv] = UL[i + UL.grid.extent[0] * (nv)];
        }

        for (int nv = 0; nv < num_aux; ++nv) {
          Ul[nv + num_vars] = UL.aux[i + UL.grid.extent[0] * (nv)];
        }

        switch_to_cons_single(Ul);

        for (int nv = 0; nv < num_vars; ++nv) {
          UL[i + UL.grid.extent[0] * (nv)] = Ul[nv];
        }
      } // for i
    }

    if (ndim == 2) {
#pragma omp parallel for
      for (int j = 0; j < UL.grid.extent[1]; ++j)
        #pragma omp simd
        for (int i = 0; i < UL.grid.extent[0]; ++i) {

          std::array<typename Tstorage::data_t,
                     num_vars + num_aux>
              Ul;

          for (int nv = 0; nv < num_vars; ++nv) {
            Ul[nv] = UL[i + UL.grid.extent[0] * (j + UL.grid.extent[1] * nv)];
          }

          for (int nv = 0; nv < num_aux; ++nv) {
            Ul[num_vars + nv] =
                UL.aux[i + UL.grid.extent[0] * (j + UL.grid.extent[1] * nv)];
          }

          switch_to_cons_single(Ul);

          for (int nv = 0; nv < num_vars; ++nv) {
            UL[i + UL.grid.extent[0] * (j + UL.grid.extent[1] * nv)] = Ul[nv];
          }

        } // for i
    };

    if (ndim == 3) {
      assert(!"Not implemented yet");
    };
    return; // No C2P needed
  };

  template <typename Tstorage>
  static inline void switch_to_prims(Tstorage &UL) {
    UL.is_primitive = true;
    if (ndim == 1) {
      //#pragma omp simd
      for (int i = 0; i < UL.grid.extent[0]; ++i) {

        std::array<typename Tstorage::data_t,
                   num_vars + num_aux>
            Ul;

        for (int nv = 0; nv < num_vars; ++nv) {
          Ul[nv] = UL[i + UL.grid.extent[0] * (nv)];
        }

        for (int nv = 0; nv < num_aux; ++nv) {
          Ul[nv + num_vars] = UL.aux[i + UL.grid.extent[0] * (nv)];
        }

        switch_to_prims_single(Ul);

        for (int nv = 0; nv < num_vars; ++nv) {
          UL[i + UL.grid.extent[0] * (nv)] = Ul[nv];
        }
      } // for i
    }

    if (ndim == 2) {
#pragma omp parallel for
      for (int j = 0; j < UL.grid.extent[1]; ++j)
        //#pragma omp simd
        for (int i = 0; i < UL.grid.extent[0]; ++i) {

          std::array<typename Tstorage::data_t,
                     num_vars + num_aux>
              Ul;

          for (int nv = 0; nv < num_vars; ++nv) {
            Ul[nv] = UL[i + UL.grid.extent[0] * (j + UL.grid.extent[1] * nv)];
          }

          for (int nv = 0; nv < num_aux; ++nv) {
            Ul[num_vars + nv] =
                UL.aux[i + UL.grid.extent[0] * (j + UL.grid.extent[1] * nv)];
          }

          switch_to_prims_single(Ul);

          for (int nv = 0; nv < num_vars; ++nv) {
            UL[i + UL.grid.extent[0] * (j + UL.grid.extent[1] * nv)] = Ul[nv];
          }

        } // for i
    };

    if (ndim == 3) {
      assert(!"Not implemented yet");
    };
    return; // No C2P needed
  };

  template <int dir, typename Td>
  static inline decltype(auto)
  compute_max_characteristics(std::array<Td, num_vars> &U) {

    auto const z2 = U[WVX] * U[WVX] + U[WVY] * U[WVY] + U[WVZ] * U[WVZ];
    auto const lorentz2 = 1. + z2;
    auto const lorentz = std::sqrt(lorentz2);
    auto const lorentzi = 1./lorentz;
    auto const lorentzi2 = 1./lorentz2;

    auto const v2 = z2 * lorentzi2;

    // Need EOS call here!!
    auto const cs2 = EOS_t::cs2__rho_eps_th(U[RHOB], U[EPS]);

    auto const z2_par = U[WVX + dir] * U[WVX + dir];

    auto const tmp = std::sqrt( 
	cs2 * lorentzi2 * (1. - (z2_par + (z2 - z2_par) * cs2) *lorentzi2)
	);

//    auto const tmp = lorentzi2*std::sqrt( 
//	cs2 * (1. +(z2 - z2_par) * (1.-cs2)));


    auto const p1 = U[WVX + dir] * lorentzi * (1. - cs2);

    auto const invden =1./ (1. - v2 * cs2);

    return std::array<T, 2>{{(p1 + tmp) * invden, ((p1 - tmp) * invden)}};
  };
};

template <int Tndim, typename T, typename EOS_t> T SRHD<Tndim,T,EOS_t>::rho_atmo = 1.e-10;

template <int Tndim, typename T,typename EOS_t> T SRHD<Tndim,T,EOS_t>::eps_atmo = 1.e-10;

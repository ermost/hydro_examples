/*
 * =====================================================================================
 *
 *       Filename:  riemann.hh
 *
 *    Description:  Riemann solvers
 *
 *        Version:  1.0
 *        Created:  12.06.2018 21:49:11
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Elias Roland Most (ERM), emost@itp.uni-frankfurt.de
 *   Organization:  Goethe University Frankfurt
 *
 * =====================================================================================
 */

#pragma once
#include <algorithm>
#include <cmath>
#include <functional>
#include <iostream>

template <typename Tsystem, bool LLF = false, bool pp=false> class HLL_RiemannSolver {

public:
  template <int dir, typename Tstorage>
  static inline void solve(Tstorage &UL, Tstorage &UR, Tstorage &flux) {

    auto hll_solve_single = [&](auto &Ul, auto &Ur) {

      auto const Fl = Tsystem::template compute_flux<dir>(Ul);
      auto const Fr = Tsystem::template compute_flux<dir>(Ur);


      typename Tstorage::data_t cp =1.;
      typename Tstorage::data_t cm =1.;

      if(!pp){
	      auto  cL = Tsystem::template compute_max_characteristics<dir>(Ul);
	      auto  cR = Tsystem::template compute_max_characteristics<dir>(Ur);

	      cp = std::max(std::max(cL[0], cR[0]), 0.);
	      cm = std::max(std::max(-cL[1], -cR[1]), 0.);
      };

      if (LLF) {
        cp = std::max(cp, cm);
        cm = cp;
      }

      Tsystem::switch_to_cons_single(Ul);
      Tsystem::switch_to_cons_single(Ur);

      std::array<typename Tstorage::data_t, Tsystem::num_vars> tot_flux;

      if(std::max(cp,cm) < 1.e-4) {cp =1.; cm=1;};

      auto ctot = 1./(cp+cm);

	for (int nv = 0; nv < Tsystem::num_vars; ++nv) {
	  tot_flux[nv] = 
	      (cm * Fr[nv] + cp * Fl[nv] - cp * cm * (Ur[nv] - Ul[nv]))*ctot;
	}

      cp = cm = 1.;
      ctot = 0.5;

//      tot_flux[Tsystem::PHI] = (cm * Fr[Tsystem::PHI] + cp * Fl[Tsystem::PHI] - cp * cm * (Ur[Tsystem::PHI] - Ul[Tsystem::PHI]))*ctot;
//      tot_flux[Tsystem::BX+dir] = (cm * Fr[Tsystem::BX+dir] + cp * Fl[Tsystem::BX+dir] - cp * cm * (Ur[Tsystem::BX+dir] - Ul[Tsystem::BX+dir]))*ctot;

      return tot_flux;
    };

    if (Tsystem::ndim == 1) {
//#pragma omp simd
      for (int i = 0; i < UL.grid.extent[0]; ++i) {

        std::array<typename Tstorage::data_t, Tsystem::num_vars+Tsystem::num_aux> Ul, Ur;

        for (int nv = 0; nv < Tsystem::num_vars; ++nv) {
          Ul[nv] = UL[i + UL.grid.extent[0] * (nv)];
          Ur[nv] = UR[i + UR.grid.extent[0] * (nv)];
        }

        for (int nv = 0; nv < Tsystem::num_aux; ++nv) {
          Ul[nv+Tsystem::num_vars] = UL.aux[i + UL.grid.extent[0] * (nv)];
          Ur[nv+Tsystem::num_vars] = UR.aux[i + UR.grid.extent[0] * (nv)];
        }

        auto tot_flux = hll_solve_single(Ul, Ur);

        for (int nv = 0; nv < Tsystem::num_vars; ++nv) {
          flux[i + UL.grid.extent[0] * (nv)] = tot_flux[nv];
        }
      } // for i
    }

    if (Tsystem::ndim == 2) {
#pragma omp parallel for
      for (int j = 0; j < UL.grid.extent[1]; ++j)
//#pragma omp simd
        for (int i = 0; i < UL.grid.extent[0]; ++i) {

          std::array<typename Tstorage::data_t, Tsystem::num_vars+Tsystem::num_aux> Ul, Ur;

          for (int nv = 0; nv < Tsystem::num_vars; ++nv) {
            Ul[nv] = UL[i + UL.grid.extent[0] * (j + UL.grid.extent[1] * nv)];
            Ur[nv] = UR[i + UR.grid.extent[0] * (j + UR.grid.extent[1] * nv)];
          }

          for (int nv = 0; nv < Tsystem::num_aux; ++nv) {
            Ul[Tsystem::num_vars+nv] = UL.aux[i + UL.grid.extent[0] * (j + UL.grid.extent[1] * nv)];
            Ur[Tsystem::num_vars+nv] = UR.aux[i + UR.grid.extent[0] * (j + UR.grid.extent[1] * nv)];
          }


          auto tot_flux = hll_solve_single(Ul, Ur);

          for (int nv = 0; nv < Tsystem::num_vars; ++nv) {
            flux[i + UL.grid.extent[0] * (j + UL.grid.extent[1] * nv)] =
                tot_flux[nv];
          }

        } // for i
    };

    if (Tsystem::ndim == 3) {
      assert(!"Not implemented yet");
    };

    flux.is_primitive = false;
  }
};

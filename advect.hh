/*
 * =====================================================================================
 *
 *       Filename:  advect.hh
 *
 *    Description:  Simple advection routine
 *
 *        Version:  1.0
 *        Created:  11.06.2018 19:18:06
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Elias Roland Most (ERM), emost@itp.uni-frankfurt.de
 *   Organization:  Goethe University Frankfurt
 *
 * =====================================================================================
 */

#pragma once

#include "storage.hh"
#include <algorithm>
#include <array>
#include <cassert>
#include <cmath>
#include <cstring>
#include <functional>
#include <memory>

template <bool high_order, typename Tsystem, typename T, typename grid_t,
          typename Treconstruct, typename TriemannSolver>
class McCorquodale_FV {

private:
  using storage_t = SimpleStorage<grid_t, Tsystem>;

  storage_t scratch, scratchL, scratchR, scratchL2, scratchR2, scratchF, scratchC, scratchP;
     // ,final_flux;

  using TR = Treconstruct;
  using TS = Tsystem;
  using TRiem = TriemannSolver;

  grid_t &grid;

  template <int dir, typename Tstorage>
  inline void calculate_flux(Tstorage &U) {
    static_assert(std::is_same<grid_t, typename Tstorage::grid_t>::value,
                  "Grids don't match");

    // Actual reconstruction step
    TR::template reconstruct<TS::ndim, dir>(U, scratchR, scratchL);
    //Need to fill aux
    TS::fill_aux(scratchR);
    TS::fill_aux(scratchL);

    if (TS::ndim == 2 && high_order) {

      auto correct_rec = [&](auto sL, auto sO) {
        std::memset(scratch.U, 0,
                    U.ndof * Tstorage::nsystem *
                        sizeof(typename Tstorage::data_t));

        // Need to correct with the transverse laplacian
        if (dir == 0) {
          for (int nv = 0; nv < Tsystem::num_vars; ++nv)
#pragma omp parallel for
            for (int j = 1; j < U.grid.extent[1] - 1; ++j)
              for (int i = 0; i < U.grid.extent[0]; ++i) {
                const int ij =
                    i + U.grid.extent[0] * (j + U.grid.extent[1] * nv);
                const int ijp =
                    i + U.grid.extent[0] * (j + 1 + U.grid.extent[1] * nv);
                const int ijm =
                    i + U.grid.extent[0] * (j - 1 + U.grid.extent[1] * nv);

                auto const tmp = sL.U[ijp] + sL.U[ijm] - 2. * sL.U[ij];
                sO.U[ij] = sL.U[ij] - 1. / 24. * tmp;
              };
        } else { // dir==1

          for (int nv = 0; nv < Tsystem::num_vars; ++nv)
#pragma omp parallel for
            for (int j = 0; j < U.grid.extent[1]; ++j)
#pragma omp simd
              for (int i = 1; i < U.grid.extent[0] - 1; ++i) {
                const int ij =
                    i + U.grid.extent[0] * (j + U.grid.extent[1] * nv);
                const int ipj =
                    i + 1 + U.grid.extent[0] * (j + U.grid.extent[1] * nv);
                const int imj =
                    i - 1 + U.grid.extent[0] * (j + U.grid.extent[1] * nv);

                auto const tmp = (sL.U[ipj] + sL.U[imj] - 2. * sL.U[ij]);
                sO.U[ij] = sL.U[ij] - 1. / 24. * tmp;
              };

        }; // switch dir

	sO.is_primitive = sL.is_primitive;
      };   // lambda correct_rec

      /// For the high order algorithm need to correct using the transverse
      // Laplacian

      correct_rec(scratchR, scratchR2);
      correct_rec(scratchL, scratchL2);

    TS::fill_aux(scratchR2);
    TS::fill_aux(scratchL2);

    }; // ndim ==2

    if (Tsystem::ndim == 3 && high_order) {
      assert(!"Not implemented yet");
    }

    TRiem::template solve<dir>(scratchL, scratchR, scratch);

    // Call Riemann solver
    if (TS::ndim > 1 && high_order){
      TRiem::template solve<dir>(scratchL2, scratchR2, scratchF);
    } 

    //We use scratchF in the rest of the routine...
    if(TS::ndim ==1 || !high_order ) scratchF = scratch;


    if (Tsystem::ndim == 2 && high_order) {
        // Need to correct with the transverse laplacian
        if (dir == 0) {
          for (int nv = 0; nv < Tsystem::num_vars; ++nv)
#pragma omp parallel for
            for (int j = 1; j < U.grid.extent[1] - 1; ++j)
              for (int i = 0; i < U.grid.extent[0]; ++i) {
                const int ij =
                    i + U.grid.extent[0] * (j + U.grid.extent[1] * nv);
                const int ijp =
                    i + U.grid.extent[0] * (j + 1 + U.grid.extent[1] * nv);
                const int ijm =
                    i + U.grid.extent[0] * (j - 1 + U.grid.extent[1] * nv);

                auto const tmp = scratch[ijp] + scratch[ijm] - 2. * scratch[ij];
                scratchF.U[ij] += 1. / 24. * tmp;
              };
        } else { // dir==1

          for (int nv = 0; nv < Tsystem::num_vars; ++nv)
#pragma omp parallel for
            for (int j = 0; j < U.grid.extent[1]; ++j)
#pragma omp simd
              for (int i = 1; i < U.grid.extent[0] - 1; ++i) {
                const int ij =
                    i + U.grid.extent[0] * (j + U.grid.extent[1] * nv);
                const int ipj =
                    i + 1 + U.grid.extent[0] * (j + U.grid.extent[1] * nv);
                const int imj =
                    i - 1 + U.grid.extent[0] * (j + U.grid.extent[1] * nv);

                auto const tmp = (scratch[ipj] + scratch[imj] - 2. * scratch[ij]);
                scratchF[ij] += 1. / 24. * tmp;
              };

        }; // switch dir
    }

    if (Tsystem::ndim == 3 && high_order) {
      assert(!"Not implemented yet");
    }

    /*  

    //Apply positivity preserving limiter (FIXME For now only to first component by default!)
  
      //U is naturally UR, but UL is U shifted by one!
      if (Tsystem::ndim == 1) {
        for (int nv = 0; nv < Tsystem::num_vars; ++nv)
          for (int ijk = 1; ijk < U.grid.extent[0] - 1; ++ijk) {
            const int i = ijk + nv * U.grid.extent[0];
            scratchL[i] = U[i-1];
          };
      }

      if (Tsystem::ndim == 2) {

	if(dir==0){
        for (int nv = 0; nv < Tsystem::num_vars; ++nv)
          for (int j = 0; j < U.grid.extent[1]; ++j)
            for (int i = 1; i < U.grid.extent[0] - 1; ++i) {
              const int ij = i + U.grid.extent[0] * (j + U.grid.extent[1] * nv);
              const int imj =
                  i - 1 + U.grid.extent[0] * (j + U.grid.extent[1] * nv);

              scratchL[ij] = U[imj];
	    }
         };

	if(dir==1){
        for (int nv = 0; nv < Tsystem::num_vars; ++nv)
          for (int j = 1; j < U.grid.extent[1]-1; ++j)
            for (int i = 0; i < U.grid.extent[0]; ++i) {
              const int ij = i + U.grid.extent[0] * (j + U.grid.extent[1] * nv);
              const int ijm =
                  i + U.grid.extent[0] * (j -1 + U.grid.extent[1] * nv);
              scratchL[ij] = U[ijm];
	    }
         };

      };

    if (Tsystem::ndim == 3 ) {
      assert(!"Not implemented yet");
    }




      HLL_RiemannSolver<TS,true,true>::template solve<dir>(scratchL,U,scratchR);

      auto const posp = [&] (auto rLFm, auto rm, auto drm, auto rLFp, auto rp,auto drp, auto epsmin){
	      T thetam = 1.;
	      T thetap = 1.;

//	      thetap = std::min( (epsmin - rLFp)/(rp - rLFp + 1.e-200), 1.0);
//	      thetam = std::min( (epsmin - rLFm)/(rm - rLFm + 1.e-200), 1.0);
//
	if(rp < epsmin)
	      thetap = std::min( (epsmin - rLFp)/(drp + 1.e-200), 1.0);
	if(rm < epsmin)
	      thetam = std::min( (epsmin - rLFm)/(drm + 1.e-200), 1.0);

	      auto result = std::max(0.,std::min(thetap,thetam));

//	      assert( result =>0. && result <=1.0);
	      return result;
      }; 

      auto const acfl = 2.*0.45; //Strict CFL condition!!

      if (Tsystem::ndim == 1) {
	  int const nv = 0; //FIXME!!
          for (int ijk = 1; ijk < U.grid.extent[0] - 1; ++ijk) {
            const int i = ijk + nv * U.grid.extent[0];
	    auto rLFm = scratchC[i-1] - acfl*scratchR[i];
	    auto rm = scratchC[i-1] - acfl*scratchF[i];
	    auto drm = - acfl*scratchF[i] + acfl*scratchR[i];
	    auto rLFp = scratchC[i] + acfl*scratchR[i];
	    auto rp = scratchC[i] + acfl*scratchF[i];
	    auto drp = acfl*scratchF[i] - acfl*scratchR[i];

              scratchL[i] = posp(rLFm, rm, drm, rLFp, rp,drp, TS::rho_atmo);
          };
      }

      if (Tsystem::ndim == 2) {
	  int const nv = 0; //FIXME!!

	if(dir==0){
          for (int j = 0; j < U.grid.extent[1]; ++j)
            for (int i = 1; i < U.grid.extent[0] - 1; ++i) {
              const int ij = i + U.grid.extent[0] * (j + U.grid.extent[1] * nv);
              const int imj =
                  i - 1 + U.grid.extent[0] * (j + U.grid.extent[1] * nv);
	      auto rLFm = scratchC[imj] - acfl*scratchR[ij];
	      auto rm = scratchC[imj] - acfl*scratchF[ij];
	      auto drm = - acfl*scratchF[ij] + acfl*scratchR[ij];
	      auto rLFp = scratchC[ij] + acfl*scratchR[ij];
	      auto rp = scratchC[ij] + acfl*scratchF[ij];
	      auto drp = acfl*scratchF[ij] - acfl*scratchR[ij];

              scratchL[ij] = posp(rLFm, rm,drm, rLFp, rp,drp, TS::rho_atmo);
	    }
         };

	if(dir==1){
          for (int j = 1; j < U.grid.extent[1]-1; ++j)
            for (int i = 0; i < U.grid.extent[0]; ++i) {
              const int ij = i + U.grid.extent[0] * (j + U.grid.extent[1] * nv);
              const int ijm =
                  i + U.grid.extent[0] * (j -1 + U.grid.extent[1] * nv);
	      auto rLFm = scratchC[ijm] - acfl*scratchR[ij];
	      auto rm = scratchC[ijm] - acfl*scratchF[ij];
	      auto drm = - acfl*scratchF[ij] + acfl*scratchR[ij];
	      auto rLFp = scratchC[ij] + acfl*scratchR[ij];
	      auto rp = scratchC[ij] + acfl*scratchF[ij];
	      auto drp = acfl*scratchF[ij] - acfl*scratchR[ij];

	      scratchL[ij] = posp(rLFm, rm,drm, rLFp, rp, drp, TS::rho_atmo);
	    }
         };

      };




    if (Tsystem::ndim == 3 ) {
      assert(!"Not implemented yet");
    }

    //Hybridize fluxes!
        for (int nv = 0; nv < Tsystem::num_vars; ++nv)
         for (int ijk = 0; ijk < U.ndof; ++ijk) {
		size_t const ijkO = ijk + U.ndof*nv;
		scratchF[ijkO] = scratchL[ijk]*scratchF[ijkO] + (1.-scratchL[ijk])*scratchR[ijkO];
          };
*/

	
/*  
    MinMod_Reconstruct<false,T>::template reconstruct<TS::ndim, dir>(scratchP, scratchR, scratchL);
      HLL_RiemannSolver<TS,true,false>::template solve<dir>(scratchL,scratchR, scratch);

      //Test PLUTO fix with 1st order
      pluto_flattener<Tsystem>(scratchP,scratchL.U,scratchR.U);

    //Hybridize fluxes!
        for (int nv = 0; nv < Tsystem::num_vars; ++nv)
         for (int ijk = 0; ijk < U.ndof; ++ijk) {
		size_t const ijkO = ijk + U.ndof*nv;
//		PLUTO works the other way!
		scratchF[ijkO] = (1.-scratchL[ijk])*scratchF[ijkO] + scratchL[ijk]*scratch[ijkO];
          };
	  */
     

  };

public:
  McCorquodale_FV<high_order, Tsystem, T, grid_t, Treconstruct, TriemannSolver>(
      grid_t &_grid)
      : grid(_grid), scratch(_grid), scratchL(_grid), scratchR(_grid),
        scratchL2(_grid), scratchR2(_grid), scratchF(_grid), scratchC(_grid), scratchP(_grid) {};
       // ,final_flux(_grid){};


  template <typename Tstorage> inline void switch_to_volume(Tstorage &U){

     if(!high_order) return;

      if (Tsystem::ndim == 1) {
        for (int nv = 0; nv < Tsystem::num_vars; ++nv)
#pragma omp simd
          for (int ijk = 1; ijk < U.grid.extent[0] - 1; ++ijk) {
            const int i = ijk + nv * U.grid.extent[0];
            scratch[i] = U[i] + 1. / 24. * (U[i + 1] + U[i - 1] - 2. * U[i]);
          };
      }

      if (Tsystem::ndim == 2) {
        for (int nv = 0; nv < Tsystem::num_vars; ++nv)
#pragma omp parallel for
          for (int j = 0; j < U.grid.extent[1]; ++j)
#pragma omp simd
            for (int i = 1; i < U.grid.extent[0] - 1; ++i) {
              const int ij = i + U.grid.extent[0] * (j + U.grid.extent[1] * nv);
              const int ipj =
                  i + 1 + U.grid.extent[0] * (j + U.grid.extent[1] * nv);
              const int imj =
                  i - 1 + U.grid.extent[0] * (j + U.grid.extent[1] * nv);

              scratch[ij] = U[ij] + 1. / 24 * (U[ipj] + U[imj] - 2. * U[ij]);
            };

        for (int nv = 0; nv < Tsystem::num_vars; ++nv)
#pragma omp parallel for
          for (int j = 1; j < U.grid.extent[1] - 1; ++j)
            for (int i = 0; i < U.grid.extent[0]; ++i) {
              const int ij = i + U.grid.extent[0] * (j + U.grid.extent[1] * nv);
              const int ijp =
                  i + U.grid.extent[0] * (j + 1 + U.grid.extent[1] * nv);
              const int ijm =
                  i + U.grid.extent[0] * (j - 1 + U.grid.extent[1] * nv);

              scratch[ij] += 1. / 24. * (U[ijp] + U[ijm] - 2. * U[ij]);
            };
      };

      if (Tsystem::ndim == 3) {
        assert(!"Not implemented yet");
      }

      //Switch pointers
      U = scratch;

  };


  template <typename Tstorage> inline void switch_to_point(Tstorage &U){

     if(!high_order) return;

      if (Tsystem::ndim == 1) {
        for (int nv = 0; nv < Tsystem::num_vars; ++nv)
#pragma omp simd
          for (int ijk = 1; ijk < U.grid.extent[0] - 1; ++ijk) {
            const int i = ijk + nv * U.grid.extent[0];
            scratch.U[i] =
                U.U[i] - 1. / 24. * (U.U[i + 1] + U.U[i - 1] - 2. * U.U[i]);
          };
      }

      if (Tsystem::ndim == 2) {
        for (int nv = 0; nv < Tsystem::num_vars; ++nv)
#pragma omp parallel for
          for (int j = 0; j < U.grid.extent[1]; ++j)
#pragma omp simd
            for (int i = 1; i < U.grid.extent[0] - 1; ++i) {
              const int ij = i + U.grid.extent[0] * (j + U.grid.extent[1] * nv);
              const int ipj =
                  i + 1 + U.grid.extent[0] * (j + U.grid.extent[1] * nv);
              const int imj =
                  i - 1 + U.grid.extent[0] * (j + U.grid.extent[1] * nv);

              scratch.U[ij] =
                  U.U[ij] - 1. / 24 * (U.U[ipj] + U.U[imj] - 2. * U.U[ij]);
            };

        for (int nv = 0; nv < Tsystem::num_vars; ++nv)
#pragma omp parallel for
          for (int j = 1; j < U.grid.extent[1] - 1; ++j)
            for (int i = 0; i < U.grid.extent[0]; ++i) {
              const int ij = i + U.grid.extent[0] * (j + U.grid.extent[1] * nv);
              const int ijp =
                  i + U.grid.extent[0] * (j + 1 + U.grid.extent[1] * nv);
              const int ijm =
                  i + U.grid.extent[0] * (j - 1 + U.grid.extent[1] * nv);

              scratch.U[ij] -= 1. / 24. * (U.U[ijp] + U.U[ijm] - 2. * U.U[ij]);
            };
      };

      if (Tsystem::ndim == 3) {
        assert(!"Not implemented yet");
      }

      //Switch pointers
      U = scratch;

  };

  template <typename Tstorage> inline decltype(auto) advect(Tstorage &U, Tstorage & final_flux) {

    static_assert(std::is_same<grid_t, typename Tstorage::grid_t>::value,
                  "Grids don't match");

    // Need some scratch space

    std::memset(scratch.U, 0,
                U.ndof * Tstorage::nsystem * sizeof(typename Tstorage::data_t));

    // Step 0. (only if Tsystem::needs_c2p)
    // Compute values at cell centre

    if (Tsystem::needs_c2p) {

    std::memcpy(scratchC.U,U.U, U.ndof * Tstorage::nsystem * sizeof(typename Tstorage::data_t));


      if(high_order){



	
	      switch_to_point(U);
	      U = scratch;
	      // Finally compute the primitive values at the centre of the cell
	      Tsystem::switch_to_prims(U);
	      Tsystem::switch_to_prims(scratch);
    	      std::memcpy(scratchP.U,U.U, U.ndof * Tstorage::nsystem * sizeof(typename Tstorage::data_t));
	      
	      // Now need to translate back
	      if (Tsystem::ndim == 1) {
		for (int nv = 0; nv < Tsystem::num_vars; ++nv)
#pragma omp simd
		  for (int ijk = 1; ijk < U.grid.extent[0] - 1; ++ijk) {
		    const int i = ijk + nv * U.grid.extent[0];
		    scratch.U[i] += 1. / 24. * (U.U[i + 1] + U.U[i - 1] - 2. * U.U[i]);
		  };
	      }

	      if (Tsystem::ndim == 2) {
		for (int nv = 0; nv < Tsystem::num_vars; ++nv)
#pragma omp parallel for
		  for (int j = 0; j < U.grid.extent[1]; ++j)
#pragma omp simd
		    for (int i = 1; i < U.grid.extent[0] - 1; ++i) {
		      const int ij = i + U.grid.extent[0] * (j + U.grid.extent[1] * nv);
		      const int ipj =
			  i + 1 + U.grid.extent[0] * (j + U.grid.extent[1] * nv);
		      const int imj =
			  i - 1 + U.grid.extent[0] * (j + U.grid.extent[1] * nv);

		      scratch.U[ij] += 1. / 24 * (U.U[ipj] + U.U[imj] - 2. * U.U[ij]);
		    };

		for (int nv = 0; nv < Tsystem::num_vars; ++nv)
#pragma omp parallel for
		  for (int j = 1; j < U.grid.extent[1] - 1; ++j)
		    for (int i = 0; i < U.grid.extent[0]; ++i) {
		      const int ij = i + U.grid.extent[0] * (j + U.grid.extent[1] * nv);
		      const int ijp =
			  i + U.grid.extent[0] * (j + 1 + U.grid.extent[1] * nv);
		      const int ijm =
			  i + U.grid.extent[0] * (j - 1 + U.grid.extent[1] * nv);

		      scratch.U[ij] += 1. / 24. * (U.U[ijp] + U.U[ijm] - 2. * U.U[ij]);
		    };
	      };

	      if (Tsystem::ndim == 3) {
		assert(!"Not implemented yet");
	      }

	      //swap pointers
	      U=scratch;


/*  
	      Tsystem::switch_to_prims(U);
    	      std::memcpy(scratchP.U,U.U, U.ndof * Tstorage::nsystem * sizeof(typename Tstorage::data_t));
	      Tsystem::switch_to_cons(U);

	      switch_to_point(U);
	      Tsystem::switch_to_prims(U);
	      switch_to_volume(U);
*/

      }else{
	      Tsystem::switch_to_prims(U);
    	      std::memcpy(scratchP.U,U.U, U.ndof * Tstorage::nsystem * sizeof(typename Tstorage::data_t));
      }//low order


    } // If high_order && needs_c2p

    calculate_flux<0>(U);

    if (Tsystem::ndim == 1) {
      for (int nv = 0; nv < Tsystem::num_vars; ++nv)
#pragma omp simd
        for (int ijk = scratchF.grid.gz[0];
             ijk < scratchF.grid.extent[0] - scratchF.grid.gz[0]; ++ijk) {
          int const i = ijk + scratchF.grid.extent[0] * nv;
          int const ip = ijk + scratchF.grid.extent[0] * nv + 1;
          final_flux[i] = (scratchF[i] - scratchF[ip]) * scratchF.grid.idx[0];
        };
    };

    if (Tsystem::ndim == 2) {

      for (int nv = 0; nv < Tsystem::num_vars; ++nv)
#pragma omp parallel for
        for (int j = scratchF.grid.gz[1];
             j < scratchF.grid.extent[1] - scratchF.grid.gz[1]; ++j)
#pragma omp simd
          for (int i = scratchF.grid.gz[0];
               i < scratchF.grid.extent[0] - scratchF.grid.gz[0]; ++i) {

            int const ij =
                i + scratchF.grid.extent[0] * j +
                scratchF.grid.extent[0] * scratchF.grid.extent[1] * nv;

            int const ipj =
                i + 1 + scratchF.grid.extent[0] * j +
                scratchF.grid.extent[0] * scratchF.grid.extent[1] * nv;

            final_flux.U[ij] =
                (scratchF[ij] - scratchF[ipj]) * scratchF.grid.idx[0];
          };

      calculate_flux<1>(U);

      for (int nv = 0; nv < Tsystem::num_vars; ++nv)
#pragma omp parallel for
        for (int j = scratchF.grid.gz[1];
             j < scratchF.grid.extent[1] - scratchF.grid.gz[1]; ++j)
          for (int i = scratchF.grid.gz[0];
               i < scratchF.grid.extent[0] - scratchF.grid.gz[0]; ++i) {

            int const ij =
                i + scratchF.grid.extent[0] * j +
                scratchF.grid.extent[0] * scratchF.grid.extent[1] * nv;

            int const ijp =
                i + scratchF.grid.extent[0] * (j + 1) +
                scratchF.grid.extent[0] * scratchF.grid.extent[1] * nv;

            final_flux.U[ij] +=
                (scratchF[ij] - scratchF[ijp]) * scratchF.grid.idx[1];
          };
    };

    if (Tsystem::ndim == 3) {
      assert(!"Not implemented yet");
    };


    //Need to switch back to conservatives!

    if (Tsystem::needs_c2p) {
      U=scratchC;
    } // If high_order && needs_c2p

//    return final_flux;
  };
};




template <bool high_order, typename Tsystem, typename T, typename grid_t,
          typename Treconstruct, typename TriemannSolver>
class DelZanna_FD {

private:
  using storage_t = SimpleStorage<grid_t, Tsystem>;

  storage_t scratch, scratchL, scratchR, scratchF,scratchC;
     // ,final_flux;

  using TR = Treconstruct;
  using TS = Tsystem;
  using TRiem = TriemannSolver;

  grid_t &grid;

  template <int dir, typename Tstorage>
  inline void calculate_flux(Tstorage &U) {
    static_assert(std::is_same<grid_t, typename Tstorage::grid_t>::value,
                  "Grids don't match");

    // Actual reconstruction step
    TR::template reconstruct<TS::ndim, dir>(U, scratchR, scratchL);

    //Need to fill aux
    TS::fill_aux(scratchR);
    TS::fill_aux(scratchL);

    //Riemann solver  
    TRiem::template solve<dir>(scratchL, scratchR, scratch);

    //We use scratchF in the rest of the routine...
    if(!high_order) scratchF = scratch;
    else{

      if (Tsystem::ndim == 1) {
        for (int nv = 0; nv < Tsystem::num_vars; ++nv)
#pragma omp simd
          for (int ijk = 1; ijk < U.grid.extent[0] - 1; ++ijk) {
            const int i = ijk + nv * U.grid.extent[0];
            scratchF[i] =
                scratch[i] - 1. / 24. * (scratch[i + 1] + scratch[i - 1] - 2. * scratch[i]);
          };
      }

      if (Tsystem::ndim == 2) {
        for (int nv = 0; nv < Tsystem::num_vars; ++nv)
#pragma omp parallel for
          for (int j = 0; j < U.grid.extent[1]; ++j)
#pragma omp simd
            for (int i = 1; i < U.grid.extent[0] - 1; ++i) {
              const int ij = i + U.grid.extent[0] * (j + U.grid.extent[1] * nv);
              const int ipj =
                  i + 1 + U.grid.extent[0] * (j + U.grid.extent[1] * nv);
              const int imj =
                  i - 1 + U.grid.extent[0] * (j + U.grid.extent[1] * nv);

              scratchF[ij] =
                  scratch[ij] - 1. / 24 * (scratch[ipj] + scratch[imj] - 2. * scratch[ij]);
            };

        for (int nv = 0; nv < Tsystem::num_vars; ++nv)
#pragma omp parallel for
          for (int j = 1; j < U.grid.extent[1] - 1; ++j)
            for (int i = 0; i < U.grid.extent[0]; ++i) {
              const int ij = i + U.grid.extent[0] * (j + U.grid.extent[1] * nv);
              const int ijp =
                  i + U.grid.extent[0] * (j + 1 + U.grid.extent[1] * nv);
              const int ijm =
                  i + U.grid.extent[0] * (j - 1 + U.grid.extent[1] * nv);

              scratchF[ij] -= 1. / 24. * (scratch[ijp] + scratch[ijm] - 2. * scratch[ij]);
            };
      };

      if (Tsystem::ndim == 3) {
        assert(!"Not implemented yet");
      }

    }


  };

public:
  DelZanna_FD<high_order, Tsystem, T, grid_t, Treconstruct, TriemannSolver>(
      grid_t &_grid)
      : grid(_grid), scratch(_grid), scratchL(_grid), scratchR(_grid),
        scratchF(_grid),scratchC(_grid) {};
       // ,final_flux(_grid){};


  template <typename Tstorage> inline void switch_to_volume(Tstorage &U){
     return;
  };

  template <typename Tstorage> inline void switch_to_point(Tstorage &U){

     return;
  };

  template <typename Tstorage> inline decltype(auto) advect(Tstorage &U, Tstorage & final_flux) {

    static_assert(std::is_same<grid_t, typename Tstorage::grid_t>::value,
                  "Grids don't match");

    // Need some scratch space

    std::memset(scratch.U, 0,
                U.ndof * Tstorage::nsystem * sizeof(typename Tstorage::data_t));

    // Step 0. (only if Tsystem::needs_c2p)
    if(Tsystem::needs_c2p){
      std::memcpy(scratchC.U,U.U, U.ndof * Tstorage::nsystem * sizeof(typename Tstorage::data_t));
      Tsystem::switch_to_prims(U);
    }

    calculate_flux<0>(U);

    if (Tsystem::ndim == 1) {
      for (int nv = 0; nv < Tsystem::num_vars; ++nv)
#pragma omp simd
        for (int ijk = scratchF.grid.gz[0];
             ijk < scratchF.grid.extent[0] - scratchF.grid.gz[0]; ++ijk) {
          int const i = ijk + scratchF.grid.extent[0] * nv;
          int const ip = ijk + scratchF.grid.extent[0] * nv + 1;
          final_flux[i] = (scratchF[i] - scratchF[ip]) * scratchF.grid.idx[0];
        };
    };

    if (Tsystem::ndim == 2) {

      for (int nv = 0; nv < Tsystem::num_vars; ++nv)
#pragma omp parallel for
        for (int j = scratchF.grid.gz[1];
             j < scratchF.grid.extent[1] - scratchF.grid.gz[1]; ++j)
#pragma omp simd
          for (int i = scratchF.grid.gz[0];
               i < scratchF.grid.extent[0] - scratchF.grid.gz[0]; ++i) {

            int const ij =
                i + scratchF.grid.extent[0] * j +
                scratchF.grid.extent[0] * scratchF.grid.extent[1] * nv;

            int const ipj =
                i + 1 + scratchF.grid.extent[0] * j +
                scratchF.grid.extent[0] * scratchF.grid.extent[1] * nv;

            final_flux.U[ij] =
                (scratchF[ij] - scratchF[ipj]) * scratchF.grid.idx[0];
          };

      calculate_flux<1>(U);

      for (int nv = 0; nv < Tsystem::num_vars; ++nv)
#pragma omp parallel for
        for (int j = scratchF.grid.gz[1];
             j < scratchF.grid.extent[1] - scratchF.grid.gz[1]; ++j)
          for (int i = scratchF.grid.gz[0];
               i < scratchF.grid.extent[0] - scratchF.grid.gz[0]; ++i) {

            int const ij =
                i + scratchF.grid.extent[0] * j +
                scratchF.grid.extent[0] * scratchF.grid.extent[1] * nv;

            int const ijp =
                i + scratchF.grid.extent[0] * (j + 1) +
                scratchF.grid.extent[0] * scratchF.grid.extent[1] * nv;

            final_flux.U[ij] +=
                (scratchF[ij] - scratchF[ijp]) * scratchF.grid.idx[1];
          };
    };

    if (Tsystem::ndim == 3) {
      assert(!"Not implemented yet");
    };

    if (Tsystem::needs_c2p) U=scratchC;

//    return final_flux;
  };
};

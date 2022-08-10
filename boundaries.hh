/*
 * =====================================================================================
 *
 *       Filename:  boundaries.hh
 *
 *    Description:  Simple boundaries
 *
 *        Version:  1.0
 *        Created:  15.06.2018 18:45:36
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Elias Roland Most (ERM), emost@itp.uni-frankfurt.de
 *   Organization:  Goethe University Frankfurt
 *
 * =====================================================================================
 */

#pragma once

#include <cassert>

template <int ndim, typename Tstorage> inline void PeriodicBC(Tstorage &U) {

  if (ndim == 1) {
    for (int nv = 0; nv < Tstorage::nsystem; ++nv)
      for (int ijk = 0; ijk < U.grid.gz[0]; ++ijk) {
        int const i = ijk + U.grid.extent[0] * nv;
        int const ip = U.grid.extent[0]*(nv+1) - 2 * U.grid.gz[0] + i;

        U[i] = U[ip];
        U[ip + U.grid.gz[0]] = U[i + U.grid.gz[0]];
      };
  };

  if (ndim == 2) {
    for (int nv = 0; nv < Tstorage::nsystem; ++nv)
#pragma omp parallel for
      for (int j = 0; j < U.grid.gz[1]; ++j)
        for (int i = 0; i < U.grid.extent[0]; ++i) {

          int const ij = i + U.grid.extent[0] * j +
                         U.grid.extent[0] * U.grid.extent[1] * nv;

          int const ijp =
              i + U.grid.extent[0] * (j + U.grid.extent[1] - 2 * U.grid.gz[1]) +
              U.grid.extent[0] * U.grid.extent[1] * nv;

          U[ij] = U[ijp];
          const int offset = U.grid.extent[0] * U.grid.gz[1];
          U[ijp + offset] = U[ij + offset];
        };

    for (int nv = 0; nv < Tstorage::nsystem; ++nv)
#pragma omp parallel for
      for (int j = 0; j < U.grid.extent[1]; ++j)
        for (int i = 0; i < U.grid.gz[0]; ++i) {

          int const ij = i + U.grid.extent[0] * j +
                         U.grid.extent[0] * U.grid.extent[1] * nv;

          int const ipj = i+ U.grid.extent[0] - 2 * U.grid.gz[0] +
                          U.grid.extent[0] * j +
                          U.grid.extent[0] * U.grid.extent[1] * nv;

          U[ij] = U[ipj];
          const int offset = U.grid.gz[0];
          U[ipj + offset] = U[ij + offset];
        };
  };

  if (ndim == 3)
    assert(!"Not implemented yet!");
};



template <int ndim, typename Tstorage> inline void FlatBC(Tstorage &U) {

  if (ndim == 1) {
    for (int nv = 0; nv < Tstorage::nsystem; ++nv){
      int const ilow = U.grid.gz[0] + U.grid.extent[0] * nv;
      int const ihigh = U.grid.extent[0] * (nv+1) - U.grid.gz[0] - 1;

      auto const Ulow = U[ilow];
      auto const Uhigh = U[ihigh];

      for (int ijk = 0; ijk < U.grid.gz[0]; ++ijk) {
        U[ilow - ijk-1] = Ulow;
        U[ihigh + ijk+1] = Uhigh;
      };
    }
  };

  if (ndim == 2) {
#pragma omp parallel for
    for (int nv = 0; nv < Tstorage::nsystem; ++nv)
      for (int j = 0; j < U.grid.gz[1]; ++j)
        for (int i = 0; i < U.grid.extent[0]; ++i) {

          int const ijL = i + U.grid.extent[0] * U.grid.gz[1] +
                         U.grid.extent[0] * U.grid.extent[1] * nv;

          int const ijH =
              i + U.grid.extent[0] * (U.grid.extent[1] - U.grid.gz[1] -1) +
              U.grid.extent[0] * U.grid.extent[1] * nv;

          U[ijL - U.grid.extent[0]*(j-1) ] = U[ijL];
          U[ijL + U.grid.extent[0]*(j+1) ] = U[ijH];
        };

    for (int nv = 0; nv < Tstorage::nsystem; ++nv)
#pragma omp parallel for
      for (int j = 0; j < U.grid.extent[1]; ++j){

	int const ilow = U.grid.gz[0] + U.grid.extent[0] *j +
                         U.grid.extent[0] * U.grid.extent[1] * nv;
	int const ihigh = U.grid.extent[0] * (j+1) - U.grid.gz[0] - 1 +
                         U.grid.extent[0] * U.grid.extent[1] * nv;

	auto const Ulow = U[ilow];
	auto const Uhigh = U[ihigh];
	for (int ijk = 0; ijk < U.grid.gz[0]; ++ijk) {

	  U[ilow - ijk-1] = Ulow;
	  U[ihigh + ijk+1] = Uhigh;
	};

      }
  };

  if (ndim == 3)
    assert(!"Not implemented yet!");
};

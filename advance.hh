/*
 * =====================================================================================
 *
 *       Filename:  advance.hh
 *
 *    Description:  Simple time stepper
 *
 *        Version:  1.0
 *        Created:  13.06.2018 18:48:11
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
#include <cassert>
#include <cmath>
#include <functional>

template <typename Tsystem, typename Tstorage> class TimeStepperRK4 {
  // Note that this class is *NOT* using low memory version of this
  // algorithm!

public:
  using data_t = typename Tstorage::data_t;
  using grid_t = typename Tstorage::grid_t;

private:
  Tstorage k1, k2, k3, k4, ktmp; // Runge Kutta coefficients

  data_t delta_t = 0;

  data_t current_time = 0;

  int iteration = 0;

public:
  TimeStepperRK4<Tsystem, Tstorage>(grid_t &grid)
      : k1(Tstorage(grid)), k2(Tstorage(grid)), k3(Tstorage(grid)),
        k4(Tstorage(grid)), ktmp(Tstorage(grid)){};

  inline data_t get_current_time() const { return current_time; };
  inline int get_iteration() const { return iteration; };

  inline void set_timestep(data_t dt){ delta_t = dt;};

  inline void set_timestep_from_CFL(data_t CFL) {

    assert(CFL < 1.);

    data_t dx = 1.e99;
    // Get min direction
    for (int i = 0; i < Tsystem::ndim; ++i)
      dx = std::min(dx, k1.grid.dx[i]);

    delta_t = CFL * dx;
  };

  template <typename Tstorage2, typename Tf>
  inline void advance(Tstorage2 &U, Tf &func) {

    auto fma = [&](Tstorage &a, Tstorage &b, data_t c, Tstorage &d) {
      if (Tsystem::ndim == 1) {

        for (int nv = 0; nv < Tsystem::num_vars; ++nv)
#pragma omp simd
          for (int ijk = 0; ijk < U.grid.extent[0];
               ++ijk) {
            int const i = ijk + U.grid.extent[0] * nv;
            a.U[i] = b.U[i] + c * d.U[i];
          };

      } // nimd==1

      if (Tsystem::ndim == 2) {

        for (int nv = 0; nv < Tsystem::num_vars; ++nv)
#pragma omp parallel for
          for (int j = 0; j < U.grid.extent[1] ; ++j)
#pragma omp simd
            for (int i = 0; i < U.grid.extent[0];
                 ++i) {

              int const ij =
                  i + U.grid.extent[0] * j
                    + U.grid.extent[0] * U.grid.extent[1] * nv;

              a.U[ij] = b.U[ij] + c * d.U[ij];
            };

      }; // ndim==2

      if (Tsystem::ndim == 3) {
        assert(!"Not implemented yet");
      }
    };

    func(U, k1);

    fma(k3, U, 0.5 * delta_t, k1);

    func(k3, k2);

    fma(k4, U, 0.5 * delta_t, k2);

    func(k4, k3);

    fma(ktmp, U, delta_t, k3);

    func(ktmp, k4);

    if (Tsystem::ndim == 1) {

      for (int nv = 0; nv < Tsystem::num_vars; ++nv)
#pragma omp simd
        for (int ijk = 0; ijk < U.grid.extent[0];
             ++ijk) {
          int const i = ijk + U.grid.extent[0] * nv;
          U[i] +=
              delta_t * (1. / 3. * (k2[i] + k3[i]) + 1. / 6. * (k1[i] + k4[i]));
        };

    } // nimd==1

    if (Tsystem::ndim == 2) {

      for (int nv = 0; nv < Tsystem::num_vars; ++nv)
#pragma omp parallel for
        for (int j = 0; j < U.grid.extent[1]; ++j)
#pragma omp simd
          for (int i = 0; i < U.grid.extent[0]; ++i) {

              int const ij =
                  i + U.grid.extent[0] * j
                    + U.grid.extent[0] * U.grid.extent[1] * nv;

            U[ij] += delta_t * (1. / 3. * (k2[ij] + k3[ij]) +
                                1. / 6. * (k1[ij] + k4[ij]));
          };

    }; // ndim==2

    if (Tsystem::ndim == 3) {
      assert(!"Not implemented yet");
    }
    iteration++;
    current_time += delta_t;
  };
};

 template <typename Tsystem, typename Tstorage> class TimeStepperRK3 {
   // Note that this class is *NOT* using low memory version of this
   // algorithm!
 
 public:
   using data_t = typename Tstorage::data_t;
   using grid_t = typename Tstorage::grid_t;
 
 private:
   Tstorage k1, k2, ktmp; // Runge Kutta coefficients
 
   data_t delta_t = 0;
 
   data_t current_time = 0;
 
   int iteration = 0;
 
 public:
   TimeStepperRK3<Tsystem, Tstorage>(grid_t &grid)
       : k1(Tstorage(grid)), k2(Tstorage(grid)) , ktmp(Tstorage(grid)){};
 
   inline data_t get_current_time() const { return current_time; };
   inline int get_iteration() const { return iteration; };
 
   inline void set_timestep(data_t dt){ delta_t = dt;};
 
   inline void set_timestep_from_CFL(data_t CFL) {
 
     assert(CFL < 1.);
 
     data_t dx = 1.e99;
     // Get min direction
     for (int i = 0; i < Tsystem::ndim; ++i)
       dx = std::min(dx, k1.grid.dx[i]);
 
     delta_t = CFL * dx;
   };
 
   template <typename Tstorage2, typename Tf>
   inline void advance(Tstorage2 &U, Tf &func) {
 
     auto fma = [&](Tstorage &a, Tstorage &b, data_t c, Tstorage &d) {
       if (Tsystem::ndim == 1) {
 
         for (int nv = 0; nv < Tsystem::num_vars; ++nv)
 #pragma omp simd
           for (int ijk = 0; ijk < U.grid.extent[0];
                ++ijk) {
             int const i = ijk + U.grid.extent[0] * nv;
             a.U[i] = b.U[i] + c * d.U[i];
           };
 
       } // nimd==1
 
       if (Tsystem::ndim == 2) {
 
         for (int nv = 0; nv < Tsystem::num_vars; ++nv)
 #pragma omp parallel for
           for (int j = 0; j < U.grid.extent[1] ; ++j)
 #pragma omp simd
             for (int i = 0; i < U.grid.extent[0];
                  ++i) {
 
               int const ij =
                   i + U.grid.extent[0] * j
                     + U.grid.extent[0] * U.grid.extent[1] * nv;
 
               a.U[ij] = b.U[ij] + c * d.U[ij];
             };
 
       }; // ndim==2
 
       if (Tsystem::ndim == 3) {
         assert(!"Not implemented yet");
       }
     };
 
     auto fma3 = [&](data_t ca, Tstorage &a, data_t cb, Tstorage &b, data_t cd, Tstorage &d, Tstorage &r) {
       if (Tsystem::ndim == 1) {
 
         for (int nv = 0; nv < Tsystem::num_vars; ++nv)
 #pragma omp simd
           for (int ijk = 0; ijk < U.grid.extent[0];
                ++ijk) {
             int const i = ijk + U.grid.extent[0] * nv;
             r[i] = ca*a[i] + cb*b[i] + cd * d[i];
           };
 
       } // nimd==1
 
       if (Tsystem::ndim == 2) {
 
         for (int nv = 0; nv < Tsystem::num_vars; ++nv)
 #pragma omp parallel for
           for (int j = 0; j < U.grid.extent[1] ; ++j)
 #pragma omp simd
             for (int i = 0; i < U.grid.extent[0];
                  ++i) {
 
               int const ij =
                   i + U.grid.extent[0] * j
                     + U.grid.extent[0] * U.grid.extent[1] * nv;
 
               r[ij] = ca*a[ij] + cb*b[ij] + cd * d[ij];
             };
 
       }; // ndim==2
 
       if (Tsystem::ndim == 3) {
         assert(!"Not implemented yet");
       }
     };
 
     func(U, k1);
 
     fma(ktmp, U, delta_t, k1);
 
     func(ktmp, k1);
 
     fma3(0.75, U, 0.25, ktmp, 0.25*delta_t, k1, k2);
 
     func(k2,k1);
 
     fma3(1./3., U, 2./3., k2, 2./3.*delta_t, k1, U);
 
     iteration++;
     current_time += delta_t;
   };
 };

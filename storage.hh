/*
 * =====================================================================================
 *
 *       Filename:  storage.hh
 *
 *    Description:  Simple storage class for advection toy
 *
 *        Version:  1.0
 *        Created:  11.06.2018 17:55:20
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Elias Roland Most (ERM), emost@itp.uni-frankfurt.de
 *   Organization:  Goethe University Frankfurt
 *
 * =====================================================================================
 */

#pragma once

template <int Tndim = 1, typename T = double> class SimpleCartesianGrid {

public:
  static constexpr int ndim = Tndim;
  std::array<int, ndim> extent;
  std::array<int, ndim> gz;
  std::array<T, ndim> dx, idx;
  std::array<T, ndim> origin;

  int ndof;

  SimpleCartesianGrid<Tndim, T>(std::array<int, ndim> &&extent,
                                std::array<int, ndim> &&gz,
                                std::array<T, ndim> &&dx,
                                std::array<T, ndim> &&origin)
      : extent(std::move(extent)), gz(std::move(gz)), dx(std::move(dx)),
        origin(std::move(origin)) {
    ndof = 1;
    for (auto &a : this->extent)
      ndof *= a;
    for (int i = 0; i < ndim; ++i)
      idx[i] = 1. / dx[i];
  };

  SimpleCartesianGrid<Tndim, T>() = default;

  inline T get_coords(int dim, int index) {
    return origin[dim] + index * dx[dim];
  };
};

template <typename Tgrid, typename Tsystem, typename T = double>
class SimpleStorage {

public:
  using this_t = SimpleStorage<Tgrid, Tsystem, T>;
  using data_t = T;
  using grid_t = Tgrid;

  bool is_primitive = false;
  T *U = nullptr;
  T *aux = nullptr;

  grid_t &grid;

  int ndof;

  static constexpr int nsystem = Tsystem::num_vars;
  static constexpr int naux = Tsystem::num_aux;

  SimpleStorage<Tgrid, Tsystem, T>(grid_t &grid) : grid(grid), ndof(grid.ndof) {
    U = new T[ndof * nsystem];

    if (naux > 0)
      aux = new T[naux * ndof];
  };

  SimpleStorage<Tgrid, Tsystem, T>(SimpleStorage<Tgrid, Tsystem, T> const &A)
      : is_primitive(A.is_primitive), U(A.U), aux(A.aux), grid(A.grid),
        ndof(A.grid.ndof){};

  inline double &operator[](int ijk) { return U[ijk]; };

  inline SimpleStorage<Tgrid, Tsystem, T> &
  operator=(SimpleStorage<Tgrid, Tsystem, T> &A) {

    std::swap(is_primitive, A.is_primitive);
    std::swap(U, A.U);
    std::swap(aux, A.aux);
    std::swap(grid, A.grid);
    std::swap(ndof, A.ndof);

    //Assignment would be really bad as the storage object would no longer
    //be independent
//    is_primitive = A.is_primitive;
//    U = A.U;
//    aux = A.aux;
//    grid = A.grid;
//    ndof = grid.ndof;

    return *this;
  };


  //   ~SimpleStorage<Tgrid,Tsystem,T>(){
  //     if(U!=nullptr) delete U;
  //     if(aux !=nullptr) delete aux;
  //   };
};

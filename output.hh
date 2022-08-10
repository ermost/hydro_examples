/*
 * =====================================================================================
 *
 *       Filename:  output.hh
 *
 *    Description:  Simple text output
 *
 *        Version:  1.0
 *        Created:  15.06.2018 22:13:37
 *       Revision:  none
 *       Compiler:  gcc
 *
 *         Author:  Elias Roland Most (ERM), emost@itp.uni-frankfurt.de
 *   Organization:  Goethe University Frankfurt
 *
 * =====================================================================================
 */

#pragma once

#include <fstream>
#include <iostream>
#include <sstream>
#include <iomanip>
#include <string>

template <typename Tstorage> inline bool output(Tstorage &U, int iteration) {
  std::stringstream ss;
  ss << "output/output_" << std::setfill('0') << std::setw(8) << iteration << ".asc";
  std::ofstream ostrm(ss.str());

  assert(ostrm.good() && "'output' folder not present");

  for (int ijk = 0; ijk < U.ndof; ++ijk) {
      int ijkN = ijk;
    for (int nn = 0; nn < Tstorage::grid_t::ndim; ++nn){
      int d = U.grid.extent[nn];
      int ijkL= ijkN % d;
      if(nn == Tstorage::grid_t::ndim -1) ijkL = ijkN;
      ijkN/=d;
      ostrm << U.grid.get_coords(nn, ijkL) << "\t";
    }
    for (int nv = 0; nv < Tstorage::nsystem; ++nv)
      ostrm << U[ijk + U.ndof * nv] << "\t";
    ostrm << std::endl;
  };

  ostrm.flush();
  ostrm.close();

  return true;
};

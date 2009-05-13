// QuCoSi - Quantum Computer Simulation
// Copyright © 2009 Frank S. Thomas <frank@thomas-alfeld.de>
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses/>.

#ifndef QUCOSI_AUX_H
#define QUCOSI_AUX_H

#include <limits>

#include <Eigen/Core>

namespace QuCoSi {

typedef double fptype;
typedef std::complex<fptype> field;
typedef Eigen::Matrix<field, Eigen::Dynamic, 1> VectorXc;
typedef Eigen::Matrix<field, Eigen::Dynamic, Eigen::Dynamic> MatrixXc;

inline bool isZero(const fptype x)
{
  return std::abs(x) <= std::numeric_limits<fptype>::epsilon();
}

inline bool isOne(const fptype x)
{
  return isZero(x-1.);
}

inline int log2(unsigned value)
{
  unsigned l = 0;
  while ((value >> l) != 0) {
    l++;
  }
  return l-1;
}

} // namespace QuCoSi

#endif // QUCOSI_AUX_H

// vim: shiftwidth=2 textwidth=78

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

#ifndef QUCOSI_BASE_H
#define QUCOSI_BASE_H

#include <Eigen/Core>

namespace QuCoSi {

typedef double fptype;
typedef std::complex<fptype> field;
typedef Eigen::Matrix<field, 2, 1> Vector;
typedef Eigen::Matrix<field, Eigen::Dynamic, 1> MultiVector;
typedef Eigen::Matrix<field, Eigen::Dynamic, 1> DynamicVector;

class Tensor : public DynamicVector {};

class Base 

} // namespace QuCoSi

#endif // QUCOSI_BASE_H

// vim: shiftwidth=2 textwidth=78

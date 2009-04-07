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

#ifndef QUCOSI_MULTIQUBIT_H
#define QUCOSI_MULTIQUBIT_H

#include <Base.h>
#include <Qubit.h>

namespace QuCoSi {

class MultiQubit : public MultiVector {
  public:
    inline MultiQubit() : MultiVector() {}

    inline MultiQubit(const int size) : MultiVector(size) {}

    inline MultiQubit& operator=(const Qubit& v)
    {
      MultiVector::operator=(v);
      return *this;
    }

    inline MultiQubit& otimes(const Qubit& v)
    {
      if (rows() == 0) {
        resize(v.rows());
        *this = v;
      } else {
        MultiQubit tmp(rows()*v.rows());

        int n = 0;
        for (int i = 0; i < rows(); i++) {
          for (int j = 0; j < v.rows(); j++, n++) {
            tmp[n] = (*this)[i] * v[j];
          }
        } 
        *this = tmp;
      }
      return *this;
    }

    inline MultiQubit& otimes(const Qubit& v1, const Qubit& v2)
    {
      return otimes(v1).otimes(v2);
    }

    inline MultiQubit& otimes(const Qubit& v1, const Qubit& v2,
                              const Qubit& v3)
    {
      return otimes(v1).otimes(v2).otimes(v3);
    }

    inline MultiQubit& otimes(const Qubit& v1, const Qubit& v2,
                              const Qubit& v3, const Qubit& v4)
    {
      return otimes(v1).otimes(v2).otimes(v3).otimes(v4);
    }
};

} // namespace QuCoSi

#endif // QUCOSI_MULTIQUBIT_H

// vim: shiftwidth=2 textwidth=78

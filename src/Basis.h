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

#ifndef QUCOSI_BASIS_H
#define QUCOSI_BASIS_H

#include <ostream>
#include <vector>

#include <Vector.h>

namespace QuCoSi {

class Basis : public std::vector<Vector>
{
  public:
    inline Basis()
    {
      setStandardBasis(2);
    }

    inline Basis(const int dim)
    {
      setStandardBasis(dim);
    }

    inline Basis& setStandardBasis(const int dim)
    {
      resize(dim);
      for (int i = 0; i < dim; i++) {
        Vector e(dim);
        e(i) = field(1,0);
        (*this)[i] = e;
      }
      return *this;
    }

    inline Basis tensorDot(const Basis& f) const
    {
      Basis g(size()*f.size());
      int k = 0;
      for (int i = 0; i < size(); i++) {
        for (int j = 0; j < f.size(); j++, k++) {
          g[k] = (*this)[i].tensorDot(f[j]);
        }
      }
      return g;
    }

    inline Basis& tensorDotSet(const Basis& f)
    {
      if (size() > 0 && f.size() > 0) {
        *this = tensorDot(f);
      }
      return *this;
    }

    inline bool isNormalized() const
    {
      for (int i = 0; i < size(); i++) {
        if (!at(i).isNormalized()) {
          return false;
        }
      }
      return true;
    }

    inline bool isOrthogonal() const
    {
      int n = size();
      for (int i = 0; i < n; i++) {
        for (int j = i + 1; j < n; j++) {
          if (!at(i).isOrthogonal(at(j))) {
            return false;
          }
        }
      }
      return true;
    }

    inline bool isOrthonormal() const
    {
      if (isNormalized()) {
        if (isOrthogonal()) {
          return true;
        }
      }
      return false;
    }
};

inline std::ostream& operator<<(std::ostream& os, const Basis& b)
{
  int n = b.size()-1;
  os << "{" << std::endl;
  for (int i = 0; i < n; i++) {
    os << b.at(i) << "," << std::endl << std::endl;
  }
  if (n >= 0) {
    os << b.at(n) << std::endl;
  }
  return os << "}";
}

} // namespace QuCoSi

#endif // QUCOSI_BASIS_H

// vim: shiftwidth=2 textwidth=78

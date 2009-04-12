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

class Basis : public std::vector<Vector> {
  public:
    inline Basis() {}

    inline Basis(const int dim)
    {
      setNaturalBasis(dim);
    }

    inline Basis& setNaturalBasis(const int dim)
    {
      clear();
      for (int i = 0; i < dim; i++) {
        Vector e(dim);
        e[i] = field(1,0);
        push_back(e);
      }
      return *this;
    }

    inline Basis otimes(const Basis& f) const
    {
      Basis g(size()*f.size());

      int k = 0;
      for (int i = 0; i < size(); i++) {
        for (int j = 0; j < f.size(); j++, k++) {
          g[k] = (*this)[i].otimes(f[j]);
        }
      }
      return g;
    }

    inline Basis& otimesSet(const Basis& f)
    {
      if (size() > 0 && f.size() > 0) {
        *this = otimes(f);
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
};

inline std::ostream& operator<<(std::ostream& os, const Basis& b)
{
  int n = b.size();
  os << "{" << std::endl;
  for (int i = 0; i < n; i++) {
    os << b.at(i);
    if (i < n-1) {
      os << "," << std::endl << std::endl;
    }
  }
  return os << std::endl << "}";
}

} // namespace QuCoSi

#endif // QUCOSI_BASIS_H

// vim: shiftwidth=2 textwidth=78

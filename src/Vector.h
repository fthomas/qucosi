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

#ifndef QUCOSI_VECTOR_H
#define QUCOSI_VECTOR_H

#include <limits>

#include <Eigen/Core>

namespace QuCoSi {

typedef double fptype;
typedef std::complex<fptype> field;
typedef Eigen::Matrix<field, Eigen::Dynamic, 1> VectorXc;

inline bool isZero(const fptype x)
{
  return std::abs(x) <= std::numeric_limits<fptype>::epsilon();
}

inline bool isOne(const fptype x)
{
  return isZero(x-1.);
}

class Vector : public VectorXc {
  public:
    inline Vector() : VectorXc(2)
    {
      (*this)[0] = field(1,0);
      (*this)[1] = field(0,0);
    }

    inline Vector(const int dim) : VectorXc(dim) {}

    inline Vector(const field& c0, const field& c1) : VectorXc(2)
    {
      (*this)[0] = c0;
      (*this)[1] = c1;
    }

    inline Vector& operator=(const VectorXc& v)
    {
      VectorXc::operator=(v);
      return *this;
    }

    inline bool isNormalized() const
    {
      return isOne(norm());
    }

    inline Vector& randomize()
    {
      setRandom().normalize();
      return *this;
    }

    inline Vector otimes(const Vector& v) const
    {
      Vector w(rows()*v.rows());

      int k = 0;
      for (int i = 0; i < rows(); i++) {
        for (int j = 0; j < v.rows(); j++, k++) {
          w[k] = (*this)[i]*v[j];
        }
      }
      return w;
    }

    inline Vector& otimesSet(const Vector& v)
    {
      if (rows() > 0 && v.rows() > 0) {
        *this = otimes(v);
      }
      return *this;
    }
};

} // namespace QuCoSi

#endif // QUCOSI_VECTOR_H

// vim: shiftwidth=2 textwidth=78

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

#ifndef QUCOSI_QUBIT_H
#define QUCOSI_QUBIT_H

#include <cstdlib>
#include <vector>

#include <Aux.h>
#include <Vector.h>

namespace QuCoSi {

class Qubit : public Vector {
  public:
    inline Qubit() : Vector(2) {}

    inline Qubit(const int dim) : Vector(dim) {}

    inline Qubit(const field& c0, const field& c1) : Vector(c0, c1) {}

    inline Qubit(const int x, const int n) : Vector(std::pow(2,n))
    {
      (*this)(x) = field(1,0);
    }

    inline Qubit& operator=(const VectorXc& v)
    {
      VectorXc::operator=(v);
      return *this;
    }

    inline Qubit& operator=(const Vector& v)
    {
      Vector::operator=(v);
      return *this;
    }

    inline bool isPureState() const
    {
      for (int i = 0; i < size(); i++) {
        if (isOne(std::pow(std::abs((*this)(i)),2))) {
          return true;
        }
      }
      return false;
    }

    inline Qubit& measure()
    {
      int n = size();
      std::vector<fptype> p(n);

      for (int i = 0; i < n; i++) {
        p[i] = std::pow(std::abs((*this)(i)),2);
        if (isOne(p[i])) {
          return *this;
        }
      }

      fptype s = 0., r = fptype(std::rand())/RAND_MAX;
      for (int j = 0; j < n; j++) {
        s += p[j];
        if (s >= r) {
          field c = (*this)(j)/std::abs((*this)(j));
          setZero();
          (*this)(j) = c;
          return *this;
        }
      }
      return *this;
    }

    inline Qubit& measurePartial(const int p)
    {
      int n = log2(size());
      int q = n-p;
      int pn = std::pow(2,p);
      int qn = std::pow(2,q);
      std::vector<fptype> pj(pn);

      for (int j = 0; j < pn; j++) {
        pj[j] = 0.;
        for (int r = 0; r < qn; r++) {
          pj[j] += std::pow(std::abs((*this)(j*qn+r)),2);
        }
      }

      fptype s = 0., t = fptype(std::rand())/RAND_MAX;
      for (int j = 0; j < pn; j++) {
        s += pj[j];
        if (s >= t) {
          Qubit rq(qn);
          for (int r = 0; r < qn; r++) {
              rq += (*this)(j*qn+r)/std::sqrt(pj[j])*Qubit(r,q);
          }
          *this = Qubit(j,p).tensorDot(rq);
          return *this;
        }
      }
      return *this;
    }
};

} // namespace QuCoSi

#endif // QUCOSI_QUBIT_H

// vim: shiftwidth=2 textwidth=78

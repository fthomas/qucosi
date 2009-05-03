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

#ifndef QUCOSI_GATE_H
#define QUCOSI_GATE_H

#include <Vector.h>

namespace QuCoSi {

typedef Eigen::Matrix<field, Eigen::Dynamic, Eigen::Dynamic> MatrixXc;

class Gate : public MatrixXc
{
  public:
    inline Gate() : MatrixXc() {}

    inline Gate(const int r, const int c) : MatrixXc(r, c) {}

    inline Gate& operator=(const MatrixXc& m)
    {
      MatrixXc::operator=(m);
      return *this;
    }

    inline Gate tensorDot(const Gate& m) const
    {
      int c1, r1, c2, r2;
      r1 = rows();
      c1 = cols();
      r2 = m.rows();
      c2 = m.cols();
      Gate x(r1*r2, c1*c2);

      for (int c = 0; c < c1; c++) {
        for (int r = 0; r < r1; r++) {
          x.block(r*r2, c*c2, r2, c2) = (*this)(r,c)*m;
        }
      }
      return x;
    }

    inline Gate& tensorDotSet(const Gate& m)
    {
      if (size() > 0 && m.size() > 0) {
        *this = tensorDot(m);
      }
      return *this;
    }

    inline Gate tensorPow(const int n) const
    {
      Gate x = *this;
      for (int i = 1; i < n; i++) {
        x.tensorDotSet(*this);
      }
      return x;
    }

    inline Gate& tensorPowSet(const int n)
    {
      *this = tensorPow(n);
      return *this;
    }

    inline Gate extendToPos(const int k, const int n) const
    {
      int l = n-k-log2(rows());
      Gate id, x = *this;

      if (l > 0) {
        id.resize(std::pow(2,l), std::pow(2,l));
        id.setIdentity();
        x = id.tensorDot(x);
      }
      if (l >= 0 && k > 0) {
        id.resize(std::pow(2,k), std::pow(2,k));
        id.setIdentity();
        x.tensorDotSet(id);
      }
      return x;
    }
};

} // namespace QuCoSi

#endif // QUCOSI_GATE_H

// vim: shiftwidth=2 textwidth=78

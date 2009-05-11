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

#include <vector>

#include <Vector.h>

namespace QuCoSi {

typedef Eigen::Matrix<field, Eigen::Dynamic, Eigen::Dynamic> MatrixXc;

class Gate : public MatrixXc
{
  public:
    inline Gate() : MatrixXc() {}

    inline Gate(const int r, const int c) : MatrixXc(r,c) {}

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

    /** \brief \b X gate (NOT gate)
      *
      * \f[\mathbf{X} =
      *   \left(\begin{array}{cc}
      *     0 & 1\\
      *     1 & 0
      *   \end{array}\right)
      * \f]
      */
    inline Gate& XGate()
    {
      resize(2,2);
      setZero();
      (*this)(0,1) = field(1,0);
      (*this)(1,0) = field(1,0);
      return *this;
    }

    /** \brief \b Y gate
      *
      * \f[\mathbf{Y} =
      *   \left(\begin{array}{cc}
      *     0 & -i\\
      *     i & 0
      *   \end{array}\right)
      * \f]
      */
    inline Gate& YGate()
    {
      resize(2,2);
      setZero();
      (*this)(0,1) = field(0,-1);
      (*this)(1,0) = field(0,1);
      return *this;
    }

    /** \brief \b Z gate
      *
      * \f[\mathbf{Z} =
      *   \left(\begin{array}{cc}
      *     1 & 0\\
      *     0 & -1
      *   \end{array}\right)
      * \f]
      */
    inline Gate& ZGate()
    {
      resize(2,2);
      setZero();
      (*this)(0,0) = field(1,0);
      (*this)(1,1) = field(-1,0);
      return *this;
    }

    /** \brief \b H gate (Hadamard gate)
      *
      * \f[\mathbf{H} = \frac{1}{\sqrt{2}}
      *   \left(\begin{array}{cc}
      *     1 & 1\\
      *     1 & -1
      *   \end{array}\right)
      * \f]
      */
    inline Gate& HGate()
    {
      resize(2,2);
      fptype c = std::sqrt(0.5);
      (*this)(0,0) = field(c,0);
      (*this)(0,1) = field(c,0);
      (*this)(1,0) = field(c,0);
      (*this)(1,1) = field(-c,0);
      return *this;
    }

    /** \brief \b P gate (phase gate)
      *
      * \f[\mathbf{P} =
      *   \left(\begin{array}{cc}
      *     1 & 0\\
      *     0 & i
      *   \end{array}\right)
      * \f]
      */
    inline Gate& PGate()
    {
      resize(2,2);
      setZero();
      (*this)(0,0) = field(1,0);
      (*this)(1,1) = field(0,1);
      return *this;
    }

    /** \brief \b T gate (pi/4 phase shift gate)
      *
      * \f[\mathbf{T} =
      *   \left(\begin{array}{cc}
      *     1 & 0\\
      *     0 & e^{\frac{i \pi}{4}}
      *   \end{array}\right)
      * \f]
      */
    inline Gate& TGate()
    {
      resize(2,2);
      setZero();
      fptype c = std::sqrt(0.5);
      (*this)(0,0) = field(1,0);
      (*this)(1,1) = field(c,c);
      return *this;
    }

    /** \brief \b CNOT gate (controlled NOT gate)
      *
      * \f[\mathbf{CNOT} =
      *   \left(\begin{array}{cccc}
      *     1 & 0 & 0 & 0\\
      *     0 & 1 & 0 & 0\\
      *     0 & 0 & 0 & 1\\
      *     0 & 0 & 1 & 0
      *   \end{array}\right)
      * \f]
      */
    inline Gate& CNOTGate()
    {
      resize(4,4);
      setZero();
      (*this)(0,0) = field(1,0);
      (*this)(1,1) = field(1,0);
      (*this)(2,3) = field(1,0);
      (*this)(3,2) = field(1,0);
      return *this;
    }

    /** \brief \b SWAP gate
      *
      * \f[\mathbf{SWAP} =
      *   \left(\begin{array}{cccc}
      *     1 & 0 & 0 & 0\\
      *     0 & 0 & 1 & 0\\
      *     0 & 1 & 0 & 0\\
      *     0 & 0 & 0 & 1
      *   \end{array}\right)
      * \f]
      */
    inline Gate& SWAPGate()
    {
      resize(4,4);
      setZero();
      (*this)(0,0) = field(1,0);
      (*this)(1,2) = field(1,0);
      (*this)(2,1) = field(1,0);
      (*this)(3,3) = field(1,0);
      return *this;
    }

    inline Gate& UfGate(const Vector& v, const std::vector<int>& f)
    {
      int s = v.size();
      resize(s,s);
      setIdentity();

      int x = 0, y = 0, z = 0, w = 0;
      for ( ; z < s; z++) {
        if (v(z) != field(0,0)) {
          x = z/2;
          y = z%2;
          w = z-y;

          if (y != (f.at(x)+y)%2) {
            (*this)(w,w) = field(0,0);
            (*this)(w+1,w+1) = field(0,0);
            (*this)(w+1,w) = field(1,0);
            (*this)(w,w+1) = field(1,0);
          }
        }
      }
      return *this;
    }

    inline Gate& UfGate2(const Vector& v, const std::vector<int>& f)
    {
      int s = v.size();
      resize(s,s);
      setIdentity();

      int x = 0, y = 0, z = 0;
      for ( ; z < s; z++) {
        if (v(z) != field(0,0)) {
          break;
        }
      }
      x = z/2;
      y = z%2;
      z -= y;

      if (z < s && y != (f.at(x)+y)%2) {
        (*this)(z,z) = field(0,0);
        (*this)(z+1,z+1) = field(0,0);
        (*this)(z+1,z) = field(1,0);
        (*this)(z,z+1) = field(1,0);
      }
      return *this;
    }
};

} // namespace QuCoSi

#endif // QUCOSI_GATE_H

// vim: shiftwidth=2 textwidth=78

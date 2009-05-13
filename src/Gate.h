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

#include <bitset>
#include <cmath>
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

    inline Gate applyToPos(const int k, const int n) const
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
      *   = \mathbf{R}(2)
      * \f]
      *
      * \sa RGate()
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
      *   = \frac{1}{\sqrt{2}} (\mathbf{X} + \mathbf{Z})
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
      *   = \mathbf{R}(4)
      * \f]
      *
      * \sa RGate()
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
      *     0 & e^{\pi i/4}
      *   \end{array}\right)
      *   = \mathbf{R}(8)
      * \f]
      *
      * \sa RGate()
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

    /** \brief <b>R</b>(\e k) gate (general phase shift gate)
      *
      * \f[\mathbf{R}(k) =
      *   \left(\begin{array}{cc}
      *     1 & 0\\
      *     0 & e^{2 \pi i/k}
      *   \end{array}\right)
      * \f]
      */
    inline Gate& RGate(const int k)
    {
      resize(2,2);
      setZero();
      (*this)(0,0) = field(1,0);
      (*this)(1,1) = std::exp(2*M_PI/k*field(0,1));
      return *this;
    }

    /** \brief <b>C</b>(\e U) gate (controlled \e U gate)
      *
      * \f[\mathbf{C}(U) =
      *   \left(\begin{array}{cccc}
      *     1 & 0 & 0 & 0\\
      *     0 & 1 & 0 & 0\\
      *     0 & 0 & U_{11} & U_{12}\\
      *     0 & 0 & U_{21} & U_{22}
      *   \end{array}\right)
      * \f]
      */
    inline Gate& CGate(const Gate& U)
    {
      resize(4,4);
      setZero();
      (*this)(0,0) = field(1,0);
      (*this)(1,1) = field(1,0);
      block(2,2,2,2) = U;
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
      *   = \mathbf{C}(\mathbf{X})
      * \f]
      *
      * \sa CGate()
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

    /** \brief Swaps Qubits according to the permutation \p sigma
      *
      * This methods constructs a tensor permutation matrix that permutes
      * qubits according to the permutation \p sigma.
      *
      * \note
      * The implementation of this method is based on proposition 6.2 in the
      * paper arXiv:math/0508053v2 by Rakotonirina Christian. Read the paper
      * if you are interested why this method constructs a permutation
      * matrix. This implementation takes advantage of the fact that the
      * dimension of single qubits is 2 so that the multiple row and column
      * indices \f$i_1 \ldots i_k\f$ and \f$j_1 \ldots j_k\f$ can be obtained
      * from the row and column indices of the permutation matrix with
      * <tt>std::bitset</tt>s.
      *
      * \sa http://arxiv.org/abs/math/0508053
      */
    inline Gate& SGate(const std::vector<int>& sigma)
    {
      int n = sigma.size();
      int dim = std::pow(2,n);
      resize(dim,dim);
      setZero();

      for (int c = 0; c < dim; c++) {
        for (int r = 0; r < dim; r++) {
          // Create bitsets of the row and column index. Since the dimension
          // of single qubits is 2, these bitsets can be interpreted as the
          // multiple row and colum indices.
          std::bitset< std::numeric_limits<unsigned long>::digits >
            br(r), bc(c), bcp(0);

          // Permute the bitset bc using the permutation sigma.
          for (int i = 0; i < n; i++) {
            bcp[n-1-i] = bc[n-1-sigma[i]];
          }

          // If the bitset of the row index is equal to the permuted bitset of
          // the column index, set a 1 at this position. This equals the
          // multiple Kronecker deltas.
          if (br == bcp) {
            (*this)(r,c) = field(1,0);
          }
        }
      }
      return *this;
    }

    /** \brief <b>S</b>\e pq gate - swaps the Qubits at position \e p and \e q
      *
      * \sa SGate()
      */
    inline Gate& SpqGate(const int p, const int q, const int n)
    {
      std::vector<int> sigma(n);
      for (int i = 0; i < n; i++) {
        sigma[i] = i;
      }
      sigma[p] = q;
      sigma[q] = p;
      SGate(sigma);
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
      *   = \mathbf{S}_{10}
      * \f]
      *
      * \sa SpqGate(), SGate()
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

      int x = 0, y = 0, z = 0;
      for (int i = 0; i < s; i++) {
        if (v(i) != field(0,0)) {
          x = i/2;
          y = i%2;
          z = i-y;

          if (y != (f.at(x)+y)%2) {
            (*this)(z,z) = field(0,0);
            (*this)(z+1,z+1) = field(0,0);
            (*this)(z+1,z) = field(1,0);
            (*this)(z,z+1) = field(1,0);
          }
        }
      }
      return *this;
    }

    /** \brief \b F gate (quantum Fourier transform)
      *
      * \f[\mathbf{F_N} = \frac{1}{\sqrt{N}}
      *   \left(\begin{array}{cccccc}
      *     1 & 1 & 1 & 1 & \cdots & 1 \\
      *     1 & \omega^1 & \omega^2 & \omega^3 & \cdots & \omega^{N-1}\\
      *     1 & \omega^2 & \omega^4 & \omega^6 & \cdots & \omega^{2(N-1)}\\
      *     1 & \omega^3 & \omega^6 & \omega^9 & \cdots & \omega^{3(N-1)}\\
      *     \vdots & \vdots & \vdots & \vdots & \ddots & \vdots \\
      *     1 & \omega^{N-1} & \omega^{2(N-1)} & \omega^{3(N-1)} & \cdots &
      *       \omega^{(N-1)(N-1)}
      *   \end{array}\right),
      *   \quad \omega = \exp\left(2 \pi i/N\right)
      * \f]
      */
    inline Gate& FGate(const int N)
    {
      resize(N,N);
      for (int x = 0; x < N; x++) {
        for (int y = x; y < N; y++) {
          (*this)(x,y) = std::exp(2*M_PI*x*y/N*field(0,1));
        }
      }
      *this *= std::sqrt(1./N);

      MatrixXc t = transpose();
      t.diagonal().setZero();
      *this += t;

      return *this;
    }
};

} // namespace QuCoSi

#endif // QUCOSI_GATE_H

// vim: shiftwidth=2 textwidth=78

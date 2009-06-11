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
#include <limits>
#include <vector>

#include <Aux.h>
#include <Vector.h>

namespace QuCoSi {

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

    /** \brief \b X gate (NOT gate, Pauli \b X matrix)
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

    /** \brief \b Y gate (Pauli \b Y matrix)
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

    /** \brief \b Z gate (Pauli \b Z matrix)
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

    /** \brief <b>R</b>(\p k) gate (general phase shift gate)
      *
      * \f[\mathbf{R}(k) =
      *   \left(\begin{array}{cc}
      *     1 & 0\\
      *     0 & e^{2 \pi i/k}
      *   \end{array}\right)
      * \f]
      */
    inline Gate& RGate(const fptype k)
    {
      resize(2,2);
      setZero();
      (*this)(0,0) = field(1,0);
      (*this)(1,1) = std::exp(2*M_PI/k*field(0,1));
      return *this;
    }

    /** \brief \b I gate (identity gate)
      *
      * \f[\mathbf{I} =
      *   \left(\begin{array}{cc}
      *     1 & 0\\
      *     0 & 1
      *   \end{array}\right)
      *   = \mathbf{R}(1)
      * \f]
      */
    inline Gate& IGate()
    {
      resize(2,2);
      setIdentity();
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
      *   = \mathbf{C}_{102}(\mathbf{X})
      * \f]
      *
      * \sa CGate(), XGate()
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

    /** \brief \b CCNOT gate (Toffoli gate, controlled \b CNOT gate)
      *
      * \f[\mathbf{CCNOT} =
      *   \left(\begin{array}{cccccccc}
      *     1 & 0 & 0 & 0 & 0 & 0 & 0 & 0\\
      *     0 & 1 & 0 & 0 & 0 & 0 & 0 & 0\\
      *     0 & 0 & 1 & 0 & 0 & 0 & 0 & 0\\
      *     0 & 0 & 0 & 1 & 0 & 0 & 0 & 0\\
      *     0 & 0 & 0 & 0 & 1 & 0 & 0 & 0\\
      *     0 & 0 & 0 & 0 & 0 & 1 & 0 & 0\\
      *     0 & 0 & 0 & 0 & 0 & 0 & 0 & 1\\
      *     0 & 0 & 0 & 0 & 0 & 0 & 1 & 0
      *   \end{array}\right)
      *   = \mathbf{C}_{103}(\mathbf{CNOT})
      * \f]
      *
      * \sa CGate(), CNOTGate()
      */
    inline Gate& CCNOTGate()
    {
      resize(8,8);
      setIdentity();
      (*this)(6,6) = field(0,0);
      (*this)(7,7) = field(0,0);
      (*this)(7,6) = field(1,0);
      (*this)(6,7) = field(1,0);
      return *this;
    }

    /** \brief \b CSWAP gate (Fredkin gate, controlled \b SWAP gate)
      *
      * \f[\mathbf{CSWAP} =
      *   \left(\begin{array}{cccccccc}
      *     1 & 0 & 0 & 0 & 0 & 0 & 0 & 0\\
      *     0 & 1 & 0 & 0 & 0 & 0 & 0 & 0\\
      *     0 & 0 & 1 & 0 & 0 & 0 & 0 & 0\\
      *     0 & 0 & 0 & 1 & 0 & 0 & 0 & 0\\
      *     0 & 0 & 0 & 0 & 1 & 0 & 0 & 0\\
      *     0 & 0 & 0 & 0 & 0 & 0 & 1 & 0\\
      *     0 & 0 & 0 & 0 & 0 & 1 & 0 & 0\\
      *     0 & 0 & 0 & 0 & 0 & 0 & 0 & 1
      *   \end{array}\right)
      *   = \mathbf{C}_{103}(\mathbf{SWAP})
      * \f]
      *
      * \sa CGate(), SWAPGate()
      */
    inline Gate& CSWAPGate()
    {
      resize(8,8);
      setIdentity();
      (*this)(5,5) = field(0,0);
      (*this)(6,6) = field(0,0);
      (*this)(6,5) = field(1,0);
      (*this)(5,6) = field(1,0);
      return *this;
    }

    /** \brief <b>C</b><sub>\p tcn</sub>(\p U) gate (controlled \p U gate)
      */
    inline Gate& CGate(const int t, const int c, const int n, const Gate& U)
    {
      // Construct the controlled U gate with the first qubit as control and
      // the second qubit as target.
      int d = U.rows();
      Gate cu(2*d,2*d);
      cu.setIdentity();
      cu.block(d,d,d,d) = U;

      // Extend the controlled U gate for n qubits if it does not already
      // operate on n qubits.
      if (std::pow(2,n) > 2*d) {
        cu = cu.applyToPos(0,n);
      }

      // Construct a cyclic permutation that displaces t to position 1.
      std::vector<int> sigma(n);
      for (int i = 0; i < n; i++) {
        sigma[(n-t+i+1)%n] = i;
      }
      // Permute c to position 0.
      int newc = (n-t+c+1)%n, tmp = sigma[0];
      sigma[0] = sigma[newc];
      sigma[newc] = tmp;

      // Construct an S gate from the permutation sigma.
      Gate s;
      s.SGate(sigma);

      // Now combine the S and the controlled U gate.
      *this = s*cu*s;
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
      *   = \mathbf{S}_{102} = \mathbf{S}_{012}
      * \f]
      *
      * \sa SGate()
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

    /** \brief <b>S</b><sub>\p pqn</sub> gate
      *
      * This method constructs a \f$2^n \times 2^n\f$ tensor permutation
      * matrix that permutes the \f$p^{\mathrm{th}}\f$ and \f$q^{\mathrm{th}}
      * \f$ qubits in a tensor product of \f$n\f$ qubits.
      */
    inline Gate& SGate(const int p, const int q, const int n)
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

    /** \brief <b>S</b>(\p sigma) gate
      *
      * This method constructs a tensor permutation matrix that permutes
      * qubits according to the permutation \p sigma. For example, this
      * matrix acts on the tensor product of the qubits \f$ q_0,\ q_1,\
      * \ldots,\ q_k\f$ as follows:
      * \f[
      *   \mathbf{S}(\sigma) \left(q_0 \otimes q_1 \otimes \ldots \otimes
      *     q_k\right) = q_{\sigma(0)} \otimes q_{\sigma(1)} \otimes \ldots
      *     \otimes q_{\sigma(k)}
      * \f]
      *
      * \note The implementation of this method is based on proposition 6.2
      * in the paper arXiv:math/0508053v2 by Rakotonirina Christian. It takes
      * advantage of the fact that the dimension of single qubits is 2 so
      * that the multiple row and column indices \f$i_1 \ldots i_k\f$ and
      * \f$j_1 \ldots j_k\f$ can be obtained from the row and column indices
      * of the permutation matrix with <tt>std::bitset</tt>s.
      *
      * \sa http://arxiv.org/abs/math/0508053
      */
    inline Gate& SGate(const std::vector<int>& sigma)
    {
      int n = sigma.size();
      int dim = std::pow(2,n);
      resize(dim,dim);
      setZero();

      // The first and last component of any tensor stays unchanged for all
      // tensor permutations, therefore we ignore the outermost rows and
      // columns.
      (*this)(0,0) = field(1,0);
      (*this)(dim-1,dim-1) = field(1,0);

      for (int c = 1; c < dim-1; c++) {
        for (int r = 1; r < dim-1; r++) {
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
            break;
          }
        }
      }
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

    /** \brief <b>F</b><sub>\p n</sub> gate (quantum Fourier transform)
      *
      * \f[\mathbf{F_n} = 2^{-n/2}
      *   \left(\begin{array}{cccccc}
      *     1 & 1 & 1 & 1 & \cdots & 1 \\
      *     1 & \omega^1 & \omega^2 & \omega^3 & \cdots & \omega^{2^n-1}\\
      *     1 & \omega^2 & \omega^4 & \omega^6 & \cdots & \omega^{2(2^n-1)}\\
      *     1 & \omega^3 & \omega^6 & \omega^9 & \cdots & \omega^{3(2^n-1)}\\
      *     \vdots & \vdots & \vdots & \vdots & \ddots & \vdots \\
      *     1 & \omega^{2^n-1} & \omega^{2(2^n-1)} & \omega^{3(2^n-1)} &
      *       \cdots & \omega^{(2^n-1)^2}
      *   \end{array}\right),
      *   \quad \omega = \exp\left(2 \pi i/2^n\right)
      * \f]
      */
    inline Gate& FGate(const int n)
    {
      int N = std::pow(2,n);
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

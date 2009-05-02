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

#include <Eigen/Array>
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

/** \class Vector
  *
  * \brief Dynamic size vector of complex numbers
  *
  * The Vector class is the base for all vectors and qubits used in QuCoSi. It
  * is of dynamic size and uses complex numbers. Besides the standard methods
  * that are inherited from the Eigen base class it offers some convenient
  * methods like isNormalized() and randomize(). The most important feature of
  * this class is the canonical tensor product tensorDot() and tensorDotSet().
  *
  * \sa Qubit
  */
class Vector : public VectorXc
{
  public:
    inline Vector() : VectorXc(2)
    {
      (*this)[0] = field(1,0);
      (*this)[1] = field(0,0);
    }

    inline Vector(const unsigned dim) : VectorXc(dim) {}

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

    /** \brief Calculates the tensor product of this vector with \a v
      *
      * This method calculates the canonical tensor product of this vector
      * with Vector \a v. The canonical tensor product \f$x \otimes y\f$ of
      * the vectors \f$x \in K^n\f$ and \f$y \in K^m\f$ is defined as:
      * \f[
      *   \left(\begin{array}{c}
      *     x_1 \\ x_2 \\ \vdots \\ x_n
      *   \end{array}\right)
      *   \otimes
      *   \left(\begin{array}{c}
      *     y_1 \\ y_2 \\ \vdots \\ y_m
      *   \end{array}\right)
      *   =
      *   \left(\begin{array}{c}
      *     x_1 y_1 \\ \vdots \\ x_1 y_m \\
      *     x_2 y_1 \\ \vdots \\ x_2 y_m \\
      *     \vdots \\ x_n y_m
      *   \end{array}\right)
      * \f]
      *
      * \param v The right hand side operand of the tensor product
      * \return The canonical tensor product of this vector with Vector \a v
      */
    inline Vector tensorDot(const Vector& v) const
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

    /** \brief Sets the tensor product of this vector and \a v as this vector
      *
      * Sets the result of the tensor product of this vector and Vector \a v
      * as this vector. For two vectors \c x and \c y \code x.tensorDotSet(y)
      * \endcode is practically identical to \code x = x.tensorDot(y) \endcode
      *
      * \sa tensorDot()
      */
    inline Vector& tensorDotSet(const Vector& v)
    {
      if (rows() > 0 && v.rows() > 0) {
        *this = tensorDot(v);
      }
      return *this;
    }
};

} // namespace QuCoSi

#endif // QUCOSI_VECTOR_H

// vim: shiftwidth=2 textwidth=78

/*
   <one line to give the program's name and a brief idea of what it does.>
   Copyright © 2009 Frank S. Thomas <frank@thomas-alfeld.de>

   This program is free software: you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation, either version 3 of the License, or
   (at your option) any later version.

   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/

#ifndef QUBIT_H
#define QUBIT_H

#include <cstdlib>
#include <limits>
#include <stdexcept>

#include <Eigen/Array>
#include <Eigen/Core>

typedef Eigen::Vector2cf Vector2D;
typedef std::complex<float> field;

class Qubit : public Vector2D {
  protected:
    Vector2D m_std_base[2];

  private:
    inline bool isOne(const float x) const
    {
      return std::abs(x-1.) < std::numeric_limits<float>::epsilon();
    }

  public:
    inline Qubit() : Vector2D()
    {
      (*this)[0] = field(1,0);
      (*this)[1] = field(0,0);
      setStdBase();
    }

    inline Qubit(const field c0, const field c1) : Vector2D()
    {
      (*this)[0] = c0;
      (*this)[1] = c1;
      setStdBase();
    }

    inline Qubit& operator=(const Vector2D& v)
    {
      Vector2D::operator=(v);
      return *this;
    }

    inline void setStdBase()
    {
      Vector2D v0, v1;

      v0[0] = field(1, 0);
      v0[1] = field(0, 0);

      v1[0] = field(0, 0);
      v1[1] = field(0, 1);

      m_std_base[0] = v0;
      m_std_base[1] = v1;
    }

    inline void setStdBase(const Vector2D& baseVector0,
                           const Vector2D& baseVector1)
    {
      if (!isOne(baseVector0.norm()) || !isOne(baseVector1.norm())) {
        throw std::logic_error("Qubit::setStdBase(): "
          "one of the alleged base vectors is not normalized");
      }
      if (!baseVector0.isOrthogonal(baseVector1)) {
        throw std::logic_error("Qubit::setStdBase(): "
          "the alleged base vectors are not orthogonal");
      }
      m_std_base[0] = baseVector0;
      m_std_base[1] = baseVector1;
    }

    inline field coefficient(const unsigned index) const
    {
      return m_std_base[index%2].dot(*this);
    }

    inline bool isNormalized() const
    {
      return isOne(norm());
    }

    inline bool isPureState() const
    {
      return isPureState(m_std_base[0], m_std_base[1]);
    }

    inline bool isPureState(const Vector2D& baseVector0,
                            const Vector2D& baseVector1) const
    {
      field c[2];
      c[0] = baseVector0.dot(*this);
      c[1] = baseVector1.dot(*this);

      for (int i = 0; i <= 1; i++) {
        if (std::pow(std::abs(c[i]),2) == 0) {
          return true;
        }
      }
      return false;
    }

    inline Qubit& setRandom()
    {
      Vector2D::setRandom().normalize();
      return *this;
    }

    inline Qubit& measure()
    {
      return measure(m_std_base[0], m_std_base[1]);
    }

    inline Qubit& measure(const Vector2D& baseVector0,
                          const Vector2D& baseVector1)
    {
      field c[2];
      c[0] = baseVector0.dot(*this);
      c[1] = baseVector1.dot(*this);

      double p[2];
      p[0] = std::pow(std::abs(c[0]),2);
      p[1] = std::pow(std::abs(c[1]),2);

      if (p[0] == 0 || p[1] == 0) {
        return *this;
      }

      double w = double(std::rand())/RAND_MAX;
      if (w <= p[0]) {
        *this = c[0]/std::abs(c[0]) * baseVector0;
      }
      else {
        *this = c[1]/std::abs(c[1]) * baseVector1;
      }
      return *this;
    }
};

#endif // QUBIT_H

// vim: shiftwidth=2

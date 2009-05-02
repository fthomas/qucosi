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
#include <stdexcept>

#include <Basis.h>
#include <Vector.h>

namespace QuCoSi {

class Qubit : public Vector {
  private:
    Basis m_std_basis;

  public:
    inline Qubit() : Vector(2)
    {
      m_std_basis = Basis(2);
    }

    inline Qubit(const field& c0, const field& c1) : Vector(c0, c1)
    {
      m_std_basis = Basis(2);
    }

    inline Qubit& operator=(const Vector& v)
    {
      Vector::operator=(v);
      return *this;
    }

    inline Basis getStdBasis() const
    {
      return m_std_basis;
    }

    inline void setStdBasis()
    {
      m_std_basis.setNaturalBasis(rows());
    }

    inline void setStdBasis(const Basis& b)
    {
      if (b.size() == rows()) {
        m_std_basis = b;
      }
    }

    inline Qubit& tensorDotSet(const Qubit& q)
    {
      Vector::tensorDotSet(q);
      m_std_basis.tensorDotSet(q.getStdBasis());
      return *this;
    }

    inline field coefficient(const unsigned index) const
    {
      return m_std_basis[index%rows()].dot(*this);
    }

    inline bool isPureState() const
    {
      return isPureState(m_std_basis);
    }

    inline bool isPureState(const Basis& b) const
    {
      field c;
      for (int i = 0; i < rows(); i++) {
        c = b.at(i).dot(*this);
        if (isOne(std::pow(std::abs(c),2))) {
          return true;
        }
      }
      return false;
    }

    inline Qubit& measure()
    {
      return measure(m_std_basis);
    }

    inline Qubit& measure(const Basis& b)
    {
      int n = rows();
      std::vector<field> c(n);
      std::vector<fptype> p(n);

      for (int i = 0; i < n; i++) {
        c[i] = b.at(i).dot(*this);
        p[i] = std::pow(std::abs(c[i]),2);
        if (isOne(p[i])) {
          return *this;
        }
      }

      fptype s = 0., r = fptype(std::rand())/RAND_MAX;
      for (int j = 0; j < n; j++) {
        s += p[j];
        if (s >= r) {
          *this = b.at(j);
          *this *= c[j]/std::abs(c[j]);
          return *this;
        }
      }
      return *this;
    }
};

} // namespace QuCoSi

#endif // QUCOSI_QUBIT_H

// vim: shiftwidth=2 textwidth=78

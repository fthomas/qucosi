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

#ifndef QUCOSI_GATETEST_H
#define QUCOSI_GATETEST_H

#include <cppunit/TestFixture.h>
#include <cppunit/extensions/HelperMacros.h>

#include <Gate.h>

namespace QuCoSi {

class GateTest : public CppUnit::TestFixture
{
  CPPUNIT_TEST_SUITE(GateTest);
  CPPUNIT_TEST(testTensorPow);
  CPPUNIT_TEST_SUITE_END();

  public:
    void setUp() {}

    void tearDown() {}

    void testTensorPow()
    {
      Gate h, h_byhand;

      h.HGate().tensorPowSet(2);
      h_byhand.resize(4,4);
      h_byhand << 1,  1,  1,  1,
                  1, -1,  1, -1,
                  1,  1, -1, -1,
                  1, -1, -1,  1;
      h_byhand *= 0.5;
      CPPUNIT_ASSERT( h.isApprox(h_byhand) );

      h.HGate().tensorPowSet(3);
      h_byhand.resize(8,8);
      h_byhand << 1,  1,  1,  1,  1,  1,  1,  1,
                  1, -1,  1, -1,  1, -1,  1, -1,
                  1,  1, -1, -1,  1,  1, -1, -1,
                  1, -1, -1,  1,  1, -1, -1,  1,
                  1,  1,  1,  1, -1, -1, -1, -1,
                  1, -1,  1, -1, -1,  1, -1,  1,
                  1,  1, -1, -1, -1, -1,  1,  1,
                  1, -1, -1,  1, -1,  1,  1, -1;
      h_byhand *= std::sqrt(1./8);
      CPPUNIT_ASSERT( h.isApprox(h_byhand) );

      Gate i, i_byhand;
      i.resize(2,2);
      i.setIdentity();

      for (int n = 2; n <= 8; n++) {
        int p = std::pow(2,n);
        i_byhand.resize(p,p);
        i_byhand.setIdentity();
        CPPUNIT_ASSERT( i_byhand.isApprox(i.tensorPow(n)) );
      }
    }
};

} // namespace QuCoSi

#endif // QUCOSI_GATETEST_H

// vim: shiftwidth=2 textwidth=78

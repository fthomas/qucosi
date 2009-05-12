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
  CPPUNIT_TEST(testApplyToPos);
  CPPUNIT_TEST(testFGate);
  CPPUNIT_TEST_SUITE_END();

  public:
    void setUp() {}

    void tearDown() {}

    void testTensorPow()
    {
      Gate h, h_byhand;

      h.HGate().tensorPowSet(1);
      h_byhand.resize(2,2);
      h_byhand << 1,  1,
                  1, -1;
      h_byhand *= std::sqrt(1./2);
      CPPUNIT_ASSERT( h.isApprox(h_byhand) );
      CPPUNIT_ASSERT( h.isUnitary() );
      CPPUNIT_ASSERT( h_byhand.isUnitary() );

      h.HGate().tensorPowSet(2);
      h_byhand.resize(4,4);
      h_byhand << 1,  1,  1,  1,
                  1, -1,  1, -1,
                  1,  1, -1, -1,
                  1, -1, -1,  1;
      h_byhand *= std::sqrt(1./4);
      CPPUNIT_ASSERT( h.isApprox(h_byhand) );
      CPPUNIT_ASSERT( h.isUnitary() );
      CPPUNIT_ASSERT( h_byhand.isUnitary() );

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
      CPPUNIT_ASSERT( h.isUnitary() );
      CPPUNIT_ASSERT( h_byhand.isUnitary() );

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

    void testApplyToPos()
    {
      Gate h, h0, h1;
      h.HGate();
      h0.resize(4,4);
      h1.resize(4,4);

      h0 << 1,  1,  0,  0,
            1, -1,  0,  0,
            0,  0,  1,  1,
            0,  0,  1, -1;
      h0 *= std::sqrt(1./2);

      h1 << 1,  0,  1,  0,
            0,  1,  0,  1,
            1,  0, -1,  0,
            0,  1,  0, -1;
      h1 *= std::sqrt(1./2);

      CPPUNIT_ASSERT( h0.isApprox(h.applyToPos(0,2)) );
      CPPUNIT_ASSERT( h1.isApprox(h.applyToPos(1,2)) );
      CPPUNIT_ASSERT( h.applyToPos(0,2).isUnitary() );
      CPPUNIT_ASSERT( h.applyToPos(1,2).isUnitary() );
    }

    void testFGate()
    {
      field i = field(0,1);
      Qubit q0(0,4),   q1(1,4),   q2(2,4),   q3(3,4),
            q4(4,4),   q5(5,4),   q6(6,4),   q7(7,4),
            q8(8,4),   q9(9,4),   q10(10,4), q11(11,4),
            q12(12,4), q13(13,4), q14(14,4), q15(15,4);

      Gate f;
      f.FGate(16);

      CPPUNIT_ASSERT( (f*(0.5*(q0 + q4 + q8 + q12))).isApprox(
                          0.5*(q0 + q4 + q8 + q12)) );

      CPPUNIT_ASSERT( (f*(0.5*(q1 + q5 + q9 + q13))).isApprox(
                          0.5*(q0 + i*q4 - q8 - i*q12)) );

      CPPUNIT_ASSERT( (f*(0.5*(q2 + q6 + q10 + q14))).isApprox(
                          0.5*(q0 - q4 + q8  - q12)) );

      CPPUNIT_ASSERT( (f*(0.5*(q3 + q7 + q11 + q15))).isApprox(
                          0.5*(q0 - i*q4 - q8 + i*q12)) );
    }
};

} // namespace QuCoSi

#endif // QUCOSI_GATETEST_H

// vim: shiftwidth=2 textwidth=78

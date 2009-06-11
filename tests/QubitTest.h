// QuCoSi - Quantum Computer Simulation
// Copyright © 2009 Frank S. Thomas <frank@blue-dwarf.de>
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

#ifndef QUCOSI_QUBITTEST_H
#define QUCOSI_QUBITTEST_H

#include <cstdlib>
#include <ctime>

#include <cppunit/TestFixture.h>
#include <cppunit/extensions/HelperMacros.h>

#include <Qubit.h>

namespace QuCoSi {

class QubitTest : public CppUnit::TestFixture
{
  CPPUNIT_TEST_SUITE(QubitTest);
  CPPUNIT_TEST(testMeasure);
  CPPUNIT_TEST(testMeasurePartial);
  CPPUNIT_TEST_SUITE_END();

  public:
    void setUp()
    {
      std::srand((unsigned)std::time(NULL) + (unsigned)std::clock());
    }

    void tearDown() {}

    void testMeasure()
    {
      Qubit v, x, q0(0,2), q1(1,2), q2(2,2);
      x = std::sqrt(0.5)*q0 - 0.5*q1 + 0.5*q2;

      int n = 0, r0 = 0, r1 = 0, r2 = 0;
      for ( ; n < 1000; n++) {
        v = x;
        v.measure();
        if (v == q0) r0++;
        if (v == -1*q1) r1++;
        if (v == q2) r2++;
      }
      CPPUNIT_ASSERT( r0 + r1 + r2 == n );
      CPPUNIT_ASSERT( r0 > 450 && r0 < 550);
      CPPUNIT_ASSERT( r1 > 200 && r1 < 300);
      CPPUNIT_ASSERT( r2 > 200 && r2 < 300);
    }

    void testMeasurePartial()
    {
      Qubit b, r1, r2, q0(0,2), q1(1,2), q2(2,2), q3(3,2);

      b = std::sqrt(0.5)*q0 - 0.5*q1 + 0.5*q2;
      b.measurePartial(1);
      r1 = std::sqrt(2./3.)*q0 - std::sqrt(1./3.)*q1;
      r2 = q2;
      CPPUNIT_ASSERT( b.isApprox(r1) || b.isApprox(r2) );

      b = std::sqrt(0.5)*(q0 + q3);
      b.measurePartial(1);
      r1 = q0;
      r2 = q3;
      CPPUNIT_ASSERT( b.isApprox(r1) || b.isApprox(r2) );
    }
};

} // namespace QuCoSi

#endif // QUCOSI_QUBITTEST_H

// vim: shiftwidth=2 textwidth=78

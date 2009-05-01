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

#ifndef QUCOSI_VECTORTEST_H
#define QUCOSI_VECTORTEST_H

#include <cppunit/TestFixture.h>
#include <cppunit/extensions/HelperMacros.h>

#include <Vector.h>

namespace QuCoSi {

class VectorTest : public CppUnit::TestFixture
{
  CPPUNIT_TEST_SUITE(VectorTest);
  CPPUNIT_TEST(testIsNormalized);
  CPPUNIT_TEST(testOtimes);
  CPPUNIT_TEST_SUITE_END();

  public:
    void setUp()
    {
    }

    void tearDown()
    {
    }

    void testIsNormalized()
    {
      Vector v1(field(1,0), field(0,0)),
             v2(field(1,0), field(1,0)),
             v3(field(1,0), field(0.001,0)),
             v4(field(0.7071067811865475,0), field(0,0.7071067811865475));

      CPPUNIT_ASSERT( v1.isNormalized() == true );
      CPPUNIT_ASSERT( v2.isNormalized() == false );
      CPPUNIT_ASSERT( v3.isNormalized() == false );
      CPPUNIT_ASSERT( v4.isNormalized() == true );

      CPPUNIT_ASSERT( v1.randomize().isNormalized() == true );
      CPPUNIT_ASSERT( v2.randomize().isNormalized() == true );
      CPPUNIT_ASSERT( v3.randomize().isNormalized() == true );
      CPPUNIT_ASSERT( v4.randomize().isNormalized() == true );

      v1.resize(10);
      v1.setZero();
      v1[9] = field(1,0);
      CPPUNIT_ASSERT( v1.isNormalized() == true );
    }

    void testOtimes()
    {
      Vector v1(field(1,0), field(0,0)),
             v2(field(1,0), field(0,0)),
             v3(4);

      v3[0] = field(1,0);
      CPPUNIT_ASSERT( v3 == v1.otimes(v1) );
      CPPUNIT_ASSERT( v3.otimes(v1) == v1.otimes(v1).otimes(v1) );

      v1[0] = field(2,0);
      v1[1] = field(3,0);

      v2[0] = field(5,0);
      v2[1] = field(7,0);

      v3[0] = field(10,0);
      v3[1] = field(14,0);
      v3[2] = field(15,0);
      v3[3] = field(21,0);

      CPPUNIT_ASSERT( v3 == v1.otimes(v2) );

      CPPUNIT_ASSERT( v1.otimes(v2).size() == 4 );
      CPPUNIT_ASSERT( v3.otimes(v2).size() == 8 );
      CPPUNIT_ASSERT( v3.otimes(v3).size() == 16 );
    }
};

} // namespace QuCoSi

#endif // QUCOSI_VECTORTEST_H

// vim: shiftwidth=2 textwidth=78

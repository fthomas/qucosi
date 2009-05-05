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

#ifndef QUCOSI_ALGORITHMSTEST_H
#define QUCOSI_ALGORITHMSTEST_H

#include <cstdlib>
#include <ctime>

#include <cppunit/TestFixture.h>
#include <cppunit/extensions/HelperMacros.h>

#include <Gate.h>
#include <Qubit.h>

namespace QuCoSi {

class AlgorithmsTest : public CppUnit::TestFixture
{
  CPPUNIT_TEST_SUITE(AlgorithmsTest);
  CPPUNIT_TEST(testCoinFlipping);
  CPPUNIT_TEST_SUITE_END();

  public:
    void setUp()
    {
      std::srand((unsigned)std::time(NULL) + (unsigned)std::clock());
    }

    void tearDown() {}

    void testCoinFlipping()
    {
      // First test the classical coin.
      fptype coin;
      int n = 0, heads = 0, tails = 0;

      for ( ; n < 48779; n++) {
        coin = fptype(std::rand())/RAND_MAX;
        if (coin <= 0.5) {
          tails++;
        }
        else {
          heads++;
        }
      }
      CPPUNIT_ASSERT( heads + tails == n );
      fptype pc = fptype(heads) / n;

      // With a maximum error of 0.01 the classical coin is fair at 99.999%
      // level of confidence.
      CPPUNIT_ASSERT( (0.5 - 0.01) <= pc );
      CPPUNIT_ASSERT( (0.5 + 0.01) >= pc );

      // Now test the quantum coin.
      Qubit q0(field(1,0), field(0,0)),
            q1(field(0,0), field(1,0)), x;
      Gate h;
      h.HGate();

      n = 0, heads = 0, tails = 0;
      for ( ; n < 48779; n++) {
        x = q0;
        x = h*x;
        x.measure();

        if (x == q0) {
          tails++;
        }
        else if (x == q1) {
          heads++;
        }
      }
      CPPUNIT_ASSERT( heads + tails == n );
      fptype pq = fptype(heads) / n;

      // With a maximum error of 0.01 the quantum coin is also fair at
      // 99.999% level of confidence.
      CPPUNIT_ASSERT( (0.5 - 0.01) <= pq );
      CPPUNIT_ASSERT( (0.5 + 0.01) >= pq );
    }
};

} // namespace QuCoSi

#endif // QUCOSI_ALGORITHMSTEST_H

// vim: shiftwidth=2 textwidth=78

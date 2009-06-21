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
  CPPUNIT_TEST(testDeutsch);
  CPPUNIT_TEST(testDeutschJozsa);
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

    void testDeutsch()
    {
      Qubit q0(0,1), q1(1,1), x, r0 ,r1;

      Gate h, u;
      h.HGate().tensorPowSet(2);
      r0 = q0.tensorDot(q1);
      r1 = q1.tensorDot(q1);
      std::vector<int> f(2);

      // The function f is constant: f(0) == f(1).
      f[0] = 0;
      f[1] = 0;
      x = h * q0.tensorDot(q1);
      x = h * u.UfGate(f) * x;
      CPPUNIT_ASSERT( x.isApprox(r0) );

      f[0] = 1;
      f[1] = 1;
      x = h * q0.tensorDot(q1);
      x = h * u.UfGate(f) * x;
      CPPUNIT_ASSERT( x.isApprox(-r0) );

      // The function f is balanced: f(0) != f(1).
      f[0] = 0;
      f[1] = 1;
      x = h * q0.tensorDot(q1);
      x = h * u.UfGate(f) * x;
      CPPUNIT_ASSERT( x.isApprox(r1) );

      f[0] = 1;
      f[1] = 0;
      x = h * q0.tensorDot(q1);
      x = h * u.UfGate(f) * x;
      CPPUNIT_ASSERT( x.isApprox(-r1) );
    }

    void testDeutschJozsa()
    {
      Qubit x, y, r;
      Gate h, u, i;
      std::vector<int> f1(2), f2(4);

      y = (Qubit(0,1) - Qubit(1,1)) * std::sqrt(0.5);

      //
      // Test with 2 qubits.
      //
      r = Qubit(0,1).tensorDot(y);

      // The function f is constant: f(0) == f(1).
      f1[0] = 0;
      f1[1] = 0;

      x = Qubit(1,2);
      x = h.HGate().tensorPow(2) * x;
      x = u.UfGate(f1) * x;
      x = h.tensorDot(i.IGate()) * x;
      x.measurePartial(1);
      CPPUNIT_ASSERT( x.isApprox(r) );

      f1[0] = 1;
      f1[1] = 1;

      x = Qubit(1,2);
      x = h.HGate().tensorPow(2) * x;
      x = u.UfGate(f1) * x;
      x = h.tensorDot(i.IGate()) * x;
      x.measurePartial(1);
      CPPUNIT_ASSERT( x.isApprox(-r) );

      // The function f is balanced: f(0) != f(1).
      f1[0] = 1;
      f1[1] = 0;

      x = Qubit(1,2);
      x = h.HGate().tensorPow(2) * x;
      x = u.UfGate(f1) * x;
      x = h.tensorDot(i.IGate()) * x;
      x.measurePartial(1);
      CPPUNIT_ASSERT( !x.isApprox(r) && !x.isApprox(-r) );

      f1[0] = 0;
      f1[1] = 1;

      x = Qubit(1,2);
      x = h.HGate().tensorPow(2) * x;
      x = u.UfGate(f1) * x;
      x = h.tensorDot(i.IGate()) * x;
      x.measurePartial(1);
      CPPUNIT_ASSERT( !x.isApprox(r) && !x.isApprox(-r) );

      //
      // Test with 3 qubits.
      //
      r = Qubit(0,2).tensorDot(y);

      // The function f is constant: f(x) = f(y) for all x, y.
      f2[0] = 0;
      f2[1] = 0;
      f2[2] = 0;
      f2[3] = 0;

      x = Qubit(1,3);
      x = h.HGate().tensorPow(3) * x;
      x = u.UfGate(f2) * x;
      x = h.HGate().tensorPow(2).tensorDot(i.IGate()) * x;
      x.measurePartial(2);
      CPPUNIT_ASSERT( x.isApprox(r) );

      f2[0] = 1;
      f2[1] = 1;
      f2[2] = 1;
      f2[3] = 1;

      x = Qubit(1,3);
      x = h.HGate().tensorPow(3) * x;
      x = u.UfGate(f2) * x;
      x = h.HGate().tensorPow(2).tensorDot(i.IGate()) * x;
      x.measurePartial(2);
      CPPUNIT_ASSERT( x.isApprox(-r) );
      
      // The function f is balanced.
      f2[0] = 0;
      f2[1] = 0;
      f2[2] = 1;
      f2[3] = 1;

      x = Qubit(1,3);
      x = h.HGate().tensorPow(3) * x;
      x = u.UfGate(f2) * x;
      x = h.HGate().tensorPow(2).tensorDot(i.IGate()) * x;
      x.measurePartial(2);
      CPPUNIT_ASSERT( !x.isApprox(r) && !x.isApprox(-r) );

      f2[0] = 1;
      f2[1] = 1;
      f2[2] = 0;
      f2[3] = 0;

      x = Qubit(1,3);
      x = h.HGate().tensorPow(3) * x;
      x = u.UfGate(f2) * x;
      x = h.HGate().tensorPow(2).tensorDot(i.IGate()) * x;
      x.measurePartial(2);
      CPPUNIT_ASSERT( !x.isApprox(r) && !x.isApprox(-r) );

      f2[0] = 1;
      f2[1] = 0;
      f2[2] = 1;
      f2[3] = 0;

      x = Qubit(1,3);
      x = h.HGate().tensorPow(3) * x;
      x = u.UfGate(f2) * x;
      x = h.HGate().tensorPow(2).tensorDot(i.IGate()) * x;
      x.measurePartial(2);
      CPPUNIT_ASSERT( !x.isApprox(r) && !x.isApprox(-r) );

      f2[0] = 0;
      f2[1] = 1;
      f2[2] = 0;
      f2[3] = 1;

      x = Qubit(1,3);
      x = h.HGate().tensorPow(3) * x;
      x = u.UfGate(f2) * x;
      x = h.HGate().tensorPow(2).tensorDot(i.IGate()) * x;
      x.measurePartial(2);
      CPPUNIT_ASSERT( !x.isApprox(r) && !x.isApprox(-r) );

      f2[0] = 1;
      f2[1] = 0;
      f2[2] = 0;
      f2[3] = 1;

      x = Qubit(1,3);
      x = h.HGate().tensorPow(3) * x;
      x = u.UfGate(f2) * x;
      x = h.HGate().tensorPow(2).tensorDot(i.IGate()) * x;
      x.measurePartial(2);
      CPPUNIT_ASSERT( !x.isApprox(r) && !x.isApprox(-r) );

      f2[0] = 0;
      f2[1] = 1;
      f2[2] = 1;
      f2[3] = 0;

      x = Qubit(1,3);
      x = h.HGate().tensorPow(3) * x;
      x = u.UfGate(f2) * x;
      x = h.HGate().tensorPow(2).tensorDot(i.IGate()) * x;
      x.measurePartial(2);
      CPPUNIT_ASSERT( !x.isApprox(r) && !x.isApprox(-r) );
    }
};

} // namespace QuCoSi

#endif // QUCOSI_ALGORITHMSTEST_H

// vim: shiftwidth=2 textwidth=78

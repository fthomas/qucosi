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

#include <QuCoSi/Gate>
#include <QuCoSi/Qubit>

namespace QuCoSi {

class AlgorithmsTest : public CppUnit::TestFixture
{
  CPPUNIT_TEST_SUITE(AlgorithmsTest);
  CPPUNIT_TEST(testCoinFlipping);
  CPPUNIT_TEST(testDeutsch);
  CPPUNIT_TEST(testDeutschJozsa);
  CPPUNIT_TEST(testBernsteinVazirani);
  CPPUNIT_TEST(testSimon);
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
      Qubit q0(0,1), q1(1,1), x;
      Gate h;
      h.H();

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
      h.H().tensorPowSet(2);
      r0 = q0.tensorDot(q1);
      r1 = q1.tensorDot(q1);
      std::vector<int> f(2);

      // The function f is constant: f(0) == f(1).
      f[0] = 0;
      f[1] = 0;
      x = h * q0.tensorDot(q1);
      x = h * u.U(f) * x;
      CPPUNIT_ASSERT( x.isApprox(r0) );

      f[0] = 1;
      f[1] = 1;
      x = h * q0.tensorDot(q1);
      x = h * u.U(f) * x;
      CPPUNIT_ASSERT( x.isApprox(-r0) );

      // The function f is balanced: f(0) != f(1).
      f[0] = 0;
      f[1] = 1;
      x = h * q0.tensorDot(q1);
      x = h * u.U(f) * x;
      CPPUNIT_ASSERT( x.isApprox(r1) );

      f[0] = 1;
      f[1] = 0;
      x = h * q0.tensorDot(q1);
      x = h * u.U(f) * x;
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
      x = h.H().tensorPow(2) * x;
      x = u.U(f1) * x;
      x = h.tensorDot(i.I()) * x;
      x.measurePartial(1);
      CPPUNIT_ASSERT( x.isApprox(r) );

      f1[0] = 1;
      f1[1] = 1;

      x = Qubit(1,2);
      x = h.H().tensorPow(2) * x;
      x = u.U(f1) * x;
      x = h.tensorDot(i.I()) * x;
      x.measurePartial(1);
      CPPUNIT_ASSERT( x.isApprox(-r) );

      // The function f is balanced: f(0) != f(1).
      f1[0] = 1;
      f1[1] = 0;

      x = Qubit(1,2);
      x = h.H().tensorPow(2) * x;
      x = u.U(f1) * x;
      x = h.tensorDot(i.I()) * x;
      x.measurePartial(1);
      CPPUNIT_ASSERT( !x.isApprox(r) && !x.isApprox(-r) );

      f1[0] = 0;
      f1[1] = 1;

      x = Qubit(1,2);
      x = h.H().tensorPow(2) * x;
      x = u.U(f1) * x;
      x = h.tensorDot(i.I()) * x;
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
      x = h.H().tensorPow(3) * x;
      x = u.U(f2) * x;
      x = h.H().tensorPow(2).tensorDot(i.I()) * x;
      x.measurePartial(2);
      CPPUNIT_ASSERT( x.isApprox(r) );

      f2[0] = 1;
      f2[1] = 1;
      f2[2] = 1;
      f2[3] = 1;

      x = Qubit(1,3);
      x = h.H().tensorPow(3) * x;
      x = u.U(f2) * x;
      x = h.H().tensorPow(2).tensorDot(i.I()) * x;
      x.measurePartial(2);
      CPPUNIT_ASSERT( x.isApprox(-r) );
      
      // The function f is balanced.
      f2[0] = 0;
      f2[1] = 0;
      f2[2] = 1;
      f2[3] = 1;

      x = Qubit(1,3);
      x = h.H().tensorPow(3) * x;
      x = u.U(f2) * x;
      x = h.H().tensorPow(2).tensorDot(i.I()) * x;
      x.measurePartial(2);
      CPPUNIT_ASSERT( !x.isApprox(r) && !x.isApprox(-r) );

      f2[0] = 1;
      f2[1] = 1;
      f2[2] = 0;
      f2[3] = 0;

      x = Qubit(1,3);
      x = h.H().tensorPow(3) * x;
      x = u.U(f2) * x;
      x = h.H().tensorPow(2).tensorDot(i.I()) * x;
      x.measurePartial(2);
      CPPUNIT_ASSERT( !x.isApprox(r) && !x.isApprox(-r) );

      f2[0] = 1;
      f2[1] = 0;
      f2[2] = 1;
      f2[3] = 0;

      x = Qubit(1,3);
      x = h.H().tensorPow(3) * x;
      x = u.U(f2) * x;
      x = h.H().tensorPow(2).tensorDot(i.I()) * x;
      x.measurePartial(2);
      CPPUNIT_ASSERT( !x.isApprox(r) && !x.isApprox(-r) );

      f2[0] = 0;
      f2[1] = 1;
      f2[2] = 0;
      f2[3] = 1;

      x = Qubit(1,3);
      x = h.H().tensorPow(3) * x;
      x = u.U(f2) * x;
      x = h.H().tensorPow(2).tensorDot(i.I()) * x;
      x.measurePartial(2);
      CPPUNIT_ASSERT( !x.isApprox(r) && !x.isApprox(-r) );

      f2[0] = 1;
      f2[1] = 0;
      f2[2] = 0;
      f2[3] = 1;

      x = Qubit(1,3);
      x = h.H().tensorPow(3) * x;
      x = u.U(f2) * x;
      x = h.H().tensorPow(2).tensorDot(i.I()) * x;
      x.measurePartial(2);
      CPPUNIT_ASSERT( !x.isApprox(r) && !x.isApprox(-r) );

      f2[0] = 0;
      f2[1] = 1;
      f2[2] = 1;
      f2[3] = 0;

      x = Qubit(1,3);
      x = h.H().tensorPow(3) * x;
      x = u.U(f2) * x;
      x = h.H().tensorPow(2).tensorDot(i.I()) * x;
      x.measurePartial(2);
      CPPUNIT_ASSERT( !x.isApprox(r) && !x.isApprox(-r) );
    }

    void testBernsteinVazirani()
    {
      Gate h, u, x, x0, x1 ,x2, x3, x4;
      h.H().tensorPowSet(6);

      x.X();
      x0.C(5,0,6,x);
      x1.C(5,1,6,x);
      x2.C(5,2,6,x);
      x3.C(5,3,6,x);
      x4.C(5,4,6,x);

      CPPUNIT_ASSERT( (h*(x0*x1*x4)*h*Qubit(1,6)).isApprox(
        Qubit(25,5).tensorDot(Qubit(1,1))) );
      CPPUNIT_ASSERT( (h*(x0*x1*x2)*h*Qubit(1,6)).isApprox(
        Qubit(28,5).tensorDot(Qubit(1,1))) );
      CPPUNIT_ASSERT( (h*(x2)*h*Qubit(1,6)).isApprox(
        Qubit(4,5).tensorDot(Qubit(1,1))) );

      std::vector<int> f(32);
      for (int i = 0; i < 32; i++) {
        f[i] = 0;
      }

      // a = 25 =  1  1  0  0  1
      //          a4 a3 a2 a1 a0
      // f[x] = a0*x0 o+ a1*x1 o+ ... , where o+ denotes the binary addition:
      // 0 o+ 0 = 0, 0 o+ 1 = 1, 1 o+ 0 = 1, 1 o+ 1 = 0

      f[1]  = 1; // 1*0 o+ 1*0 o+ 0*0 o+ 0*0 o+ 1*1 = 1
      f[3]  = 1; // 1*0 o+ 1*0 o+ 0*0 o+ 0*1 o+ 1*1 = 1
      f[5]  = 1; // 1*0 o+ 1*0 o+ 0*1 o+ 0*0 o+ 1*1 = 1
      f[7]  = 1; // 1*0 o+ 1*0 o+ 0*1 o+ 0*1 o+ 1*1 = 1
      f[8]  = 1; // 1*0 o+ 1*1 o+ 0*0 o+ 0*0 o+ 1*0 = 1
      f[10] = 1; // 1*0 o+ 1*1 o+ 0*0 o+ 0*1 o+ 1*0 = 1
      f[12] = 1; // 1*0 o+ 1*1 o+ 0*1 o+ 0*0 o+ 1*0 = 1
      f[14] = 1; // 1*0 o+ 1*1 o+ 0*1 o+ 0*1 o+ 1*0 = 1
      // and so on and so forth
      f[16] = 1; f[18] = 1; f[20] = 1; f[22] = 1;
      f[25] = 1; f[27] = 1; f[29] = 1; f[31] = 1;

      u.U(f);
      CPPUNIT_ASSERT( u == x0*x1*x4 );
      CPPUNIT_ASSERT( (h*u*h*Qubit(1,6)).isApprox(
        Qubit(25,5).tensorDot(Qubit(1,1))) );

      // Automated way to create the appropriate function for U().
      std::vector<int> g(32);
      for (int i = 0; i < 32; i++) {
        g[i] = bwise_bin_dot(i,25);
      }
      CPPUNIT_ASSERT( u == Gate().U(g) );
    }

    void testSimon()
    {
      std::vector<int> f(4);
      // f: {0, 1, 2, 3} -> {0, 1, 2, 3}
      // f[x] = f[x o+ a], with a = 2 ^= 10
      int a = 2;
      f[0] = 1; f[2] = 1;
      f[1] = 2; f[3] = 2;

      Qubit r(0,4), t(0,4);
      Gate h, u;

      r = h.H().tensorPowSet(2).applyToSet(0,4) * r;

      // Check that our qubits are intialized correctly.
      t.setZero();
      for (int x = 0; x <= 3; x++) {
        t += 0.5 * Qubit(x,2).tensorDot(Qubit(0,2));
      }
      CPPUNIT_ASSERT( r.isApprox(t) );

      r = u.U(f,2) * r;

      // Check the call of the oracle U_f.
      t.setZero();
      for (int x = 0; x <= 3; x++) {
        t += 0.5 * Qubit(x,2).tensorDot(Qubit(f[x],2));
      }
      CPPUNIT_ASSERT( r.isApprox(t) );

      // Measure the output qubits.
      std::vector<int> p(4);
      p[0] = 2; p[1] = 3; p[2] = 0; p[3] = 1;
      Gate s; s.S(p);
      r = s * r; r.measurePartial(2); r = s.transpose() * r;

      // Check the qubit state after measuring the output qubits.
      Qubit t1[4];
      int x0 = 0;
      for (int x = 0; x <= 3; x++) {
        t1[x] = std::sqrt(0.5)*(Qubit(x,2) + Qubit(x^a,2));
        t1[x].tensorDotSet(Qubit(f[x],2));
        if (r.isApprox(t1[x])) {
          x0 = x^a;
        }
      }
      CPPUNIT_ASSERT( r.isApprox(t1[0]) || r.isApprox(t1[1]) ||
                      r.isApprox(t1[2]) || r.isApprox(t1[3]) );

      r = h * r;
      r = r.first(2);

      // Check the input qubit state after applying the Hadamard gate.
      t.resize(4); t.setZero();
      for (int y = 0; y <= 1; y++) {
        if (bwise_bin_dot(a,y) == 0) {
          t += 1./std::pow(2, 1./2.) * std::pow(-1,bwise_bin_dot(x0,y)) *
               Qubit(y,2);
        }
      }
      CPPUNIT_ASSERT( r.isApprox(t) );

      r.measure();

      // Our input a = 2. Our result y = 1.
      // => a_1*0 o+ a_0*1 = 0. Because a != 0 it follows that a = 2.
      CPPUNIT_ASSERT( r.isApprox(Qubit(0,2)) || r.isApprox(Qubit(1,2)) ||
                      r.isApprox(-Qubit(1,2)) );
    }
};

} // namespace QuCoSi

#endif // QUCOSI_ALGORITHMSTEST_H

// vim: shiftwidth=2 textwidth=78

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
  CPPUNIT_TEST(testSGate);
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

    void testSGate()
    {
      Gate s,t;
      Qubit a(field(2,0), field(3,0)),
            b(field(4,0), field(5,0)),
            c(field(6,0), field(7,0)),
            d(field(8,0), field(9,0));

      // Permute two qubits.
      Qubit ab, ba;
      ab = a.tensorDot(b);
      ba = b.tensorDot(a);

      CPPUNIT_ASSERT( s.SpqGate(0,1,2)*ab == ba );
      CPPUNIT_ASSERT( s.SpqGate(0,1,2) == t.SWAPGate() );

      // Permute three qubits.
      Qubit abc, acb, bac, bca, cab, cba;
      abc = a.tensorDot(b.tensorDot(c));
      acb = a.tensorDot(c.tensorDot(b));
      bac = b.tensorDot(a.tensorDot(c));
      bca = b.tensorDot(c.tensorDot(a));
      cab = c.tensorDot(a.tensorDot(b));
      cba = c.tensorDot(b.tensorDot(a));

      CPPUNIT_ASSERT( s.SpqGate(0,0,3)*abc == abc );
      CPPUNIT_ASSERT( s.SpqGate(0,1,3)*abc == bac );
      CPPUNIT_ASSERT( s.SpqGate(1,2,3)*abc == acb );
      CPPUNIT_ASSERT( s.SpqGate(0,2,3)*abc == cba );
      CPPUNIT_ASSERT( s.SpqGate(1,2,3)*t.SpqGate(0,1,3)*abc == bca );
      CPPUNIT_ASSERT( s.SpqGate(0,1,3)*t.SpqGate(1,2,3)*abc == cab );

      std::vector<int> p1(3), p2(3);
      p1[0] = 1; p2[0] = 2;
      p1[1] = 2; p2[1] = 0;
      p1[2] = 0; p2[2] = 1;

      CPPUNIT_ASSERT( s.SGate(p1)*abc == bca );
      CPPUNIT_ASSERT( s.SGate(p2)*abc == cab );

      // Permute four qubits.
      Qubit abcd, abdc, acbd, acdb, adbc, adcb,
            bacd, badc, bcad, bcda, bdac, bdca,
            cabd, cadb, cbad, cbda, cdab, cdba,
            dacb, dabc, dbac, dbca, dcab, dcba;
      abcd = a.tensorDot(b.tensorDot(c.tensorDot(d)));
      abdc = a.tensorDot(b.tensorDot(d.tensorDot(c)));
      acbd = a.tensorDot(c.tensorDot(b.tensorDot(d)));
      acdb = a.tensorDot(c.tensorDot(d.tensorDot(b)));
      adbc = a.tensorDot(d.tensorDot(b.tensorDot(c)));
      adcb = a.tensorDot(d.tensorDot(c.tensorDot(b)));

      bacd = b.tensorDot(a.tensorDot(c.tensorDot(d)));
      badc = b.tensorDot(a.tensorDot(d.tensorDot(c)));
      bcad = b.tensorDot(c.tensorDot(a.tensorDot(d)));
      bcda = b.tensorDot(c.tensorDot(d.tensorDot(a)));
      bdac = b.tensorDot(d.tensorDot(a.tensorDot(c)));
      bdca = b.tensorDot(d.tensorDot(c.tensorDot(a)));

      cabd = c.tensorDot(a.tensorDot(b.tensorDot(d)));
      cadb = c.tensorDot(a.tensorDot(d.tensorDot(b)));
      cbad = c.tensorDot(b.tensorDot(a.tensorDot(d)));
      cbda = c.tensorDot(b.tensorDot(d.tensorDot(a)));
      cdab = c.tensorDot(d.tensorDot(a.tensorDot(b)));
      cdba = c.tensorDot(d.tensorDot(b.tensorDot(a)));

      dabc = d.tensorDot(a.tensorDot(b.tensorDot(c)));
      dacb = d.tensorDot(a.tensorDot(c.tensorDot(b)));
      dbac = d.tensorDot(b.tensorDot(a.tensorDot(c)));
      dbca = d.tensorDot(b.tensorDot(c.tensorDot(a)));
      dcab = d.tensorDot(c.tensorDot(a.tensorDot(b)));
      dcba = d.tensorDot(c.tensorDot(b.tensorDot(a)));

      CPPUNIT_ASSERT( s.SpqGate(0,1,4)*abcd == bacd );
      CPPUNIT_ASSERT( s.SpqGate(0,2,4)*abcd == cbad );
      CPPUNIT_ASSERT( s.SpqGate(0,3,4)*abcd == dbca );
      CPPUNIT_ASSERT( s.SpqGate(1,2,4)*abcd == acbd );
      CPPUNIT_ASSERT( s.SpqGate(1,3,4)*abcd == adcb );
      CPPUNIT_ASSERT( s.SpqGate(2,3,4)*abcd == abdc );

      std::vector<int> p3(4), p4(4), p5(4), p6(4);
      p3[0] = 3; p4[0] = 1; p5[0] = 3; p6[0] = 3;
      p3[1] = 2; p4[1] = 2; p5[1] = 0; p6[1] = 2;
      p3[2] = 1; p4[2] = 3; p5[2] = 1; p6[2] = 0;
      p3[3] = 0; p4[3] = 0; p5[3] = 2; p6[3] = 1;

      CPPUNIT_ASSERT( s.SGate(p3)*abcd == dcba );
      CPPUNIT_ASSERT( s.SGate(p3).transpose()*dcba == abcd );
      CPPUNIT_ASSERT( s.SGate(p4)*abcd == bcda );
      CPPUNIT_ASSERT( s.SGate(p4).transpose()*bcda == abcd );
      CPPUNIT_ASSERT( s.SGate(p5)*abcd == dabc );
      CPPUNIT_ASSERT( s.SGate(p5).transpose()*dabc == abcd );
      CPPUNIT_ASSERT( s.SGate(p6)*abcd == dcab );
      CPPUNIT_ASSERT( s.SGate(p6).transpose()*dcab == abcd );
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

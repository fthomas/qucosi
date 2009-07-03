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

#ifndef QUCOSI_GATETEST_H
#define QUCOSI_GATETEST_H

#include <cppunit/TestFixture.h>
#include <cppunit/extensions/HelperMacros.h>

#include <QuCoSi/Gate>

namespace QuCoSi {

class GateTest : public CppUnit::TestFixture
{
  CPPUNIT_TEST_SUITE(GateTest);
  CPPUNIT_TEST(testTensorPow);
  CPPUNIT_TEST(testApplyToPos);
  CPPUNIT_TEST(testCGate);
  CPPUNIT_TEST(testSGate);
  CPPUNIT_TEST(testUfGate);
  CPPUNIT_TEST(testFGate);
  CPPUNIT_TEST(testCircuitIdentities);
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

      h0 << 1,  0,  1,  0,
            0,  1,  0,  1,
            1,  0, -1,  0,
            0,  1,  0, -1;
      h0 *= std::sqrt(1./2);

      h1 << 1,  1,  0,  0,
            1, -1,  0,  0,
            0,  0,  1,  1,
            0,  0,  1, -1;
      h1 *= std::sqrt(1./2);

      CPPUNIT_ASSERT( h0.isApprox(h.applyToPos(0,2)) );
      CPPUNIT_ASSERT( h1.isApprox(h.applyToPos(1,2)) );
      CPPUNIT_ASSERT( h.applyToPos(0,2).isUnitary() );
      CPPUNIT_ASSERT( h.applyToPos(1,2).isUnitary() );

      h0 = h0.HGate().tensorPowSet(2).applyToPos(0,3);
      h1 = h1.HGate().tensorPowSet(2).tensorDotSet(Gate().IGate());
      CPPUNIT_ASSERT( h0.isApprox(h1) );
    }

    void testCGate()
    {
      Gate cnot, ccnot, cswap, c, u, g;
      cnot.CNOTGate();
      ccnot.CCNOTGate();
      cswap.CSWAPGate();

      CPPUNIT_ASSERT( c.CGate(1,0,2,u.XGate()) == cnot );
      CPPUNIT_ASSERT( c.CGate(1,0,3,u.CGate(1,0,2,g.XGate())) == ccnot );
      CPPUNIT_ASSERT( c.CGate(1,0,3,u.CNOTGate()) == ccnot );
      CPPUNIT_ASSERT( c.CGate(1,0,3,u.SWAPGate()) == cswap );
      CPPUNIT_ASSERT( c.CGate(2,1,3,u.XGate()) == cnot.applyToPos(1,3) );

      g.resize(4,4);
      g << 1, 0, 0, 0,
           0, 0, 0, 1,
           0, 0, 1, 0,
           0, 1, 0, 0;
      CPPUNIT_ASSERT( c.CGate(0,1,2,u.XGate()) == g );
      u.SWAPGate();
      CPPUNIT_ASSERT( u*cnot*u == g );
    }

    void testSGate()
    {
      Gate s,t;
      Qubit a(field(2,0), field(3,0)),
            b(field(4,0), field(5,0)),
            c(field(6,0), field(7,0)),
            d(field(8,0), field(9,0)),
            e(field(10,0), field(11,0));

      // Permute two qubits.
      Qubit ab, ba;
      ab = a.tensorDot(b);
      ba = b.tensorDot(a);

      CPPUNIT_ASSERT( s.SGate(0,1,2)*ab == ba );
      CPPUNIT_ASSERT( s.SGate(0,1,2) == t.SWAPGate() );

      // Permute three qubits.
      Qubit abc, acb, bac, bca, cab, cba;
      abc = a.tensorDot(b.tensorDot(c));
      acb = a.tensorDot(c.tensorDot(b));
      bac = b.tensorDot(a.tensorDot(c));
      bca = b.tensorDot(c.tensorDot(a));
      cab = c.tensorDot(a.tensorDot(b));
      cba = c.tensorDot(b.tensorDot(a));

      CPPUNIT_ASSERT( s.SGate(0,0,3)*abc == abc );
      CPPUNIT_ASSERT( s.SGate(0,1,3)*abc == bac );
      CPPUNIT_ASSERT( s.SGate(1,2,3)*abc == acb );
      CPPUNIT_ASSERT( s.SGate(0,2,3)*abc == cba );
      CPPUNIT_ASSERT( s.SGate(1,2,3)*t.SGate(0,1,3)*abc == bca );
      CPPUNIT_ASSERT( s.SGate(0,1,3)*t.SGate(1,2,3)*abc == cab );

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

      CPPUNIT_ASSERT( s.SGate(0,1,4)*abcd == bacd );
      CPPUNIT_ASSERT( s.SGate(0,2,4)*abcd == cbad );
      CPPUNIT_ASSERT( s.SGate(0,3,4)*abcd == dbca );
      CPPUNIT_ASSERT( s.SGate(1,2,4)*abcd == acbd );
      CPPUNIT_ASSERT( s.SGate(1,3,4)*abcd == adcb );
      CPPUNIT_ASSERT( s.SGate(2,3,4)*abcd == abdc );

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

      // Permute five qubits.
      Qubit abcde, abced, abdce, abedc, acbde, acbed, adcbe, aecdb, bacde,
            cbade, dbcae, ebcda, edcab;
      abcde = a.tensorDot(b).tensorDot(c).tensorDot(d).tensorDot(e);
      abced = a.tensorDot(b).tensorDot(c).tensorDot(e).tensorDot(d);
      abdce = a.tensorDot(b).tensorDot(d).tensorDot(c).tensorDot(e);
      abedc = a.tensorDot(b).tensorDot(e).tensorDot(d).tensorDot(c);
      acbde = a.tensorDot(c).tensorDot(b).tensorDot(d).tensorDot(e);
      acbed = a.tensorDot(c).tensorDot(b).tensorDot(e).tensorDot(d);
      adcbe = a.tensorDot(d).tensorDot(c).tensorDot(b).tensorDot(e);
      aecdb = a.tensorDot(e).tensorDot(c).tensorDot(d).tensorDot(b);
      bacde = b.tensorDot(a).tensorDot(c).tensorDot(d).tensorDot(e);
      cbade = c.tensorDot(b).tensorDot(a).tensorDot(d).tensorDot(e);
      dbcae = d.tensorDot(b).tensorDot(c).tensorDot(a).tensorDot(e);
      ebcda = e.tensorDot(b).tensorDot(c).tensorDot(d).tensorDot(a);
      edcab = e.tensorDot(d).tensorDot(c).tensorDot(a).tensorDot(b);

      std::vector<int> r1(5), r2(5);
      r1[0] = 0; r2[0] = 4;
      r1[1] = 2; r2[1] = 3;
      r1[2] = 1; r2[2] = 2;
      r1[3] = 4; r2[3] = 0;
      r1[4] = 3; r2[4] = 1;

      CPPUNIT_ASSERT( s.SGate(r1)*abcde == acbed );
      CPPUNIT_ASSERT( s.SGate(r2)*abcde == edcab );

      CPPUNIT_ASSERT( s.SGate(0,1,5)*abcde == bacde );
      CPPUNIT_ASSERT( s.SGate(0,2,5)*abcde == cbade );
      CPPUNIT_ASSERT( s.SGate(0,3,5)*abcde == dbcae );
      CPPUNIT_ASSERT( s.SGate(0,4,5)*abcde == ebcda );

      CPPUNIT_ASSERT( s.SGate(1,2,5)*abcde == acbde );
      CPPUNIT_ASSERT( s.SGate(1,3,5)*abcde == adcbe );
      CPPUNIT_ASSERT( s.SGate(1,4,5)*abcde == aecdb );

      CPPUNIT_ASSERT( s.SGate(2,3,5)*abcde == abdce );
      CPPUNIT_ASSERT( s.SGate(2,4,5)*abcde == abedc );

      CPPUNIT_ASSERT( s.SGate(3,4,5)*abcde == abced );
    }

    void testUfGate()
    {
      Gate uf, uc, uh, a1, a2, a3, a4, a5;

      //
      // Test all 4x4 UfGates.
      //
      std::vector<int> f(2);

      f[0] = 0; f[1] = 0;
      uh.resize(4,4);
      uh << 1, 0, 0, 0,
            0, 1, 0, 0,
            0, 0, 1, 0,
            0, 0, 0, 1;
      uc.IGate().tensorPowSet(2);
      uf.UfGate(f);
      CPPUNIT_ASSERT( uh == uc );
      CPPUNIT_ASSERT( uh == uf );

      f[0] = 0; f[1] = 1;
      uh.resize(4,4);
      uh << 1, 0, 0, 0,
            0, 1, 0, 0,
            0, 0, 0, 1,
            0, 0, 1, 0;
      uc.CGate(1,0,2,a1.XGate());
      uf.UfGate(f);
      CPPUNIT_ASSERT( uh == uc );
      CPPUNIT_ASSERT( uh == uf );

      f[0] = 1; f[1] = 0;
      uh.resize(4,4);
      uh << 0, 1, 0, 0,
            1, 0, 0, 0,
            0, 0, 1, 0,
            0, 0, 0, 1;
      uc = a1.XGate().applyToPos(1,2) * a2.CGate(1,0,2,a3.XGate());
      uf.UfGate(f);
      CPPUNIT_ASSERT( uh == uc );
      CPPUNIT_ASSERT( uh == uf );

      f[0] = 1; f[1] = 1;
      uh.resize(4,4);
      uh << 0, 1, 0, 0,
            1, 0, 0, 0,
            0, 0, 0, 1,
            0, 0, 1, 0;
      uc = uc.XGate().applyToPos(1,2);
      uf.UfGate(f);
      CPPUNIT_ASSERT( uh == uc );
      CPPUNIT_ASSERT( uh == uf );

      //
      // Test all 8x8 UfGates.
      //
      f.resize(4);

      // binary 0
      f[0] = 0; f[1] = 0; f[2] = 0; f[3] = 0;
      uh.resize(8,8);
      uh.setIdentity();
      uc.IGate().tensorPowSet(3);
      uf.UfGate(f);
      CPPUNIT_ASSERT( uh == uc );
      CPPUNIT_ASSERT( uh == uf );

      // binary 15
      f[0] = 1; f[1] = 1; f[2] = 1; f[3] = 1;
      uh.resize(8,8);
      uh << 0, 1, 0, 0, 0, 0, 0, 0,
            1, 0, 0, 0, 0, 0, 0, 0,
            0, 0, 0, 1, 0, 0, 0, 0,
            0, 0, 1, 0, 0, 0, 0, 0,
            0, 0, 0, 0, 0, 1, 0, 0,
            0, 0, 0, 0, 1, 0, 0, 0,
            0, 0, 0, 0, 0, 0, 0, 1,
            0, 0, 0, 0, 0, 0, 1, 0;
      uc = uc.XGate().applyToPos(2,3);
      uf.UfGate(f);
      CPPUNIT_ASSERT( uh == uc );
      CPPUNIT_ASSERT( uh == uf );

      // binary 1
      f[0] = 0; f[1] = 0; f[2] = 0; f[3] = 1;
      uh.resize(8,8);
      uh << 1, 0, 0, 0, 0, 0, 0, 0,
            0, 1, 0, 0, 0, 0, 0, 0,
            0, 0, 1, 0, 0, 0, 0, 0,
            0, 0, 0, 1, 0, 0, 0, 0,
            0, 0, 0, 0, 1, 0, 0, 0,
            0, 0, 0, 0, 0, 1, 0, 0,
            0, 0, 0, 0, 0, 0, 0, 1,
            0, 0, 0, 0, 0, 0, 1, 0;
      uc = uc.CGate(1,0,3,a1.CNOTGate());
      uf.UfGate(f);
      CPPUNIT_ASSERT( uh == uc );
      CPPUNIT_ASSERT( uh == uf );

      // binary 14
      f[0] = 1; f[1] = 1; f[2] = 1; f[3] = 0;
      uh.resize(8,8);
      uh << 0, 1, 0, 0, 0, 0, 0, 0,
            1, 0, 0, 0, 0, 0, 0, 0,
            0, 0, 0, 1, 0, 0, 0, 0,
            0, 0, 1, 0, 0, 0, 0, 0,
            0, 0, 0, 0, 0, 1, 0, 0,
            0, 0, 0, 0, 1, 0, 0, 0,
            0, 0, 0, 0, 0, 0, 1, 0,
            0, 0, 0, 0, 0, 0, 0, 1;
      uc = a1.XGate().applyToPos(2,3) * uc.CGate(1,0,3,a2.CNOTGate());
      uf.UfGate(f);
      CPPUNIT_ASSERT( uh == uc );
      CPPUNIT_ASSERT( uh == uf );

      // binary 2
      f[0] = 0; f[1] = 0; f[2] = 1; f[3] = 0;
      uh.resize(8,8);
      uh << 1, 0, 0, 0, 0, 0, 0, 0,
            0, 1, 0, 0, 0, 0, 0, 0,
            0, 0, 1, 0, 0, 0, 0, 0,
            0, 0, 0, 1, 0, 0, 0, 0,
            0, 0, 0, 0, 0, 1, 0, 0,
            0, 0, 0, 0, 1, 0, 0, 0,
            0, 0, 0, 0, 0, 0, 1, 0,
            0, 0, 0, 0, 0, 0, 0, 1;
      uc = a1.XGate().applyToPos(1,3) * uc.CGate(1,0,3,a2.CNOTGate()) *
           a3.XGate().applyToPos(1,3);
      uf.UfGate(f);
      CPPUNIT_ASSERT( uh == uc );
      CPPUNIT_ASSERT( uh == uf );

      // binary 13
      f[0] = 1; f[1] = 1; f[2] = 0; f[3] = 1;
      uh.resize(8,8);
      uh << 0, 1, 0, 0, 0, 0, 0, 0,
            1, 0, 0, 0, 0, 0, 0, 0,
            0, 0, 0, 1, 0, 0, 0, 0,
            0, 0, 1, 0, 0, 0, 0, 0,
            0, 0, 0, 0, 1, 0, 0, 0,
            0, 0, 0, 0, 0, 1, 0, 0,
            0, 0, 0, 0, 0, 0, 0, 1,
            0, 0, 0, 0, 0, 0, 1, 0;
      uc = a1.XGate().applyToPos(1,3) * a2.XGate().applyToPos(2,3) *
           uc.CGate(1,0,3,a3.CNOTGate()) * a4.XGate().applyToPos(1,3);
      uf.UfGate(f);
      CPPUNIT_ASSERT( uh == uc );
      CPPUNIT_ASSERT( uh == uf );

      // binary 3
      f[0] = 0; f[1] = 0; f[2] = 1; f[3] = 1;
      uh.resize(8,8);
      uh << 1, 0, 0, 0, 0, 0, 0, 0,
            0, 1, 0, 0, 0, 0, 0, 0,
            0, 0, 1, 0, 0, 0, 0, 0,
            0, 0, 0, 1, 0, 0, 0, 0,
            0, 0, 0, 0, 0, 1, 0, 0,
            0, 0, 0, 0, 1, 0, 0, 0,
            0, 0, 0, 0, 0, 0, 0, 1,
            0, 0, 0, 0, 0, 0, 1, 0;
      uc = uc.CGate(2,0,3,a1.XGate());
      uf.UfGate(f);
      CPPUNIT_ASSERT( uh == uc );
      CPPUNIT_ASSERT( uh == uf );

      // binary 12
      f[0] = 1; f[1] = 1; f[2] = 0; f[3] = 0;
      uh.resize(8,8);
      uh << 0, 1, 0, 0, 0, 0, 0, 0,
            1, 0, 0, 0, 0, 0, 0, 0,
            0, 0, 0, 1, 0, 0, 0, 0,
            0, 0, 1, 0, 0, 0, 0, 0,
            0, 0, 0, 0, 1, 0, 0, 0,
            0, 0, 0, 0, 0, 1, 0, 0,
            0, 0, 0, 0, 0, 0, 1, 0,
            0, 0, 0, 0, 0, 0, 0, 1;
      uc = a1.XGate().applyToPos(0,3) * uc.CGate(2,0,3,a2.XGate()) *
           a3.XGate().applyToPos(0,3);
      uf.UfGate(f);
      CPPUNIT_ASSERT( uh == uc );
      CPPUNIT_ASSERT( uh == uf );

      // binary 4
      f[0] = 0; f[1] = 1; f[2] = 0; f[3] = 0;
      uh.resize(8,8);
      uh << 1, 0, 0, 0, 0, 0, 0, 0,
            0, 1, 0, 0, 0, 0, 0, 0,
            0, 0, 0, 1, 0, 0, 0, 0,
            0, 0, 1, 0, 0, 0, 0, 0,
            0, 0, 0, 0, 1, 0, 0, 0,
            0, 0, 0, 0, 0, 1, 0, 0,
            0, 0, 0, 0, 0, 0, 1, 0,
            0, 0, 0, 0, 0, 0, 0, 1;
      uc = a1.XGate().applyToPos(0,3) * uc.CCNOTGate() *
           a2.XGate().applyToPos(0,3);
      uf.UfGate(f);
      CPPUNIT_ASSERT( uh == uc );
      CPPUNIT_ASSERT( uh == uf );

      // binary 11
      f[0] = 1; f[1] = 0; f[2] = 1; f[3] = 1;
      uh.resize(8,8);
      uh << 0, 1, 0, 0, 0, 0, 0, 0,
            1, 0, 0, 0, 0, 0, 0, 0,
            0, 0, 1, 0, 0, 0, 0, 0,
            0, 0, 0, 1, 0, 0, 0, 0,
            0, 0, 0, 0, 0, 1, 0, 0,
            0, 0, 0, 0, 1, 0, 0, 0,
            0, 0, 0, 0, 0, 0, 0, 1,
            0, 0, 0, 0, 0, 0, 1, 0;
      uc = a1.XGate().applyToPos(0,3) * a2.XGate().applyToPos(2,3) *
           uc.CCNOTGate() * a3.XGate().applyToPos(0,3);
      uf.UfGate(f);
      CPPUNIT_ASSERT( uh == uc );
      CPPUNIT_ASSERT( uh == uf );

      // binary 5
      f[0] = 0; f[1] = 1; f[2] = 0; f[3] = 1;
      uh.resize(8,8);
      uh << 1, 0, 0, 0, 0, 0, 0, 0,
            0, 1, 0, 0, 0, 0, 0, 0,
            0, 0, 0, 1, 0, 0, 0, 0,
            0, 0, 1, 0, 0, 0, 0, 0,
            0, 0, 0, 0, 1, 0, 0, 0,
            0, 0, 0, 0, 0, 1, 0, 0,
            0, 0, 0, 0, 0, 0, 0, 1,
            0, 0, 0, 0, 0, 0, 1, 0;
      uc = uc.CNOTGate().applyToPos(1,3);
      uf.UfGate(f);
      CPPUNIT_ASSERT( uh == uc );
      CPPUNIT_ASSERT( uh == uf );

      // binary 10
      f[0] = 1; f[1] = 0; f[2] = 1; f[3] = 0;
      uh.resize(8,8);
      uh << 0, 1, 0, 0, 0, 0, 0, 0,
            1, 0, 0, 0, 0, 0, 0, 0,
            0, 0, 1, 0, 0, 0, 0, 0,
            0, 0, 0, 1, 0, 0, 0, 0,
            0, 0, 0, 0, 0, 1, 0, 0,
            0, 0, 0, 0, 1, 0, 0, 0,
            0, 0, 0, 0, 0, 0, 1, 0,
            0, 0, 0, 0, 0, 0, 0, 1;
      uc = a1.XGate().applyToPos(1,3) * uc.CNOTGate().applyToPos(1,3) *
           a2.XGate().applyToPos(1,3);
      uf.UfGate(f);
      CPPUNIT_ASSERT( uh == uc );
      CPPUNIT_ASSERT( uh == uf );

      // binary 6
      f[0] = 0; f[1] = 1; f[2] = 1; f[3] = 0;
      uh.resize(8,8);
      uh << 1, 0, 0, 0, 0, 0, 0, 0,
            0, 1, 0, 0, 0, 0, 0, 0,
            0, 0, 0, 1, 0, 0, 0, 0,
            0, 0, 1, 0, 0, 0, 0, 0,
            0, 0, 0, 0, 0, 1, 0, 0,
            0, 0, 0, 0, 1, 0, 0, 0,
            0, 0, 0, 0, 0, 0, 1, 0,
            0, 0, 0, 0, 0, 0, 0, 1;
      uc = a1.CGate(2,1,3,a3.XGate()) * a2.CGate(2,0,3,a4.XGate());
      uf.UfGate(f);
      CPPUNIT_ASSERT( uh == uc );
      CPPUNIT_ASSERT( uh == uf );

      // binary 9
      f[0] = 1; f[1] = 0; f[2] = 0; f[3] = 1;
      uh.resize(8,8);
      uh << 0, 1, 0, 0, 0, 0, 0, 0,
            1, 0, 0, 0, 0, 0, 0, 0,
            0, 0, 1, 0, 0, 0, 0, 0,
            0, 0, 0, 1, 0, 0, 0, 0,
            0, 0, 0, 0, 1, 0, 0, 0,
            0, 0, 0, 0, 0, 1, 0, 0,
            0, 0, 0, 0, 0, 0, 0, 1,
            0, 0, 0, 0, 0, 0, 1, 0;
      uc = uc.XGate().applyToPos(2,3) * a1.CGate(2,1,3,a2.XGate()) *
           a3.CGate(2,0,3,a4.XGate());
      uf.UfGate(f);
      CPPUNIT_ASSERT( uh == uc );
      CPPUNIT_ASSERT( uh == uf );

      // binary 7
      f[0] = 0; f[1] = 1; f[2] = 1; f[3] = 1;
      uh.resize(8,8);
      uh << 1, 0, 0, 0, 0, 0, 0, 0,
            0, 1, 0, 0, 0, 0, 0, 0,
            0, 0, 0, 1, 0, 0, 0, 0,
            0, 0, 1, 0, 0, 0, 0, 0,
            0, 0, 0, 0, 0, 1, 0, 0,
            0, 0, 0, 0, 1, 0, 0, 0,
            0, 0, 0, 0, 0, 0, 0, 1,
            0, 0, 0, 0, 0, 0, 1, 0;
      uc = a1.XGate().applyToPos(0,3) * a2.XGate().applyToPos(1,3) *
           a3.XGate().applyToPos(2,3) * uc.CCNOTGate() *
           a4.XGate().applyToPos(0,3) * a5.XGate().applyToPos(1,3);
      uf.UfGate(f);
      CPPUNIT_ASSERT( uh == uc );
      CPPUNIT_ASSERT( uh == uf );

      // binary 8
      f[0] = 1; f[1] = 0; f[2] = 0; f[3] = 0;
      uh.resize(8,8);
      uh << 0, 1, 0, 0, 0, 0, 0, 0,
            1, 0, 0, 0, 0, 0, 0, 0,
            0, 0, 1, 0, 0, 0, 0, 0,
            0, 0, 0, 1, 0, 0, 0, 0,
            0, 0, 0, 0, 1, 0, 0, 0,
            0, 0, 0, 0, 0, 1, 0, 0,
            0, 0, 0, 0, 0, 0, 1, 0,
            0, 0, 0, 0, 0, 0, 0, 1;
      uc = a1.XGate().applyToPos(0,3) * a2.XGate().applyToPos(1,3) *
           uc.CCNOTGate() * a3.XGate().applyToPos(0,3) *
           a4.XGate().applyToPos(1,3);
      uf.UfGate(f);
      CPPUNIT_ASSERT( uh == uc );
      CPPUNIT_ASSERT( uh == uf );
    }

    void testFGate()
    {
      field i = field(0,1);
      Qubit q0(0,4),   q1(1,4),   q2(2,4),   q3(3,4),
            q4(4,4),   q5(5,4),   q6(6,4),   q7(7,4),
            q8(8,4),   q9(9,4),   q10(10,4), q11(11,4),
            q12(12,4), q13(13,4), q14(14,4), q15(15,4);

      Gate f;
      f.FGate(4);

      CPPUNIT_ASSERT( (f*(0.5*(q0 + q4 + q8 + q12))).isApprox(
                          0.5*(q0 + q4 + q8 + q12)) );

      CPPUNIT_ASSERT( (f*(0.5*(q1 + q5 + q9 + q13))).isApprox(
                          0.5*(q0 + i*q4 - q8 - i*q12)) );

      CPPUNIT_ASSERT( (f*(0.5*(q2 + q6 + q10 + q14))).isApprox(
                          0.5*(q0 - q4 + q8  - q12)) );

      CPPUNIT_ASSERT( (f*(0.5*(q3 + q7 + q11 + q15))).isApprox(
                          0.5*(q0 - i*q4 - q8 + i*q12)) );
    }

    void testCircuitIdentities()
    {
      Gate g1, g2, g3, g4, g5;

      CPPUNIT_ASSERT( (g1.HGate() * g2.HGate()).isApprox(g3.IGate()) );

      CPPUNIT_ASSERT( (g1.HGate() * g2.XGate() * g3.HGate()).isApprox(
        g4.ZGate()) );

      CPPUNIT_ASSERT( (g1.HGate() * g2.ZGate() * g3.HGate()).isApprox(
        g4.XGate()) );

      CPPUNIT_ASSERT( (g1.HGate().applyToPos(1,2) *
        g2.CGate(1,0,2,g3.ZGate()) * g4.HGate().applyToPos(1,2)).isApprox(
        g5.CNOTGate()) );

      CPPUNIT_ASSERT( g1.CGate(1,0,2,g2.ZGate()).isApprox(
        g3.CGate(0,1,2,g4.ZGate())) );

      CPPUNIT_ASSERT( (g1.HGate().tensorPow(2) * g2.CNOTGate() *
        g3.HGate().tensorPow(2)).isApprox(g4.CGate(0,1,2,g5.XGate())) );
    }
};

} // namespace QuCoSi

#endif // QUCOSI_GATETEST_H

// vim: shiftwidth=2 textwidth=78

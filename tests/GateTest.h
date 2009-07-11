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
  CPPUNIT_TEST(testC);
  CPPUNIT_TEST(testS);
  CPPUNIT_TEST(testU);
  CPPUNIT_TEST(testF);
  CPPUNIT_TEST(testCircuitIdentities);
  CPPUNIT_TEST_SUITE_END();

  public:
    void setUp() {}

    void tearDown() {}

    void testTensorPow()
    {
      Gate h, h_byhand;

      h.H().tensorPowSet(1);
      h_byhand.resize(2,2);
      h_byhand << 1,  1,
                  1, -1;
      h_byhand *= std::sqrt(1./2);
      CPPUNIT_ASSERT( h.isApprox(h_byhand) );
      CPPUNIT_ASSERT( h.isUnitary() );
      CPPUNIT_ASSERT( h_byhand.isUnitary() );

      h.H().tensorPowSet(2);
      h_byhand.resize(4,4);
      h_byhand << 1,  1,  1,  1,
                  1, -1,  1, -1,
                  1,  1, -1, -1,
                  1, -1, -1,  1;
      h_byhand *= std::sqrt(1./4);
      CPPUNIT_ASSERT( h.isApprox(h_byhand) );
      CPPUNIT_ASSERT( h.isUnitary() );
      CPPUNIT_ASSERT( h_byhand.isUnitary() );

      h.H().tensorPowSet(3);
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
      h.H();
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

      CPPUNIT_ASSERT( h0.isApprox(h.applyTo(0,2)) );
      CPPUNIT_ASSERT( h1.isApprox(h.applyTo(1,2)) );
      CPPUNIT_ASSERT( h.applyTo(0,2).isUnitary() );
      CPPUNIT_ASSERT( h.applyTo(1,2).isUnitary() );

      h0 = h0.H().tensorPowSet(2).applyTo(0,3);
      h1 = h1.H().tensorPowSet(2).tensorDotSet(Gate().I());
      CPPUNIT_ASSERT( h0.isApprox(h1) );
    }

    void testC()
    {
      Gate cnot, ccnot, cswap, c, u, g;
      cnot.CNOT();
      ccnot.CCNOT();
      cswap.CSWAP();

      CPPUNIT_ASSERT( c.C(1,0,2,u.X()) == cnot );
      CPPUNIT_ASSERT( c.C(1,0,3,u.C(1,0,2,g.X())) == ccnot );
      CPPUNIT_ASSERT( c.C(1,0,3,u.CNOT()) == ccnot );
      CPPUNIT_ASSERT( c.C(1,0,3,u.SWAP()) == cswap );
      CPPUNIT_ASSERT( c.C(2,1,3,u.X()) == cnot.applyTo(1,3) );

      g.resize(4,4);
      g << 1, 0, 0, 0,
           0, 0, 0, 1,
           0, 0, 1, 0,
           0, 1, 0, 0;
      CPPUNIT_ASSERT( c.C(0,1,2,u.X()) == g );
      u.SWAP();
      CPPUNIT_ASSERT( u*cnot*u == g );
    }

    void testS()
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

      CPPUNIT_ASSERT( s.S(0,1,2)*ab == ba );
      CPPUNIT_ASSERT( s.S(0,1,2) == t.SWAP() );

      // Permute three qubits.
      Qubit abc, acb, bac, bca, cab, cba;
      abc = a.tensorDot(b.tensorDot(c));
      acb = a.tensorDot(c.tensorDot(b));
      bac = b.tensorDot(a.tensorDot(c));
      bca = b.tensorDot(c.tensorDot(a));
      cab = c.tensorDot(a.tensorDot(b));
      cba = c.tensorDot(b.tensorDot(a));

      CPPUNIT_ASSERT( s.S(0,0,3)*abc == abc );
      CPPUNIT_ASSERT( s.S(0,1,3)*abc == bac );
      CPPUNIT_ASSERT( s.S(1,2,3)*abc == acb );
      CPPUNIT_ASSERT( s.S(0,2,3)*abc == cba );
      CPPUNIT_ASSERT( s.S(1,2,3)*t.S(0,1,3)*abc == bca );
      CPPUNIT_ASSERT( s.S(0,1,3)*t.S(1,2,3)*abc == cab );

      std::vector<int> p1(3), p2(3);
      p1[0] = 1; p2[0] = 2;
      p1[1] = 2; p2[1] = 0;
      p1[2] = 0; p2[2] = 1;

      CPPUNIT_ASSERT( s.S(p1)*abc == bca );
      CPPUNIT_ASSERT( s.S(p2)*abc == cab );

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

      CPPUNIT_ASSERT( s.S(0,1,4)*abcd == bacd );
      CPPUNIT_ASSERT( s.S(0,2,4)*abcd == cbad );
      CPPUNIT_ASSERT( s.S(0,3,4)*abcd == dbca );
      CPPUNIT_ASSERT( s.S(1,2,4)*abcd == acbd );
      CPPUNIT_ASSERT( s.S(1,3,4)*abcd == adcb );
      CPPUNIT_ASSERT( s.S(2,3,4)*abcd == abdc );

      std::vector<int> p3(4), p4(4), p5(4), p6(4);
      p3[0] = 3; p4[0] = 1; p5[0] = 3; p6[0] = 3;
      p3[1] = 2; p4[1] = 2; p5[1] = 0; p6[1] = 2;
      p3[2] = 1; p4[2] = 3; p5[2] = 1; p6[2] = 0;
      p3[3] = 0; p4[3] = 0; p5[3] = 2; p6[3] = 1;

      CPPUNIT_ASSERT( s.S(p3)*abcd == dcba );
      CPPUNIT_ASSERT( s.S(p3).transpose()*dcba == abcd );
      CPPUNIT_ASSERT( s.S(p4)*abcd == bcda );
      CPPUNIT_ASSERT( s.S(p4).transpose()*bcda == abcd );
      CPPUNIT_ASSERT( s.S(p5)*abcd == dabc );
      CPPUNIT_ASSERT( s.S(p5).transpose()*dabc == abcd );
      CPPUNIT_ASSERT( s.S(p6)*abcd == dcab );
      CPPUNIT_ASSERT( s.S(p6).transpose()*dcab == abcd );

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

      CPPUNIT_ASSERT( s.S(r1)*abcde == acbed );
      CPPUNIT_ASSERT( s.S(r2)*abcde == edcab );

      CPPUNIT_ASSERT( s.S(0,1,5)*abcde == bacde );
      CPPUNIT_ASSERT( s.S(0,2,5)*abcde == cbade );
      CPPUNIT_ASSERT( s.S(0,3,5)*abcde == dbcae );
      CPPUNIT_ASSERT( s.S(0,4,5)*abcde == ebcda );

      CPPUNIT_ASSERT( s.S(1,2,5)*abcde == acbde );
      CPPUNIT_ASSERT( s.S(1,3,5)*abcde == adcbe );
      CPPUNIT_ASSERT( s.S(1,4,5)*abcde == aecdb );

      CPPUNIT_ASSERT( s.S(2,3,5)*abcde == abdce );
      CPPUNIT_ASSERT( s.S(2,4,5)*abcde == abedc );

      CPPUNIT_ASSERT( s.S(3,4,5)*abcde == abced );
    }

    void testU()
    {
      Gate uf, uc, uh, g[16];

      //
      // Test all 4x4 U gates.
      //
      std::vector<int> f(2);

      f[0] = 0; f[1] = 0;
      uh.resize(4,4);
      uh << 1, 0, 0, 0,
            0, 1, 0, 0,
            0, 0, 1, 0,
            0, 0, 0, 1;
      uc.I().tensorPowSet(2);
      uf.U(f);
      CPPUNIT_ASSERT( uh == uc );
      CPPUNIT_ASSERT( uh == uf );
      CPPUNIT_ASSERT( uh == uf.U(f,1) );

      f[0] = 0; f[1] = 1;
      uh.resize(4,4);
      uh << 1, 0, 0, 0,
            0, 1, 0, 0,
            0, 0, 0, 1,
            0, 0, 1, 0;
      uc.C(1,0,2,g[0].X());
      uf.U(f);
      CPPUNIT_ASSERT( uh == uc );
      CPPUNIT_ASSERT( uh == uf );
      CPPUNIT_ASSERT( uh == uf.U(f,1) );

      f[0] = 1; f[1] = 0;
      uh.resize(4,4);
      uh << 0, 1, 0, 0,
            1, 0, 0, 0,
            0, 0, 1, 0,
            0, 0, 0, 1;
      uc = g[0].X().applyTo(1,2) * g[1].C(1,0,2,g[2].X());
      uf.U(f);
      CPPUNIT_ASSERT( uh == uc );
      CPPUNIT_ASSERT( uh == uf );
      CPPUNIT_ASSERT( uh == uf.U(f,1) );

      f[0] = 1; f[1] = 1;
      uh.resize(4,4);
      uh << 0, 1, 0, 0,
            1, 0, 0, 0,
            0, 0, 0, 1,
            0, 0, 1, 0;
      uc.X().applyToSet(1,2);
      uf.U(f);
      CPPUNIT_ASSERT( uh == uc );
      CPPUNIT_ASSERT( uh == uf );
      CPPUNIT_ASSERT( uh == uf.U(f,1) );

      //
      // Test all 8x8 U gates.
      //
      f.resize(4);

      // binary 0
      f[0] = 0; f[1] = 0; f[2] = 0; f[3] = 0;
      uh.resize(8,8);
      uh.setIdentity();
      uc.I().tensorPowSet(3);
      uf.U(f);
      CPPUNIT_ASSERT( uh == uc );
      CPPUNIT_ASSERT( uh == uf );
      CPPUNIT_ASSERT( uh == uf.U(f,1) );

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
      uc.X().applyToSet(2,3);
      uf.U(f);
      CPPUNIT_ASSERT( uh == uc );
      CPPUNIT_ASSERT( uh == uf );
      CPPUNIT_ASSERT( uh == uf.U(f,1) );

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
      uc.C(1,0,3,g[0].CNOT());
      uf.U(f);
      CPPUNIT_ASSERT( uh == uc );
      CPPUNIT_ASSERT( uh == uf );
      CPPUNIT_ASSERT( uh == uf.U(f,1) );

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
      uc = g[0].X().applyTo(2,3) * g[1].C(1,0,3,g[2].CNOT());
      uf.U(f);
      CPPUNIT_ASSERT( uh == uc );
      CPPUNIT_ASSERT( uh == uf );
      CPPUNIT_ASSERT( uh == uf.U(f,1) );

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
      uc = g[0].X().applyTo(1,3) * g[1].C(1,0,3,g[2].CNOT()) *
           g[3].X().applyTo(1,3);
      uf.U(f);
      CPPUNIT_ASSERT( uh == uc );
      CPPUNIT_ASSERT( uh == uf );
      CPPUNIT_ASSERT( uh == uf.U(f,1) );

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
      uc = g[0].X().applyTo(1,3) * g[1].X().applyTo(2,3) *
           g[2].C(1,0,3,g[3].CNOT()) * g[4].X().applyTo(1,3);
      uf.U(f);
      CPPUNIT_ASSERT( uh == uc );
      CPPUNIT_ASSERT( uh == uf );
      CPPUNIT_ASSERT( uh == uf.U(f,1) );

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
      uc = g[0].C(2,0,3,g[1].X());
      uf.U(f);
      CPPUNIT_ASSERT( uh == uc );
      CPPUNIT_ASSERT( uh == uf );
      CPPUNIT_ASSERT( uh == uf.U(f,1) );

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
      uc = g[0].X().applyTo(0,3) * g[1].C(2,0,3,g[2].X()) *
           g[3].X().applyTo(0,3);
      uf.U(f);
      CPPUNIT_ASSERT( uh == uc );
      CPPUNIT_ASSERT( uh == uf );
      CPPUNIT_ASSERT( uh == uf.U(f,1) );

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
      uc = g[0].X().applyTo(0,3) * g[1].CCNOT() * g[2].X().applyTo(0,3);
      uf.U(f);
      CPPUNIT_ASSERT( uh == uc );
      CPPUNIT_ASSERT( uh == uf );
      CPPUNIT_ASSERT( uh == uf.U(f,1) );

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
      uc = g[0].X().applyTo(0,3) * g[1].X().applyTo(2,3) * g[2].CCNOT() *
           g[3].X().applyTo(0,3);
      uf.U(f);
      CPPUNIT_ASSERT( uh == uc );
      CPPUNIT_ASSERT( uh == uf );
      CPPUNIT_ASSERT( uh == uf.U(f,1) );

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
      uc.CNOT().applyToSet(1,3);
      uf.U(f);
      CPPUNIT_ASSERT( uh == uc );
      CPPUNIT_ASSERT( uh == uf );
      CPPUNIT_ASSERT( uh == uf.U(f,1) );

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
      uc = g[0].X().applyTo(1,3) * g[1].CNOT().applyTo(1,3) *
           g[2].X().applyTo(1,3);
      uf.U(f);
      CPPUNIT_ASSERT( uh == uc );
      CPPUNIT_ASSERT( uh == uf );
      CPPUNIT_ASSERT( uh == uf.U(f,1) );

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
      uc = g[0].C(2,1,3,g[1].X()) * g[2].C(2,0,3,g[3].X());
      uf.U(f);
      CPPUNIT_ASSERT( uh == uc );
      CPPUNIT_ASSERT( uh == uf );
      CPPUNIT_ASSERT( uh == uf.U(f,1) );

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
      uc = g[0].X().applyTo(2,3) * g[1].C(2,1,3,g[2].X()) *
           g[3].C(2,0,3,g[4].X());
      uf.U(f);
      CPPUNIT_ASSERT( uh == uc );
      CPPUNIT_ASSERT( uh == uf );
      CPPUNIT_ASSERT( uh == uf.U(f,1) );

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
      uc = g[0].X().applyTo(0,3) * g[1].X().applyTo(1,3) *
           g[2].X().applyTo(2,3) * g[3].CCNOT() * g[4].X().applyTo(0,3) *
           g[5].X().applyTo(1,3);
      uf.U(f);
      CPPUNIT_ASSERT( uh == uc );
      CPPUNIT_ASSERT( uh == uf );
      CPPUNIT_ASSERT( uh == uf.U(f,1) );

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
      uc = g[0].X().applyTo(0,3) * g[1].X().applyTo(1,3) * g[2].CCNOT() *
           g[3].X().applyTo(0,3) * g[4].X().applyTo(1,3);
      uf.U(f);
      CPPUNIT_ASSERT( uh == uc );
      CPPUNIT_ASSERT( uh == uf );
      CPPUNIT_ASSERT( uh == uf.U(f,1) );

      //
      // Test all 8x8 U gates for two output qubits.
      //
      f.resize(2);

      f[0] = 0; f[1] = 0;
      uh.resize(8,8);
      uh.setIdentity();
      uc.I().tensorPowSet(3);
      uf.U(f,2);
      CPPUNIT_ASSERT( uh == uc );
      CPPUNIT_ASSERT( uh == uf );

      f[0] = 0; f[1] = 1;
      uh.resize(8,8);
      uh << 1, 0, 0, 0, 0, 0, 0, 0,
            0, 1, 0, 0, 0, 0, 0, 0,
            0, 0, 1, 0, 0, 0, 0, 0,
            0, 0, 0, 1, 0, 0, 0, 0,
            0, 0, 0, 0, 0, 1, 0, 0,
            0, 0, 0, 0, 1, 0, 0, 0,
            0, 0, 0, 0, 0, 0, 0, 1,
            0, 0, 0, 0, 0, 0, 1, 0;
      uc = g[0].C(2,0,3,g[1].X());
      uf.U(f,2);
      CPPUNIT_ASSERT( uh == uc );
      CPPUNIT_ASSERT( uh == uf );

      f[0] = 0; f[1] = 2;
      uh.resize(8,8);
      uh << 1, 0, 0, 0, 0, 0, 0, 0,
            0, 1, 0, 0, 0, 0, 0, 0,
            0, 0, 1, 0, 0, 0, 0, 0,
            0, 0, 0, 1, 0, 0, 0, 0,
            0, 0, 0, 0, 0, 0, 1, 0,
            0, 0, 0, 0, 0, 0, 0, 1,
            0, 0, 0, 0, 1, 0, 0, 0,
            0, 0, 0, 0, 0, 1, 0, 0;
      uc = g[0].C(1,0,3,g[1].X());
      uf.U(f,2);
      CPPUNIT_ASSERT( uh == uc );
      CPPUNIT_ASSERT( uh == uf );

      f[0] = 0; f[1] = 3;
      uh.resize(8,8);
      uh << 1, 0, 0, 0, 0, 0, 0, 0,
            0, 1, 0, 0, 0, 0, 0, 0,
            0, 0, 1, 0, 0, 0, 0, 0,
            0, 0, 0, 1, 0, 0, 0, 0,
            0, 0, 0, 0, 0, 0, 0, 1,
            0, 0, 0, 0, 0, 0, 1, 0,
            0, 0, 0, 0, 0, 1, 0, 0,
            0, 0, 0, 0, 1, 0, 0, 0;
      uc = g[0].C(2,0,3,g[1].X()) * g[2].C(1,0,3,g[3].X());
      uf.U(f,2);
      CPPUNIT_ASSERT( uh == uc );
      CPPUNIT_ASSERT( uh == uf );

      f[0] = 1; f[1] = 0;
      uh.resize(8,8);
      uh << 0, 1, 0, 0, 0, 0, 0, 0,
            1, 0, 0, 0, 0, 0, 0, 0,
            0, 0, 0, 1, 0, 0, 0, 0,
            0, 0, 1, 0, 0, 0, 0, 0,
            0, 0, 0, 0, 1, 0, 0, 0,
            0, 0, 0, 0, 0, 1, 0, 0,
            0, 0, 0, 0, 0, 0, 1, 0,
            0, 0, 0, 0, 0, 0, 0, 1;
      uc = g[0].X().applyTo(0,3) * g[1].C(2,0,3,g[2].X()) *
           g[3].X().applyTo(0,3);
      uf.U(f,2);
      CPPUNIT_ASSERT( uh == uc );
      CPPUNIT_ASSERT( uh == uf );

      f[0] = 1; f[1] = 1;
      uh.resize(8,8);
      uh << 0, 1, 0, 0, 0, 0, 0, 0,
            1, 0, 0, 0, 0, 0, 0, 0,
            0, 0, 0, 1, 0, 0, 0, 0,
            0, 0, 1, 0, 0, 0, 0, 0,
            0, 0, 0, 0, 0, 1, 0, 0,
            0, 0, 0, 0, 1, 0, 0, 0,
            0, 0, 0, 0, 0, 0, 0, 1,
            0, 0, 0, 0, 0, 0, 1, 0;
      uc = g[0].C(2,0,3,g[1].X()) * g[2].X().applyTo(0,3) *
           g[3].C(2,0,3,g[4].X()) * g[5].X().applyTo(0,3);
      uf.U(f,2);
      CPPUNIT_ASSERT( uh == uc );
      CPPUNIT_ASSERT( uh == g[6].X().applyTo(2,3) );
      CPPUNIT_ASSERT( uh == uf );

      f[0] = 1; f[1] = 2;
      uh.resize(8,8);
      uh << 0, 1, 0, 0, 0, 0, 0, 0,
            1, 0, 0, 0, 0, 0, 0, 0,
            0, 0, 0, 1, 0, 0, 0, 0,
            0, 0, 1, 0, 0, 0, 0, 0,
            0, 0, 0, 0, 0, 0, 1, 0,
            0, 0, 0, 0, 0, 0, 0, 1,
            0, 0, 0, 0, 1, 0, 0, 0,
            0, 0, 0, 0, 0, 1, 0, 0;
      uc = g[0].C(1,0,3,g[1].X()) * g[2].X().applyTo(0,3) *
           g[3].C(2,0,3,g[4].X()) * g[5].X().applyTo(0,3);
      uf.U(f,2);
      CPPUNIT_ASSERT( uh == uc );
      CPPUNIT_ASSERT( uh == uf );

      f[0] = 1; f[1] = 3;
      uh.resize(8,8);
      uh << 0, 1, 0, 0, 0, 0, 0, 0,
            1, 0, 0, 0, 0, 0, 0, 0,
            0, 0, 0, 1, 0, 0, 0, 0,
            0, 0, 1, 0, 0, 0, 0, 0,
            0, 0, 0, 0, 0, 0, 0, 1,
            0, 0, 0, 0, 0, 0, 1, 0,
            0, 0, 0, 0, 0, 1, 0, 0,
            0, 0, 0, 0, 1, 0, 0, 0;
      uc = g[0].C(2,0,3,g[1].X()) * g[2].C(1,0,3,g[3].X()) *
           g[4].X().applyTo(0,3) * g[5].C(2,0,3,g[6].X()) *
           g[7].X().applyTo(0,3);
      uf.U(f,2);
      CPPUNIT_ASSERT( uh == uc );
      CPPUNIT_ASSERT( uh == uf );

      f[0] = 2; f[1] = 0;
      uh.resize(8,8);
      uh << 0, 0, 1, 0, 0, 0, 0, 0,
            0, 0, 0, 1, 0, 0, 0, 0,
            1, 0, 0, 0, 0, 0, 0, 0,
            0, 1, 0, 0, 0, 0, 0, 0,
            0, 0, 0, 0, 1, 0, 0, 0,
            0, 0, 0, 0, 0, 1, 0, 0,
            0, 0, 0, 0, 0, 0, 1, 0,
            0, 0, 0, 0, 0, 0, 0, 1;
      uc = g[0].X().applyTo(0,3) * g[1].C(1,0,3,g[2].X()) *
           g[3].X().applyTo(0,3);
      uf.U(f,2);
      CPPUNIT_ASSERT( uh == uc );
      CPPUNIT_ASSERT( uh == uf );

      f[0] = 2; f[1] = 1;
      uh.resize(8,8);
      uh << 0, 0, 1, 0, 0, 0, 0, 0,
            0, 0, 0, 1, 0, 0, 0, 0,
            1, 0, 0, 0, 0, 0, 0, 0,
            0, 1, 0, 0, 0, 0, 0, 0,
            0, 0, 0, 0, 0, 1, 0, 0,
            0, 0, 0, 0, 1, 0, 0, 0,
            0, 0, 0, 0, 0, 0, 0, 1,
            0, 0, 0, 0, 0, 0, 1, 0;
      uc = g[0].C(2,0,3,g[1].X()) * g[2].X().applyTo(0,3) *
           g[3].C(1,0,3,g[4].X()) * g[5].X().applyTo(0,3);
      uf.U(f,2);
      CPPUNIT_ASSERT( uh == uc );
      CPPUNIT_ASSERT( uh == uf );

      f[0] = 2; f[1] = 2;
      uh.resize(8,8);
      uh << 0, 0, 1, 0, 0, 0, 0, 0,
            0, 0, 0, 1, 0, 0, 0, 0,
            1, 0, 0, 0, 0, 0, 0, 0,
            0, 1, 0, 0, 0, 0, 0, 0,
            0, 0, 0, 0, 0, 0, 1, 0,
            0, 0, 0, 0, 0, 0, 0, 1,
            0, 0, 0, 0, 1, 0, 0, 0,
            0, 0, 0, 0, 0, 1, 0, 0;
      uc = g[0].C(1,0,3,g[1].X()) * g[2].X().applyTo(0,3) *
           g[3].C(1,0,3,g[4].X()) * g[5].X().applyTo(0,3);
      uf.U(f,2);
      CPPUNIT_ASSERT( uh == uc );
      CPPUNIT_ASSERT( uh == uf );

      f[0] = 2; f[1] = 3;
      uh.resize(8,8);
      uh << 0, 0, 1, 0, 0, 0, 0, 0,
            0, 0, 0, 1, 0, 0, 0, 0,
            1, 0, 0, 0, 0, 0, 0, 0,
            0, 1, 0, 0, 0, 0, 0, 0,
            0, 0, 0, 0, 0, 0, 0, 1,
            0, 0, 0, 0, 0, 0, 1, 0,
            0, 0, 0, 0, 0, 1, 0, 0,
            0, 0, 0, 0, 1, 0, 0, 0;
      uc = g[0].C(1,0,3,g[1].X()) * g[2].C(2,0,3,g[3].X()) *
           g[4].X().applyTo(0,3) * g[5].C(1,0,3,g[6].X()) *
           g[7].X().applyTo(0,3);
      uf.U(f,2);
      CPPUNIT_ASSERT( uh == uc );
      CPPUNIT_ASSERT( uh == uf );

      f[0] = 3; f[1] = 0;
      uh.resize(8,8);
      uh << 0, 0, 0, 1, 0, 0, 0, 0,
            0, 0, 1, 0, 0, 0, 0, 0,
            0, 1, 0, 0, 0, 0, 0, 0,
            1, 0, 0, 0, 0, 0, 0, 0,
            0, 0, 0, 0, 1, 0, 0, 0,
            0, 0, 0, 0, 0, 1, 0, 0,
            0, 0, 0, 0, 0, 0, 1, 0,
            0, 0, 0, 0, 0, 0, 0, 1;
      uc = g[0].X().applyTo(0,3) * g[1].C(1,0,3,g[2].X()) *
           g[3].C(2,0,3,g[4].X()) * g[5].X().applyTo(0,3);
      uf.U(f,2);
      CPPUNIT_ASSERT( uh == uc );
      CPPUNIT_ASSERT( uh == uf );

      f[0] = 3; f[1] = 1;
      uh.resize(8,8);
      uh << 0, 0, 0, 1, 0, 0, 0, 0,
            0, 0, 1, 0, 0, 0, 0, 0,
            0, 1, 0, 0, 0, 0, 0, 0,
            1, 0, 0, 0, 0, 0, 0, 0,
            0, 0, 0, 0, 0, 1, 0, 0,
            0, 0, 0, 0, 1, 0, 0, 0,
            0, 0, 0, 0, 0, 0, 0, 1,
            0, 0, 0, 0, 0, 0, 1, 0;
      uc = g[0].C(2,0,3,g[1].X()) * g[2].X().applyTo(0,3) *
           g[3].C(1,0,3,g[4].X()) * g[5].C(2,0,3,g[6].X()) *
           g[7].X().applyTo(0,3);
      uf.U(f,2);
      CPPUNIT_ASSERT( uh == uc );
      CPPUNIT_ASSERT( uh == uf );

      f[0] = 3; f[1] = 2;
      uh.resize(8,8);
      uh << 0, 0, 0, 1, 0, 0, 0, 0,
            0, 0, 1, 0, 0, 0, 0, 0,
            0, 1, 0, 0, 0, 0, 0, 0,
            1, 0, 0, 0, 0, 0, 0, 0,
            0, 0, 0, 0, 0, 0, 1, 0,
            0, 0, 0, 0, 0, 0, 0, 1,
            0, 0, 0, 0, 1, 0, 0, 0,
            0, 0, 0, 0, 0, 1, 0, 0;
      uc = g[0].C(1,0,3,g[1].X()) * g[2].X().applyTo(0,3) *
           g[3].C(1,0,3,g[4].X()) * g[5].C(2,0,3,g[6].X()) *
           g[7].X().applyTo(0,3);
      uf.U(f,2);
      CPPUNIT_ASSERT( uh == uc );
      CPPUNIT_ASSERT( uh == uf );

      f[0] = 3; f[1] = 3;
      uh.resize(8,8);
      uh << 0, 0, 0, 1, 0, 0, 0, 0,
            0, 0, 1, 0, 0, 0, 0, 0,
            0, 1, 0, 0, 0, 0, 0, 0,
            1, 0, 0, 0, 0, 0, 0, 0,
            0, 0, 0, 0, 0, 0, 0, 1,
            0, 0, 0, 0, 0, 0, 1, 0,
            0, 0, 0, 0, 0, 1, 0, 0,
            0, 0, 0, 0, 1, 0, 0, 0;
      uc = g[0].C(1,0,3,g[1].X()) * g[2].C(2,0,3,g[3].X()) *
           g[4].X().applyTo(0,3) * g[5].C(1,0,3,g[6].X()) *
           g[7].C(2,0,3,g[8].X()) * g[9].X().applyTo(0,3);
      uf.U(f,2);
      CPPUNIT_ASSERT( uh == uc );
      CPPUNIT_ASSERT( uh == uf );

      //
      // Test some 16x16 U gates for two input and output qubits.
      //
      f.resize(4);

      f[0] = 1; f[1] = 2; f[2] = 0; f[3] = 0;
      uh.resize(16,16);
      uh << 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
            1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
            0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
            0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
            0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0,
            0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0,
            0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
            0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
            0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0,
            0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0,
            0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0,
            0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0,
            0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0,
            0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0,
            0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0,
            0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1;
      // --X--o-----o--X--
      //      |     |
      // -----o--X--o--X--
      //      |     |
      // -----X-----|-----
      //            |
      // -----------X-----
      uc = g[0].X().applyTo(0,4) * g[1].C(1,0,4,g[2].CNOT()) *
           g[3].X().applyTo(1,4) * g[4].C(1,0,4,g[5].C(2,0,3,g[6].X())) *
           g[7].X().tensorPow(2).applyTo(0,4);
      uf.U(f,2);
      CPPUNIT_ASSERT( uh == uc );
      CPPUNIT_ASSERT( uh == uf );
    }

    void testF()
    {
      field i = field(0,1);
      Qubit q0(0,4),   q1(1,4),   q2(2,4),   q3(3,4),
            q4(4,4),   q5(5,4),   q6(6,4),   q7(7,4),
            q8(8,4),   q9(9,4),   q10(10,4), q11(11,4),
            q12(12,4), q13(13,4), q14(14,4), q15(15,4);

      Gate f;
      f.F(4);

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
      Gate g[8];

      CPPUNIT_ASSERT( ( g[0].H() * g[1].H() ).isApprox(g[2].I()) );

      CPPUNIT_ASSERT( ( g[0].H() * g[1].X() * g[2].H() ).isApprox(g[3].Z()) );

      CPPUNIT_ASSERT( ( g[0].H() * g[1].Z() * g[2].H() ).isApprox(g[3].X()) );

      CPPUNIT_ASSERT( ( g[0].H().applyTo(1,2) * g[1].C(1,0,2,g[2].Z()) *
        g[3].H().applyTo(1,2) ).isApprox(g[4].CNOT()) );

      CPPUNIT_ASSERT( g[0].C(1,0,2,g[1].Z()).isApprox(g[2].C(0,1,2,g[3].Z())) );

      CPPUNIT_ASSERT( ( g[0].H().tensorPow(2) * g[1].CNOT() *
        g[2].H().tensorPow(2) ).isApprox(g[3].C(0,1,2,g[4].X())) );
    }
};

} // namespace QuCoSi

#endif // QUCOSI_GATETEST_H

// vim: shiftwidth=2 textwidth=78

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

#include <bitset>
#include <cstdlib>
#include <ctime>
#include <limits>

#include <QuCoSi/Aux>
#include <QuCoSi/Basis>
#include <QuCoSi/Gate>
#include <QuCoSi/Qubit>
#include <QuCoSi/Vector>

using namespace std;
using namespace Eigen;
using namespace QuCoSi;

int main(int, char *[])
{
  std::srand((unsigned)std::time(NULL) + (unsigned)std::clock());

  Gate x;
  for(int i = 0; i < 5000000; i++) {
    x.CNOT();
  }

/*
  Gate u;
  vector<int> f(4);
  f[0] = 1; f[1] = 1; f[2] = 1; f[3] = 1;
  for(int i = 0; i < 1000000; i++) {
    u.U2(f);
  }
*/

/*
  cout << bwise_bin_add(0,1) << (1^0) << endl;
  cout << bwise_bin_add(1,1) << (1^1) << endl;
  cout << bwise_bin_add(2,1) << (1^2) << endl;
  cout << bwise_bin_add(3,1) << (1^3) << endl;

  cout << bwise_bin_dot(25,3) << endl;
*/

/*
  Gate c,h;
  h.XGate();
  cout << c.CGate(0,1,2,h);
*/

/*
  cout << (0+3)%4 << endl;
  cout << (1+3)%4 << endl;
  cout << (2+3)%4 << endl;
  cout << (3+3)%4 << endl;
*/

/*
  Gate s,c;
  s.SWAPGate();
  c.CNOTGate();

  cout << s*c*s << endl;
*/

/*
  bitset<std::numeric_limits<unsigned long>::digits>
    c(10), d(0);
  cout << c << endl;

  vector<int> p(4);
  int n = p.size();
  p[0] = 1;
  p[1] = 2;
  p[2] = 3;
  p[3] = 0;

  for(int i = 0; i < n; i++) {
    int o = n-1-i;
    int u = n-1-p[i];

    cout << o <<  " " << u << endl;
    d[o] = c[u];
  }
  cout << d << endl;
*/

/*
  // passt fuer alle permutationen
  vector<int> p(3);
  int n = p.size();
  p[0] = 2;
  p[1] = 0;
  p[2] = 1;

  for(int i = 0; i < n; i++) {
    int o = n-1-i;
    int u = n-1-p[i];

    cout << o <<  " " << u << endl;
    d[o] = c[u];
  }
  cout << d << endl;
*/

/*
  //Gate s;
  //cout << s.SGate(p) << endl << endl;

  Gate c0,c1,x;
  c0.CNOTGate();
  c1.CGate(x.XGate());
  cout << c0 << endl;
  cout <<( c1 == c0) << endl;
*/

/*
  Gate f;
  f.FGate(2);
  cout << f << endl;
  cout << (f*f.conjugate()).isIdentity() << endl;
*/

/*
  Qubit q0(field(1,0), field(0,0));
  Qubit q1(field(0,0), field(1,0));
  Qubit xy;

  xy = q0.tensorDot(q1);
  Gate h; h.HGate().tensorPowSet(2);
  xy = h*xy;

  vector<int> f(2);
  f[0] = 1;
  f[1] = 0;
  Gate u; u.UfGate(xy,f);
  xy = u*xy;
  xy = h*xy;
  cout << xy << endl;
*/

/*
  Qubit q0(field(1,0), field(0,0));
  Qubit q1(field(0,0), field(1,0));
  Qubit xy, xy2;
  Gate h;
  h.HGate().tensorPowSet(2);
  Gate u;

  xy.resize(4);
  xy2.resize(4);
  xy2.setZero();

  vector<int> f(2);
  f[0] = 0;
  f[1] = 0;

  //xy = q0.tensorDot(q1);
  //cout << xy << endl << endl;
  //cout << h << endl << endl;
  //xy = h*xy;
  xy(0) = 2;
  xy(1) = 3;
  xy(2) = 5;
  xy(3) = 7;
  cout << xy << endl << endl;

  Basis b(xy.size());
  //cout << b << b.size()<< endl;
  Gate uu,uu2;
  uu2.resize(4,4);

  cout << u.UfGate(xy,f) << endl;


  for (int i = 0; i < b.size(); i++) {
    Qubit tmp;
    tmp = b[i].dot(xy)*b[i];
    //cout << tmp << endl;
    uu.UfGate2(tmp,f);
    uu2.col(i) = uu.col(i);
    xy2 = xy2 + uu*tmp;
  }
  cout << xy2 << endl << endl;
  cout << u*xy << endl;

  //xy = u.UfGate(xy,f)*xy;
  //cout << xy << endl << endl;

  //xy2 = h*xy2;
  //cout << xy2 << endl << endl;
*/

/*
  cout << (1+1)%2 << endl;
  cout << (1+0)%2 << endl;
  cout << (0+1)%2 << endl;
  cout << (0+0)%2 << endl;
*/

/*
  Qubit x(field(0,0), field(1,0));
  Gate h;
  h.HGate();
  x = h*x;
  cout << x.measure() << endl;
*/

/*
  Gate m;
  cout << m.CNOTGate() << endl;
*/

/*
  Qubit q0(field(1,0), field(0,0)),
        q1(field(0,0), field(1,0));
  XGate x = XGate();
  YGate y = YGate();
  ZGate z = ZGate();
  HGate h = HGate();
  PGate p = PGate();
  TGate t = TGate();

  cout << h.extendToPos(1,2) << endl;

  cout << x*q0 << endl << endl;
  cout << x*q1 << endl << "---" << endl;

  cout << y*q0 << endl << endl;
  cout << y*q1 << endl << "---" << endl;

  cout << z*q0 << endl << endl;
  cout << z*q1 << endl << "---" << endl;

  cout << h*q0 << endl << endl;
  cout << h*q1 << endl << "---" << endl;

  cout << p*q0 << endl << endl;
  cout << p*q1 << endl << "---" << endl;

  cout << t*q0 << endl << endl;
  cout << t*q1 << endl << "---" << endl;
*/

/*
  Gate m(4,4);

  m.setZero();
  m(0,0) = 1;
  m(2,1) = 1;
  m(1,2) = 1;
  m(3,3) = 1;

  cout << m.extendToPos(1,3);
*/

/*
  Gate m(4,4), i(2,2);

  m.setZero();
  m(0,0) = 1;
  m(2,1) = 1;
  m(1,2) = 1;
  m(3,3) = 1;
  i.setIdentity();

  cout << i.tensorPow(0) << endl << endl;
  cout << m << endl << endl;

  cout << m.tensorDot(i);
*/

/*
    Gate m1(2,2), m2(2,2);
    m1 << 1,2,3,4;
    m2 << 5,6,7,8;
    cout << m1 << endl << endl << m2 << endl << endl;
    cout << m1.tensorDot(m2);
*/

/*
    MatrixXf m1(2,2), m2(2,2);
    m1 << 1,2,3,4;
    m2 << 5,6,7,8;
    cout << m1 << endl << endl << m2 << endl << endl;

    MatrixXf m3(m1.rows()*m2.rows(), m1.cols()*m2.cols());

    for (int i = 0; i < m1.cols(); i++) {
        for (int j = 0; j < m1.rows(); j++) {
            m3.block(i*m2.rows(), j*m2.cols(), m2.rows(), m2.cols()) =  m1(i,j)*m2;
        }
    }

    cout << m3 << endl;
*/

/*
    Qubit
        q1(field(0,0), field(1,0)),
        q2(field(0,0), field(1,0)),
        q3(field(1,0), field(0,0));

    cout << q1 << endl << endl;
    cout << q2 << endl << endl;
    cout << q3 << endl << endl;
    cout << q1.otimesSet(q2).otimesSet(q3) << endl << endl;
    //cout << q1.getStdBasis() << endl << endl;
*/

/*
    Vector v1(field(1,0), field(0,1)),v2(2);
    Qubit q;
    Basis b(2);
    q=v1;
    cout << b << endl << endl;
    cout << q.getStdBasis() << endl;
    v2=(field(3,4)* b[0]);
    cout << field(55,0)*b[1] << endl;
*/

/*
    Basis b1(2), b2(2);
    cout << b1 << endl << endl;
    cout << b2 << endl << endl;
    cout << b1.otimesSet(b2) << endl;
    cout << b1.setNaturalBasis(0) << endl;
*/

/*
    Vector v(field(0.9999999999999999,0), field(0,0));
    cout << v.isNormalized() << endl;
*/

/*
    Qubit q;
    MultiQubit mq;

    mq.otimes(q);

    cout << q << endl<<endl;
    cout << mq << endl<<endl;
    mq.otimes(q);
    cout << mq << endl<<endl;
    mq.otimes(mq);
    cout << mq << endl<<endl;
*/

/*
    std::srand((unsigned)std::time(NULL) + (unsigned)std::clock());

    Qubit k, q(std::complex<float>(0.707106781187,0),
               std::complex<float>(0,0.707106781187));


    k.setRandom();
    cout << k << endl;
    cout << k.measure() << endl;
    cout << k.setRandom() << endl;
    cout << k.measure() << endl<< endl;
    cout <<  std::pow(std::numeric_limits<float>::epsilon(),2);
*/

/*
    cout << q[0] <<endl;
    cout << k.isPureState() << endl;
    cout << q.isPureState() << endl;
    cout << k << endl << q << endl;
    cout << pow(abs(k.dot(q)),2) << endl;
    cout << q.coefficient(1) << endl;
    cout << (q.norm()-1) << ( (q.norm()-1) == 0 ) << endl;
*/

/*
    cout << 0%2 << 1%2 << 2%2 << 3%2 << 4%2<< 5%2 <<6%2;
    q.setRandom();
    std::cout << q << endl <<  abs(q.norm()-1)  << abs(-1.23)<<endl;
    std::cout << q << endl << ( (q.norm()-1)  < numeric_limits<float>::epsilon()) << endl;
*/

/*
    Eigen::Vector2cf v1;
    v1[0] = std::complex<float>(1,0);
    v1[1] = std::complex<float>(0,1);

    std::cout << v1 << std::endl << v2.norm() << v1[1].real() << std::endl;
    std::cout << v1.isOrthogonal(v2) << v1.isOnes();
    std::cout << q.isNormalized() << std::endl << std::endl;
    q.setStdBase(v1, v2);
*/

  return 0;
}

// vim: shiftwidth=2 textwidth=78

#include <cstdlib>
#include <ctime>
#include <limits>

#include <Qubit.h>
#include <Vector.h>
#include <Basis.h>

using namespace std;
using namespace QuCoSi;

int main(int, char *[])
{
    Qubit 
        q1(field(0,0), field(1,0)),
        q2(field(0,0), field(1,0)),
        q3(field(1,0), field(0,0));

    cout << q1 << endl << endl;
    cout << q2 << endl << endl;
    cout << q3 << endl << endl;
    cout << q1.otimesSet(q2).otimesSet(q3) << endl << endl;
    //cout << q1.getStdBasis() << endl << endl;

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
//    cout << 0%2 << 1%2 << 2%2 << 3%2 << 4%2<< 5%2 <<6%2;
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

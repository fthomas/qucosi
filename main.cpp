#include <Qubit.h>
#include <cstdlib>
#include <ctime>
#include <limits>

using namespace std;

int main(int, char *[])
{
    std::srand((unsigned)std::time(NULL) + (unsigned)std::clock());

    Qubit k, q(std::complex<float>(0.707106781187,0),
               std::complex<float>(0,0.707106781187));


    k.setRandom();
    cout << k << endl;
    cout << k.measure() << endl;
    cout << k.setRandom() << endl;
    cout << k.measure() << endl<< endl;
    cout <<  std::pow(std::numeric_limits<float>::epsilon(),2);

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

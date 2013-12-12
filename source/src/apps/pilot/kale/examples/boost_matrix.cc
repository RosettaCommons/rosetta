#include <iostream>

#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/io.hpp>

#include <numeric/types.hh>
#include <numeric/ublas_extras/utilities.hh>

using namespace std;
using namespace numeric;
using namespace boost::numeric;
using namespace numeric::ublas_extras;

int main (int argc, char** argv) {

	ublas::matrix<Real> a (3, 3);
	ublas::matrix<Real> b (4, 4);

	a(0,0) = 3;  a(0,1) = 5;  a(0,2) = 5;
	a(1,0) = 2;  a(1,1) = 3;  a(1,2) = 4;
	a(2,0) = 3;  a(2,1) = 4;  a(2,2) = 1;

	b(0,0) = 3;  b(0,1) = 4;  b(0,2) = 5;  b(0,3) = 4;
	b(1,0) = 3;  b(1,1) = 4;  b(1,2) = 2;  b(1,3) = 2;
	b(2,0) = 1;  b(2,1) = 3;  b(2,2) = 4;  b(2,3) = 2;
	b(3,0) = 4;  b(3,1) = 3;  b(3,2) = 3;  b(3,3) = 2;

	/*
	a(0,0) =  2;   a(0,1) = -3;   a(0,2) = -2;
	a(1,0) = -6;   a(1,1) =  3;   a(1,2) =  3;
	a(2,0) = -2;   a(2,1) = -3;   a(2,2) = -2;

	b(0,0) = -4;   b(0,1) =  5;   b(0,2) =  2;
	b(1,0) = -3;   b(1,1) =  4;   b(1,2) =  2;
	b(2,0) = -1;   b(2,1) =  2;   b(2,2) =  5;

	cout << "Original Matrices" << endl;
	cout << a << endl;
	cout << b << endl << endl;

	cout << "Arithmetic Operators:" << endl;
	cout << a + b << " (a + b)" << endl;
	cout << a - b << " (a - b)" << endl;
	cout << a * 2 << " (a * b)" << endl;
	cout << a / 2 << " (a / b)" << endl << endl;
	*/

	cout << "Linear Algebra:" << endl;
	cout << determinant(a) << " (det |A| = 12)" << endl;
	cout << determinant(b) << " (det |B| = -3)" << endl;
}

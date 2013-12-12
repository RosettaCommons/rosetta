#include <iostream>

#include <utility/vector1.hh>
#include <numeric/xyzVector.hh>
#include <numeric/MathVector.hh>
#include <numeric/MathVector_operations.hh>

using namespace std;
using namespace numeric;

void print_vector(MathVector<Real> const &vector) {
	cout << "( ";
	cout << vector(0) << ", ";
	cout << vector(1) << ", ";
	cout << vector(2) << ", ";
	cout << vector(3) << ", ";
	cout << vector(4);
	cout << " )" << endl;
}

void print_vector(xyzVector<Real> const &vector) {
	cout << "< ";
	cout << vector.x() << ", ";
	cout << vector.y() << ", ";
	cout << vector.z();
	cout << " >" << endl;
}

void print_vector(utility::vector1<Real> const &vector) {
	cout << "[ ";
	cout << vector[0] << ", ";
	cout << vector[1] << ", ";
	cout << vector[2];
	cout << " ]" << endl;
}

int main(int argc, char** argv) {

	MathVector<Real> a(5), b(5), c(5);
	xyzVector<Real> r, s;
	utility::vector1<Real> m;

	// Initialize MathVectors.

	a(0) = Real(1.0);
	a(1) = Real(2.0);
	a(2) = Real(3.0);
	a(3) = Real(4.0);
	a(4) = Real(5.0);

	b(0) = Real(2.0);
	b(1) = Real(2.0);
	b(2) = Real(4.0);
	b(3) = Real(4.0);
	b(4) = Real(6.0);

	// Initialize vector1s.

	m.push_back(1.0);
	m.push_back(3.0);
	m.push_back(5.0);

	// Perform mathematical operations on MathVectors.

	print_vector(a);
	print_vector(b);
	print_vector(b - a);

	cout << a.norm() << endl;
	cout << a.square_norm() << endl << endl;

	// Construct xyzVectors from other vector types.
	
	r = a.begin();
	s = &(*m.begin());
	s = &(m.front());

	print_vector(r);
	print_vector(s);
	print_vector(r.cross(s));
}

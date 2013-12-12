#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/io.hpp>

using namespace std;

int main(int argc, char** argv) {

	namespace boost = boost::numeric::ublas;

	boost::vector<double> a (3);
	boost::vector<double> b (3);

	a(0) = 1; a(1) = 2; a(2) = 3;
	b(0) = 1; b(1) = 3; b(2) = 5;

	cout << "Original Vectors" << endl;
	cout << a << endl;
	cout << b << endl << endl;

	cout << "After a Swap:" << endl;
	a.swap(b);
	cout << a << endl;
	cout << b << endl << endl;

	cout << "Arithmetic Operators:" << endl;
	cout << a + b << " (a + b)" << endl;
	cout << a - b << " (a - b)" << endl;
	cout << a * 2 << " (a * b)" << endl;
	cout << a / 2 << " (a / b)" << endl << endl;

	cout << "Type Conversions:" << endl;

	// Note that attempting to extract a c-style array from a vector is a very 
	// dangerous thing to do.  It assumes that I know how the vector's data 
	// stored is arranged in memory (which I don't necessarily) and that the 
	// vector will never change how its data is stored (which it will, when it 
	// resizes).  But this approach is much faster than copying the data and is
	// probably good enough for short term conversions.

	double *data = &a(0);
	cout << data << " (data address)" << endl;
	cout << "<" << data[0] << ", " << data[1] << ", " << data[2] << "> ";
	cout << "(double array)" << endl;

}

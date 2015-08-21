#include <devel/optest/B2.hh>
//#include <devel/optest/A.hh>

#include <iostream>

namespace devel {
namespace optest {

void c2_func1() {
	std::cout << "BEGIN c2_func1: standard constructor / destructor" << std::endl;
	B2 b;
	std::cout << "b: "; b.status();
	std::cout << "END c2_func1" << std::endl << std::endl;
}

void c2_func2() {
	std::cout << "BEGIN c2_func2: set A" << std::endl;
	B2 b;
	std::cout << "b after default ctor: "; b.status();

	b.set_default();
	std::cout << "b after set_default(): "; b.status();
	std::cout << "END c2_func2" << std::endl << std::endl;
}

void c2_func3() {
	std::cout << "BEGIN c2_func3: copy constructor" << std::endl;
	B2 b1;
	std::cout << "b1: "; b1.status();

	B2 b2( b1 );
	std::cout << "b2: "; b2.status();

	b1.increment_a();
	std::cout << "b1: "; b1.status();
	std::cout << "b2: "; b2.status();
	std::cout << "END c2_func3" << std::endl << std::endl;
}

void c2_func4() {
	std::cout << "BEGIN c2_func4: assignment operator" << std::endl;
	B2 b1;
	std::cout << "b1: "; b1.status();

	B2 b2;
	std::cout << "b2: "; b2.status();

	b2 = b1;
	std::cout << "b2: "; b2.status();

	b1.increment_a();
	std::cout << "b1: "; b1.status();
	std::cout << "b2: "; b2.status();

	std::cout << "END c2_func4" << std::endl << std::endl;
}


void c2_func5() {
	std::cout << "BEGIN c2_func5: clone" << std::endl;
	B2 b1;
	std::cout << "b1: "; b1.status();

	B2OP b2op = b1.clone();
	std::cout << "b2op: "; b2op->status();

	b1.increment_a();
	std::cout << "b1: "; b1.status();
	std::cout << "b2op: "; b2op->status();

	std::cout << "END c2_func5" << std::endl << std::endl;
}


}
}

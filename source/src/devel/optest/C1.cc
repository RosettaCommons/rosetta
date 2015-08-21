#include <devel/optest/B1.hh>
//#include <devel/optest/A.hh>

#include <iostream>

namespace devel {
namespace optest {

void c1_func1() {
	std::cout << "BEGIN c1_func1: standard constructor / destructor" << std::endl;
	B1 b;
	std::cout << "b: "; b.status();
	std::cout << "END c1_func1" << std::endl << std::endl;
}

void c1_func2() {
	std::cout << "BEGIN c1_func2: set A pointer" << std::endl;
	B1 b;
	//b.set_aptr( new A );
	b.set_default_A();
	std::cout << "b: "; b.status();
	std::cout << "END c1_func2" << std::endl << std::endl;
}

void c1_func3() {
	std::cout << "BEGIN c1_func3: copy constructor" << std::endl;
	B1 b1;
	//b1.set_aptr( new A );
	b1.set_default_A();
	std::cout << "b1: "; b1.status();

	B1 b2( b1 );
	std::cout << "b2: "; b2.status();
	std::cout << "END c1_func3" << std::endl << std::endl;
}

void c1_func4() {
	std::cout << "BEGIN c1_func4: assignment operator" << std::endl;
	B1 b1;
	//b1.set_aptr( new A );
	b1.set_default_A();
	std::cout << "b1: "; b1.status();

	B1 b2;
	std::cout << "b2: "; b2.status();

	b2 = b1;
	std::cout << "b2: "; b2.status();

	std::cout << "END c1_func4" << std::endl << std::endl;
}


void c1_func5() {
	std::cout << "BEGIN c1_func5: clone" << std::endl;
	B1 b1;
	//b1.set_aptr( new A );
	b1.set_default_A();
	std::cout << "b1: "; b1.status();

	B1OP b2op = b1.clone();
	std::cout << "b2op: "; b2op->status();

	std::cout << "END c1_func5" << std::endl << std::endl;
}


}
}

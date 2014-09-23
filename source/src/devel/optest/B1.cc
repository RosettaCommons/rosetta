
#include <devel/optest/B1.hh>
#include <devel/optest/A.hh>


#include <iostream>

namespace devel {
namespace optest {

/// @details Auto-generated virtual destructor
B1::~B1() {}


/*
B1::B1() : a_pointer_( 0 ) {}
B1::~B1() {}
B1::B1( B1 const & src ) : a_pointer_( src.a_pointer_ ) {}
B1 const & B1::operator = ( B1 const & src )
{
	if ( this != & src ) {
		a_pointer_ = src.a_pointer_;
	}
	return *this;
}
*/


void B1::set_default_A() {
	a_pointer_ = AOP( new A );
}


void B1::set_aptr( AOP aptr ) {
	a_pointer_ = aptr;
}

void B1::status() {
   if ( a_pointer_ ) {
		std::cout << "B1::status -- a_pointer_ active, a_pointer_->my_int() == " << a_pointer_->my_int() << std::endl;
	} else {
		std::cout << "B1::status -- a_pointer_ points at 0" << std::endl;
	}
}

/*B1OP B1::clone() const {
	return new B1( *this );
}*/

AOP B1::get() {
	return a_pointer_;
}

}
}

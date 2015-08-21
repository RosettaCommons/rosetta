#include <devel/optest/B2.hh>

#include <iostream>

namespace devel {
namespace optest {

/// @details Auto-generated virtual destructor
B2::~B2() {}


/*
WARNING.  Changing from an instance to a pointer changes the semantics of these
functions.  What do you have to do in order to duplicate the old behavior?
Three of these four functions have to be changed.
B2::B2() : a_instance_() {}
B2::~B2() {}
B2::B2( B2 const & src ) : a_instance_( src.a_instance_ ) {}
B2 const & B2::operator = ( B2 const & rhs ) { a_instance_ = rhs.a_instance_; }
*/


/// this function has to change when a_instance_ is replaced by a_pointer_. What's the right way to do that?
void B2::set_default() {
	A anotherA;
	a_instance_ = anotherA;
}

/// this function has to change when a_instance_ is replaced by a_pointer_. What's the right way to do that?
void B2::set_a( A const & a ) {
	a_instance_ = a;
}

/// @brief Increment the value of integer that the a_instance_ holds
void B2::increment_a() {
	a_instance_.my_int( a_instance_.my_int() + 1 );
}

void B2::status() {

	std::cout << "B2::status -- My A has value " << a_instance_.my_int() << std::endl;

	/// When you've replaced a_instance_ with a_pointer_, use the code below to report
	/// B2's status instead of the code above

	/*if ( a_instance_ )  { /// NULL POINTER CHECK
	std::cout << "B2::status -- My A has value " << a_instance_->my_int() << std::endl;
	} else {
	std::cout << "B2::status -- My A is null " << std::endl;
	}*/
}

/*B2OP B2::clone() const {
return new B2( *this );
}*/

}
}

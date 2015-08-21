#ifndef INCLUDED_devel_optest_B1_hh
#define INCLUDED_devel_optest_B1_hh

#include <devel/optest/B1.fwd.hh>
#include <devel/optest/A.hh>  /// <--- Try to replace with devel/optest/A.fwd.hh

#include <utility/pointer/ReferenceCount.hh>

namespace devel {
namespace optest {

/// GOAL: replace A.hh with A.fwd.hh
///
/// What to do:
/// Uncomment the magic four functions below:
/// The default construtor, the destructor, the copy
/// constructor, and the assignment operator.
/// Uncomment these functions in the .cc file as well
/// Replace the #inclusion of A.hh with A.fwd.hh in this file
/// Try to compile C1.cc
///
/// ALSO TEST:
/// First:  Uncomment the four functions below from here and the .cc file.
/// and replace the #inclusion of A.hh with A.fwd.hh.  C1.cc will compile.
/// Second: Comment out B1's destructor both here and in the .cc file.
/// Third:  Try to compile C1.cc. Can you explain why C1.cc will not compile?
///
/// Try this one-at-a-time removal and readdition with the other
/// three functions.
///
/// ALSO TEST:
/// First:  Uncomment the four functions below from here and the .cc file
/// and replace the #inclusion of A.hh with A.fwd.hh.  C1.cc will compile.
/// Second: Restore the implementation of the "get()" function in B1 to
/// this header, and comment it out from B1.cc.
/// Third:  Try to compile C1.cc.  Can you explain why C1.cc will not compile?

class B1 : public utility::pointer::ReferenceCount {
public:
	/// @brief Automatically generated virtual destructor for class deriving directly from ReferenceCount
	virtual ~B1();

	/*B1();
	~B1();
	B1( B1 const & );
	B1 const & operator = ( B1 const & );*/

	void set_default_A();
	void set_aptr( AOP );
	void status();

	B1OP clone() const { return B1OP( new B1( *this ) ); }

	AOP get();// { return a_pointer_; }

private:
	AOP a_pointer_;

};

}
}

#endif

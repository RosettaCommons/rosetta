#ifndef INCLUDED_devel_optest_B2_hh
#define INCLUDED_devel_optest_B2_hh

#include <devel/optest/B2.fwd.hh>
#include <devel/optest/A.hh>   /// <---- Try to replace A.hh with A.fwd.hh

#include <utility/pointer/ReferenceCount.hh>

/// GOAL:  Replace the instance of class A in this class
/// with an owning pointer to class A, and remove A.hh from
/// this header.
///
/// What to do:
/// First, you need to keep save information about the old behavior.
/// Call this the "before" behavior.
/// Changing an instance of an object to a pointer of an object can
/// cause the class's behavior to change.
/// Build the pilot app andrew/testop.cc and the corresponding
/// exectuable once and save the output.
/// (e.g. ./bin/optest.macosgccdebug > optest_before.log )
///
/// Now you're ready to start modifying this class.
/// In this file, replace "A a_instance_;" with "AOP a_pointer_;"
/// In B2.cc, replace "a_instance_." with "a_pointer_->"
/// CAUTION: what happens if a_pointer_ is NULL?
/// Uncomment from the header the magic four functions:
/// the default constructor, the destructor, the copy
/// constructor, and the assignment operator.
/// Implement these four functions in the .cc file.
///
/// Now verify that the new behavior is the same as the old behavior.
/// Recompile the code and rerun the executable, saving the output
/// to a new file.  Call this the "after" behavior.
/// (e.g. ./bin/optest.macosgccdebug > optest_after.log).
/// Is the after behavior the same as the before behavior?
/// (e.g. diff optest_before.log optest_after.log)
///
/// Figure out what you have to do to make the outputs identical.
/// Email me if you'd like a hint (aleaverfay@gmail.com)
///
/// WARNING: if you're unable to make the after behavior match
/// the before behavior, then do not attempt to replace an object
/// with a pointer in real code.  Ask for help, but do not break
/// working code to improve compile time.

namespace devel {
namespace optest {

class B2 : public utility::pointer::ReferenceCount {
public:
	/// @brief Automatically generated virtual destructor for class deriving directly from ReferenceCount
	virtual ~B2();

	/* B2();
	~B2();
	B2( B2 const & );
	B2 const & operator = ( B2 );
	*/

	void set_default();
	void set_a( A const & a );
	void status();
	void increment_a();
	B2OP clone() const { return B2OP( new B2( *this ) ); }

private:
	int some_int_;
	A a_instance_;

};

}
}

#endif

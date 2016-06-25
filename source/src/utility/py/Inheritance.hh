// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   utility/py/Inheritance.hh
/// @author Sergey Lyskov
///
/// @note Defining various help classes to test C++ SubClassing in PyRosetta

#ifndef INCLUDED_utility_py_Inheritance_hh
#define INCLUDED_utility_py_Inheritance_hh

#include <utility/pointer/owning_ptr.hh>
#include <utility/exit.hh>

#include <string>
#include <iostream>

namespace utility { namespace py {

/** @details Dummy Value class to test if we can pass argument by-reference to PyRosetta sub-classes **/
class Value
{
public:
	Value(std::string s="unknown") { s_=s; }
	Value( Value const & src ) { std::cout << "Value Copy constructor! ["<< src.s_ << "]" << std::endl; utility_exit_with_message("Got Value copy-constructor call.. this really should not be happening... most likely is that PyRosetta sub-classing got broken!"); }

	~Value() { std::cout << "Value destructor!  ["<< s_ << "]" << std::endl; }

	std::string s() const { return s_; }
	void s(std::string const & s) { s_=s; }

private:
	std::string s_;
};

inline std::ostream & operator << ( std::ostream & os, Value const & v) { os << "Value('" << v.s() << "') "; return os; }

typedef utility::pointer::shared_ptr< Value > ValueOP;
typedef utility::pointer::shared_ptr< Value const > ValueCOP;


/** @details Dummy Base class to test if inheritance work properly in PyRosetta **/
class Base
{
public:
	Base() : int_value_(0) { std::cout << "Base::Base()" << std::endl; }

	int int_value_;
	std::string string_value_;

	virtual void foo_int(int i)            {  std::cout << "Base::foo_int( " << i << " )" << std::endl; }
	virtual void foo_string(std::string s) {  std::cout << "Base::foo_string( '" << s << "' )" << std::endl; }

	virtual void foo_value(Value& v    )   {  std::cout << "Base::foo_value( " << v << " )" << std::endl; }
	virtual void foo_value_sp(ValueOP v)   {  std::cout << "Base::foo_value_sp( " << *v << " )" << std::endl; }
};

class OverloadTest
{
public:
	Value value_;

	void self_test_virtual_ref() {
		std::cout << "Calling test(Value &v)..." << std::endl;
		test_ref(value_);
		std::cout << "Calling test(Value &v)... Done!" << std::endl;
	}

	void self_pure_test_virtual_ref() {
		std::cout << "Calling pure_test(Value &v)..." << std::endl;
		pure_test_ref(value_);
		std::cout << "Calling pure_test(Value &v)... Done!" << std::endl;
	}

	void self_test_virtual_p() {
		std::cout << "Calling test(Value *v)..." << std::endl;
		test_p(&value_);
		std::cout << "Calling test(Value *v)... Done!" << std::endl;
	}

	void self_pure_test_virtual_p() {
		std::cout << "Calling pure_test(Value *v)..." << std::endl;
		pure_test_p(&value_);
		std::cout << "Calling pure_test(Value *v)... Done!" << std::endl;
	}

	virtual void test_p(Value *v) {}
	virtual void pure_test_p(Value *v) = 0;

	virtual void test_ref(Value &v) {}
	virtual void pure_test_ref(Value &v) = 0;
};


typedef utility::pointer::shared_ptr< Base > BaseOP;
typedef utility::pointer::shared_ptr< Base const > BaseCOP;


inline void foo_all(Base &b, int i, Value &v, std::string s)     { b.foo_int(i);   b.foo_value(v);   b.foo_string(s); }
inline void foo_all_sp(BaseOP b, int i, Value &v, std::string s) { b->foo_int(i);  b->foo_value(v);  b->foo_string(s); }


} } // namespace py namespace utility

#endif // INCLUDED_utility_py_Inheritance_hh

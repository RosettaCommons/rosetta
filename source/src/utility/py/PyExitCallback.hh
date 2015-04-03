// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
// :noTabs=false:tabSize=4:indentSize=4:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   src/utility/PyExitCallback.hh
/// @brief  Class to specify utility_exit callback in PyRosetta.
/// @author Sergey Lyskov


#ifndef INCLUDED_utility_PyExitCallback_hh
#define INCLUDED_utility_PyExitCallback_hh

#include <utility/py/PyExitCallback.fwd.hh>

#include <utility/pointer/ReferenceCount.hh>

//#include <boost/utility.hpp>

#include <string>

namespace utility { namespace py {


/// This class for holding callback function.
class PyExitCallback : public utility::pointer::ReferenceCount //, boost::noncopyable
{
public:
	virtual ~PyExitCallback() {}

	virtual void exit_callback();

	static void set_PyExitCallBack(PyExitCallbackOP exit_callback_object);


	// Bindings Testing/Demo functions to test Windows DLL static inits
	static std::string       get_static_string() { return static_string_; }
	static std::string const get_static_string_const() { return static_string_const_; }
	static int               get_static_int() { return static_int_; }

private:
	static std::string       static_string_;
	static std::string const static_string_const_;
	static int               static_int_;

private:
	static void global_exit_callback(void);

	static PyExitCallbackOP & current_callback_object();
};


} } // namespace py namespace utility

#endif // INCLUDED_utility_PyExitCallback_hh

// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
// :noTabs=false:tabSize=4:indentSize=4:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   src/utility/PyExitCallback.cc
/// @brief  Class to specify utility_exit callback in PyRosetta.
/// @author Sergey Lyskov

#include <utility/PyExitCallback.hh>

#include <utility/exit.hh>

//#include <iostream>

namespace utility {

PyExitCallbackOP & PyExitCallback::current_callback_object()
{
	static PyExitCallbackOP callback;
	return callback;
}

std::string       PyExitCallback::static_string_ = "Just static string...";
std::string const PyExitCallback::static_string_const_ = "Just static string const...";
int               PyExitCallback::static_int_ = 42;

void PyExitCallback::exit_callback()
{
 	//std::cout << "PyExitCallback::exit_callback() default handler!" << std::endl;
	throw "PyRosetta Exception!";
}

void PyExitCallback::set_PyExitCallBack(PyExitCallbackOP exit_callback_object)
{
	current_callback_object() = exit_callback_object;
	set_main_exit_callback(global_exit_callback);
}

void PyExitCallback::global_exit_callback(void)
{
	//std::cout << "PyExitCallback::global_exit_callback" << std::endl;
	current_callback_object()->exit_callback();
}


} // namespace utility

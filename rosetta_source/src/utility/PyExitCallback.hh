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
///

#ifndef INCLUDED_utility_PyExitCallback_hh
#define INCLUDED_utility_PyExitCallback_hh

#include <utility/PyExitCallback.fwd.hh>

#include <utility/pointer/ReferenceCount.hh>

namespace utility {

/// This class for holding callback function.
class PyExitCallback : public utility::pointer::ReferenceCount
{
public:
	virtual ~PyExitCallback() {}

	virtual void exit_callback();

	static void set_PyExitCallBack(PyExitCallbackOP exit_callback_object);

private:
	static void global_exit_callback(void);

	static PyExitCallbackOP current_callback_object_;
};


} // namespace utility

#endif // INCLUDED_utility_PyExitCallback_hh

// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   utility/py/PyAssert.hh
/// @author Sergey Lyskov
///
/// @note Defining special assert version for PyRosetta build. It will allow to perform some additional
///       verification for PyRosetta build.


#ifndef INCLUDED_utility_py_PyAssert_hh
#define INCLUDED_utility_py_PyAssert_hh

#include <utility/exit.hh>

/// @brief PyAssert macro. Check if statement is true and call utility_exit otherwise, which in turn
///        will throw boost::python exception that python interpreter will catch.
///        This macro works only when PYROSETTA macro was defined.
///
/// Example of usage: PyAssert( (residue_num > 0) && (residue_num<=pose.n_residue()), "Residue not found in pose!"  )
///
#if defined(PYROSETTA) || ! defined(NDEBUG)  // PyRosetta asserts, we also enable then in our DEBUG mode as well
#define PyAssert(x, m) { if(!(x)) utility_exit_with_message(m); }
#else
	#define PyAssert(x, m)
#endif

namespace utility {
namespace py {
} // namespace py
} // namespace utility

#endif // INCLUDED_utility_py_PyAssert_HH

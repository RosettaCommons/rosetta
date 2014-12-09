// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   utility/Show.hh
///
/// @brief  uniformed printing functionality for data types.
/// @author Sergey Lyskov (Sergey.Lyskov@jhu.edu)

#ifndef INCLUDED_utility_Show_hh
#define INCLUDED_utility_Show_hh

// C++ headers
#include <istream>
#include <vector>

namespace utility {

//
// Base class for any showable object. To add ability to print to your object - just inherit
// this class, and overload 'show' method.
//
class Show {
public:
	virtual ~Show() {};

	virtual void show(std::ostream &) const {};
};


} // namespace utility


#endif // INCLUDED_utility_Show_HH

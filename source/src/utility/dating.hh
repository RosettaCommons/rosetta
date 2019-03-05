// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file    utility/dating.hh
/// @brief   Declarations for utility functions involving calendar dates, not for
///          finding new romantic partners.
/// @author  Labonte <JWLabonte@jhu.edu>
/// @note    This file was created at the request of Vikram.  Currently, the only
///          code using this (currently) single function is .pdb file output
///          code, but perhaps other needs for dating output might be needed in
///          the future....


#ifndef INCLUDED_utility_dating_HH
#define INCLUDED_utility_dating_HH

// C++ Headers
#include <string>


namespace utility {

/// @brief  Enumeration of acceptable date formats.
enum DateFormat {
	PDB_FORMAT = 1  // dd-MMM-yy
};


/// @brief  Return current date in the requested format.
std::string get_current_date( DateFormat const format=PDB_FORMAT );

}  // namespace utility

#endif  // INCLUDED_utility_dating_HH

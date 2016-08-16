// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/fragment/picking_old/vall/vall_io.hh
/// @brief  reading/writing of Vall libraries
/// @author Yih-En Andrew Ban (yab@u.washington.edu)

#ifndef INCLUDED_core_fragment_picking_old_vall_vall_io_hh
#define INCLUDED_core_fragment_picking_old_vall_vall_io_hh

// type headers
#include <core/types.hh>

// package headers
#include <core/fragment/picking_old/vall/VallLibrary.fwd.hh>

// C++ headers
#include <string>


namespace core {
namespace fragment {
namespace picking_old {
namespace vall {


/// @brief load standard Vall library from file
/// @param[in] filename
/// @param[out] vall
void vall_library_from_file( std::string const & filename, VallLibrary & library, core::Size const preallocate = 0 );


} // vall
} // picking_old
} // fragment
} // core


#endif /* INCLUDED_core_fragment_picking_old_vall_vall_io_HH */

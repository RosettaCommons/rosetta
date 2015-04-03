// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file protocols/forge/methods/fragment_functions.hh
/// @brief methods for manipulating fragments
/// @author Yih-En Andrew Ban (yab@u.washington.edu)

#ifndef INCLUDED_protocols_forge_methods_fragment_functions_hh
#define INCLUDED_protocols_forge_methods_fragment_functions_hh

// type headers
#include <core/types.hh>

// project headers
#include <core/fragment/ConstantLengthFragSet.fwd.hh>

#include <core/fragment/FrameIterator.fwd.hh>


namespace protocols {
namespace forge {
namespace methods {


/// @brief create small-mers from large-mers
/// @param[in] all_possible_smallmers Default false.  If true, grab all
///  possible smallmers using every possible starting position in a largemer
///  (you could end up with a *lot* of fragments per position).  If false,
///  grab smallmers using only the first position of a largemer as the starting
///  position.
core::fragment::ConstantLengthFragSetOP
smallmer_from_largemer(
	core::fragment::ConstFrameIterator begin,
	core::fragment::ConstFrameIterator end,
	core::Size const smallmer_size = 1,
	bool const all_possible_smallmers = false
);


} // methods
} // forge
} // protocols


#endif /* INCLUDED_protocols_forge_methods_fragment_functions_HH */

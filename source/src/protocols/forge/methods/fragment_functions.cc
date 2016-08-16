// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/forge/methods/fragment_functions.cc
/// @brief methods for manipulating fragments
/// @author Yih-En Andrew Ban (yab@u.washington.edu)

// unit headers
#include <protocols/forge/methods/fragment_functions.hh>

// project headers
#include <core/fragment/ConstantLengthFragSet.hh>
#include <core/fragment/FrameIteratorWorker_.hh>

#include <core/fragment/FrameIterator.hh>
#include <utility/vector1.hh>

//Auto Headers
#ifdef WIN32
#include <core/fragment/FragID.hh>
#endif


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
	core::Size const smallmer_size,
	bool const all_possible_smallmers
)
{
	using namespace core;
	using namespace core::fragment;

	assert( smallmer_size > 0 );

	ConstantLengthFragSetOP small( new ConstantLengthFragSet( smallmer_size ) );

	for ( ConstFrameIterator f = begin; f != end; ++f ) {
		Size const ie = all_possible_smallmers ? f->length() - smallmer_size + 1 : 1;
		for ( Size i = 1; i <= ie; ++i ) {
			small->add( f->generate_sub_frame( smallmer_size, i ) );
		}
	}

	return small;
}


} // methods
} // forge
} // protocols

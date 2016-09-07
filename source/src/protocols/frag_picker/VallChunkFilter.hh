// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   protocols/frag_picker/VallChunkFilter.hh
/// @brief  says whether a given chunk is interesting or not
/// @author Dominik Gront (dgront@chem.uw.edu.pl)


#ifndef INCLUDED_protocols_frag_picker_VallChunkFilter_hh
#define INCLUDED_protocols_frag_picker_VallChunkFilter_hh

// unit headers
#include <protocols/frag_picker/VallChunkFilter.fwd.hh>

// package headers
#include <protocols/frag_picker/VallChunk.fwd.hh>

// utility headers
#include <utility/pointer/ReferenceCount.hh>

namespace protocols {
namespace frag_picker {

/// @brief  a base class for a chunk filtering mechanism
/// @details Chunk filtering is used to screen chunks before any fragment is evaluated
/// Therefore it is the fastest way to excluded unwanted proteins
/// @see AllowPdbIdFilter and DenyPdbIdFilter for implementations
class VallChunkFilter: public utility::pointer::ReferenceCount {
public:
	~VallChunkFilter() override = default;
	/// @brief if a chunk seems to be interesting, return true. Otherwise say false
	virtual bool test_chunk(VallChunkOP) = 0;
};

} // frag_picker
} // protocols


#endif /* INCLUDED_protocols_frag_picker_VallChunkFilter_HH */

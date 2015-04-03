// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/fragments/FrameList.hh
/// @brief  set of fragments for a certain alignment frame
/// @author Oliver Lange (olange@u.washington.edu)
/// @date   Wed Oct 20 12:08:31 2007
///
#ifndef INCLUDED_core_fragment_FrameList_HH
#define INCLUDED_core_fragment_FrameList_HH

// Unit Headers
#include <core/fragment/FrameList.fwd.hh>

// Package Headers
#include <core/fragment/Frame.hh>

// Utility headers
#include <utility/vector1.fwd.hh>

// C++ STL Headers

#include <core/fragment/FragID.fwd.hh>
#include <core/types.hh>

#include <utility/vector1_bool.hh>

namespace core {
namespace fragment {

class FrameList : public utility::vector1< FrameOP > {
public:
	// access physical fragments in the FrameList by their "flat" nr,
	// i.e., if the FrameList would be iterated with a FragID_Iterator
	// this allows uniform sampling over all fragments in the FrameList.
	FragID fragID ( Size flat_nr );
	core::Size flat_size() const;

	utility::vector1<FrameOP> frame_vector();
};

extern std::ostream& operator<< ( std::ostream& out, FrameList const& );
}
}

#endif

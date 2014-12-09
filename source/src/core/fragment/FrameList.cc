// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/fragments/FragData.cc
/// @brief  a collection classes of the FragData and SingleResidueFragData class hirachy
/// @author Oliver Lange (olange@u.washington.edu)
/// @date   Wed Oct 20 12:08:31 2007
///

// Unit Headers
#include <core/fragment/FrameList.hh>

// Package Headers
// AUTO-REMOVED #include <core/fragment/FragData.hh>
#include <core/fragment/Frame.hh>

// Project Headers
// AUTO-REMOVED #include <core/pose/Pose.hh>
// AUTO-REMOVED #include <core/conformation/ResidueFactory.hh>

// ObjexxFCL Headers

// Utility headers
#include <basic/Tracer.hh>

#include <core/fragment/FragID.hh>
#include <utility/exit.hh>
#include <utility/vector1.hh>


namespace core {
namespace fragment {

static thread_local basic::Tracer tr( "core.fragment" );

FragID FrameList::fragID ( Size flat_nr ) {
  Size passed_frags( 0 );
  iterator it = begin(), eit=end();

	while (it!=eit && (passed_frags + (*it)->nr_frags()) < flat_nr)
	{
    passed_frags += (*it)->nr_frags();
		++it;
	}

  if ( it!=eit )
	{
    return FragID( *it, flat_nr - passed_frags );
  }
	else
	{
		runtime_assert( 0 ); //out of bounds
		return FragID();
	}
}

Size FrameList::flat_size() const {
  Size frags( 0 );
  for (const_iterator it=begin(), eit=end(); it!=eit ; ++it) {
    frags += (*it)->nr_frags();
  }
  return frags;
}

utility::vector1<FrameOP> FrameList::frame_vector()
{
	return utility::vector1<FrameOP>(*this);
}

std::ostream& operator<< ( std::ostream& out, FrameList const& frags) {
  for ( FrameList::const_iterator it=frags.begin(), eit=frags.end(); it!=eit; ++it ) {
    (*it)->show( out );
  }
	return out;
}

} //fragment
} //core

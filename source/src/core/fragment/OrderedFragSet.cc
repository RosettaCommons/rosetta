// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/fragments/FragSet.cc
/// @brief  set of fragments for a certain alignment frame
/// @author Oliver Lange (olange@u.washington.edu)
/// @author James Thompson (tex@u.washington.edu)
/// @date   Wed Oct 20 12:08:31 2007


// Unit Headers
#include <core/fragment/OrderedFragSet.hh>

// Package Headers
//#include <core/fragment/BBTorsionSRFD.hh>
#include <core/fragment/OrderedFragSetIterator_.hh>
#include <core/fragment/Frame.hh>

// Project Headers
#include <core/types.hh>


// ObjexxFCL Headers

// Utility headers
#include <utility/vector1.fwd.hh>
#include <utility/pointer/owning_ptr.hh>
#include <basic/Tracer.hh>
#include <ostream>

#include <core/fragment/FrameIterator.hh>
#include <utility/vector1.hh>

#ifdef WIN32
#include <iterator>
#endif

namespace core {
namespace fragment {

using namespace kinematics;

static THREAD_LOCAL basic::Tracer tr( "core.fragments" );
// preliminary reader method --- reads classic rosetta++ frag files


OrderedFragSet::OrderedFragSet() {}
OrderedFragSet::~OrderedFragSet() = default;

FragSetOP OrderedFragSet::clone() const {
	return FragSetOP( new OrderedFragSet( *this ) );
}
FragSetOP OrderedFragSet::empty_clone() const {
	return FragSetOP( new OrderedFragSet() );
}

/// @brief get fragments that start somewhere between start and end
Size OrderedFragSet::region(
	MoveMap const&,
	core::Size start,
	core::Size end, //not used
	core::Size, //min_overlap not used
	core::Size, //min_length not used
	FrameList &frame_list
) const {
	Size count( 0 );
	for ( Size pos=start; pos<=end; pos++ ) {
		count += frames( pos, frame_list );
	}
	return count;
}


/// @brief Accessor for the Frame at the specified insertion position. Returns false if
/// there is no frame at the specified position.
Size OrderedFragSet::frames( Size pos, FrameList &out_frames ) const
{
	auto it = frames_.find(pos);
	if ( it == frames_.end() ) return 0;
	if ( it->second.begin() != it->second.end() ) {
		copy( it->second.begin(), it->second.end(), back_inserter( out_frames ) ); // should append frames
		return it->second.size();
	}

	return 0;
}

ConstFrameIterator OrderedFragSet::begin() const {
	return ConstFrameIterator( FrameIteratorWorker_OP( new OrderedFragSetIterator_( frames_.begin(), frames_.end() ) ) );
}

ConstFrameIterator OrderedFragSet::end() const {
	return ConstFrameIterator( FrameIteratorWorker_OP( new OrderedFragSetIterator_( frames_.end(), frames_.end() ) ) );
}

FrameIterator OrderedFragSet::nonconst_begin() {
	return FrameIterator( FrameIteratorWorker_OP( new OrderedFragSetIterator_( frames_.begin(), frames_.end() ) ) );
}

FrameIterator OrderedFragSet::nonconst_end() {
	return FrameIterator( FrameIteratorWorker_OP( new OrderedFragSetIterator_( frames_.end(), frames_.end() ) ) );
}

bool OrderedFragSet::empty() const {
	return frames_.empty();
}


void OrderedFragSet::add_( FrameOP aframe )
{
	Size seqpos( aframe->start() );
	frames_[ seqpos ].push_back( aframe );
}

}//fragment
}// core

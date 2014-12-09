// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/fragments/FragSet.hh
/// @brief  set of fragments
/// @author Oliver Lange ( olange@u.washington.edu)
/// @date   Wed Aug 22 12:08:31 2007
///

#ifndef INCLUDED_core_fragment_FrameIteratorWorker__HH
#define INCLUDED_core_fragment_FrameIteratorWorker__HH

// Unit Headers
#include <core/fragment/FrameIteratorWorker_.fwd.hh>

// type headers
// AUTO-REMOVED #include <core/types.hh>

// Package Headers
#include <core/fragment/Frame.hh>
// AUTO-REMOVED #include <core/fragment/FragID.hh>

// Project Headers
// AUTO-REMOVED #include <core/pose/Pose.fwd.hh>

// AUTO-REMOVED #include <core/fragment/FrameIterator.fwd.hh>
// AUTO-REMOVED #include <core/fragment/FragID_Iterator.fwd.hh>


// Utility headers
// AUTO-REMOVED #include <utility/vector1.hh>
#include <utility/pointer/ReferenceCount.hh>

// std Headers
// AUTO-REMOVED #include <iterator>

#include <core/fragment/Frame.fwd.hh>
#include <iostream>


/* Just a mad thought: with fragments becoming ever more "Residue" like one might want to use the
	 packer to choose a combination of good fragments instead of makeing independent choices.
	 I guess, it is only a question of keeping the combinatorics in control...
	 maybe it makes sense to pack with only "unconfident" regions of the backbone flexible ..
*/

namespace core {
namespace fragment {


// virtual Base class for FrameIterator
class FrameIteratorWorker_ : public utility::pointer::ReferenceCount {
	friend class FrameIterator;
	friend class ConstFrameIterator;
	friend class FragID_Iterator;
protected:
	//	FrameIteratorWorker_ ( BaseFragSet& bfs );
	virtual bool operator != ( FrameIteratorWorker_ const& ) const = 0;
	virtual bool operator == ( FrameIteratorWorker_ const& fit) const {
		return !( (*this)!= fit );
	};

	virtual FrameIteratorWorker_& operator++ ()=0;

	virtual FrameIteratorWorker_& operator+ ( Size offset ) {
		for ( Size i = 1; i<=offset; i++ ) ++(*this);
		return *this;
	}

	virtual FrameIteratorWorker_& operator = ( FrameIteratorWorker_ const& itr ) = 0;

	FrameOP operator* () {
		return frame_ptr();
	}

	FrameCOP operator* () const {
		return frame_ptr();
	}

	FrameOP operator-> () {
		return ( frame_ptr() );
	}

	FrameCOP operator-> () const {
		return ( frame_ptr() );
	}

	virtual FrameOP frame_ptr() = 0;
	virtual FrameCOP frame_ptr() const = 0;
};

}
}
#endif

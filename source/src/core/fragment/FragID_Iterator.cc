// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/fragments/FragID_Iterator.cc
/// @brief
/// @author Oliver Lange ( olange@u.washington.edu)
/// @date

// Unit Headers
#include <core/fragment/FragID_Iterator.hh>

// Package Headers
#include <core/fragment/FragID.hh>
#include <core/fragment/Frame.hh>
#include <core/fragment/FrameIteratorWorker_.hh>
#include <core/fragment/FrameList.hh>
#include <core/fragment/FrameListIterator_.hh>

// Project Headers
// AUTO-REMOVED #include <core/kinematics/MoveMap.hh>
#include <core/types.hh>

// ObjexxFCL Headers

// Utility headers
// AUTO-REMOVED #include <utility/vector1.hh>
#include <utility/pointer/ReferenceCount.hh>

// std Headers
// AUTO-REMOVED #include <iterator>

#include <core/fragment/FrameIterator.hh>
#include <utility/vector1.hh>


/* Just a mad thought: with fragments becoming ever more "Residue" like one might want to use the
	 packer to choose a combination of good fragments instead of makeing independent choices.
	 I guess, it is only a question of keeping the combinatorics in control...
	 maybe it makes sense to pack with only "unconfident" regions of the backbone flexible ..
*/

namespace core {
namespace fragment {

FragID_Iterator::FragID_Iterator( FrameIterator it ) : it_( it.it_ ), ipos_(1) {}
FragID_Iterator::FragID_Iterator( ConstFrameIterator it ) : it_( it.it_ ), ipos_(1) {}
FragID_Iterator::FragID_Iterator( FrameIteratorWorker_OP it ) : it_( it ), ipos_(1) {}
FragID_Iterator::FragID_Iterator( FrameList::iterator it ) : it_( FrameIteratorWorker_OP( new FrameListIterator_( it ) ) ), ipos_( 1 ) {}

FragID_Iterator::~FragID_Iterator() {}

FragID_Iterator::FragID_Iterator() : it_( /* NULL */ ) {}

bool FragID_Iterator::operator != ( FragID_Iterator const& fi) const {
	return (*it_) != (*fi.it_) || ipos_ != fi.ipos_;
}

bool FragID_Iterator::operator == ( FragID_Iterator const& fi) const {
	return !operator!=( fi);
}

FragID_Iterator& FragID_Iterator::operator++ () {
	//		std::cout << "it - nr_frags " << (*it_)->nr_frags() << std::endl;
	//		std::cout << "ipos_ " << ipos_ << std::endl;
	if ( (*it_)->nr_frags() > ipos_ ) {
		++ipos_;
	} else {
		// if we assume that NEVER an empty frame is in the fragset we can avoid the use of eit
		// this assumption is TRUE now: we require the FrameIterator to always show to a valid frame
		//			while ( ++it != eit ) {
		//				if (it_->nr_frags() ) break;
		//			}
		//			ipos_ = 1;
		++(*it_);
		ipos_ = 1;
	}
	return *this;
}

FragID_Iterator& FragID_Iterator::operator+ ( Size offset ) {
	for ( Size i = 1; i<=offset ; i++ )	operator++();
	return *this;
}

FragID_Iterator& FragID_Iterator::operator = ( FragID_Iterator const& itr ) {
	it_=itr.it_; //copy the pointers to the real iterators
	ipos_= itr.ipos_;
	return *this;
}

FragID FragID_Iterator::frag_id() {
	return FragID( it_->frame_ptr() , (*it_)->frag_id( ipos_) );
}

FragID FragID_Iterator::operator* () {
	return frag_id();
}

// can't provide that operator, can I ? Could point to a member instance of FragID
// normally this operator is never used to actually asked for the pointer to something...
FragID* FragID_Iterator::operator-> () {
		my_frag_id_ = frag_id();
		return &my_frag_id_;
}


}
}

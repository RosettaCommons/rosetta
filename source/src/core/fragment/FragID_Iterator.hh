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


#ifndef INCLUDED_core_fragment_FragID_Iterator_HH
#define INCLUDED_core_fragment_FragID_Iterator_HH

// Unit Headers

// Package Headers
#include <core/fragment/FragID.fwd.hh>
#include <core/fragment/FrameIterator.fwd.hh>
#include <core/fragment/FrameIteratorWorker_.fwd.hh>
#include <core/fragment/FrameList.hh>

// Project Headers
#include <core/types.hh>

// Utility headers

// std Headers

#include <core/fragment/FragID.hh>
#include <utility/vector1.hh>


/* Just a mad thought: with fragments becoming ever more "Residue" like one might want to use the
	 packer to choose a combination of good fragments instead of makeing independent choices.
	 I guess, it is only a question of keeping the combinatorics in control...
	 maybe it makes sense to pack with only "unconfident" regions of the backbone flexible ..
*/

namespace core {
namespace fragment {

class FragID_Iterator : std::iterator< std::forward_iterator_tag, FragID > {
public:
	FragID_Iterator( ConstFrameIterator it );
	FragID_Iterator( FrameIterator it );
	FragID_Iterator( FrameIteratorWorker_OP it );
	FragID_Iterator( FrameList::iterator it );
	~FragID_Iterator();

	FragID_Iterator();

	bool operator != ( FragID_Iterator const& fi) const;

	bool operator == ( FragID_Iterator const& fi) const;

	FragID_Iterator& operator++ ();

	FragID_Iterator& operator+ ( Size offset );

	FragID_Iterator& operator = ( FragID_Iterator const& itr );

	FragID frag_id();

	FragID operator* ();

	// can't provide that operator, can I ? Could point to a member instance of FragID
	// normally this operator is never used to actually asked for the pointer to something...
	FragID* operator-> ();

protected:
	//BaseFragSet& fragset_;
	FrameIteratorWorker_OP it_;
	Size ipos_; //intra_frame_pos
	FragID my_frag_id_;

};

}
}

#endif

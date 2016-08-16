// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/fragments/ConstantLengthFragSet.hh
/// @brief  yields a simple implementation of a fragset
/// @author Oliver Lange ( olange@u.washington.edu)
/// @date   Wed Aug 22 12:08:31 2007


#ifndef INCLUDED_core_fragment_ConstantLengthFragSetIterator__HH
#define INCLUDED_core_fragment_ConstantLengthFragSetIterator__HH

// Unit Headers
//#include <core/fragment/ConstantLengthFragSet.fwd.hh>

// Package Headers
#include <core/fragment/FrameList.hh>
#include <core/fragment/FrameIteratorWorker_.hh>


// Project Headers
#include <core/types.hh>

// Utility headers
#include <utility/vector1.fwd.hh>
#include <utility/pointer/ReferenceCount.hh>

// std Headers

#include <core/fragment/Frame.hh>
#include <utility/vector1.hh>


/* Just a mad thought: with fragments becoming ever more "Residue" like one might want to use the
packer to choose a combination of good fragments instead of makeing independent choices.
I guess, it is only a question of keeping the combinatorics in control...
maybe it makes sense to pack with only "unconfident" regions of the backbone flexible ..
*/

namespace core {
namespace fragment {


class ConstantLengthFragSetIterator_ : public FrameIteratorWorker_ {
	friend class ConstantLengthFragSet;

	ConstantLengthFragSetIterator_ & operator = (ConstantLengthFragSetIterator_ const&);
protected:
	ConstantLengthFragSetIterator_( FrameList::const_iterator it, FrameList::const_iterator eit ) : it_( it ), eit_( eit ) {
		if ( it != eit ) {
			if ( *it ) {
				if ( (*it)->nr_frags() ) return;
			}
			++(*this); //if not already pointing to some valid stuff ... increment until it is
		}
	};

	bool operator != ( FrameIteratorWorker_ const& fiw ) const {
		ConstantLengthFragSetIterator_ const& fsit ( dynamic_cast< ConstantLengthFragSetIterator_ const& > ( fiw ) );
		return it_!=fsit.it_;
	};

	FrameIteratorWorker_& operator++ () {
		while ( ++it_ != eit_ ) {
			if ( *it_ ) {
				if ( (*it_)->nr_frags() ) return *this;
			}
		};
		return *this;
	}

	FrameIteratorWorker_& operator = ( FrameIteratorWorker_ const& fiw ) {
		ConstantLengthFragSetIterator_ const& fsit ( dynamic_cast< ConstantLengthFragSetIterator_ const& > ( fiw ) );
		it_= fsit.it_;
		return *this;
	}

	FrameOP frame_ptr() {
		return *it_;
	}

	FrameCOP frame_ptr() const {
		return *it_;
	}

private:
	FrameList::const_iterator it_;
	FrameList::const_iterator eit_;
};

}
}

#endif

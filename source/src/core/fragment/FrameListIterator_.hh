// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/fragments/ConstantLengthFragSet.hh
/// @brief  yields a simple implementation of a fragset
/// @author Oliver Lange ( olange@u.washington.edu)
/// @date   Wed Aug 22 12:08:31 2007


#ifndef INCLUDED_core_fragment_FramelistIterator__HH
#define INCLUDED_core_fragment_FramelistIterator__HH

// Package Headers
#include <core/fragment/FrameList.hh>

// Project Headers

// Package Headers

#include <core/fragment/FrameIteratorWorker_.hh>

/* Just a mad thought: with fragments becoming ever more "Residue" like one might want to use the
packer to choose a combination of good fragments instead of makeing independent choices.
I guess, it is only a question of keeping the combinatorics in control...
maybe it makes sense to pack with only "unconfident" regions of the backbone flexible ..
*/

namespace core {
namespace fragment {


class FrameListIterator_ : public FrameIteratorWorker_ {
	friend class FrameList;
	friend class FragID_Iterator;
protected:
	FrameListIterator_( FrameList::iterator it ) : it_( it ) {};

	bool operator != ( FrameIteratorWorker_ const& fiw ) const {
		FrameListIterator_ const& fsit ( dynamic_cast< FrameListIterator_ const& > ( fiw ) );
		return it_!=fsit.it_;
	};

	FrameIteratorWorker_& operator++ () {
		++it_;
		return *this;
	}

	FrameIteratorWorker_& operator = ( FrameIteratorWorker_ const& fiw ) {
		FrameListIterator_ const& fsit ( dynamic_cast< FrameListIterator_ const& > ( fiw ) );
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
	FrameList::iterator it_;
};

}
}

#endif

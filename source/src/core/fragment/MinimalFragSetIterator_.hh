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
/// @author David E Kim ( dekim@u.washington.edu)
/// @date   Wed Aug 22 12:08:31 2007


#ifndef INCLUDED_core_fragment_MinimalFragSetIterator__HH
#define INCLUDED_core_fragment_MinimalFragSetIterator__HH

// Unit Headers
//#include <core/fragment/ConstantLengthFragSet.fwd.hh>

// Package Headers
#include <core/fragment/FrameList.hh>
#include <core/fragment/FrameIteratorWorker_.hh>
//#include <core/fragment/MinimalFragSet.hh>

//#include <core/fragment/FragSet.hh>
//#include <core/fragment/Frame.hh>
//#include <core/fragment/Frame.fwd.hh>


// Project Headers
//#include <core/pose/Pose.hh>
//#include <core/kinematics/MoveMap.hh>
#include <core/types.hh>

// ObjexxFCL Headers

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

/* iterate over a map of vectors such that if feels like a flat data-structure */
class MinimalFragSetIterator_ : public FrameIteratorWorker_ {
	friend class MinimalFragSet;
	typedef std::map< Size, FrameList > FrameMap;
	typedef FrameMap::const_iterator OuterIterator;
	typedef FrameList::const_iterator InnerIterator;

	MinimalFragSetIterator_ & operator = (MinimalFragSetIterator_ const&);
protected:
	MinimalFragSetIterator_( OuterIterator it, OuterIterator eit ) : outer_( it ), outer_end_( eit ) {
		if ( outer_ != outer_end_ ) {
			inner_ = outer_->second.begin();
			inner_end_ = outer_->second.end();
			if ( !(inner_ != inner_end_) ) increment_outer();
			else if ( (*this)->nr_frags() == 0 ) ++(*this); //move forward to first real frame
		}
	};

	bool operator != ( FrameIteratorWorker_ const& fiw ) const override {
		MinimalFragSetIterator_ const& fsit ( dynamic_cast< MinimalFragSetIterator_ const& > ( fiw ) );
		bool bOut ( outer_!=fsit.outer_ );
		if ( !bOut && fsit.outer_ != fsit.outer_end_ && outer_ != outer_end_ ) {
			return inner_ != fsit.inner_;
		} else return bOut;
	};

	FrameIteratorWorker_& operator++ () override {
		if ( increment_inner() ) return *this; //increment inner-loop to the end
		increment_outer(); //then increment outer-loop
		return *this;
	}

	bool increment_outer () {
		while ( ++outer_ != outer_end_ ) {
			inner_ = outer_->second.begin();
			inner_end_ = outer_->second.end();
			if ( inner_ != inner_end_ ) {
				if ( (*inner_)->nr_frags() ) return true;
				else if ( increment_inner() ) return true;
			}
		}
		return false;
	}

	bool increment_inner () {
		while ( ++inner_ != inner_end_ ) {
			if ( (*inner_)->nr_frags() ) return true;
		}
		return false;
	}


	FrameIteratorWorker_& operator = ( FrameIteratorWorker_ const& fiw ) override {
		MinimalFragSetIterator_ const& fsit ( dynamic_cast< MinimalFragSetIterator_ const& > ( fiw ) );
		inner_ = fsit.inner_;
		inner_end_ = fsit.inner_end_;
		outer_ = fsit.outer_;
		outer_end_ = fsit.outer_end_;
		return *this;
	}

	FrameOP frame_ptr() override {
		return *inner_; //call get() of owning_ptr
	}

	FrameCOP frame_ptr() const override {
		return *inner_; //call get() of owning_ptr
	}

private:
	OuterIterator outer_;
	OuterIterator outer_end_;
	InnerIterator inner_;
	InnerIterator inner_end_;
};

}
}

#endif

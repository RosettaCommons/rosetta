// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
// :noTabs=false:tabSize=4:indentSize=4:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/fragments/FrameIterator.cc
/// @brief
/// @author Oliver Lange ( olange@u.washington.edu)
/// @date

// Unit Headers
#include <core/fragment/FrameIterator.hh>

// Package Headers
#include <core/fragment/FrameIteratorWorker_.hh>
// AUTO-REMOVED #include <core/fragment/FragID_Iterator.hh>

// Project Headers
#include <core/types.hh>

// ObjexxFCL Headers

// Utility headers
// AUTO-REMOVED #include <utility/vector1.hh>
#include <utility/pointer/ReferenceCount.hh>

#include <utility/vector1.hh>


// std Headers
// AUTO-REMOVED #include <iterator>

/*
Might change the FrameIterator to return FrameOP instead of Frame&...
This has the advantage of bein closer to the std:: implementations...
for instance one could do a
FrameList my_list;
FrameIterator it=fragset.begin(), eit=fragset.end();
copy( it, eit,back_inserter( my_list ));

the it-> wouldn't change, it still returns a Frame*
the *it would change and return a Frame* ...

*/

namespace core {
namespace fragment {

FrameIterator::FrameIterator( FrameIteratorWorker_OP it ) : it_( it ) {}
FrameIterator::FrameIterator() : it_( NULL ) {}
FrameIterator::~FrameIterator() {}

bool FrameIterator::operator != ( FrameIterator const& fi) const {
	return (*it_)!=(*fi.it_);
}

bool FrameIterator::operator == ( FrameIterator const& fi) const {
	return (*it_)==(*fi.it_);
}

FrameIterator& FrameIterator::operator++ () {
	++(*it_);
	return *this;
}

FrameIterator& FrameIterator::operator+ ( Size offset ) {
	(*it_)+( offset);
	return *this;
}

FrameIterator const & FrameIterator::operator = ( FrameIterator const& itr ) {
	it_=itr.it_; //copy the pointers to the real iterators
	return *this;
}

Frame* FrameIterator::operator* () {
	return frame_ptr();
}

Frame const* FrameIterator::operator* () const {
	return frame_ptr();
}

Frame* FrameIterator::operator-> () {
	return frame_ptr();
}

Frame const* FrameIterator::operator-> () const {
	return frame_ptr();
}

Frame* FrameIterator::frame_ptr() {
	return it_->frame_ptr();
}

Frame const* FrameIterator::frame_ptr() const{
	return it_->frame_ptr();
}



}
}

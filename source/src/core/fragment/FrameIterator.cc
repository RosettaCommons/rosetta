// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
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

ConstFrameIterator::ConstFrameIterator( FrameIteratorWorker_OP it ) : it_( it ) {}
ConstFrameIterator::ConstFrameIterator() : it_( /* NULL */ ) {}
ConstFrameIterator::~ConstFrameIterator() {}

bool ConstFrameIterator::operator != ( ConstFrameIterator const& fi) const {
	return (*it_)!=(*fi.it_);
}

bool ConstFrameIterator::operator == ( ConstFrameIterator const& fi) const {
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

ConstFrameIterator& ConstFrameIterator::operator++ () {
	++(*it_);
	return *this;
}

ConstFrameIterator& ConstFrameIterator::operator+ ( Size offset ) {
	(*it_)+( offset);
	return *this;
}


FrameIterator const & FrameIterator::operator = ( FrameIterator const& itr ) {
	it_=itr.it_; //copy the pointers to the real iterators
	return *this;
}

ConstFrameIterator const & ConstFrameIterator::operator = ( ConstFrameIterator const& itr ) {
	it_=itr.it_; //copy the pointers to the real iterators
	return *this;
}

FrameOP FrameIterator::operator* () {
	return frame_ptr();
}

FrameCOP ConstFrameIterator::operator* () const {
	return frame_ptr();
}

FrameOP FrameIterator::operator-> () {
	return frame_ptr();
}

FrameCOP ConstFrameIterator::operator-> () const {
	return frame_ptr();
}

FrameOP FrameIterator::frame_ptr() {
	return it_->frame_ptr();
}

FrameCOP ConstFrameIterator::frame_ptr() const{
	return it_->frame_ptr();
}



}
}

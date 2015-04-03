// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/fragments/FragSet.cc
/// @brief  set of fragments for a certain alignment frame
/// @author Oliver Lange (olange@u.washington.edu)
/// @author James Thompson (tex@u.washington.edu)
/// @date   Wed Oct 20 12:08:31 2007


// Unit Headers
#include <core/fragment/FragSetCollection.hh>
#ifdef WIN32
#include <core/fragment/Frame.hh>
#endif

// Package Headers
#include <core/fragment/FragSet.hh>
#include <core/fragment/FrameIteratorWorker_.hh>

// Project Headers
#include <core/types.hh>

// Utility headers
#include <utility/vector1.fwd.hh>
#include <utility/exit.hh>

//#include <basic/Tracer.hh>


#include <core/fragment/FrameIterator.hh>
#include <utility/vector1.hh>


namespace core {
namespace fragment {

//static basic::Tracer tr("core.fragment");
using namespace kinematics;

FragSetCollection::FragSetCollection() {}
FragSetCollection::~FragSetCollection() {}
FragSetCollection::FragSetCollection( FragSetCollection const & src ) :
	Parent( src ),
	fragset_list_( src.fragset_list_ )
{}

FragSetOP FragSetCollection::clone() const {
	return FragSetOP( new FragSetCollection( *this ) );
}
FragSetOP FragSetCollection::empty_clone() const {
	return FragSetOP( new FragSetCollection() );
}


Size
FragSetCollection::region(
	MoveMap const& move_map,
	Size start,
	Size end,
	Size min_overlap,
	Size min_length,
	FrameList &frames
) const {
  Size count ( 0 );
  for ( FragSetList::const_iterator it = fragset_list_.begin(), eit = fragset_list_.end();
        it != eit; ++it ) {
    count += (*it)->region(  move_map, start, end, min_overlap, min_length, frames );
  }
  return count;
}

void FragSetCollection::add_fragset( FragSetOP fragset ) {
  if ( fragset->max_frag_length() > max_frag_length() ) set_max_frag_length( fragset->max_frag_length() );
  if ( fragset->min_pos() < min_pos() ) set_min_pos( fragset->min_pos() );
  if ( fragset->max_pos() > max_pos() ) set_max_pos( fragset->max_pos() );
  fragset_list_.push_back( fragset );
}

ConstFrameIterator FragSetCollection::begin() const {
	std::cout << "iterator of FragSetCollection has stubbed out " << std::endl;
debug_assert( 0 );
	utility_exit_with_message( "iterator of FragSetCollection has stubbed out " );
	return ConstFrameIterator(); //too make compiler happy
}

ConstFrameIterator FragSetCollection::end() const {
	std::cout << "iterator of FragSetCollection has stubbed out " << std::endl;
debug_assert( 0 );
	utility_exit_with_message( "iterator of FragSetCollection has stubbed out " );
	return ConstFrameIterator(); //too make compiler happy
}

FrameIterator FragSetCollection::nonconst_begin() {
	std::cout << "iterator of FragSetCollection has stubbed out " << std::endl;
debug_assert( 0 );
	utility_exit_with_message( "iterator of FragSetCollection has stubbed out " );
	return FrameIterator(); //too make compiler happy
}

FrameIterator FragSetCollection::nonconst_end() {
	std::cout << "iterator of FragSetCollection has stubbed out " << std::endl;
debug_assert( 0 );
	utility_exit_with_message( "iterator of FragSetCollection has stubbed out " );
	return FrameIterator(); //too make compiler happy
}

bool FragSetCollection::empty() const {
	for ( FragSetList::const_iterator it = fragset_list_.begin(); it != fragset_list_.end(); ++it ) {
		if ( !(*it)->empty() ) return false;
	}
	return true;
}

void FragSetCollection::add_( FrameOP ) {
	//tricky which FragSet should it add the frame to? needs a way to determine this
	std::cout << "add Frame to FragSetCollection has stubbed out " << std::endl;
debug_assert( 0 );
	utility_exit_with_message( "add Frame to FragSetCollection has stubbed out " );
}


} //fragment
} //core

// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/fragments/FragID.cc
/// @brief  set of fragments for a certain alignment frame
/// @author Oliver Lange (olange@u.washington.edu)
/// @author James Thompson
/// @date   Wed Oct 20 12:08:31 2007

// Unit Headers
#include <core/fragment/FragID.hh>

// Package Headers
#include <core/fragment/Frame.hh>

#include <core/fragment/FragData.hh>
#include <utility/vector1.hh>


namespace core {
namespace fragment {


FragID::FragID() : first( /* 0 */ ), second( 0 ) {}

FragID::FragID( FrameCOP frame, Size frag_id )
	: first( frame ), second( frag_id )
{
	assert( frame->nr_frags() >= frag_id );
}

FragID::FragID( FragID const & src ) :
	first( src.first ),
	second( src.second )
{}

FragID::~FragID() {}

FragID const & FragID::operator = ( FragID const & rhs )
{
	if ( this != &rhs ) {
		first = rhs.first;
		second = rhs.second;
	}
	return *this;
}


FrameCOP
FragID::frame_ptr() const { return first; }

Frame const&
FragID::frame() const { return *first; }

//Frame&
//FragID::frame() { return *first; }

Size
FragID::id() const { return second; } // stores the nr you get by asking frame.frag_id ( nr );

FragData const &
FragID::fragment() const { return frame().fragment( id() ); }

FragDataCOP
FragID::fragment_ptr() const { return frame().fragment_ptr( id() ); }

Size
FragID::apply( kinematics::MoveMap const& mm, pose::Pose& pose) const
{
	return frame().apply( mm, id(), pose );
}

Size FragID::apply( pose::Pose & pose ) const
{
	return frame().apply( id(), pose );
}
// if we enable id != nr frame will need a map that can do frag_id --> nr
bool
FragID::is_valid() const
{
	return frame_ptr() && ( id() > 0 ) && fragment().is_valid();
}

Size
FragID::apply_ss( kinematics::MoveMap const& mm, std::string& ss ) const {
	return frame().apply_ss( mm, id(), ss );
}

bool FragID::operator == ( FragID const & other ) const
{
	return first == other.first && second == other.second;
}

bool FragID::operator <  ( FragID const & other ) const
{
	if ( first < other.first ) return true;
	else if ( first == other.first && second < other.second ) return true;

	return false;
}


} //fragment
} //core

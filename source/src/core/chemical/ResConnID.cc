// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file
/// @brief
/// @author Andrew Leaver-Fay


// Unit Headers
#include <core/chemical/ResConnID.hh>

namespace core {
namespace chemical {

ResConnID::ResConnID() : res_id_( 0 ), conn_id_( 0 ) {}

ResConnID::ResConnID( ResConnID const & rhs ) : res_id_( rhs.res_id_ ), conn_id_( rhs.conn_id_ ) {}

ResConnID::ResConnID( Size resid, Size connid ) : res_id_( resid ), conn_id_( connid ) {}

ResConnID & ResConnID::operator = ( ResConnID const & rhs )
{
	res_id_ = rhs.res_id_;
	conn_id_ = rhs.conn_id_;
	return *this;
}

bool operator < ( ResConnID const & lhs, ResConnID const & rhs )
{
	return ( lhs.res_id_ == rhs.res_id_ ? lhs.conn_id_ < rhs.conn_id_ : lhs.res_id_ < rhs.res_id_ );
}


bool operator == ( ResConnID const & lhs, ResConnID const & rhs )
{
	return ( lhs.res_id_ == rhs.res_id_ && lhs.conn_id_ == rhs.conn_id_ );
}

bool operator != ( ResConnID const & lhs, ResConnID const & rhs )
{
	return !( lhs == rhs );
}


Size
ResConnID::resid() const
{
	return res_id_;
}

void
ResConnID::resid( Size setting )
{
	res_id_ = setting;
}

Size ResConnID::connid() const
{
	return conn_id_;
}

void ResConnID::connid( Size setting )
{
	conn_id_ = setting;
}

bool ResConnID::incomplete() const
{
	return res_id_ == 0 || conn_id_ == 0;
}

void ResConnID::mark_incomplete()
{
	res_id_ = conn_id_ = 0;
}

}
}

// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/sewing/data_storage/Basis.cc
/// @brief Class for storing residue information needed to generate alignments
/// @author Minnie Langlois (minnie@email.unc.edu)

#include <protocols/sewing/data_storage/Basis.hh>
#include <basic/Tracer.hh>


//External headers
#include <boost/unordered_map.hpp>
#include <boost/functional/hash.hpp>
#include <boost/tuple/tuple.hpp>
#include <boost/tuple/tuple_comparison.hpp>





static basic::Tracer TR( "protocols.sewing.data_storage.Basis" );


namespace protocols {
namespace sewing {
namespace data_storage {

Basis::Basis():
	utility::pointer::ReferenceCount(),
	segment_id_(0),
	resnum_(0)
{}

Basis::Basis( core::Size segment_id, core::Size resnum ){
	segment_id_ = segment_id;
	resnum_ = resnum;
}

Basis::~Basis(){}

Basis::Basis( Basis const & src ) {
	segment_id_ = src.segment_id_;
	resnum_  = src.resnum_;
}

BasisOP
Basis::clone() const {
	return BasisOP( new Basis( *this ) );
}

bool
operator< ( Basis const & a, Basis const & b ){
	return boost::tie( a.segment_id_, a.resnum_) < boost::tie(b.segment_id_, b.resnum_);
}

void
Basis::segment_id( core::Size segment_id ){
	segment_id_ = segment_id;
}

void
Basis::resnum( core::Size resnum ){
	resnum_ = resnum;
}

Basis &
Basis::operator=( Basis const & other ){
	segment_id_ = other.segment_id();
	resnum_ = other.resnum();
	return *this;
}


} //data_storage
} //protocols
} //sewing







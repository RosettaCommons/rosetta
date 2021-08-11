// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/pose_sewing/data_storage/SegmentEnvelope.cc
/// @brief a descriptor of properties a SewAnything secondar structral segment may have
/// @author frankdt (frankdt@email.unc.edu)

#include <protocols/pose_sewing/data_storage/SegmentEnvelope.hh>
#include <basic/Tracer.hh>

static basic::Tracer TR( "protocols.pose_sewing.data_storage.SegmentEnvelope" );


namespace protocols {
namespace pose_sewing {
namespace data_storage {

SegmentEnvelope::SegmentEnvelope():
	utility::VirtualBase()
{

}

SegmentEnvelope::~SegmentEnvelope(){}

SegmentEnvelope::SegmentEnvelope( SegmentEnvelope const & ) = default;



SegmentEnvelopeOP
SegmentEnvelope::clone() const {
	return SegmentEnvelopeOP( new SegmentEnvelope( *this ) );
}

void
SegmentEnvelope::set_minimum_length(core::Size new_min){
	minimum_length_ = new_min;
}

void
SegmentEnvelope::set_maximum_length(core::Size new_max){
	maximum_length_ = new_max;
}

void
SegmentEnvelope::set_permissible_secondary_structures(std::string & new_perms){
	permissible_secondary_structures_ = new_perms;
}

core::Size
SegmentEnvelope::get_minimum_length(){
	return minimum_length_;
}

core::Size
SegmentEnvelope::get_maximum_length(){
	return maximum_length_;
}

std::string
SegmentEnvelope::get_permissible_secondary_structures(){
	return permissible_secondary_structures_;
}

bool
SegmentEnvelope::is_permissible_secondary_structure(char DSSP){
	if ( permissible_secondary_structures_.find(DSSP)!=std::string::npos ) {
		return true;
	} else {
		return false;
	}
}
bool
SegmentEnvelope::is_permissible_length(core::Size length){
	if ( length >= minimum_length_ && length <= maximum_length_ ) {
		return true;
	} else {
		return false;
	}
}
bool
SegmentEnvelope::is_valid(char DSSP, core::Size length){
	if ( this->is_permissible_secondary_structure(DSSP) && this->is_permissible_length(length) ) {
		return true;
	} else {
		return false;
	}
}

bool
SegmentEnvelope::is_valid(PoseSegmentOP pose_segment){
	return this->is_valid(pose_segment->get_dssp_code(),pose_segment->get_length());
}

} //protocols
} //pose_sewing
} //data_storage







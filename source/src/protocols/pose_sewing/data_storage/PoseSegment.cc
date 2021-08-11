// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file protocols/pose_sewing/data_storage/PoseSegment.cc
/// @brief a region of a Pose with contiguous secondary structure
/// @author frankdt (frankdt@email.unc.edu)

#include <protocols/pose_sewing/data_storage/PoseSegment.hh>
#include <basic/Tracer.hh>

static basic::Tracer TR( "protocols.pose_sewing.data_storage.PoseSegment" );


namespace protocols {
namespace pose_sewing {
namespace data_storage {

PoseSegment::PoseSegment():
	utility::VirtualBase() , source_pose_(nullptr)
{
	dssp_code_ = 'X';
}

PoseSegment::PoseSegment(core::Size start,core::Size end, char dssp, core::pose::PoseCOP posecop):
	utility::VirtualBase(), source_pose_(posecop) {
	starting_residue_=start;
	ending_residue_=end;
	dssp_code_ = dssp;
	//source_pose_=posecop;
}

PoseSegment::~PoseSegment(){}

PoseSegment::PoseSegment( PoseSegment const & ) = default;

PoseSegment::PoseSegment( PoseSegment const & src, core::Size start,core::Size end) : source_pose_(src.get_source_pose()) {
	//source_pose_ = src.get_source_pose();
	dssp_code_ = src.get_dssp_code();
	starting_residue_ = start;
	ending_residue_ = end;
}



PoseSegmentOP
PoseSegment::clone() const {
	return PoseSegmentOP( new PoseSegment( *this ) );
}

void
PoseSegment::set_starting_residue(core::Size new_starting){
	starting_residue_ = new_starting;
}
void
PoseSegment::set_ending_residue(core::Size new_ending){
	ending_residue_ = new_ending;
}
void
PoseSegment::set_dssp_code(char new_code){
	dssp_code_ = new_code;
}
void
PoseSegment::set_source_pose(core::pose::PoseCOP new_posecop){
	source_pose_ = new_posecop;
}

core::Size
PoseSegment::get_starting_residue() const {
	return starting_residue_;
}
core::Size
PoseSegment::get_ending_residue() const {
	return ending_residue_;
}
char
PoseSegment::get_dssp_code() const {
	return dssp_code_;
}
core::pose::PoseCOP
PoseSegment::get_source_pose() const {
	return source_pose_;
}
core::Size
PoseSegment::get_length() const {
	return ending_residue_ - starting_residue_ + 1;
}

} //protocols
} //pose_sewing
} //data_storage







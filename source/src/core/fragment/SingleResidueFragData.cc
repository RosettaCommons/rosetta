// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/fragments/FragData.cc
/// @brief  a collection classes of the FragData and SingleResidueFragData class hirachy
/// @author Oliver Lange (olange@u.washington.edu)
/// @date   Wed Oct 20 12:08:31 2007
///

// Unit Headers
#include <core/fragment/SingleResidueFragData.hh>

// Package Headers
//#include <core/fragment/FragData.hh>
#include <core/fragment/Frame.hh>

// Project Headers
// AUTO-REMOVED #include <core/kinematics/MoveMap.hh>
#include <core/pose/Pose.hh>
#include <core/conformation/Residue.hh>

#include <utility/vector1.hh>


// ObjexxFCL Headers
//#include <ObjexxFCL/format.hh>

// Utility headers
//#include <utility/vector1.fwd.hh>
//#include <basic/Tracer.hh>

namespace core {
namespace fragment {

/// @details Auto-generated virtual destructor
SingleResidueFragData::~SingleResidueFragData() {}

//static basic::Tracer tr("core.fragment");

bool SingleResidueFragData::steal(pose::Pose const& pose,  Size seq_pos) {
  sequence_ = oneletter_code_from_aa( pose.residue( seq_pos ).aa() );
  return true;
}

// insert fragment_data into pose at position given by Frame.seqpos( intra_frame_pos );
bool SingleResidueFragData::apply(pose::Pose &pose, Size const intra_frame_pos, Frame const& frame) const {
	return apply( pose, frame.seqpos( intra_frame_pos ) );
}

/// @brief insert fragment_data into pose at position given by Frame.seqpos( intra_frame_pos )
///  for dofs that are allowed to move in the MoveMap
bool SingleResidueFragData::apply( kinematics::MoveMap const & movemap, pose::Pose & pose, Size const intra_frame_pos, Frame const & frame ) const {
	return apply( movemap, pose, frame.seqpos( intra_frame_pos ) );
}

// insert fragment_data sec-struct into ss-string at position seq_pos
bool SingleResidueFragData::apply_ss( std::string& ss, Size intra_frame_pos, Frame const& frame) const {
  return apply_ss( ss, frame.seqpos( intra_frame_pos ) );
}

// insert fragment_data into pose at position seq_pos
bool SingleResidueFragData::steal(pose::Pose const& pose, Size intra_frame_pos, Frame const& frame) {
  return steal( pose, frame.seqpos( intra_frame_pos ) );
}

// check weather dofs can be moved
bool SingleResidueFragData::is_applicable( kinematics::MoveMap const& mm, Size intra_frame_pos, Frame const& frame) const {
  return is_applicable( mm, frame.seqpos( intra_frame_pos ) );
}

void SingleResidueFragData::show( std::ostream &out ) const {
	out << type() << " ";
	//out << "SINGLE RESIDUE FRAG DATA -- show() not overloaded for class " << typeid(*this).name() << std::endl;
}


void SingleResidueFragData::read_data( std::istream& ) {
	//nothing to read
}


std::string
SingleResidueFragData::type() const {
	return "_SRFD_Base_ERROR_";
}

std::string
SingleResidueFragData::_static_type_name() {
	return "SingleResidueFragData";
}


}
}

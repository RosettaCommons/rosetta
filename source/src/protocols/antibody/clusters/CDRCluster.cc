// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available
// (c) under license. The Rosetta software is developed by the contributing
// (c) members of the Rosetta Commons. For more information, see
// (c) http://www.rosettacommons.org. Questions about this can be addressed to
// (c) University of Washington UW TechTransfer,email:license@u.washington.edu.

/// @file protocols/antibody/clusters/CDRCluster.cc
/// @brief Simple class to hold and access CDRCluster information of a region in a pose.
/// @author Jared Adolf-Bryfogle (jadolfbr@gmail.com)

#include <protocols/antibody/clusters/CDRCluster.hh>
#include <core/pose/PDBInfo.hh>
#include <math.h>

namespace protocols {
namespace antibody {
namespace clusters {
	using namespace protocols::antibody;


CDRCluster::CDRCluster(core::pose::Pose const & pose, CDRNameEnum const cdr, core::Size const cdr_length,
		CDRClusterEnum const cluster, core::Size const start, core::Real const distance):
		cdr_(cdr),
		cluster_(cluster),
		distance_(distance),
		start_(start),
		length_(cdr_length)
{
	end_ = start+cdr_length-1;
	normalized_distance_ = distance/(cdr_length*2);
	set_pdb_numbering(pose, start, end_);
}

CDRCluster::~CDRCluster(){}

core::Real
CDRCluster::normalized_distance_in_degrees() const {
	
	core::Real result = acos(1-(distance_/(2*length_)));
	return result;
}

void
CDRCluster::set_pdb_numbering(core::pose::Pose const & pose, core::Size start, core::Size end){
	pdb_start_ = pose.pdb_info()->number(start);
	pdb_end_ = pose.pdb_info()->number(end);
	
	pdb_start_insertion_code_ = pose.pdb_info()->icode(start);
	pdb_end_insertion_code_ = pose.pdb_info()->icode(end);
	chain_ = pose.pdb_info()->chain(start);
	
}

} //design
} //antibody
} //protocols

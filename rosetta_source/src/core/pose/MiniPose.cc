// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file   core/pose/MiniPose.cc
/// @brief  MiniPose class
/// @details minimal class with xyz and fold_tree but not the atom_tree + energies machinery...
/// @author Rhiju Das

// Unit headers
#include <core/pose/Pose.hh>
#include <core/pose/MiniPose.hh>
#include <core/types.hh>

// Project headers
#include <core/conformation/Residue.hh>
#include <core/id/AtomID.hh>
#include <core/kinematics/FoldTree.hh>
#include <utility/vector1.hh>
#include <numeric/xyzVector.hh>

#include <string>

namespace core {
namespace pose {

	///////////////////////////////////////////////////////////////////////
	MiniPose::MiniPose( core::pose::Pose const & pose )
	{
		for ( Size i = 1; i <= pose.total_residue(); i++ ) {
			utility::vector1< PointPosition > xyz;
			for ( Size j = 1; j <= pose.residue(i).natoms(); j++ ) {
				xyz.push_back( pose.residue(i).xyz( j ) ) ;
			}
			coords_.push_back( xyz );
		}
		sequence_ = pose.sequence();
		fold_tree_ = pose.fold_tree();
	}


	core::kinematics::FoldTree const &
	MiniPose::fold_tree() const{
		return fold_tree_;
	}

	utility::vector1< utility::vector1< PointPosition > > const &
	MiniPose::coords() const{
		return coords_;
	}

	Size
	MiniPose::size() const {
		return coords_.size();
	}

	std::string const &
	MiniPose::sequence() const {
		return sequence_;
	}

	PointPosition const &
	MiniPose::xyz( core::id::AtomID atom_id ) const{
		return coords_[ atom_id.rsd() ][ atom_id.atomno() ];
	}


} // pose
} // core

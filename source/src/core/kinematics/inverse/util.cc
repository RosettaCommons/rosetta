// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

///@file   core/kinematics/inverse/jump.cc
///@brief  Utility functions for calculating jumps by knowing desired atom positions
///@author Jack Maguire

// Unit headers
#include <core/kinematics/inverse/jump.hh>
#include <core/kinematics/inverse/AlignmentAtom.hh>

#include <basic/Tracer.hh>
#include <numeric/HomogeneousTransform.hh>
#include <core/kinematics/Jump.hh>
#include <core/kinematics/FoldTree.hh>
#include <core/kinematics/Edge.hh>

#include <utility/excn/Exceptions.hh>
//#include <stdexcept>

namespace core {
namespace kinematics {
namespace inverse {

static basic::Tracer TR( "core.kinematics.inverse.util" );

void
assert_atoms_are_downstream_of_jump(
	conformation::Conformation const & conformation,
	core::Size const jump_id,
	AlignmentAtomArray const & atom_arr
){
	utility::vector1< bool > is_downstream =
		conformation.fold_tree().partition_by_jump( jump_id );
	core::Size const resid_that_should_be_downstream =
		conformation.fold_tree().jump_edge( jump_id ).stop();

	//Flip partition if needed
	if ( !is_downstream[ resid_that_should_be_downstream ] ) {
		for ( core::Size ii = 1; ii <= is_downstream.size(); ++ii ) {
			is_downstream[ ii ] = !is_downstream[ ii ];
		}
	}

	for ( AlignmentAtom const & aa : atom_arr.atoms ) {
		if ( ! is_downstream[ aa.id.rsd() ] ) {
			throw CREATE_EXCEPTION( utility::excn::BadInput, "Residue " + std::to_string(aa.id.rsd()) + " is being used to calculate the destination of jump " + std::to_string(jump_id) + " even though this residue is upstream of this jump." );
		}
	}
}

void
assert_atoms_are_upstream_of_jump(
	conformation::Conformation const & conformation,
	core::Size const jump_id,
	AlignmentAtomArray const & atom_arr
){
	utility::vector1< bool > is_downstream =
		conformation.fold_tree().partition_by_jump( jump_id );
	core::Size const resid_that_should_be_upstream =
		conformation.fold_tree().jump_edge( jump_id ).start();

	//Flip partition if needed
	if ( is_downstream[ resid_that_should_be_upstream ] ) {
		for ( core::Size ii = 1; ii <= is_downstream.size(); ++ii ) {
			is_downstream[ ii ] = !is_downstream[ ii ];
		}
	}

	for ( AlignmentAtom const & aa : atom_arr.atoms ) {
		if ( is_downstream[ aa.id.rsd() ] ) {
			throw CREATE_EXCEPTION( utility::excn::BadInput, "Residue " + std::to_string(aa.id.rsd()) + " is being used to calculate the base of jump " + std::to_string(jump_id) + " even though this residue is downstream of this jump." );
		}
	}
}

} // namespace inverse
} // namespace kinematics
} // namespace core

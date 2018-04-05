// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file src/core/pack/rotamer_set/bb_independent_rotamers.hh
/// @brief a bunch of utility functions used in enzdes
/// @author Florian Richter, floric@u.washington.edu

#include <core/pack/rotamer_set/bb_independent_rotamers.hh>

#include <utility/vector1.hh>
#include <core/conformation/Residue.hh>
#include <core/pack/task/PackerTask.hh>
#include <core/pose/Pose.hh>
#include <core/pose/variant_util.hh>
#include <core/scoring/ScoreFunction.hh>

#include <utility/graph/Graph.hh>

#include <core/pack/task/TaskFactory.hh>
#include <core/pack/packer_neighbors.hh>
#include <core/pack/rotamer_set/RotamerSetOperation.hh>
#include <core/pack/rotamer_set/RotamerSetFactory.hh>
#include <core/pack/rotamer_set/RotamerSet.hh>

namespace core {
namespace pack {
namespace rotamer_set {


utility::vector1< core::conformation::ResidueCOP >
bb_independent_rotamers(
	core::chemical::ResidueTypeCOP rot_restype,
	bool ignore_cmdline
)
{
	core::conformation::Residue firstres( *rot_restype, true );
	core::pose::Pose dummy_pose;
	dummy_pose.append_residue_by_jump( firstres, (core::Size) 0 );
	if ( rot_restype->is_polymer() ) {
		core::pose::add_lower_terminus_type_to_pose_residue( dummy_pose, 1 ); //prolly critical so that the dunbrack library uses neutral phi
		core::pose::add_upper_terminus_type_to_pose_residue( dummy_pose, 1 ); //prolly critical so that the dunbrack library uses neutral psi
	}
	core::scoring::ScoreFunction dummy_sfxn;
	dummy_sfxn( dummy_pose );
	core::pack::task::PackerTaskOP dummy_task = core::pack::task::TaskFactory::create_packer_task( dummy_pose );
	if ( !ignore_cmdline ) dummy_task->initialize_from_command_line();
	dummy_task->nonconst_residue_task( 1 ).restrict_to_repacking();
	dummy_task->nonconst_residue_task( 1 ).or_include_current( false ); //need to do this because the residue was built from internal coords and is probably crumpled up
	dummy_task->nonconst_residue_task( 1 ).or_fix_his_tautomer( true ); //since we only want rotamers for the specified restype
	utility::graph::GraphOP dummy_png = core::pack::create_packer_graph( dummy_pose, dummy_sfxn, dummy_task );

	core::pack::rotamer_set::RotamerSetOP rotset( RotamerSetFactory::create_rotamer_set( dummy_pose ) ); // build_rotamers should be symmetry independent
	rotset->set_resid( 1 );
	rotset->build_rotamers( dummy_pose, dummy_sfxn, *dummy_task, dummy_png );

	utility::vector1< core::conformation::ResidueCOP > to_return;

	//now when creating the rotamers, we have to make sure we don't sneak in the additional variant types
	for ( core::Size i = 1; i <= rotset->num_rotamers(); ++i ) {
		core::conformation::ResidueOP rot( firstres.clone() );
		for ( core::Size j =1; j <= firstres.nchi(); ++j ) rot->set_chi( j, rotset->rotamer( i )->chi( j ) );
		to_return.push_back( rot );
	}

	return to_return;
}

}
}
}

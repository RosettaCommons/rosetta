// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/pack/guidance_scoreterms/approximate_buried_unsat_penalty/util.cc
/// @brief  Utility functions for ApproximateBuriedUnsatPenalty
/// @author Brian Coventry (bcov@uw.edu)

// #include <core/pack/guidance_scoreterms/approximate_buried_unsat_penalty/util.hh>

#include <core/id/AtomID.hh>
#include <core/conformation/Residue.hh>
#include <core/pack/rotamer_set/RotamerSet_.hh>
#include <core/pack/task/TaskFactory.hh>
#include <core/pack/hbonds/HBondGraph_util.hh>
#include <core/scoring/atomic_depth/AtomicDepth.hh>
#include <core/scoring/ScoreFunction.hh>

#include <utility/graph/Graph.hh>
#include <boost/format.hpp>

namespace core {
namespace pack {
namespace guidance_scoreterms {
namespace approximate_buried_unsat_penalty {



/// @brief Construct an AtomLevelHBondGraph from a partial rotsets (like you might see during packing)
scoring::hbonds::graph::AtomLevelHBondGraphOP
hbond_graph_from_partial_rotsets(
	pose::Pose const & pose_in,
	pack::rotamer_set::RotamerSetsOP const & original_rotsets,
	scoring::ScoreFunctionOP const & scorefxn,
	pack::rotamer_set::RotamerSetsOP & complete_rotsets_out,
	utility::vector1<bool> & position_had_rotset,
	float minimum_hb_cut /* =0 */
) {
	// OK, so first we need to make a RotamerSets that has an entry for every single seqpos.
	//     We'll reuse whatever we can from original_rotsets and then start adding the current rotamers.

	// We need to score the pose. So we should probably make a copy so that we don't confuse the packer.

	pose::Pose pose = pose_in;

	scorefxn->score( pose );

	// Pretend like we're going to design every position
	utility::vector1<bool> true_vect( pose.size(), true );
	scorefxn->setup_for_packing( pose, true_vect, true_vect );

	// Blank packer task that says we're going to pack every position
	pack::task::PackerTaskOP task = pack::task::TaskFactory::create_packer_task( pose );

	// Our internal rotsets that the AtomLevelHBondGraph needs
	pack::rotamer_set::RotamerSetsOP rotsets( new pack::rotamer_set::RotamerSets() );
	rotsets->set_task( task );

	position_had_rotset.clear();
	position_had_rotset.resize( pose.size() );

	// Fill in the rotsets with either the RotamerSet from original_rotsets or make a new one from the current res
	for ( Size resnum = 1; resnum <= pose.size(); resnum++ ) {

		// We need to make a new rotset in either case
		//  Either because there was no rotset for this position
		//  Or because we'll confuse the scorefunction machinery if we don't make a new one
		pack::rotamer_set::RotamerSetOP rotset( new pack::rotamer_set::RotamerSet_() );
		rotset->set_resid( resnum );

		if ( original_rotsets->has_rotamer_set_for_residue( resnum ) ) {
			pack::rotamer_set::RotamerSetCOP original_rotset = original_rotsets->rotamer_set_for_residue( resnum );
			for ( Size irot = 1; irot <= original_rotset->num_rotamers(); irot++ ) {
				rotset->add_rotamer( *original_rotset->rotamer( irot ) );
			}
			position_had_rotset[resnum] = true;
		} else {
			rotset->add_rotamer( *pose.residue(resnum).create_rotamer() );
			position_had_rotset[resnum] = false;
		}

		scorefxn->prepare_rotamers_for_packing( pose, *rotset );
		rotsets->set_explicit_rotamers( rotsets->resid_2_moltenres( resnum ), rotset );
	}
	rotsets->update_offset_data();

	// Use some of the HBNet machinery to load the AtomLevelHBondGraph
	scoring::hbonds::graph::AtomLevelHBondGraphOP hb_graph;
	hb_graph = pack::hbonds::create_init_and_create_edges_for_atom_level_hbond_graph(
		rotsets,    // our temporary rotsets
		* scorefxn, pose,
		0,              // allow all hbonds during interation graph generation
		1e6,            // allow all clashes
		minimum_hb_cut, // only allow good hbonds to make it into the final graph
		true            // include backbone atoms
	);

	complete_rotsets_out = rotsets;

	return hb_graph;
}




}
}
}
}

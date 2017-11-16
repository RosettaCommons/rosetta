// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   core/pack/rotamer_trials.cc
/// @brief  rotamer trials module
/// @author Andrew Leaver-Fay (leaverfa@email.unc.edu)

// Unit headers
#include <core/pack/rotamer_trials.hh>

// Package headers

#include <core/pack/packer_neighbors.hh>
#include <core/pack/rotamer_set/RotamerSet.hh>
#include <core/pack/rotamer_set/RotamerSetFactory.hh>
#include <core/pack/rotamer_set/symmetry/SymmetricRotamerSetFactory.hh>
#include <core/conformation/symmetry/SymmetryInfo.hh>
#include <core/pack/task/PackerTask.hh>
#include <core/pack/make_symmetric_task.hh>

// Project headers
#include <core/types.hh>

#include <utility/graph/Graph.hh>

#include <core/pose/Pose.hh>

#include <core/scoring/ScoreFunction.hh>

#include <core/conformation/Residue.hh>

#include <basic/prof.hh> // should this guy go into utilities?
#include <basic/Tracer.hh>

// Utility headers
#include <utility/vector1.hh>
#include <utility/vector1.functions.hh>

// Numeric Headers
#include <numeric/random/random.hh>
#include <numeric/random/random_permutation.hh>

#include <core/pose/symmetry/util.hh>


// STL headers
#ifdef WIN32
#include <ctime>
#endif

#include <core/conformation/symmetry/SymmetricConformation.hh>
#include <utility/vector0.hh>


namespace core {
namespace pack {

typedef conformation::symmetry::SymmetricConformation SymmetricConformation;
typedef conformation::symmetry::SymmetryInfo SymmetryInfo;


static basic::Tracer TR( "core.pack.rotamer_trials" );

utility::vector1< uint >
symmetric_repackable_residues(
	task::PackerTask const & the_task,
	pose::Pose & pose
);

void
rotamer_trials(
	pose::Pose & pose,
	scoring::ScoreFunction const & scfxn,
	task::PackerTaskCOP input_task
) {
	using namespace numeric::random;

	//fpd safety check for symmetry
	if ( core::pose::symmetry::is_symmetric( pose ) ) {
		symmetric_rotamer_trials( pose, scfxn, input_task );
		return;
	}

	//clock_t starttime = clock();
	PROF_START( basic::ROTAMER_TRIALS );

	pack_scorefxn_pose_handshake( pose, scfxn);

	pose.update_residue_neighbors();

	utility::vector1< uint > residues_for_trials( repackable_residues( *input_task ));
	// Replace this random shuffle with one based on rosetta's RNG
	//random__shuffle(residues_for_trials.begin(), residues_for_trials.end() );
	random_permutation( residues_for_trials, numeric::random::rg() );

	task::PackerTaskOP rottrial_task( input_task->clone() );

	rottrial_task->set_bump_check( false );
	rottrial_task->or_include_current( true );
	rottrial_task->temporarily_fix_everything();

	// this will call setup fxns for each scoring method, eg HBondEnergy will
	// compute backbone hbonds to prepare for hbchecking,
	// PairEnergy will update actcoords...
	scfxn.setup_for_packing( pose, rottrial_task->repacking_residues(), rottrial_task->designing_residues() );

	rotamer_set::RotamerSetFactory rsf;

	//Energy last_energy( 0.0 );
	utility::graph::GraphOP packer_neighbor_graph = create_packer_graph( pose, scfxn, input_task );

	bool replaced_residue( false );

	Size const num_in_trials = residues_for_trials.size();
	TR.Trace << "Performing rotamer trails on " << num_in_trials << " residues..." << std::endl;
	for ( Size ii = 1; ii <= num_in_trials; ++ii ) {
		pose.update_residue_neighbors(); // will return if uptodate

		int resid = residues_for_trials[ ii ];
		conformation::Residue const & trial_res = pose.residue( resid );

		//pretend this is a repacking and only this residue is being repacked
		//while all other residues are being held fixed.
		rottrial_task->temporarily_set_pack_residue( resid, true );

		rotamer_set::RotamerSetOP rotset = rsf.create_rotamer_set( trial_res );
		rotset->set_resid( resid );
		rotset->build_rotamers( pose, scfxn, *rottrial_task, packer_neighbor_graph );

		scfxn.prepare_rotamers_for_packing( pose, *rotset );

		utility::vector1< core::PackerEnergy > one_body_energies( rotset->num_rotamers() );

		rotset->compute_one_body_energies( pose, scfxn, *rottrial_task, packer_neighbor_graph, one_body_energies );

		//select the best rotamer
		Size bestrot = utility::arg_min( one_body_energies );

		//don't replace if the best rotamer is the one that's already assigned
		TR.Trace << "rottrial at position " << resid << " nrot= " << rotset->num_rotamers() << std::endl;
		if ( bestrot != rotset->id_for_current_rotamer() ) {
			replaced_residue = true;
			conformation::ResidueOP newresidue( rotset->rotamer( bestrot )->clone() );//create_residue() );
			pose.replace_residue ( resid, *newresidue, false );
			scfxn.update_residue_for_packing( pose, resid );
			TR.Trace << "rottrial accept: " << resid << " bestrot: " << bestrot << ' ' << one_body_energies[ bestrot ];
			if ( rotset->id_for_current_rotamer() != 0 ) { //more output in this case
				//sml this output is protected because id_for_current_rotamer is incompatible with forced mutations (eg PIKAA)
				TR.Trace << ' ' << one_body_energies[ rotset->id_for_current_rotamer() ] << ' ' <<
					one_body_energies[ bestrot ] - one_body_energies[ rotset->id_for_current_rotamer() ];
			}
			TR.Trace << std::endl;
		}

		rottrial_task->temporarily_set_pack_residue( resid, false );

	}
	//clock_t stoptime = clock();
	PROF_STOP ( basic::ROTAMER_TRIALS );
	//std::cout << "Rotamer trials took " << ((double) stoptime - starttime ) / CLOCKS_PER_SEC << " seconds" << std::endl;
	if ( replaced_residue ) {
		// rescore here to make Energies GOOD
		scfxn( pose );
	}
}

utility::vector1< uint >
repackable_residues( task::PackerTask const & the_task )
{
	utility::vector1< int > to_be_packed( the_task.num_to_be_packed() );
	uint count = 0;
	for ( uint ii = 1; ii <= the_task.total_residue(); ++ii ) {
		if ( the_task.pack_residue( ii ) ) {
			++count;
			to_be_packed[ count ] = ii;
		}
	}
	return to_be_packed;
}

void
symmetric_rotamer_trials(
	pose::Pose & pose,
	scoring::ScoreFunction const & scfxn,
	task::PackerTaskCOP non_symmetric_task
) {
	using namespace numeric::random;
	using namespace conformation::symmetry;

	//clock_t starttime = clock();
	PROF_START( basic::ROTAMER_TRIALS );

	task::PackerTaskCOP input_task = make_new_symmetric_PackerTask_by_requested_method( pose, non_symmetric_task );

	pack_scorefxn_pose_handshake( pose, scfxn);

	pose.update_residue_neighbors();

	utility::vector1< uint > residues_for_trials( symmetric_repackable_residues( *input_task, pose ));
	// Replace this random shuffle with one based on rosetta's RNG
	//random__shuffle(residues_for_trials.begin(), residues_for_trials.end() );
	random_permutation( residues_for_trials, numeric::random::rg() );

	task::PackerTaskOP rottrial_task( input_task->clone() );

	rottrial_task->set_bump_check( false );
	rottrial_task->or_include_current( true );
	rottrial_task->temporarily_fix_everything();

	// this will call setup fxns for each scoring method, eg HBondEnergy will
	// compute backbone hbonds to prepare for hbchecking,
	// PairEnergy will update actcoords...
	scfxn.setup_for_packing( pose, rottrial_task->repacking_residues(), rottrial_task->designing_residues() );

	rotamer_set::symmetry::SymmetricRotamerSetFactory rsf;

	//Energy last_energy( 0.0 );
	utility::graph::GraphOP packer_neighbor_graph = create_packer_graph( pose, scfxn, input_task );

	bool replaced_residue( false );

	Size const num_in_trials = residues_for_trials.size();
	for ( Size ii = 1; ii <= num_in_trials; ++ii ) {
		pose.update_residue_neighbors(); // will return if uptodate

		int resid = residues_for_trials[ ii ];
		conformation::Residue const & trial_res = pose.residue( resid );

		//pretend this is a repacking and only this residue is being repacked
		//while all other residues are being held fixed.
		rottrial_task->temporarily_set_pack_residue( resid, true );

		rotamer_set::RotamerSetOP rotset = rsf.create_rotamer_set( trial_res );
		rotset->set_resid( resid );
		rotset->initialize_pose_for_rotset_creation( pose );
		rotset->build_rotamers( pose, scfxn, *rottrial_task, packer_neighbor_graph );

		scfxn.prepare_rotamers_for_packing( pose, *rotset );

		utility::vector1< core::PackerEnergy > one_body_energies( rotset->num_rotamers() );

		rotset->compute_one_body_energies( pose, scfxn, *rottrial_task, packer_neighbor_graph, one_body_energies );

		//select the best rotamer
		Size bestrot = utility::arg_min( one_body_energies );

		//don't replace if the best rotamer is the one that's already assigned
		//  TR.Trace << "rottrial at position " << resid << " nrot= " << rotset->num_rotamers() << std::endl;
		if ( bestrot != rotset->id_for_current_rotamer() ) {
			replaced_residue = true;
			conformation::ResidueOP newresidue(  rotset->rotamer( bestrot )->clone() );//create_residue() );
			pose.replace_residue ( resid, *newresidue, false );
			scfxn.update_residue_for_packing( pose, resid );
			// clone to symmetrical positions
			SymmetricConformation const & SymmConf ( dynamic_cast<SymmetricConformation const &> ( pose.conformation()) );
			SymmetryInfoCOP symm_info( SymmConf.Symmetry_Info() );

			//fpd replace_residue is symmetric now ... still need to update scorefxn though
			for ( core::Size clone : symm_info->bb_clones( resid ) ) {
				//    conformation::ResidueOP sym_rsd = newresidue->clone();
				//     sym_rsd->orient_onto_residue(pose.residue( *clone) );
				//    pose.replace_residue ( *clone, *sym_rsd, false );
				scfxn.update_residue_for_packing( pose, clone );
			}
			//   TR.Trace << "rottrial accept: " << resid << " bestrot: " << bestrot << ' ' << one_body_energies[ bestrot ];
			if ( rotset->id_for_current_rotamer() != 0 ) { //more output in this case
				//sml this output is protected because id_for_current_rotamer is incompatible with forced mutations (eg PIKAA)
				// TR.Trace << ' ' << one_body_energies[ rotset->id_for_current_rotamer() ] << ' ' <<
				//     one_body_energies[ bestrot ] - one_body_energies[ rotset->id_for_current_rotamer() ];

			}
			//TR.Trace << std::endl;
		}
		rottrial_task->temporarily_set_pack_residue( resid, false );
	}
	//clock_t stoptime = clock();
	PROF_STOP ( basic::ROTAMER_TRIALS );
	//std::cout << "Rotamer trials took " << ((double) stoptime - starttime ) / CLOCKS_PER_SEC << " seconds" << std::endl;
	if ( replaced_residue ) {
		// rescore here to make Energies GOOD
		scfxn( pose );
	}

}

utility::vector1< uint >
symmetric_repackable_residues(
	task::PackerTask const & the_task,
	pose::Pose & pose
)
{
	using namespace conformation::symmetry;

	// find SymmInfo
	SymmetricConformation const & SymmConf (
		dynamic_cast<SymmetricConformation const &> ( pose.conformation() ) );
	SymmetryInfoCOP symm_info( SymmConf.Symmetry_Info() );

	// First find out how many residues to pack. Silly, lets come up with a better way of sizing this array...
	uint num_molten = 0;
	for ( uint res = 1; res <= the_task.total_residue(); ++res ) {
		if ( the_task.pack_residue( res ) && symm_info->fa_is_independent( res )
				&& res <= symm_info->num_total_residues() ) {
			++num_molten;
		}
	}
	utility::vector1< int > to_be_packed( num_molten );
	uint count = 0;
	for ( uint ii = 1; ii <= the_task.total_residue(); ++ii ) {
		if ( the_task.pack_residue( ii ) && symm_info->fa_is_independent( ii )
				&& ii <= symm_info->num_total_residues() ) {
			++count;
			to_be_packed[ count ] = ii;
		}
	}
	return to_be_packed;
}


} //end namespace core
} //end namespace pack

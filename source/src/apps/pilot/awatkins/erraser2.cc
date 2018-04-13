// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file apps/pilot/awatkins/erraser2.cc
/// @brief a single erraser app. no python!
/// @author Andy Watkins, amw579@stanford.edu


// protocols
#include <protocols/viewer/viewers.hh>
#include <protocols/stepwise/monte_carlo/util.hh>
#include <protocols/stepwise/monte_carlo/options/StepWiseMonteCarloOptions.hh>
#include <protocols/stepwise/monte_carlo/mover/StepWiseMasterMover.hh>
#include <protocols/stepwise/modeler/rna/util.hh>
#include <protocols/stepwise/modeler/util.hh>


// core
#include <core/types.hh>
#include <devel/init.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/Energies.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/rms_util.hh>
#include <core/scoring/rna/RNA_ScoringInfo.hh>
#include <core/chemical/ChemicalManager.hh>
#include <core/io/silent/SilentFileData.hh>
#include <core/io/silent/BinarySilentStruct.hh>
#include <core/import_pose/import_pose.hh>
#include <core/import_pose/pose_stream/PDBPoseInputStream.hh>
#include <core/pose/Pose.hh>
#include <core/pose/PDBInfo.hh>
#include <core/pose/util.hh>
#include <core/pose/rna/util.hh>
#include <core/pose/full_model_info/FullModelInfo.hh>
#include <core/pose/full_model_info/util.hh>
#include <core/import_pose/FullModelPoseBuilder.hh>
#include <protocols/rna/movers/ErraserMinimizerMover.hh>
#include <core/pose/rna/RNA_SuiteName.hh>

// basic
#include <utility/file/file_sys_util.hh>
#include <basic/Tracer.hh>
#include <basic/options/option.hh>
#include <basic/options/option_macros.hh>
#include <basic/options/keys/OptionKeys.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
//#include <basic/options/keys/file.OptionKeys.gen.hh>
#include <basic/options/keys/out.OptionKeys.gen.hh>
#include <basic/options/keys/chemical.OptionKeys.gen.hh>
#include <basic/options/keys/score.OptionKeys.gen.hh>

// utility
#include <utility/vector1.hh>
#include <utility/vector1.functions.hh>
#include <utility/excn/Exceptions.hh>
#include <utility/io/izstream.hh>

// C++ headers
#include <iostream>
#include <string>

using namespace protocols;
using namespace basic::options;
using namespace basic::options::OptionKeys;
using namespace core;
using namespace core::chemical;
using namespace core::pose;
using namespace core::scoring;
using namespace protocols::rna::movers;
using namespace core::pose::rna;
using namespace utility::file;
using namespace core::import_pose::pose_stream;

OPT_KEY( Integer, rounds )
OPT_KEY( Boolean, minimize_protein )

static basic::Tracer TR( "apps.pilot.awatkins.erraser2" );




std::string first_preminimized_namer( std::string const & s ) {
	std::stringstream minimized_name;
	minimized_name << s << "_preminimized_round_1.pdb";
	return minimized_name.str();
}
std::string second_preminimized_namer( std::string const & s ) {
	std::stringstream minimized_name;
	minimized_name << s << "_preminimized_round_2.pdb";
	return minimized_name.str();
}
std::string third_preminimized_namer( std::string const & s ) {
	std::stringstream minimized_name;
	minimized_name << s << "_preminimized_round_3.pdb";
	return minimized_name.str();
}
std::string round_minimized_namer( std::string const & s, Size const ii ) {
	std::stringstream minimized_name;
	minimized_name << s << "_minimized_" << ii << ".pdb";
	return minimized_name.str();
}
std::string round_resampled_namer( std::string const & s, Size const ii ) {
	std::stringstream minimized_name;
	minimized_name << s << "_resampled_" << ii << ".pdb";
	return minimized_name.str();
}


void fill_input_poses( utility::vector1< PoseOP > & input_poses, std::string const & s ) {
	TR << "Filling input pose from " << s << "..." << std::endl;

	// setup residue types
	ResidueTypeSetCOP rsd_set = ChemicalManager::get_instance()->residue_type_set( FA_STANDARD );
	PDBPoseInputStream input( s );
	// iterate over input stream
	PoseOP start_pose( new Pose );
	input.fill_pose( *start_pose, *rsd_set );
	input_poses.emplace_back( start_pose );
}

/// @details This function looks in the current working directory to see if
/// there are any checkpoint files. (These are based off the name of
/// option[ in::file::s ](). In fact, we have now abstracted those naming
/// systems into functions. It sets the 'checkpoints_to_pass' counter to
/// the number of checkpointing opportunities that can be crossed now that
/// we have recovered a file with this type of name...
void fill_input_poses_based_on_possible_checkpoints(
	utility::vector1< PoseOP > & input_poses,
	Size & checkpoints_to_pass
) {

	// Looks for 'latest possible' first.
	for ( Size round = option[ rounds ].value(); round >= 1; --round ) {
		if ( file_exists( round_resampled_namer( option[ in::file::s ]()[1], round ) ) ) {
			fill_input_poses( input_poses, round_resampled_namer( option[ in::file::s ]()[1], round ) );
			checkpoints_to_pass = 3 + 2 * round;
			return;
		}
		if ( file_exists( round_minimized_namer( option[ in::file::s ]()[1], round ) ) ) {
			fill_input_poses( input_poses, round_minimized_namer( option[ in::file::s ]()[1], round ) );
			checkpoints_to_pass = 3 + 2 * round - 1;
			return;
		}
	}

	// This won't happen, probably -- it's been disabled down below.
	if ( file_exists( third_preminimized_namer( option[ in::file::s ]()[1] ) ) ) {
		fill_input_poses( input_poses, third_preminimized_namer( option[ in::file::s ]()[1] ) );
		checkpoints_to_pass = 3;
		return;
	}
	if ( file_exists( second_preminimized_namer( option[ in::file::s ]()[1] ) ) ) {
		fill_input_poses( input_poses, second_preminimized_namer( option[ in::file::s ]()[1] ) );
		checkpoints_to_pass = 2;
		return;
	}
	if ( file_exists( first_preminimized_namer( option[ in::file::s ]()[1] ) ) ) {
		fill_input_poses( input_poses, first_preminimized_namer( option[ in::file::s ]()[1] ) );
		checkpoints_to_pass = 1;
		return;
	}

	fill_input_poses( input_poses, option[ in::file::s ]()[1] );
	checkpoints_to_pass = 0;
}

void resample_full_model( pose::Pose & start_pose, ScoreFunctionOP const & scorefxn,  utility::vector1< Size > const & definite_residues  );

utility::vector1< Size >
all_pose_residues( core::pose::Pose const & pose ) {
	utility::vector1< Size > foo;
	for ( Size ii = 1; ii <= pose.size(); ++ii ) foo.emplace_back( ii );
	return foo;
}

void show_accuracy_report( pose::Pose const & start_pose, std::string const & tag, Size const round ) {
	RNA_SuiteName suite_namer;

	utility::vector1< Size > bad_suites;
	for ( Size ii = 1; ii <= start_pose.size(); ++ii ) {
		// This should pick up bad suites
		auto suite_assignment = suite_namer.assign( start_pose, ii );
		if ( suite_assignment.name == "!!" ) bad_suites.push_back( ii );
	}

	// TODO: don't output round if # is zero
	// TODO: don't output bad suites if empty vector
	// TODO: output other things too.
	TR << tag << " " << round << ": " << bad_suites << std::endl;
}

utility::vector1< Size >
analyze_poses_for_convergence( utility::vector1< pose::Pose > const & poses ) {
	// Assume poses already aligned.
	utility::vector1< Size > res;
	utility::vector1< Real > RMSFs( poses[ 1 ].size() );

	for ( Size ii = 1; ii <= poses[ 1 ].size(); ++ii ) {
		utility::vector1< Real > dev_this_rsd_all_poses( poses.size(), 0 );
		for ( Size jj = 1; jj <= poses[ 1 ].residue_type( ii ).natoms(); ++jj ) {
			numeric::xyzVector< Real > avg_xyz( 0 );
			for ( auto const & pose : poses ) {
				avg_xyz += pose.residue( ii ).atom( jj ).xyz();
			}
			avg_xyz /= poses.size();

			for ( Size kk = 1; kk <= poses.size(); ++kk ) {
				dev_this_rsd_all_poses[ kk ] += poses[ kk ].residue( ii ).atom( jj ).xyz().distance_squared( avg_xyz );
			}
		}
		Real tot = 0;
		for ( Real const dev : dev_this_rsd_all_poses ) tot += dev;
		RMSFs[ ii ] = std::sqrt( 1.0/poses.size() /  poses[ 1 ].residue_type( ii ).natoms() * tot );
		TR.Trace << "RSD " << ii << " RMSF " << RMSFs[ ii ] << std::endl;
	}


	// OK which indices have the biggest RMSF?
	res.resize( 10 );
	arg_greatest_several( RMSFs, res );
	TR <<" Evaluating residues: " << res << std::endl;
	return res;
}

///////////////////////////////////////////////////////////////////////////////
void
erraser2_test()
{
	// Outline:
	// 1. Read in a pose.
	// 2. Minimize it (using the ErraserMinimizerMover)
	// 3. Rebuild residues that moved a lot plus geometric outliers
	//   a. implement via resample_full_model with a custom sample res
	// 4. [wash again if desired]
	// 5. Minimize it

	using namespace core::scoring;
	using namespace core::chemical;
	using namespace core::conformation;
	using namespace core::io::silent;
	using namespace core::import_pose;
	using namespace core::import_pose::pose_stream;
	using namespace core::pose;
	using namespace core::pose::rna;
	using namespace core::pose::full_model_info;

	// setup residue types
	ResidueTypeSetCOP rsd_set = ChemicalManager::get_instance()->residue_type_set( FA_STANDARD );

	// setup score function
	core::scoring::ScoreFunctionOP scorefxn;
	if ( option[ score::weights ].user () ) {
		scorefxn = get_score_function();
	} else {
		// Don't allow scorefunction to default.
		utility_exit_with_message( "You must provide a scoring function via -score:weights. Try -score:weights stepwise/rna/rna_res_level_energy4.wts with -set_weights elec_dens_fast 10.0." );
	}

	if ( option[ in::file::s ]().empty() ) {
		utility_exit_with_message( "You must provide a starting model via -in:file:s in PDB format." );
	}


	utility::vector1< PoseOP > input_poses;
	PoseOP start_pose( new Pose );

	utility::vector1< Size > unconverged_res;

	Size checkpoints_to_pass = 0;
	if ( option[ in::file::s ]().size() > 1 ) {
		// We interpret multiple provided files as an attempt to say "hey, I need to know which residues to rebuild!"
		PDBPoseInputStream input( option[ in::file::s ]() );
		Size ii = 1;
		utility::vector1< Pose > to_be_analyzed;
		// Still just put the first one into the system. The rest aren't other_poses.
		while ( input.has_another_pose() ) {
			if ( ii == 1 ) {
				input.fill_pose( *start_pose, *rsd_set );
				input_poses.emplace_back( start_pose );
				to_be_analyzed.push_back( *start_pose );
			}

			Pose new_pose;
			input.fill_pose( new_pose, *rsd_set );
			to_be_analyzed.emplace_back( new_pose );

			++ii;
		}

		unconverged_res = analyze_poses_for_convergence( to_be_analyzed );
	} else {

		fill_input_poses_based_on_possible_checkpoints( input_poses, checkpoints_to_pass );

	}

	// setup poses
	Pose seq_rebuild_pose;

	FullModelPoseBuilder builder;
	builder.set_input_poses( input_poses );
	builder.set_options( option );
	builder.initialize_further_from_options();
	builder.build(); // hope this will update original_poses[ 1 ]
	start_pose = input_poses[ 1 ];
	(*scorefxn)(*start_pose);

	// Need to switch after minimizer!
	core::kinematics::FoldTree resample_fold_tree = start_pose->fold_tree();
	core::kinematics::FoldTree emm_ft = start_pose->fold_tree();
	if ( option[ in::file::fold_tree ].user() ) {
		// Especially important in denovo from density cases
		utility::io::izstream foo( option[ in::file::fold_tree ].value() );
		foo >> emm_ft;
		foo.close();
	}

	show_accuracy_report( *start_pose, "Start", 0/*, TR*/ );

	core::Size const nrounds = option[ rounds ].value();

	ErraserMinimizerMover erraser_minimizer;
	erraser_minimizer.scorefxn( scorefxn );
	erraser_minimizer.edens_scorefxn( scorefxn );
	erraser_minimizer.minimize_protein( option[ minimize_protein ] );

	if ( checkpoints_to_pass < 1 ) { // First pass: constrain P, initial FT
		TR << "First pass preminimization..." << std::endl;
		//if ( !utility::file::file_exists( minimized_name.str() ) ) {
		erraser_minimizer.constrain_phosphate( true );
		erraser_minimizer.apply( *start_pose );
		start_pose->dump_pdb( first_preminimized_namer( option[ in::file::s ]()[1]  ) );
	}

	if ( checkpoints_to_pass < 2 ) { // Second pass: constrain P, fun FT
		if ( emm_ft.nres() == start_pose->size() ) start_pose->fold_tree( emm_ft ); // no-op if none specified
		erraser_minimizer.constrain_phosphate( true );
		erraser_minimizer.apply( *start_pose );
		start_pose->dump_pdb( second_preminimized_namer( option[ in::file::s ]()[1]  ) );
		if ( resample_fold_tree.nres() == start_pose->size() ) start_pose->fold_tree( resample_fold_tree ); // no-op if none specified
	}

	/*
	* We only want this to be active when we are adapting models from other methods
	* (like rna_denovo or coarse grained)
	if ( checkpoints_to_pass < 3 ) { // First pass: no constrain P, fun FT
	std::stringstream minimized_name;
	minimized_name << option[ in::file::s ]()[1] << "_preminimized_round_3.pdb";
	//if ( !utility::file::file_exists( minimized_name.str() ) ) {
	start_pose->fold_tree( emm_ft ); // no-op if none specified
	erraser_minimizer.constrain_phosphate( false );
	erraser_minimizer.apply( *start_pose );
	start_pose->dump_pdb( minimized_name.str() );
	start_pose->fold_tree( resample_fold_tree ); // no-op if none specified
	//}
	}
	*/


	// Each round eats up checkpoints 4,5; 6,7; 8,9: or, 3 + 2*ii - 1, 3 + 2*ii
	for ( Size ii = 1; ii <= nrounds; ++ii ) {
		if ( checkpoints_to_pass < 3 + 2 * ii - 1 ) {
			if ( emm_ft.nres() == start_pose->size() ) start_pose->fold_tree( emm_ft ); // no-op if none specified
			erraser_minimizer.apply( *start_pose );
			show_accuracy_report( *start_pose, "Minimized", ii/*, TR*/ );
			start_pose->dump_pdb( round_minimized_namer( option[ in::file::s ]()[1], ii ) );
		}

		// Let's be REALLY SURE
		//auto ft = start_pose.fold_tree();
		//ft.reorder(start_pose.size());
		//start_pose.fold_tree( ft );
		if ( checkpoints_to_pass < 3 + 2 * ii ) {
			if ( resample_fold_tree.nres() == start_pose->size() ) start_pose->fold_tree( resample_fold_tree ); // no-op if none specified
			resample_full_model( *start_pose, scorefxn, unconverged_res );
			show_accuracy_report( *start_pose, "Resampled", ii/*, TR*/ );
			start_pose->dump_pdb( round_resampled_namer( option[ in::file::s ]()[1], ii ) );
		}
	}

	if ( emm_ft.nres() == start_pose->size() ) start_pose->fold_tree( emm_ft ); // no-op if none specified
	erraser_minimizer.apply( *start_pose );
	show_accuracy_report( *start_pose, "FINAL", 0/*, TR*/ );

	std::stringstream filename_stream;
	filename_stream << option[ in::file::s ]()[1] << "FINISHED.pdb";
	start_pose->dump_pdb( filename_stream.str() );
}

bool atoms_have_bond_to_bonded_atoms( pose::Pose const & pose, Size const ai, Size const ii, Size const aj, Size const jj ) {
	// Loop through atoms in both residues
	// Accumulate set of atom IDs bonded to each atom.
	// Check if there are bonds between any pairs.

	utility::vector1< core::id::AtomID > bonded_to_ai, bonded_to_aj;

	for ( Size ak = 1; ak <= pose.residue( ii ).natoms(); ++ak ) {
		if ( pose.conformation().atoms_are_bonded( core::id::AtomID( ai, ii ), core::id::AtomID( ak, ii ) ) ) {
			bonded_to_ai.emplace_back( ak, ii );
		}
		if (  pose.conformation().atoms_are_bonded( core::id::AtomID( aj, jj ), core::id::AtomID( ak, ii ) ) ) {
			bonded_to_aj.emplace_back( ak, ii );
		}
	}
	for ( Size ak = 1; ak <= pose.residue( jj ).natoms(); ++ak ) {
		if ( pose.conformation().atoms_are_bonded( core::id::AtomID( ai, ii ), core::id::AtomID( ak, jj ) ) ) {
			bonded_to_ai.emplace_back( ak, jj );
		}
		if ( pose.conformation().atoms_are_bonded( core::id::AtomID( aj, jj ), core::id::AtomID( ak, jj ) ) ) {
			bonded_to_aj.emplace_back( ak, jj );
		}
	}

	for ( auto const & b_to_ai : bonded_to_ai ) {
		for ( auto const & b_to_aj : bonded_to_aj ) {
			if ( pose.conformation().atoms_are_bonded( b_to_ai, b_to_aj ) ) return true;
		}
	}

	return false;
}

bool atoms_have_mutual_bond_to_atom( pose::Pose const & pose, Size const ai, Size const ii, Size const aj, Size const jj ) {
	// Loop through atoms in both residues
	// Assumes no 1-atom polymeric residues. That's ok.
	for ( Size ak = 1; ak <= pose.residue( ii ).natoms(); ++ak ) {
		if ( pose.conformation().atoms_are_bonded( core::id::AtomID( ai, ii ), core::id::AtomID( ak, ii ) )
				&& pose.conformation().atoms_are_bonded( core::id::AtomID( aj, jj ), core::id::AtomID( ak, ii ) ) ) return true;
	}
	for ( Size ak = 1; ak <= pose.residue( jj ).natoms(); ++ak ) {
		if ( pose.conformation().atoms_are_bonded( core::id::AtomID( ai, ii ), core::id::AtomID( ak, jj ) )
				&& pose.conformation().atoms_are_bonded( core::id::AtomID( aj, jj ), core::id::AtomID( ak, jj ) ) ) return true;
	}

	return false;
}

bool
bump_check( pose::Pose const & pose, Size const ii ) {
	for ( Size jj = 1; jj <= pose.size(); ++jj ) {
		if ( jj == ii ) continue;
		if ( jj == ii + 1 ) continue;
		if ( jj == ii - 1 ) continue;

		// If first atoms are more than 10A apart no clashes
		if ( pose.residue( ii ).xyz( 1 ).distance_squared( pose.residue( jj ).xyz( 1 ) ) > 100 ) continue;

		for ( Size ai = 1; ai <= pose.residue( ii ).natoms(); ++ai ) {
			for ( Size aj = 1; aj <= pose.residue( jj ).natoms(); ++aj ) {
				// If the two atoms are bonded ignore.
				if ( pose.conformation().atoms_are_bonded( core::id::AtomID( ai, ii ), core::id::AtomID( aj, jj ) ) ) continue;
				if ( atoms_have_mutual_bond_to_atom( pose, ai, ii, aj, jj ) ) continue;
				if ( atoms_have_bond_to_bonded_atoms( pose, ai, ii, aj, jj ) ) continue;

				Real ai_r = pose.residue_type( ii ).atom_type( ai ).lj_radius();
				Real aj_r = pose.residue_type( jj ).atom_type( aj ).lj_radius();

				if ( pose.residue( ii ).xyz( ai ).distance_squared( pose.residue( jj ).xyz( aj ) ) <
						( aj_r + ai_r - 0.35 ) * ( aj_r + ai_r - 0.35 ) ) {
					TR.Trace << "Bump check failed because residue " << ii << " atom " << pose.residue( ii ).atom_name( ai ) << " (lj_radius " << ai_r << ") is within "
						<< ( aj_r + ai_r - 0.35 ) * ( aj_r + ai_r - 0.35 ) << " of residue " << jj << " atom " << pose.residue( jj ).atom_name( aj ) << " (lj_radius " << aj_r << ")" << std::endl;
					return true;
				}
			}
		}
	}
	return false;
}

void resample_full_model( pose::Pose & start_pose, ScoreFunctionOP const & scorefxn, utility::vector1< Size > const & definite_residues ) {

	/*// Remove any single-residue checkpoint indicators if present.
	#include <glob.h>
	glob_t *pglob;
	glob("seq_rebuild_temp_*pdb", );
	seq_rebuild_temp_A:14.pdb
	utility::delete_file( aaa)
	*/
	RNA_SuiteName suite_namer;

	// Rebuilds to be based on torsion score.
	// (For now... later, do an actual PHENIX integration or something.)
	std::set< Size > residues_to_rebuild( definite_residues.begin(), definite_residues.end() );
	for ( Size ii = 1; ii <= start_pose.size(); ++ii ) {
		if ( !option[ minimize_protein ].value() && start_pose.residue_type( ii ).is_protein() ) continue;

		// This should pick up bad torsions
		Real const torsion_score = start_pose.energies().residue_total_energies( ii )[ rna_torsion ];
		if ( torsion_score > 1.2 ) {
			TR << "Rebuilding residue " << ii <<" for a bad torsion score of " << torsion_score << std::endl;
			residues_to_rebuild.insert( ii );
		}

		// Bad suites, if it hasn't been added already
		// TODO: chain aware
		auto suite_assignment = suite_namer.assign( start_pose, ii );
		if ( suite_assignment.name == "!!" ) {
			TR << "Rebuilding residue " << ii << " with its neighbors as possible for a bad suite" << std::endl;
			if ( ii > 1 ) {
				residues_to_rebuild.insert( ii - 1 );
			}
			residues_to_rebuild.insert( ii );
			if ( ii < start_pose.size() ) {
				residues_to_rebuild.insert( ii + 1 );
			}
		}

		// This should pick up bad bond lengths and angles.
		// (If cart_bonded is on...)
		Real const cart = start_pose.energies().residue_total_energies( ii )[ cart_bonded ];
		if ( cart > 2 ) {
			TR << "Rebuilding residue " << ii <<" for a bad cart_bonded score of " << cart << std::endl;
			residues_to_rebuild.insert( ii );
		}


		// Clashes
		if ( bump_check( start_pose, ii ) ) {
			TR << "Rebuilding residue " << ii <<" for a bad bump check" << std::endl;
			residues_to_rebuild.insert( ii );
		}

		// Make sure bad glycosylation sites get picked up here. Need samplers.

		// Make sure branched nucleotides are always resampled.

		// Make sure ions get refined... have to think about how to handle them, water coordination, etc.

		// (at first just Mg?)

		// Make sure cyclic dinucleotide ligands are redocked. Also JUMP_DOCK.
		if ( !start_pose.residue( ii ).is_polymer() ) residues_to_rebuild.insert( ii );
	}
	TR << "Total to resample this round: " << residues_to_rebuild << std::endl;

	// OK, what residues might these be? What might be connected?
	for ( auto const elem : residues_to_rebuild ) {
		TR << elem << ": " << " " << start_pose.pdb_info()->chain( elem ) << start_pose.pdb_info()->number( elem ) << " " << start_pose.residue_type( elem ).name() << std::endl;
	}
	//return;

	using namespace protocols::stepwise::monte_carlo;
	pose::Pose seq_rebuild_pose = start_pose;
	protocols::viewer::add_conformation_viewer(seq_rebuild_pose.conformation(), "current", 500, 500);


	// initialize options
	options::StepWiseMonteCarloOptionsOP options( new options::StepWiseMonteCarloOptions );
	options->initialize_from_command_line();
	options->set_erraser( true );
	options->set_enumerate( true );
	options->set_skip_deletions( true );
	options->set_output_minimized_pose_list( false );
	options->set_force_moving_res_for_erraser( true );

	// run StepWiseMasterMover::resample_full_model
	mover::StepWiseMasterMover master_mover( scorefxn, options );
	// This is essential for rmsd_screen apparently
	master_mover.set_native_pose( pose::PoseCOP( new pose::Pose( start_pose ) ) );

	//master_mover.resample_full_model( start_pose, seq_rebuild_pose, true /*checkpointing_breadcrumbs*/ );

	//return;
	master_mover.resample_full_model( start_pose, seq_rebuild_pose, true /*checkpointing_breadcrumbs*/, utility::vector1< Size >( residues_to_rebuild.begin(), residues_to_rebuild.end() ) );
	// DEBUG 82
	//master_mover.resample_full_model( start_pose, seq_rebuild_pose, true /*checkpointing_breadcrumbs*/, utility::vector1< Size >( 1, 82 ) );
	//master_mover.resample_full_model( start_pose, seq_rebuild_pose, true /*checkpointing_breadcrumbs*/, utility::vector1< Size >( 1, 153 ) );

	// score seq_rebuild pose
	(*scorefxn)(seq_rebuild_pose);

	// show scores of start_pose and full_model_pose
	//if ( option[ show_scores ]() ) {
	TR << "\n Score before seq_rebuild:" << std::endl;
	scorefxn->show( std::cout, start_pose );
	TR << "\n Score after seq_rebuild:" << std::endl;
	scorefxn->show( std::cout, seq_rebuild_pose );
	//}

	// dump pdb -- should figure out file name based on inputs
	seq_rebuild_pose.dump_pdb( "SEQ_REBUILD.pdb" );

	// write out to silent file
	//std::string tag = "S_0";
	//stepwise::monte_carlo::output_to_silent_file( tag, silent_file_out, full_model_pose );
	start_pose = seq_rebuild_pose;
	protocols::viewer::clear_conformation_viewers();
}


///////////////////////////////////////////////////////////////
void*
my_main( void* )
{
	erraser2_test();
	exit( 0 );
}

///////////////////////////////////////////////////////////////////////////////
int
main( int argc, char * argv [] )
{

	try {

		NEW_OPT( rounds, "Number of total rounds", 3 );
		NEW_OPT( minimize_protein, "Minimize protein", false );
		// help
		std::cout << std::endl << "Basic usage:  " << argv[0] << "  -in::file::silent <silent file>" << std::endl;
		std::cout << std::endl << " Type -help for full slate of options." << std::endl << std::endl;

		// options

		// setup
		devel::init(argc, argv);

		// run
		protocols::viewer::viewer_main( my_main );
		exit( 0 );

	} catch ( utility::excn::Exception const & e ) {
		std::cout << "caught exception " << e.msg() << std::endl;
		return -1;
	}

}

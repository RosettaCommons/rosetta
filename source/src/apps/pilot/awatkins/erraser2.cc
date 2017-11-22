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
#include <protocols/stepwise/setup/FullModelInfoSetupFromCommandLine.hh>


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
using namespace core::scoring;
using namespace protocols::rna::movers;
using namespace core::pose::rna;


static basic::Tracer TR( "apps.pilot.awatkins.erraser2" );

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
	ResidueTypeSetCOP rsd_set;
	rsd_set = ChemicalManager::get_instance()->residue_type_set( FA_STANDARD );

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
		PDBPoseInputStream input( option[ in::file::s ]() );
		// iterate over input stream
		input.fill_pose( *start_pose, *rsd_set );
		input_poses.emplace_back( start_pose );
	}

	// setup poses
	Pose seq_rebuild_pose;

	FullModelPoseBuilder builder;
	builder.set_input_poses( input_poses );
	builder.set_options( option );
	builder.initialize_further_from_options();
	builder.build(); // hope this will update original_poses[ 1 ]
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

	core::Size const nrounds = 1;

	ErraserMinimizerMover erraser_minimizer;
	erraser_minimizer.scorefxn( scorefxn );
	erraser_minimizer.edens_scorefxn( scorefxn );

	for ( Size ii = 1; ii <= nrounds; ++ii ) {
		std::stringstream filename_stream;
		filename_stream << "minimized_" << ii << ".pdb";
		std::string const minimized_name = filename_stream.str();
		filename_stream.str( "" );
		filename_stream << "resampled_" << ii << ".pdb";
		std::string const resampled_name = filename_stream.str();

		start_pose->fold_tree( emm_ft ); // no-op if none specified
		erraser_minimizer.apply( *start_pose );
		show_accuracy_report( *start_pose, "Minimized", ii/*, TR*/ );
		start_pose->dump_pdb( minimized_name );

		// Let's be REALLY SURE
		//auto ft = start_pose.fold_tree();
		//ft.reorder(start_pose.size());
		//start_pose.fold_tree( ft );
		start_pose->fold_tree( resample_fold_tree ); // no-op if none specified
		resample_full_model( *start_pose, scorefxn, unconverged_res );
		show_accuracy_report( *start_pose, "Resampled", ii/*, TR*/ );
		start_pose->dump_pdb( resampled_name );
	}

	start_pose->fold_tree( emm_ft ); // no-op if none specified
	erraser_minimizer.apply( *start_pose );
	show_accuracy_report( *start_pose, "FINAL", 0/*, TR*/ );
	start_pose->dump_pdb( "FINISHED.pdb" );
}

void resample_full_model( pose::Pose & start_pose, ScoreFunctionOP const & scorefxn, utility::vector1< Size > const & definite_residues ) {

	RNA_SuiteName suite_namer;

	// Rebuilds to be based on torsion score.
	// (For now... later, do an actual PHENIX integration or something.)
	std::set< Size > residues_to_rebuild( definite_residues.begin(), definite_residues.end() );
	for ( Size ii = 1; ii <= start_pose.size(); ++ii ) {
		// This should pick up bad torsions
		Real const torsion_score = start_pose.energies().residue_total_energies( ii )[ rna_torsion ];
		TR.Trace << "Torsion score " << ii << " " << torsion_score << std::endl;
		if ( torsion_score > 1.2 ) residues_to_rebuild.insert( ii );

		// Bad suites, if it hasn't been added already
		auto suite_assignment = suite_namer.assign( start_pose, ii );
		if ( suite_assignment.name == "!!" )  residues_to_rebuild.insert( ii );

		// This should pick up bad bond lengths and angles.
		// (If cart_bonded is on...)
		Real const cart = start_pose.energies().residue_total_energies( ii )[ cart_bonded ];
		TR.Trace << "Cart-bonded score " << ii << " " << cart << std::endl;
		if ( cart > 12 ) residues_to_rebuild.insert( ii );

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
	//master_mover.resample_full_model( start_pose, seq_rebuild_pose, true /*checkpointing_breadcrumbs*/ );
	// debug -- residues with odd jump behavior in use case
	residues_to_rebuild.erase(131);

	TR << "Post deletion: " << std::endl;

	// OK, what residues might these be? What might be connected?
	for ( auto const elem : residues_to_rebuild ) {
		TR << elem << ": " << " " << start_pose.pdb_info()->chain( elem ) << start_pose.pdb_info()->number( elem ) << " " << start_pose.residue_type( elem ).name() << std::endl;
	}
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

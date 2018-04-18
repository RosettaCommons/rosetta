// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file apps/pilot/kkappel/score_rnp.cc
/// @brief Align RNA structure using command line input residues
/// @author Kalli Kappel kappel@stanford.edu


//Unit Headers
//Package Headers
//Project Headers
#include <core/pose/Pose.hh>
#include <core/pose/util.hh>
#include <core/import_pose/import_pose.hh>
#include <core/import_pose/pose_stream/PoseInputStream.hh>
#include <core/import_pose/pose_stream/PoseInputStream.fwd.hh>
#include <core/import_pose/pose_stream/PDBPoseInputStream.hh>
#include <core/import_pose/pose_stream/SilentFilePoseInputStream.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/Energies.hh>
#include <core/chemical/ChemicalManager.hh>
#include <core/io/silent/BinarySilentStruct.hh>
#include <protocols/rigid/RigidBodyMover.hh>

//Utility Headers
#include <basic/options/option.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <utility/pointer/owning_ptr.hh>
#include <basic/Tracer.hh>
//Numeric Headers
#include <numeric/random/random.hh>
//C++ Headers
#include <iostream>
#include <fstream>

#include <devel/init.hh>

#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <basic/options/keys/out.OptionKeys.gen.hh>
#include <basic/options/keys/score.OptionKeys.gen.hh>
#include <basic/options/option_macros.hh>

// New options for this application
using namespace basic::options::OptionKeys;

static basic::Tracer TR( "apps.pilot.kkappel.move_jump" );

/////////////////////////////////////////////////////////////////////////////////
void score_rnp() {
	using namespace core::scoring;
	using namespace core::chemical;
	using namespace basic::options;
	using namespace core::io::silent;
	using namespace core::import_pose::pose_stream;
	using namespace basic::options::OptionKeys;
	core::pose::Pose pose;

	// Setup the score function, it will default to flicker rnp sfxn
	ScoreFunctionOP sfxn;
	if ( option[ OptionKeys::score::weights ].user() ) {
		sfxn = get_score_function();
	} else { //default to flicker
		sfxn = ScoreFunctionFactory::create_score_function( "flicker2015_sol2_fa_elec3.wts" );
	}

	core::Real min_score( 100000000.0 );

	core::Real partition_function( 0.0 );
	core::Real kT( 0.59 ); //Boltzmann factor room temp

	// Now get the input poses and show their scores
	PoseInputStreamOP input;
	if ( option[ in::file::silent ].user() ) {
		input = PoseInputStreamOP( new SilentFilePoseInputStream( option[ in::file::silent ]() ) );
	} else {
		input = PoseInputStreamOP( new PDBPoseInputStream( option[ in::file::s ]() ) );
	}

	ResidueTypeSetCOP rsd_set;
	rsd_set = core::chemical::ChemicalManager::get_instance()->residue_type_set( FA_STANDARD );


	// utility::vector1< std::string > input_structs = option[ OptionKeys::in::file::s ]();
	// for ( core::Size i = 1; i <= input_structs.size(); ++i ) {
	while ( input->has_another_pose() ) {
		//core::import_pose::pose_from_pdb( pose, input_structs[ i ]);
		input->fill_pose( pose, *rsd_set );
		// super hacky, change the fold tree
		core::kinematics::FoldTree separation_foldtree;
		separation_foldtree.add_edge(1,258, core::kinematics::Edge::PEPTIDE );
		separation_foldtree.add_edge(259,277, core::kinematics::Edge::PEPTIDE );
		separation_foldtree.add_edge(1,259, 1);
		separation_foldtree.reorder( 1 );
		pose.fold_tree( separation_foldtree );
		// Do a super hacky move of the jump
		protocols::rigid::RigidBodyTransMoverOP separate_complex( new protocols::rigid::RigidBodyTransMover( pose, 1 /* jump # */) );
		separate_complex->step_size( -1000.0 );
		separate_complex->apply( pose );
		pose.dump_pdb( "aligned.pdb" );
		//pose.dump_pdb( 'aligned' + option[ in::file::s ]() );

		//TR << "Score for " << core::pose::tag_from_pose( pose ) << ":" << std::endl;
		////TR << "Score for " << input_structs[ i ] << ":" << std::endl;
		//sfxn->show( pose );
		//core::Real score = (*sfxn)(pose);
		//if ( (*sfxn)(pose) < min_score ) {
		// min_score = score;
		//}
		//partition_function += std::exp( -1.0 * score / kT );

	}

	// core::Real boltzmann_wtd_avg = -1.0 * kT * std::log( partition_function );

	//std::cout << "Minimum score: " << min_score << std::endl;
	//std::cout << "Boltzmann weighted average: " << boltzmann_wtd_avg << std::endl;
}

/////////////////////////////////////////////////////////////////////////////////
int main( int argc, char ** argv ) {

	try {
		using namespace basic::options;

		devel::init( argc, argv );
		score_rnp();
	} catch ( utility::excn::EXCN_Base const & e ) {
		std::cout << "caught exception " << e.msg() << std::endl;
		return -1;
	}
	return 0;
}

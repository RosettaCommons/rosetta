// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   test_glycan_linkages.cc
/// @brief  This is simple pilot app for testing carbohydrates.
/// @author Labonte <JWLabonte@jhu.edu>


// Project headers
#include <devel/init.hh>

#include <core/types.hh>
#include <core/chemical/ChemicalManager.hh>
#include <core/chemical/ResidueTypeSet.hh>
#include <core/pose/Pose.hh>
#include <core/pose/annotated_sequence.hh>
#include <core/pose/carbohydrates/util.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/pack/task/PackerTask.hh>
#include <core/pack/task/TaskFactory.hh>

#include <protocols/minimization_packing/PackRotamersMover.hh>
#include <protocols/minimization_packing/RotamerTrialsMover.hh>

// Utility headers
#include <utility/excn/Exceptions.hh>

// C++ header
#include <iostream>


// Constants
int const SUCCESS( 0 );
int const FAILURE( -1 );

std::string const PATH = "/home/labonte/Workspace/Carbohydrates/testing/";


int
main( int argc, char *argv[] )
{
	using namespace std;

	try {
		using namespace core;
		using namespace chemical;
		using namespace pose;
		//using namespace scoring;
		//using namespace pack::task;
		//using namespace protocols::minimization_packing;

		// Initialize core.
		devel::init( argc, argv );

		// Declare variables.
		Pose starting_pose, pose;
		ResidueTypeSetCOP residue_set( ChemicalManager::get_instance()->residue_type_set( "fa_standard" ) );


		// Create initial test pose.
		make_pose_from_sequence( starting_pose, "ANA", *residue_set );

		cout << "Created pose (peptide only): " << starting_pose << endl;

		// Extend peptide.
		for ( core::uint i( 1 ); i <=2 ; ++i ) {
			starting_pose.set_omega( i, 180.0 );
		}

		cout << "Extended peptide." << endl;

		pose::carbohydrates::glycosylate_pose( starting_pose, 2,
			"b-D-GlcpNAc-(1->4)-[a-L-Fucp-(1->6)]-b-D-GlcpNAc-" );

		cout << "Glycosylated pose: " << starting_pose << endl;

		// Idealize glycan.
		//starting_pose.set_phi( 6, -75.0 );
		//starting_pose.set_phi( 7, -75.0 );
		//starting_pose.set_psi( 7, 115.0 );
		//starting_pose.set_phi( 8, -75.0 );
		//starting_pose.set_psi( 8, 115.0 );
		//starting_pose.set_omega( 8, 180.0 );

		cout << "Idealized glycan." << endl;

		starting_pose.dump_pdb( PATH + "N-linked_test.start.pdb" );

		cout << "Output starting pose." << endl;

		/*
		// Set up ScoreFunction.
		ScoreFunctionOP sf( get_score_function() );

		// Set up PackerTask.
		PackerTaskOP task( TaskFactory::create_packer_task( starting_pose ) );
		task->restrict_to_repacking();

		// Set up two test Movers.
		PackRotamersMover pack_rotamers_mover( sf, task );
		RotamerTrialsMover rotamer_trials_mover( sf, *task );

		pose = starting_pose;
		pack_rotamers_mover.apply( pose );
		pose.dump_pdb( PATH + "N-linked_test.prm.pdb" );

		cout << "Output pose packed with PackRotamersMover." << endl;

		pose = starting_pose;
		rotamer_trials_mover.apply( pose );
		pose.dump_pdb( PATH + "N-linked_test.rtm.pdb" );

		cout << "Output pose packed with RotamerTrialsMover." << endl;
		*/

	} catch ( utility::excn::Exception const & e ) {
		cerr << "Caught exception: " << e.msg() << endl;
		return FAILURE;
	}
	return SUCCESS;
}

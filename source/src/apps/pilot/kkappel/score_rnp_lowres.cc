// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file apps/pilot/kkappel/calculate_distances.cc
/// @brief score RNA/protein complex with low-res score function
/// @author Kalli Kappel kappel@stanford.edu


//Unit Headers
//Package Headers
//Project Headers
#include <core/pose/Pose.hh>
#include <core/chemical/AtomType.hh>
#include <core/chemical/ChemicalManager.hh>
#include <core/import_pose/pose_stream/PoseInputStream.hh>
#include <core/import_pose/pose_stream/PoseInputStream.fwd.hh>
#include <core/import_pose/pose_stream/SilentFilePoseInputStream.hh>
#include <core/import_pose/pose_stream/PDBPoseInputStream.hh>
#include <utility/vector1.hh>
#include <core/pose/extra_pose_info_util.hh>
#include <core/pose/rna/util.hh>
#include <core/import_pose/import_pose.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/Energies.hh>
#include <core/scoring/dssp/Dssp.hh>
#include <core/pack/task/PackerTask.hh>
#include <core/pack/task/TaskFactory.hh>
#include <core/pack/pack_rotamers.hh>
#include <core/kinematics/MoveMap.hh>
#include <core/optimization/MinimizerOptions.hh>
#include <core/optimization/AtomTreeMinimizer.hh>
#include <core/io/silent/BinarySilentStruct.hh>
#include <core/io/silent/SilentFileData.hh>
#include <core/io/silent/SilentFileOptions.hh>
#include <core/pose/full_model_info/FullModelInfo.hh>
#include <core/conformation/ResidueFactory.hh>
#include <core/util/SwitchResidueTypeSet.hh>
//Utility Headers
#include <basic/options/option.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>
#include <utility/pointer/owning_ptr.hh>
#include <basic/Tracer.hh>
//Numeric Headers
#include <numeric/random/random.hh>
#include <numeric/xyzMatrix.hh>
//C++ Headers
#include <iostream>
#include <fstream>

#include <devel/init.hh>
#include <protocols/moves/MonteCarlo.hh>
#include <protocols/moves/PyMOLMover.hh>
#include <protocols/jd2/JobDistributor.hh>
#include <protocols/rigid/RigidBodyMover.hh>

#include <core/scoring/constraints/ConstraintIO.hh>
#include <core/scoring/constraints/ConstraintSet.hh>
#include <core/scoring/constraints/AtomPairConstraint.hh>

#include <basic/options/keys/out.OptionKeys.gen.hh>
#include <basic/options/keys/score.OptionKeys.gen.hh>
#include <basic/options/keys/constraints.OptionKeys.gen.hh>
#include <basic/options/option_macros.hh>

// New options for this application
using namespace basic::options::OptionKeys;

OPT_KEY( Boolean, dump )
OPT_KEY( Boolean, low_res )

static basic::Tracer TR( "apps.pilot.kkappel.score_rnp_lowres" );

void check_structures() {

	using namespace core::import_pose::pose_stream;
	using namespace basic::options;
	using namespace basic::options::OptionKeys;
	using namespace core::chemical;
	using namespace core::io::silent;
	using namespace core::scoring;

	// input stream
	PoseInputStreamOP input;
	if ( option[ in::file::silent ].user() ) {
		if ( option[ in::file::tags ].user() ) {
			input = PoseInputStreamOP( new SilentFilePoseInputStream(
				option[ in::file::silent ](),
				option[ in::file::tags ]()
				) );
		} else {
			input = PoseInputStreamOP( new SilentFilePoseInputStream( option[ in::file::silent ]() ) );
		}
	} else {
		input = PoseInputStreamOP( new PDBPoseInputStream( option[ in::file::s ]() ) );
	}

	TR << TR.Blue << "Done loading the pose." << std::endl;

	ResidueTypeSetCOP rsd_set;
	rsd_set = core::chemical::ChemicalManager::get_instance()->residue_type_set( FA_STANDARD );
	core::pose::Pose pose;

	//core::scoring::ScoreFunctionOP sfxn;
	//sfxn = core::scoring::get_score_function();


	core::scoring::ScoreFunctionOP sfxn;
	//sfxn = core::scoring::get_score_function( "rnp_lowres" );
	//sfxn = core::scoring::ScoreFunctionFactory::create_score_function( "rnp_lores" );
	//sfxn = core::scoring::get_score_function();
	if ( option[ score::weights ].user() ) {
		sfxn = get_score_function();
	} else {
		sfxn = core::scoring::ScoreFunctionFactory::create_score_function( "rnp_lores" );
	}


	// Silent file output setup
	std::string const silent_file = option[ out::file::silent  ]();
	SilentFileOptions opts; // initialized from the command line
	SilentFileData silent_file_data(opts);

	while ( input->has_another_pose() ) {

		input->fill_pose( pose, *rsd_set );

		//std::cout << "FOLD TREE: " << std::endl;
		//pose.fold_tree().show( std::cout );
		//pose.dump_pdb( "pose_dump.pdb" );

		// tag
		std::string tag = tag_from_pose( pose );

		if ( option[ low_res ]() ) {
			core::util::switch_to_residue_type_set( pose, core::chemical::CENTROID, false /* no sloppy match */, true /* only switch protein residues */, true /* keep energies! */ );
		}

		if ( basic::options::option[ dump ]() ) {
			std::string out_pdb_file = tag + "_dump.pdb";
			pose.dump_pdb( out_pdb_file );
		}

		// score it
		(*sfxn)( pose );

		// write it to the silent file
		BinarySilentStruct s( opts, pose, tag );
		silent_file_data.write_silent_struct( s, silent_file, true /*write score only*/ );
	}

}

int main( int argc, char ** argv ) {

	try {
		using namespace basic::options;
		option.add_relevant( score::weights );
		NEW_OPT( dump, "dump the final structure", false );
		NEW_OPT( low_res, "switch protein to centroid", true );

		devel::init( argc, argv );
		check_structures();
	} catch (utility::excn::Exception const & e ) {
		std::cout << "caught exception " << e.msg() << std::endl;
		return -1;
	}
	return 0;
}

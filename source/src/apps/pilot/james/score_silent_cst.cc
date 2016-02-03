// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington UW TechTransfer, email: license@u.washington.edu.

/// @file
/// @brief
/// @author James Thompson

// libRosetta headers


#include <core/types.hh>

#include <core/chemical/ChemicalManager.hh>
#include <core/chemical/ResidueTypeSet.fwd.hh>
#include <core/scoring/rms_util.hh>
#include <core/scoring/Energies.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/ScoreFunction.fwd.hh>
#include <core/scoring/constraints/util.hh>
#include <core/scoring/constraints/ConstraintIO.hh>
#include <core/scoring/constraints/ConstraintSet.hh>
#include <core/scoring/constraints/ConstraintSet.fwd.hh>

#include <core/pose/Pose.hh>
#include <core/pose/Pose.fwd.hh>

#include <devel/init.hh>

#include <basic/options/option.hh>
#include <core/io/pdb/pdb_writer.hh>

#include <core/io/silent/silent.fwd.hh>
#include <core/io/silent/ProteinSilentStruct.hh>
#include <core/io/silent/SilentFileData.hh>

#include <protocols/simple_filters/RmsdEvaluator.hh>

#include <utility/vector1.hh>
#include <utility/io/ozstream.hh>

// C++ headers
#include <fstream>
#include <iostream>
#include <string>

// option key includes

#include <basic/options/keys/out.OptionKeys.gen.hh>
#include <basic/options/keys/in.OptionKeys.gen.hh>

//Auto Headers
#include <core/import_pose/import_pose.hh>

#include <utility/excn/Exceptions.hh>


int
main( int argc, char* argv[] )
{
	try {

	// options, random initialization
	devel::init( argc, argv );

	using namespace core::scoring;
	using namespace basic::options;
	using namespace basic::options::OptionKeys;
	using namespace core::io::silent;
	using std::string;
	using utility::vector1;

	// setup residue types
	core::chemical::ResidueTypeSetCAP rsd_set;
	rsd_set = core::chemical::ChemicalManager::get_instance()->residue_type_set(
		option[ in::file::residue_type_set ]()
	);

	// configure score function
	core::scoring::ScoreFunctionOP scorefxn = core::scoring::get_score_function();
	core::scoring::EnergyMap emap;
	emap.zero();

	core::io::silent::SilentFileData sfd;
	core::scoring::constraints::ConstraintSetOP cstset = NULL;

	std::string infile	= *(option[ in::file::silent ]().begin());
	std::string outfile = option[ out::file::silent ]();
	utility::io::ozstream output;

	if ( option[ in::file::silent ].user() ) {
		sfd.read_file( infile );
	}

	core::pose::Pose native_pose, pose;
	FArray2D< core::Real > native_points, decoy_points;
	utility::vector1< bool > subset; // positions with constraints
	if ( option[ in::file::native ].user() ) {
		// read in pdb and constraints if necessary
		core::import_pose::pose_from_file(
			native_pose, *rsd_set, option[ in::file::native ]()
		);
		if ( !cstset ) {
			cstset = core::scoring::constraints::ConstraintIO::read_constraints(
				core::scoring::constraints::get_cst_file_option(),
				new core::scoring::constraints::ConstraintSet,
			 	native_pose
			);
			native_pose.constraint_set( cstset );
		}

		subset.reserve( native_pose.total_residue() );
		scorefxn->set_weight( atom_pair_constraint, 1.0 );
		// score and output native score.
		core::Real score = (*scorefxn)(native_pose);
		emap = native_pose.energies().total_energies();
		core::Real cst_score = emap[ core::scoring::atom_pair_constraint ];
		score = score - cst_score; // substract out cst_score!
		emap.zero();
		//score = 0; // for now.
		output.open( outfile.c_str() );
		//native_points = get_subset_CA_coords( native_pose, subset );

		output << "tag alignment score cst_score rmsd maxsub gdtmm" << "\n";
		output 	<< "native native " << score << ' ' << cst_score << ' ' << 0.0 << ' '
						<< native_pose.total_residue() << ' ' << 1.0 << "\n";

		//output << "tag score cst_score rmsd maxsub gdtmm" << "\n";
		//output 	<< "native " << score << ' ' << cst_score << ' ' << 0.0 << ' '
		//				<< native_pose.total_residue() << ' ' << 1.0 << "\n";
		//std::cerr << "tag alignment score cst_score rmsd maxsub gdtmm" << "\n";
		//std::cerr << "native native " << score << ' ' << cst_score << ' ' << 0.0 << ' '
		//				<< native_pose.total_residue() << ' ' << 1.0 << "\n";
		std::cerr << "native cst_score = " << cst_score << std::endl;
	}
	core::io::silent::SilentFileData sfd_out;

	vector1< string > user_tags;
	if ( option[ in::file::user_tags ].user() ) {
		user_tags = option[ in::file::user_tags ]();
	}

	for ( SilentFileData::iterator iter = sfd.begin(), end = sfd.end();
				iter != end; ++iter
	) {
		std::string alignment_id = "empty";
		if ( user_tags.size() > 0 ) {
			std::string temp_id( iter->get_comment( "user_tag" ) );
			if ( temp_id == "" ) temp_id = iter->get_comment( "user_ta" );
			if ( temp_id == "" ) continue;
			typedef vector1< string >::const_iterator iter;
			bool valid_tag( false );
			for ( iter it = user_tags.begin(), end = user_tags.end();
						it != end && !valid_tag; ++it
			) {
				if ( temp_id.find( *it ) != std::string::npos ) valid_tag = true;
				alignment_id = temp_id;
			} // for user_tags
			if ( !valid_tag ) continue;
		} // user_tags.size() > 0

		iter->fill_pose( pose, *rsd_set );

		if ( !cstset ) {
			cstset = core::scoring::constraints::ConstraintIO::read_constraints(
				core::scoring::constraints::get_cst_file_option(),
				new core::scoring::constraints::ConstraintSet,
			 	pose
			);
		}
		pose.constraint_set( cstset );

		core::Real score = (*scorefxn)(pose);
		core::io::silent::ProteinSilentStruct ss(
			pose, iter->decoy_tag(), option[ in::file::fullatom]()
		);
		emap.zero();
		emap = pose.energies().total_energies();

		core::Real cst_score = emap[ atom_pair_constraint ];
		score = score - cst_score;
		if ( option[ in::file::native ].user() ) {
			core::Real rmsd	 = protocols::simple_filters::native_CA_rmsd( pose, native_pose );
			core::Real gdtmm = core::scoring::CA_gdtmm ( pose, native_pose );
			int maxsub_nres	 = core::scoring::CA_maxsub( pose, native_pose );
			output 	<< iter->decoy_tag()
							<< ' ' << alignment_id
							<< ' ' << score
							<< ' ' << cst_score
							<< ' ' << rmsd
							<< ' ' << maxsub_nres
							<< ' ' << gdtmm << std::endl;
			//std::cerr 	<< iter->decoy_tag()
			//				//<< ' ' << alignment_id
			//				<< ' ' << score
			//				<< ' ' << cst_score
			//				<< ' ' << rmsd
			//				<< ' ' << maxsub_nres
			//				<< ' ' << gdtmm << std::endl;
		} else {
			//output << iter->decoy_tag() << ' ' << score << ' ' << cst_score << "\n";
			//core::io::silent::SilentStructOP pss(
			//	new core::io::silent::ScoreFileSilentStruct( pose )
			//);
			//pss->add_energy( "atompair_constraint", cst_score );
			//pss->decoy_tag( iter->decoy_tag() );
			//sfd_out.write_silent_struct( *pss, outfile, true );
			iter->add_energy( "atompair_constraint", cst_score );
			//sfd_out.write_silent_struct( *(*iter), outfile, true );
		}
	} // for ss in SilentFileData
	output.close();

	} catch ( utility::excn::EXCN_Base const & e ) {
		std::cout << "caught exception " << e.msg() << std::endl;
		return -1;
	}

	return 0;
} // main
